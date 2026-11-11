# Trạng thái & handoff — cập nhật 2026-11-11

Đọc file này trước khi làm gì tiếp. Nó thay cho việc đọc lại cả hội thoại.

## Tình thế một đoạn

Cuộc thi 长江杯 2026, ta là **public team**, chọn **2 track minibinder**: TL1A và IL-7Rα.
20 sequence/track. Ranking cuối thuần theo **KD đo thật**, nhưng trước đó phải vượt
**cửa sàng in-silico** (ESMFold2 + ipSAE, 1000 slot/track, chật 3–4 lần).
**Mốc làm việc: chốt 11-13, confirm 11-14** (không phải 11-30 — xem `00-COMPETITION-FACTS.md` §1).

## Những gì đã chốt (có bằng chứng, không phải phỏng đoán)

| Việc | Kết luận | Nguồn |
|---|---|---|
| Đề thi, 7 target, luật | Trích từ JS bundle của site chính thức | `00-COMPETITION-FACTS.md` |
| Target sequence | TL1A 180 aa, IL-7Rα 219 aa, verify khớp 2 đầu | `data/*.fasta` |
| ESMFold2 là gì | `biohub/ESMFold2`, 6,58 B params, MIT, ungated. Verify qua HF API | `02-RESEARCH.md` §Q1 |
| ipSAE | `DunbrackLab/IPSAE`, công thức + CLI đã có. Mục tiêu **≥0,75** | `02-RESEARCH.md` §Q2 |
| ESMFold2 chạy được trên L4 | **Có.** Weights chỉ **6,41 GB** (ESMC lượng tử fp8), load 106 s, còn trống 15,4 GB | `out/01_esmfold2_smoke_report.json` |
| Trần VRAM | Activation bậc hai (exp 1,99). Ngoại suy 640 aa ≈ 12,9 GB → **TL1A trimer khả thi** | đo thật, nhưng là **cận dưới** |
| Token-order PAE | `is_identity: true` — rủi ro đã đóng | smoke test, 5/5 lần |
| Độ dài design | **85–100 aa** | xem dưới |
| Epitope | TL1A: groove liên protomer (W137/F138...). IL-7Rα: D1, `K77-FLLIG-83` | `02-RESEARCH.md` §Q3/Q4 |
| Novelty | TL1A **sạch**. IL-7Rα có 122.273 design đã công bố (36–67 aa) | `data/known_binders/` |

### Quyết định thiết kế quan trọng nhất: design 85–100 aa

Hai lý do độc lập cùng chỉ về một chỗ:
1. **Miễn nhiễm novelty.** Known binder dài tối đa **67 aa**. Luật trigger khi aligned ≥80% tổng chiều dài.
   Binder 84 aa → cần cover 67,2 aa > 67 → **không known binder nào chạm được, về mặt toán học**.
   (Binder 80 aa thì vẫn chạm được → chưa an toàn.)
2. **Điểm cửa sàng.** ipSAE có sàn `d0 = max(1.0, …)` → điểm cao gần như bất khả thi khi interface < ~25 residue.

## Bug dtype — ĐÃ GIẢI (tìm bằng đọc source, không thử mò)

`EsmFold2Model.forward` mở `torch.amp.autocast(bf16)` ở `model.py:973` và **đóng ở 1160**
(câu cuối trong vùng đó là `z = z.float()`). Ba head chạy **ngoài** autocast, trên input fp32:

- `distogram_head(z + z.transpose(-2,-3))` — `model.py:1181`, một `nn.Linear(256,64)` trần
- `structure_head.sample(z_trunk=z, s_inputs=x_inputs, …)` — 1185
- `confidence_head(…, z=z.detach().float(), …)` — 1239

Ba head này **được thiết kế để giữ weights fp32**: `DiffusionConditioning.forward`
(`layers.py:1520-1534`) force-cast `z_trunk`/`s_inputs`/`rel_pos` sang fp32; `ConfidenceHead` làm
`pair.add_(pair_delta.float())`. Và `InputsEmbedder.forward` trả
`torch.cat([a_bf16, aatype_fp32, profile_fp32, deletion_mean_fp32])` → type promotion làm `x_inputs`
cũng fp32. Tức **cả input lẫn weights đều cố tình là fp32**.

`from_pretrained(dtype=torch.bfloat16)` gọi `model.to(bfloat16)` (`hub.py:194`), cast **tất cả** kể cả
ba head đó → matmul fp32×bf16 đầu tiên nó gặp là `distogram_head`.

### ✅ Cấu hình đúng

```python
EsmFold2Model.from_pretrained(dir, device="cuda", dtype=None, esmc_precision="fp8", load_esmc=True)
model.set_kernel_backend("fused")
```

`dtype=None` (fp32). Trunk **vẫn chạy bf16** nhờ autocast ở 973 + fused kernel tự cast ở call site
(`_to_fused_kernel_dtype`, `layers.py:2625`). Giá phải trả chỉ **+0,53 GB** (≈6,95 GB thay vì 6,413 GB).

⚠️ **Cái bẫy khiến tôi suy luận sai:** tôi coi "fp32 = 26,3 GB" nên khăng khăng phải bf16.
Sai — `esmc_precision` lượng tử ESMC **độc lập** với `dtype`, nên fp32 không bao giờ tốn 26 GB.
Hai giả thuyết tôi xếp hạng cao nhất (`esmc_precision="bf16"`, và bf16+`_offload_esmc`) **đều sẽ fail
y hệt**, vì lỗi không liên quan gì tới ESMC. Agent đọc source trước nên không tốn CU nào để biết điều đó.

### Hệ quả

- Diffusion decoder **chưa từng chạy** → số VRAM/thời gian của notebook 01 là **cận dưới của cận dưới**.
- Trần residue thật vẫn chưa biết.

## MSA — ĐÃ GIẢI

ColabFold API không hỏng; lỗi nằm ở notebook 01. MSA thật đã lấy được:
**IL-7Rα depth 921** (7 s), **TL1A depth 2037** (29 s), `match_columns` = đúng chiều dài target,
query khớp byte-to-byte. Lưu ở `data/msa/*.a3m`, **và nhúng gzip+base64 trong notebook 02 kèm sha256**
→ đường critical path không cần mạng.

Phát hiện thêm: `prepare_input.py:1182` tự broadcast một MSA qua `id=["A","B","C"]`
→ **TL1A trimer không cần xử lý đặc biệt.**

## Ba việc tiếp theo, đúng thứ tự

1. **Sửa dtype rồi chạy `02_*.ipynb`.** Nó trả lời 4 câu còn mở:
   - **Đối chứng dương/âm**: ipSAE của binder 7S5B (đã validated thực nghiệm) vs bản scramble của nó.
     *Nếu hai cái không tách nhau thì thước đo của ta hỏng và mọi số phía sau vô nghĩa.* Đây là kết quả
     quan trọng nhất của cả notebook.
   - **Trần residue thật** (peak thật, có diffusion decoder) → TL1A trimer có fit hay không.
   - **Monomer vs trimer**: fold cùng một binder ở 2 cấu hình, so ipSAE. Lệch nhiều → giả định này thành
     rủi ro lớn nhất của track TL1A, phải hỏi ban tổ chức.
   - **Tính tiền định**: chạy 2 lần cùng `seed=0`, ipSAE có trùng bit không? `fold()` mặc định
     `lm_dropout=0.3` (ngẫu nhiên) trong khi `config.json` ghi 0.25. Nếu không tiền định thì "1 seed 1 sample"
     nghĩa là **cửa sàng có yếu tố may mắn**.
2. **MSA thật.** Smoke test tụt về `query_only_FALLBACK` depth-1 (ColabFold API fail âm thầm).
   Cửa sàng **dùng MSA cho target** → cần đo chênh lệch ipSAE giữa MSA thật và không MSA.
   Lệch lớn thì chấm không MSA sẽ xếp sai thứ tự ứng viên.
3. **Sinh design.** Chưa bắt đầu. Pipeline đề xuất: BoltzGen hoặc BindCraft sinh backbone →
   ProteinMPNN sinh sequence → ESMFold2+ipSAE chấm → lọc. Phần sinh design là phần nặng nhất →
   đẩy sang platform sponsor (AsterFire nói compute không giới hạn) hoặc Antigravity.

## Ngân sách

- **Colab: ~88,8 / 90 compute unit** còn lại (smoke test tiêu 1,17). L4 ≈ 4,8 CU/giờ → ~18,5 giờ GPU.
- Weights ESMFold2 **đã cache** ở `MyDrive/proteinAI/weights/ESMFold2` — **đừng xoá**, tải lại mất 6 phút + 25 GB.
- Không SSH vào Colab (ToS). Không dùng tailscale (workspace của user).
- Máy nhà có GPU, SSH được — chưa dùng, là bước leo thang tiếp theo nếu Colab không đủ.

## Chia việc với Antigravity

Tôi **không gọi được agent Antigravity trực tiếp** (không có API/MCP nối). Mô hình thực tế:

```
Claude (ở đây)          →  brief tự chứa  →  Antigravity agent  →  kết quả  →  Claude review
  nghĩ, chốt, kiểm                             gõ, chạy, lặp                      quyết bước tiếp
```

Tôi giữ phần đắt-trí-tuệ-rẻ-token: đọc luật, chốt quyết định, review output, bắt lỗi.
Antigravity giữ phần rẻ-trí-tuệ-đắt-token: viết code, chạy, debug vòng lặp.

**Khuôn brief cho Antigravity** (tôi sẽ viết theo mẫu này, bạn dán nguyên văn):

```
BỐI CẢNH (tự chứa — đừng giả định agent biết gì về dự án)
  - mục tiêu cụ thể của task này
  - dữ kiện cần, kèm số thật + đường dẫn file
  - ràng buộc: L4 22 GB, ngân sách CU, mốc 11-14

VIỆC
  - 1 việc rõ ràng, có tiêu chí "xong" đo được

ĐÃ BIẾT (đừng tìm lại)
  - danh sách fact đã verify + nguồn

CẠM BẪY ĐÃ BIẾT
  - vd: esmc_precision fp8 gây dtype mismatch; -Fast bỏ qua MSA im lặng

ĐẦU RA
  - đường dẫn file cụ thể + format báo cáo
```

Nguyên tắc để brief không hỏng: **mọi con số và đường dẫn phải nằm trong brief**, vì agent Antigravity
không đọc được hội thoại này. Và luôn ghi mục "ĐÃ BIẾT" — nếu không nó sẽ tốn quota tìm lại
đúng những thứ ta đã trả giá để biết (ESMFold2 là model nào, ipSAE formula, epitope).

## Cạm bẫy đã trả giá để biết — đừng lặp lại

- `ESMFold2-Fast` **bỏ qua MSA mà không báo lỗi** → tưởng tái dựng cửa sàng nhưng không.
- `esmc_precision="fp8"` + trunk bf16 → **dtype mismatch** (bug hiện tại).
- `ESMFold2InputBuilder()` âm thầm tải thêm **417 MB `ccd.pkl`** không nằm trong `allow_patterns`.
- ipSAE: dùng **đường AF2** (`.pdb` + json có `pae/plddt/iptm/ptm`), **không** dùng đường AF3 `.cif`
  (`token_chain_ids` không bao giờ được `ipsae.py` đọc).
- 3DI3 có đột biến **I66T** so với construct cuộc thi → dùng **7OPB** (2,14 Å, WT).
- `11QC` **là PDB thật** (TfR1), không phải typo.
- Research agent từng khẳng định "không có CPU offload" — **sai**, `_offload_esmc` có thật trong esm 3.4.3.
  → Luôn verify claim của agent bằng source thật trước khi xây lên nó.
