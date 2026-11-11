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
| ESMFold2 chạy được trên L4 | **Có.** Weights chỉ **6,41 GB** (ESMC lượng tử fp8), load 108,9 s, VRAM weights 6,85 GB | `notebooks/output_01.json` (Notebook 02) |
| Bug Dtype | **Đã giải & verify.** `dtype=None, esmc_precision="fp8", kernel_backend="fused"` | `notebooks/output_01.json` |
| Thước đo ipSAE | **Hoạt động xuất sắc.** Dương 0,7728 (min 0,7488) vs Âm 0,0 (Gap 0,77, verdict: PASS vs REVIEW) | `03_screening_report.json` & Notebook 02 |
| Trần VRAM & Residues | **775 aa OK.** Trimer 716 aa peak 15,71 GB (200 steps) / 9,79 GB (50 steps) an toàn trên L4 22 GB | đo thật Notebook 02 & 03 |
| TL1A Monomer vs Trimer | Monomer ipSAE 0,267 vs Trimer 0,545 (delta 0,277) → **Bắt buộc fold dạng Trimer** | đo thật Notebook 02 |
| MSA Delta | Lệch rất nhỏ: -0,01326 (có MSA 0,7483 vs không MSA 0,7616). Đã nạp từ `data/msa/*.a3m` | đo thật Notebook 02 |
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

### ✅ Kết quả thực nghiệm Notebook 02 (`notebooks/output_01.json` — ĐÃ CHẠY THÀNH CÔNG)

Notebook 02 đã chạy hoàn tất trên NVIDIA L4 (22,03 GB VRAM, PyTorch 2.14.1+cu130, esm 3.4.3). Toàn bộ 4 câu hỏi còn mở và giả thuyết đo lường đã được trả lời bằng số liệu thực nghiệm chính xác:

1. **Bug Dtype & Decoder đã giải quyết hoàn toàn:**
   - Cấu hình chuẩn: `dtype=None, esmc_precision="fp8", kernel_backend="fused"`.
   - VRAM load weights: **6,85 GB** (load trong 108,9 s).
   - Diffusion decoder và các heads chạy bình thường, tính PAE và structure hoàn chỉnh không còn lỗi `Float and BFloat16`.

2. **Đối chứng dương/âm — Thước đo sàng lọc HOẠT ĐỘNG XUẤT SẮC:**
   - **Positive control** (7S5B binder 55 aa + IL-7Rα 219 aa):
     - `ipSAE_min = 0.748344` (asym A->B: 0.748344, B->A: 0.769589)
     - `ipTM = 0.93714`, `pLDDT = 0.87899`
   - **Negative control** (scrambled sequence của 7S5B + IL-7Rα 219 aa):
     - `ipSAE_min = 0.0` (asym A->B: 0.0, B->A: 0.0)
     - `ipTM = 0.2368`, `pLDDT = 0.75839`
   - **Gap = 0.74834** (vượt xa ngưỡng dương 0,5 và ngưỡng gap 0,15).
   - **Verdict: `THUOC_DO_HOAT_DONG`**. Thước đo hoàn toàn phân tách được binder thật picomolar và chuỗi ngẫu nhiên.

3. **Trần residue thật & Khả năng chạy TL1A Trimer trên L4:**
   - Đã test probe trần: **775 aa OK không OOM** (`first_oom_residues: null`).
   - TL1A trimer (3 × 180 aa = 540 aa) + binder 85–100 aa = **625–640 aa** hoàn toàn nằm trong vùng an toàn trên GPU L4.
   - Thử nghiệm thực tế TL1A trimer (540 aa) + DcR3 probe (176 aa) = **716 aa**: peak VRAM chỉ **9,785 GB**, thời gian chạy 102,67 s. Không hề chạm ngưỡng OOM của L4 (22 GB).

4. **TL1A Monomer vs Trimer — BẮT BUỘC mô hình hóa dạng Trimer:**
   - Monomer (TL1A 180 aa + DcR3 176 aa = 356 aa): `ipSAE_binder_best = 0.267342` (ipTM 0.57335, peak VRAM 7.98 GB, 23.37 s).
   - Trimer (TL1A 540 aa + DcR3 176 aa = 716 aa): `ipSAE_binder_best = 0.545034` (ipTM 0.80854, peak VRAM 9.785 GB, 102.67 s).
   - **Delta: +0.27769** (>0.277) — chênh lệch cực lớn!
   - Kết luận: Đúng với sinh học cấu trúc (TL1A tạo trimer chức năng và groove nằm ở giao diện giữa các protomer). **Bắt buộc phải fold đánh giá candidate dạng trimer (3x180 aa + binder)**. Chấm trên monomer sẽ làm mất oan các binder tốt.

5. **Tác động của MSA (MSA Delta):**
   - Chạy kèm MSA (depth 921): `ipSAE = 0.748344`, `ipTM = 0.93714`.
   - Chạy không MSA (single sequence mode): `ipSAE = 0.761601`, `ipTM = 0.94157`.
   - **Delta: -0.01326** (rất nhỏ, ~1,3%).
   - Target MSA đã có sẵn và nạp từ `data/msa/*.a3m` (IL-7Rα: `IL7Ra_P16871_E21-D239.uniref.a3m` depth 921; TL1A: `TL1A_O95150_L72-L251.uniref.a3m` depth 2037).

6. **Tính tiền định (Determinism):**
   - 2 lần chạy cùng seed: `ipSAE_abs_diff = 0.000703` (0.748344 vs 0.749047), `pae_max_abs_diff = 3.879`.
   - Do `lm_dropout=0.3` (ngẫu nhiên) trong fold default.

### ✅ Kết quả thực nghiệm Notebook 03 Screening Pipeline trên Colab L4 (`output_03_esmfold2_screening_colab/03_screening_report.json`)

Notebook 03 đã chạy trọn vẹn trên Colab GPU NVIDIA L4 (19,45 phút, tiêu tốn ~1,56 CU), thực thi chuẩn xác gate đánh giá 20 loops, 200 diffusion steps, seed=0, 1 sample, no-template:

1. **Positive control 7S5B trên IL-7Rα (phức 274 aa):**
   - `ipSAE_binder = 0.772791` (asym B->A: 0.772791, A->B: 0.74884), `ipSAE_min = 0.74884`.
   - `ipTM = 0.9374`, `pLDDT_mean = 87.86%`.
   - Thời gian fold: **35,6 s**, peak VRAM: **8,27 GB**.
   - Verdict: **`PASS`** (vượt xa ngưỡng cạnh tranh).
2. **Negative control scrambled trên IL-7Rα (phức 274 aa):**
   - `ipSAE_binder = 0.0`, `ipSAE_min = 0.0` (asym A->B: 0.0, B->A: 0.0).
   - `ipTM = 0.2071`, `pLDDT_mean = 78.22%`.
   - Thời gian fold: **22,2 s**, peak VRAM: **8,27 GB**.
   - Verdict: **`REVIEW`** (bị loại hoàn toàn).
3. **TL1A Trimer probe (DcR3 176 aa + TL1A 3×180 aa = 716 aa):**
   - `ipSAE_binder = 0.561002` (cặp Z->C: 0.561002, Z->B: 0.545559), `ipSAE_min = 0.343129`.
   - `ipTM = 0.8009`, `pLDDT_mean = 79.57%`.
   - Thời gian fold: **196,6 s** (~3,3 phút), peak VRAM: **15,71 GB**.
   - **Xác nhận trần VRAM**: 15,71 GB nằm an toàn dưới trần 22,03 GB của GPU L4; TL1A trimer + binder 85–100 aa (625–640 aa) sẽ chỉ tốn ~11–13 GB, chạy hoàn toàn ổn định không lo OOM.
   - Verdict: **`PASS`**.
4. **Xử lý Novelty check:**
   - Trong phiên chạy Colab, thư mục `data/known_binders/` chưa có trên Drive nên notebook báo `NO_DB` và fallback an toàn không crash.
   - Công cụ cục bộ `tools/novelty_check.py` trong repo đã có sẵn toàn bộ dữ liệu 122k design Cao 2022, 7S5B và kháng thể mẫu để thực hiện novelty check độc lập hoặc offline.

## MSA — ĐÃ GIẢI

ColabFold API không hỏng; lỗi nằm ở notebook 01. MSA thật đã lấy được:
**IL-7Rα depth 921** (7 s), **TL1A depth 2037** (29 s), `match_columns` = đúng chiều dài target,
query khớp byte-to-byte. Lưu ở `data/msa/*.a3m`, **và nhúng gzip+base64 trong notebook 02 kèm sha256**
→ đường critical path không cần mạng.

Phát hiện thêm: `prepare_input.py:1182` tự broadcast một MSA qua `id=["A","B","C"]`
→ **TL1A trimer không cần xử lý đặc biệt.**

## Các việc tiếp theo & Thứ tự ưu tiên hiện tại

1. **[X] ĐÃ XONG: Sửa dtype & chạy Notebook 02.** Thước đo xác nhận hoạt động (`THUOC_DO_HOAT_DONG`), trần residue 775 aa không OOM, TL1A trimer khả thi và bắt buộc dùng trimer.
2. **[X] ĐÃ XONG: MSA thật & đo delta.** Lệch chỉ -0.01326, MSA đã sẵn sàng trong `data/msa/*.a3m`.
3. **[X] ĐÃ XONG: Công cụ kiểm tra submission & novelty.** `tools/novelty_check.py` và `tools/build_submission.py` cùng test suite đã sẵn sàng, khóa chặt luật dài 85–100 aa.
4. **[X] ĐÃ XONG: Pipeline sàng lọc in-silico production (`notebooks/03_esmfold2_screening.ipynb`).**
   - Tái dựng chuẩn xác cửa sàng: 1 seed, 1 sample, no-template, nạp MSA target từ Drive cache (`IL7Ra` depth 921, `TL1A` depth 2037).
   - Mô hình hóa IL-7Rα dạng monomer (chain A 219 aa + binder chain B) và TL1A dạng trimer (chains A,B,C 540 aa broadcast MSA + binder chain Z).
   - Đánh giá ipSAE (AF2 path, cutoffs 10/10) + tích hợp Smith-Waterman parasail novelty filter (FAIL nếu cov ≥ 80% & ID > 70%).
   - Tự động lưu `03_screening_results.csv`, `03_screening_top_candidates.csv`, `03_screening_report.json` và khối `PASTE CHO CLAUDE`.
   - **Cơ chế Checkpoint Resume & Fast Reconnect**: `STAGE_TO_LOCAL = False` mặc định đọc mmap từ Drive (~100s, tránh copy 25GB mỗi lần đứt mạng); Cell 9 tự động scan checkpoint CSV/JSON cũ để bỏ qua các candidate đã fold xong khi bị ngắt kết nối.
5. **[X] ĐÃ XONG: Hạ tầng sinh candidates & Đã tạo trọn vẹn 100 candidates cho cả 2 track (`notebooks/candidate_pool_for_screening.csv`).**
   - **Track TL1A**: Đã chạy sinh thành công 50 design TL1A bằng Notebook 04 (`bindcraft_scaffold` + ColabDesign / ProteinMPNN):
     - 50/50 chuỗi hợp lệ tuyệt đối: độ dài 85–100 aa (min 85, max 100 aa), 20 amino acid chuẩn, không Cys tự do.
     - Phân bổ epitope: 37 design Nhóm A (liên diện rãnh trimer groove) + 13 design Nhóm B (bề mặt monomer).
     - Kiểm tra tính mới: 50/50 PASS (zero hit đối với cơ sở dữ liệu known binders của TL1A).
   - **Track IL-7Rα**: Đã sinh thành công 50 design IL-7Rα (`notebooks/candidate_pool_il7ra.csv`):
     - 50/50 chuỗi hợp lệ: độ dài 85–100 aa, 20 amino acid chuẩn, không Cys tự do, không motif khử amin (NG/DG/DP).
     - Phân bổ epitope: 38 design Nhóm A (rãnh chức năng D1 KKFLLIG 77–83) + 12 design Nhóm B (mặt bên/topologies thay thế).
     - Kiểm tra tính mới với `tools/novelty_check.py`: **50/50 PASS** (chạy 6.114.250 phép so khớp Smith-Waterman đối chiếu toàn bộ 122.285 known binders Cao 2022 + kháng thể; identity tối đa chỉ ~20–34% << 70%, kết quả ghi tại `notebooks/il7ra_novelty_report.json`).
   - **Manifest tổng hợp**: `notebooks/candidate_pool_for_screening.csv` đã hợp nhất đủ 100 design (50 TL1A + 50 IL-7Rα).
   - Prompt MatvenUS tiếng Trung chuẩn cho IL-7Rα (`docs/briefs/MATVENUS-02-il7ra-design.md`) với chiến lược miễn nhiễm novelty 85–100 aa và epitope D1 KKFLLIG (77–83).
   - **Kết quả MatvenUS S1 Pilot (IL-7Rα)**:
     - 40/40 RFdiffusion backbones nằm chuẩn xác trong khoảng [86, 100] aa.
     - 38/40 (95%) đạt geometric QC, khoảng cách tới glycan N49 đều ≥ 11,7 Å (vượt xa yêu cầu an toàn 8 Å).
     - Đã chọn Top 12 backbones tối ưu theo số tiếp xúc và hotspot (A8: 30 contacts, B8: 29 contacts, A17: 28 contacts, A24: 26 contacts, tiếp xúc đủ 7/7 hotspot).
     - **Đang chạy MatvenUS Task (`call_00_i4tlmhhrt9zcefhbd7otkgqz`)**: sinh 72 sequences (12 backbones × 6 seqs, `omit_AAs: "C"`, không motif NG/DG/DP, pI tránh 6.5–7.5, ưu tiên GRAVY > -1.2) và tự động sàng lọc tiếp qua ESMFold2 multimer.
   - **Kết quả sàng lọc Top 20 IL-7Rα qua Notebook 03 (`topcandidate.csv`)**:
     - Đã hoàn tất sàng lọc toàn bộ pool IL-7Rα trên Colab L4 (20 loops, 200 diffusion steps).
     - Candidate dẫn đầu: `IL7RA_cand_0002` (89 aa, ipSAE = 0.1522, ipTM = 0.5543, pTM = 0.7065, pLDDT = 70.30%, thời gian fold 28.8s, peak VRAM 8.61 GB, PASS novelty).
     - Ranks 2–5: `IL7RA_cand_0041` (ipTM 0.4473), `IL7RA_cand_0001` (ipTM 0.2213), `IL7RA_cand_0045` (ipTM 0.3719), `IL7RA_cand_0047` (ipTM 0.2645).
     - 20/20 candidate đều đạt độ dài 85–98 aa và 100% PASS novelty filter. Peak VRAM đồng đều ở mức 8.57–8.71 GB, thời gian fold ~27s/candidate.
   - Script hỗ trợ `tools/candidate_generator.py` và test suite cú pháp/unit tests đều xanh.
6. **[ ] ƯU TIÊN CAO NHẤT HIỆN TẠI: Chạy sàng lọc qua Notebook 03 và chọn Top 20 mỗi track.**
   - **Track TL1A**: Đang triển khai sàng lọc 50 candidate tự hành trên Colab L4 qua `scripts/run_screening_colab.py` (ESMFold2 trimer 540 aa + candidate Z 85–100 aa, checkpoint resume trên Drive).
   - **Track IL-7Rα**: Đã có Top 20 sơ bộ từ Notebook 03 (`topcandidate.csv`), đang chờ batch 72 sequence từ MatvenUS S1 pilot để sàng lọc bổ sung nâng cao ipSAE.
   - Trích xuất top 20 candidates mỗi track có ipSAE cao nhất và Novelty PASS.
   - Dùng `tools/build_submission.py` để đóng gói 2 file CSV nộp bài trước hạn chót.

## Ngân sách

- **Colab: 81,94 compute unit** còn lại sau các đợt test. L4 ≈ 4,8 CU/giờ → ~17 giờ GPU. Đã xác thực CLI (`colab auth login`) cho phép chạy tự hành headless session.
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
