# 长江杯 2026 — Fact base (nguồn: official site + group chat)

> Nguồn: JS bundle của `https://aiprotein.sjtu.edu.cn` (tải 2026-11-11, file `assets/index-BZIjDIfc.js`)
> + `instruction.md` (log group chat WeChat, tới 2026-11-11 08:42).
> Raw data target: `data/targets_raw_from_bundle.txt`.
> **Mọi con số trong file này lấy từ nguồn chính thức — không suy diễn.** Chỗ nào suy diễn sẽ ghi rõ `[SUY LUẬN]`.

## 1. Mốc thời gian

| Mốc | Ngày |
|---|---|
| Công bố đề (7 target) | 2026-11-01 |
| Hạn submit theo rules doc | 2026-11-30 23:59:59 |
| 🔴 **MỐC LÀM VIỆC THỰC TẾ CỦA TA** | **chốt sequence 11-13, confirm 11-14** |
| Công bố shortlist vào wet-lab | 2026-12-01 → giữa 12/2026 |
| Wet-lab validation | 12/2026 |
| Công bố kết quả | **2/2027** |

Hôm nay **2026-11-11**.

🔴 **Tại sao làm theo mốc 11-14 chứ không phải 11-30** (verify agent đã kiểm nguồn):
- `/api/pages/rules` (sửa 2026-11-06) và trang schedule đều ghi **2026-11-30 23:59:59**.
- NHƯNG `/api/public/meta` — chính config điều khiển `submission_open`, tức nút submit — **vẫn ghi hết hạn
  2026-11-15**, và backend có sẵn lỗi `已过提交截止时间`.
- Bundle đã redeploy **2026-11-10** mà **chưa sửa** bug này, dù đã bị báo từ 2026-11-08.
- **Draft đã lưu có thể KHÔNG được tính**: luật nói bản *confirmed* cuối cùng mới tính.
→ Giả định xấu nhất: nút submit chết sau 11-15. Không đặt cược vào việc họ kịp sửa.

⚠️ Sửa lỗi của bản fact base đầu: `11QC` **là PDB thật** (TfR1 ectodomain, cryo-EM 2.4 Å, công bố
2026-05-27), không phải typo như tôi đoán ban đầu.

## 2. Bảy track (mỗi target = 1 track, xếp hạng/giải độc lập)

| # | Target | Loại phân tử | UniProt | Range | Ref structure | Oligomer | Tag | Yêu cầu thêm |
|---|---|---|---|---|---|---|---|---|
| 1 | EGFR | **VHH** | P00533-1 | Met1–Ser645 | 6ARU chain A | — | C-term His | **đo affinity ở pH 6.0** |
| 2 | TNF-α | **VHH** | P01375 | Val77–Leu233 | 1TNF chain A,B,C | Trimer | C-term His | — |
| 3 | GLP-1R | **VHH** | P43220 | Arg24–Ser432 | 6X18 chain R | — | N-term Flag + C-term His | **Agonist** |
| 4 | TfR1 (CD71) | **VHH** | P02786 | Cys89–Phe760 | 11QC chain A,B | Dimer | N-term His | — |
| 5 | **TL1A** (TNFSF15) | **minibinder** | O95150 | Leu72–Leu251 | 2RE9 chain A,B,C | Trimer | N-term His | — |
| 6 | TNFR1 (TNFRSF1A) | **VHH** | P19438-1 | Leu30–Thr211 | 9S3I chain A | — | C-term His | — |
| 7 | **IL-7Rα** (CD127) | **minibinder** | P16871-1 | Glu21–Asp239 | 3DI3 chain B | — | C-term His | — |

Ghi chú chính thức: *"靶点相关信息以序列（UniProt编号）和序列范围为准，参考结构仅作展示"* → **sequence + range là chuẩn, structure chỉ minh hoạ**.
`11QC` có vẻ là typo của PDB (PDB id 4 ký tự: khả năng là `1IQC`? cần verify) — ghi nhận để kiểm tra.

→ **5 track VHH, 2 track minibinder.**

## 3. Yêu cầu thiết kế

**Minibinder track** (TL1A, IL-7Rα):
- Dài **40–100 aa** (bounds **inclusive**, 40 và 100 đều hợp lệ)
- Dạng: 1 chuỗi polypeptide, 1 module binding
- Novelty: so với known minibinder của cùng target, khi aligned region ≥80% full length → **sequence identity ≤70%**

🟢 **QUYẾT ĐỊNH THIẾT KẾ: design 85–100 aa.** Hai lý do độc lập chỉ về cùng một chỗ:

1. **Miễn nhiễm luật novelty.** Known de novo IL-7Rα minibinder (122.273 design, Cao et al. 2022 IPD archive,
   file `data/known_binders/IL7Ra_denovo_minibinder_designs_IPD_Cao2022.fasta`) có độ dài **36–67 aa, median 63,
   không cái nào >67 aa** (tôi đã đếm trực tiếp). Luật trigger khi aligned ≥80% tổng chiều dài:
   - binder 80 aa → cần cover 64 aa → known binder 67 aa **vẫn chạm được** ⚠️
   - binder **84 aa** → cần cover 67,2 aa > 67 aa → **không một known binder nào chạm được về mặt toán học** ✅

   → ≥84 aa thì an toàn *do cấu trúc bài toán*, không cần may mắn. Lấy 85 aa làm sàn cho chắc.
   (Verify agent đo: pipeline thường vượt 70% identity **9,2%** số lần ≈ 2/20 sequence bị loại nếu design ngắn.)

2. **Tối ưu điểm cửa sàng.** ipSAE có sàn `d0 = max(1.0, ...)` khiến điểm cao gần như không đạt được khi
   interface < ~25 residue → binder dài, footprint rộng cho điểm cao hơn. Chi tiết `02-RESEARCH.md` §Q6.

**Độ cạnh tranh cửa sàng** (verify agent ước lượng): ~800 slot public vs ~2.900 sequence public → **chật 3–4 lần**,
tỉ lệ nhận ~28%, cần khoảng **p72** (bi quan: p88). → **Mục tiêu ipSAE ≥0.75.**
AsterFire quảng cáo 0,63 là "đáng tin" ⇒ 0,63 chỉ là trung vị của đám đông, không phải mốc an toàn.

**VHH track** (EGFR, TNF-α, GLP-1R, TfR1, TNFR1):
- Dài **100–150 aa**
- Dạng: 1 VHH variable domain, **đủ framework + CDR** (không được submit riêng CDR)
- Novelty: so với known VHH cùng target, khi aligned region ≥80% vùng CDR → **CDR identity ≤70%**
- Xác nhận từ organizer (2026-11-10): **được phép giữ FR của VHH đã biết và design de novo CDR**. Novelty chỉ xét **CDR**, không so FR. Nhưng full-length alignment vẫn phải thoả.

**Chung:**
- Chỉ 20 amino acid chuẩn
- Chỉ submit **mature protein**: không signal peptide, không tag, không Fc, không fusion
- Phải là **de novo**; cấm submit sequence đã biết dương tính với cùng target, hoặc lấy nó làm điểm khởi đầu để optimize affinity
- **Tối đa 20 sequence / track**. Được tham gia nhiều track.

## 4. Luật vào vòng wet-lab (QUAN TRỌNG NHẤT với team public)

- Tổng ~**7000** sequence được wet-lab miễn phí, **1000 slot / track**.
- **Seed team** (20 team, do organizer mời): submit → compliance check → **vào hết**, không lọc chất lượng.
- **Public team** (chúng ta): score **từng sequence**, sort giảm dần, lấy top:
  - Structure prediction: **ESMFold2**
  - Score: **ipSAE**
  - Config: **1 seed, 1 sample, no template, MSA CHỈ cho target sequence, KHÔNG MSA cho protein được design**
- Top 50 affinity mỗi track còn được test thêm developability + cell activity (không ảnh hưởng ranking chính).

⚠️ Cảnh báo từ group (2026-11-05): **ESMFold2 ≠ ESMFold**, hai model khác nhau rất nhiều, cuộc thi dùng **ESMFold2**. → Phải xác định chính xác ESMFold2 là model nào. **Việc này đã giao cho research agent.**

## 5. Xếp hạng & giải thưởng

- Ranking **hoàn toàn** theo KD đo thực tế (nhỏ → lớn). Không có đánh giá chủ quan.
- Mỗi team trong 1 track chỉ lấy **1 sequence tốt nhất** để xếp hạng; không trùng giải trong cùng track.
- Mỗi track top 6 có giải: **nhất 60k RMB (1), nhì 25k RMB (2), ba 10k RMB (3)**.
- Bonus trong 6 giải: developability pass +5k, cell activity dương +5k (cộng dồn tối đa +10k).
- Tổng prize pool ~1.4M RMB.

## 6. Cách submit

1. Login `aiprotein.sjtu.edu.cn` bằng **email captain + số điện thoại captain** (password ban đầu = số điện thoại; phải đổi password trước khi tạo team / submit).
2. Sau khi đăng ký thành công phải **chờ 1 ngày làm việc** mới login được.
3. Vào trang submit, chọn track.
4. Upload:
   - **Material 1**: sequence CSV — header **phải có `name` và `sequence`**; tối đa 20 dòng; file ≤10MB; hỗ trợ CSV/XLSX/XLS; **mỗi target 1 file riêng**.
   - **Material 2**: design summary — text (≤4000 từ khi confirm; ≤8000 ở chỗ khác) hoặc **PDF ≤20MB, mỗi target 1 bản**. *Không bắt buộc công khai model weights / training set / source code.*
5. Save được ghi đè nhiều lần; **confirm rồi là không sửa được track đó nữa**. Re-upload phải ghi lý do (≤2000 từ).
6. Chỉ **captain** submit được.
7. Trùng sequence giữa các team: coverage >90% và identity =100% → duplicate; chỉ giữ bản submit đầu tiên (kháng thể chỉ so CDR).

## 7. Tool / platform sponsor (miễn phí trong thời gian thi)

Không giới hạn model/platform/technical route. Platform sponsor (xem `docs/01-SPONSOR-SIGNUP.md` để biết cách đăng ký):

| Platform | URL | Ưu đãi |
|---|---|---|
| MatwingsVenus® 晓鹜 | https://matvenus.com/zh | 100k điểm + 100 RMB token; có **ESMFold2 + ipSAE**, RFdiffusion, Rfantibody, BindCraft, BoltzGen, Boltz2, AntiBMPNN, Mol-1 |
| AsterFire™ 星使智算 / AsterFire Go | https://sidereus-ai.com/ | compute không giới hạn trong thời gian thi; Protenix, Mol-1, MD, agent-driven |
| SPONGE Mokda™ | https://sponge-mokda.cn/ | RFdiffusion + SPONGE MD miễn phí |
| 书生·端砚 (InternAI Discovery) | https://discovery.intern-ai.org.cn/chat | compute ("墨点值"), dataset, lit search |
| 元生·开物 OriBench | https://oribench.lglab.ac.cn/ | 300+ drug-discovery skills, token miễn phí |

## 8. Câu hỏi mở cần verify

- [ ] Deadline submit thật là 11-30 hay 11-15? (frontend báo lỗi) → verify agent đang chốt
- [x] **ESMFold2 = `biohub/ESMFold2`** — đã xác định và tôi đã verify độc lập qua HF API:
      6,579,755,843 params, MIT, ungated, 352k downloads, tạo 2026-05-13, lưu **F32 → 26.3 GB**.
      Frozen ESMC-6B + 48 folding layer + diffusion all-atom decoder. **Nhận MSA per-chain** (khớp đúng
      config cửa sàng), trả `result.pae/plddt/ptm/iptm`. **Không phải** ESMFold 2023, **không phải** ESM3.
      ⚠️ Checkpoint `-Fast` **bỏ qua MSA mà không báo lỗi** → cửa sàng phải đang dùng model full.
      Chi tiết: `docs/02-RESEARCH.md` §Q1.
- [x] **ipSAE** = `DunbrackLab/IPSAE`. `1/(1+(PAE/d0)²)`, `d0=max(1.0, 1.24(L−15)^⅓−1.8)`, max trên residue aligned.
      Ngưỡng hiệu chuẩn: **≥0.61** là mốc dự báo thành công thực nghiệm tốt nhất (meta-analysis 3766 binder);
      nhắm **≥0.70** cho minibinder. Chi tiết §Q2.
- [ ] `11QC` là PDB id gì? (không thuộc track ta chọn, ưu tiên thấp)
- [x] Known binder DB → `data/known_binders/`. **Bất đối xứng quan trọng**: TL1A **chưa có** de novo
      minibinder công bố; IL-7Rα **đã có** (PDB 7S5B, 55 aa) + bộ design Cao 2022.
- [x] Target sequence đã cắt + verify khớp hai đầu → `data/*.fasta`

## 9. Ràng buộc tài nguyên thật (đo trên máy user, 2026-11-11)

- GPU **L4, 21.5 GB cấp phát được**. torch 2.14.1+cu130, jax 0.11.1. 52 GB RAM, 12 core.
- Drive 106 GB trống, ghi 208 MB/s. Mạng tải weights OK (các lỗi 403/404 lúc probe là do test root path).
- **Compute: 90 unit, Colab Pro. L4 ≈ 4.8 CU/giờ → chỉ còn ~18.7 giờ GPU cho toàn bộ dự án.**
  Đây là ràng buộc chặt nhất của cả kế hoạch.
- Hệ quả VRAM: ESMFold2 bf16 ≈ 13 GB chỉ riêng weights; maintainer nói 24 GB kham **~500 residue**.
  → IL-7Rα + binder (259–319 aa) **fit**. TL1A **trimer** + binder (580–640 aa) **không fit**.
