# AGY-01 — Novelty checker

> Brief cho agent Antigravity. **Tự chứa** — làm đúng file này, không cần context khác.
> Làm xong ghi báo cáo vào `docs/briefs/AGY-01-RESULT.md`.

## BỐI CẢNH

Dự án: thi thiết kế protein (长江杯 2026). Ta sẽ submit tối đa **20 sequence minibinder** cho mỗi target.
Ban tổ chức có **luật novelty** — vi phạm là **bị loại sequence đó**:

> So sánh với các minibinder đã biết của cùng target: khi vùng aligned **cover ≥80% tổng chiều dài**
> thì **sequence identity phải ≤70%**.
>
> Nguyên văn: 与每个靶点的已知minibinder进行比对，比对区域覆盖总长度的至少 80%时，序列一致性须 ≤70%

Rủi ro thật: với target IL-7Rα có **122.273 minibinder đã công bố** (Cao et al., Nature 2022).
Pipeline thiết kế của ta có thể tình cờ ra sequence quá giống → mất slot.

Cần một công cụ chạy **local trên laptop** (không dùng GPU, không dùng credit cloud) để tự kiểm **trước khi** submit.

## VIỆC

Viết `tools/novelty_check.py`.

**Input:** một file CSV/FASTA chứa các sequence ứng viên của ta, và một target (`tl1a` hoặc `il7ra`).
**Output:** với **mỗi** sequence ứng viên, báo:

- `max_identity_pct` — identity cao nhất tìm được trong số các known binder **mà thoả điều kiện coverage ≥80%**
- `worst_hit_id` — tên known binder gây ra con số đó
- `coverage_pct` của hit đó
- `verdict` — `PASS` nếu không vi phạm, `FAIL` nếu có hit nào vừa cover ≥80% vừa identity >70%
- `n_hits_over_80pct_coverage` — số known binder đạt ngưỡng coverage

Thêm một bảng tổng hợp cuối: bao nhiêu PASS / bao nhiêu FAIL.

### Cách tính, định nghĩa chính xác — đọc kỹ

- Dùng **local alignment Smith-Waterman** (không phải global), vì luật nói "vùng aligned".
- `coverage_pct` = (số residue của **sequence ứng viên của ta** nằm trong vùng aligned) / (tổng chiều dài sequence
  ứng viên của ta) × 100.
  → **Mẫu số là chiều dài sequence CỦA TA**, không phải của known binder. Đây là chỗ dễ hiểu sai nhất.
- `identity_pct` = (số vị trí khớp đúng trong vùng aligned) / (chiều dài vùng aligned) × 100.
- Vi phạm khi: `coverage_pct >= 80` **VÀ** `identity_pct > 70`. Cả hai cùng lúc.
- Vì luật mơ hồ ở chỗ mẫu số coverage, hãy **tính và báo cả hai cách**: cột
  `coverage_ours_pct` và `coverage_theirs_pct` (mẫu số là chiều dài known binder).
  Verdict chính dùng `coverage_ours_pct`, nhưng thêm cột `verdict_strict` dùng
  `max(coverage_ours_pct, coverage_theirs_pct)` để ta biết trường hợp xấu nhất.

### Hiệu năng

122.273 × 20 = ~2,4 triệu alignment. Pure Python sẽ mất hàng giờ → **không chấp nhận**.
Dùng một trong hai, thứ tự ưu tiên:

1. `pip install parasail` — SIMD Smith-Waterman, cho ra score + traceback. Nhanh nhất, đơn giản nhất.
2. `mmseqs2` (AUR: `yay -S mmseqs2`) với `easy-search`, nếu parasail không lắp được.

Mục tiêu: **chạy xong 20 sequence dưới 5 phút**. Báo thời gian thật trong kết quả.
Nếu dùng mmseqs, chú ý nó là tool heuristic — phải đặt sensitivity cao (`-s 7.5`) và `--max-seqs` lớn,
nếu không nó bỏ sót hit. Và ghi rõ trong báo cáo là đã dùng heuristic.

## ĐÃ BIẾT — ĐỪNG TÌM LẠI

Dữ liệu **đã có sẵn trên đĩa**, đừng tải lại, đừng đi tìm nguồn khác:

```
data/known_binders/
  IL7Ra_denovo_minibinder_designs_IPD_Cao2022.fasta   122.273 design, dài 36-67 aa, median 63  ← DB chính cho il7ra
  IL7Ra_denovo_minibinders_VALIDATED.fasta            2 design, cả hai dài 55 aa
  IL7Ra_therapeutic_antibodies.fasta                  6 seq, 107-122 aa
  IL7Ra_PDB_antibody_Fabs.fasta                       4 seq, 213-225 aa
  TL1A_therapeutic_antibodies.fasta                   6 seq, 106-124 aa
  TL1A_natural_binder_DcR3.fasta                      2 seq, 176 aa và 300 aa
```

Dữ kiện đã xác nhận, dùng luôn:

- Known **de novo minibinder** dài nhất của IL-7Rα là **67 aa**. Không có cái nào dài hơn.
- **TL1A không có de novo minibinder nào được công bố** — chỉ có kháng thể và receptor tự nhiên DcR3.
  Nên với target `tl1a`, DB chỉ gồm kháng thể + DcR3, và kỳ vọng mọi ứng viên đều PASS dễ dàng.
- Sequence ứng viên của ta sẽ dài **85-100 aa** (đã chốt). Suy ra: binder 84 aa trở lên thì
  `0.8 × 84 = 67.2 > 67` → về mặt toán học **không known binder nào chạm được ngưỡng coverage**.
  **Công cụ phải tái tạo đúng kết luận này** — xem self-test case 4.

## SELF-TEST BẮT BUỘC

Viết `tools/test_novelty_check.py`, chạy được bằng `python3 tools/test_novelty_check.py`, in PASS/FAIL từng case.
**Không báo cáo xong việc nếu còn case nào đỏ.**

| # | Input | Kết quả bắt buộc |
|---|---|---|
| 1 | Lấy nguyên một sequence từ `IL7Ra_denovo_minibinders_VALIDATED.fasta` làm ứng viên, target `il7ra` | `max_identity_pct == 100.0`, `coverage_ours_pct == 100.0`, `verdict == FAIL` |
| 2 | Một sequence random 90 aa (seed cố định, in ra seed) , target `il7ra` | `verdict == PASS`, và `max_identity_pct < 70` |
| 3 | Lấy một known binder 55 aa rồi **đột biến 40% số vị trí** (seed cố định) , target `il7ra` | identity rơi vào khoảng 55-70% → `verdict == PASS`. In số thật ra để kiểm bằng mắt |
| 4 | Lấy một known binder 60 aa, **nối thêm 30 aa random** vào cuối (tổng 90 aa) , target `il7ra` | `coverage_ours_pct ≈ 66.7` (60/90) → **dưới 80 → `verdict == PASS`** dù identity vùng đó rất cao. Đây là case chứng minh logic coverage đúng; nếu nó FAIL thì bạn đã lấy sai mẫu số |
| 5 | Một sequence bất kỳ với target `tl1a` | chạy không lỗi, `n_hits_over_80pct_coverage` báo ra được (kỳ vọng 0) |
| 6 | File input rỗng, và file input có ký tự không hợp lệ (vd `X`, `*`) | báo lỗi rõ ràng, không crash bằng traceback thô |

Case 4 là case quan trọng nhất. Nếu nó sai, cả công cụ vô dụng.

## ĐẦU RA

```
tools/novelty_check.py          — CLI: python3 tools/novelty_check.py --input <file> --target il7ra|tl1a [--out <csv>]
tools/test_novelty_check.py     — self-test, in PASS/FAIL từng case
docs/briefs/AGY-01-RESULT.md    — báo cáo
```

Báo cáo `AGY-01-RESULT.md` phải có:

1. Dùng engine alignment nào, phiên bản bao nhiêu, vì sao chọn
2. **Kết quả self-test, dán nguyên output**, 6/6 case
3. Thời gian chạy thật cho 20 sequence × 122.273 DB
4. Cách cài dependency (lệnh chính xác, để lần sau dựng lại được)
5. Những chỗ bạn phải tự quyết định vì brief chưa nói rõ — ghi rõ đã quyết thế nào
6. Hạn chế đã biết của công cụ

## RÀNG BUỘC

- **TUYỆT ĐỐI KHÔNG bịa sequence protein.** Không chắc thì ghi `UNCONFIRMED`.
- **KHÔNG sửa, không xoá, không ghi đè** bất cứ gì trong `data/`. Chỉ đọc.
  Đặc biệt `data/known_binders/` và `data/msa/` là dữ liệu đã trả giá để thu thập.
- Chỉ chạy CPU trên laptop. **Không** dùng GPU, không gọi API cloud, không cần mạng
  (trừ lúc `pip install`).
- Máy: Arch Linux, python 3.14, 8 core, RAM 15 GB nhưng **chỉ còn ~3 GB khả dụng** →
  đừng nạp cả 122k sequence vào RAM dạng object nặng; stream hoặc xử lý theo batch.
- Code tiếng Anh, báo cáo tiếng Việt.
