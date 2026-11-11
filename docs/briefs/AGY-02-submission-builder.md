# AGY-02 — Submission builder + validator

> Brief cho agent Antigravity. **Tự chứa.** Làm xong ghi báo cáo vào `docs/briefs/AGY-02-RESULT.md`.
> Làm **sau** AGY-01 (nó dùng kết quả novelty check), nhưng code được song song.

## BỐI CẢNH

Cuối cuộc thi ta phải upload **một CSV cho mỗi target** lên website ban tổ chức.
Nếu CSV sai format, upload bị từ chối — và **deadline có thể là 14/11**, không có cơ hội sửa lại.
Nên cần một công cụ build CSV và **tự validate y hệt validator của website** trước khi upload.

May là ta đã trích được logic validator thật từ JS bundle của website.
Bản tham chiếu nằm ở **`data/validate_submission_csv.py`** — đọc nó trước khi viết gì.

## VIỆC

Viết `tools/build_submission.py`: nhận danh sách sequence ứng viên + metadata, xuất ra CSV đúng chuẩn,
rồi tự validate và **từ chối xuất file nếu không pass**.

## ĐÃ BIẾT — ĐỪNG TÌM LẠI

Luật format, đã trích từ bundle của site, **chính xác, đừng đi tìm lại**:

- Header **phải chứa** `name` và `sequence`. **Không phân biệt hoa/thường**, **không quan trọng thứ tự**,
  **cột lạ bị bỏ qua** (không gây lỗi).
- Tối đa **20 sequence** → tức tối đa **21 dòng** kể cả header.
- `name`: dài **1–64 ký tự** (inclusive cả hai đầu).
- `sequence`: bị **strip whitespace** và **uppercase** trước khi kiểm.
- Chỉ cho phép **20 amino acid chuẩn**: `ACDEFGHIKLMNPQRSTVWY`. Ký tự khác (`X`, `*`, `-`, số…) → lỗi.
- Độ dài sequence cho **track minibinder**: **40–100 aa**, bounds **inclusive** (40 và 100 đều hợp lệ).
  (Track VHH là 100–150 aa, nhưng ta **không** thi track VHH — chỉ minibinder.)
- Encoding **UTF-8**, **BOM phải được strip** (`﻿`).
- File ≤ **10 MB**.
- Mỗi target **một file riêng**.

Ta thi **2 track, cả hai đều minibinder**:

| target key | tên | file target sequence |
|---|---|---|
| `tl1a` | TL1A / TNFSF15 | `data/TL1A_O95150_L72-L251.fasta` (180 aa) |
| `il7ra` | IL-7Rα / CD127 | `data/IL7Ra_P16871_E21-D239.fasta` (219 aa) |

Dữ kiện khác đã chốt:

- Sequence ta submit sẽ dài **85–100 aa** (trong khoảng hợp lệ 40–100).
- Chỉ submit **mature protein**: không signal peptide, không tag, không Fc, không fusion.
  → Công cụ nên **cảnh báo** nếu phát hiện dấu hiệu tag còn sót, ví dụ chuỗi `HHHHHH`, `GGGGS` ở đầu/cuối,
  hoặc các motif tag phổ biến (`DYKDDDDK` = FLAG, `EQKLISEEDL` = myc). Cảnh báo, không tự xoá —
  **không bao giờ tự sửa sequence của người dùng**.

## YÊU CẦU CHỨC NĂNG

1. **Build**: `python3 tools/build_submission.py --target il7ra --input <candidates> --out out/il7ra_submission.csv`
2. **Validate trước khi ghi.** Không pass → in lỗi theo từng dòng (số dòng + tên + lý do), **không ghi file**.
3. **Đối chiếu chéo** với `data/validate_submission_csv.py`: sau khi ghi, chạy luôn validator tham chiếu đó
   lên file vừa tạo và khẳng định nó cũng pass. Hai bên không đồng ý → **coi là lỗi**, in ra chỗ lệch.
   *(Đây là điểm cốt lõi của task: validator của bạn phải khớp bản tham chiếu, không phải chỉ "trông có vẻ đúng".)*
4. **Kiểm trùng nội bộ**: nếu hai sequence trong cùng file giống nhau 100%, cảnh báo —
   luật thi coi coverage >90% + identity 100% là trùng, và trùng thì chỉ bản submit đầu được giữ.
5. **Nối với novelty check**: có cờ `--novelty-report <csv>` nhận output của `tools/novelty_check.py` (AGY-01).
   Nếu sequence nào `verdict == FAIL`, **từ chối build** kèm thông báo rõ, trừ khi có `--force`.
6. In **bảng tóm tắt** trước khi ghi: số sequence, min/max/median độ dài, có cảnh báo gì.

## SELF-TEST BẮT BUỘC

`tools/test_build_submission.py`, in PASS/FAIL từng case. Không báo xong nếu còn case đỏ.

| # | Input | Kết quả bắt buộc |
|---|---|---|
| 1 | 20 sequence hợp lệ, dài 85–100 aa | build thành công, validator tham chiếu cũng pass |
| 2 | **21** sequence | từ chối, báo vượt giới hạn 20 |
| 3 | Một sequence dài **39 aa**, và một dài **101 aa** | cả hai bị từ chối (ngoài khoảng 40–100) |
| 4 | Một sequence dài đúng **40 aa**, một đúng **100 aa** | **cả hai PASS** (bounds inclusive) |
| 5 | Sequence chứa `X` | từ chối, báo ký tự không hợp lệ |
| 6 | `name` dài 65 ký tự; `name` rỗng | cả hai bị từ chối. `name` dài đúng 64 → pass |
| 7 | Header viết `NAME,SEQUENCE`; header thứ tự `sequence,name`; header có thêm cột `note` | **cả ba đều pass** |
| 8 | File có BOM ở đầu | pass (BOM được strip) |
| 9 | Sequence có khoảng trắng giữa (`"ABC DEF"`) và chữ thường | pass sau khi strip + uppercase |
| 10 | Hai sequence trùng nhau 100% | build được nhưng **có cảnh báo trùng** |
| 11 | Sequence kết thúc bằng `HHHHHH` | build được nhưng **có cảnh báo nghi tag** |
| 12 | `--novelty-report` có một dòng `FAIL` | từ chối build; thêm `--force` thì build được |

Case 4 và case 7 là hai case dễ làm sai nhất — bounds inclusive, và header linh hoạt.

## ĐẦU RA

```
tools/build_submission.py
tools/test_build_submission.py
docs/briefs/AGY-02-RESULT.md
```

Báo cáo phải có: output self-test dán nguyên (12/12), những chỗ validator của bạn từng lệch với bản
tham chiếu và đã sửa thế nào, các quyết định bạn tự đưa ra vì brief chưa nói rõ, hạn chế còn lại.

## RÀNG BUỘC

- **KHÔNG tự sửa sequence của người dùng.** Chỉ cảnh báo. Việc im lặng "làm sạch" một sequence rồi submit
  là kiểu lỗi tệ nhất có thể xảy ra ở đây.
- **KHÔNG sửa/xoá gì trong `data/`.** Chỉ đọc. `data/validate_submission_csv.py` là bản tham chiếu, giữ nguyên.
- **KHÔNG bịa sequence protein.** Self-test cần sequence thì sinh random có seed cố định và ghi rõ là random.
- Không cần mạng, không cần GPU.
- Code tiếng Anh, báo cáo tiếng Việt.
