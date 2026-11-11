
---

# ✅ KIỂM CHỨNG (do Claude chạy, 2026-11-11)

## Self-test: 12/12 PASS, ngay vòng đầu

```
.venv/bin/python3 tools/test_build_submission.py
→ Results: 12/12 PASS   0/12 FAIL — All cases passed.
```

Hai case dễ sai nhất đều đúng: **Case 4** (bounds inclusive — 40 aa và 100 aa đều hợp lệ) và
**Case 7** (header linh hoạt — `NAME,SEQUENCE` hoa, đảo thứ tự, thêm cột lạ, cả ba pass).

## Fuzz đối chiếu độc lập: 0/22 lệch

Không tin test của agy, tôi tự dựng fuzz so **trực tiếp** validator của agy với bản tham chiếu
`data/validate_submission_csv.py` (trích từ JS bundle của website) trên 22 case biên:

```
REF vs AGY._cross_validate lệch: 0/22
```

Case đã thử: header hoa/đảo/thêm cột/thiếu, biên 39–40–100–101 aa, ký tự `X`, chữ thường,
khoảng trắng giữa sequence, BOM, name rỗng/64/65 ký tự, 20 vs 21 sequence, dòng trống giữa file,
CRLF, file rỗng, dấu `"` trong CSV.

**Kết luận: validator khớp bản tham chiếu trên toàn bộ case biên đã thử.**

## Ghi chú về quá trình kiểm

Lần probe đầu của tôi báo lệch 21/21 — nhưng đó là **lỗi harness của tôi**, không phải của agy:
tôi đoán sai tên hàm (`validate_csv_text`) nên nhận `AttributeError` ở mọi case. Tên thật là
`_validate_our` / `_validate_reference` / `_cross_validate`. Dấu hiệu nhận biết: lệch *100%* ở mọi
case, kể cả case đơn giản nhất, thường có nghĩa là harness sai chứ không phải code sai.

## Cách dùng

```bash
.venv/bin/python3 tools/build_submission.py --target il7ra|tl1a \
    --input <candidates> --out out/<target>_submission.csv \
    [--novelty-report <novelty_check output csv>] [--force]
```

Chặn build nếu có sequence `verdict == FAIL` trong novelty report (trừ khi `--force`).
Cảnh báo (không tự sửa) nếu nghi còn tag: `HHHHHH`, `GGGGS`, `DYKDDDDK`, `EQKLISEEDL`.
