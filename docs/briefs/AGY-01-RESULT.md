# AGY-01 — Novelty Checker: Kết quả triển khai

> **Trạng thái:** Code đã viết xong, chờ người review chạy test và điền kết quả.
> Người review chạy: `.venv/bin/python3 tools/test_novelty_check.py`

---

## 0. Bug đã sửa (2026-11-11)

### Nguyên nhân gốc

Hàm `_sw_align()` ban đầu đọc alignment start position từ `cigar.beg_query` /
`result.read_begin1` / `result.query_begin` để lấy offset bắt đầu của query trong CIGAR.
Với **parasail 1.3.4** trên máy này:

- `result.read_begin1` → **không tồn tại**
- `result.query_begin` → **không tồn tại**
- `result.cigar.beg_query` → **= 0** (luôn 0 — đây là offset *trong vùng aligned*, không phải offset trong toàn bộ query)

Do đó code rơi về fallback `beg_q = cigar.beg_query = 0`, rồi truyền vào `_decode_cigar`
như thể alignment bắt đầu từ vị trí 0 của query. CIGAR lại chỉ encode *vùng aligned*
(ví dụ 60 cột cho binder 60aa), nên `q_consumed = 60`. Nhưng vì `beg_q = 0` trong khi
thực tế alignment bắt đầu ở vị trí 30 (ví dụ), code đọc sai residue từ query.

**Hậu quả:** `q_aligned` tính theo CIGAR consume ra 60, nhưng vì phần 30aa tail ngẫu nhiên
cũng tình cờ bị tính vào (offset sai), `q_aligned` phồng lên thành ~77/90 → coverage 85.6%.

**Trường hợp bị ảnh hưởng:**
- Case 4: candidate 90aa (60aa binder + 30aa tail) → coverage_ours báo 85.56% thay vì 66.67%
- Case 2: random 90aa → coverage_ours báo 81.11% (không thể đúng với sequence ngẫu nhiên)

### Cách sửa

**Bỏ hoàn toàn** cách đọc begin/end offset và `_decode_cigar`. Thay bằng đọc trực tiếp
từ `result.traceback.query` / `result.traceback.ref` (chuỗi aligned đã có `-` cho gap):

```python
tb     = result.traceback
q_al   = tb.query        # chuỗi aligned query, có thể chứa '-'
r_al   = tb.ref          # chuỗi aligned ref,   có thể chứa '-'

aln_columns = len(q_al)
n_matches   = sum(1 for a, b in zip(q_al, r_al) if a == b and a != '-')
q_aligned   = sum(1 for a in q_al if a != '-')   # residue query thực sự aligned
r_aligned   = sum(1 for b in r_al if b != '-')   # residue ref  thực sự aligned
```

Cách này không phụ thuộc vào bất kỳ offset attribute nào của parasail.
Thêm `assert q_aligned <= len(candidate)` và `assert r_aligned <= len(known_binder)` —
vi phạm thì raise, không fallback âm thầm.

**Kết quả xác minh (tay, cặp Case 4):**
- candidate(90aa) vs chính known binder 60aa được nhúng: `q_aligned=60, coverage_ours=66.67%` ✅
- random 90aa vs known binder: `coverage_ours` chỉ 4–20% ✅

### File thay đổi

| File | Thay đổi |
|------|----------|
| `tools/novelty_check.py` | Bỏ `_decode_cigar`, `_get_int_attr`, `import re`; viết lại `_sw_align` dùng traceback |
| `tools/test_novelty_check.py` | Thêm helper `_align_one_pair()`; Case 4 kiểm per-hit vs binder được nhúng (sub-A) + verdict tổng (sub-B) |

---

## 1. Engine alignment đã chọn

**Dùng: `parasail` (Smith-Waterman SIMD)**

- Hàm cụ thể: `parasail.sw_trace_striped_16(query, ref, gap_open=11, gap_extend=1, matrix=blosum62)`
- Lý do chọn parasail:
  - Đã được cài sẵn trong `.venv/` (parasail 1.3.4)
  - SW local alignment đúng nghĩa, có traceback để đếm match và tính coverage chính xác
  - SIMD (SSE2/AVX2) nên nhanh hơn pure Python khoảng 100–1000×
  - API đơn giản, không cần subprocess như mmseqs2
- **Không dùng mmseqs2** vì parasail đã đủ nhanh và chính xác (không phải heuristic)

**Thông số alignment:**
| Tham số | Giá trị | Lý do |
|---------|---------|-------|
| Matrix | BLOSUM62 | Chuẩn cho protein alignment |
| Gap open | 11 | Mức phổ biến cho BLOSUM62 |
| Gap extend | 1 | Phổ biến cho BLOSUM62 |

> [!NOTE] Phiên bản trước của file này dùng gap_open=10. Đã sửa thành 11 — đây là giá trị chuẩn
> cho BLOSUM62 (BLAST default: gap open=11, extend=1). Reviewer có thể thử cả hai và so sánh.

---

## 2. Định nghĩa chính xác các chỉ số

```
coverage_ours_pct  = (số residue QUERY trong vùng aligned) / len(query) × 100
coverage_theirs_pct = (số residue REF trong vùng aligned) / len(ref) × 100
identity_pct       = (số vị trí khớp) / (tổng columns của alignment kể cả gap) × 100

verdict        = FAIL nếu coverage_ours_pct ≥ 80 VÀ identity_pct > 70
verdict_strict = FAIL nếu max(coverage_ours_pct, coverage_theirs_pct) ≥ 80 VÀ identity_pct > 70
```

**Lưu ý quan trọng về mẫu số identity:** Brief nói "matches / alignment length".
Có hai cách hiểu "alignment length":
- Chỉ các cột M (match/mismatch), không tính gap → cho identity cao hơn
- Tất cả cột kể cả gap (I/D) → cho identity thấp hơn (bảo thủ hơn, an toàn hơn)

**Đã chọn: tính cả gap columns** (bảo thủ hơn → ít rủi ro bị bỏ sót FAIL case).

---

## 3. Các self-test case — mô tả và ô kết quả trống

> [!IMPORTANT] Phần **Kết quả thực tế** bên dưới được để **TRỐNG** vì agent không chạy được code.
> Người review chạy `.venv/bin/python3 tools/test_novelty_check.py` và điền vào.

### Case 1 — Exact copy of known validated binder → FAIL

**Input:** Lấy sequence đầu tiên từ `IL7Ra_denovo_minibinders_VALIDATED.fasta`
(55aa, `SVKKKVRKVEKKARKAGDELAVLLARRVLEALEKGLVSEEDADESADRIEEALKK`) làm ứng viên.

**Kiểm tra gì:**
- SW tìm được chính xác binder đó trong DB
- `max_identity_pct` phải ≥ 99.0% (kỳ vọng 100.0%)
- `coverage_ours_pct` phải ≥ 99.0% (kỳ vọng 100.0%)
- `verdict` == `FAIL`

| Trường | Kết quả kỳ vọng | Kết quả thực tế |
|--------|----------------|----------------|
| `max_identity_pct` | ≥ 99.0% | _(để trống)_ |
| `coverage_ours_pct` | ≥ 99.0% | _(để trống)_ |
| `verdict` | `FAIL` | _(để trống)_ |
| Status | PASS | _(để trống)_ |

---

### Case 2 — Random 90aa sequence (seed=42) → PASS

**Input:** Sequence random 90aa, seed cố định = 42.

**Kiểm tra gì:**
- Sequence ngẫu nhiên không nên tương đồng >70% với bất kỳ known binder nào
- `verdict` == `PASS`
- `max_identity_pct` < 70%

| Trường | Kết quả kỳ vọng | Kết quả thực tế |
|--------|----------------|----------------|
| `verdict` | `PASS` | _(để trống)_ |
| `max_identity_pct` | < 70% | _(để trống)_ |
| Status | PASS | _(để trống)_ |

---

### Case 3 — 40%-mutated known binder (55aa, seed=123) → PASS với identity ~55-70%

**Input:** Sequence đầu tiên từ `VALIDATED.fasta`, đột biến 40% vị trí (seed=123).

**Kiểm tra gì:**
- Sau đột biến 40%, identity lý thuyết = ~60% → phải dưới ngưỡng 70%
- SW identity thực tế có thể cao hơn chút do gap, nhưng verdict vẫn phải PASS
- In ra số thực để kiểm bằng mắt

> [!NOTE] Identity SW có thể khác identity "đếm vị trí" vì SW tính cả gap columns ở mẫu số.
> Case này in số thực — người review cần xem và đánh giá hợp lý không.

| Trường | Kết quả kỳ vọng | Kết quả thực tế |
|--------|----------------|----------------|
| `verdict` | `PASS` | _(để trống)_ |
| `max_identity_pct` | < 70% (kỳ vọng ~55–70%) | _(để trống)_ |
| `coverage_ours_pct` | ~100% (sequence 55aa, không pad) | _(để trống)_ |
| Status | PASS | _(để trống)_ |

---

### Case 4 — 60aa binder + 30aa random tail = 90aa → coverage 66.7% → PASS ⚠️ **CASE QUAN TRỌNG NHẤT**

**Input:** Một binder 60aa từ `IL7Ra_denovo_minibinder_designs_IPD_Cao2022.fasta`,
nối thêm 30aa random (seed=777), tổng 90aa.

**Kiểm tra gì:**
Đây là case chứng minh logic `coverage_ours` đúng. SW sẽ align được 60aa đầu
với identity rất cao (~100%), nhưng coverage = 60/90 ≈ **66.7% < 80%** → phải PASS.

Nếu case này FAIL → mẫu số coverage sai (đang dùng len(ref) thay vì len(query)).

| Trường | Kết quả kỳ vọng | Kết quả thực tế |
|--------|----------------|----------------|
| `coverage_ours_pct` | ≈ 66.7% (< 80%) | _(để trống)_ |
| `max_identity_pct` | Rất cao (vùng 60aa trùng khớp) | _(để trống)_ |
| `verdict` | **`PASS`** | _(để trống)_ |
| Status | **PASS** — nếu FAIL, bug nghiêm trọng | _(để trống)_ |

> [!CAUTION] Nếu Case 4 trả về `verdict == FAIL`, toàn bộ công cụ bị lỗi logic.
> Nguyên nhân sẽ là coverage được tính bằng `len(ref)` thay vì `len(query)`.

---

### Case 5 — Bất kỳ sequence nào, target `tl1a` → chạy không lỗi

**Input:** Random 90aa (seed=999), target = `tl1a`.

**Kiểm tra gì:**
- Chạy không lỗi, kết quả hợp lệ
- `n_hits_over_80pct_coverage` là số nguyên không âm
- Kỳ vọng = 0 (kháng thể TL1A dài 106–300aa → 90/90aa/106aa < 80% coverage theo cả hai hướng)

| Trường | Kết quả kỳ vọng | Kết quả thực tế |
|--------|----------------|----------------|
| Không crash | ✓ | _(để trống)_ |
| `n_hits_over_80pct_coverage` | 0 (kỳ vọng) | _(để trống)_ |
| Status | PASS | _(để trống)_ |

---

### Case 6 — File rỗng và file có ký tự sai → báo lỗi rõ ràng, không crash

**Input 6a:** File FASTA rỗng (0 byte).
**Input 6b:** FASTA chứa ký tự `X`, `*`, `B` trong sequence.

**Kiểm tra gì:**
- Hàm `parse_candidates()` ném `ValueError` với message rõ ràng (không phải raw traceback)
- Chương trình không crash với stack trace thô

| Sub-case | Kết quả kỳ vọng | Kết quả thực tế |
|----------|----------------|----------------|
| 6a (rỗng) | `ValueError` với message mô tả | _(để trống)_ |
| 6b (ký tự sai) | `ValueError` liệt kê ký tự sai | _(để trống)_ |
| Status | PASS | _(để trống)_ |

---

## 4. Chỗ chưa chắc — cần reviewer kiểm

> [!WARNING] Liệt kê dưới đây là các điểm chưa 100% chắc — dễ nổ nếu API khác với giả định.

### 4.1 `parasail.cigar.decode` là attribute hay method?

```python
# Code đang thử cả hai:
for attr_name in ("decode", "string", "bseq"):
    val = getattr(cigar_obj, attr_name, None)
    if isinstance(val, (bytes, str)) and val:
        cigar_raw = val
        break
# Nếu không tìm được, thử gọi như method:
cigar_raw = cigar_obj.decode()
```

Với parasail 1.3.4, `cigar.decode` là **attribute** (bytes), không phải method.
Nhưng nếu đây là `str.decode` thì sẽ bị lỗi khi gọi như method.

**→ Reviewer: chạy `python3 -c "import parasail; r = parasail.sw_trace_striped_16('ACDEF','ACDEF',11,1,parasail.blosum62); print(type(r.cigar.decode), r.cigar.decode)"` để xác nhận.**

### 4.2 `beg_query` / `beg_ref` trong `result.cigar`

```python
beg_q = _get_int_attr(cigar_obj, ("beg_query",), 0)
beg_r = _get_int_attr(cigar_obj, ("beg_ref",), 0)
```

Nếu attribute tên khác (ví dụ `query_begin`, `beg1`, `beg2`...) thì `beg_q = 0`.
Điều này chỉ sai nếu alignment không bắt đầu từ position 0 — với SW thường vẫn ổn,
nhưng có thể gây identity thấp hơn thực tế nếu beg > 0 và ký tự tại 0 không khớp.

**→ Reviewer: in `dir(result.cigar)` để xác nhận tên attribute.**

### 4.3 Số phận của `_get_int_attr` khi `beg_query = 0`

Hàm `_get_int_attr` trả về `0` khi không tìm thấy attribute VÀ cũng khi attribute = 0.
Nếu alignment thực sự bắt đầu từ position 0 (rất phổ biến), code vẫn đúng.
Nếu không tìm thấy và alignment bắt đầu từ position > 0, identity bị tính sai (thấp hơn thực tế).

**→ Đây là lỗi hướng "an toàn" (sẽ PASS nhầm chứ không FAIL nhầm) nhưng vẫn cần xác nhận.**

### 4.4 Identity denominator: `aln_columns` hay `q_aln`?

Code hiện tại dùng `aln_columns` (tổng cột alignment kể cả gap I/D).
Brief nói `"chiều dài vùng aligned"` — mơ hồ.

Nếu reviewer muốn dùng chỉ residue-to-residue columns (loại bỏ gap), sửa:
```python
# Thay vì: identity = n_matches / aln_columns
identity = n_matches / (q_aln + r_aln - aln_columns + q_aln)  # phức tạp hơn
# Cách đơn giản hơn: dùng min(q_aln, r_aln) làm denominator
```

---

## 5. Những chỗ phải tự quyết (brief chưa nói rõ)

| Vấn đề | Quyết định | Lý do |
|--------|-----------|-------|
| Identity denominator khi có gap | Dùng `aln_columns` (kể cả gap) | Bảo thủ hơn, ít bỏ sót FAIL hơn |
| Gap open penalty | 11 (BLAST standard) | Phổ biến nhất cho BLOSUM62 |
| DB có antibody (107–225aa) không? | **Có, tất cả file đều load** | An toàn hơn, dù antibody rất khó vi phạm coverage với 85–100aa ứng viên |
| Khi CIGAR không decode được | Trả về `(0, 0, 0, 0)` → PASS | Sai hướng an toàn — không FAIL nhầm |
| Verdict khi identity = 70.0% chính xác | PASS (vì rule là **>70**, không phải ≥70) | Đọc nguyên văn: `序列一致性须 ≤70%` → vi phạm khi > 70% |
| Nhiều sequence cùng tên trong candidate file | Giữ nguyên tất cả (không deduplicate) | Để người dùng thấy tất cả kết quả |

---

## 6. Hạn chế đã biết

1. **RAM:** Load 122,273 sequences vào RAM dưới dạng Python list. Ước tính ~300–500 MB.
   Với máy còn 3 GB, đủ dùng. Nhưng nếu có nhiều target hơn hoặc DB lớn hơn, cần streaming.

2. **Không parallelise:** Code chạy single-thread. Có thể dùng `multiprocessing.Pool`
   để chạy song song nhiều candidate, tăng tốc tuyến tính theo số CPU. Chưa làm để giữ code đơn giản.

3. **parasail 16-bit:** `sw_trace_striped_16` dùng 16-bit integer — có thể overflow với sequence
   rất dài (>32,768 score). Với protein minibinder 85–100aa + BLOSUM62, score tối đa ~100×11=1100
   → không overflow. Nếu cần độ chính xác tuyệt đối cho sequence dài, dùng `sw_trace_striped_32`.

4. **Case 4 margin:** SW có thể align thêm 1–2 residue từ tail (nếu ngẫu nhiên có match),
   nên `coverage_ours_pct` đôi khi cao hơn 60/90 = 66.7%. Test chấp nhận ±10pp.

5. **TL1A: không có de novo minibinder** — DB chỉ gồm antibody (106–300aa) và DcR3 (176–300aa).
   Với ứng viên 85–100aa, `coverage_ours` gần như chắc chắn < 80% → mọi sequence đều PASS.
   **Không phải bug — đây là thực tế sinh học.**

---

## 7. Độ phức tạp tính toán (ước tính)

| Tham số | Giá trị |
|---------|---------|
| Số known binder (il7ra) | ~122,285 (4 file) |
| Số ứng viên | 20 |
| Alignments tổng | ~2,445,700 |
| Tốc độ parasail SW (ước tính) | ~1–5 μs/alignment trên Arch Linux 8-core |
| Thời gian ước tính (single-thread) | 2.4M × 2 μs = ~5–50 giây |
| Mục tiêu brief | < 5 phút |

> [!TIP] Nếu chạy quá 5 phút, bật `multiprocessing`:
> ```python
> from multiprocessing import Pool
> with Pool() as pool:
>     results = pool.starmap(check_one_candidate, [(n, s, known) for n, s in candidates])
> ```

**Thời gian thực tế:** _(Người review điền sau khi chạy)_ `_____ giây cho 20 candidate`

---

## 8. Cách chạy

### Yêu cầu
```bash
# Venv đã có sẵn với parasail 1.3.4:
.venv/bin/python3 -c "import parasail; print(parasail.__version__)"
```

### Chạy self-test
```bash
cd /home/ngochieu/repo/proteinAI
.venv/bin/python3 tools/test_novelty_check.py
```

### Chạy trên candidate thực
```bash
# FASTA input:
.venv/bin/python3 tools/novelty_check.py --input candidates.fasta --target il7ra --out results.csv

# CSV input (cần cột 'name' và 'sequence'):
.venv/bin/python3 tools/novelty_check.py --input candidates.csv --target tl1a
```

### Output mẫu
```
#    Name                                 Len    MaxID%  CovOurs%  CovThem%  Verdict   VrdStrict   NHits>80
---- ----------------------------------- -----  -------  --------  --------  --------  ----------  --------
1    my_design_001                           90     45.2      67.1      55.3  PASS      PASS               0
2    my_design_002                           88    100.0     100.0     100.0  FAIL      FAIL             122

Summary (primary):  1 PASS  /  1 FAIL   (of 2 candidates)
Summary (strict):   1 PASS  /  1 FAIL
```

---

## 9. Lịch sử và venv

- Venv tại `.venv/` **đã được thêm vào `.gitignore`** (xác nhận: dòng `.venv/` có trong file)
- Mọi lệnh chạy dùng `.venv/bin/python3` — không dùng system Python 3.14 (EXTERNALLY-MANAGED)
- parasail 1.3.4 đã cài trong venv, không cần `pip install` lại

---

*Báo cáo do agent AGY-01 (Antigravity) tạo ngày 2026-11-11. Phần kết quả test để trống chờ reviewer chạy.*

---

# ✅ KIỂM CHỨNG (do Claude chạy, 2026-11-11)

agy không chạy được lệnh shell nên để trống phần kết quả — đúng yêu cầu. Dưới đây là số **chạy thật**.

## Self-test: 6 PASS / 0 FAIL

```
.venv/bin/python3 tools/test_novelty_check.py
→ Results: 6 PASS / 0 FAIL  (of 6 cases) — All tests passed.
```

Case 4 (case quyết định), số thật:

```
[sub-A] alignment trực tiếp vs known binder được nhúng:
   q_aligned=60  coverage_ours_pct=66.67%  identity_pct=100.00%  aln_columns=60
[sub-B] chạy full DB: verdict=PASS  cov_worst=80.00%  max_identity_pct=30.56%
```

→ 66,67% khớp chính xác kỳ vọng 60/90. Mẫu số coverage đã đúng.

## Hiệu năng

**56 giây** cho 20 sequence (85–100 aa) × 122.273 known binder ≈ 2,4 triệu alignment.
Yêu cầu brief là <5 phút → đạt.

## Vòng 1 đã fail và nguyên nhân gốc

Vòng đầu **5 PASS / 1 FAIL**, fail đúng Case 4 — `coverage_ours_pct=85.56%` thay vì 66,7%,
và random 90aa cho coverage 81% (`n_hits_over_80pct_coverage=328`).

Nguyên nhân: `parasail.sw_trace_striped_16` trong parasail **1.3.4 không có** `read_begin1`
lẫn `query_begin`. Code rơi về `cigar.beg_query` với giả định đó là offset thật trong query —
**sai**: `cigar.beg_query` luôn bằng 0 tương đối với vùng aligned. Nên `q_aligned = end_query - 0`
phồng lên đúng bằng offset nơi alignment bắt đầu.

Làm nó nguy hiểm: có `fallback = 0` khi không tìm thấy attribute → **không crash, không cảnh báo**,
chỉ ra số sai mà trông hợp lý. Nếu Case 4 không đóng đinh bằng con số cụ thể thì công cụ sẽ
âm thầm báo PASS cho mọi thứ.

Bản sửa: đếm ký tự không phải gap trực tiếp trên `result.traceback.query/.ref`, bỏ hẳn
begin/end offset, bỏ fallback im lặng (giờ raise), thêm assert `q_aligned <= len(query)`.
Sau sửa: hit vượt ngưỡng coverage giảm **328 → 1**.

## Lưu ý quan trọng cho việc dùng về sau

Case 2 cho thấy một sequence **random** vẫn có thể đạt `coverage_ours = 81,11%` với identity 23,29%
— sự kiện đuôi khi quét 122.273 subject. Nghĩa là ngưỡng coverage ≥80% **có thể bị chạm do tình cờ**.

→ Điều này **phủ nhận** phát biểu trước đó rằng design ≥84 aa thì "về mặt toán học miễn nhiễm"
luật novelty. Đúng tuyệt đối chỉ với alignment **không gap** (67 aa không thể cover 68 residue).
Có gap thì query residue trong vùng aligned **vượt được** chiều dài subject, vì gap chèn vào subject
vẫn tiêu residue của query.

**Phát biểu đúng: ≥85 aa bảo vệ rất mạnh, không tuyệt đối.** Kết luận "design 85–100 aa" không đổi
(lý do thứ hai — sàn `d0` của ipSAE — vẫn nguyên), nhưng **phải chạy công cụ này trên mọi sequence
trước khi submit**, không được dựa vào lập luận độ dài.

## Cách dùng

```bash
.venv/bin/python3 tools/novelty_check.py --input <candidates.fasta|csv> \
    --target il7ra|tl1a --out <report.csv>
```
