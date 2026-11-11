# ═══════════════════════════════════════════════════════════════════════
# PATCH cho notebook 02 — lỗi "A3M match-column count mismatch: 219 vs 220"
#
# CÁCH DÙNG: thêm cell MỚI ngay TRƯỚC cell 10 (cell MSA), dán nguyên file này,
#            chạy nó, rồi chạy lại cell 10 và tiếp tục bình thường.
#            KHÔNG cần chạy lại từ đầu — model đã load, weights đã cache.
#
# NGUYÊN NHÂN (đã chẩn đoán trên bản a3m ở local, không đoán):
#   api.colabfold.com kết thúc file uniref.a3m bằng MỘT null byte \x00.
#   Nó dính vào dòng CUỐI của MSA, khiến dòng đó đếm ra 220 match column
#   thay vì 219 → MSA.from_a3m raise ValueError.
#
#   Đo thật trên 2 file:
#     IL7Ra_P16871_E21-D239.uniref.a3m : 921 record, đúng 1 dòng sai (idx 920)
#        đuôi dòng đó: ...'WSDWSRVEEYTTLEA--------\x00'
#     TL1A_O95150_L72-L251.uniref.a3m  : 2037 record, cũng đúng 1 null byte
#
#   Vì sao validation của notebook không bắt được: hàm a3m_stats() chỉ đếm
#   match-column của dòng QUERY (dòng đầu). Null byte nằm ở dòng cuối nên lọt.
# ═══════════════════════════════════════════════════════════════════════

from esm.utils.msa import MSA as _MSA

_ALLOWED_CTRL = (9, 10, 13)          # tab, LF, CR
_orig_from_a3m = _MSA.from_a3m


def _sanitize_a3m(path):
    """Bỏ mọi byte điều khiển khỏi file a3m. Trả về số byte đã bỏ."""
    with open(path, "rb") as f:
        raw = f.read()
    clean = bytes(b for b in raw if b >= 32 or b in _ALLOWED_CTRL)
    if clean != raw:
        with open(path, "wb") as f:
            f.write(clean)
    return len(raw) - len(clean)


def _audit_a3m(path):
    """Đếm match-column của MỌI dòng, không chỉ query. Trả về (depth, {count: n})."""
    from collections import Counter
    recs, buf, started = [], [], False
    with open(path, "r", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("#"):
                continue
            if line.startswith(">"):
                if started:
                    recs.append("".join(buf))
                buf, started = [], True
            elif line.strip():
                buf.append(line.strip())
    if started:
        recs.append("".join(buf))
    counts = Counter(sum(1 for c in s if not c.islower()) for s in recs)
    return len(recs), counts


def _patched_from_a3m(*args, **kwargs):
    path = kwargs.get("path") or (args[0] if args else None)
    if isinstance(path, str):
        removed = _sanitize_a3m(path)
        depth, counts = _audit_a3m(path)
        tag = path.split("/")[-1]
        if removed:
            print(f"  [PATCH] {tag}: bỏ {removed} byte điều khiển")
        if len(counts) == 1:
            mc = next(iter(counts))
            print(f"  [PATCH] {tag}: depth={depth}, match_columns={mc} (đồng nhất) ✅")
        else:
            print(f"  [PATCH] {tag}: ⚠️ match_columns KHÔNG đồng nhất: {dict(counts)}")
            print(f"  [PATCH] {tag}: vẫn còn dòng lệch — báo lại cho Claude, đừng chạy tiếp")
    return _orig_from_a3m(*args, **kwargs)


_MSA.from_a3m = _patched_from_a3m

print("=" * 68)
print("PATCH đã áp dụng: MSA.from_a3m giờ tự làm sạch a3m trước khi parse")
print("và audit match-column của MỌI dòng (không chỉ query).")
print()
print("→ Giờ chạy lại CELL 10. Không cần chạy lại từ đầu.")
print("=" * 68)
