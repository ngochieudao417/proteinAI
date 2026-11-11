#!/usr/bin/env python3
"""
build_submission.py — Submission builder + validator for Changjiang Cup 2026.

Builds a submission CSV for a given target, validates it against the site's
reference validator logic, and refuses to write the file if validation fails.

Usage:
    python3 tools/build_submission.py --target il7ra --input <candidates> --out out/il7ra_submission.csv
    python3 tools/build_submission.py --target tl1a  --input <candidates> --out out/tl1a_submission.csv \
        [--novelty-report novelty_results.csv] [--force]

Targets supported: il7ra, tl1a  (both minibinder track, 40–100 aa)
"""

from __future__ import annotations

import argparse
import csv
import re
import statistics
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# Constants — mirrored from data/validate_submission_csv.py (reference)
# ---------------------------------------------------------------------------

ALLOWED_AA: set[str] = set("ACDEFGHIKLMNPQRSTVWY")

# Minibinder track: 40–100 aa, bounds inclusive
MIN_LEN = 40
MAX_LEN = 100
MAX_SEQUENCES = 20

TARGETS = {
    "il7ra": "IL-7Rα / CD127",
    "tl1a":  "TL1A / TNFSF15",
}

# Tag motifs that suggest the sequence is not a mature protein (warn, never auto-remove)
TAG_MOTIFS: list[tuple[str, str]] = [
    (r"H{6,}",        "His-tag (6+ His)"),
    (r"DYKDDDDK",     "FLAG tag"),
    (r"EQKLISEEDL",   "Myc tag"),
    (r"(GGGGS){2,}",  "GS linker (GGGGS repeat, ≥2×)"),
    (r"GGGGS",        "GS linker (GGGGS at terminus region)"),
    (r"RSHHHHHH",     "His-tag with RS"),
    (r"HHHHHHHH",     "His-tag (8+ His)"),
    (r"HHHHHH",       "His-tag (6 His)"),
]

# How many characters at each end to scan for tag motifs (N-term / C-term)
TAG_SCAN_WINDOW = 30

# ---------------------------------------------------------------------------
# Reference-compatible CSV splitter
# (faithful copy of split_csv from data/validate_submission_csv.py)
# ---------------------------------------------------------------------------

def _split_csv_line(line: str) -> list[str]:
    """Split a single CSV line, handling quoted fields (RFC 4180 style)."""
    out: list[str] = []
    cur = ""
    q = False
    i = 0
    while i < len(line):
        c = line[i]
        if q:
            if c == '"':
                if i + 1 < len(line) and line[i + 1] == '"':
                    cur += '"'
                    i += 1
                else:
                    q = False
            else:
                cur += c
        else:
            if c == '"':
                q = True
            elif c == ',':
                out.append(cur)
                cur = ''
            else:
                cur += c
        i += 1
    out.append(cur)
    return out


# ---------------------------------------------------------------------------
# Reference-compatible validator
# (faithful re-implementation of validate() from data/validate_submission_csv.py)
# ---------------------------------------------------------------------------

def _validate_reference(text: str) -> tuple[int, list[tuple[int, str, str]]]:
    """
    Validate CSV text using the EXACT same logic as the site's reference validator.

    Args:
        text: Full CSV file contents as a string (UTF-8, BOM already stripped by caller).

    Returns:
        (n_sequences_parsed, list_of_errors)
        Each error is (line_number_1indexed, name, reason_string).
        Empty error list means the file would be accepted by the site.

    NOTE: reason strings are kept in Chinese to exactly match the reference validator's
    output format (data/validate_submission_csv.py uses Chinese error strings).
    """
    maxrows = MAX_SEQUENCES
    mn = MIN_LEN
    mx = MAX_LEN
    req = "minibinder 为 40—100 aa"

    text = text.lstrip('\ufeff')           # strip BOM if any (reference does this)
    lines = re.split(r'\r?\n', text)
    nonblank = [l for l in lines if l.strip()]
    errs: list[tuple[int, str, str]] = []

    if len(nonblank) > maxrows + 1:
        return len(nonblank) - 1, [(len(nonblank), '', f'最多提交 {maxrows} 条序列（含表头共 {maxrows + 1} 行）')]

    hdr = [h.strip().lower() for h in _split_csv_line(lines[0] if lines else '')]
    if 'name' not in hdr or 'sequence' not in hdr:
        return 0, [(1, '', '表头必须包含 name 和 sequence')]

    oi = hdr.index('name')
    ui = hdr.index('sequence')
    seen: set[str] = set()
    m = 0

    for g in range(1, len(lines)):
        y = lines[g]
        if not y.strip():
            continue
        j = _split_csv_line(y)
        x = (j[oi] if oi < len(j) else '').strip()
        S = re.sub(r'\s+', '', (j[ui] if ui < len(j) else '')).upper()
        p = g + 1
        if not x and not S:
            continue
        m += 1
        if m > maxrows:
            return m, [(p, '', f'最多提交 {maxrows} 条序列（含表头共 {maxrows + 1} 行）')]
        if not x or len(x) > 64:
            errs.append((p, x, '名称须为 1 到 64 个字符'))
            continue
        if x in seen:
            errs.append((p, x, '名称重复'))
            continue
        seen.add(x)
        if not S:
            errs.append((p, x, '序列为空'))
            continue
        if any(ch not in ALLOWED_AA for ch in S):
            errs.append((p, x, '含有非法氨基酸'))
            continue
        if len(S) < mn or len(S) > mx:
            errs.append((p, x, f'序列长度不符合要求：{req}'))

    if m == 0:
        errs.append((1, '', '没有序列'))

    return m, errs


# ---------------------------------------------------------------------------
# Our own validator — same logic, English messages for user-friendliness
# ---------------------------------------------------------------------------

def _validate_our(
    records: list[tuple[str, str]],
) -> list[tuple[int, str, str]]:
    """
    Validate a list of (name, sequence) pairs using our own English-language messages.
    Returns list of (record_index_1based, name, reason) for each error found.

    This runs BEFORE writing the CSV. It operates on the in-memory list, not on
    the serialised text, so it catches problems earlier and gives cleaner messages.
    """
    errs: list[tuple[int, str, str]] = []

    if len(records) > MAX_SEQUENCES:
        errs.append((0, '', f'Too many sequences: {len(records)} submitted, max is {MAX_SEQUENCES}.'))
        return errs   # no point continuing

    if not records:
        errs.append((0, '', 'No sequences provided.'))
        return errs

    seen_names: set[str] = set()

    for idx, (name, raw_seq) in enumerate(records, 1):
        # Normalise sequence exactly as the reference does: strip whitespace, uppercase
        seq = re.sub(r'\s+', '', raw_seq).upper()

        # Validate name length: 1–64 chars
        if not name or len(name.strip()) == 0:
            errs.append((idx, name, 'Name is empty (must be 1–64 characters).'))
            continue
        name_stripped = name.strip()
        if len(name_stripped) > 64:
            errs.append((idx, name, f'Name is too long ({len(name_stripped)} chars); max is 64.'))
            continue

        # Duplicate name check
        if name_stripped in seen_names:
            errs.append((idx, name_stripped, 'Duplicate name.'))
            continue
        seen_names.add(name_stripped)

        # Empty sequence
        if not seq:
            errs.append((idx, name_stripped, 'Sequence is empty after stripping whitespace.'))
            continue

        # Invalid amino acid characters
        bad_chars = sorted(set(seq) - ALLOWED_AA)
        if bad_chars:
            errs.append((idx, name_stripped,
                         f'Invalid amino acid character(s): {bad_chars}. '
                         f'Only standard 20 AA allowed: ACDEFGHIKLMNPQRSTVWY.'))
            continue

        # Length check: 40–100 aa (inclusive, minibinder track)
        if len(seq) < MIN_LEN or len(seq) > MAX_LEN:
            errs.append((idx, name_stripped,
                         f'Sequence length {len(seq)} aa is outside allowed range '
                         f'[{MIN_LEN}, {MAX_LEN}] (minibinder track, bounds inclusive).'))

    return errs


# ---------------------------------------------------------------------------
# Tag / contamination warnings
# ---------------------------------------------------------------------------

def _check_tags(name: str, seq: str) -> list[str]:
    """
    Return a list of warning messages if the sequence appears to contain
    purification tags or fusion motifs.

    We only check a window at the N-terminus and C-terminus (first/last
    TAG_SCAN_WINDOW residues) for some motifs, and the full sequence for others.
    We NEVER modify the sequence — warnings only.
    """
    warnings: list[str] = []
    seq_upper = seq.upper()
    n_term = seq_upper[:TAG_SCAN_WINDOW]
    c_term = seq_upper[-TAG_SCAN_WINDOW:]

    for pattern, label in TAG_MOTIFS:
        # Check full sequence for all motifs (a tag anywhere is suspicious)
        if re.search(pattern, seq_upper):
            warnings.append(f'  ⚠  [{name}] Possible {label} detected in sequence.')

    return warnings


# ---------------------------------------------------------------------------
# Duplicate sequence warnings
# ---------------------------------------------------------------------------

def _check_duplicates(records: list[tuple[str, str]]) -> list[str]:
    """Return warning messages for any pairs of identical sequences (100% identity)."""
    warnings: list[str] = []
    # Normalise sequences for comparison
    normed = [(name, re.sub(r'\s+', '', seq).upper()) for name, seq in records]

    seen: dict[str, str] = {}  # seq -> first name
    for name, seq in normed:
        if seq in seen:
            warnings.append(
                f'  ⚠  Duplicate sequence: [{name}] is 100% identical to [{seen[seq]}]. '
                f'Competition rules: only the first submission is kept for identical sequences.'
            )
        else:
            seen[seq] = name

    return warnings


# ---------------------------------------------------------------------------
# Summary table
# ---------------------------------------------------------------------------

def _print_summary(records: list[tuple[str, str]], warnings: list[str]) -> None:
    """Print a summary table of the sequences to be submitted."""
    normed_seqs = [re.sub(r'\s+', '', seq).upper() for _, seq in records]
    lengths = [len(s) for s in normed_seqs]

    print()
    print("=" * 60)
    print("  SUBMISSION SUMMARY")
    print("=" * 60)
    print(f"  Sequences    : {len(records)}")
    if lengths:
        print(f"  Length min   : {min(lengths)} aa")
        print(f"  Length max   : {max(lengths)} aa")
        if len(lengths) >= 2:
            print(f"  Length median: {statistics.median(lengths):.1f} aa")
        else:
            print(f"  Length median: {lengths[0]:.1f} aa")
    if warnings:
        print(f"  Warnings     : {len(warnings)}")
        for w in warnings:
            print(w)
    else:
        print("  Warnings     : none")
    print("=" * 60)
    print()


# ---------------------------------------------------------------------------
# CSV serialisation
# ---------------------------------------------------------------------------

def _build_csv_text(records: list[tuple[str, str]]) -> str:
    """
    Serialise records to a CSV string with header 'name,sequence'.
    Sequences are stripped of whitespace and uppercased before writing.
    UTF-8, no BOM, CRLF line endings (safest for cross-platform validators).
    """
    lines: list[str] = ["name,sequence"]
    for name, raw_seq in records:
        seq = re.sub(r'\s+', '', raw_seq).upper()
        # Quote fields that contain commas or quotes
        name_field = _csv_quote(name)
        seq_field  = _csv_quote(seq)
        lines.append(f"{name_field},{seq_field}")
    return "\r\n".join(lines) + "\r\n"


def _csv_quote(value: str) -> str:
    """Minimal CSV quoting: only quote when necessary."""
    if ',' in value or '"' in value or '\n' in value or '\r' in value:
        return '"' + value.replace('"', '""') + '"'
    return value


# ---------------------------------------------------------------------------
# Novelty report parsing
# ---------------------------------------------------------------------------

def _load_novelty_report(path: Path) -> dict[str, str]:
    """
    Parse a novelty_check.py CSV output.
    Returns a dict mapping sequence name -> verdict ('PASS' or 'FAIL').
    """
    result: dict[str, str] = {}
    with open(path, newline='', encoding='utf-8-sig') as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames is None:
            raise ValueError(f"Novelty report '{path}' appears to be empty.")
        fields_lower = [f.strip().lower() for f in reader.fieldnames]
        if 'name' not in fields_lower or 'verdict' not in fields_lower:
            raise ValueError(
                f"Novelty report '{path}' must contain 'name' and 'verdict' columns. "
                f"Found: {list(reader.fieldnames)}"
            )
        # Map normalised field name back to original
        orig_name    = reader.fieldnames[fields_lower.index('name')]
        orig_verdict = reader.fieldnames[fields_lower.index('verdict')]
        for row in reader:
            n = row.get(orig_name, '').strip()
            v = row.get(orig_verdict, '').strip().upper()
            if n:
                result[n] = v
    return result


def _check_novelty(
    records: list[tuple[str, str]],
    novelty_map: dict[str, str],
) -> list[str]:
    """
    Return a list of FAIL names for sequences whose verdict is FAIL in the novelty report.
    """
    fails: list[str] = []
    for name, _ in records:
        verdict = novelty_map.get(name)
        if verdict == 'FAIL':
            fails.append(name)
    return fails


# ---------------------------------------------------------------------------
# Input parsing (CSV or FASTA)
# ---------------------------------------------------------------------------

def _parse_fasta(path: Path) -> list[tuple[str, str]]:
    """Parse a FASTA file. Returns list of (header, sequence)."""
    records: list[tuple[str, str]] = []
    header: str | None = None
    seq_parts: list[str] = []
    with open(path, 'r', encoding='utf-8-sig') as fh:
        for raw in fh:
            line = raw.rstrip('\r\n')
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    records.append((header, ''.join(seq_parts)))
                header = line[1:].strip()
                seq_parts = []
            else:
                seq_parts.append(line.strip())
    if header is not None:
        records.append((header, ''.join(seq_parts)))
    return records


def _parse_input_csv(path: Path) -> list[tuple[str, str]]:
    """
    Parse a CSV input file with flexible column name matching:
        name column:     'name', 'id', 'seq_id', 'sequence_id'
        sequence column: 'sequence', 'seq', 'aa', 'protein'
    Returns list of (name, raw_sequence).
    """
    with open(path, newline='', encoding='utf-8-sig') as fh:
        sample = fh.read(4096)
        fh.seek(0)
        delimiter = '\t' if '\t' in sample else ','
        reader = csv.DictReader(fh, delimiter=delimiter)
        if reader.fieldnames is None:
            raise ValueError(f"CSV '{path}' appears to be empty.")

        fields_lower = [f.strip().lower() for f in reader.fieldnames]

        name_key = _find_col(fields_lower, ('name', 'id', 'seq_id', 'sequence_id'))
        seq_key  = _find_col(fields_lower, ('sequence', 'seq', 'aa', 'protein'))

        if name_key is None or seq_key is None:
            raise ValueError(
                f"CSV '{path}' must have name and sequence columns. "
                f"Recognised name cols: name/id/seq_id/sequence_id. "
                f"Recognised seq cols: sequence/seq/aa/protein. "
                f"Found: {list(reader.fieldnames)}"
            )

        orig_name = reader.fieldnames[fields_lower.index(name_key)]
        orig_seq  = reader.fieldnames[fields_lower.index(seq_key)]

        return [
            (row[orig_name].strip(), row[orig_seq].strip())
            for row in reader
            if row.get(orig_name, '').strip() or row.get(orig_seq, '').strip()
        ]


def _find_col(fields: list[str], candidates: tuple[str, ...]) -> str | None:
    for c in candidates:
        if c in fields:
            return c
    return None


def parse_input(path: Path) -> list[tuple[str, str]]:
    """Auto-detect format and parse candidate sequences."""
    suffix = path.suffix.lower()
    if suffix in ('.fa', '.fasta', '.faa', '.fna'):
        return _parse_fasta(path)
    elif suffix in ('.csv', '.tsv'):
        return _parse_input_csv(path)
    else:
        with open(path, 'r', encoding='utf-8-sig') as fh:
            first = fh.read(1)
        if first == '>':
            return _parse_fasta(path)
        return _parse_input_csv(path)


# ---------------------------------------------------------------------------
# Cross-check: run reference validator on the generated CSV text
# ---------------------------------------------------------------------------

def _cross_validate(csv_text: str) -> list[tuple[int, str, str]]:
    """
    Run the reference validator (from data/validate_submission_csv.py) on the CSV text.
    Returns the error list from the reference validator.
    """
    # Import inline so this file can be used without the data/ directory on PYTHONPATH
    repo_root = Path(__file__).resolve().parent.parent
    sys.path.insert(0, str(repo_root / 'data'))
    try:
        import validate_submission_csv as ref_validator  # type: ignore
        _, ref_errs = ref_validator.validate(csv_text, 20, 40, 100, 'minibinder 为 40—100 aa')
        return ref_errs
    except ImportError as e:
        print(f"[WARN] Could not import reference validator: {e}", file=sys.stderr)
        return []
    finally:
        if str(repo_root / 'data') in sys.path:
            sys.path.remove(str(repo_root / 'data'))


# ---------------------------------------------------------------------------
# Main build function
# ---------------------------------------------------------------------------

def build_submission(
    target: str,
    records: list[tuple[str, str]],
    out_path: Path,
    novelty_map: dict[str, str] | None = None,
    force: bool = False,
) -> bool:
    """
    Build and validate a submission CSV.

    Args:
        target:      Target key ('il7ra' or 'tl1a').
        records:     List of (name, raw_sequence) pairs from the input file.
        out_path:    Path to write the output CSV.
        novelty_map: Dict of name -> verdict from novelty_check.py output (optional).
        force:       If True, build even if novelty check FAILs.

    Returns:
        True if the file was written successfully, False otherwise.
    """
    target_label = TARGETS.get(target, target.upper())
    print(f"\n[BUILD] Target: {target_label} ({target})")
    print(f"[BUILD] Input sequences: {len(records)}")

    # ── Step 1: Novelty check ──────────────────────────────────────────────
    if novelty_map is not None:
        fails = _check_novelty(records, novelty_map)
        if fails:
            print("\n[ERROR] The following sequences have FAIL verdict in the novelty report:")
            for name in fails:
                print(f"  ✗  {name}")
            if not force:
                print("\n[ABORT] Build refused. Use --force to override novelty check failures.")
                return False
            else:
                print("\n[WARN] --force flag is set; proceeding despite novelty FAIL(s).")

    # ── Step 2: Our own validation ─────────────────────────────────────────
    our_errs = _validate_our(records)
    if our_errs:
        print("\n[ERROR] Validation failed:")
        for (idx, name, reason) in our_errs:
            if idx == 0:
                print(f"  ✗  (file-level) {reason}")
            else:
                print(f"  ✗  row {idx} [{name}]: {reason}")
        print("\n[ABORT] File not written.")
        return False

    # ── Step 3: Collect warnings (non-fatal) ──────────────────────────────
    warnings: list[str] = []

    # Duplicate sequence warning
    warnings.extend(_check_duplicates(records))

    # Tag motif warnings
    for name, raw_seq in records:
        seq = re.sub(r'\s+', '', raw_seq).upper()
        warnings.extend(_check_tags(name, seq))

    # ── Step 4: Summary table ──────────────────────────────────────────────
    _print_summary(records, warnings)

    # ── Step 5: Serialise CSV ─────────────────────────────────────────────
    csv_text = _build_csv_text(records)

    # Size check (< 10 MB)
    size_bytes = len(csv_text.encode('utf-8'))
    if size_bytes > 10 * 1024 * 1024:
        print(f"[ERROR] CSV size {size_bytes / 1024 / 1024:.2f} MB exceeds 10 MB limit.")
        print("[ABORT] File not written.")
        return False

    # ── Step 6: Cross-validate with reference validator ───────────────────
    ref_errs = _cross_validate(csv_text)
    if ref_errs:
        print("\n[ERROR] Reference validator (data/validate_submission_csv.py) REJECTED the output:")
        for (line_no, name, reason) in ref_errs:
            print(f"  ✗  line {line_no} [{name}]: {reason}")
        print(
            "\n[ERROR] Our validator and the reference validator DISAGREE — this is a bug. "
            "The file has NOT been written. Please report this discrepancy."
        )
        return False

    print("[OK] Reference validator (data/validate_submission_csv.py) also PASSES. ✓")

    # ── Step 7: Write file ────────────────────────────────────────────────
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(csv_text, encoding='utf-8')

    print(f"[OK] Written: {out_path}  ({size_bytes} bytes, {len(records)} sequences)")
    return True


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=(
            "Build and validate a competition submission CSV for the Changjiang Cup 2026.\n"
            "Validates using the same logic as the site's reference validator before writing."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python3 tools/build_submission.py --target il7ra --input candidates.fasta --out out/il7ra_submission.csv
  python3 tools/build_submission.py --target tl1a  --input candidates.csv   --out out/tl1a_submission.csv
  python3 tools/build_submission.py --target il7ra --input candidates.csv   \\
      --novelty-report novelty_results.csv --out out/il7ra_submission.csv
        """,
    )
    p.add_argument('--target', '-t', required=True, choices=list(TARGETS),
                   help='Target protein (il7ra or tl1a).')
    p.add_argument('--input', '-i', required=True,
                   help='Input file with candidate sequences (.fasta or .csv).')
    p.add_argument('--out', '-o', required=True,
                   help='Output CSV path.')
    p.add_argument('--novelty-report', metavar='CSV',
                   help='Output CSV from tools/novelty_check.py. Sequences with verdict=FAIL are rejected.')
    p.add_argument('--force', action='store_true',
                   help='Build even if some sequences have novelty verdict=FAIL.')
    return p


def main() -> None:
    args = _build_parser().parse_args()

    input_path = Path(args.input)
    if not input_path.exists():
        sys.exit(f"[ERROR] Input file not found: {input_path}")

    out_path = Path(args.out)

    # Parse input
    try:
        records = parse_input(input_path)
    except Exception as e:
        sys.exit(f"[ERROR] Failed to parse input file: {e}")

    if not records:
        sys.exit(f"[ERROR] No sequences found in {input_path}.")

    # Load novelty report if provided
    novelty_map: dict[str, str] | None = None
    if args.novelty_report:
        novelty_path = Path(args.novelty_report)
        if not novelty_path.exists():
            sys.exit(f"[ERROR] Novelty report not found: {novelty_path}")
        try:
            novelty_map = _load_novelty_report(novelty_path)
            print(f"[INFO] Loaded novelty report: {len(novelty_map)} entries from {novelty_path}")
        except Exception as e:
            sys.exit(f"[ERROR] Failed to load novelty report: {e}")

    success = build_submission(
        target=args.target,
        records=records,
        out_path=out_path,
        novelty_map=novelty_map,
        force=args.force,
    )
    sys.exit(0 if success else 1)


if __name__ == '__main__':
    main()
