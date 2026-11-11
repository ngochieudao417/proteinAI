#!/usr/bin/env python3
"""
novelty_check.py — Novelty checker for protein minibinder design submissions.

Checks candidate sequences against known binders using Smith-Waterman local
alignment (via parasail). Reports coverage and identity metrics, and verdicts
based on competition rules (Changjiang Cup 2026):

  FAIL if: coverage_ours_pct >= 80  AND  identity_pct > 70

where:
  coverage_ours_pct = (query residues in aligned region) / (query total length) × 100
  identity_pct      = (identical positions in aligned region) / (aligned columns) × 100

Usage:
  .venv/bin/python3 tools/novelty_check.py \\
      --input <file.csv|file.fasta> --target il7ra|tl1a [--out <output.csv>]

IMPORTANT: Run with .venv/bin/python3 (parasail must be installed in venv).
"""

import argparse
import csv
import sys
import time
from pathlib import Path

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

REPO_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = REPO_ROOT / "data" / "known_binders"

# Known binder database files, keyed by target.
KNOWN_BINDER_FILES: dict[str, list[tuple[str, str]]] = {
    "il7ra": [
        # (filename, label)
        ("IL7Ra_denovo_minibinder_designs_IPD_Cao2022.fasta", "IPD_Cao2022_122k"),
        ("IL7Ra_denovo_minibinders_VALIDATED.fasta",          "VALIDATED"),
        ("IL7Ra_therapeutic_antibodies.fasta",                "therapeutics"),
        ("IL7Ra_PDB_antibody_Fabs.fasta",                     "PDB_Fabs"),
    ],
    "tl1a": [
        ("TL1A_therapeutic_antibodies.fasta", "therapeutics"),
        ("TL1A_natural_binder_DcR3.fasta",    "DcR3_natural"),
    ],
}

COVERAGE_THRESHOLD = 80.0  # percent
IDENTITY_THRESHOLD = 70.0  # percent — FAIL when STRICTLY GREATER than this

# Standard 20 amino acids (uppercase only; sequences are uppercased before checks)
VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")

# Smith-Waterman parameters: gap open 11, gap extend 1 — standard for protein SW with BLOSUM62
_GAP_OPEN   = 11
_GAP_EXTEND = 1

# ---------------------------------------------------------------------------
# parasail import — guard against missing installation
# ---------------------------------------------------------------------------

_parasail_import_error_msg = ""  # set below on failure

try:
    import parasail  # type: ignore

    # Verify that the specific symbols we use are present (paranoia check for version mismatches)
    if not hasattr(parasail, "sw_trace_striped_16"):
        raise ImportError("parasail.sw_trace_striped_16 not found — unexpected parasail version.")
    if not hasattr(parasail, "blosum62"):
        raise ImportError("parasail.blosum62 not found — unexpected parasail version.")

    _PARASAIL_OK = True

except ImportError as _parasail_import_error:
    _PARASAIL_OK = False
    _parasail_import_error_msg = str(_parasail_import_error)


def _assert_parasail() -> None:
    if not _PARASAIL_OK:
        sys.exit(
            f"[ERROR] parasail is not available: {_parasail_import_error_msg}\n"
            "Install it:  .venv/bin/pip install parasail\n"
            "Then run:    .venv/bin/python3 tools/novelty_check.py ..."
        )


# ---------------------------------------------------------------------------
# FASTA parsing — returns list (loads whole file; for 122k× we need that in RAM)
# ---------------------------------------------------------------------------

def parse_fasta(path: Path) -> list[tuple[str, str]]:
    """Parse a (possibly multi-line) FASTA file.

    Returns list of (header, sequence) tuples.
    Sequences are uppercased. Header is the full string after '>'.

    Raises FileNotFoundError if path doesn't exist.
    Raises ValueError on format errors.
    """
    records: list[tuple[str, str]] = []
    header: str | None = None
    seq_parts: list[str] = []

    with open(path, "r", encoding="utf-8") as fh:
        for raw in fh:
            line = raw.rstrip("\r\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_parts).upper()))
                header = line[1:].strip()
                seq_parts = []
            else:
                # Strip whitespace within sequence lines (should not appear in valid FASTA,
                # but some files include spaces or tabs)
                seq_parts.append(line.strip())

    if header is not None:
        records.append((header, "".join(seq_parts).upper()))

    return records


# ---------------------------------------------------------------------------
# CSV / FASTA candidate input parsing
# ---------------------------------------------------------------------------

def parse_candidates(path: Path) -> list[tuple[str, str]]:
    """Parse candidate sequences from a CSV or FASTA file.

    CSV must have columns for name and sequence; flexible column name matching:
      name column:     'name', 'id', 'seq_id', 'sequence_id'  (first match wins)
      sequence column: 'sequence', 'seq', 'aa', 'protein'     (first match wins)

    Returns list of (name, sequence).

    Raises ValueError with a descriptive message on:
      - empty file
      - unrecognised format
      - invalid/empty sequences
    """
    suffix = path.suffix.lower()

    if suffix in (".fa", ".fasta", ".faa", ".fna"):
        records = parse_fasta(path)
    elif suffix in (".csv", ".tsv"):
        records = _parse_csv(path)
    else:
        # Try FASTA heuristic (presence of '>' at line start), then CSV
        with open(path, "r", encoding="utf-8") as fh:
            first_char = fh.read(1)
        if first_char == ">":
            records = parse_fasta(path)
        else:
            records = _parse_csv(path)

    if not records:
        raise ValueError(
            f"Input file '{path}' contains no sequences. "
            "Check that the file is non-empty and in FASTA or CSV format."
        )

    # Validate sequences — fail fast with a clear message before alignment
    cleaned: list[tuple[str, str]] = []
    errors: list[str] = []
    for name, seq in records:
        if not seq:
            errors.append(f"  '{name}': empty sequence body.")
            continue
        bad = sorted(set(seq) - VALID_AA)
        if bad:
            errors.append(
                f"  '{name}': invalid amino acid characters {bad!r}. "
                "Only standard 20 amino acids (A-Y, no B/J/O/U/X/Z) are accepted."
            )
            continue
        cleaned.append((name, seq))

    if errors:
        raise ValueError(
            f"Sequence validation failed for {len(errors)} sequence(s):\n"
            + "\n".join(errors)
        )

    return cleaned


def _parse_csv(path: Path) -> list[tuple[str, str]]:
    with open(path, newline="", encoding="utf-8") as fh:
        sample = fh.read(4096)
        fh.seek(0)
        delimiter = "\t" if "\t" in sample else ","
        reader = csv.DictReader(fh, delimiter=delimiter)
        if reader.fieldnames is None:
            return []

        fields_lower = [f.strip().lower() for f in reader.fieldnames]

        name_key = _find_col(fields_lower, ("name", "id", "seq_id", "sequence_id"))
        seq_key  = _find_col(fields_lower, ("sequence", "seq", "aa", "protein"))

        if name_key is None or seq_key is None:
            raise ValueError(
                f"CSV '{path}' must contain name and sequence columns. "
                f"Recognised name cols: name/id/seq_id. "
                f"Recognised seq cols: sequence/seq/aa/protein. "
                f"Found: {list(reader.fieldnames)}"
            )

        orig_name = reader.fieldnames[fields_lower.index(name_key)]
        orig_seq  = reader.fieldnames[fields_lower.index(seq_key)]

        return [
            (row[orig_name].strip(), row[orig_seq].strip().upper())
            for row in reader
            if row.get(orig_name, "").strip() and row.get(orig_seq, "").strip()
        ]


def _find_col(fields: list[str], candidates: tuple[str, ...]) -> str | None:
    for c in candidates:
        if c in fields:
            return c
    return None


# ---------------------------------------------------------------------------
# Smith-Waterman alignment via parasail
# ---------------------------------------------------------------------------

def _sw_align(query: str, ref: str) -> tuple[int, int, int, int]:
    """Run SW local alignment (query vs ref) using parasail.sw_trace_striped_16.

    Returns:
        (n_matches, q_aligned, r_aligned, aln_columns)

    where:
        n_matches    — identical residue pairs in the aligned region
        q_aligned    — number of QUERY residues in the aligned region (gaps excluded)
        r_aligned    — number of REF residues in the aligned region (gaps excluded)
        aln_columns  — total alignment columns (including gap columns)
                        → identity = n_matches / aln_columns

    IMPLEMENTATION NOTE — traceback-based counting:
      We use result.traceback.query / result.traceback.ref (the actual aligned
      strings with '-' for gaps) to count residues directly.  This avoids any
      reliance on beg_query / beg_ref / cigar offset attributes, which vary
      between parasail versions and caused inflated coverage numbers in
      parasail 1.3.4 (cigar.beg_query is always 0 relative to the aligned
      region, not an absolute offset into the full query string).

    DESIGN NOTE — identity denominator:
      The brief says "identity = matches / alignment_length" where alignment_length
      is the length of the aligned region *including gaps*.  Using aln_columns (which
      includes I/D gap columns) gives a more conservative (lower) identity than using
      only the residue-to-residue matched columns.  We use aln_columns to match the
      conventional definition.
    """
    result = parasail.sw_trace_striped_16(query, ref, _GAP_OPEN, _GAP_EXTEND, parasail.blosum62)

    # Use the traceback strings directly — most robust across parasail versions.
    tb = result.traceback
    q_al = tb.query   # aligned query string (may contain '-')
    r_al = tb.ref     # aligned ref   string (may contain '-')

    if not q_al or not r_al:
        raise RuntimeError(
            f"parasail traceback returned empty strings for query len={len(query)} "
            f"ref len={len(ref)}. This should not happen with valid non-empty inputs."
        )

    if len(q_al) != len(r_al):
        raise RuntimeError(
            f"parasail traceback strings have different lengths: "
            f"q_al={len(q_al)}, r_al={len(r_al)}. Unexpected parasail output."
        )

    aln_columns = len(q_al)
    n_matches   = sum(1 for a, b in zip(q_al, r_al) if a == b and a != '-')
    q_aligned   = sum(1 for a in q_al if a != '-')   # query residues in aligned region
    r_aligned   = sum(1 for b in r_al if b != '-')   # ref   residues in aligned region

    # Defensive assertions — these should never fire with valid inputs
    if q_aligned > len(query):
        raise RuntimeError(
            f"BUG: q_aligned={q_aligned} > len(query)={len(query)}. "
            "Coverage numerator exceeds sequence length."
        )
    if r_aligned > len(ref):
        raise RuntimeError(
            f"BUG: r_aligned={r_aligned} > len(ref)={len(ref)}. "
            "Coverage numerator exceeds sequence length."
        )

    return n_matches, q_aligned, r_aligned, aln_columns


# ---------------------------------------------------------------------------
# Database loading
# ---------------------------------------------------------------------------

def load_known_binders(target: str, verbose: bool = True) -> list[tuple[str, str]]:
    """Load all known binder sequences for the given target into a list.

    Returns list of (name, sequence).

    NOTE ON RAM: 122,273 sequences × ~60 chars each ≈ 7–8 MB of raw string data.
    Python overhead adds ~300–500 MB for the list of (str, str) tuples.
    With 3 GB available, this is acceptable.
    """
    files = KNOWN_BINDER_FILES.get(target.lower())
    if files is None:
        raise ValueError(f"Unknown target '{target}'. Choose from: {list(KNOWN_BINDER_FILES)}")

    all_records: list[tuple[str, str]] = []
    for filename, label in files:
        fpath = DATA_DIR / filename
        if not fpath.exists():
            print(f"  [WARN] File not found, skipping: {fpath}", file=sys.stderr)
            continue
        records = parse_fasta(fpath)
        if verbose:
            print(f"  Loaded {len(records):>7,} seqs from {filename}  ({label})", file=sys.stderr)
        all_records.extend(records)

    if not all_records:
        raise RuntimeError(
            f"No known binder sequences loaded for target '{target}'. "
            "Check that data/known_binders/ files exist."
        )

    return all_records


# ---------------------------------------------------------------------------
# Core novelty check — single candidate
# ---------------------------------------------------------------------------

def check_one_candidate(
    cand_name: str,
    cand_seq: str,
    known_binders: list[tuple[str, str]],
) -> dict:
    """Run novelty check for a single candidate against all known binders.

    Does ONE pass over known_binders, computing both primary and strict verdicts.

    Returns a dict with all output columns.
    """
    cand_len = len(cand_seq)

    # --- primary tracking: best identity among hits with coverage_ours >= 80 ---
    best_identity        = 0.0
    best_cov_ours        = 0.0
    best_cov_theirs      = 0.0
    best_hit_id          = ""
    n_hits_cov_ours_80   = 0   # number of known binders where coverage_ours >= 80

    # --- strict tracking: best identity among hits with max(cov_ours, cov_theirs) >= 80 ---
    strict_best_identity = 0.0
    strict_best_cov      = 0.0
    strict_best_hit_id   = ""

    for kb_name, kb_seq in known_binders:
        kb_len = len(kb_seq)
        if kb_len == 0:
            continue

        try:
            n_matches, q_aln, r_aln, aln_cols = _sw_align(cand_seq, kb_seq)
        except Exception as exc:
            print(
                f"  [WARN] Alignment failed for '{kb_name}': {exc}",
                file=sys.stderr,
            )
            continue

        if aln_cols == 0:
            continue

        # --- coverage ---
        cov_ours   = (q_aln / cand_len) * 100.0
        cov_theirs = (r_aln / kb_len)   * 100.0
        cov_strict = max(cov_ours, cov_theirs)

        # --- identity (denominator = alignment columns, includes gaps) ---
        identity = (n_matches / aln_cols) * 100.0

        # --- primary verdict tracking ---
        if cov_ours >= COVERAGE_THRESHOLD:
            n_hits_cov_ours_80 += 1
            if identity > best_identity:
                best_identity   = identity
                best_cov_ours   = cov_ours
                best_cov_theirs = cov_theirs
                best_hit_id     = kb_name

        # --- strict verdict tracking ---
        if cov_strict >= COVERAGE_THRESHOLD:
            if identity > strict_best_identity:
                strict_best_identity = identity
                strict_best_cov      = cov_strict
                strict_best_hit_id   = kb_name

    # --- verdicts ---
    # Primary: FAIL if best identity > 70 for any hit with coverage_ours >= 80
    verdict = "FAIL" if best_identity > IDENTITY_THRESHOLD else "PASS"

    # Strict: FAIL if best identity > 70 for any hit with max(cov_ours, cov_theirs) >= 80
    verdict_strict = "FAIL" if strict_best_identity > IDENTITY_THRESHOLD else "PASS"

    return {
        "name":                       cand_name,
        "length":                     cand_len,
        "max_identity_pct":           round(best_identity, 2),
        "worst_hit_id":               best_hit_id,
        "coverage_ours_pct":          round(best_cov_ours, 2),
        "coverage_theirs_pct":        round(best_cov_theirs, 2),
        "verdict":                    verdict,
        "verdict_strict":             verdict_strict,
        "n_hits_over_80pct_coverage": n_hits_cov_ours_80,
    }


# ---------------------------------------------------------------------------
# Batch run with progress display
# ---------------------------------------------------------------------------

FIELDNAMES = [
    "name", "length",
    "max_identity_pct", "worst_hit_id",
    "coverage_ours_pct", "coverage_theirs_pct",
    "verdict", "verdict_strict",
    "n_hits_over_80pct_coverage",
]


def run_novelty_check(
    candidates: list[tuple[str, str]],
    known_binders: list[tuple[str, str]],
    verbose: bool = True,
) -> list[dict]:
    """Run novelty check for all candidates. Returns list of result dicts."""
    results: list[dict] = []
    n = len(candidates)
    t0 = time.perf_counter()

    for i, (name, seq) in enumerate(candidates, 1):
        t_cand = time.perf_counter()
        row = check_one_candidate(name, seq, known_binders)
        elapsed = time.perf_counter() - t_cand

        if verbose:
            verdict_tag = "FAIL" if row["verdict"] == "FAIL" else "PASS"
            print(
                f"  [{i:>3}/{n}] {name[:40]:<40}  "
                f"id={row['max_identity_pct']:5.1f}%  cov={row['coverage_ours_pct']:5.1f}%  "
                f"→ {verdict_tag}  ({elapsed:.1f}s)",
                file=sys.stderr,
            )
        results.append(row)

    total_elapsed = time.perf_counter() - t0
    if verbose:
        n_pass = sum(1 for r in results if r["verdict"] == "PASS")
        n_fail = n - n_pass
        print(
            f"\n  Total: {total_elapsed:.1f}s  |  {n} candidates vs {len(known_binders):,} known binders"
            f"  |  PASS {n_pass}  FAIL {n_fail}",
            file=sys.stderr,
        )

    return results


# ---------------------------------------------------------------------------
# Output formatting
# ---------------------------------------------------------------------------

def print_table(results: list[dict]) -> None:
    """Print results as a formatted ASCII table to stdout."""
    hdr = (
        f"{'#':<4} {'Name':<35} {'Len':>5}  {'MaxID%':>7}  {'CovOurs%':>9}  "
        f"{'CovThem%':>9}  {'Verdict':<8}  {'VrdStrict':<10}  {'NHits>80':>8}"
    )
    print(hdr)
    print("-" * len(hdr))
    for i, r in enumerate(results, 1):
        print(
            f"{i:<4} {r['name'][:34]:<35} {r['length']:>5}  "
            f"{r['max_identity_pct']:>7.1f}  {r['coverage_ours_pct']:>9.1f}  "
            f"{r['coverage_theirs_pct']:>9.1f}  {r['verdict']:<8}  "
            f"{r['verdict_strict']:<10}  {r['n_hits_over_80pct_coverage']:>8}"
        )
    print()
    n_pass = sum(1 for r in results if r["verdict"] == "PASS")
    n_fail = len(results) - n_pass
    n_pass_s = sum(1 for r in results if r["verdict_strict"] == "PASS")
    n_fail_s = len(results) - n_pass_s
    print(f"Summary (primary):  {n_pass} PASS  /  {n_fail} FAIL   (of {len(results)} candidates)")
    print(f"Summary (strict):   {n_pass_s} PASS  /  {n_fail_s} FAIL")


def write_csv(results: list[dict], out_path: Path) -> None:
    with open(out_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=FIELDNAMES, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(results)
    print(f"[INFO] Results written to {out_path}", file=sys.stderr)


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Novelty checker for protein minibinder design candidates (Changjiang Cup 2026).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  .venv/bin/python3 tools/novelty_check.py --input candidates.fasta --target il7ra
  .venv/bin/python3 tools/novelty_check.py --input candidates.csv   --target tl1a --out results.csv
        """,
    )
    p.add_argument("--input",  "-i", required=True,  help="Candidate sequences (.fasta or .csv)")
    p.add_argument(
        "--target", "-t", required=True,
        choices=["il7ra", "tl1a"],
        help="Target protein (determines which known-binder DB to use)",
    )
    p.add_argument("--out",   "-o", default=None,       help="Output CSV path (optional)")
    p.add_argument("--quiet", "-q", action="store_true", help="Suppress per-candidate progress")
    return p


def main() -> None:
    _assert_parasail()
    args = _build_parser().parse_args()

    input_path = Path(args.input)
    if not input_path.exists():
        sys.exit(f"[ERROR] Input file not found: {input_path}")

    # Parse and validate candidates
    try:
        candidates = parse_candidates(input_path)
    except ValueError as e:
        sys.exit(f"[ERROR] Invalid input: {e}")

    print(f"[INFO] Loaded {len(candidates)} candidate(s) from {input_path}", file=sys.stderr)

    # Load known binders
    print(f"[INFO] Loading known binders for target '{args.target}'...", file=sys.stderr)
    try:
        known_binders = load_known_binders(args.target, verbose=True)
    except (ValueError, RuntimeError) as e:
        sys.exit(f"[ERROR] {e}")

    print(f"[INFO] Total known binders: {len(known_binders):,}", file=sys.stderr)
    n_align = len(candidates) * len(known_binders)
    print(f"[INFO] Alignments to run: {n_align:,}", file=sys.stderr)
    print(f"[INFO] Running Smith-Waterman via parasail...\n", file=sys.stderr)

    # Run checks
    results = run_novelty_check(candidates, known_binders, verbose=not args.quiet)

    # Print table
    print_table(results)

    # Write CSV if requested
    if args.out:
        write_csv(results, Path(args.out))


if __name__ == "__main__":
    main()
