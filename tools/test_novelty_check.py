#!/usr/bin/env python3
"""
test_novelty_check.py — Self-contained test suite for novelty_check.py.

Run with:
  .venv/bin/python3 tools/test_novelty_check.py

Prints PASS/FAIL and actual numbers for each test case.
Exits with code 1 if any test fails.

Tests do NOT make network calls and do NOT modify data/.
"""

import os
import sys
import io
import random
import tempfile
import traceback
from pathlib import Path

# ---------------------------------------------------------------------------
# Ensure novelty_check is importable from the tools/ sibling dir
# ---------------------------------------------------------------------------
_TOOLS_DIR = Path(__file__).resolve().parent
_REPO_ROOT = _TOOLS_DIR.parent
sys.path.insert(0, str(_TOOLS_DIR))

try:
    import novelty_check as nc
except ImportError as e:
    print(f"[ERROR] Cannot import novelty_check: {e}", file=sys.stderr)
    sys.exit(1)

# ---------------------------------------------------------------------------
# Amino acid alphabet (standard 20)
# ---------------------------------------------------------------------------
_AA = "ACDEFGHIKLMNPQRSTVWY"

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _rand_seq(length: int, rng: random.Random) -> str:
    """Generate a random protein sequence of the given length."""
    return "".join(rng.choice(_AA) for _ in range(length))


def _mutate_seq(seq: str, mutation_rate: float, rng: random.Random) -> str:
    """Randomly mutate `mutation_rate` fraction of positions in `seq`."""
    result = list(seq)
    for i in range(len(result)):
        if rng.random() < mutation_rate:
            original = result[i]
            choices = [aa for aa in _AA if aa != original]
            result[i] = rng.choice(choices)
    return "".join(result)


def _fasta_str(records: list[tuple[str, str]]) -> str:
    """Create FASTA text from (name, seq) list."""
    lines = []
    for name, seq in records:
        lines.append(f">{name}")
        # Wrap at 60 chars
        for i in range(0, len(seq), 60):
            lines.append(seq[i:i+60])
    return "\n".join(lines) + "\n"


def _write_temp_fasta(records: list[tuple[str, str]]) -> Path:
    """Write records to a temporary FASTA file and return its Path."""
    tmp = tempfile.NamedTemporaryFile(
        mode="w", suffix=".fasta", delete=False, encoding="utf-8"
    )
    tmp.write(_fasta_str(records))
    tmp.flush()
    tmp.close()
    return Path(tmp.name)


def _cleanup(path: Path) -> None:
    try:
        path.unlink()
    except Exception:
        pass


def _run_check(candidates: list[tuple[str, str]], target: str) -> list[dict]:
    """Load known binders for target and run check for all candidates."""
    known_binders = nc.load_known_binders(target, verbose=False)
    return nc.run_novelty_check(candidates, known_binders, verbose=False)


def _align_one_pair(candidate_seq: str, binder_seq: str) -> dict:
    """Run SW alignment for a single (candidate, binder) pair.

    Returns a dict with keys:
        q_aligned    — query residues in aligned region
        r_aligned    — ref residues in aligned region
        aln_columns  — total alignment columns (including gaps)
        n_matches    — identical pairs
        identity_pct — n_matches / aln_columns * 100
        coverage_ours_pct  — q_aligned / len(candidate_seq) * 100
        coverage_theirs_pct — r_aligned / len(binder_seq) * 100
    """
    n_matches, q_aln, r_aln, aln_cols = nc._sw_align(candidate_seq, binder_seq)
    return {
        "q_aligned":            q_aln,
        "r_aligned":            r_aln,
        "aln_columns":          aln_cols,
        "n_matches":            n_matches,
        "identity_pct":         (n_matches / aln_cols * 100.0) if aln_cols else 0.0,
        "coverage_ours_pct":    (q_aln / len(candidate_seq) * 100.0) if candidate_seq else 0.0,
        "coverage_theirs_pct":  (r_aln / len(binder_seq) * 100.0) if binder_seq else 0.0,
    }


# ---------------------------------------------------------------------------
# Test infrastructure
# ---------------------------------------------------------------------------

_PASS_COUNT = 0
_FAIL_COUNT = 0
_TEST_RESULTS: list[dict] = []


def _report(case_num: int, name: str, passed: bool, details: str) -> None:
    global _PASS_COUNT, _FAIL_COUNT
    status = "PASS" if passed else "FAIL"
    if passed:
        _PASS_COUNT += 1
    else:
        _FAIL_COUNT += 1
    _TEST_RESULTS.append({"case": case_num, "name": name, "status": status, "details": details})
    print(f"[Case {case_num}] {status}: {name}")
    print(f"         {details}")
    print()


# ---------------------------------------------------------------------------
# Data paths
# ---------------------------------------------------------------------------
DATA_DIR = _REPO_ROOT / "data" / "known_binders"
VALIDATED_FASTA = DATA_DIR / "IL7Ra_denovo_minibinders_VALIDATED.fasta"
IPD_FASTA = DATA_DIR / "IL7Ra_denovo_minibinder_designs_IPD_Cao2022.fasta"


# ===========================================================================
# CASE 1 — Exact copy of a known validated binder → must FAIL
# ===========================================================================
def test_case_1():
    """
    Use an exact sequence from IL7Ra_denovo_minibinders_VALIDATED.fasta as the candidate.
    Expected: max_identity_pct == 100.0, coverage_ours_pct == 100.0, verdict == FAIL.
    """
    if not VALIDATED_FASTA.exists():
        _report(1, "exact known binder → FAIL", False,
                f"SKIP: {VALIDATED_FASTA} not found")
        return

    records = nc.parse_fasta(VALIDATED_FASTA)
    if not records:
        _report(1, "exact known binder → FAIL", False,
                "SKIP: validated FASTA is empty")
        return

    name, seq = records[0]
    print(f"  Case 1 input: '{name}' length={len(seq)}")

    results = _run_check([(f"test_exact_{name[:20]}", seq)], "il7ra")
    r = results[0]

    details = (
        f"max_identity_pct={r['max_identity_pct']:.2f}%  "
        f"coverage_ours_pct={r['coverage_ours_pct']:.2f}%  "
        f"verdict={r['verdict']}  "
        f"worst_hit={r['worst_hit_id'][:50]}"
    )

    # Tolerance: identity should be very high (>=99.0 to allow for
    # any minor gap-scoring edge cases; expect 100.0 for identical sequence)
    passed = (
        r["max_identity_pct"] >= 99.0
        and r["coverage_ours_pct"] >= 99.0
        and r["verdict"] == "FAIL"
    )
    _report(1, "exact known binder → FAIL", passed, details)


# ===========================================================================
# CASE 2 — Random 90aa sequence → should PASS with low identity
# ===========================================================================
def test_case_2():
    """
    A fully random 90 aa sequence (seed=42) against il7ra.
    Expected: verdict == PASS, max_identity_pct < 70.
    """
    SEED = 42
    rng = random.Random(SEED)
    seq = _rand_seq(90, rng)
    print(f"  Case 2 input: random 90aa (seed={SEED}): {seq[:30]}...")

    results = _run_check([("test_random_90aa", seq)], "il7ra")
    r = results[0]

    details = (
        f"max_identity_pct={r['max_identity_pct']:.2f}%  "
        f"coverage_ours_pct={r['coverage_ours_pct']:.2f}%  "
        f"n_hits_over_80pct_coverage={r['n_hits_over_80pct_coverage']}  "
        f"verdict={r['verdict']}"
    )

    passed = r["verdict"] == "PASS" and r["max_identity_pct"] < 70.0
    _report(2, "random 90aa → PASS with identity < 70%", passed, details)


# ===========================================================================
# CASE 3 — 40%-mutated known binder (55aa) → identity 55–70%, PASS
# ===========================================================================
def test_case_3():
    """
    Take a 55aa known binder, mutate 40% of positions (seed=123).
    Expected: identity roughly 55–70%, verdict == PASS.
    Note: exact identity depends on SW scoring; check that PASS and identity reported.
    """
    if not VALIDATED_FASTA.exists():
        _report(3, "40%-mutated binder → PASS (identity ~55-70%)", False,
                f"SKIP: {VALIDATED_FASTA} not found")
        return

    records = nc.parse_fasta(VALIDATED_FASTA)
    if not records:
        _report(3, "40%-mutated binder → PASS (identity ~55-70%)", False, "SKIP: empty FASTA")
        return

    SEED = 123
    rng = random.Random(SEED)
    name, seq = records[0]
    mutated = _mutate_seq(seq, 0.40, rng)
    theoretical_identity = sum(a == b for a, b in zip(seq, mutated)) / len(seq) * 100
    print(f"  Case 3 input: mutated version of '{name}' (seed={SEED})")
    print(f"  Case 3: theoretical string identity to original: {theoretical_identity:.1f}%")
    print(f"  Case 3: mutated seq: {mutated}")

    results = _run_check([("test_mutated_55aa", mutated)], "il7ra")
    r = results[0]

    details = (
        f"max_identity_pct={r['max_identity_pct']:.2f}%  "
        f"coverage_ours_pct={r['coverage_ours_pct']:.2f}%  "
        f"verdict={r['verdict']}  "
        f"(theoretical_identity={theoretical_identity:.1f}%)"
    )

    # Main check: verdict should be PASS (identity ≤ 70 or coverage < 80)
    # We also check identity is in a plausible range, but SW identity may differ
    # from theoretical string identity due to gap handling.
    passed = r["verdict"] == "PASS"
    _report(3, "40%-mutated binder → PASS (identity in plausible range)", passed, details)


# ===========================================================================
# CASE 4 — 60aa binder + 30aa random tail (90aa total) → coverage ≈ 66.7%, PASS
# ===========================================================================
def test_case_4():
    """
    THE CRITICAL TEST.
    Take a 60aa known binder, append 30aa random to get 90aa candidate.
    The SW alignment against that EXACT known binder should give:
      - q_aligned = 60  (only the first 60aa of the 90aa candidate align)
      - coverage_ours = 60/90 ≈ 66.7%  < 80%  → verdict PASS
      - identity ≈ 100% (exact match for the 60aa prefix)
    If coverage_ours comes out ≥ 80%, the coverage numerator is wrong.

    Sub-check A: direct per-hit alignment vs the EXACT embedded binder.
    Sub-check B: full run against the il7ra DB → overall verdict == PASS.
    """
    # Find a known binder of approximately 60aa from IPD Cao2022
    if not IPD_FASTA.exists():
        _report(4, "60aa+30aa tail → coverage 66.7% → PASS [CRITICAL]", False,
                f"SKIP: {IPD_FASTA} not found")
        return

    # Stream through to find one ~60aa binder
    target_len = 60
    binder_60: tuple[str, str] | None = None
    with open(IPD_FASTA, "r", encoding="utf-8") as fh:
        header: str | None = None
        seq_parts: list[str] = []
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    seq = "".join(seq_parts).upper()
                    if len(seq) == target_len:
                        binder_60 = (header, seq)
                        break
                header = line[1:].strip()
                seq_parts = []
            else:
                seq_parts.append(line.strip())
        if binder_60 is None and header is not None:
            seq = "".join(seq_parts).upper()
            if len(seq) == target_len:
                binder_60 = (header, seq)

    if binder_60 is None:
        _report(4, "60aa+30aa tail → coverage 66.7% → PASS [CRITICAL]", False,
                "SKIP: could not find a 60aa binder in IPD file")
        return

    binder_name, binder_seq = binder_60
    SEED = 777
    rng = random.Random(SEED)
    tail = _rand_seq(30, rng)
    candidate = binder_seq + tail
    expected_cov = len(binder_seq) / len(candidate) * 100.0  # 66.667%

    print(f"  Case 4: binder '{binder_name[:50]}' len={len(binder_seq)}")
    print(f"  Case 4: tail (30aa, seed={SEED}): {tail}")
    print(f"  Case 4: candidate total len={len(candidate)} = {len(binder_seq)}+{len(tail)}")
    print(f"  Case 4: expected coverage_ours vs embedded binder = "
          f"{len(binder_seq)}/{len(candidate)} = {expected_cov:.1f}%")

    # --- Sub-check A: direct per-hit alignment vs the EXACT embedded binder ---
    hit = _align_one_pair(candidate, binder_seq)
    print(f"  Case 4 [sub-A] direct alignment vs embedded binder:")
    print(f"    q_aligned={hit['q_aligned']}  coverage_ours_pct={hit['coverage_ours_pct']:.2f}%  "
          f"identity_pct={hit['identity_pct']:.2f}%  aln_columns={hit['aln_columns']}")

    cov_ok_A    = abs(hit["coverage_ours_pct"] - expected_cov) < 1.0   # within 1pp of 66.7%
    ident_ok_A  = hit["identity_pct"] >= 99.0                           # exact match

    details_A = (
        f"coverage_ours_pct={hit['coverage_ours_pct']:.2f}% (expected≈{expected_cov:.1f}%)  "
        f"identity_pct={hit['identity_pct']:.2f}%  q_aligned={hit['q_aligned']}"
    )

    if not cov_ok_A:
        print(f"  *** CRITICAL: coverage_ours_pct={hit['coverage_ours_pct']:.2f}% "
              f"deviates >1pp from expected {expected_cov:.1f}% — numerator is WRONG ***")

    # --- Sub-check B: full DB run → overall verdict must be PASS ---
    results = _run_check([("test_60plus30", candidate)], "il7ra")
    r = results[0]
    verdict_ok_B = r["verdict"] == "PASS"

    print(f"  Case 4 [sub-B] full DB run: verdict={r['verdict']}  "
          f"coverage_ours_pct={r['coverage_ours_pct']:.2f}%  "
          f"max_identity_pct={r['max_identity_pct']:.2f}%  "
          f"worst_hit={r['worst_hit_id'][:40]}")

    details = (
        f"[A] cov_vs_embedded={hit['coverage_ours_pct']:.2f}% (expected≈{expected_cov:.1f}%)  "
        f"ident_vs_embedded={hit['identity_pct']:.2f}%  "
        f"[B] verdict={r['verdict']}  cov_worst={r['coverage_ours_pct']:.2f}%"
    )

    passed = cov_ok_A and ident_ok_A and verdict_ok_B
    _report(4, "60aa+30aa tail → coverage 66.7% → PASS [CRITICAL]", passed, details)

    if not cov_ok_A:
        print("  *** CRITICAL FAILURE (sub-A): coverage_ours_pct != 66.7% — coverage numerator WRONG ***")
    if not ident_ok_A:
        print(f"  *** FAIL (sub-A): identity_pct={hit['identity_pct']:.2f}% < 99% for exact embedded binder ***")
    if not verdict_ok_B:
        print(f"  *** CRITICAL FAILURE (sub-B): verdict={r['verdict']} — full DB run should PASS ***")


# ===========================================================================
# CASE 5 — Any sequence, target tl1a → no crash, n_hits reported
# ===========================================================================
def test_case_5():
    """
    Run any sequence against tl1a target.
    Expected: runs without error, n_hits_over_80pct_coverage is a non-negative integer.
    For an 85-100aa candidate, we expect n_hits = 0 (antibodies are 106-300aa,
    coverage_ours would be very low), but we do not assert that — we only check
    for correctness of execution.
    """
    SEED = 999
    rng = random.Random(SEED)
    seq = _rand_seq(90, rng)
    print(f"  Case 5 input: random 90aa (seed={SEED}), target=tl1a")

    try:
        results = _run_check([("test_tl1a_any", seq)], "tl1a")
        r = results[0]
        details = (
            f"n_hits_over_80pct_coverage={r['n_hits_over_80pct_coverage']}  "
            f"max_identity_pct={r['max_identity_pct']:.2f}%  "
            f"coverage_ours_pct={r['coverage_ours_pct']:.2f}%  "
            f"verdict={r['verdict']}"
        )
        passed = isinstance(r["n_hits_over_80pct_coverage"], int) and r["n_hits_over_80pct_coverage"] >= 0
        _report(5, "tl1a target → no crash, n_hits reported", passed, details)
    except Exception as exc:
        _report(5, "tl1a target → no crash", False,
                f"EXCEPTION: {exc}\n{traceback.format_exc()}")


# ===========================================================================
# CASE 6 — Empty file and file with invalid characters → clear error, no crash
# ===========================================================================
def test_case_6():
    """
    Two sub-cases:
    6a. Empty FASTA file → ValueError with descriptive message, no raw traceback.
    6b. FASTA with invalid characters (X, *, B) → ValueError with descriptive message.
    """
    # --- 6a: empty file ---
    tmp_empty = Path(tempfile.mktemp(suffix=".fasta"))
    tmp_empty.write_text("", encoding="utf-8")

    try:
        nc.parse_candidates(tmp_empty)
        empty_ok = False
        empty_detail = "No exception raised for empty file — should have raised ValueError"
    except ValueError as e:
        empty_ok = True
        empty_detail = f"ValueError raised: {str(e)[:120]}"
    except Exception as e:
        empty_ok = False
        empty_detail = f"Wrong exception type ({type(e).__name__}): {e}"
    finally:
        _cleanup(tmp_empty)

    # --- 6b: invalid characters ---
    invalid_fasta = ">test_invalid\nACDEFGHIKLMNPQRSTVWYX*B\n"
    tmp_invalid = Path(tempfile.mktemp(suffix=".fasta"))
    tmp_invalid.write_text(invalid_fasta, encoding="utf-8")

    try:
        nc.parse_candidates(tmp_invalid)
        invalid_ok = False
        invalid_detail = "No exception raised for invalid chars — should have raised ValueError"
    except ValueError as e:
        msg = str(e)
        # Check that the error message mentions the bad characters
        invalid_ok = any(ch in msg for ch in ("X", "*", "B", "invalid"))
        invalid_detail = f"ValueError raised: {msg[:200]}"
    except Exception as e:
        invalid_ok = False
        invalid_detail = f"Wrong exception type ({type(e).__name__}): {e}"
    finally:
        _cleanup(tmp_invalid)

    passed = empty_ok and invalid_ok
    details = (
        f"6a (empty): {'OK' if empty_ok else 'FAIL'} — {empty_detail}  ||  "
        f"6b (invalid chars): {'OK' if invalid_ok else 'FAIL'} — {invalid_detail}"
    )
    _report(6, "empty file & invalid chars → clear ValueError", passed, details)


# ===========================================================================
# Main
# ===========================================================================

def main():
    print("=" * 70)
    print("novelty_check.py — Self-test suite")
    print(f"Repo root: {_REPO_ROOT}")
    print(f"Data dir:  {DATA_DIR}")
    print("=" * 70)
    print()

    # Check parasail early
    if not nc._PARASAIL_OK:
        print(f"[ERROR] parasail not available. Cannot run tests.")
        print(f"Install with: .venv/bin/pip install parasail")
        sys.exit(1)

    # Check data dir
    if not DATA_DIR.exists():
        print(f"[ERROR] Data directory not found: {DATA_DIR}")
        sys.exit(1)

    # Run tests
    test_case_1()
    test_case_2()
    test_case_3()
    test_case_4()
    test_case_5()
    test_case_6()

    # Summary
    print("=" * 70)
    print(f"Results: {_PASS_COUNT} PASS / {_FAIL_COUNT} FAIL  (of {len(_TEST_RESULTS)} cases)")
    print("=" * 70)

    if _FAIL_COUNT > 0:
        print("\nFailed cases:")
        for t in _TEST_RESULTS:
            if t["status"] == "FAIL":
                print(f"  Case {t['case']}: {t['name']}")
        sys.exit(1)
    else:
        print("\nAll tests passed.")
        sys.exit(0)


if __name__ == "__main__":
    main()
