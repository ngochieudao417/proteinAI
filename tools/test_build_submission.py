#!/usr/bin/env python3
"""
test_build_submission.py — Self-test for build_submission.py

Runs all 12 required test cases and prints PASS/FAIL for each.

Usage:
    python3 tools/test_build_submission.py

NOTE: This test file generates synthetic random sequences with a fixed seed.
      It does NOT use real protein sequences from the competition.
      Sequences are labelled "RANDOM_SEED42_..." to make this explicit.
"""

from __future__ import annotations

import io
import os
import random
import re
import sys
import tempfile
from pathlib import Path

# ---------------------------------------------------------------------------
# Make sure we can import build_submission from the tools/ directory
# ---------------------------------------------------------------------------

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / 'tools'))

from build_submission import (
    _validate_our,
    _validate_reference,
    _build_csv_text,
    _check_duplicates,
    _check_tags,
    _load_novelty_report,
    _check_novelty,
    build_submission,
    parse_input,
)

# ---------------------------------------------------------------------------
# Helpers — synthetic sequence generation
# ---------------------------------------------------------------------------

STANDARD_AA = list("ACDEFGHIKLMNPQRSTVWY")


def _rand_seq(length: int, rng: random.Random) -> str:
    """Generate a random amino acid sequence of exactly `length` residues."""
    return ''.join(rng.choices(STANDARD_AA, k=length))


def _make_valid_records(
    n: int,
    min_len: int = 85,
    max_len: int = 100,
    seed: int = 42,
) -> list[tuple[str, str]]:
    """
    Generate `n` valid (name, sequence) records.
    Lengths are drawn uniformly from [min_len, max_len].
    Uses a fixed seed for reproducibility.
    """
    rng = random.Random(seed)
    records = []
    for i in range(n):
        length = rng.randint(min_len, max_len)
        seq = _rand_seq(length, rng)
        name = f"RANDOM_SEED42_{i + 1:03d}_len{length}"
        records.append((name, seq))
    return records


# ---------------------------------------------------------------------------
# Test infrastructure
# ---------------------------------------------------------------------------

_PASS = 0
_FAIL = 0
_results: list[tuple[int, str, str]] = []  # (case_no, label, status)


def _report(case_no: int, label: str, passed: bool, details: str = '') -> None:
    global _PASS, _FAIL
    status = 'PASS' if passed else 'FAIL'
    if passed:
        _PASS += 1
    else:
        _FAIL += 1
    _results.append((case_no, label, status))
    mark = '✓' if passed else '✗'
    print(f"  [{mark}] Case {case_no:>2}: {label}  → {status}")
    if details:
        for line in details.strip().splitlines():
            print(f"           {line}")


def _build_to_tempfile(records, **kwargs) -> tuple[bool, str]:
    """
    Call build_submission() writing to a temp file. Returns (success, out_path_str).
    Captures stdout/stderr to suppress noise during tests unless VERBOSE=1.
    """
    verbose = os.environ.get('VERBOSE', '0') == '1'
    with tempfile.NamedTemporaryFile(suffix='.csv', delete=False) as tmp:
        tmp_path = Path(tmp.name)

    # Redirect stdout if not verbose
    if not verbose:
        old_stdout = sys.stdout
        sys.stdout = io.StringIO()

    try:
        success = build_submission(
            target='il7ra',
            records=records,
            out_path=tmp_path,
            **kwargs,
        )
    finally:
        if not verbose:
            sys.stdout = old_stdout

    return success, str(tmp_path)


def _write_csv_to_temp(csv_text: str, encoding='utf-8', bom=False) -> str:
    """Write raw CSV text to a temp file and return the path."""
    with tempfile.NamedTemporaryFile(
        suffix='.csv', delete=False, mode='wb'
    ) as tmp:
        data = csv_text.encode(encoding)
        if bom:
            data = b'\xef\xbb\xbf' + data
        tmp.write(data)
        return tmp.name


def _write_novelty_csv(entries: list[tuple[str, str]]) -> str:
    """Write a minimal novelty_check.py style CSV and return path."""
    lines = ['name,verdict']
    for name, verdict in entries:
        lines.append(f'{name},{verdict}')
    csv_text = '\n'.join(lines) + '\n'
    with tempfile.NamedTemporaryFile(
        suffix='.csv', delete=False, mode='w', encoding='utf-8'
    ) as tmp:
        tmp.write(csv_text)
        return tmp.name


# ---------------------------------------------------------------------------
# Individual test cases
# ---------------------------------------------------------------------------

def run_case_01():
    """20 valid sequences, 85–100 aa → build succeeds, reference validator passes."""
    records = _make_valid_records(20)
    assert len(records) == 20

    # Our validator: should have zero errors
    our_errs = _validate_our(records)

    # Build CSV text and cross-check with reference validator
    csv_text = _build_csv_text(records)
    _, ref_errs = _validate_reference(csv_text)

    # Also try building to file
    success, tmp_path = _build_to_tempfile(records)
    try:
        os.unlink(tmp_path)
    except Exception:
        pass

    passed = (len(our_errs) == 0 and len(ref_errs) == 0 and success)
    details = ''
    if not passed:
        details = f'our_errs={our_errs}, ref_errs={ref_errs}, success={success}'
    _report(1, '20 valid sequences (85–100 aa) → build OK', passed, details)


def run_case_02():
    """21 sequences → refused, over-limit reported."""
    records = _make_valid_records(21)
    our_errs = _validate_our(records)
    # Must have exactly one error about exceeding the limit
    has_limit_err = any('20' in reason or 'Too many' in reason for _, _, reason in our_errs)
    success, tmp_path = _build_to_tempfile(records)
    try:
        os.unlink(tmp_path)
    except Exception:
        pass
    passed = (len(our_errs) > 0 and has_limit_err and not success)
    details = '' if passed else f'our_errs={our_errs}, success={success}'
    _report(2, '21 sequences → refused (over limit)', passed, details)


def run_case_03():
    """One sequence 39 aa, one 101 aa → both rejected (out of [40,100])."""
    rng = random.Random(42)
    records = [
        ('short_39', _rand_seq(39, rng)),
        ('long_101', _rand_seq(101, rng)),
    ]
    our_errs = _validate_our(records)
    # Both names should appear in errors
    err_names = {name for _, name, _ in our_errs}
    passed = ('short_39' in err_names and 'long_101' in err_names)
    details = '' if passed else f'our_errs={our_errs}'
    _report(3, '39 aa and 101 aa sequences → both rejected', passed, details)


def run_case_04():
    """Sequence of exactly 40 aa and exactly 100 aa → both PASS (bounds inclusive)."""
    rng = random.Random(99)
    records = [
        ('exactly_40', _rand_seq(40, rng)),
        ('exactly_100', _rand_seq(100, rng)),
    ]
    our_errs = _validate_our(records)
    csv_text = _build_csv_text(records)
    _, ref_errs = _validate_reference(csv_text)
    success, tmp_path = _build_to_tempfile(records)
    try:
        os.unlink(tmp_path)
    except Exception:
        pass
    passed = (len(our_errs) == 0 and len(ref_errs) == 0 and success)
    details = '' if passed else f'our_errs={our_errs}, ref_errs={ref_errs}, success={success}'
    _report(4, 'Exactly 40 aa and 100 aa → both PASS (inclusive bounds)', passed, details)


def run_case_05():
    """Sequence containing 'X' → rejected, invalid character reported."""
    rng = random.Random(7)
    seq_with_X = _rand_seq(85, rng)[:42] + 'X' + _rand_seq(85, rng)[43:]
    records = [('has_X', seq_with_X)]
    our_errs = _validate_our(records)
    has_invalid_err = any('Invalid amino acid' in reason or 'invalid' in reason.lower()
                          for _, _, reason in our_errs)
    passed = (len(our_errs) > 0 and has_invalid_err)
    details = '' if passed else f'our_errs={our_errs}'
    _report(5, "Sequence with 'X' → rejected (invalid character)", passed, details)


def run_case_06():
    """name=65 chars rejected; name=0 chars rejected; name=64 chars passes."""
    rng = random.Random(13)
    seq = _rand_seq(85, rng)
    name_65 = 'A' * 65
    name_0  = ''
    name_64 = 'B' * 64

    errs_65  = _validate_our([(name_65,  seq)])
    errs_0   = _validate_our([(name_0,   seq)])
    errs_64  = _validate_our([(name_64,  seq)])

    # 65-char name: must have an error
    has_err_65 = len(errs_65) > 0 and any('64' in r or 'long' in r.lower() for _, _, r in errs_65)
    # empty name: must have an error
    has_err_0 = len(errs_0) > 0 and any('empty' in r.lower() or '1' in r for _, _, r in errs_0)
    # 64-char name: must have NO errors
    has_no_err_64 = (len(errs_64) == 0)

    passed = has_err_65 and has_err_0 and has_no_err_64
    details = '' if passed else (
        f'errs_65={errs_65}, errs_0={errs_0}, errs_64={errs_64}'
    )
    _report(6, 'name 65 chars → rejected; empty → rejected; 64 chars → pass', passed, details)


def run_case_07():
    """Headers: 'NAME,SEQUENCE', 'sequence,name', extra column → all 3 pass."""
    rng = random.Random(21)
    seq = _rand_seq(85, rng)
    name = 'test_header_flexibility'

    # Build a CSV with the given header and validate via reference validator
    def check_header_csv(header_line: str, row_line: str) -> bool:
        csv_text = f"{header_line}\r\n{row_line}\r\n"
        _, errs = _validate_reference(csv_text)
        return len(errs) == 0

    # Case A: uppercase header
    row = f'{name},{seq}'
    ok_a = check_header_csv('NAME,SEQUENCE', row)

    # Case B: reversed column order
    row_b = f'{seq},{name}'
    ok_b = check_header_csv('sequence,name', row_b)

    # Case C: extra column 'note'
    row_c = f'{name},{seq},some_note'
    ok_c = check_header_csv('name,sequence,note', row_c)

    passed = ok_a and ok_b and ok_c
    details = '' if passed else f'ok_a={ok_a}, ok_b={ok_b}, ok_c={ok_c}'
    _report(7, 'Flexible headers (uppercase / reversed / extra col) → all pass', passed, details)


def run_case_08():
    """CSV with UTF-8 BOM → passes (BOM stripped)."""
    rng = random.Random(31)
    seq = _rand_seq(85, rng)
    csv_text_no_bom = f"name,sequence\r\nbom_test,{seq}\r\n"

    # Validate text with BOM prepended
    csv_text_with_bom = '\ufeff' + csv_text_no_bom
    _, errs = _validate_reference(csv_text_with_bom)

    passed = len(errs) == 0
    details = '' if passed else f'ref_errs={errs}'
    _report(8, 'CSV with UTF-8 BOM → passes (BOM stripped)', passed, details)


def run_case_09():
    """Sequence with internal spaces and lowercase → passes after strip+uppercase."""
    # Construct a sequence string with spaces and lowercase
    rng = random.Random(37)
    raw_seq_upper = _rand_seq(85, rng)
    # Introduce lowercase and spaces
    half = len(raw_seq_upper) // 2
    raw_seq_mixed = raw_seq_upper[:half].lower() + ' ' + raw_seq_upper[half:]

    records = [('mixed_case_spaces', raw_seq_mixed)]
    our_errs = _validate_our(records)

    # Also check via reference validator (write CSV manually with the raw sequence)
    csv_text = f"name,sequence\r\nmixed_case_spaces,{raw_seq_mixed}\r\n"
    _, ref_errs = _validate_reference(csv_text)

    passed = (len(our_errs) == 0 and len(ref_errs) == 0)
    details = '' if passed else f'our_errs={our_errs}, ref_errs={ref_errs}'
    _report(9, 'Sequence with spaces and lowercase → passes after normalisation', passed, details)


def run_case_10():
    """Two identical sequences → build succeeds but warns about duplicates."""
    rng = random.Random(43)
    seq = _rand_seq(85, rng)
    records = [
        ('dup_A', seq),
        ('dup_B', seq),
    ]
    dup_warnings = _check_duplicates(records)
    # Validation should pass (no hard error for duplicates — only a warning)
    our_errs = _validate_our(records)
    success, tmp_path = _build_to_tempfile(records)
    try:
        os.unlink(tmp_path)
    except Exception:
        pass

    has_dup_warning = len(dup_warnings) > 0
    passed = (has_dup_warning and len(our_errs) == 0 and success)
    details = '' if passed else f'dup_warnings={dup_warnings}, our_errs={our_errs}, success={success}'
    _report(10, 'Two identical sequences → build OK with duplicate warning', passed, details)


def run_case_11():
    """Sequence ending in HHHHHH → build succeeds but warns about His-tag."""
    rng = random.Random(55)
    # Build a 79-aa random sequence + 6 His = 85 aa total
    seq = _rand_seq(79, rng) + 'HHHHHH'
    assert len(seq) == 85
    records = [('his_tag_test', seq)]
    tag_warnings = _check_tags('his_tag_test', seq)
    our_errs = _validate_our(records)
    success, tmp_path = _build_to_tempfile(records)
    try:
        os.unlink(tmp_path)
    except Exception:
        pass

    has_tag_warning = len(tag_warnings) > 0
    passed = (has_tag_warning and len(our_errs) == 0 and success)
    details = '' if passed else f'tag_warnings={tag_warnings}, our_errs={our_errs}, success={success}'
    _report(11, 'Sequence ending with HHHHHH → build OK with tag warning', passed, details)


def run_case_12():
    """--novelty-report with one FAIL → refused; with --force → succeeds."""
    rng = random.Random(61)
    records = [
        ('novelty_ok',   _rand_seq(85, rng)),
        ('novelty_fail', _rand_seq(87, rng)),
    ]

    # Write a fake novelty report
    novelty_path = _write_novelty_csv([
        ('novelty_ok',   'PASS'),
        ('novelty_fail', 'FAIL'),
    ])

    try:
        novelty_map = _load_novelty_report(Path(novelty_path))

        # Without force: should be refused
        fails_without_force = _check_novelty(records, novelty_map)
        success_without, tmp1 = _build_to_tempfile(records, novelty_map=novelty_map, force=False)
        try:
            os.unlink(tmp1)
        except Exception:
            pass

        # With force: should succeed
        success_with_force, tmp2 = _build_to_tempfile(records, novelty_map=novelty_map, force=True)
        try:
            os.unlink(tmp2)
        except Exception:
            pass

    finally:
        try:
            os.unlink(novelty_path)
        except Exception:
            pass

    has_fail = len(fails_without_force) > 0 and 'novelty_fail' in fails_without_force
    passed = (has_fail and not success_without and success_with_force)
    details = '' if passed else (
        f'fails={fails_without_force}, '
        f'success_without={success_without}, success_with_force={success_with_force}'
    )
    _report(12, '--novelty-report FAIL → refused; --force → succeeds', passed, details)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> None:
    print()
    print("=" * 65)
    print("  build_submission.py — Self-test  (12 required cases)")
    print("  Sequences: random, seed=42, labelled RANDOM_SEED42_...")
    print("=" * 65)
    print()

    run_case_01()
    run_case_02()
    run_case_03()
    run_case_04()
    run_case_05()
    run_case_06()
    run_case_07()
    run_case_08()
    run_case_09()
    run_case_10()
    run_case_11()
    run_case_12()

    print()
    print("=" * 65)
    total = _PASS + _FAIL
    print(f"  Results: {_PASS}/{total} PASS   {_FAIL}/{total} FAIL")
    if _FAIL == 0:
        print("  ✓  All cases passed.")
    else:
        print("  ✗  Some cases FAILED — see details above.")
        for no, label, status in _results:
            if status == 'FAIL':
                print(f"      Case {no}: {label}")
    print("=" * 65)
    print()

    sys.exit(0 if _FAIL == 0 else 1)


if __name__ == '__main__':
    main()
