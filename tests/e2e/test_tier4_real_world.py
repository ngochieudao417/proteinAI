"""
tests/e2e/test_tier4_real_world.py - Tier 4: Real-World Application Workload Scenarios (S1-S5).
Covers production workflows for TL1A, IL-7Ra, live/mock Quandela cloud, dual submission,
and adversarial fault injection & recovery (5 tests total).
"""

import os
import re
from pathlib import Path
import pandas as pd
import pytest
import requests
import perceval as pcvl

from quantum_folding.lattice import (
    encode_sequence,
    compute_coordinates,
    is_self_avoiding,
)
from quantum_folding.hamiltonian import (
    build_folding_hamiltonian,
)
from quantum_folding.circuit import (
    build_photonic_ansatz,
    build_input_state,
)
from quantum_folding.quandela_backend import (
    QuandelaCloudClient,
)
from quantum_folding.pipeline import (
    QuantumFoldingPipeline,
)
from tools.submit_sjtu import (
    SJTUPlatformClient,
    validate_sequence,
    validate_submission_csv,
)
from tests.e2e.conftest import (
    TL1A_TARGET_ID,
    IL7RA_TARGET_ID,
    VALID_QUANDELA_TOKEN,
    SJTU_EMAIL,
    SJTU_PASSWORD,
)


def test_s1_full_minibatch_pipeline_tl1a(mock_cloud_client, temp_output_dir):
    """
    S1: Full Mini-Batch Pipeline for TL1A Target.
    Realistic workflow: extracts TL1A candidates, computes quantum stability score,
    exports validated CSV, submits to SJTU platform, and verifies batch record.
    """
    tl1a_candidates = [
        {"name": f"TL1A_CAND_{i:02d}", "sequence": "AALEELLKALEELLKALEELLKALEPWHETVMSINLFMIEFIIAHPSPKEIQITWELKDFLKKNFEPVSKSFHDRLFPFGDYKWRSPLEIGK"}
        for i in range(1, 6)
    ]

    pipeline = QuantumFoldingPipeline(backend=mock_cloud_client)
    scored = pipeline.run_batch(tl1a_candidates)
    assert len(scored) == 5

    out_csv = Path(temp_output_dir) / "tl1a_submission.csv"
    lines = ["name,sequence"]
    for c in scored:
        lines.append(f"{c['name']},{c['sequence']}")
    out_csv.write_text("\n".join(lines))

    # Pre-submission verification
    count, errs = validate_submission_csv(out_csv.read_text(), modality="minibinder")
    assert count == 5
    assert len(errs) == 0

    # Submission to SJTU platform
    client = SJTUPlatformClient()
    client.login()
    sub_res = client.submit_candidates(
        target_id=TL1A_TARGET_ID,
        csv_path=str(out_csv),
        action="confirm",
        agreed=True,
        abstract="TL1A quantum VQE folding stabilization pipeline",
    )
    assert sub_res["status"] == "confirmed"
    assert sub_res["target_code"] == "tl1a"
    assert sub_res["count"] == 5
    assert "id" in sub_res


def test_s2_full_minibatch_pipeline_il7ra(mock_cloud_client, temp_output_dir):
    """
    S2: Full Mini-Batch Pipeline for IL-7Ra Target.
    Evaluates binding motif contacts, quantum stability, generates platform CSV,
    and captures submission confirmation.
    """
    il7ra_candidates = [
        {"name": f"IL7RA_CAND_{i:02d}", "sequence": "AALEELLKALEELLKALEKLAPKRTQHADWTNSQIIQERVHEQDYRSDINIIFNIQYWPGLFGEVAQANVENDTTYHY"}
        for i in range(1, 6)
    ]

    pipeline = QuantumFoldingPipeline(backend=mock_cloud_client)
    scored = pipeline.run_batch(il7ra_candidates)
    assert len(scored) == 5

    out_csv = Path(temp_output_dir) / "il7ra_submission.csv"
    lines = ["name,sequence"]
    for c in scored:
        lines.append(f"{c['name']},{c['sequence']}")
    out_csv.write_text("\n".join(lines))

    # Validate against minibinder track rules
    count, errs = validate_submission_csv(out_csv.read_text(), modality="minibinder")
    assert count == 5
    assert len(errs) == 0

    client = SJTUPlatformClient()
    client.login()
    sub_res = client.submit_candidates(
        target_id=IL7RA_TARGET_ID,
        csv_path=str(out_csv),
        action="confirm",
        agreed=True,
        note="IL-7Ra quantum-optimized candidate batch",
    )
    assert sub_res["status"] == "confirmed"
    assert sub_res["target_code"] == "il7ra"
    assert sub_res["count"] == 5


def test_s3_quandela_cloud_smoke_test(mock_quandela_token):
    """
    S3: Quandela Cloud Authentication & Execution Smoke Test.
    Verifies authentication token formatting, platform availability, and job execution.
    """
    client = QuandelaCloudClient(token=mock_quandela_token, platform="sim:slos")
    assert client.token.startswith("_T_")

    avail = client.check_platform_status()
    assert avail["status"] == "available"

    # Circuit submission smoke test
    c = pcvl.Circuit(2) // pcvl.BS()
    inp = pcvl.BasicState([1, 1])
    sub = client.submit_job(c, inp, shots=50)

    job_id = sub["job_id"]
    status = client.get_job_status(job_id)
    assert status == "completed"

    res = client.get_job_results(job_id)
    assert "results" in res
    assert sum(res["raw_counts"].values()) == 50


def test_s4_multi_target_dual_submission_batch_verification(tmp_path):
    """
    S4: Multi-Target Dual Submission & Batch Verification.
    Simultaneous submission across both TL1A and IL-7Ra competition tracks.
    """
    client = SJTUPlatformClient()
    client.login()

    # TL1A submission
    csv_tl1a = tmp_path / "sub_tl1a.csv"
    csv_tl1a.write_text("name,sequence\nTL1A_TOP_01," + "A" * 70 + "\nTL1A_TOP_02," + "A" * 75 + "\n")
    res_tl1a = client.submit_candidates(target_id=TL1A_TARGET_ID, csv_path=str(csv_tl1a), agreed=True)

    # IL-7Ra submission
    csv_il7ra = tmp_path / "sub_il7ra.csv"
    csv_il7ra.write_text("name,sequence\nIL7RA_TOP_01," + "A" * 80 + "\nIL7RA_TOP_02," + "A" * 85 + "\n")
    res_il7ra = client.submit_candidates(target_id=IL7RA_TARGET_ID, csv_path=str(csv_il7ra), agreed=True)

    # Verify both batch records are distinct and verified
    assert res_tl1a["id"] != res_il7ra["id"]
    assert res_tl1a["target_code"] == "tl1a"
    assert res_il7ra["target_code"] == "il7ra"

    # Query all submissions
    all_subs = client.get_submissions(grouped=True)
    items = all_subs["groups"][0]["items"]
    submitted_targets = {item["target_id"] for item in items}
    assert TL1A_TARGET_ID in submitted_targets
    assert IL7RA_TARGET_ID in submitted_targets


def test_s5_adversarial_fault_injection_and_end_to_end_recovery(tmp_path, mock_cloud_client, monkeypatch):
    """
    S5: Adversarial Fault Injection & End-to-End Recovery.
    Injects invalid amino acids, out-of-bound lengths, duplicates, and transient network glitches.
    Verifies automatic filtering, network retry, and 100% clean submission.
    """
    # 1. Dirty candidate pool with errors
    raw_candidates = [
        {"name": "valid_01", "sequence": "AALEELLKALEELLKALEELLKALEPWHETVMSINLFMIEFIIAHPSPKEIQITWELKDFLKKNFEPVSKSFHDRLFPFGDYKWRSPLEIGK"},
        {"name": "valid_02", "sequence": "AALEELLKALEELLKALEKLAPKRTQHADWTNSQIIQERVHEQDYRSDINIIFNIQYWPGLFGEVAQANVENDTTYHY"},
        {"name": "corrupt_noncanonical", "sequence": "AALEELLKALEELLK" + "X" * 10 + "A" * 30},
        {"name": "corrupt_underflow", "sequence": "A" * 25},  # < 40 aa
        {"name": "corrupt_overflow", "sequence": "A" * 115},  # > 100 aa
        {"name": "valid_01", "sequence": "A" * 60},  # Duplicate name
    ]

    # SUT Candidate Filtering: Pipeline autonomously drops invalid/corrupt rows and deduplicates
    pipeline = QuantumFoldingPipeline(backend=mock_cloud_client)
    scored = pipeline.run_batch(raw_candidates, filter_invalid=True)

    # Assert SUT dropped corrupt candidates and deduplicated
    assert len(scored) == 2, f"Expected 2 valid candidates after SUT filtering, got {len(scored)}"
    assert {c["name"] for c in scored} == {"valid_01", "valid_02"}

    # Export SUT-evaluated candidates to submission CSV
    clean_csv = tmp_path / "recovered_submission.csv"
    lines = ["name,sequence"]
    for c in scored:
        lines.append(f"{c['name']},{c['sequence']}")
    clean_csv.write_text("\n".join(lines))

    # 2. SUT Network Recovery: Platform client autonomously retries on transient 500 error
    client = SJTUPlatformClient()
    client.login()

    attempt_count = 0
    orig_post = client.session.post

    def flaking_post(url, *args, **kwargs):
        nonlocal attempt_count
        attempt_count += 1
        if attempt_count == 1 and "/api/submissions/csv" in url:
            # Simulate transient database timeout on 1st attempt
            resp = requests.Response()
            resp.status_code = 500
            resp._content = b'{"detail": "Transient database timeout"}'
            resp.url = url
            return resp
        return orig_post(url, *args, **kwargs)

    monkeypatch.setattr(client.session, "post", flaking_post)

    # Invoke SUT with retry configuration - NO test-level retry loops
    sub_res = client.submit_candidates(
        target_id=TL1A_TARGET_ID,
        csv_path=str(clean_csv),
        agreed=True,
        max_retries=3,
    )

    # Verification of SUT autonomous retry success
    assert sub_res is not None
    assert sub_res["status"] == "confirmed"
    assert sub_res["count"] == 2
    assert attempt_count == 2, "SUT client must have autonomously retried and succeeded on attempt 2"
