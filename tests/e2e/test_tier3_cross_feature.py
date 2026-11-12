"""
tests/e2e/test_tier3_cross_feature.py - Tier 3: Cross-Feature Combinations & Data Flow (C1-C5).
Covers end-to-end component interactions, fallbacks, and multi-stage pipelines (5 tests total).
"""

import os
from pathlib import Path
import pandas as pd
import pytest
import perceval as pcvl

from quantum_folding.lattice import (
    encode_sequence,
    compute_coordinates,
    is_self_avoiding,
    get_contacts,
)
from quantum_folding.hamiltonian import (
    build_folding_hamiltonian,
)
from quantum_folding.circuit import (
    build_photonic_ansatz,
    build_input_state,
    state_to_conformation,
    decode_distribution,
)
from quantum_folding.quandela_backend import (
    QuandelaCloudClient,
)
from quantum_folding.vqe import (
    evaluate_vqe,
)
from quantum_folding.pipeline import (
    QuantumFoldingPipeline,
)
from tools.submit_sjtu import (
    SJTUPlatformClient,
    validate_submission_csv,
)
from tests.e2e.conftest import (
    TL1A_TARGET_ID,
    IL7RA_TARGET_ID,
    VALID_QUANDELA_TOKEN,
)


def test_c1_sequence_to_lattice_to_circuit_to_quandela_execution(mock_cloud_client, sample_peptides):
    """
    C1: Sequence -> Lattice -> Circuit -> Quandela Execution.
    Verifies full quantum computational path from amino acid string to energy expectation.
    """
    seq = sample_peptides["short_mj"]  # "ACDEFY"
    protein = encode_sequence(seq, model="MJ")
    hamiltonian = build_folding_hamiltonian(protein)

    num_turns = len(seq) - 1
    modes = num_turns * 2
    circuit, params = build_photonic_ansatz(num_modes=modes, num_photons=2)
    input_state = build_input_state(num_modes=modes, num_photons=2)

    # Set parameters and submit to Quandela client
    for i, p in enumerate(params):
        p.set_value(float(i) * 0.1)

    sub = mock_cloud_client.submit_job(circuit, input_state, shots=100)
    job_id = sub["job_id"]
    status = mock_cloud_client.get_job_status(job_id)
    assert status == "completed"

    res = mock_cloud_client.get_job_results(job_id)
    distribution = decode_distribution(res["results"], num_turns)
    assert len(distribution) > 0
    assert sum(distribution.values()) == pytest.approx(1.0, abs=1e-5)

    # Compute expected energy over decoded conformations
    expected_energy = sum(hamiltonian.energy(conf) * prob for conf, prob in distribution.items())
    assert isinstance(expected_energy, float)
    assert not pd.isna(expected_energy)


def test_c2_candidate_pool_to_quantum_scoring_to_ranking(mock_cloud_client, temp_output_dir):
    """
    C2: Candidate Pool -> Quantum Scoring -> Ranking & Filtering.
    Verifies processing candidate designs into ranked quantum folding stability output.
    """
    pool_file = Path("notebooks/candidate_pool_for_screening.csv")
    if not pool_file.exists():
        # Fallback to local representative pool
        pool_file = Path(temp_output_dir) / "pool.csv"
        pool_file.write_text(
            "name,sequence,target\n"
            "C1,AALEELLKALEELLKALEELLKALEPWHETVMSINLFMIEFIIAHPSPKEIQITWELKDFLKKNFEPVSKSFHDRLFPFGDYKWRSPLEIGK,TL1A\n"
            "C2,EKLRKLEKQARKQGDEVLVMLARMVLEWSETVMSINLFMIEFIIAHPSPKEIQITWELKDFLKKNFEPVSKSFHDRLFPFGDYKWRSPLEIGK,TL1A\n"
            "C3,AALEELLKALEELLKALEKLAPKRTQHADWTNSQIIQERVHEQDYRSDINIIFNIQYWPGLFGEVAQANVENDTTYHY,IL7RA\n"
        )

    pipeline = QuantumFoldingPipeline(backend=mock_cloud_client)
    out_scores = Path(temp_output_dir) / "quantum_folding_scores.csv"

    # Evaluate mini-slice of pool
    df = pipeline.evaluate_pool(str(pool_file), max_candidates=5, output_csv=str(out_scores))
    assert out_scores.exists()
    assert "candidate_id" in df.columns
    assert "quantum_score" in df.columns
    assert "rank" in df.columns

    # Ranking is strictly descending by quantum_score
    scores = df["quantum_score"].tolist()
    assert all(scores[i] >= scores[i + 1] for i in range(len(scores) - 1))
    assert df["rank"].tolist() == list(range(1, len(df) + 1))


def test_c3_top_ranked_to_csv_validation_to_submission_payload(mock_cloud_client, temp_output_dir):
    """
    C3: Top Ranked Sequences -> CSV Validation & Formatting -> Submission Payload.
    Verifies export of top 20 candidates into fully compliant competition CSV.
    """
    # Create 20 top-scoring candidates
    rows = ["name,sequence"]
    for i in range(20):
        # 85-95 aa minibinder sequence
        seq = "AALEELLK" + ("A" * (77 + (i % 10)))
        rows.append(f"Q_TOP_{i+1:02d},{seq}")

    submission_csv = Path(temp_output_dir) / "tl1a_submission_ready.csv"
    submission_csv.write_text("\n".join(rows))

    # Validate against competition validator
    csv_text = submission_csv.read_text()
    count, errs = validate_submission_csv(csv_text, modality="minibinder")
    assert count == 20
    assert len(errs) == 0, f"Validation errors found: {errs}"

    # Package into client submission structure
    client = SJTUPlatformClient()
    client.login()
    res = client.submit_candidates(
        target_id=TL1A_TARGET_ID,
        csv_path=str(submission_csv),
        action="confirm",
        agreed=True,
        abstract="Quantum-classical hybrid VQE folding pipeline",
    )
    assert res["status"] == "confirmed"
    assert res["count"] == 20
    assert res["target_code"] == "tl1a"


def test_c4_sjtu_login_to_target_query_to_submission_to_verification(tmp_path):
    """
    C4: SJTU Client Login -> Target Query -> Submission -> Batch Verification.
    Verifies end-to-end platform interaction lifecycle.
    """
    client = SJTUPlatformClient()
    # 1. Login
    auth = client.login()
    assert client.access_token is not None
    assert auth["user"]["team"]["name"] == "QAAA Lab"

    # 2. Query targets and verify IL-7Ra track
    targets = client.get_targets()
    il7ra_target = next(t for t in targets if t["id"] == IL7RA_TARGET_ID)
    assert il7ra_target["code"] == "il7ra"
    assert il7ra_target["modality"] == "minibinder"

    # 3. Create valid CSV and submit
    csv_file = tmp_path / "il7ra_batch.csv"
    csv_file.write_text("name,sequence\nIL7RA_C1," + "A" * 80 + "\n")

    sub_res = client.submit_candidates(
        target_id=IL7RA_TARGET_ID,
        csv_path=str(csv_file),
        agreed=True,
        note="Initial IL-7Ra quantum submission",
    )
    batch_id = sub_res["id"]
    assert batch_id.startswith("batch-")

    # 4. Verification of latest submission
    latest = client.get_latest_submission(IL7RA_TARGET_ID)
    assert latest is not None
    assert latest["id"] == batch_id
    assert latest["status"] == "confirmed"
    assert latest["count"] == 1


def test_c5_quantum_cloud_failure_fallback_to_local_simulation_uninterrupted_submission(tmp_path):
    """
    C5: Quantum Cloud Failure -> Fallback to Local Simulation -> Uninterrupted Submission.
    Verifies system resilience: when remote QPU is in maintenance, pipeline falls back
    to local sim:slos and completes platform submission without error.
    """
    # Initialize client targeting maintenance QPU
    client_fallback = QuandelaCloudClient(VALID_QUANDELA_TOKEN, platform="qpu:ascella")
    status = client_fallback.check_platform_status()
    assert status["status"] == "maintenance"

    # Pipeline uses client, which automatically executes via local SLOS
    pipeline = QuantumFoldingPipeline(backend=client_fallback)
    cand_eval = pipeline.evaluate_candidate("AALEELLKALEELLK", candidate_id="cand_resilient")
    assert "quantum_score" in cand_eval
    assert not pd.isna(cand_eval["quantum_score"])

    # Proceed to SJTU submission with evaluated candidates
    csv_file = tmp_path / "resilient_submission.csv"
    csv_file.write_text(f"name,sequence\n{cand_eval['candidate_id']}," + "A" * 70 + "\n")

    sjtu_client = SJTUPlatformClient()
    sjtu_client.login()
    res = sjtu_client.submit_candidates(
        target_id=TL1A_TARGET_ID,
        csv_path=str(csv_file),
        agreed=True,
        note="Resilient submission via fallback simulation",
    )
    assert res["status"] == "confirmed"
    assert res["count"] == 1
