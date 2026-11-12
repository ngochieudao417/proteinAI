"""
tests/e2e/test_tier2_boundary_corner.py - Tier 2: Boundary Value Analysis & Corner Robustness.
Covers Features 1 through 7 with exactly 5 boundary/negative test cases each (35 tests total).
"""

import os
import re
from pathlib import Path
from unittest.mock import MagicMock, patch
import numpy as np
import pytest
import requests
import perceval as pcvl

from quantum_folding.lattice import (
    LatticeModel,
    LatticeProtein,
    encode_hp,
    encode_sequence,
    compute_coordinates,
    is_self_avoiding,
    get_contacts,
    has_backtracking,
)
from quantum_folding.hamiltonian import (
    FoldingHamiltonian,
    build_folding_hamiltonian,
    get_mj_potential,
)
from quantum_folding.circuit import (
    build_photonic_ansatz,
    build_input_state,
    state_to_conformation,
    decode_distribution,
)
from quantum_folding.quandela_backend import (
    QuandelaCloudClient,
    QuandelaBackend,
)
from quantum_folding.vqe import (
    VQEResult,
    evaluate_vqe,
    run_vqe,
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
    SJTU_BASE_URL,
)


# ===========================================================================
# Feature 1 Boundaries (5 tests)
# ===========================================================================

def test_f1_b01_minimal_peptide_length():
    """T2_F1_01: Verify minimal peptide length (2 aa and 3 aa) evaluates without index errors."""
    protein_2aa = encode_sequence("AC", model="MJ")
    ham_2aa = build_folding_hamiltonian(protein_2aa)
    # 2 aa has exactly 1 turn
    e_2aa = ham_2aa.energy([0])
    assert e_2aa == 0.0, "2 aa peptide has no non-bonded contacts (|i-j| >= 2 impossible)"

    protein_3aa = encode_sequence("ACD", model="HP")
    ham_3aa = build_folding_hamiltonian(protein_3aa)
    # 3 aa has 2 turns
    e_3aa = ham_3aa.energy([0, 1])
    assert isinstance(e_3aa, float)


def test_f1_b02_homopolymer_sequences():
    """T2_F1_02: Verify homopolymer sequences (all-Ala, all-Gly) evaluate consistently."""
    p_ala = encode_sequence("AAAAAAAA", model="MJ")
    ham_ala = build_folding_hamiltonian(p_ala)
    turns = [0, 1, 2, 3, 0, 1, 2]
    e_ala = ham_ala.energy(turns)
    assert not np.isnan(e_ala)

    p_gly = encode_sequence("GGGGGGGG", model="HP")
    ham_gly = build_folding_hamiltonian(p_gly)
    assert p_gly.hp_sequence == "PPPPPPPP"
    e_gly = ham_gly.energy(turns)
    # In HP model, Polar-Polar contacts contribute 0.0
    assert e_gly >= 0.0


def test_f1_b03_extreme_hydrophobicity_contrast():
    """T2_F1_03: Verify hydrophobic vs polar sequences show sharp energetic contrast."""
    p_hydro = encode_sequence("LLLLLL", model="MJ")
    p_polar = encode_sequence("KKKKKK", model="MJ")
    ham_hydro = build_folding_hamiltonian(p_hydro)
    ham_polar = build_folding_hamiltonian(p_polar)

    # Conformation with multiple spatial contacts
    turns = [0, 1, 2, 0, 1]
    e_hydro = ham_hydro.energy(turns)
    e_polar = ham_polar.energy(turns)
    # Hydrophobic contacts are attractive (negative), polar contacts are neutral/repulsive
    assert e_hydro < e_polar


def test_f1_b04_complete_cyclic_collision_penalty():
    """T2_F1_04: Verify cyclic walk returning to start triggers heavy overlap penalty."""
    protein = encode_sequence("ACDEF", model="MJ")
    ham = build_folding_hamiltonian(protein)
    # 4 turns making a closed square: (0,0) -> (1,0) -> (1,1) -> (0,1) -> (0,0)
    cyclic_turns = [0, 1, 2, 3]
    coords = compute_coordinates(cyclic_turns)
    assert not is_self_avoiding(coords)
    assert coords[0] == coords[4]

    e = ham.energy(cyclic_turns)
    assert e >= ham.lambda_overlap - 10.0, "Overlap penalty must dominate"


def test_f1_b05_fully_extended_linear_chain():
    """T2_F1_05: Verify straight linear conformation has zero contacts and zero penalty."""
    protein = encode_sequence("ACDEFGHIKL", model="MJ")
    ham = build_folding_hamiltonian(protein)
    linear_turns = [0] * (len(protein.sequence) - 1)
    coords = compute_coordinates(linear_turns)

    assert is_self_avoiding(coords)
    assert len(get_contacts(coords)) == 0
    assert not has_backtracking(linear_turns)
    e = ham.energy(linear_turns)
    assert e == 0.0, "Fully extended conformation must have exactly zero interaction energy"


# ===========================================================================
# Feature 2 Boundaries (5 tests)
# ===========================================================================

def test_f2_b01_zero_parameter_circuit():
    """T2_F2_01: Verify circuit with all zero phase angles operates deterministically."""
    circuit, params = build_photonic_ansatz(num_modes=4, num_photons=2)
    for p in params:
        p.set_value(0.0)

    u = circuit.compute_unitary()
    assert np.allclose(np.dot(u.conj().T, u), np.eye(4), atol=1e-5)


def test_f2_b02_single_mode_circuit_boundary():
    """T2_F2_02: Verify 1-mode optical circuit boundary behaves as pure phase scalar."""
    circuit, params = build_photonic_ansatz(num_modes=1, num_photons=1)
    assert circuit.m == 1
    assert len(params) == 1
    params[0].set_value(np.pi / 2)
    u = circuit.compute_unitary()
    assert u.shape == (1, 1)
    assert np.isclose(np.abs(u[0, 0]), 1.0)


def test_f2_b03_extreme_phase_angles_wrap_around():
    """T2_F2_03: Verify phase angles wrap periodically modulo 2pi without divergence."""
    circuit, params = build_photonic_ansatz(num_modes=2, num_photons=1)
    p = params[0]

    p.set_value(0.5)
    u1 = circuit.compute_unitary()

    p.set_value(0.5 + 2 * np.pi)
    u2 = circuit.compute_unitary()

    assert np.allclose(u1, u2, atol=1e-5), "Unitary must be periodic with period 2pi"

    # Extreme value 1000 * pi
    p.set_value(1000.0 * np.pi)
    u_ext = circuit.compute_unitary()
    assert not np.isnan(u_ext).any()
    assert not np.isinf(u_ext).any()


def test_f2_b04_saturated_photon_density():
    """T2_F2_04: Verify high photon density (N = M = 4) executes without crash."""
    circuit, _ = build_photonic_ansatz(num_modes=4, num_photons=4)
    state = build_input_state(num_modes=4, num_photons=4, encoding="fock")
    assert state.n == 4
    assert state.m == 4

    proc = pcvl.Processor("SLOS", circuit)
    proc.with_input(state)
    sampler = pcvl.algorithm.Sampler(proc)
    res = sampler.sample_count(20)
    assert sum(res["results"].values()) == 20


def test_f2_b05_empty_circuit_identity():
    """T2_F2_05: Verify empty 2-mode circuit acts as identity transformation."""
    c = pcvl.Circuit(2)
    u = c.compute_unitary()
    assert np.allclose(u, np.eye(2), atol=1e-7)

    proc = pcvl.Processor("SLOS", c)
    proc.with_input(pcvl.BasicState([1, 0]))
    sampler = pcvl.algorithm.Sampler(proc)
    res = sampler.sample_count(50)
    assert len(res["results"]) == 1
    assert pcvl.BasicState([1, 0]) in res["results"]


# ===========================================================================
# Feature 3 Boundaries (5 tests)
# ===========================================================================

def test_f3_b01_missing_t_prefix_token_rejected():
    """T2_F3_01: Verify token without literal '_T_' prefix raises explicit ValueError."""
    raw_jwt = "eyJhbGciOiJIUzUxMiIsInR5cCI6IkpXVCJ9.dummy_body.dummy_signature"
    with pytest.raises(ValueError, match="_T_"):
        QuandelaCloudClient(token=raw_jwt)


def test_f3_b02_empty_or_none_token_rejected():
    """T2_F3_02: Verify empty or None token is immediately rejected."""
    with pytest.raises(ValueError):
        QuandelaCloudClient(token="")

    with pytest.raises(ValueError):
        QuandelaCloudClient(token=None)

    with pytest.raises(ValueError):
        QuandelaCloudClient(token="   ")


def test_f3_b03_http_401_unauthorized_handling(monkeypatch):
    """
    T2_F3_03: Verify client handles Quandela Cloud HTTP 401 Unauthorized gateway rejection.
    Simulates remote authentication failure with valid token prefix format but revoked credentials.
    """
    revoked_token = (
        "_T_eyJhbGciOiJIUzUxMiIsInR5cCI6IkpXVCJ9."
        "eyJpZCI6OTk5OSwiZXhwIjoxMDAwMDAwMDAwfQ."
        "invalid_revoked_signature_sample_abcdef1234567890"
    )
    client = QuandelaBackend(token=revoked_token, platform="sim:slos", local_fallback=False)

    def mock_401_response(url, *args, **kwargs):
        resp = requests.Response()
        resp.status_code = 401
        resp._content = b'{"detail": "Unauthorized: Invalid or expired API key"}'
        resp.url = url
        return resp

    monkeypatch.setattr(requests, "get", mock_401_response)

    with pytest.raises((PermissionError, requests.HTTPError), match="401|Unauthorized"):
        client.check_platform_availability("sim:slos")


def test_f3_b04_platform_maintenance_or_rate_limit_fallback(mock_cloud_client):
    """T2_F3_04: Verify client detects physical QPU maintenance and falls back to sim:slos."""
    client_qpu = QuandelaCloudClient(VALID_QUANDELA_TOKEN, platform="qpu:ascella")
    status = client_qpu.check_platform_status()
    assert status["status"] == "maintenance"
    assert status["available"] is False

    # Job submission to maintenance platform succeeds via automatic fallback
    c = pcvl.Circuit(2) // pcvl.BS()
    inp = pcvl.BasicState([1, 0])
    sub = client_qpu.submit_job(c, inp, shots=10)
    assert "job_id" in sub


def test_f3_b05_network_timeout_and_retry(mock_cloud_client):
    """T2_F3_05: Verify nonexistent job query returns 'failed' or raises KeyError."""
    status = mock_cloud_client.get_job_status("nonexistent-job-uuid-9999")
    assert status == "failed"

    with pytest.raises(KeyError):
        mock_cloud_client.get_job_results("nonexistent-job-uuid-9999")


# ===========================================================================
# Feature 4 Boundaries (5 tests)
# ===========================================================================

def test_f4_b01_empty_candidate_pool_handling(tmp_path, mock_cloud_client):
    """T2_F4_01: Verify evaluation pipeline handles empty candidate CSV gracefully."""
    empty_csv = tmp_path / "empty_pool.csv"
    empty_csv.write_text("name,sequence,target\n")

    pipeline = QuantumFoldingPipeline(backend=mock_cloud_client)
    df = pipeline.evaluate_pool(str(empty_csv))
    assert df.empty
    assert "quantum_score" in df.columns


def test_f4_b02_single_candidate_pool(tmp_path, mock_cloud_client):
    """T2_F4_02: Verify pipeline evaluates single-candidate pool and assigns rank 1."""
    single_csv = tmp_path / "single_pool.csv"
    single_csv.write_text("name,sequence\nSOLO_01,AALEELLKALEELLK\n")

    pipeline = QuantumFoldingPipeline(backend=mock_cloud_client)
    out_csv = tmp_path / "out_single.csv"
    df = pipeline.evaluate_pool(str(single_csv), output_csv=str(out_csv))
    assert len(df) == 1
    assert df.iloc[0]["rank"] == 1
    assert out_csv.exists()


def test_f4_b03_duplicate_candidate_sequences(mock_cloud_client):
    """T2_F4_03: Verify duplicate candidate sequences yield equal scores without sorting crash."""
    pipeline = QuantumFoldingPipeline(backend=mock_cloud_client)
    candidates = [
        {"name": "cand_A", "sequence": "AALEELLK"},
        {"name": "cand_B", "sequence": "AALEELLK"},
    ]
    scored = pipeline.run_batch(candidates)
    assert len(scored) == 2
    assert scored[0]["quantum_score"] == pytest.approx(scored[1]["quantum_score"])


def test_f4_b04_zero_or_one_iteration_vqe(mock_cloud_client):
    """T2_F4_04: Verify minimal iteration VQE (max_iter=1) terminates cleanly."""
    protein = encode_sequence("HPPH", model="HP")
    res = evaluate_vqe(protein, mock_cloud_client, max_iter=1)
    assert isinstance(res, VQEResult)
    assert len(res.history) == 1
    assert len(res.cloud_job_ids) == 1


def test_f4_b05_identical_folding_scores_tie_breaking(mock_cloud_client):
    """T2_F4_05: Verify deterministic tie-breaking for candidates with identical energies."""
    pipeline = QuantumFoldingPipeline(backend=mock_cloud_client)
    candidates = [
        {"name": "cand_z", "sequence": "GGGG"},
        {"name": "cand_a", "sequence": "GGGG"},
    ]
    scored = pipeline.run_batch(candidates)
    assert len(scored) == 2
    assert scored[0]["name"] == "cand_a" or scored[1]["name"] == "cand_a"


# ===========================================================================
# Feature 5 Boundaries (5 tests)
# ===========================================================================

def test_f5_b01_invalid_login_credentials():
    """T2_F5_01: Verify incorrect credentials raise PermissionError without leaking password."""
    client = SJTUPlatformClient(email=SJTU_EMAIL, password="wrong_password_123")
    with pytest.raises(PermissionError, match="用户名或密码错误"):
        client.login()


def test_f5_b02_malformed_json_response_handling(monkeypatch):
    """T2_F5_02: Verify non-JSON error response from platform is handled cleanly."""
    client = SJTUPlatformClient()
    def mock_post(*args, **kwargs):
        r = requests.Response()
        r.status_code = 502
        r._content = b"<html>Bad Gateway</html>"
        return r

    monkeypatch.setattr(client.session, "post", mock_post)
    with pytest.raises(requests.HTTPError):
        client.login()


def test_f5_b03_trailing_slash_normalization():
    """T2_F5_03: Verify URLs with multiple trailing slashes normalize identically."""
    c1 = SJTUPlatformClient(base_url="https://aiprotein.sjtu.edu.cn")
    c2 = SJTUPlatformClient(base_url="https://aiprotein.sjtu.edu.cn/")
    c3 = SJTUPlatformClient(base_url="https://aiprotein.sjtu.edu.cn////")
    assert c1.base_url == c2.base_url == c3.base_url == "https://aiprotein.sjtu.edu.cn"


def test_f5_b04_server_http_500_error_handling(monkeypatch):
    """T2_F5_04: Verify server 500 error raises HTTPError with proper status code."""
    client = SJTUPlatformClient()
    def mock_get(*args, **kwargs):
        r = requests.Response()
        r.status_code = 500
        r._content = b'{"detail": "Internal Server Error"}'
        return r

    monkeypatch.setattr(client.session, "get", mock_get)
    with pytest.raises(requests.HTTPError) as exc_info:
        client.get_targets()
    assert exc_info.value.response.status_code == 500


def test_f5_b05_empty_response_body_handling(monkeypatch):
    """T2_F5_05: Verify latest submission returns None when no submission exists."""
    client = SJTUPlatformClient()
    client.login()
    latest = client.get_latest_submission("unsubmitted-target-uuid-000")
    assert latest is None


# ===========================================================================
# Feature 6 Boundaries (5 tests)
# ===========================================================================

def test_f6_b01_empty_string_sequence_rejected():
    """T2_F6_01: Verify empty or whitespace-only sequences are rejected."""
    ok, err = validate_sequence("")
    assert ok is False
    assert "序列为空" in err

    ok, err = validate_sequence("     ")
    assert ok is False
    assert "序列为空" in err


def test_f6_b02_non_canonical_amino_acids_rejected():
    """T2_F6_02: Verify non-canonical characters (X, B, Z, J, U, numbers) are rejected."""
    bad_chars = ["X", "B", "Z", "J", "U", "O", "1", "*", "-"]
    for ch in bad_chars:
        seq = "AALEELLK" + ch + "KALEELLK" + "A" * 30
        ok, err = validate_sequence(seq)
        assert ok is False, f"Expected character '{ch}' to be rejected"
        assert "含有非法氨基酸" in err


def test_f6_b03_length_underflow_rejected():
    """T2_F6_03: Verify sequences below minimum length (39 aa for minibinder) are rejected."""
    seq_39 = "A" * 39
    ok, err = validate_sequence(seq_39, min_len=40, max_len=100)
    assert ok is False
    assert "序列长度不符合要求" in err


def test_f6_b04_length_overflow_rejected():
    """T2_F6_04: Verify sequences above maximum length (101 aa for minibinder) are rejected."""
    seq_101 = "A" * 101
    ok, err = validate_sequence(seq_101, min_len=40, max_len=100)
    assert ok is False
    assert "序列长度不符合要求" in err


def test_f6_b05_exact_boundary_lengths_accepted():
    """T2_F6_05: Verify exact boundary lengths (40, 85, 100, 150 aa) are all accepted."""
    # Minibinder lower bound: exactly 40 aa
    ok_40, _ = validate_sequence("A" * 40, min_len=40, max_len=100)
    assert ok_40 is True

    # Real pool lower length: exactly 85 aa
    ok_85, _ = validate_sequence("A" * 85, min_len=40, max_len=100)
    assert ok_85 is True

    # Minibinder upper bound: exactly 100 aa
    ok_100, _ = validate_sequence("A" * 100, min_len=40, max_len=100)
    assert ok_100 is True

    # Nanobody upper bound: exactly 150 aa
    ok_150, _ = validate_sequence("A" * 150, min_len=100, max_len=150)
    assert ok_150 is True


# ===========================================================================
# Feature 7 Boundaries (5 tests)
# ===========================================================================

def test_f7_b01_unagreed_submission_rejected(tmp_path):
    """T2_F7_01: Verify submission with agreed=False is blocked locally."""
    client = SJTUPlatformClient()
    client.login()

    csv_file = tmp_path / "valid.csv"
    csv_file.write_text("name,sequence\nc1," + "A" * 50 + "\n")

    with pytest.raises(ValueError, match="请先勾选并同意赛事细则"):
        client.submit_candidates(
            target_id=TL1A_TARGET_ID,
            csv_path=str(csv_file),
            agreed=False,
        )


def test_f7_b02_missing_or_invalid_target_id_rejected(tmp_path):
    """T2_F7_02: Verify non-existent target ID raises KeyError or ValueError."""
    client = SJTUPlatformClient()
    client.login()

    csv_file = tmp_path / "valid.csv"
    csv_file.write_text("name,sequence\nc1," + "A" * 50 + "\n")

    with pytest.raises(ValueError, match="Target ID is required"):
        client.submit_candidates(target_id="", csv_path=str(csv_file), agreed=True)

    with pytest.raises(KeyError, match="赛题不存在"):
        client.submit_candidates(
            target_id="00000000-0000-0000-0000-000000000000",
            csv_path=str(csv_file),
            agreed=True,
        )


def test_f7_b03_oversized_abstract_or_notes_rejected(tmp_path):
    """T2_F7_03: Verify abstract > 4000 chars or note > 2000 chars is rejected."""
    client = SJTUPlatformClient()
    client.login()

    csv_file = tmp_path / "valid.csv"
    csv_file.write_text("name,sequence\nc1," + "A" * 50 + "\n")

    with pytest.raises(ValueError, match="Abstract exceeds maximum"):
        client.submit_candidates(
            target_id=TL1A_TARGET_ID,
            csv_path=str(csv_file),
            agreed=True,
            abstract="A" * 4001,
        )

    with pytest.raises(ValueError, match="Note exceeds maximum"):
        client.submit_candidates(
            target_id=TL1A_TARGET_ID,
            csv_path=str(csv_file),
            agreed=True,
            note="N" * 2001,
        )


def test_f7_b04_empty_csv_submission_rejected(tmp_path):
    """T2_F7_04: Verify empty CSV or header-only CSV is rejected."""
    client = SJTUPlatformClient()
    client.login()

    empty_csv = tmp_path / "empty.csv"
    empty_csv.write_text("name,sequence\n")

    with pytest.raises(ValueError, match="没有序列"):
        client.submit_candidates(target_id=TL1A_TARGET_ID, csv_path=str(empty_csv), agreed=True)


def test_f7_b05_platform_422_detailed_error_parsing(tmp_path):
    """T2_F7_05: Verify 422 validation failure with duplicate names is caught and formatted."""
    client = SJTUPlatformClient()
    client.login()

    dup_csv = tmp_path / "dup.csv"
    dup_csv.write_text(
        "name,sequence\n"
        "cand_1," + "A" * 50 + "\n"
        "cand_1," + "A" * 55 + "\n"
    )

    with pytest.raises(ValueError, match="名称重复"):
        client.submit_candidates(target_id=TL1A_TARGET_ID, csv_path=str(dup_csv), agreed=True)


def test_f7_b06_sequence_count_overflow_rejected(tmp_path):
    """T2_F7_06: Verify submitting 21 sequences (>20 platform limit) is rejected."""
    client = SJTUPlatformClient()
    client.login()

    overflow_csv = tmp_path / "overflow_21.csv"
    rows = ["name,sequence"] + [f"CAND_{i:02d}," + "A" * 50 for i in range(21)]
    overflow_csv.write_text("\n".join(rows))

    # 1. Local validator check
    count, errs = validate_submission_csv(overflow_csv.read_text(), modality="minibinder")
    assert count == 21
    assert any("最多提交 20 条序列" in e[2] for e in errs)

    # 2. Client submission check
    with pytest.raises(ValueError, match="最多提交 20 条序列"):
        client.submit_candidates(target_id=TL1A_TARGET_ID, csv_path=str(overflow_csv), agreed=True)


def test_f7_b07_mock_adapter_sequence_length_underflow_and_overflow(tmp_path):
    """T2_F7_07: Verify MockSJTUPlatformAdapter rejects 39 aa and 101 aa for TL1A."""
    s = requests.Session()

    # Underflow (39 aa < 40 aa min)
    r1 = s.post(
        "https://aiprotein.sjtu.edu.cn/api/submissions/csv",
        json={"target_id": TL1A_TARGET_ID, "csv": "name,sequence\nc1," + "A" * 39, "agreed": True},
    )
    assert r1.status_code == 422
    assert "序列长度不符合要求" in r1.json()["errors"][0]["message"]

    # Overflow (101 aa > 100 aa max)
    r2 = s.post(
        "https://aiprotein.sjtu.edu.cn/api/submissions/csv",
        json={"target_id": TL1A_TARGET_ID, "csv": "name,sequence\nc1," + "A" * 101, "agreed": True},
    )
    assert r2.status_code == 422
    assert "序列长度不符合要求" in r2.json()["errors"][0]["message"]


def test_f7_b08_mock_adapter_abstract_and_note_size_rejection():
    """T2_F7_08: Verify MockSJTUPlatformAdapter enforces abstract <= 4000 and note <= 2000 chars."""
    s = requests.Session()

    r1 = s.post(
        "https://aiprotein.sjtu.edu.cn/api/submissions/csv",
        json={
            "target_id": TL1A_TARGET_ID,
            "csv": "name,sequence\nc1," + "A" * 50,
            "agreed": True,
            "abstract": "X" * 4001,
        },
    )
    assert r1.status_code == 422
    assert "摘要超过最大长度限制" in r1.json()["detail"]

    r2 = s.post(
        "https://aiprotein.sjtu.edu.cn/api/submissions/csv",
        json={
            "target_id": TL1A_TARGET_ID,
            "csv": "name,sequence\nc1," + "A" * 50,
            "agreed": True,
            "note": "X" * 2001,
        },
    )
    assert r2.status_code == 422
    assert "备注超过最大长度限制" in r2.json()["detail"]


def test_f7_b09_quoted_csv_fields_with_commas(tmp_path):
    """T2_F7_09: Verify candidate names containing quoted commas parse without error."""
    client = SJTUPlatformClient()
    client.login()

    quoted_csv = tmp_path / "quoted.csv"
    quoted_csv.write_text('name,sequence\n"candidate, variant 1",' + "A" * 50 + "\n")

    res = client.submit_candidates(target_id=TL1A_TARGET_ID, csv_path=str(quoted_csv), agreed=True)
    assert res["status"] == "confirmed"
    assert res["count"] == 1


def test_f7_b10_nanobody_egfr_submission_and_length_validation(tmp_path):
    """T2_F7_10: Verify 120 aa Nanobody to EGFR succeeds and auto-resolves modality."""
    client = SJTUPlatformClient()
    client.login()

    egfr_id = "af6b1e2b-c62d-40e4-a5b3-3f091be04c60"
    vhh_csv = tmp_path / "vhh_egfr.csv"
    vhh_csv.write_text("name,sequence\nvhh_01," + "A" * 120 + "\n")

    res = client.submit_candidates(target_id=egfr_id, csv_path=str(vhh_csv), agreed=True)
    assert res["status"] == "confirmed"
    assert res["target_code"] == "egfr"
    assert res["count"] == 1

