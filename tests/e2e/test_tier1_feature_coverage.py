"""
tests/e2e/test_tier1_feature_coverage.py - Tier 1: Happy-Path Feature Coverage.
Covers Features 1 through 7 with exactly 5 isolated test cases each (35 tests total).
"""

import os
import re
from pathlib import Path
import numpy as np
import pytest
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
    MJ_AA_ORDER,
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
# Feature 1: Quantum Hamiltonian Encoding (5 tests)
# ===========================================================================

def test_f1_01_hp_sequence_encoding():
    """T1_F1_01: Verify HP sequence encoding maps canonical amino acids correctly."""
    seq = "ACDEFGHIKLMNPQRSTVWY"
    hp = encode_hp(seq)
    assert len(hp) == len(seq)
    # Check known hydrophobic vs polar
    assert hp[0] == "H"  # Alanine
    assert hp[1] == "H"  # Cysteine
    assert hp[2] == "P"  # Aspartic Acid
    assert hp[3] == "P"  # Glutamic Acid
    assert set(hp).issubset({"H", "P"})

    protein = encode_sequence("HPPH", model="HP")
    assert isinstance(protein, LatticeProtein)
    assert protein.hp_sequence == "HPPH"


def test_f1_02_mj_contact_matrix_potentials():
    """T1_F1_02: Verify Miyazawa-Jernigan 20x20 contact potentials symmetry and ordering."""
    # Symmetry check
    for a1 in "ACDEF":
        for a2 in "GHIKL":
            e1 = get_mj_potential(a1, a2)
            e2 = get_mj_potential(a2, a1)
            assert e1 == pytest.approx(e2), f"Asymmetry between {a1} and {a2}"

    # Cys-Cys disulfide contact is significantly more attractive than Lys-Lys
    e_cys = get_mj_potential("C", "C")
    e_lys = get_mj_potential("K", "K")
    assert e_cys < e_lys
    assert e_cys < 0.0, "Cys-Cys contact should be strongly negative"


def test_f1_03_cost_hamiltonian_scalar_energy():
    """T1_F1_03: Verify FoldingHamiltonian computes scalar float energy for conformations."""
    protein = encode_sequence("ACDEFGHIKL", model="MJ")
    hamiltonian = build_folding_hamiltonian(protein)

    turns = (0, 1, 2, 0, 1, 2, 3, 0, 1)
    energy = hamiltonian.energy(turns)
    assert isinstance(energy, float)
    assert not np.isnan(energy)
    assert not np.isinf(energy)


def test_f1_04_coordinate_self_avoidance_penalty():
    """T1_F1_04: Verify coordinate self-avoidance penalty enforcement on collisions."""
    protein = encode_sequence("HPPH", model="HP")
    hamiltonian = build_folding_hamiltonian(protein)

    # Turns [0, 1, 2, 3]: Right, Up, Left, Down -> 4th step returns to (0, 0)
    cyclic_turns = (0, 1, 2, 3)
    coords_cyclic = compute_coordinates(cyclic_turns)
    assert not is_self_avoiding(coords_cyclic), "Cyclic turns must detect duplicate coordinates"

    # Straight line turns [0, 0, 0, 0]
    straight_turns = (0, 0, 0, 0)
    coords_straight = compute_coordinates(straight_turns)
    assert is_self_avoiding(coords_straight), "Straight line must be self-avoiding"

    e_cyclic = hamiltonian.energy(cyclic_turns)
    e_straight = hamiltonian.energy(straight_turns)
    assert e_cyclic >= hamiltonian.lambda_overlap - 1.0
    assert e_cyclic > e_straight


def test_f1_05_direction_vector_translation():
    """T1_F1_05: Verify 2D lattice coordinate generation and Euclidean unit step lengths."""
    turns = [0, 1, 2, 3]
    coords = compute_coordinates(turns, lattice=LatticeModel.SQUARE_2D)
    expected = [(0, 0), (1, 0), (1, 1), (0, 1), (0, 0)]
    assert coords == expected

    # Step length between successive residues must be exactly 1
    for k in range(len(coords) - 1):
        c1 = coords[k]
        c2 = coords[k + 1]
        dist = abs(c1[0] - c2[0]) + abs(c1[1] - c2[1])
        assert dist == 1


# ===========================================================================
# Feature 2: Perceval Photonic Circuit Generation (5 tests)
# ===========================================================================

def test_f2_01_interferometer_generation():
    """T1_F2_01: Verify Clements linear optical interferometer synthesis for N modes."""
    circuit_4m, params_4m = build_photonic_ansatz(num_modes=4, num_photons=2)
    assert circuit_4m.m == 4
    assert len(params_4m) > 0

    circuit_6m, params_6m = build_photonic_ansatz(num_modes=6, num_photons=2)
    assert circuit_6m.m == 6
    assert len(params_6m) > len(params_4m)


def test_f2_02_parameterized_phase_shifters_update():
    """T1_F2_02: Verify parameterized phase shifters accept numerical value updates."""
    circuit, params = build_photonic_ansatz(num_modes=4, num_photons=2)
    p0 = params[0]
    p0.set_value(1.234)
    assert float(p0) == pytest.approx(1.234)

    # Unitary reflects parameter modification
    u1 = circuit.compute_unitary()
    p0.set_value(2.567)
    u2 = circuit.compute_unitary()
    assert not np.allclose(u1, u2), "Circuit unitary must update when parameter changes"


def test_f2_03_fock_state_input_creation():
    """T1_F2_03: Verify input Fock state generation in dual-rail and concentrated forms."""
    state_dr = build_input_state(num_modes=4, num_photons=2, encoding="dual_rail")
    assert state_dr.m == 4
    assert state_dr.n == 2
    assert state_dr[0] == 1 and state_dr[2] == 1

    state_fock = build_input_state(num_modes=4, num_photons=2, encoding="fock")
    assert state_fock.m == 4
    assert state_fock.n == 2
    assert state_fock[0] == 1 and state_fock[1] == 1


def test_f2_04_local_simulator_execution():
    """T1_F2_04: Verify local SLOS simulator executes circuit and exhibits Hong-Ou-Mandel bunching."""
    c = pcvl.Circuit(2) // pcvl.BS()
    proc = pcvl.Processor("SLOS", c)
    proc.with_input(pcvl.BasicState([1, 1]))
    sampler = pcvl.algorithm.Sampler(proc)
    res = sampler.sample_count(100)
    counts = res["results"]

    total = sum(counts.values())
    assert total == 100
    # In ideal 50:50 beam splitter with |1, 1> input, Hong-Ou-Mandel effect outputs |2, 0> and |0, 2>
    state_20 = pcvl.BasicState([2, 0])
    state_02 = pcvl.BasicState([0, 2])
    assert state_20 in counts or state_02 in counts


def test_f2_05_circuit_unitary_and_photon_conservation():
    """T1_F2_05: Verify interferometer matrix is unitary and preserves photon count."""
    circuit, _ = build_photonic_ansatz(num_modes=4, num_photons=2)
    u = circuit.compute_unitary()
    # Check U^dagger * U == I
    identity = np.eye(4, dtype=complex)
    product = np.dot(u.conj().T, u)
    assert np.allclose(product, identity, atol=1e-5), "Interferometer must be unitary"

    # Sample and verify photon conservation
    proc = pcvl.Processor("SLOS", circuit)
    proc.with_input(pcvl.BasicState([1, 0, 1, 0]))
    sampler = pcvl.algorithm.Sampler(proc)
    res = sampler.sample_count(50)
    for state, cnt in res["results"].items():
        assert state.n == 2, f"State {state} violated photon number conservation"


# ===========================================================================
# Feature 3: Quandela Cloud Backend Client (5 tests)
# ===========================================================================

def test_f3_01_client_init_token_validation(mock_quandela_token):
    """T1_F3_01: Verify client initialization succeeds with valid '_T_' prefixed token."""
    client = QuandelaCloudClient(token=mock_quandela_token, platform="sim:slos")
    assert client.token == mock_quandela_token
    assert client.platform == "sim:slos"
    assert "quandela.com" in client.base_url


def test_f3_02_platform_configuration(mock_cloud_client):
    """T1_F3_02: Verify client platform configuration and status querying."""
    status_sim = mock_cloud_client.check_platform_status("sim:slos")
    assert status_sim["status"] == "available"
    assert status_sim["available"] is True

    status_qpu = mock_cloud_client.check_platform_status("qpu:ascella")
    assert status_qpu["status"] == "maintenance"
    assert status_qpu["available"] is False


def test_f3_03_remote_job_submission_response(mock_cloud_client):
    """T1_F3_03: Verify remote job submission returns valid job tracking structure."""
    circuit, _ = build_photonic_ansatz(num_modes=4, num_photons=2)
    input_state = build_input_state(num_modes=4, num_photons=2)
    sub = mock_cloud_client.submit_job(circuit, input_state, shots=50)

    assert "job_id" in sub
    assert isinstance(sub["job_id"], str)
    assert sub["status"] in ("completed", "waiting", "running")


def test_f3_04_job_status_polling(mock_cloud_client):
    """T1_F3_04: Verify client job status polling endpoint parses standard states."""
    circuit, _ = build_photonic_ansatz(num_modes=2, num_photons=1)
    input_state = pcvl.BasicState([1, 0])
    sub = mock_cloud_client.submit_job(circuit, input_state, shots=20)
    job_id = sub["job_id"]

    status = mock_cloud_client.get_job_status(job_id)
    assert status == "completed"

    unknown_status = mock_cloud_client.get_job_status("nonexistent-job-uuid")
    assert unknown_status == "failed"


def test_f3_05_job_results_retrieval(mock_cloud_client):
    """T1_F3_05: Verify job result retrieval extracts sample counts and performance metrics."""
    circuit, _ = build_photonic_ansatz(num_modes=2, num_photons=1)
    input_state = pcvl.BasicState([1, 0])
    sub = mock_cloud_client.submit_job(circuit, input_state, shots=50)
    job_id = sub["job_id"]

    res = mock_cloud_client.get_job_results(job_id)
    assert "results" in res
    assert "global_perf" in res
    assert isinstance(res["results"], dict)
    assert len(res["results"]) > 0


# ===========================================================================
# Feature 4: Quantum Folding Evaluation Pipeline (5 tests)
# ===========================================================================

def test_f4_01_evaluate_vqe_single_peptide(mock_cloud_client):
    """T1_F4_01: Verify VQE execution on single peptide returns complete VQEResult."""
    protein = encode_sequence("HPPH", model="HP")
    result = evaluate_vqe(protein, mock_cloud_client, max_iter=5)

    assert isinstance(result, VQEResult)
    assert isinstance(result.optimal_energy, float)
    assert len(result.optimal_turns) > 0
    assert len(result.cloud_job_ids) == 5
    assert len(result.history) == 5


def test_f4_02_batch_candidate_evaluation(mock_cloud_client):
    """T1_F4_02: Verify pipeline evaluates candidate batch and returns scored items."""
    pipeline = QuantumFoldingPipeline(backend=mock_cloud_client)
    candidates = [
        {"name": "c1", "sequence": "AALEELLK"},
        {"name": "c2", "sequence": "KKKKLLLL"},
        {"name": "c3", "sequence": "GGGGAAAA"},
    ]
    scored = pipeline.run_batch(candidates)
    assert len(scored) == 3
    assert all("quantum_score" in c for c in scored)
    # Ranked descending by quantum_score
    assert scored[0]["quantum_score"] >= scored[1]["quantum_score"] >= scored[2]["quantum_score"]


def test_f4_03_csv_output_generation(mock_cloud_client, temp_candidate_csv, temp_output_dir):
    """T1_F4_03: Verify evaluation pipeline reads pool and writes ranked output CSV."""
    pipeline = QuantumFoldingPipeline(backend=mock_cloud_client)
    out_csv = Path(temp_output_dir) / "quantum_folding_scores.csv"

    df = pipeline.evaluate_pool(temp_candidate_csv, output_csv=str(out_csv))
    assert out_csv.exists()
    assert "quantum_score" in df.columns
    assert "candidate_id" in df.columns
    assert "rank" in df.columns
    assert len(df) == 5


def test_f4_04_hybrid_scoring_weight_integration(mock_cloud_client):
    """T1_F4_04: Verify hybrid multi-objective scoring formula strictly respects weights."""
    pipeline = QuantumFoldingPipeline(backend=mock_cloud_client, weight_quantum=0.8, weight_classical=0.2)
    res = pipeline.evaluate_candidate("ACDEFY", candidate_id="test_cand")
    expected_combined = 0.8 * res["quantum_score"] + 0.2 * 1.0
    assert res["combined_score"] == pytest.approx(expected_combined, abs=1e-4)


def test_f4_05_vqe_optimizer_convergence(mock_cloud_client):
    """T1_F4_05: Verify VQE optimization history is non-increasing."""
    protein = encode_sequence("ACDEFY", model="MJ")
    res = evaluate_vqe(protein, mock_cloud_client, max_iter=8)

    # History tracks running best energy
    for i in range(len(res.history) - 1):
        assert res.history[i + 1] <= res.history[i] + 1e-9


# ===========================================================================
# Feature 5: SJTU Platform Client & Auth (5 tests)
# ===========================================================================

def test_f5_01_client_initialization_defaults_and_custom_base_url():
    """T1_F5_01: Verify SJTU client initializes and normalizes trailing slashes in base URL."""
    client = SJTUPlatformClient(email=SJTU_EMAIL, password=SJTU_PASSWORD, base_url="https://aiprotein.sjtu.edu.cn///")
    assert client.base_url == "https://aiprotein.sjtu.edu.cn"
    assert client.email == SJTU_EMAIL


def test_f5_02_login_authentication_success_and_token_storage():
    """T1_F5_02: Verify login authenticates contestant, stores tokens, and sets auth headers."""
    client = SJTUPlatformClient()
    auth_data = client.login()

    assert client.access_token is not None
    assert client.refresh_token is not None
    assert client.session.headers["Authorization"] == f"Bearer {client.access_token}"
    assert auth_data["user"]["team"]["name"] == "QAAA Lab"
    assert auth_data["user"]["name"] == "Mai Duc Khang"


def test_f5_03_get_targets_list_and_schema():
    """T1_F5_03: Verify targets query discovers all 7 tracks including TL1A and IL-7Ra."""
    client = SJTUPlatformClient()
    client.login()
    targets = client.get_targets()

    assert len(targets) == 7
    codes = {t["code"]: t for t in targets}
    assert "tl1a" in codes
    assert "il7ra" in codes

    tl1a = codes["tl1a"]
    assert tl1a["id"] == TL1A_TARGET_ID
    assert tl1a["modality"] == "minibinder"
    assert tl1a["min_length"] == 40
    assert tl1a["max_length"] == 100
    assert tl1a["max_sequences"] == 20


def test_f5_04_session_token_refresh():
    """T1_F5_04: Verify session token refresh updates bearer authorization header."""
    client = SJTUPlatformClient()
    client.login()
    old_token = client.access_token

    refresh_data = client.refresh_token_session()
    assert refresh_data["access_token"] != old_token
    assert client.access_token == refresh_data["access_token"]
    assert client.session.headers["Authorization"] == f"Bearer {client.access_token}"


def test_f5_05_get_latest_submission_and_history():
    """T1_F5_05: Verify queries for latest target submission and overall team history."""
    client = SJTUPlatformClient()
    client.login()

    # Query before any submission
    latest = client.get_latest_submission(TL1A_TARGET_ID)
    assert latest is None or isinstance(latest, dict)

    submissions = client.get_submissions(grouped=True)
    assert "groups" in submissions


# ===========================================================================
# Feature 6: Candidate Sequence Validation & Formatting (5 tests)
# ===========================================================================

def test_f6_01_canonical_20_amino_acids_validation():
    """T1_F6_01: Verify canonical 20 amino acid sequence passes validation without errors."""
    valid_seq = "ACDEFGHIKLMNPQRSTVWY" * 3  # 60 aa
    ok, err = validate_sequence(valid_seq, min_len=40, max_len=100)
    assert ok is True
    assert err is None


def test_f6_02_minibinder_length_validation_bounds(sample_peptides):
    """T1_F6_02: Verify minibinder length validation accepts candidates within 40-100 aa."""
    tl1a = sample_peptides["tl1a_sample"]
    assert 40 <= len(tl1a) <= 100
    ok, err = validate_sequence(tl1a, min_len=40, max_len=100)
    assert ok is True

    il7ra = sample_peptides["il7ra_sample"]
    assert 40 <= len(il7ra) <= 100
    ok, err = validate_sequence(il7ra, min_len=40, max_len=100)
    assert ok is True


def test_f6_03_nanobody_vhh_length_validation_bounds():
    """T1_F6_03: Verify nanobody/VHH track validates lengths within 100-150 aa."""
    vhh_seq = "A" * 125
    ok, err = validate_sequence(vhh_seq, min_len=100, max_len=150)
    assert ok is True

    short_vhh = "A" * 80
    ok, err = validate_sequence(short_vhh, min_len=100, max_len=150)
    assert ok is False
    assert "序列长度不符合要求" in err


def test_f6_04_csv_header_and_structure_validation():
    """T1_F6_04: Verify CSV header case-insensitivity and extra column tolerance."""
    csv_upper = "NAME,SEQUENCE\ncand_1," + "A" * 50 + "\n"
    count, errs = validate_submission_csv(csv_upper)
    assert count == 1
    assert len(errs) == 0

    csv_extra = "name,sequence,quantum_score,notes\ncand_1," + "A" * 50 + ",0.95,top\n"
    count, errs = validate_submission_csv(csv_extra)
    assert count == 1
    assert len(errs) == 0


def test_f6_05_sequence_sanitization_and_normalization():
    """T1_F6_05: Verify sanitization removes whitespace and handles UTF-8 BOM."""
    raw_csv = "\ufeffname,sequence\r\n cand_1 ,  aa leel  lkal  eell  kale  ellk  alep  whet  vmsi  nlfm  iefi  iahp  spke  iqit  welk  dflk  knfe  pvsk  sfhd  rlfp  fgdy  kwrs  plei  gk \r\n"
    count, errs = validate_submission_csv(raw_csv)
    assert count == 1
    assert len(errs) == 0


# ===========================================================================
# Feature 7: Candidate Sequence Submission & Feedback (5 tests)
# ===========================================================================

def test_f7_01_submission_payload_construction_and_contract(tmp_path):
    """T1_F7_01: Verify submission payload structure conforms to platform requirements."""
    client = SJTUPlatformClient()
    client.login()

    csv_file = tmp_path / "valid_tl1a.csv"
    csv_file.write_text("name,sequence\ncand_01," + "A" * 50 + "\n")

    res = client.submit_candidates(
        target_id=TL1A_TARGET_ID,
        csv_path=str(csv_file),
        action="confirm",
        agreed=True,
        abstract="Quantum folding optimization pipeline methodology",
        note="Initial candidate submission",
    )
    assert "id" in res
    assert res["target_id"] == TL1A_TARGET_ID
    assert res["count"] == 1


def test_f7_02_successful_submission_returns_batch_id_and_confirmation(tmp_path):
    """T1_F7_02: Verify successful submission returns batch ID and confirmed status."""
    client = SJTUPlatformClient()
    client.login()

    csv_file = tmp_path / "submission_batch.csv"
    rows = ["name,sequence"]
    for i in range(5):
        rows.append(f"CAND_{i:02d}," + "A" * (45 + i))
    csv_file.write_text("\n".join(rows))

    res = client.submit_candidates(
        target_id=TL1A_TARGET_ID,
        csv_path=str(csv_file),
        agreed=True,
    )
    assert res["status"] == "confirmed"
    assert res["count"] == 5
    assert res["id"].startswith("batch-")
    assert "confirmed_at" in res


def test_f7_03_target_id_binding_and_validation(tmp_path):
    """T1_F7_03: Verify correct target UUID binding for TL1A and IL-7Ra."""
    client = SJTUPlatformClient()
    client.login()

    csv_tl1a = tmp_path / "tl1a.csv"
    csv_tl1a.write_text("name,sequence\ntl1a_c1," + "A" * 50 + "\n")
    res_tl1a = client.submit_candidates(target_id=TL1A_TARGET_ID, csv_path=str(csv_tl1a), agreed=True)
    assert res_tl1a["target_code"] == "tl1a"

    csv_il7ra = tmp_path / "il7ra.csv"
    csv_il7ra.write_text("name,sequence\nil7ra_c1," + "A" * 50 + "\n")
    res_il7ra = client.submit_candidates(target_id=IL7RA_TARGET_ID, csv_path=str(csv_il7ra), agreed=True)
    assert res_il7ra["target_code"] == "il7ra"


def test_f7_04_submission_audit_trail_logging(tmp_path):
    """T1_F7_04: Verify client submission history queries retrieve newly submitted batch."""
    client = SJTUPlatformClient()
    client.login()

    csv_file = tmp_path / "audit_test.csv"
    csv_file.write_text("name,sequence\naudit_cand," + "A" * 60 + "\n")

    res = client.submit_candidates(target_id=TL1A_TARGET_ID, csv_path=str(csv_file), agreed=True)
    batch_id = res["id"]

    latest = client.get_latest_submission(TL1A_TARGET_ID)
    assert latest is not None
    assert latest["id"] == batch_id


def test_f7_05_re_upload_with_explanation_note(tmp_path):
    """T1_F7_05: Verify re-upload with explanatory note updates target submission."""
    client = SJTUPlatformClient()
    client.login()

    csv_file = tmp_path / "reupload.csv"
    csv_file.write_text("name,sequence\nreup_cand," + "A" * 65 + "\n")

    res = client.submit_candidates(
        target_id=IL7RA_TARGET_ID,
        csv_path=str(csv_file),
        agreed=True,
        note="Updated with refined quantum-scored conformations",
    )
    assert res["status"] == "confirmed"
    assert res["note"] == "Updated with refined quantum-scored conformations"
