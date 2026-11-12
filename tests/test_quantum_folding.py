"""
Unit and Integration Test Suite for Quantum Folding Package.
Covers lattice modeling, Hamiltonian formulation, circuit synthesis,
backend execution (local SLOS and cloud mock), VQE loop, and screening pipeline.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import perceval as pcvl

# Ensure project root is in sys.path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from quantum_folding.circuit import (
    build_ansatz_circuit,
    build_fock_state,
    decode_distribution,
    state_to_conformation,
)
from quantum_folding.hamiltonian import (
    FoldingHamiltonian,
    build_qubo_matrix,
    get_hp_energy,
    get_mj_energy,
)
from quantum_folding.lattice import (
    Direction2D,
    Direction3D,
    LatticeModel,
    LatticeProtein,
    count_backtracks,
    count_overlaps,
    encode_hp,
    has_backtracking,
    is_self_avoiding,
    turns_to_coordinates,
)
from quantum_folding.pipeline import (
    QuantumFoldingPipeline,
    calculate_stability_score,
    compute_random_walk_baseline,
    extract_hotspot_window,
)
from quantum_folding.quandela_backend import ExecutionResult, QuandelaBackend
from quantum_folding.vqe import VQEResult, compute_cvar, evaluate_vqe, run_vqe


# ============================================================================
# Section 1: Lattice & Encoding Tests
# ============================================================================

class TestLatticeModeling:
    def test_canonical_hp_encoding(self):
        """Verify 20 canonical amino acids correctly split into H and P."""
        hydrophobic = "ACFILMVWY"
        polar = "DEGHKNPQRST"
        assert encode_hp(hydrophobic) == "HHHHHHHHH"
        assert encode_hp(polar) == "PPPPPPPPPPP"

    def test_invalid_amino_acid_raises(self):
        """Invalid characters raise ValueError."""
        with pytest.raises(ValueError):
            encode_hp("ACFZ123")

    def test_2d_square_coordinates_straight(self):
        """Straight turns generate linear bead coordinates."""
        prot = LatticeProtein("HHHH", lattice_model=LatticeModel.SQUARE_2D)
        coords = prot.compute_coordinates((0, 0, 0))
        assert coords == [(0, 0), (1, 0), (2, 0), (3, 0)]
        assert prot.is_self_avoiding(coords) is True

    def test_self_avoiding_walk_detection(self):
        """Closed loop should be flagged as overlapping."""
        prot = LatticeProtein("HHHHH", lattice_model=LatticeModel.SQUARE_2D)
        loop_turns = (0, 1, 2, 3)
        coords = prot.compute_coordinates(loop_turns)
        assert coords[0] == coords[-1]
        assert prot.is_self_avoiding(coords) is False
        assert count_overlaps(coords) >= 1

    def test_contact_calculation(self):
        """Verify non-bonded topological contact detection (|i - j| >= 2, dist=1)."""
        prot = LatticeProtein("HHHH", lattice_model=LatticeModel.SQUARE_2D)
        u_turns = (0, 1, 2)
        coords = prot.compute_coordinates(u_turns)
        contacts = prot.get_contacts(coords)
        assert len(contacts) == 1
        assert (0, 3) in contacts or (3, 0) in contacts

    def test_backtracking_detection(self):
        """Immediate 180-degree reversal is detected."""
        prot = LatticeProtein("HHH", lattice_model=LatticeModel.SQUARE_2D)
        # Right then Left -> (0, 0) -> (1, 0) -> (0, 0)
        backtrack_turns = (0, 2)
        coords = prot.compute_coordinates(backtrack_turns)
        assert has_backtracking(coords) is True
        assert count_backtracks(coords) == 1

    def test_3d_cubic_coordinates(self):
        """3D cubic coordinates are generated correctly."""
        prot = LatticeProtein("HHHH", lattice_model=LatticeModel.CUBIC_3D)
        turns = (Direction3D.POS_X, Direction3D.POS_Y, Direction3D.POS_Z)
        coords = prot.compute_coordinates(turns)
        assert len(coords) == 4
        assert coords[0] == (0, 0, 0)
        assert coords[1] == (1, 0, 0)
        assert coords[2] == (1, 1, 0)
        assert coords[3] == (1, 1, 1)


# ============================================================================
# Section 2: Hamiltonian & Energy Tests
# ============================================================================

class TestFoldingHamiltonian:
    def test_mj_matrix_symmetry(self):
        """Miyazawa-Jernigan contact matrix must be strictly symmetric."""
        aas = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
        for a1 in aas:
            for a2 in aas:
                assert get_mj_energy(a1, a2) == pytest.approx(get_mj_energy(a2, a1), abs=1e-5)

    def test_mj_favorable_aromatic_contact(self):
        """Aromatic and hydrophobic pairs have negative contact potential."""
        assert get_mj_energy("C", "C") < -3.0
        assert get_mj_energy("F", "F") < -1.5
        assert get_mj_energy("L", "I") < -1.5

    def test_hp_energy_model(self):
        """HP energy model yields -1.0 for H-H and 0.0 otherwise."""
        assert get_hp_energy("H", "H") == -1.0
        assert get_hp_energy("H", "P") == 0.0
        assert get_hp_energy("P", "P") == 0.0

    def test_overlap_penalty_dominates(self):
        """Overlapping conformations must incur massive energy penalty."""
        prot = LatticeProtein("HHHHH", lattice_model=LatticeModel.SQUARE_2D)
        h = FoldingHamiltonian(protein=prot, model="HP", lambda_overlap=100.0)
        valid_energy = h.energy((0, 0, 0, 0))
        overlap_energy = h.energy((0, 1, 2, 3))
        assert overlap_energy > valid_energy + 50.0

    def test_qubo_matrix_structure(self):
        """QUBO matrix matches variable dimensions."""
        seq = "HPPHP"
        q, offset = build_qubo_matrix(seq, model="HP")
        num_turns = len(seq) - 2  # 3 turns
        expected_vars = num_turns * 4  # 12 variables
        assert q.shape == (expected_vars, expected_vars)
        assert isinstance(offset, float)


# ============================================================================
# Section 3: Photonic Circuit & Dual-Rail Tests
# ============================================================================

class TestPhotonicCircuit:
    def test_ansatz_circuit_creation(self):
        """GenericInterferometer should build with correct mode and parameter count."""
        circuit, params = build_ansatz_circuit(num_modes=4, shape="rectangle")
        assert circuit.m == 4
        assert len(params) >= 6
        # All parameters defined
        assert circuit.defined

    def test_reck_interferometer_shape(self):
        """Reck (triangular) decomposition is supported."""
        circuit, params = build_ansatz_circuit(num_modes=4, design="reck")
        assert circuit.m == 4
        assert len(params) >= 6

    def test_fock_state_generation(self):
        """Dual-rail Fock state generation places photons in even modes."""
        state = build_fock_state(num_modes=6, num_photons=3, pattern="dual_rail")
        assert state[0] == 1 and state[2] == 1 and state[4] == 1
        assert state[1] == 0 and state[3] == 0 and state[5] == 0

    def test_fock_state_to_conformation_mapping(self):
        """Dual-rail Fock states map correctly to discrete turns."""
        state = pcvl.BasicState([1, 0, 0, 1])
        turns = state_to_conformation(state, num_turns=2)
        assert len(turns) == 2
        assert all(isinstance(t, int) for t in turns)

    def test_decode_distribution(self):
        """Distribution is correctly decoded and normalized."""
        state1 = pcvl.BasicState([1, 0, 1, 0])
        state2 = pcvl.BasicState([0, 1, 0, 1])
        raw = {state1: 30, state2: 70}
        dist = decode_distribution(raw, num_turns=2)
        assert sum(dist.values()) == pytest.approx(1.0)
        assert len(dist) == 2


# ============================================================================
# Section 4: Quandela Backend Client Tests
# ============================================================================

class TestQuandelaBackend:
    def test_token_validation_literal_prefix(self):
        """Token lacking literal '_T_' prefix raises ValueError."""
        with pytest.raises(ValueError, match="must start with literal '_T_'"):
            QuandelaBackend(token="eyJhbGciOiJIUzUxMiIsInR5cCI6IkpXVCJ9...", local_fallback=False)

    def test_empty_token_error(self):
        """Empty string token raises ValueError."""
        with pytest.raises(ValueError):
            QuandelaBackend(token="", local_fallback=False)

    def test_local_slos_execution(self):
        """Local SLOS simulation mode executes synchronously without network."""
        backend = QuandelaBackend(platform="local:slos")
        c = pcvl.Circuit(2) // pcvl.BS()
        state = pcvl.BasicState([1, 1])
        res = backend.run(circuit=c, input_state=state, shots=50)

        assert isinstance(res, ExecutionResult)
        assert res.status == "SUCCESS"
        assert res.platform == "local:slos"
        assert sum(res.counts.values()) == 50
        # HOM bunching verification
        assert res.counts.get("|2,0>", 0) > 0
        assert res.counts.get("|0,2>", 0) > 0
        assert res.counts.get("|1,1>", 0) == 0

    def test_platform_availability_offline_safe(self):
        """Availability check handles local simulator mode safely."""
        backend = QuandelaBackend(platform="local:slos")
        info = backend.check_platform_availability()
        assert info["status"] == "available"
        assert info["type"] == "local_simulator"

    def test_run_circuit_dict_return(self):
        """run_circuit returns expected dictionary schema."""
        backend = QuandelaBackend(platform="local:slos")
        c = pcvl.Circuit(2) // pcvl.BS()
        state = pcvl.BasicState([1, 1])
        res = backend.run_circuit(c, state, shots=30)
        assert "results" in res
        assert "counts" in res
        assert "job_id" in res
        assert res["status"] == "SUCCESS"


# ============================================================================
# Section 5: VQE & Screening Pipeline Tests
# ============================================================================

class TestVQEAndPipeline:
    def test_hotspot_window_extraction(self):
        """Hotspot window extraction captures hydrophobic core."""
        seq = "DDDDDDIILLFFDDDDDD"
        win, idx = extract_hotspot_window(seq, window_size=6)
        assert win == "IILLFF"
        assert idx == 6

    def test_random_walk_baseline(self):
        """Baseline computes deterministic mean and std across valid SAWs."""
        mean_e, std_e, n_saws = compute_random_walk_baseline("INLFMI", model="HP")
        assert n_saws == 71  # Exactly 71 valid 6-step SAWs in 2D square
        assert isinstance(mean_e, float)
        assert std_e >= 0.0

    def test_stability_score_bounds(self):
        """Normalized stability score is strictly within [0.0, 1.0]."""
        s_high = calculate_stability_score(vqe_energy=-10.0, baseline_mean=0.0, baseline_std=2.0)
        s_low = calculate_stability_score(vqe_energy=5.0, baseline_mean=0.0, baseline_std=2.0)
        assert 0.0 <= s_high <= 1.0
        assert 0.0 <= s_low <= 1.0
        assert s_high > s_low

    def test_vqe_fast_execution(self):
        """Fast VQE optimization completes cleanly in local mode."""
        prot = LatticeProtein("HPPH")
        res = run_vqe(prot, max_iter=5, execution_mode="local")
        assert isinstance(res, VQEResult)
        assert len(res.best_conformation) == 2
        assert isinstance(res.best_energy, float)
        assert res.num_iterations >= 1

    def test_mini_candidate_screening(self, tmp_path):
        """Evaluate a mini pool and verify exported CSV schema."""
        backend = QuandelaBackend(platform="local:slos")
        pipeline = QuantumFoldingPipeline(backend=backend, window_size=6, vqe_max_iter=5)

        sample_candidates = [
            {"name": "TL1A_test_1", "target": "TL1A", "sequence": "AALEELLKALEELLKALEELLKALEPWHETVMSINLFMIEFIIAHPSPKEIQIT"},
            {"name": "IL7R_test_1", "target": "IL7RA", "sequence": "EKLRKLEKQARKQGDEVLVMLARMVLEYLEKGWVYSPFRHHETAEGRTQEMPNWW"},
        ]

        out_csv = tmp_path / "scores.csv"
        df = pipeline.evaluate_pool(sample_candidates, output_csv=out_csv)

        assert len(df) == 2
        assert out_csv.exists()
        expected_cols = {
            "name", "target", "full_length", "window_seq", "vqe_energy",
            "random_baseline_mean", "quantum_folding_stability_score", "cloud_job_id"
        }
        assert expected_cols.issubset(df.columns)
        assert all(0.0 <= s <= 1.0 for s in df["quantum_folding_stability_score"])


# ============================================================================
# Section 6: Live Cloud Opt-In Integration Test
# ============================================================================

@pytest.mark.cloud
@pytest.mark.skipif(not os.getenv("QUANDELA_TOKEN"), reason="Live cloud test requires QUANDELA_TOKEN")
def test_live_quandela_cloud_execution():
    """Live integration test against Quandela Cloud sim:slos."""
    token = os.getenv("QUANDELA_TOKEN")
    backend = QuandelaBackend(token=token, platform="sim:slos")
    c = pcvl.Circuit(2) // pcvl.BS()
    state = pcvl.BasicState([1, 1])
    res = backend.run(circuit=c, input_state=state, shots=20, timeout=60)
    assert res.status == "SUCCESS"
    assert res.job_id != ""
    assert res.counts.get("|2,0>", 0) + res.counts.get("|0,2>", 0) > 0
