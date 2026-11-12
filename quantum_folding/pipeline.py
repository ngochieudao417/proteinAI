"""
quantum_folding.pipeline
========================
Quantum Protein Folding Candidate Screening Pipeline.
Extracts binding epitope hotspots, evaluates folding stability via VQE,
computes baseline energy gaps, and exports multi-candidate rankings.
"""

from __future__ import annotations

import itertools
import logging
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd

from quantum_folding.circuit import build_ansatz_circuit
from quantum_folding.hamiltonian import FoldingHamiltonian
from quantum_folding.lattice import DIR_2D, LatticeModel, LatticeProtein
from quantum_folding.quandela_backend import QuandelaBackend
from quantum_folding.vqe import VQEResult, evaluate_vqe

logger = logging.getLogger("QuantumFoldingPipeline")

HYDROPHOBIC_RESIDUES = frozenset("ACFILMVWY")


def extract_hotspot_window(
    sequence: str,
    window_size: int = 6,
    strategy: str = "hydrophobic_density",
) -> Tuple[str, int]:
    """
    Extract critical binding hotspot / epitope loop from candidate sequence.

    Args:
        sequence: Full amino acid sequence (85-100 aa).
        window_size: Length of window to extract (default 6).
        strategy: 'hydrophobic_density' or 'aromatic_core'.

    Returns:
        (window_sequence, start_index)
    """
    if len(sequence) <= window_size:
        return sequence, 0

    best_score = -1.0
    best_idx = 0
    n = len(sequence)

    for i in range(n - window_size + 1):
        win = sequence[i : i + window_size]
        if strategy == "hydrophobic_density":
            score = sum(1.0 for aa in win if aa in HYDROPHOBIC_RESIDUES)
        elif strategy == "aromatic_core":
            weights = {"W": 3.0, "F": 2.5, "Y": 2.5, "L": 2.0, "I": 2.0, "M": 1.5, "V": 1.5}
            score = sum(weights.get(aa, 0.5 if aa in HYDROPHOBIC_RESIDUES else 0.0) for aa in win)
        else:
            score = 0.0

        if score > best_score:
            best_score = score
            best_idx = i

    return sequence[best_idx : best_idx + window_size], best_idx


def compute_random_walk_baseline(
    window_seq: str,
    model: str = "MJ",
    lattice_type: Union[LatticeModel, str] = LatticeModel.SQUARE_2D,
) -> Tuple[float, float, int]:
    """
    Compute exact or sampled random self-avoiding walk (SAW) energy baseline.

    For window_size=6 in 2D square lattice, enumerates all 71 valid SAWs in <1ms.
    """
    num_residues = len(window_seq)
    num_turns = num_residues - 2
    lattice = LatticeModel(lattice_type) if isinstance(lattice_type, str) else lattice_type

    protein = LatticeProtein(sequence=window_seq, lattice_model=lattice)
    hamiltonian = FoldingHamiltonian(protein=protein, model=model)

    energies: List[float] = []

    if num_turns <= 4:
        # Full enumeration: 4^4 = 256 configurations
        for turns in itertools.product(range(4), repeat=num_turns):
            full_turns = (0,) + turns
            coords = protein.compute_coordinates(full_turns)
            if protein.is_self_avoiding(coords):
                energy = hamiltonian.energy(full_turns)
                energies.append(energy)
    else:
        # Monte Carlo sampling for larger windows
        rng = np.random.default_rng(42)
        samples = 0
        while len(energies) < 100 and samples < 5000:
            samples += 1
            turns = tuple(rng.integers(0, 4, size=num_turns))
            full_turns = (0,) + turns
            coords = protein.compute_coordinates(full_turns)
            if protein.is_self_avoiding(coords):
                energies.append(hamiltonian.energy(full_turns))

    if not energies:
        return 0.0, 1.0, 0

    mean_e = float(np.mean(energies))
    std_e = float(np.std(energies)) or 1.0
    return mean_e, std_e, len(energies)


def calculate_stability_score(
    vqe_energy: float,
    baseline_mean: float,
    baseline_std: float,
) -> float:
    """
    Compute normalized folding stability score in [0.0, 1.0].
    Higher score indicates greater stability relative to random ensemble.
    """
    delta_e = baseline_mean - vqe_energy
    scale = max(baseline_std, 1.0)
    # Logistic sigmoid
    return float(1.0 / (1.0 + np.exp(-delta_e / scale)))


@dataclass
class CandidateEvaluation:
    """Evaluation output for a single protein candidate."""
    name: str
    target: str
    full_length: int
    window_seq: str
    window_start: int
    vqe_energy: float
    random_baseline_mean: float
    random_baseline_std: float
    energy_gap: float
    quantum_folding_stability_score: float
    optimal_turns: str
    is_self_avoiding: bool
    num_contacts: int
    cloud_job_id: str
    execution_platform: str
    timestamp: str


class QuantumFoldingPipeline:
    """High-level batch candidate folding evaluator."""

    def __init__(
        self,
        backend: Optional[QuandelaBackend] = None,
        window_size: int = 6,
        model: str = "MJ",
        vqe_max_iter: int = 15,
    ) -> None:
        self.backend = backend or QuandelaBackend(platform="sim:slos", local_fallback=True)
        self.window_size = window_size
        self.model = model
        self.vqe_max_iter = vqe_max_iter

    def evaluate_candidate(self, cand: Dict[str, Any]) -> CandidateEvaluation:
        """Run sliding window extraction, VQE, and stability scoring on one candidate."""
        name = cand.get("name", "unknown")
        target = cand.get("target", "unknown")
        seq = cand.get("sequence", "")

        # Step 1: Extract Hotspot Window
        win_seq, start_idx = extract_hotspot_window(seq, window_size=self.window_size)

        # Step 2: Compute Baseline
        base_mean, base_std, _ = compute_random_walk_baseline(win_seq, model=self.model)

        # Step 3: Run Quantum VQE
        protein = LatticeProtein(sequence=win_seq, lattice_model=LatticeModel.SQUARE_2D)
        vqe_res = evaluate_vqe(
            protein=protein,
            backend=self.backend,
            max_iter=self.vqe_max_iter,
            model=self.model,
        )

        opt_energy = vqe_res.optimal_energy
        opt_turns = vqe_res.optimal_turns
        coords = protein.compute_coordinates(opt_turns)
        saw = protein.is_self_avoiding(coords)
        contacts = len(protein.get_contacts(coords))
        gap = base_mean - opt_energy
        score = calculate_stability_score(opt_energy, base_mean, base_std)

        job_ids = ",".join(vqe_res.cloud_job_ids) if vqe_res.cloud_job_ids else "local"

        return CandidateEvaluation(
            name=name,
            target=target,
            full_length=len(seq),
            window_seq=win_seq,
            window_start=start_idx,
            vqe_energy=round(opt_energy, 4),
            random_baseline_mean=round(base_mean, 4),
            random_baseline_std=round(base_std, 4),
            energy_gap=round(gap, 4),
            quantum_folding_stability_score=round(score, 4),
            optimal_turns=str(list(opt_turns)),
            is_self_avoiding=saw,
            num_contacts=contacts,
            cloud_job_id=job_ids,
            execution_platform=self.backend.platform,
            timestamp=datetime.now(timezone.utc).isoformat(),
        )

    def evaluate_pool(
        self,
        candidates: Union[str, Path, pd.DataFrame, List[Dict[str, Any]]],
        max_candidates: Optional[int] = None,
        output_csv: Optional[Union[str, Path]] = "data/quantum_folding_scores.csv",
    ) -> pd.DataFrame:
        """Evaluate a batch of candidate sequences and export results."""
        if isinstance(candidates, (str, Path)):
            df_in = pd.read_csv(candidates)
            records = df_in.to_dict(orient="records")
        elif isinstance(candidates, pd.DataFrame):
            records = candidates.to_dict(orient="records")
        else:
            records = list(candidates)

        if max_candidates is not None:
            records = records[:max_candidates]

        logger.info("Starting quantum folding evaluation on %d candidates...", len(records))
        results: List[Dict[str, Any]] = []

        for idx, rec in enumerate(records, 1):
            logger.info("[%d/%d] Evaluating %s (%s)...", idx, len(records), rec.get("name"), rec.get("target"))
            res = self.evaluate_candidate(rec)
            results.append(asdict(res))

        df_out = pd.DataFrame(results)

        if output_csv:
            out_p = Path(output_csv)
            out_p.parent.mkdir(parents=True, exist_ok=True)
            df_out.to_csv(out_p, index=False)
            logger.info("Saved quantum folding scores to %s", out_p)

        return df_out
