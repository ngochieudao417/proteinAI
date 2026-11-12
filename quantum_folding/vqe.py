"""
quantum_folding.vqe
===================
Variational Photonic Eigensolver (VQE) for Lattice Protein Folding.
Perceval SDK 1.2.4 & Scipy Optimize compatible.
"""

from __future__ import annotations

import logging
import time
from dataclasses import dataclass, field
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np
import perceval as pcvl
from perceval.simulators import SimulatorFactory
import scipy.optimize as opt

from quantum_folding.circuit import (
    build_ansatz_circuit,
    build_fock_state,
    decode_distribution,
)

logger = logging.getLogger(__name__)


@dataclass
class VQEResult:
    """
    Standardized result data container for Photonic VQE optimization.

    Provides dual-interface compatibility:
    - Primary: best_energy, best_conformation, best_phases, history, cloud_job_ids
    - Interface aliases: optimal_energy, optimal_turns
    """
    best_energy: float
    best_conformation: Tuple[int, ...]
    best_phases: np.ndarray
    history: List[Dict[str, Any]] = field(default_factory=list)
    cloud_job_ids: List[str] = field(default_factory=list)
    optimal_distribution: Dict[Tuple[int, ...], float] = field(default_factory=dict)
    num_iterations: int = 0
    converged: bool = False
    execution_time_seconds: float = 0.0

    @property
    def optimal_energy(self) -> float:
        """Alias for PROJECT.md interface contract compliance."""
        return self.best_energy

    @property
    def optimal_turns(self) -> List[int]:
        """Alias for PROJECT.md interface contract compliance."""
        return list(self.best_conformation)


def compute_cvar(
    conf_distribution: Dict[Tuple[int, ...], float],
    energy_fn: Callable[[Tuple[int, ...]], float],
    alpha: float = 0.25,
) -> Tuple[float, float, Tuple[int, ...]]:
    """
    Computes the Conditional Value-at-Risk (CVaR) metric over lowest-energy conformations.

    :param conf_distribution: Normalized mapping from conformation tuple to probability.
    :param energy_fn: Evaluates Hamiltonian conformational energy E(conformation).
    :param alpha: Percentile tail fraction in (0.0, 1.0].
                  alpha=1.0 computes standard expected energy <H>.
                  alpha->0 computes ground-state minimum energy.
    :return: Tuple of (cvar_energy, best_sample_energy, best_sample_conformation).
    """
    if not conf_distribution:
        return 0.0, 0.0, ()

    # Evaluate energy for all sampled conformations
    confs_with_energy: List[Tuple[Tuple[int, ...], float, float]] = []
    best_e = float("inf")
    best_c: Tuple[int, ...] = ()

    for conf, prob in conf_distribution.items():
        e = energy_fn(conf)
        confs_with_energy.append((conf, e, prob))
        if e < best_e:
            best_e = e
            best_c = conf

    # Sort conformations ascending by energy (lowest energy first)
    confs_with_energy.sort(key=lambda item: item[1])

    # Bound alpha
    alpha_clamped = max(1e-6, min(1.0, float(alpha)))

    accumulated_prob = 0.0
    cvar_weighted_sum = 0.0

    for conf, e, prob in confs_with_energy:
        fraction_to_take = min(prob, alpha_clamped - accumulated_prob)
        cvar_weighted_sum += fraction_to_take * e
        accumulated_prob += fraction_to_take
        if accumulated_prob >= alpha_clamped - 1e-9:
            break

    cvar_value = cvar_weighted_sum / alpha_clamped
    return float(cvar_value), float(best_e), best_c


def run_vqe(
    protein: Any,
    hamiltonian: Optional[Any] = None,
    backend: Optional[Any] = None,
    max_iter: int = 30,
    method: str = "COBYLA",
    alpha: float = 0.25,
    shots: int = 100,
    execution_mode: str = "auto",
    bits_per_turn: int = 2,
    design: str = "clements",
    initial_phases: Optional[Sequence[float]] = None,
    seed: Optional[int] = None,
    model: str = "MJ",
) -> VQEResult:
    """
    Executes Variational Photonic Eigensolver (VQE) to fold a lattice protein.

    Supports:
    - Pure local simulation (Perceval SLOS backend)
    - Quandela Cloud remote execution (sim:slos / qpu:ascella)
    - Hybrid optimization: rapid local parameter convergence followed by cloud shot validation.

    :param protein: LatticeProtein instance or object with .sequence and .num_turns.
    :param hamiltonian: Optional FoldingHamiltonian instance (constructed if None).
    :param backend: Optional QuandelaBackend instance. If None, uses local Perceval Simulator.
    :param max_iter: Maximum classical optimizer iterations.
    :param method: Classical optimization algorithm ('COBYLA', 'Nelder-Mead', 'Powell').
    :param alpha: CVaR tail parameter in (0, 1]. alpha=1.0 is standard expected value.
    :param shots: Number of measurement shots for count sampling.
    :param execution_mode: 'auto' (hybrid if backend provided, else local), 'local', or 'cloud'.
    :param bits_per_turn: Number of decision bits per turn (default: 2).
    :param design: 'clements' or 'reck'.
    :param initial_phases: Optional initial phase parameters.
    :param seed: Random seed for reproducibility.
    :param model: Energy model ('MJ' or 'HP') if hamiltonian is auto-constructed.
    :return: VQEResult dataclass.
    """
    start_time = time.time()
    if seed is not None:
        np.random.seed(seed)

    # 1. Resolve protein turns and sequence
    if hasattr(protein, "num_turns"):
        num_turns = protein.num_turns
    elif hasattr(protein, "sequence"):
        num_turns = max(1, len(protein.sequence) - 2)
    else:
        raise ValueError("protein must have .num_turns or .sequence attribute")

    # 2. Resolve energy function
    if hamiltonian is not None and hasattr(hamiltonian, "energy"):
        energy_fn = hamiltonian.energy
    elif hasattr(protein, "energy"):
        energy_fn = protein.energy
    else:
        from quantum_folding.hamiltonian import FoldingHamiltonian
        ham_inst = FoldingHamiltonian(protein=protein, model=model)
        energy_fn = ham_inst.energy

    # 3. Setup Photonic Circuit
    # Adapt bits_per_turn and cap modes to 12 to respect Quandela Ascella hardware limit
    if bits_per_turn == 2 and num_turns > 2 and (num_turns * 2 * 2 > 12):
        effective_bpt = 1
    else:
        effective_bpt = bits_per_turn

    num_qubits = min(6, max(1, num_turns * effective_bpt))
    num_modes = num_qubits * 2
    num_photons = num_qubits

    circuit, params = build_ansatz_circuit(
        num_modes=num_modes,
        num_photons=num_photons,
        design=design,
        phase_shifter_layer=True,
        initial_phases=initial_phases,
    )
    input_state = build_fock_state(num_modes, num_photons, pattern="dual_rail")

    # 4. Determine execution strategy
    has_remote_backend = (
        backend is not None
        and hasattr(backend, "run_circuit")
        and not getattr(backend, "is_local", False)
    )

    if execution_mode == "auto":
        is_cloud_loop = False
        run_cloud_validation = has_remote_backend
    elif execution_mode == "cloud":
        if backend is None or not hasattr(backend, "run_circuit"):
            raise ValueError("execution_mode='cloud' requested but backend has no .run_circuit()")
        is_cloud_loop = True
        run_cloud_validation = False
    else:  # 'local'
        is_cloud_loop = False
        run_cloud_validation = False

    # Instantiate local SLOS simulator
    local_sim = SimulatorFactory.build(circuit)

    history: List[Dict[str, Any]] = []
    cloud_job_ids: List[str] = []
    best_overall_energy = float("inf")
    best_overall_conformation: Tuple[int, ...] = ()

    # 5. Define objective function
    iteration_counter = 0

    def objective_fn(phases: np.ndarray) -> float:
        nonlocal iteration_counter, best_overall_energy, best_overall_conformation
        iteration_counter += 1

        if iteration_counter > max_iter:
            return history[-1]["cvar_cost"] if history else 0.0

        # Bind parameters
        for p, val in zip(params, phases):
            p.set_value(float(val))

        # Query distribution
        if is_cloud_loop:
            res = backend.run_circuit(circuit, input_state, shots=shots)
            raw_dist = res.get("results", {})
            if "job_id" in res and res["job_id"]:
                cloud_job_ids.append(str(res["job_id"]))
        else:
            # Local exact simulation: mandatory re-set circuit to invalidate cached unitary
            local_sim.set_circuit(circuit)
            bsd = local_sim.probs(input_state)
            raw_dist = dict(bsd.items())

        # Decode distribution to conformations
        conf_dist = decode_distribution(
            raw_dist,
            num_turns=num_turns,
            bits_per_turn=effective_bpt,
            post_select_strict=False,
        )

        # Compute CVaR cost
        cvar_val, step_best_e, step_best_c = compute_cvar(conf_dist, energy_fn, alpha=alpha)

        if step_best_e < best_overall_energy:
            best_overall_energy = step_best_e
            best_overall_conformation = step_best_c

        history.append({
            "iteration": iteration_counter,
            "cvar_cost": cvar_val,
            "best_energy": best_overall_energy,
            "best_conformation": best_overall_conformation,
            "phases": phases.copy(),
        })
        return cvar_val

    # 6. Run Classical Optimizer
    init_x = np.array([float(p) for p in params], dtype=float)
    if np.allclose(init_x, 0.0):
        init_x += np.random.uniform(-0.05, 0.05, size=len(init_x))

    opt_method = method.upper()
    max_evals = max(max_iter, len(params) + 15) if opt_method == "COBYLA" else max_iter
    opt_options: Dict[str, Any] = {"maxiter": max_evals}

    logger.debug(
        "Starting VQE optimization (%s) for %d turns, %d modes, %d params. Mode: %s",
        opt_method,
        num_turns,
        num_modes,
        len(params),
        execution_mode,
    )

    opt_result = opt.minimize(
        objective_fn,
        init_x,
        method=opt_method,
        options=opt_options,
    )

    best_phases = opt_result.x

    # 7. Apply optimal phases and compute final distribution
    for p, val in zip(params, best_phases):
        p.set_value(float(val))

    local_sim.set_circuit(circuit)
    final_raw = dict(local_sim.probs(input_state).items())
    final_conf_dist = decode_distribution(
        final_raw,
        num_turns=num_turns,
        bits_per_turn=effective_bpt,
        post_select_strict=False,
    )

    # 8. Run final cloud validation shot if requested and available
    if run_cloud_validation and has_remote_backend:
        try:
            logger.info("Executing final validation circuit on Quandela Cloud...")
            cloud_res = backend.run_circuit(circuit, input_state, shots=shots)
            if "job_id" in cloud_res and cloud_res["job_id"]:
                cloud_job_ids.append(str(cloud_res["job_id"]))
            cloud_raw = cloud_res.get("results", {})
            if cloud_raw:
                cloud_conf_dist = decode_distribution(
                    cloud_raw,
                    num_turns=num_turns,
                    bits_per_turn=effective_bpt,
                    post_select_strict=False,
                )
                final_conf_dist = cloud_conf_dist
        except Exception as e:
            logger.warning("Cloud validation job failed (%s). Retaining local simulation results.", e)

    elapsed = time.time() - start_time
    logger.debug(
        "VQE finished in %.2fs. Best energy: %.4f, Conformation: %s, Cloud Jobs: %d",
        elapsed,
        best_overall_energy,
        best_overall_conformation,
        len(cloud_job_ids),
    )

    return VQEResult(
        best_energy=float(best_overall_energy),
        best_conformation=best_overall_conformation,
        best_phases=best_phases,
        history=history,
        cloud_job_ids=cloud_job_ids,
        optimal_distribution=final_conf_dist,
        num_iterations=iteration_counter,
        converged=bool(opt_result.success or iteration_counter >= max_iter),
        execution_time_seconds=elapsed,
    )


def evaluate_vqe(
    protein: Any,
    backend: Optional[Any] = None,
    max_iter: int = 50,
    model: str = "MJ",
    **kwargs: Any,
) -> VQEResult:
    """
    Interface wrapper for PROJECT.md and SCOPE.md contract compliance.
    """
    return run_vqe(protein=protein, backend=backend, max_iter=max_iter, model=model, **kwargs)
