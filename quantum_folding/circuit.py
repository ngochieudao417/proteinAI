"""
quantum_folding.circuit
=======================
Photonic Circuit Construction and Dual-Rail Conformation Mapping.
Perceval SDK 1.2.4 compatible.
"""

from __future__ import annotations

import logging
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np
import perceval as pcvl
from perceval.components import BS, PS, Circuit, GenericInterferometer
from perceval.utils import BasicState, Parameter
from perceval.utils._enums import InterferometerShape

logger = logging.getLogger(__name__)


def build_ansatz_circuit(
    num_modes: int,
    num_photons: Optional[int] = None,
    design: str = "clements",
    shape: Optional[str] = None,
    phase_shifter_layer: bool = True,
    param_prefix: str = "phi",
    initial_phases: Optional[Sequence[float]] = None,
) -> Tuple[Circuit, List[Parameter]]:
    """
    Constructs a universal parameterized linear optical circuit using Perceval.

    Implements Clements (rectangular) or Reck (triangular) decomposition with tunable
    Mach-Zehnder Interferometers (MZIs) and an optional output phase-shifter layer.

    All parameters are automatically initialized to avoid 'AssertionError: All parameters
    must be defined to compute numeric unitary matrix'.

    :param num_modes: Number of optical spatial modes (m >= 2).
    :param num_photons: (Optional) Expected photon count.
    :param design: 'clements' (rectangular mesh, default) or 'reck' (triangular mesh).
    :param shape: Alternative keyword for design ('rectangle' or 'triangle').
    :param phase_shifter_layer: If True, adds an output phase shifter layer for U(m) universality.
    :param param_prefix: Prefix string for named parameters.
    :param initial_phases: Optional initial numerical values for parameters.
    :return: Tuple of (pcvl.Circuit, List[pcvl.Parameter]).
    """
    if num_modes < 2:
        raise ValueError(f"num_modes must be at least 2, got {num_modes}")

    chosen_design = shape if shape is not None else design
    norm_design = chosen_design.lower().strip()

    if norm_design in ["clements", "rectangle", "rect"]:
        inter_shape = InterferometerShape.RECTANGLE
    elif norm_design in ["reck", "triangle", "tri"]:
        inter_shape = InterferometerShape.TRIANGLE
    else:
        raise ValueError(f"Unsupported interferometer design: '{chosen_design}'. Choose 'clements' or 'reck'.")

    # Define MZI unit cell generator
    def mzi_cell_gen(idx: int) -> Circuit:
        p = pcvl.P(f"{param_prefix}_{idx}")
        return Circuit(2, name=f"MZI_{idx}") // BS() // (0, PS(p)) // BS()

    # Define optional output phase shifter generator
    ps_layer_gen = (lambda i: PS(pcvl.P(f"{param_prefix}_out_{i}"))) if phase_shifter_layer else None

    # Instantiate generic interferometer
    interferometer = GenericInterferometer(
        num_modes,
        mzi_cell_gen,
        shape=inter_shape,
        phase_shifter_fun_gen=ps_layer_gen,
        phase_at_output=True,
    )

    params = interferometer.get_parameters()

    # Crucial: Initialize all parameters to prevent undefined parameter assertion
    if initial_phases is not None:
        if len(initial_phases) != len(params):
            raise ValueError(
                f"Mismatch in initial_phases length: expected {len(params)}, got {len(initial_phases)}"
            )
        for p, val in zip(params, initial_phases):
            p.set_value(float(val))
    else:
        for p in params:
            p.set_value(0.0)

    logger.debug(
        "Built %s ansatz with %d modes, %d parameters. Output phase layer: %s",
        chosen_design,
        num_modes,
        len(params),
        phase_shifter_layer,
    )
    return interferometer, params


# Alias for SCOPE.md and PROJECT.md contract compliance
build_photonic_ansatz = build_ansatz_circuit


def build_fock_state(
    num_modes: int,
    num_photons: int,
    pattern: str = "dual_rail",
) -> BasicState:
    """
    Generates standard input Fock states for linear optical quantum computation.

    :param num_modes: Total number of spatial modes.
    :param num_photons: Total number of indistinguishable photons.
    :param pattern:
        - 'dual_rail': |1, 0, 1, 0, ...> (1 photon per 2 modes, up to num_photons).
        - 'consecutive' / 'dense': |1, 1, ..., 0, 0> (first num_photons modes populated).
        - 'single_mode': |num_photons, 0, 0, ...>.
    :return: pcvl.BasicState instance.
    """
    if num_photons < 1:
        raise ValueError(f"num_photons must be at least 1, got {num_photons}")
    if num_photons > num_modes and pattern in ["dual_rail", "consecutive"]:
        raise ValueError(
            f"num_photons ({num_photons}) exceeds num_modes ({num_modes}) for pattern '{pattern}'"
        )

    state_list = [0] * num_modes

    if pattern == "dual_rail":
        max_dual_rail_photons = num_modes // 2
        actual_photons = min(num_photons, max_dual_rail_photons)
        for i in range(actual_photons):
            state_list[2 * i] = 1
        for i in range(actual_photons, num_photons):
            state_list[2 * (i - actual_photons) + 1] = 1

    elif pattern in ["consecutive", "dense"]:
        for i in range(num_photons):
            state_list[i] = 1

    elif pattern == "single_mode":
        state_list[0] = num_photons

    else:
        raise ValueError(f"Unknown Fock state pattern: '{pattern}'. Choose 'dual_rail', 'consecutive', or 'single_mode'.")

    return BasicState(state_list)


def state_to_conformation(
    fock_state: Union[BasicState, Sequence[int]],
    num_turns: int,
    bits_per_turn: int = 2,
    strict: bool = False,
) -> Optional[Tuple[int, ...]]:
    """
    Decodes a detected optical Fock state into a sequence of lattice turns.

    Under dual-rail encoding, each decision bit is represented by 2 spatial modes.
    For 2D lattices, each turn typically uses 2 bits (4 directions: Up=0, Right=1, Down=2, Left=3).

    :param fock_state: Measured Fock state (BasicState or tuple/list of mode occupation numbers).
    :param num_turns: Number of turns to decode.
    :param bits_per_turn: Number of bits per turn (default: 2; adaptively falls to 1 if mode count matches).
    :param strict:
        If True: Requires strictly valid dual-rail states (n_0 + n_1 == 1 per qubit). Returns None on bunching.
        If False: Uses soft argmax projection (bit = 1 if n_1 > n_0 else 0). Always returns a valid conformation.
    :return: Tuple of turn indices (e.g. (0, 3)), or None if strict=True and state is invalid.
    """
    s = list(fock_state)
    if num_turns <= 0:
        return ()

    # Adapt bits_per_turn if modes strictly correspond to 1 bit per turn
    actual_bpt = bits_per_turn
    if len(s) < num_turns * bits_per_turn * 2 and len(s) >= num_turns * 2:
        actual_bpt = len(s) // (2 * num_turns)
    if actual_bpt < 1:
        actual_bpt = 1

    num_qubits = num_turns * actual_bpt
    num_required_modes = num_qubits * 2

    # Pad with zeros if state length is shorter than required modes
    if len(s) < num_required_modes:
        s = s + [0] * (num_required_modes - len(s))

    bits: List[int] = []
    for q in range(num_qubits):
        m0 = s[2 * q]
        m1 = s[2 * q + 1]

        if strict:
            if m0 == 1 and m1 == 0:
                bits.append(0)
            elif m0 == 0 and m1 == 1:
                bits.append(1)
            else:
                return None
        else:
            # Argmax projection: photon in odd mode indicates bit 1
            bits.append(1 if m1 > m0 else 0)

    # Group bits into turn integers
    turns: List[int] = []
    for t in range(num_turns):
        turn_val = 0
        for b in range(actual_bpt):
            bit_idx = t * actual_bpt + b
            turn_val = (turn_val << 1) | bits[bit_idx]
        turns.append(turn_val)

    return tuple(turns)


def decode_distribution(
    raw_distribution: Dict[Any, Union[float, int]],
    num_turns: int,
    bits_per_turn: int = 2,
    post_select_strict: bool = False,
) -> Dict[Tuple[int, ...], float]:
    """
    Transforms a raw Fock-state probability or count distribution into a normalized
    conformation probability distribution.

    :param raw_distribution: Mapping from FockState/BasicState to float (probability) or int (count).
    :param num_turns: Number of turns per conformation.
    :param bits_per_turn: Number of bits per turn (default 2).
    :param post_select_strict: If True, filters out states with photon bunching.
                               Automatically falls back to non-strict if all states are filtered.
    :return: Dict mapping conformation tuple (e.g. (0, 3)) to normalized probability summing to 1.0.
    """
    if not raw_distribution:
        return {}

    conf_probs: Dict[Tuple[int, ...], float] = {}
    total_valid_weight = 0.0

    # Primary pass: decode states
    for state, weight in raw_distribution.items():
        w = float(weight)
        if w <= 0.0:
            continue

        conf = state_to_conformation(
            state,
            num_turns=num_turns,
            bits_per_turn=bits_per_turn,
            strict=post_select_strict,
        )
        if conf is not None:
            conf_probs[conf] = conf_probs.get(conf, 0.0) + w
            total_valid_weight += w

    # Fallback: if post_select_strict discarded 100% of counts, use argmax
    if total_valid_weight == 0.0 and post_select_strict:
        logger.warning("Strict post-selection discarded all states. Falling back to argmax projection.")
        for state, weight in raw_distribution.items():
            w = float(weight)
            if w <= 0.0:
                continue
            conf = state_to_conformation(
                state,
                num_turns=num_turns,
                bits_per_turn=bits_per_turn,
                strict=False,
            )
            if conf is not None:
                conf_probs[conf] = conf_probs.get(conf, 0.0) + w
                total_valid_weight += w

    # Normalize distribution
    if total_valid_weight > 0.0:
        for conf in conf_probs:
            conf_probs[conf] /= total_valid_weight

    return conf_probs
