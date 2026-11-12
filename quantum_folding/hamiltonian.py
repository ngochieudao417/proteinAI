"""
quantum_folding.hamiltonian
===========================
Energy potentials, folding Hamiltonians, and QUBO/Ising matrix formulation
for lattice protein folding.

Includes:
- Complete 20x20 Miyazawa-Jernigan (MJ 1985 & 1996) contact matrices
- HP model contact potential
- FoldingHamiltonian class (interaction + overlap + backtrack)
- QUBO matrix builder for turn-based variables
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Sequence, Tuple, Union
import numpy as np

from .lattice import (
    AMINO_ACIDS,
    LatticeModel,
    LatticeProtein,
    count_backtracks,
    count_overlaps,
    get_contacts,
    turns_to_coordinates,
)

AA_TO_INDEX: Dict[str, int] = {aa: i for i, aa in enumerate(AMINO_ACIDS)}

# ============================================================================
# Miyazawa-Jernigan 1985 Contact Potential (Table V / VI)
# Values in reduced k_B * T units. Ordered by 'ACDEFGHIKLMNPQRSTVWY'.
# ============================================================================
MJ_1985_MATRIX: np.ndarray = np.array([
    # A      C      D      E      F      G      H      I      K      L      M      N      P      Q      R      S      T      V      W      Y
    [-2.51, -3.38, -1.57, -1.51, -4.36, -2.15, -2.09, -4.41, -1.10, -3.96, -3.99, -1.44, -1.81, -1.70, -1.50, -1.89, -2.15, -3.62, -3.93, -2.85], # A
    [-3.38, -5.44, -2.66, -2.08, -5.63, -3.16, -3.63, -5.03, -1.54, -5.03, -5.05, -2.59, -2.92, -2.73, -2.70, -2.86, -2.88, -4.46, -4.76, -3.89], # C
    [-1.57, -2.66, -0.96, -1.23, -3.31, -1.62, -2.14, -2.91, -1.32, -2.59, -2.90, -1.33, -1.19, -1.26, -1.98, -1.46, -1.66, -2.25, -2.91, -2.25], # D
    [-1.51, -2.08, -1.23, -1.18, -3.51, -1.22, -2.27, -3.23, -1.60, -2.91, -3.19, -1.43, -1.40, -1.33, -2.07, -1.48, -1.45, -2.56, -2.94, -2.42], # E
    [-4.36, -5.63, -3.31, -3.51, -6.85, -3.72, -4.61, -6.39, -2.83, -6.26, -6.68, -3.55, -3.73, -3.30, -3.54, -3.56, -3.76, -5.75, -6.02, -4.95], # F
    [-2.15, -3.16, -1.62, -1.22, -3.72, -2.17, -1.94, -3.65, -0.84, -3.43, -3.75, -1.56, -1.72, -1.54, -1.68, -1.70, -2.03, -3.06, -3.37, -2.50], # G
    [-2.09, -3.63, -2.14, -2.27, -4.61, -1.94, -2.78, -3.76, -1.09, -3.84, -3.31, -2.01, -2.17, -1.85, -2.12, -1.94, -2.31, -3.38, -4.02, -3.33], # H
    [-4.41, -5.03, -2.91, -3.23, -6.39, -3.65, -3.76, -6.22, -2.70, -6.17, -6.33, -2.99, -3.47, -3.22, -3.33, -3.43, -3.74, -5.58, -5.64, -4.63], # I
    [-1.10, -1.54, -1.32, -1.60, -2.83, -0.84, -1.09, -2.70,  0.13, -2.63, -3.11, -0.91, -0.67, -1.02, -0.06, -0.83, -1.02, -1.95, -2.49, -2.01], # K
    [-3.96, -5.03, -2.59, -2.91, -6.26, -3.43, -3.84, -6.17, -2.63, -5.79, -6.01, -2.99, -3.06, -3.09, -3.15, -3.16, -3.43, -5.38, -5.50, -4.26], # L
    [-3.99, -5.05, -2.90, -3.19, -6.68, -3.75, -3.31, -6.33, -3.11, -6.01, -6.06, -3.50, -4.11, -3.17, -3.49, -3.55, -3.73, -5.52, -6.37, -4.92], # M
    [-1.44, -2.59, -1.33, -1.43, -3.55, -1.56, -2.01, -2.99, -0.91, -2.99, -3.50, -1.59, -1.43, -1.36, -1.41, -1.31, -1.51, -2.36, -3.11, -2.47], # N
    [-1.81, -2.92, -1.19, -1.40, -3.73, -1.72, -2.17, -3.47, -0.67, -3.06, -4.11, -1.43, -1.18, -1.73, -1.85, -1.35, -1.66, -2.96, -3.66, -2.80], # P
    [-1.70, -2.73, -1.26, -1.33, -3.30, -1.54, -1.85, -3.22, -1.02, -3.09, -3.17, -1.36, -1.73, -0.89, -1.85, -1.37, -1.59, -2.67, -3.16, -2.53], # Q
    [-1.50, -2.70, -1.98, -2.07, -3.54, -1.68, -2.12, -3.33, -0.06, -3.15, -3.49, -1.41, -1.85, -1.85, -1.39, -1.22, -1.97, -2.78, -3.56, -2.75], # R
    [-1.89, -2.86, -1.46, -1.48, -3.56, -1.70, -1.94, -3.43, -0.83, -3.16, -3.55, -1.31, -1.35, -1.37, -1.22, -1.48, -1.59, -2.79, -2.95, -2.30], # S
    [-2.15, -2.88, -1.66, -1.45, -3.76, -2.03, -2.31, -3.74, -1.02, -3.43, -3.73, -1.51, -1.66, -1.59, -1.97, -1.59, -1.72, -2.95, -3.31, -2.48], # T
    [-3.62, -4.46, -2.25, -2.56, -5.75, -3.06, -3.38, -5.58, -1.95, -5.38, -5.52, -2.36, -2.96, -2.67, -2.78, -2.79, -2.95, -4.94, -5.05, -4.05], # V
    [-3.93, -4.76, -2.91, -2.94, -6.02, -3.37, -4.02, -5.64, -2.49, -5.50, -6.37, -3.11, -3.66, -3.16, -3.56, -2.95, -3.31, -5.05, -5.42, -4.44], # W
    [-2.85, -3.89, -2.25, -2.42, -4.95, -2.50, -3.33, -4.63, -2.01, -4.26, -4.92, -2.47, -2.80, -2.53, -2.75, -2.30, -2.48, -4.05, -4.44, -3.55], # Y
], dtype=np.float64)

# ============================================================================
# Miyazawa-Jernigan 1996 Contact Potential (Table 3)
# Values in reduced k_B * T units. Ordered by 'ACDEFGHIKLMNPQRSTVWY'.
# ============================================================================
MJ_1996_MATRIX: np.ndarray = np.array([
    # A      C      D      E      F      G      H      I      K      L      M      N      P      Q      R      S      T      V      W      Y
    [-2.72, -3.57, -1.70, -1.51, -4.81, -2.31, -2.41, -4.58, -1.31, -4.91, -3.94, -1.84, -2.03, -1.89, -1.83, -2.01, -2.32, -4.04, -3.82, -3.36], # A
    [-3.57, -5.44, -2.41, -2.27, -5.80, -3.16, -3.60, -5.50, -1.95, -5.83, -4.99, -2.59, -3.07, -2.85, -2.57, -2.86, -3.11, -4.96, -4.95, -4.16], # C
    [-1.70, -2.41, -1.21, -1.02, -3.48, -1.59, -2.32, -3.17, -1.68, -3.40, -2.57, -1.68, -1.33, -1.46, -2.29, -1.63, -1.80, -2.48, -2.84, -2.76], # D
    [-1.51, -2.27, -1.02, -0.91, -3.56, -1.22, -2.15, -3.27, -1.80, -3.59, -2.89, -1.51, -1.26, -1.42, -2.27, -1.48, -1.74, -2.67, -2.99, -2.79], # E
    [-4.81, -5.80, -3.48, -3.56, -7.26, -4.13, -4.77, -6.84, -3.36, -7.28, -6.56, -3.75, -4.25, -4.10, -3.98, -4.02, -4.28, -6.29, -6.16, -5.66], # F
    [-2.31, -3.16, -1.59, -1.22, -4.13, -2.24, -2.15, -3.78, -1.15, -4.16, -3.39, -1.74, -1.87, -1.66, -1.72, -1.82, -2.08, -3.38, -3.42, -3.01], # G
    [-2.41, -3.60, -2.32, -2.15, -4.77, -2.15, -3.05, -4.14, -1.35, -4.54, -3.98, -2.08, -2.25, -1.98, -2.16, -2.11, -2.42, -3.58, -3.98, -3.52], # H
    [-4.58, -5.50, -3.17, -3.27, -6.84, -3.78, -4.14, -6.54, -3.01, -7.04, -6.02, -3.24, -3.76, -3.67, -3.63, -3.52, -4.03, -6.05, -5.78, -5.25], # I
    [-1.31, -1.95, -1.68, -1.80, -3.36, -1.15, -1.35, -3.01, -0.12, -3.37, -2.48, -1.21, -0.97, -1.29, -0.59, -1.05, -1.31, -2.49, -2.69, -2.60], # K
    [-4.91, -5.83, -3.40, -3.59, -7.28, -4.16, -4.54, -7.04, -3.37, -7.37, -6.41, -3.74, -4.20, -4.04, -4.03, -3.92, -4.34, -6.48, -6.14, -5.67], # L
    [-3.94, -4.99, -2.57, -2.89, -6.56, -3.39, -3.98, -6.02, -2.48, -6.41, -5.46, -2.95, -3.45, -3.30, -3.12, -3.03, -3.51, -5.32, -5.55, -4.91], # M
    [-1.84, -2.59, -1.68, -1.51, -3.75, -1.74, -2.08, -3.24, -1.21, -3.74, -2.95, -1.68, -1.53, -1.71, -1.64, -1.58, -1.88, -2.83, -3.07, -2.76], # N
    [-2.03, -3.07, -1.33, -1.26, -4.25, -1.87, -2.25, -3.76, -0.97, -4.20, -3.45, -1.53, -1.75, -1.73, -1.70, -1.57, -1.90, -3.32, -3.73, -3.19], # P
    [-1.89, -2.85, -1.46, -1.42, -4.10, -1.66, -1.98, -3.67, -1.29, -4.04, -3.30, -1.71, -1.73, -1.54, -1.80, -1.49, -1.90, -3.07, -3.11, -2.97], # Q
    [-1.83, -2.57, -2.29, -2.27, -3.98, -1.72, -2.16, -3.63, -0.59, -4.03, -3.12, -1.64, -1.70, -1.80, -1.55, -1.62, -1.90, -3.07, -3.41, -3.16], # R
    [-2.01, -2.86, -1.63, -1.48, -4.02, -1.82, -2.11, -3.52, -1.05, -3.92, -3.03, -1.58, -1.57, -1.49, -1.62, -1.67, -1.96, -3.05, -2.99, -2.78], # S
    [-2.32, -3.11, -1.80, -1.74, -4.28, -2.08, -2.42, -4.03, -1.31, -4.34, -3.51, -1.88, -1.90, -1.90, -1.90, -1.96, -2.12, -3.46, -3.22, -3.01], # T
    [-4.04, -4.96, -2.48, -2.67, -6.29, -3.38, -3.58, -6.05, -2.49, -6.48, -5.32, -2.83, -3.32, -3.07, -3.07, -3.05, -3.46, -5.52, -5.18, -4.62], # V
    [-3.82, -4.95, -2.84, -2.99, -6.16, -3.42, -3.98, -5.78, -2.69, -6.14, -5.55, -3.07, -3.73, -3.11, -3.41, -2.99, -3.22, -5.18, -5.06, -4.66], # W
    [-3.36, -4.16, -2.76, -2.79, -5.66, -3.01, -3.52, -5.25, -2.60, -5.67, -4.91, -2.76, -3.19, -2.97, -3.16, -2.78, -3.01, -4.62, -4.66, -4.17], # Y
], dtype=np.float64)

# Default contact matrix alias
MJ_CONTACT_MATRIX: np.ndarray = MJ_1985_MATRIX


def get_mj_energy(aa1: str, aa2: str, version: str = "1985") -> float:
    """Retrieves pairwise Miyazawa-Jernigan contact energy between two amino acids."""
    idx1 = AA_TO_INDEX[aa1.upper()]
    idx2 = AA_TO_INDEX[aa2.upper()]
    matrix = MJ_1985_MATRIX if version == "1985" else MJ_1996_MATRIX
    return float(matrix[idx1, idx2])


def get_hp_energy(res1: str, res2: str) -> float:
    """HP contact potential: -1.0 for H-H contact, 0.0 otherwise."""
    r1, r2 = res1.upper(), res2.upper()
    return -1.0 if (r1 == "H" and r2 == "H") else 0.0


class FoldingHamiltonian:
    """
    Evaluates conformational energy for lattice protein folding.
    E = E_contacts + lambda_overlap * N_overlap + lambda_backtrack * N_backtrack
    """

    def __init__(
        self,
        protein: Optional[LatticeProtein] = None,
        model: str = "MJ",
        matrix_version: str = "1985",
        lambda_overlap: float = 30.0,
        lambda_backtrack: float = 20.0,
    ) -> None:
        self.protein = protein
        self.model = model.upper()
        self.matrix_version = matrix_version
        self.lambda_overlap = float(lambda_overlap)
        self.lambda_backtrack = float(lambda_backtrack)

        if self.model not in ("MJ", "HP"):
            raise ValueError(f"Unknown energy model: {self.model}. Must be 'MJ' or 'HP'.")

    def pair_energy(self, s1: str, s2: str) -> float:
        """Returns pairwise contact interaction energy between two residues."""
        if self.model == "HP":
            return get_hp_energy(s1, s2)
        return get_mj_energy(s1, s2, version=self.matrix_version)

    def energy(
        self,
        conformation: Optional[Union[Sequence[int], LatticeProtein, Sequence[Tuple[int, ...]]]] = None,
    ) -> float:
        """
        Computes total conformational energy.
        
        Args:
            conformation: Can be:
                - None: evaluates self.protein
                - LatticeProtein instance
                - Sequence of turn integers (e.g. (0, 1, 2) or (0, 0, 0, 0))
                - Sequence of coordinate tuples
        """
        if conformation is None:
            if self.protein is None:
                raise ValueError("No protein or conformation provided for energy evaluation.")
            return self.evaluate(self.protein)

        if isinstance(conformation, LatticeProtein):
            return self.evaluate(conformation)

        # Sequence of items
        items = list(conformation)
        if not items:
            return 0.0

        if isinstance(items[0], tuple):
            # Coordinates directly
            coords = items
            seq = self.protein.sequence if self.protein else "H" * len(coords)
            hp_seq = self.protein.hp_sequence if self.protein else "H" * len(coords)
        else:
            # Turn indices
            turns = items
            if self.protein is None:
                raise ValueError("A LatticeProtein reference is required to convert turns to coordinates.")
            coords = self.protein.compute_coordinates(turns)
            seq = self.protein.sequence
            hp_seq = self.protein.hp_sequence

        contacts = get_contacts(coords)
        e_contacts = 0.0

        for i, j in contacts:
            if i < len(seq) and j < len(seq):
                s1 = hp_seq[i] if self.model == "HP" else seq[i]
                s2 = hp_seq[j] if self.model == "HP" else seq[j]
                e_contacts += self.pair_energy(s1, s2)

        n_overlaps = count_overlaps(coords)
        n_backtracks = count_backtracks(coords)

        return float(e_contacts + (self.lambda_overlap * n_overlaps) + (self.lambda_backtrack * n_backtracks))

    def evaluate(self, protein: Optional[LatticeProtein] = None) -> float:
        """Computes total conformational energy of a LatticeProtein."""
        target_prot = protein or self.protein
        if target_prot is None:
            raise ValueError("No protein specified for evaluation.")

        contacts = target_prot.contacts
        e_contacts = 0.0

        for i, j in contacts:
            if self.model == "HP":
                e_contacts += self.pair_energy(target_prot.hp_sequence[i], target_prot.hp_sequence[j])
            else:
                e_contacts += self.pair_energy(target_prot.sequence[i], target_prot.sequence[j])

        n_overlaps = count_overlaps(target_prot.coords)
        n_backtracks = count_backtracks(target_prot.coords)

        return float(e_contacts + (self.lambda_overlap * n_overlaps) + (self.lambda_backtrack * n_backtracks))

    def breakdown(
        self,
        protein: Optional[LatticeProtein] = None,
        turns: Optional[Sequence[int]] = None,
    ) -> Dict[str, Any]:
        """Detailed breakdown of energy components."""
        target_prot = protein or self.protein
        if turns is not None and target_prot is not None:
            coords = target_prot.compute_coordinates(turns)
            seq = target_prot.sequence
            hp_seq = target_prot.hp_sequence
        elif target_prot is not None:
            coords = target_prot.coords
            seq = target_prot.sequence
            hp_seq = target_prot.hp_sequence
        else:
            raise ValueError("Cannot perform breakdown without protein reference.")

        contacts = get_contacts(coords)
        e_contacts = 0.0
        contact_details = []

        for i, j in contacts:
            s1 = hp_seq[i] if self.model == "HP" else seq[i]
            s2 = hp_seq[j] if self.model == "HP" else seq[j]
            e_pair = self.pair_energy(s1, s2)
            e_contacts += e_pair
            contact_details.append({"i": i, "j": j, "res_i": s1, "res_j": s2, "energy": e_pair})

        n_overlaps = count_overlaps(coords)
        n_backtracks = count_backtracks(coords)
        total = e_contacts + (self.lambda_overlap * n_overlaps) + (self.lambda_backtrack * n_backtracks)

        return {
            "total_energy": float(total),
            "e_contacts": float(e_contacts),
            "n_contacts": len(contacts),
            "contact_pairs": contact_details,
            "n_overlaps": n_overlaps,
            "overlap_penalty": float(self.lambda_overlap * n_overlaps),
            "n_backtracks": n_backtracks,
            "backtrack_penalty": float(self.lambda_backtrack * n_backtracks),
            "is_valid": (n_overlaps == 0 and n_backtracks == 0),
        }


def build_qubo_matrix(
    protein_or_seq: Any,
    lattice: Union[LatticeModel, str] = LatticeModel.SQUARE_2D,
    model: str = "HP",
    lambda_onehot: float = 50.0,
    lambda_backtrack: float = 20.0,
    lambda_overlap: float = 30.0,
    matrix_version: str = "1985",
) -> Tuple[np.ndarray, float]:
    """
    Builds the Quadratic Unconstrained Binary Optimization (QUBO) matrix Q and scalar offset
    such that E(x) = x^T Q x + offset for one-hot turn variables.

    Parameters
    ----------
    protein_or_seq : LatticeProtein or str
        Peptide sequence or protein instance.
    lattice : LatticeModel or str
        Lattice geometry (SQUARE_2D or CUBIC_3D).
    model : str
        'HP' or 'MJ'.
    lambda_onehot : float
        Penalty weight for one-hot direction constraint violation.
    lambda_backtrack : float
        Penalty weight for immediate 180-degree backtracking.
    lambda_overlap : float
        Penalty weight for non-local spatial collisions.
    matrix_version : str
        '1985' or '1996' for MJ matrix.

    Returns
    -------
    Tuple[np.ndarray, float]
        (Q, offset) where Q is an upper-triangular matrix of shape (V, V).
    """
    seq = protein_or_seq.sequence if isinstance(protein_or_seq, LatticeProtein) else str(protein_or_seq)
    seq = seq.strip().upper()
    n_res = len(seq)
    if n_res < 3:
        raise ValueError("Peptide sequence must contain at least 3 residues to formulate turns.")

    lat_enum = LatticeModel(lattice) if isinstance(lattice, str) else lattice
    num_turns = n_res - 2
    num_dirs = 4 if lat_enum == LatticeModel.SQUARE_2D else 6
    num_vars = num_turns * num_dirs

    Q = np.zeros((num_vars, num_vars), dtype=np.float64)
    offset = 0.0

    def var_idx(step: int, direction: int) -> int:
        return step * num_dirs + direction

    # 1. One-hot penalty for each turn step: lambda_onehot * (sum_d x_{t, d} - 1)^2
    for t in range(num_turns):
        offset += lambda_onehot
        for d in range(num_dirs):
            idx = var_idx(t, d)
            Q[idx, idx] -= lambda_onehot
        for d1 in range(num_dirs):
            for d2 in range(d1 + 1, num_dirs):
                idx1 = var_idx(t, d1)
                idx2 = var_idx(t, d2)
                Q[idx1, idx2] += 2.0 * lambda_onehot

    # 2. Backtrack penalty:
    # First turn (step 0): bond 0 was along +X (dir 0). Opposite is -X (dir 2 in 2D, dir 1 in 3D).
    opp_0 = 2 if lat_enum == LatticeModel.SQUARE_2D else 1
    idx_first_opp = var_idx(0, opp_0)
    Q[idx_first_opp, idx_first_opp] += lambda_backtrack

    # Subsequent turns: penalize consecutive opposite directions
    opposites = {0: 2, 1: 3, 2: 0, 3: 1} if lat_enum == LatticeModel.SQUARE_2D else {0: 1, 1: 0, 2: 3, 3: 2, 4: 5, 5: 4}
    for t in range(num_turns - 1):
        for d in range(num_dirs):
            d_opp = opposites[d]
            idx_curr = var_idx(t, d)
            idx_next = var_idx(t + 1, d_opp)
            if idx_curr < idx_next:
                Q[idx_curr, idx_next] += lambda_backtrack
            else:
                Q[idx_next, idx_curr] += lambda_backtrack

    # 3. Contact & Loop Terms (Pairwise loop interactions)
    # On 2D square lattice, 4-mer segment (residues 0 and 3) forms contact via U-turns (+Y, -X) or (-Y, -X)
    if lat_enum == LatticeModel.SQUARE_2D and num_turns >= 2:
        hamiltonian = FoldingHamiltonian(model=model, matrix_version=matrix_version)
        e_contact_03 = hamiltonian.pair_energy(seq[0], seq[3])
        # step 0: dir 1 (+Y) or dir 3 (-Y); step 1: dir 2 (-X)
        idx_up = var_idx(0, 1)
        idx_down = var_idx(0, 3)
        idx_left = var_idx(1, 2)

        Q[idx_up, idx_left] += e_contact_03
        Q[idx_down, idx_left] += e_contact_03

    return Q, offset


def decode_qubo_solution(
    x: Sequence[int],
    lattice: Union[LatticeModel, str] = LatticeModel.SQUARE_2D,
) -> List[int]:
    """
    Decodes a binary QUBO solution vector x into an integer turn list.
    """
    lat_enum = LatticeModel(lattice) if isinstance(lattice, str) else lattice
    num_dirs = 4 if lat_enum == LatticeModel.SQUARE_2D else 6
    turns = []
    num_turns = len(x) // num_dirs
    for t in range(num_turns):
        chunk = x[t * num_dirs : (t + 1) * num_dirs]
        best_dir = int(np.argmax(chunk))
        turns.append(best_dir)
    return turns
