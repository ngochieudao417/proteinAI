"""
quantum_folding.lattice
=======================
Lattice representations, HP encoding, coordinate generation, and spatial
verification metrics for quantum protein folding simulations.

Supported geometries:
- 2D Square Lattice (4 directions)
- 3D Cubic Lattice (6 directions)
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum, IntEnum
from typing import Any, Dict, List, Optional, Sequence, Set, Tuple, Union
import numpy as np

# Standard IUPAC 20 Canonical Amino Acids
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"

# Prompt-mandated HP Classification
# Hydrophobic (H): A, C, F, I, L, M, V, W, Y (9 residues)
# Polar (P): D, E, G, H, K, N, P, Q, R, S, T (11 residues)
HYDROPHOBIC_RESIDUES: frozenset[str] = frozenset("ACFILMVWY")
POLAR_RESIDUES: frozenset[str] = frozenset("DEGHKNPQRST")

assert set(HYDROPHOBIC_RESIDUES) | set(POLAR_RESIDUES) == set(AMINO_ACIDS), "HP sets must span all 20 AAs"
assert set(HYDROPHOBIC_RESIDUES) & set(POLAR_RESIDUES) == set(), "HP sets must be disjoint"


class LatticeModel(str, Enum):
    """Supported lattice models."""
    SQUARE_2D = "2D"
    CUBIC_3D = "3D"


class Direction2D(IntEnum):
    """
    2D Square Lattice Absolute Directions.
    Right (+X), Up (+Y), Left (-X), Down (-Y).
    """
    RIGHT = 0
    UP = 1
    LEFT = 2
    DOWN = 3


class Direction3D(IntEnum):
    """
    3D Cubic Lattice Absolute Directions.
    +X, -X, +Y, -Y, +Z, -Z.
    """
    POS_X = 0
    NEG_X = 1
    POS_Y = 2
    NEG_Y = 3
    POS_Z = 4
    NEG_Z = 5


# Unit spatial displacement vectors
DIR_2D: List[Tuple[int, int]] = [(1, 0), (0, 1), (-1, 0), (0, -1)]
DIR_3D: List[Tuple[int, int, int]] = [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]

VECTORS_2D: Dict[int, Tuple[int, int]] = {
    Direction2D.RIGHT: (1, 0),
    Direction2D.UP: (0, 1),
    Direction2D.LEFT: (-1, 0),
    Direction2D.DOWN: (0, -1),
}

VECTORS_3D: Dict[int, Tuple[int, int, int]] = {
    Direction3D.POS_X: (1, 0, 0),
    Direction3D.NEG_X: (-1, 0, 0),
    Direction3D.POS_Y: (0, 1, 0),
    Direction3D.NEG_Y: (0, -1, 0),
    Direction3D.POS_Z: (0, 0, 1),
    Direction3D.NEG_Z: (0, 0, -1),
}

# Opposite direction mappings for backtrack detection
OPPOSITES_2D: Dict[int, int] = {
    Direction2D.RIGHT: Direction2D.LEFT,
    Direction2D.UP: Direction2D.DOWN,
    Direction2D.LEFT: Direction2D.RIGHT,
    Direction2D.DOWN: Direction2D.UP,
}

OPPOSITES_3D: Dict[int, int] = {
    Direction3D.POS_X: Direction3D.NEG_X,
    Direction3D.NEG_X: Direction3D.POS_X,
    Direction3D.POS_Y: Direction3D.NEG_Y,
    Direction3D.NEG_Y: Direction3D.POS_Y,
    Direction3D.POS_Z: Direction3D.NEG_Z,
    Direction3D.NEG_Z: Direction3D.POS_Z,
}


def encode_hp(sequence: str, is_already_hp: bool = False) -> str:
    """
    Encodes an amino acid sequence into Hydrophobic (H) and Polar (P) symbols.

    Parameters
    ----------
    sequence : str
        Peptide sequence (e.g. 'ACDEFGHIKLMNPQRSTVWY' or 'HPPH').
    is_already_hp : bool
        If True, validates and returns the string as HP directly.

    Returns
    -------
    str
        String containing only 'H' and 'P'.

    Raises
    ------
    ValueError
        If sequence contains invalid characters or is empty.
    """
    cleaned = sequence.strip().upper()
    if not cleaned:
        raise ValueError("Sequence cannot be empty.")

    if is_already_hp or set(cleaned).issubset({"H", "P"}):
        if not set(cleaned).issubset({"H", "P"}):
            raise ValueError(f"Invalid HP sequence: {cleaned}")
        return cleaned

    invalid = set(cleaned) - set(AMINO_ACIDS)
    if invalid:
        raise ValueError(f"Invalid amino acid characters in sequence: {invalid}")

    return "".join("H" if aa in HYDROPHOBIC_RESIDUES else "P" for aa in cleaned)


def turns_to_coordinates(
    turns: Sequence[int],
    lattice: Union[LatticeModel, str] = LatticeModel.SQUARE_2D,
    relative: bool = False,
    fix_first_step: bool = False,
    origin: Optional[Tuple[int, ...]] = None,
) -> List[Tuple[int, ...]]:
    """
    Generates spatial coordinates from a sequence of directional turns.

    Parameters
    ----------
    turns : Sequence[int]
        Sequence of direction indices.
    lattice : LatticeModel or str
        2D square or 3D cubic lattice.
    relative : bool
        If True, interprets turns relative to current heading (2D: 0=F, 1=L, 2=R, 3=B).
    fix_first_step : bool
        If True, fixes bead 0 at origin and bead 1 at (origin + +X),
        then applies the turns.
    origin : Optional[Tuple[int, ...]]
        Starting coordinate for bead 0.

    Returns
    -------
    List[Tuple[int, ...]]
        List of spatial coordinates for all beads.
    """
    lat_enum = LatticeModel(lattice) if isinstance(lattice, str) else lattice

    if lat_enum == LatticeModel.SQUARE_2D:
        start_coord = origin if origin is not None else (0, 0)
        coords: List[Tuple[int, ...]] = [start_coord]
        curr_x, curr_y = start_coord
        curr_dir = int(Direction2D.RIGHT)

        if fix_first_step:
            dx, dy = VECTORS_2D[curr_dir]
            curr_x += dx
            curr_y += dy
            coords.append((curr_x, curr_y))

        for t in turns:
            if relative:
                # 0: straight, 1: turn left (+1 mod 4), 2: turn right (-1 mod 4), 3: backtrack (+2 mod 4)
                if t == 0:
                    pass
                elif t == 1:
                    curr_dir = (curr_dir + 1) % 4
                elif t == 2:
                    curr_dir = (curr_dir - 1) % 4
                elif t == 3:
                    curr_dir = (curr_dir + 2) % 4
                else:
                    raise ValueError(f"Invalid 2D relative turn: {t}. Must be 0, 1, 2, or 3.")
            else:
                curr_dir = int(t) % 4

            dx, dy = VECTORS_2D[curr_dir]
            curr_x += dx
            curr_y += dy
            coords.append((curr_x, curr_y))

    elif lat_enum == LatticeModel.CUBIC_3D:
        start_coord = origin if origin is not None else (0, 0, 0)
        coords = [start_coord]
        curr_x, curr_y, curr_z = start_coord
        curr_dir = int(Direction3D.POS_X)

        if fix_first_step:
            dx, dy, dz = VECTORS_3D[curr_dir]
            curr_x += dx
            curr_y += dy
            curr_z += dz
            coords.append((curr_x, curr_y, curr_z))

        for t in turns:
            curr_dir = int(t) % 6
            dx, dy, dz = VECTORS_3D[curr_dir]
            curr_x += dx
            curr_y += dy
            curr_z += dz
            coords.append((curr_x, curr_y, curr_z))
    else:
        raise ValueError(f"Unsupported lattice model: {lattice}")

    return coords


def compute_coordinates(
    turns: Sequence[int],
    lattice: Union[LatticeModel, str] = LatticeModel.SQUARE_2D,
    relative: bool = False,
    fix_first_step: bool = False,
    origin: Optional[Tuple[int, ...]] = None,
) -> List[Tuple[int, ...]]:
    """Alias for turns_to_coordinates for interface compliance."""
    return turns_to_coordinates(
        turns, lattice=lattice, relative=relative, fix_first_step=fix_first_step, origin=origin
    )


def is_self_avoiding(coords: Sequence[Tuple[int, ...]]) -> bool:
    """Checks whether the coordinate walk is strictly self-avoiding (no spatial overlaps)."""
    return len(coords) == len(set(coords))


def count_overlaps(coords: Sequence[Tuple[int, ...]]) -> int:
    """Computes total pairwise spatial collisions sum_{i < j} delta(r_i, r_j)."""
    n = len(coords)
    collisions = 0
    for i in range(n):
        for j in range(i + 1, n):
            if coords[i] == coords[j]:
                collisions += 1
    return collisions


def has_backtracking(coords: Sequence[Tuple[int, ...]]) -> bool:
    """Checks if consecutive bond vectors reverse 180 degrees (r_{k+2} == r_k)."""
    if len(coords) < 3:
        return False
    for k in range(len(coords) - 2):
        if coords[k] == coords[k + 2]:
            return True
    return False


def count_backtracks(coords: Sequence[Tuple[int, ...]]) -> int:
    """Counts the number of immediate 180-degree backtracking steps."""
    if len(coords) < 3:
        return 0
    return sum(1 for k in range(len(coords) - 2) if coords[k] == coords[k + 2])


def get_contacts(coords: Sequence[Tuple[int, ...]]) -> List[Tuple[int, int]]:
    """
    Identifies all topological contacts: pairs (i, j) with |i - j| >= 2 at Euclidean distance 1.
    """
    n = len(coords)
    contacts: List[Tuple[int, int]] = []
    for i in range(n):
        for j in range(i + 2, n):
            sq_dist = sum((coords[i][d] - coords[j][d]) ** 2 for d in range(len(coords[i])))
            if sq_dist == 1:
                contacts.append((i, j))
    return contacts


class LatticeProtein:
    """
    Representation of a folded protein conformation on a discrete spatial lattice.
    """

    def __init__(
        self,
        sequence: str,
        coords: Optional[Sequence[Tuple[int, ...]]] = None,
        turns: Optional[Sequence[int]] = None,
        lattice_model: Union[LatticeModel, str] = LatticeModel.SQUARE_2D,
        lattice: Optional[Union[LatticeModel, str]] = None,
        hp_sequence: Optional[str] = None,
    ) -> None:
        self.sequence = sequence.strip().upper()
        if not self.sequence:
            raise ValueError("Sequence cannot be empty.")

        lat = lattice if lattice is not None else lattice_model
        self.lattice = LatticeModel(lat) if isinstance(lat, str) else lat
        self.lattice_model = self.lattice

        self.hp_sequence = hp_sequence if hp_sequence else encode_hp(self.sequence)

        if turns is not None:
            self.turns = list(turns)
        else:
            self.turns = []

        if coords is not None:
            if len(coords) != len(self.sequence):
                raise ValueError(
                    f"Coordinate length ({len(coords)}) must match sequence length ({len(self.sequence)})."
                )
            self.coords = list(coords)
        elif self.turns:
            self.coords = self.compute_coordinates(self.turns)
        else:
            # Default to straight coordinates
            straight_turns = [0] * (len(self.sequence) - 1)
            self.coords = self.compute_coordinates(straight_turns)

    @classmethod
    def from_turns(
        cls,
        sequence: str,
        turns: Sequence[int],
        lattice_model: Union[LatticeModel, str] = LatticeModel.SQUARE_2D,
        relative: bool = False,
    ) -> "LatticeProtein":
        """
        Constructs a LatticeProtein from turn sequence with first step fixed.
        """
        seq_clean = sequence.strip().upper()
        expected_turns = len(seq_clean) - 2
        if len(turns) != expected_turns:
            raise ValueError(
                f"Expected {expected_turns} turns for sequence of length {len(seq_clean)}, got {len(turns)}."
            )
        coords = turns_to_coordinates(
            turns, lattice=lattice_model, relative=relative, fix_first_step=True
        )
        return cls(
            sequence=seq_clean,
            coords=coords,
            turns=list(turns),
            lattice_model=lattice_model,
        )

    def compute_coordinates(
        self,
        turns: Sequence[int],
        relative: bool = False,
    ) -> List[Tuple[int, ...]]:
        """
        Generates coordinates for this protein from a turn sequence.
        
        Handles both:
        - len(turns) == len(sequence) - 1: full turns starting from origin.
        - len(turns) == len(sequence) - 2: remaining turns with first step (+X) fixed.
        """
        n_res = len(self.sequence)
        t_list = list(turns)
        if len(t_list) == n_res - 2:
            return turns_to_coordinates(
                t_list, lattice=self.lattice, relative=relative, fix_first_step=True
            )
        elif len(t_list) == n_res - 1:
            return turns_to_coordinates(
                t_list, lattice=self.lattice, relative=relative, fix_first_step=False
            )
        else:
            # Flexible: generate as many steps as turns provided
            return turns_to_coordinates(
                t_list, lattice=self.lattice, relative=relative, fix_first_step=False
            )

    def is_self_avoiding(self, coords: Optional[Sequence[Tuple[int, ...]]] = None) -> bool:
        """Checks if coords or self.coords is strictly self-avoiding."""
        c = coords if coords is not None else self.coords
        return is_self_avoiding(c)

    def has_backtracking(self, coords: Optional[Sequence[Tuple[int, ...]]] = None) -> bool:
        """Checks if coords or self.coords has immediate 180-degree backtracking."""
        c = coords if coords is not None else self.coords
        return has_backtracking(c)

    def get_contacts(self, coords: Optional[Sequence[Tuple[int, ...]]] = None) -> List[Tuple[int, int]]:
        """Returns list of contact pairs (i, j) for coords or self.coords."""
        c = coords if coords is not None else self.coords
        return get_contacts(c)

    @property
    def contacts(self) -> List[Tuple[int, int]]:
        return self.get_contacts(self.coords)

    @property
    def length(self) -> int:
        return len(self.sequence)

    @property
    def num_turns(self) -> int:
        return max(0, len(self.sequence) - 2)

    @property
    def is_valid(self) -> bool:
        return self.is_self_avoiding() and not self.has_backtracking()

    def to_dict(self) -> Dict[str, Any]:
        return {
            "sequence": self.sequence,
            "hp_sequence": self.hp_sequence,
            "length": self.length,
            "lattice": self.lattice.value,
            "turns": self.turns,
            "coords": self.coords,
            "n_contacts": len(self.contacts),
            "contacts": self.contacts,
            "is_valid": self.is_valid,
        }
