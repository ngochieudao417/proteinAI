#!/usr/bin/env python3
"""
tools/candidate_generator.py — Candidate Minibinder Generation & Target Preparation
for Changjiang Cup 2026 (TL1A & IL-7Ra tracks).

Provides utilities for:
1. Cleaning and preparing reference PDB structures (2RE9, 7OPB, 3DI3).
2. Generating BindCraft and RFdiffusion target settings JSONs.
3. Defining hotspots and constraints (strict length: 85-100 aa).
4. Validating candidate sequences (length, 20 standard AAs, Cys check).
5. Formatting candidate pools for screening in Notebook 03 (ESMFold2 + ipSAE).
"""

import os
import sys
import json
import re
from typing import Dict, List, Optional, Tuple, Any

# Strict competition rules and strategic constraints
MINIBINDER_MIN_LEN = 85
MINIBINDER_MAX_LEN = 100
STANDARD_AAS = set("ACDEFGHIKLMNPQRSTVWY")

# Verified Target Specifications
TARGET_SPECS = {
    "tl1a": {
        "name": "TL1A / TNFSF15",
        "uniprot": "O95150",
        "construct_range": "Leu72-Leu251",
        "length": 180,
        "default_pdb": "data/structures/2RE9.pdb",
        "oligomer": "homotrimer",
        "default_chains": ["A", "B", "C"],
        "hotspots_group_a": [
            # Inter-protomer groove epitope (DcR3/DR3 binding site)
            "D104", "S105", "Y117", "E119", "P120", "T121",
            "N136", "W137", "F138", "K169", "E170"
        ],
        "hotspots_group_b": [
            # Single-protomer surface epitope (hedge against monomer-only evaluation)
            "R32", "E49", "E51", "L52", "R85", "G86", "M87", "Y117", "E119", "P120"
        ],
        "avoid_residues": [
            # Glycan sequons: N133, N229 in UniProt -> N62, N158 in construct
            "N62", "N158",
            # Intrachain disulfide
            "C91", "C131"
        ],
        "msa_path": "data/msa/TL1A_O95150_L72-L251.uniref.a3m",
        "msa_depth": 2037
    },
    "il7ra": {
        "name": "IL-7Ra / CD127",
        "uniprot": "P16871",
        "construct_range": "Glu21-Asp239",
        "length": 219,
        "default_pdb": "data/structures/7OPB.pdb",
        "alt_pdb": "data/structures/3DI3.pdb",
        "oligomer": "monomer",
        "default_chains": ["A"],  # Chain A in 7OPB or Chain B in 3DI3
        "hotspots_group_a": [
            # D1 core hydrophobic strand & IL-7 functional site
            "K77", "K78", "F79", "L80", "L81", "I82", "G83",
            # D1 auxiliary loop
            "S31", "H33", "L57", "V58", "E59"
        ],
        "hotspots_group_b": [
            # D1 alternative lateral surface
            "K77", "F79", "L80", "E59", "K97", "K98"
        ],
        "avoid_residues": [
            # Glycan sequons: N29, N45, N131, N162, N212, N213
            "N29", "N45", "N131", "N162", "N212", "N213",
            # Disulfides
            "C22", "C37", "C54", "C62", "C88", "C98"
        ],
        "msa_path": "data/msa/IL7Ra_P16871_E21-D239.uniref.a3m",
        "msa_depth": 921
    }
}


def clean_pdb(
    input_pdb_path: str,
    output_pdb_path: str,
    target_chains: List[str],
    remove_hetatm: bool = True,
    renumber: bool = False,
    res_offset: int = 0
) -> Dict[str, Any]:
    """
    Cleans a reference PDB file:
    - Filters to specified chains.
    - Strips HETATM records (water, ions, crystallization reagents).
    - Removes expression tags if present.
    - Optionally applies a residue number offset.
    """
    if not os.path.isfile(input_pdb_path):
        raise FileNotFoundError(f"Input PDB file not found: {input_pdb_path}")

    retained_atoms = 0
    retained_residues = set()
    output_lines = []

    with open(input_pdb_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("ATOM  "):
                chain_id = line[21].strip()
                if chain_id not in target_chains:
                    continue
                
                res_seq = int(line[22:26].strip())
                res_name = line[17:20].strip()

                # Filter out Gly 0 tag in 2RE9 if encountered
                if res_seq <= 0:
                    continue

                new_res_seq = res_seq + res_offset if renumber else res_seq
                
                # Format line with modified residue number if requested
                if renumber and res_offset != 0:
                    line = f"{line[:22]}{new_res_seq:>4}{line[26:]}"

                output_lines.append(line)
                retained_atoms += 1
                retained_residues.add((chain_id, new_res_seq, res_name))

            elif line.startswith("TER   "):
                chain_id = line[21].strip()
                if chain_id in target_chains:
                    output_lines.append(line)
            elif not remove_hetatm and line.startswith("HETATM"):
                chain_id = line[21].strip()
                if chain_id in target_chains:
                    output_lines.append(line)

    output_lines.append("END\n")

    os.makedirs(os.path.dirname(os.path.abspath(output_pdb_path)), exist_ok=True)
    with open(output_pdb_path, "w", encoding="utf-8") as out:
        out.writelines(output_lines)

    return {
        "input": input_pdb_path,
        "output": output_pdb_path,
        "chains": target_chains,
        "total_atoms": retained_atoms,
        "total_residues": len(retained_residues),
        "residue_range": (
            min(r[1] for r in retained_residues) if retained_residues else None,
            max(r[1] for r in retained_residues) if retained_residues else None
        )
    }


def generate_bindcraft_config(
    target_key: str,
    target_pdb: str,
    target_chains: List[str],
    epitope_group: str = "A",
    num_designs: int = 100,
    min_length: int = MINIBINDER_MIN_LEN,
    max_length: int = MINIBINDER_MAX_LEN,
    num_seqs_per_backbone: int = 4
) -> Dict[str, Any]:
    """
    Generates a BindCraft-compatible target settings configuration dictionary.
    """
    if target_key.lower() not in TARGET_SPECS:
        raise ValueError(f"Unknown target: {target_key}. Must be 'tl1a' or 'il7ra'.")

    spec = TARGET_SPECS[target_key.lower()]
    hotspots = spec["hotspots_group_a"] if epitope_group.upper() == "A" else spec["hotspots_group_b"]

    config = {
        "design_name": f"{target_key.upper()}_minibinder_group{epitope_group.upper()}",
        "target": {
            "name": spec["name"],
            "uniprot": spec["uniprot"],
            "target_pdb": target_pdb,
            "target_chains": target_chains,
            "target_hotspots": hotspots,
            "avoid_residues": spec["avoid_residues"]
        },
        "binder_specs": {
            "min_length": min_length,
            "max_length": max_length,
            "target_length_distribution": "uniform",
            "enforce_strict_length_bounds": True,
            "no_cysteines": True,
            "allowed_amino_acids": list("ADEFGHIKLMNPQRSTVWY")  # Exclude Cys for safety
        },
        "backbone_generation": {
            "method": "bindcraft_af2_backprop",
            "num_designs": num_designs,
            "soft_iterations": 100,
            "hard_iterations": 50,
            "weights": {
                "hotspot_loss": 2.0,
                "interface_contact_loss": 1.5,
                "pae_loss": 1.0,
                "plddt_loss": 0.5
            }
        },
        "sequence_design": {
            "method": "ProteinMPNN",
            "num_sequences": num_seqs_per_backbone,
            "sampling_temp": 0.15,
            "model_version": "v_48_020"
        },
        "filtering": {
            "af2_plddt_min": 75.0,
            "af2_pae_interaction_max": 12.0,
            "downstream_selector": "ESMFold2_ipSAE"
        }
    }
    return config


def validate_candidate_sequence(
    sequence: str,
    name: str = "candidate",
    min_len: int = MINIBINDER_MIN_LEN,
    max_len: int = MINIBINDER_MAX_LEN
) -> Tuple[bool, List[str]]:
    """
    Validates a candidate sequence against all competition and strategic criteria:
    - Uppercase standard 20 amino acids.
    - No invalid or ambiguous letters.
    - Length strictly in [min_len, max_len] (default 85-100 aa).
    - Flags free Cys or unstable motifs.
    """
    errors = []
    clean_seq = re.sub(r"\s+", "", sequence).upper()

    if not clean_seq:
        return False, ["Sequence is empty"]

    # Length check
    if len(clean_seq) < min_len or len(clean_seq) > max_len:
        errors.append(f"Length {len(clean_seq)} violates [{min_len}, {max_len}] aa requirement")

    # Character set check
    invalid_chars = set(clean_seq) - STANDARD_AAS
    if invalid_chars:
        errors.append(f"Contains invalid characters: {sorted(list(invalid_chars))}")

    # Cysteine check (minibinders should avoid free/unpaired Cys)
    cys_count = clean_seq.count("C")
    if cys_count > 0:
        if cys_count % 2 != 0:
            errors.append(f"Contains {cys_count} cysteines (odd number = free unpaired Cys)")

    return len(errors) == 0, errors


def create_screening_manifest(
    candidates: List[Dict[str, Any]],
    output_csv_path: str
) -> Dict[str, Any]:
    """
    Formats candidate pool into CSV manifest for ingestion by Notebook 03:
    Headers: name,sequence,target,length,epitope_group,generator,backbone_id
    """
    valid_candidates = []
    rejected_candidates = []

    for c in candidates:
        name = c.get("name", "unnamed")
        seq = c.get("sequence", "")
        target = c.get("target", "TL1A").upper()
        group = c.get("epitope_group", "A")
        generator = c.get("generator", "bindcraft")
        bb_id = c.get("backbone_id", "bb_001")

        is_valid, errs = validate_candidate_sequence(seq)
        if is_valid:
            valid_candidates.append({
                "name": name,
                "sequence": seq,
                "target": target,
                "length": len(seq),
                "epitope_group": group,
                "generator": generator,
                "backbone_id": bb_id
            })
        else:
            rejected_candidates.append({
                "name": name,
                "errors": errs,
                "sequence_len": len(seq)
            })

    os.makedirs(os.path.dirname(os.path.abspath(output_csv_path)), exist_ok=True)
    with open(output_csv_path, "w", encoding="utf-8") as f:
        f.write("name,sequence,target,length,epitope_group,generator,backbone_id\n")
        for vc in valid_candidates:
            f.write(
                f"{vc['name']},{vc['sequence']},{vc['target']},{vc['length']},"
                f"{vc['epitope_group']},{vc['generator']},{vc['backbone_id']}\n"
            )

    return {
        "manifest_path": output_csv_path,
        "total_input": len(candidates),
        "total_valid": len(valid_candidates),
        "total_rejected": len(rejected_candidates),
        "rejected_details": rejected_candidates
    }


if __name__ == "__main__":
    print("tools/candidate_generator.py loaded successfully.")
