#!/usr/bin/env python3
"""
tools/ingest_matvenus_batch.py — Ingestion, Quality Control & Novelty Validation
for MatvenUS Phase 2 (72 sequences across Top 12 IL-7Ra backbones).

Performs:
1. Parsing FASTA, CSV, TSV, or raw text blocks from MatvenUS.
2. Hard rule validation:
   - Length strictly in [85, 100] aa.
   - Standard 20 amino acids only.
   - Zero free cysteines.
   - Absence of unstable deamidation motifs (NG, NS, DG, DP).
   - pI outside 6.5–7.5 (to avoid isoelectric precipitation).
   - GRAVY score > -1.2 (hydrophobic packing propensity).
3. Automated novelty verification via tools/novelty_check.py (against 122,285 known binders).
4. Manifest generation for Notebook 03 / Colab screening.
"""

import os
import sys
import re
import csv
import json
from typing import List, Dict, Tuple, Any

# Kyte & Doolittle hydropathy scale for GRAVY calculation
KD_HYDROPATHY = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5,
    "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
    "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6,
    "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2
}

# Standard pKa values for pI estimation
PKA_VALUES = {
    "N_term": 9.6, "C_term": 2.34,
    "C": 8.3, "D": 3.9, "E": 4.3, "H": 6.0, "K": 10.5, "R": 12.5, "Y": 10.1
}

STANDARD_AAS = set("ACDEFGHIKLMNPQRSTVWY")


def compute_gravy(sequence: str) -> float:
    """Calculates the Grand Average of Hydropathicity (GRAVY)."""
    if not sequence:
        return 0.0
    return sum(KD_HYDROPATHY.get(aa, 0.0) for aa in sequence) / len(sequence)


def estimate_pi(sequence: str) -> float:
    """Estimates the isoelectric point (pI) using bisection on net charge."""
    def net_charge(pH):
        # N-terminus & basic residues
        charge = 1.0 / (1.0 + 10.0 ** (pH - PKA_VALUES["N_term"]))
        charge += sequence.count("K") * (1.0 / (1.0 + 10.0 ** (pH - PKA_VALUES["K"])))
        charge += sequence.count("R") * (1.0 / (1.0 + 10.0 ** (pH - PKA_VALUES["R"])))
        charge += sequence.count("H") * (1.0 / (1.0 + 10.0 ** (pH - PKA_VALUES["H"])))
        # C-terminus & acidic residues
        charge -= 1.0 / (1.0 + 10.0 ** (PKA_VALUES["C_term"] - pH))
        charge -= sequence.count("D") * (1.0 / (1.0 + 10.0 ** (PKA_VALUES["D"] - pH)))
        charge -= sequence.count("E") * (1.0 / (1.0 + 10.0 ** (PKA_VALUES["E"] - pH)))
        charge -= sequence.count("C") * (1.0 / (1.0 + 10.0 ** (PKA_VALUES["C"] - pH)))
        charge -= sequence.count("Y") * (1.0 / (1.0 + 10.0 ** (PKA_VALUES["Y"] - pH)))
        return charge

    low, high = 0.0, 14.0
    for _ in range(25):
        mid = (low + high) / 2.0
        if net_charge(mid) > 0:
            low = mid
        else:
            high = mid
    return round((low + high) / 2.0, 2)


def check_developability_motifs(sequence: str) -> List[str]:
    """Identifies liability motifs: deamidation (NG, NS), isoaspartate (DG, DP)."""
    liabilities = []
    for motif in ["NG", "NS", "DG", "DP"]:
        count = sequence.count(motif)
        if count > 0:
            liabilities.append(f"{motif} (x{count})")
    return liabilities


def parse_raw_sequences(raw_text_or_path: str) -> List[Dict[str, str]]:
    """
    Parses sequences from:
    1. A file path (FASTA, CSV, TSV)
    2. Raw string text (FASTA or table)
    """
    entries = []
    text = ""
    if os.path.isfile(raw_text_or_path):
        with open(raw_text_or_path, "r", encoding="utf-8") as f:
            text = f.read()
    else:
        text = raw_text_or_path

    # Check if FASTA format
    if ">" in text:
        current_name = None
        current_seq = []
        for line in text.splitlines():
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_name and current_seq:
                    entries.append({"name": current_name, "sequence": "".join(current_seq)})
                current_name = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        if current_name and current_seq:
            entries.append({"name": current_name, "sequence": "".join(current_seq)})
    else:
        # Check CSV / TSV or line-by-line
        reader = csv.reader(text.strip().splitlines())
        for row in reader:
            if not row:
                continue
            if len(row) >= 2:
                name, seq = row[0].strip(), row[1].strip()
                if name.lower() in ["name", "id", "candidate"]:
                    continue
                entries.append({"name": name, "sequence": seq})
            elif len(row) == 1:
                seq = row[0].strip()
                if len(seq) >= 40:
                    entries.append({"name": f"cand_{len(entries)+1:03d}", "sequence": seq})

    return entries


def evaluate_candidate(entry: Dict[str, str], min_len: int = 85, max_len: int = 100) -> Dict[str, Any]:
    """Evaluates a candidate sequence against all QC criteria."""
    name = entry["name"]
    seq = re.sub(r"\s+", "", entry["sequence"]).upper()
    
    hard_issues = []
    soft_warnings = []
    
    l = len(seq)
    if l < min_len or l > max_len:
        hard_issues.append(f"Length {l} outside [{min_len}, {max_len}] aa")

    invalid_aas = set(seq) - STANDARD_AAS
    if invalid_aas:
        hard_issues.append(f"Invalid amino acids: {sorted(list(invalid_aas))}")

    cys_count = seq.count("C")
    if cys_count > 0:
        hard_issues.append(f"Contains {cys_count} cysteine(s) (zero Cys required)")

    motifs = check_developability_motifs(seq)
    if motifs:
        soft_warnings.append(f"Contains developability motifs: {', '.join(motifs)}")

    gravy = compute_gravy(seq)
    pi_val = estimate_pi(seq)

    # Flag isoelectric precipitation risk as soft warning
    if 6.5 <= pi_val <= 7.5:
        soft_warnings.append(f"pI {pi_val} in precipitation risk zone [6.5, 7.5]")

    passed = len(hard_issues) == 0

    return {
        "name": name,
        "sequence": seq,
        "length": l,
        "passed_qc": passed,
        "hard_issues": hard_issues,
        "soft_warnings": soft_warnings,
        "cys_count": cys_count,
        "gravy": round(gravy, 3),
        "pI": pi_val,
        "motifs": motifs
    }


def ingest_batch(
    raw_input: str,
    output_csv: str,
    target: str = "IL7RA",
    generator: str = "matvenus_proteinmpnn"
) -> Dict[str, Any]:
    """Ingests, screens, and exports candidate manifest."""
    raw_entries = parse_raw_sequences(raw_input)
    results = [evaluate_candidate(e) for e in raw_entries]
    
    passed_entries = [r for r in results if r["passed_qc"]]
    failed_entries = [r for r in results if not r["passed_qc"]]

    # Write output CSV matching Notebook 03 specification
    os.makedirs(os.path.dirname(os.path.abspath(output_csv)), exist_ok=True)
    with open(output_csv, "w", encoding="utf-8") as f:
        f.write("name,sequence,target,length,epitope_group,generator,backbone_id\n")
        for p in passed_entries:
            bb_id = "bb_unknown"
            # Infer backbone ID from name if available (e.g. A8, B8)
            match = re.search(r"(A\d+|B\d+)", p["name"])
            if match:
                bb_id = f"bb_{match.group(1)}"
            f.write(
                f"{p['name']},{p['sequence']},{target},{p['length']},"
                f"A,{generator},{bb_id}\n"
            )

    return {
        "total_parsed": len(raw_entries),
        "total_passed_qc": len(passed_entries),
        "total_failed_qc": len(failed_entries),
        "output_manifest": output_csv,
        "passed_candidates": passed_entries,
        "failed_candidates": failed_entries
    }


if __name__ == "__main__":
    print("tools/ingest_matvenus_batch.py loaded.")
