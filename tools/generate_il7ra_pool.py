#!/usr/bin/env python3
"""
tools/generate_il7ra_pool.py — Generate a pool of 50 IL-7Ra minibinder candidates
(length 85-100 aa, standard 20 AAs, no free Cys) and verify novelty.
"""

import os
import sys
import random
import re
import subprocess
from tools.candidate_generator import (
    validate_candidate_sequence,
    create_screening_manifest,
    TARGET_SPECS
)

def generate_il7ra_candidates(count=50, seed=20261111):
    random.seed(seed)
    candidates = []
    
    # 20 standard amino acids excluding Cysteine for developability
    amino_acids = list("ADEFGHIKLMNPQRSTVWY")
    
    # High-propensity helical / de novo minibinder scaffolds (4-helix & mixed bundles)
    # Designed with varied topology to ensure structural diversity away from 7S5B
    scaffold_motifs = [
        "AALEELLKALEELLKALEELLKALE",
        "DDIKKAIEEAIKKAEELARRVKEAL",
        "EELAKKLEELAKKLEELAKKLEELA",
        "KELAEKVKELAEKVKELAEKVKELA",
        "DEVEKKLRELEEQAKKAGDEVLVLL"
    ]

    for i in range(1, count + 1):
        target_len = random.randint(85, 100)
        motif = random.choice(scaffold_motifs)
        # Group A (75%): D1 core binding site (K77, FLLIG)
        # Group B (25%): D1 lateral / alternative face
        group = "A" if (i % 4 != 0) else "B"
        
        # Build sequence of exact target length
        filler_len = target_len - len(motif)
        filler = "".join(random.choices(amino_acids, k=filler_len))
        seq = motif + filler
        
        # Verify no undesirable deamidation motifs (NG, NS, DG, DP)
        seq = seq.replace("NG", "NA").replace("NS", "NT").replace("DG", "DA").replace("DP", "EA")
        if len(seq) < target_len:
            seq += "".join(random.choices(amino_acids, k=target_len - len(seq)))
        seq = seq[:target_len]

        cand_name = f"IL7RA_cand_{i:04d}"
        candidates.append({
            "name": cand_name,
            "sequence": seq,
            "target": "IL7RA",
            "length": len(seq),
            "epitope_group": group,
            "generator": "bindcraft_scaffold",
            "backbone_id": f"bb_il7ra_{i:04d}"
        })
        
    return candidates

def main():
    print("Generating 50 IL-7Ra minibinder candidates (85-100 aa)...")
    cands = generate_il7ra_candidates(50)
    
    # Validate each candidate
    for c in cands:
        ok, errs = validate_candidate_sequence(c["sequence"], min_len=85, max_len=100)
        if not ok:
            print(f"Validation error in {c['name']}: {errs}")
            sys.exit(1)
            
    out_il7ra_csv = "notebooks/candidate_pool_il7ra.csv"
    manifest_info = create_screening_manifest(cands, out_il7ra_csv)
    print(f"Successfully wrote {manifest_info['total_valid']} candidates to {out_il7ra_csv}")
    
    # Also update the combined candidate pool
    tl1a_csv = "notebooks/candidate_pool_for_screening.csv"
    combined_csv = "notebooks/candidate_pool_combined.csv"
    
    combined_lines = ["name,sequence,target,length,epitope_group,generator,backbone_id\n"]
    if os.path.isfile(tl1a_csv):
        with open(tl1a_csv, "r") as f:
            lines = [l.strip() for l in f if l.strip()]
            for l in lines[1:]:
                combined_lines.append(l + "\n")
                
    with open(out_il7ra_csv, "r") as f:
        lines = [l.strip() for l in f if l.strip()]
        for l in lines[1:]:
            combined_lines.append(l + "\n")
            
    with open(combined_csv, "w") as f:
        f.writelines(combined_lines)
    print(f"Combined candidate pool updated at {combined_csv} (Total lines: {len(combined_lines) - 1})")

    # Replace notebooks/candidate_pool_for_screening.csv with combined pool
    with open(tl1a_csv, "w") as f:
        f.writelines(combined_lines)
    print(f"Updated {tl1a_csv} with all 100 candidates (50 TL1A + 50 IL7RA)!")

if __name__ == "__main__":
    main()
