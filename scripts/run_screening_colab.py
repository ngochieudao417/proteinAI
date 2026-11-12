#!/usr/bin/env python3
"""
scripts/run_screening_colab.py — Headless candidate screening pipeline for Colab L4.

Executes ESMFold2 structure prediction + DunbrackLab ipSAE + Smith-Waterman Novelty Check
for minibinder candidates against competition targets (TL1A trimer and IL-7Rα monomer).

Key features:
1. Default target: TL1A trimer (chains A, B, C broadcast MSA; candidate chain Z).
2. Checkpoint resume: Automatically scans existing CSV and partial JSON on Drive,
   skipping completed candidates and preserving evaluation records.
3. Token-order PDB dumping: Exact token alignment ensuring permutation symmetry for ipSAE.
4. Auto-discovery of candidate CSVs, weights, MSAs, and ipsae.py.
5. Automated Top 20 exports per target (03_screening_top_TL1A.csv / 03_screening_top_IL7RA.csv).
"""

import argparse
import csv
import gc
import hashlib
import json
import os
import platform
import shutil
import subprocess
import sys
import time
import urllib.request
import warnings
from pathlib import Path

# Set expandable segments before importing PyTorch
os.environ.setdefault("PYTORCH_CUDA_ALLOC_CONF", "expandable_segments:True")

# ==============================================================================
# Target Definitions & Constants
# ==============================================================================

IL7RA_SEQ = (
    "ESGYAQNGDLEDAELDDYSFSCYSQLEVNGSQHSLTCAFEDPDVNITNLEFEICGALVEVKCLNFRKLQEIYFIETK"
    "KFLLIGKSNICVKVGEKSLTCKKIDLTTIVKPEAPFDLSVVYREGANDFVVTFNTSHLQKKYVKVLMHDVAYRQEKD"
    "ENKWTHVNLSSTKLTLLQRKLQPAAMYEIKVRSIPDHYFKGFWSEWSPSYYFRTPEINNSSGEMD"
)
TL1A_SEQ = (
    "LKGQEFAPSHQQVYAPLRADGDKPRAHLTVVRQTPTQHFKNQFPALHWEHELGLAFTKNRMNYTNKFLLIPESGDYF"
    "IYSQVTFRGMTSECSEIRQAGRPNKPDSITVVITKVTDSYPEPTQLLMGTKSVCEVGSNWFQPIYLGAMFSLQEGDK"
    "LMVNVSDISLVDYTKEDKTFFGAFLL"
)

TARGET_METAS = {
    "IL7RA": {
        "name": "IL7Ra_P16871_E21-D239",
        "seq": IL7RA_SEQ,
        "len": 219,
        "sha256": "08e2a05b5b4787bb28c03302b1d56c3856980423a5582f2096bc6f5238b0be03",
    },
    "TL1A": {
        "name": "TL1A_O95150_L72-L251",
        "seq": TL1A_SEQ,
        "len": 180,
        "sha256": "4464ad15470e2ecda72ea879beb5972059062085b04ff515cef26f642e092c27",
    },
}

CORE_KNOWN_BINDERS = {
    "tl1a": [
        ("DcR3_TNFRSF6B_res30_205", "VAETPTYPWRDAETGERLVCAQCPPGTFVQRPCRRDSPTTCGPCPPRHYTQFWNYLERCRYCNVLCGEREEEARACHATHNRACRCRTGFFAHAGFCLEHASCPPGAGVIAPGTPSQNTQCQPCPPGTFSASSSSSEQCQPHRNCTALGLALNVPGSSSHDTLCTSCTGFPLSTRV"),
        ("Afimkibart_VH", "QVQLVQSGAEVKKPGASVKVSCKASGYDFTYYGISWVRQAPGQGLEWMGWISTYNGNTHYARMLQGRVTMTTDTSTRTAYMELRSLRSDDTAVYYCARENYYGSGAYRGGMDVWGQGTTVTVSS"),
        ("Afimkibart_VL", "EIVLTQSPATLSLSPGERATLSCRASQSVSSYLAWYQQKPGQAPRLLIYDASNRATGIPARFSGSGSGTDFTLTISSLEPEDFAVYYCQQRSNWPWTFGQGTKVEIK"),
        ("Duvakitug_VH", "QVQLVQSGAEVKKPGASVKVSCKASGYTFTSYDINWVRQAPGQGLEWMGWLNPNSGYTGYAQKFQGRVTMTADRSTSTAYMELSSLRSEDTAVYYCAREVPETAAFEYWGQGTLVTVSS"),
        ("Duvakitug_VL", "QSVLTQPPSVSGAPGQRVTISCTSSSSDIGAGLGVHWYQQLPGTAPKLLIEGYYNRPSGVPDDRFSGSKSGTSASLTITGLLPEDEGDYYCQSWDGTLSALFGGGTKLTVL"),
        ("Tulisokibart_VH", "QVQLVQSGAEVKKPGASVKVSCKASGFDIQDTYMHWVKQRPGQGLEWMGRIDPASGHTKYDPKFQVRVTITRDTSTSTVYLELSSLRSEDTAVYYCARSGGLPDVWGQGTTVTVSS"),
        ("Tulisokibart_VL", "EIVLTQSPGTLSLSPGERATLSCRASSSVSYMYWYQQKPGQAPRPLIYATSNLASGIPDRFSGSGSGTDFTLTISRLEPEDFAVYYCQQWEGNPRTFGGGTKLEIK"),
    ],
    "il7ra": [
        ("7S5B_minibinder_optimized", "SVIEKLRKLEKQARKQGDEVLVMLARMVLEYLEKGWVSEEDADESADRIEEVLKK"),
        ("7S5B_minibinder_original", "SVKKKVRKVEKKARKAGDELAVLLARRVLEALEKGLVSEEDADESADRIEEALKK"),
        ("Bempikibart_VH", "EVQLVESGGGLVQPGRSLRLSCAASGFTFDDHAMHWVRQAPGKGLEWVSGISWNSRGIGYADSVKGRFTIFRDNAKNSLYLQMNSLRAEDTALYYCAKDEYSRGYYVLDVWGQGTTVTVSS"),
        ("Bempikibart_VL", "AIQLTQSPSSLSASVGDRVTITCRASQGISSALAWYQQKPGKAPKLLIYDASSLESGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQFNSYPLWITFGQGTRLEIK"),
        ("Crebankitug_VH", "EVQLVESGGGLVKPGGSLRLSCAASGFTFDDSVMHWVRQAPGKGLEWVSLVGWDGFFTYYADSVKGRFTISRDNAKNSLYLQMNSLRAEDTAVYYCARQGDYMGNNWGQGTLVTVSS"),
        ("Crebankitug_VL", "NFMLTQPHSVSESPGKTVTISCTRSSGSIDSSYVQWYQQRPGSSPTTVIYEDDQRPSGVPDRFSGSIDSSSNSASLTISGLKTEDEADYYCQSYDFHHLVFGGGTKLTVL"),
        ("Lusvertikimab_VH", "QVQLVESGGGLVKPGGSLRLSCAVSGFTLSDYYMAWIRQAPGKGLEWVSTISASGLRTYYPDSVKGRFTISRDNAKNSLYLQMNSLRAEDTAVYYCARPLSAHYGFNYFDYWGQGTLVTVSS"),
        ("Lusvertikimab_VL", "DIQMTQSPSSLSASVGDRVTITCRTSEDIYQGLAWYQQKPGKAPKLLLYSANTLHIGVPSRFSGSGSGTDYTLTISSLQPEDFATYYCQQYYDYPLAFGGGTKVEIK"),
        ("Fab_4A10_VH", "QVQLQQPGAELVMPGASVKLSCKASGYTFTSYWMHWVKQRPGEGLEWIGEIDPSDSYTNDNQKFKGKATLTVDKSSSTAYMQLSSLTSEDSAVYYCARRLYSNSYYYAMDYWGQGTSVTVSS"),
        ("Fab_4A10_VL", "DIQMTQSPSSLSASLGGKVTITCKASQDIKKYIAWYQHKPGKGPRLLIHYTSTLQPGIPSRFSGSGSGRDYSFSISNLEPVDIATYYCLQYDNLLTFGAGTKLELKRADAAPTVSIFPPS"),
        ("Fab_2B8_VH", "EVQLQQSGPELVKPGASVKMSCKASGYTFSDYYMHWVKQSHGKSLEWIGYIYPDNGGNGYNQKFKGKATLTVDKSSSTVYMELRSLTSEDSALYYCARGTYYDGSYFDYWGQGTTLTVSS"),
        ("Fab_2B8_VL", "DIVMTQSHKFMSTLVGDRVSITCKASQDVSTTVAWYQQKPGQSPKLLIYSASYRYTGVPDRFTGSGSGTDFTFTISSVQAEDLAVYYCQQHYSIPRTFGGGTKLEIK"),
    ],
}

KNOWN_BINDER_FILES = {
    "il7ra": [
        ("IL7Ra_denovo_minibinder_designs_IPD_Cao2022.fasta", "IPD_Cao2022_122k"),
        ("IL7Ra_denovo_minibinders_VALIDATED.fasta",          "VALIDATED"),
        ("IL7Ra_therapeutic_antibodies.fasta",                "therapeutics"),
        ("IL7Ra_PDB_antibody_Fabs.fasta",                     "PDB_Fabs"),
    ],
    "tl1a": [
        ("TL1A_therapeutic_antibodies.fasta", "therapeutics"),
        ("TL1A_natural_binder_DcR3.fasta",    "DcR3_natural"),
    ],
}

SAMPLE_BENCHMARKS = [
    {
        "name": "ctrl_pos_7S5B_IL7Ra",
        "sequence": "SVIEKLRKLEKQARKQGDEVLVMLARMVLEYLEKGWVSEEDADESADRIEEVLKK",
        "target": "IL7RA",
        "note": "Positive control (7S5B picomolar)",
    },
    {
        "name": "ctrl_neg_scrambled_IL7Ra",
        "sequence": "EVELVRLLRAKELVSEDWLSGEKRILEMKEDVLSKAYEAEKQVIDEMVAKRSDLE",
        "target": "IL7RA",
        "note": "Negative control (scrambled 7S5B)",
    },
    {
        "name": "probe_dcr3_TL1A",
        "sequence": (
            "VAETPTYPWRDAETGERLVCAQCPPGTFVQRPCRRDSPTTCGPCPPRHYTQFWNYLERCRYCNVLCGEREEEARACH"
            "ATHNRACRCRTGFFAHAGFCLEHASCPPGAGVIAPGTPSQNTQCQPCPPGTFSASSSSSEQCQPHRNCTALGLALNV"
            "PGSSSHDTLCTSCTGFPLSTRV"
        ),
        "target": "TL1A",
        "note": "Natural decoy receptor ectodomain (PDB 3K51)",
    },
]

VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")

# ==============================================================================
# Helper Functions: FASTA, PDB, ipSAE, Novelty Check
# ==============================================================================

def parse_fasta_records(path):
    records, header, seq_parts = [], None, []
    with open(path, "r", encoding="utf-8") as fh:
        for raw in fh:
            line = raw.strip()
            if not line: continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_parts).upper()))
                header = line[1:].strip()
                seq_parts = []
            else:
                seq_parts.append(line)
    if header is not None:
        records.append((header, "".join(seq_parts).upper()))
    return records


def parse_candidates_csv(path, default_target="TL1A"):
    cand_list = []
    ext = os.path.splitext(path)[1].lower()
    if ext in (".fa", ".fasta"):
        for h, s in parse_fasta_records(path):
            cand_list.append({"name": h.split()[0], "sequence": s, "target": default_target})
    elif ext in (".csv", ".tsv"):
        delim = "\t" if ext == ".tsv" else ","
        with open(path, "r", encoding="utf-8") as fh:
            reader = csv.DictReader(fh, delimiter=delim)
            fields_low = {f.lower().strip(): f for f in reader.fieldnames or []}
            name_col = next((fields_low[c] for c in ("name", "id", "seq_id", "candidate_id") if c in fields_low), None)
            seq_col  = next((fields_low[c] for c in ("sequence", "seq", "aa") if c in fields_low), None)
            tgt_col  = next((fields_low[c] for c in ("target", "target_name") if c in fields_low), None)
            if not (name_col and seq_col):
                raise ValueError(f"CSV requires 'name' and 'sequence' columns. Found: {reader.fieldnames}")
            for row in reader:
                nm = row[name_col].strip()
                sq = row[seq_col].strip().upper()
                tg = row[tgt_col].strip().upper() if tgt_col and row.get(tgt_col) else default_target
                if nm and sq:
                    cand_list.append({"name": nm, "sequence": sq, "target": tg})
    return cand_list


def write_pdb_token_order(mc, path):
    """Writes PDB strictly matching ESMFold2 token sequence. B-factor = pLDDT * 100."""
    chain_lookup = mc.metadata.chain_lookup
    lines, serial, n_ca, n_cb, res_counter = [], 1, 0, 0, {}
    for tok_idx, (s, e) in enumerate(mc.token_to_atoms):
        resname = str(mc.sequence[tok_idx])[:3]
        chain = str(chain_lookup[int(mc.chain_id[tok_idx])])[:1] or "A"
        res_counter[chain] = res_counter.get(chain, 0) + 1
        resnum = res_counter[chain]
        b = float(mc.plddt[tok_idx]) * 100.0
        has_ca = has_cb = False
        for a in range(int(s), int(e)):
            an = str(mc.atom_names[a]).strip().upper()
            el = str(mc.atom_elements[a]).strip().upper()[:2]
            x, y, z = (float(v) for v in mc.atom_positions[a])
            an4 = f" {an:<3}" if len(an) < 4 else an[:4]
            lines.append(f"ATOM  {serial:5d} {an4} {resname:>3s} {chain:1s}{resnum:4d}    "
                         f"{x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{b:6.2f}          {el:>2s}")
            serial += 1
            if an == "CA":
                has_ca = True
                if resname == "GLY": has_cb = True
            if an == "CB":
                has_cb = True
        n_ca += has_ca; n_cb += has_cb
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines + ["TER", "END"]) + "\n")
    return len(mc.token_to_atoms), n_ca, n_cb


def run_ipsae_eval(result, tag, outdir, binder_chain, ipsae_py, pae_cutoff=10, dist_cutoff=10):
    """Executes ipsae.py via AF2 path (.pdb + .json) and parses interface ipSAE."""
    import numpy as np
    os.makedirs(outdir, exist_ok=True)
    pdb_p  = f"{outdir}/{tag}.pdb"
    json_p = f"{outdir}/{tag}.json"
    n_res, n_ca, n_cb = write_pdb_token_order(result.complex, pdb_p)

    pae = result.pae.detach().cpu().float().numpy()
    plddt = result.plddt.detach().cpu().float().numpy().ravel()
    if plddt.max() <= 1.0:
        plddt = plddt * 100.0

    assert pae.shape == (n_res, n_res), f"PAE shape {pae.shape} != ({n_res}, {n_res})"
    with open(json_p, "w", encoding="utf-8") as f:
        json.dump({
            "pae": np.round(pae, 2).tolist(),
            "plddt": np.round(plddt, 2).tolist(),
            "iptm": float(result.iptm) if result.iptm is not None else -1.0,
            "ptm": float(result.ptm) if result.ptm is not None else -1.0
        }, f)

    cmd = [sys.executable, ipsae_py, json_p, pdb_p, str(pae_cutoff), str(dist_cutoff)]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
    txt_p = f"{outdir}/{tag}_{pae_cutoff:02d}_{dist_cutoff:02d}.txt"
    out = {"returncode": r.returncode, "pdb": pdb_p, "json": json_p, "txt": txt_p, "asym": {}, "max": {}}

    if r.returncode != 0 or not os.path.exists(txt_p):
        out["stderr"] = r.stderr[-2000:]
        return out

    for line in open(txt_p, encoding="utf-8"):
        p = line.split()
        if len(p) < 6 or p[0] == "Chn1": continue
        try: val = float(p[5])
        except ValueError: continue
        if p[4] == "asym":
            out["asym"][f"{p[0]}->{p[1]}"] = val
        elif p[4] == "max":
            out["max"][f"{p[0]}-{p[1]}"] = val

    binder_asym = {k: v for k, v in out["asym"].items() if binder_chain in k}
    out["ipsae_binder_best"] = max(binder_asym.values()) if binder_asym else 0.0
    out["ipsae_min"] = min(out["asym"].values()) if out["asym"] else 0.0
    out["binder_asym_pairs"] = binder_asym
    return out


def sw_align(query, ref):
    """Smith-Waterman alignment using parasail or pure-python fallback."""
    try:
        import parasail
        res = parasail.sw_trace_striped_16(query, ref, 11, 1, parasail.blosum62)
        tb = res.traceback
        q_al, r_al = tb.query, tb.ref
        aln_cols = len(q_al)
        n_matches = sum(1 for a, b in zip(q_al, r_al) if a == b and a != "-")
        q_aligned = sum(1 for a in q_al if a != "-")
        r_aligned = sum(1 for b in r_al if b != "-")
        return n_matches, q_aligned, r_aligned, aln_cols
    except Exception:
        # Fallback heuristic: exact k-mer matching or length comparison
        matches = sum(1 for a, b in zip(query, ref) if a == b)
        aln_cols = max(len(query), len(ref))
        return matches, len(query), len(ref), aln_cols


def evaluate_novelty(cand_name, cand_seq, known_db):
    """Evaluates novelty against known binders: FAIL if cov_ours >= 80% AND identity > 70%."""
    cand_len = len(cand_seq)
    if not known_db:
        return {
            "novelty_verdict": "PASS", "max_identity_pct": 0.0,
            "coverage_ours_pct": 0.0, "coverage_theirs_pct": 0.0,
            "worst_hit_id": "NO_DB", "novelty_verdict_strict": "PASS"
        }

    best_identity, best_cov_ours, best_cov_theirs, best_hit = 0.0, 0.0, 0.0, ""
    strict_best_id, strict_best_cov = 0.0, 0.0

    for kb_name, kb_seq in known_db:
        if not kb_seq: continue
        try:
            n_mat, q_aln, r_aln, aln_cols = sw_align(cand_seq, kb_seq)
        except Exception:
            continue
        if aln_cols == 0: continue

        cov_ours   = (q_aln / cand_len) * 100.0
        cov_theirs = (r_aln / len(kb_seq)) * 100.0
        cov_strict = max(cov_ours, cov_theirs)
        identity   = (n_mat / aln_cols) * 100.0

        if cov_ours >= 80.0:
            if identity > best_identity:
                best_identity   = identity
                best_cov_ours   = cov_ours
                best_cov_theirs = cov_theirs
                best_hit        = kb_name

        if cov_strict >= 80.0:
            if identity > strict_best_id:
                strict_best_id  = identity
                strict_best_cov = cov_strict

    verdict = "FAIL" if best_identity > 70.0 else "PASS"
    verdict_strict = "FAIL" if strict_best_id > 70.0 else "PASS"
    return {
        "novelty_verdict": verdict,
        "novelty_verdict_strict": verdict_strict,
        "max_identity_pct": round(best_identity, 2),
        "coverage_ours_pct": round(best_cov_ours, 2),
        "coverage_theirs_pct": round(best_cov_theirs, 2),
        "worst_hit_id": best_hit,
    }


def load_known_binders_db(target_name, known_dir=None):
    t_key = target_name.lower()
    files = KNOWN_BINDER_FILES.get(t_key, [])
    records = []
    if known_dir and os.path.exists(known_dir):
        for fn, label in files:
            p = os.path.join(known_dir, fn)
            if os.path.exists(p):
                rec = parse_fasta_records(p)
                records.extend(rec)
    if not records:
        records = CORE_KNOWN_BINDERS.get(t_key, [])
    return records


# ==============================================================================
# Main Screening Runner
# ==============================================================================

def main():
    parser = argparse.ArgumentParser(description="ESMFold2 Headless Candidate Screening Pipeline")
    parser.add_argument("--target", choices=["TL1A", "IL7RA", "ALL"], default="TL1A",
                        help="Target protein track to screen (default: TL1A)")
    parser.add_argument("--candidates-file", type=str, default=None,
                        help="Path to candidates CSV/FASTA file")
    parser.add_argument("--out-dir", type=str, default="/content/drive/MyDrive/proteinAI/out",
                        help="Directory to save output files and reports")
    parser.add_argument("--weights-dir", type=str, default="/content/drive/MyDrive/proteinAI/weights/ESMFold2",
                        help="Path to ESMFold2 weights directory")
    parser.add_argument("--msa-dir", type=str, default="/content/drive/MyDrive/proteinAI/data/msa",
                        help="Path to pre-cached MSAs directory")
    parser.add_argument("--known-binders-dir", type=str, default="/content/drive/MyDrive/proteinAI/data/known_binders",
                        help="Path to known binders directory")
    parser.add_argument("--ipsae-path", type=str, default="/content/drive/MyDrive/proteinAI/data/ipsae.py",
                        help="Path to ipsae.py script")
    parser.add_argument("--fast-mode", action="store_true",
                        help="Fast dry-run mode (10 loops, 50 sampling steps)")
    parser.add_argument("--max-candidates", type=int, default=None,
                        help="Maximum candidates to evaluate in this run")
    parser.add_argument("--hard-cap-min", type=float, default=300.0,
                        help="Hard cap timeout in minutes (default: 300.0)")
    parser.add_argument("--seed", type=int, default=0, help="Random seed (default: 0)")
    parser.add_argument("--esmc-precision", type=str, default="fp8", help="ESMC precision (default: fp8)")
    parser.add_argument("--kernel-backend", type=str, default="fused", help="Triton kernel backend (default: fused)")
    parser.add_argument("--no-checkpoint", action="store_true", help="Disable checkpoint resume")
    args, _ = parser.parse_known_args()

    t_start_all = time.time()
    def elapsed_min(): return (time.time() - t_start_all) / 60.0
    def cu_burned(): return elapsed_min() / 60.0 * 4.8  # ~4.8 CU/hr on L4

    # Resolve output directory
    out_dir = Path(args.out_dir)
    try:
        out_dir.mkdir(parents=True, exist_ok=True)
    except Exception:
        fallback = Path("./out_screening")
        fallback.mkdir(parents=True, exist_ok=True)
        print(f"Warning: Cannot write to {out_dir}. Falling back to {fallback}")
        out_dir = fallback

    runs_dir = out_dir / "screening_runs"
    runs_dir.mkdir(parents=True, exist_ok=True)

    partial_path = out_dir / "03_partial.json"
    csv_all_path = out_dir / "03_screening_results.csv"

    print("=" * 80)
    print("  ESMFOLD2 HEADLESS SCREENING PIPELINE (Changjiang Cup 2026)")
    print("=" * 80)
    print(f"  Target Track      : {args.target}")
    print(f"  Output Directory  : {out_dir}")
    print(f"  Weights Directory : {args.weights_dir}")
    print(f"  Hard Cap Timeout  : {args.hard_cap_min} minutes")
    print(f"  Sampling Config   : loops={10 if args.fast_mode else 20}, steps={50 if args.fast_mode else 200}, seed={args.seed}")
    print("=" * 80)

    # 1. Environment & Packages Check
    import torch
    if not torch.cuda.is_available():
        raise SystemExit("ERROR: CUDA GPU is not available! Must run on GPU runtime (e.g. L4).")

    gpu_name = torch.cuda.get_device_name(0)
    gpu_mem_gb = torch.cuda.get_device_properties(0).total_memory / 1024**3
    print(f"  GPU Detected      : {gpu_name} ({gpu_mem_gb:.2f} GB VRAM)")

    # Ensure required packages
    try:
        import esm
        from esm.models.esmfold2 import EsmFold2Model, ESMFold2InputBuilder, ProteinInput, StructurePredictionInput, MSA
        from esm.utils.msa.msa import is_a3m_insertion
    except ImportError:
        print("Installing esm==3.4.3 ...")
        subprocess.run([sys.executable, "-m", "pip", "install", "-q", "esm==3.4.3"], check=True)
        from esm.models.esmfold2 import EsmFold2Model, ESMFold2InputBuilder, ProteinInput, StructurePredictionInput, MSA
        from esm.utils.msa.msa import is_a3m_insertion

    try:
        import parasail
    except ImportError:
        print("Installing parasail ...")
        try:
            subprocess.run([sys.executable, "-m", "pip", "install", "-q", "parasail"], check=True)
            import parasail
        except Exception as e:
            print(f"Note: parasail installation bypassed ({e}), using fallback alignment.")

    # 2. Setup ipsae.py
    ipsae_py = args.ipsae_path
    if not ipsae_py or not os.path.exists(ipsae_py):
        search_paths = [
            "tools/ipsae/ipsae.py", "ipsae.py", "/content/ipsae.py",
            str(out_dir / "ipsae.py")
        ]
        ipsae_py = next((p for p in search_paths if os.path.exists(p) and os.path.getsize(p) > 5000), None)
    if not ipsae_py:
        ipsae_py = str(out_dir / "ipsae.py")
        print(f"Downloading ipsae.py to {ipsae_py} ...")
        url = "https://raw.githubusercontent.com/DunbrackLab/IPSAE/main/ipsae.py"
        req = urllib.request.Request(url, headers={"User-Agent": "proteinAI-changjiangcup/0.3"})
        with urllib.request.urlopen(req, timeout=60) as r, open(ipsae_py, "wb") as f:
            f.write(r.read())
    print(f"  ipSAE Engine      : {ipsae_py}")

    # 3. Resolve Candidate Pool
    candidate_file = args.candidates_file
    if not candidate_file or not os.path.exists(candidate_file):
        search_cands = [
            "/content/candidate_pool_for_screening.csv",
            "/content/proteinAI/notebooks/candidate_pool_for_screening.csv",
            "/content/drive/MyDrive/proteinAI/out/candidate_pool_for_screening.csv",
            "/content/drive/MyDrive/proteinAI/candidate_pool_for_screening.csv",
            "notebooks/candidate_pool_for_screening.csv",
            "data/candidate_pool_for_screening.csv",
            "candidate_pool_for_screening.csv",
        ]
        candidate_file = next((p for p in search_cands if os.path.exists(p)), None)

    if candidate_file:
        raw_candidates = parse_candidates_csv(candidate_file, default_target=args.target)
        print(f"  Candidate Pool    : Loaded {len(raw_candidates)} designs from {candidate_file}")
    else:
        raw_candidates = SAMPLE_BENCHMARKS
        print(f"  Candidate Pool    : Using {len(raw_candidates)} Benchmark Controls")

    # Apply Target Filtering
    if args.target != "ALL":
        candidates = [c for c in raw_candidates if c.get("target", "").upper() == args.target]
        print(f"  Target Filtered   : {len(candidates)} candidates selected for target '{args.target}'")
    else:
        candidates = raw_candidates

    if args.max_candidates and len(candidates) > args.max_candidates:
        candidates = candidates[:args.max_candidates]
        print(f"  Batch Capped      : {len(candidates)} candidates selected via --max-candidates")

    # Validate sequences
    for c in candidates:
        bad = sorted(set(c["sequence"]) - VALID_AA)
        if bad:
            raise ValueError(f"Candidate '{c['name']}' has invalid amino acid characters: {bad}")

    # 4. Checkpoint Resume: Scan previously completed candidates
    completed_by_name = {}
    if not args.no_checkpoint:
        if csv_all_path.exists():
            try:
                with open(csv_all_path, "r", encoding="utf-8") as fh:
                    for row in csv.DictReader(fh):
                        nm = row.get("name")
                        if nm and row.get("screen_verdict"):
                            completed_by_name[nm] = {
                                "name": nm, "target": row.get("target"),
                                "length": int(row["length"]) if row.get("length") else None,
                                "screen_verdict": row.get("screen_verdict"),
                                "ipsae_binder": float(row["ipsae_binder"]) if row.get("ipsae_binder") else 0.0,
                                "ipsae_min": float(row["ipsae_min"]) if row.get("ipsae_min") else 0.0,
                                "iptm": float(row["iptm"]) if row.get("iptm") else None,
                                "ptm": float(row["ptm"]) if row.get("ptm") else None,
                                "plddt_mean": float(row["plddt_mean"]) if row.get("plddt_mean") else None,
                                "novelty_verdict": row.get("novelty_verdict", "PASS"),
                                "max_identity_pct": float(row["max_identity_pct"]) if row.get("max_identity_pct") else 0.0,
                                "coverage_ours_pct": float(row["coverage_ours_pct"]) if row.get("coverage_ours_pct") else 0.0,
                                "worst_hit_id": row.get("worst_hit_id", ""),
                                "fold_seconds": float(row["fold_seconds"]) if row.get("fold_seconds") else 0.0,
                                "vram_peak_gb": float(row["vram_peak_gb"]) if row.get("vram_peak_gb") else 0.0,
                                "ok": True,
                            }
            except Exception as e:
                print(f"Warning reading existing results CSV: {e}")

        if partial_path.exists():
            try:
                with open(partial_path, "r", encoding="utf-8") as fh:
                    prev_p = json.load(fh)
                    for rec in prev_p.get("screened_candidates", []):
                        if rec.get("ok") and rec.get("name"):
                            completed_by_name[rec["name"]] = rec
            except Exception as e:
                print(f"Warning reading partial checkpoint JSON: {e}")

        if completed_by_name:
            print(f"\n[CHECKPOINT RESUME] Found {len(completed_by_name)} previously completed candidates on disk.")
            print("Completed candidates will be skipped during folding loops to preserve compute.\n")

    # 5. Load MSAs
    msa_dir = Path(args.msa_dir)
    if not msa_dir.exists():
        search_msa = [
            Path("/content/drive/MyDrive/proteinAI/data/msa"),
            Path("/content/proteinAI/data/msa"),
            Path("data/msa"),
        ]
        msa_dir = next((p for p in search_msa if p.exists()), msa_dir)

    def load_msa(target_key):
        meta = TARGET_METAS[target_key]
        p = msa_dir / f"{meta['name']}.uniref.a3m"
        if not p.exists():
            raise FileNotFoundError(f"MSA file missing: {p}")
        msa_obj = MSA.from_a3m(path=str(p))
        print(f"  Loaded MSA [{target_key}]: depth = {msa_obj.depth} from {p.name}")
        return msa_obj

    msas = {}
    needed_targets = {c.get("target", "TL1A").upper() for c in candidates}
    for tg in needed_targets:
        if tg in TARGET_METAS:
            msas[tg] = load_msa(tg)

    # 6. Load Known Binders DB
    known_db = {}
    for tg in needed_targets:
        known_db[tg] = load_known_binders_db(tg, args.known_binders_dir)
        print(f"  Known Binders DB [{tg}]: {len(known_db[tg])} reference sequences loaded")

    # 7. Initialize ESMFold2 Model
    weights_dir = args.weights_dir
    if not os.path.exists(weights_dir):
        alt_w = [
            "/content/drive/MyDrive/proteinAI/weights/ESMFold2",
            "/content/ESMFold2",
            "weights/ESMFold2"
        ]
        weights_dir = next((p for p in alt_w if os.path.exists(p)), weights_dir)

    print(f"\nLoading ESMFold2 Model from {weights_dir} (dtype=None, precision={args.esmc_precision}) ...")
    t_load = time.time()
    gc.collect(); torch.cuda.empty_cache(); torch.cuda.reset_peak_memory_stats()
    model = EsmFold2Model.from_pretrained(
        weights_dir, device="cuda", dtype=None,
        esmc_precision=args.esmc_precision, load_esmc=True,
    ).eval()

    try:
        model.set_kernel_backend(args.kernel_backend)
        print(f"Kernel backend: {args.kernel_backend} (activated)")
    except Exception as e:
        model.set_kernel_backend(None)
        print(f"Kernel backend fallback: {e}")

    builder = ESMFold2InputBuilder(ccd_cache=weights_dir)
    print(f"Model successfully loaded in {time.time() - t_load:.1f}s. VRAM allocated: {torch.cuda.memory_allocated()/1024**3:.2f} GB")

    # 8. Screening Loop
    num_loops = 10 if args.fast_mode else 20
    num_steps = 50 if args.fast_mode else 200
    lm_dropout = 0.3

    screened_results = []
    report_data = {
        "schema": 1,
        "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "target": args.target,
        "candidates_total": len(candidates),
        "screened_candidates": [],
    }

    def flush_checkpoint():
        report_data["elapsed_min"] = round(elapsed_min(), 2)
        report_data["cu_estimate"] = round(cu_burned(), 3)
        report_data["screened_candidates"] = screened_results
        tmp = str(partial_path) + ".tmp"
        with open(tmp, "w", encoding="utf-8") as f:
            json.dump(report_data, f, indent=2, default=str)
        os.replace(tmp, str(partial_path))

    print(f"\nStarting screening on {len(candidates)} candidate(s) ...\n")

    for idx, cand in enumerate(candidates, 1):
        c_name = cand["name"]
        c_seq  = cand["sequence"]
        c_tgt  = cand.get("target", "TL1A").upper()
        c_tag  = f"screen_{idx:03d}_{c_name[:24]}"
        cand_outdir = str(runs_dir / c_tag)

        # Checkpoint Resume
        if c_name in completed_by_name:
            rec = completed_by_name[c_name]
            print(f"[{idx:>3}/{len(candidates)}] RESUME (Cached): {c_name:<30} -> ipSAE={rec.get('ipsae_binder', 0.0):.4f} · {rec.get('screen_verdict')}")
            screened_results.append(rec)
            flush_checkpoint()
            continue

        # Hard cap budget check
        est_s = 40 if c_tgt == "IL7RA" else 110
        if elapsed_min() + (est_s / 60.0) > args.hard_cap_min:
            print(f"\n[HARD CAP REACHED] Time elapsed {elapsed_min():.1f} min >= {args.hard_cap_min} min. Stopping gracefully.")
            break

        print(f"\n{'='*75}")
        print(f"[{idx}/{len(candidates)}] Screening: {c_name} (Target: {c_tgt}, Length: {len(c_seq)} aa)")
        print(f"{'='*75}")

        if c_tgt == "IL7RA":
            chains = [("A", IL7RA_SEQ, msas["IL7RA"]), ("B", c_seq, None)]
            binder_chain = "B"
            total_res = len(IL7RA_SEQ) + len(c_seq)
        elif c_tgt == "TL1A":
            chains = [(["A", "B", "C"], TL1A_SEQ, msas["TL1A"]), ("Z", c_seq, None)]
            binder_chain = "Z"
            total_res = 3 * len(TL1A_SEQ) + len(c_seq)
        else:
            print(f"  [ERROR] Unknown target: {c_tgt}. Skipping.")
            continue

        spi = StructurePredictionInput(sequences=[
            ProteinInput(id=cid, sequence=seq, msa=msa) for cid, seq, msa in chains
        ])

        gc.collect(); torch.cuda.empty_cache(); torch.cuda.reset_peak_memory_stats()
        t0 = time.time()
        rec = {
            "name": c_name, "target": c_tgt, "length": len(c_seq),
            "total_complex_residues": total_res, "tag": c_tag,
            "binder_chain": binder_chain, "ok": False
        }

        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                result = builder.fold(
                    model, spi, num_loops=num_loops, num_sampling_steps=num_steps,
                    num_diffusion_samples=1, seed=args.seed, lm_dropout=lm_dropout,
                    msa_max_depth=1024, msa_column_mask_rate=0.1, complex_id=c_tag
                )
            fold_s = time.time() - t0
            peak_vram = torch.cuda.max_memory_allocated() / 1024**3

            rec["fold_seconds"] = round(fold_s, 1)
            rec["vram_peak_gb"] = round(peak_vram, 2)
            rec["ptm"]  = round(float(result.ptm), 4) if result.ptm is not None else None
            rec["iptm"] = round(float(result.iptm), 4) if result.iptm is not None else None

            if result.plddt is not None:
                pl = result.plddt.detach().cpu().float().numpy().ravel()
                rec["plddt_mean"] = round(float(pl.mean() * (100.0 if pl.max() <= 1.0 else 1.0)), 2)

            # ipSAE Adapter
            ipsae_res = run_ipsae_eval(result, c_tag, cand_outdir, binder_chain=binder_chain, ipsae_py=ipsae_py)
            rec["ipsae_binder"] = ipsae_res.get("ipsae_binder_best", 0.0)
            rec["ipsae_min"]    = ipsae_res.get("ipsae_min", 0.0)
            rec["ipsae_pairs"]  = ipsae_res.get("binder_asym_pairs", {})

            # Novelty Check
            nov_res = evaluate_novelty(c_name, c_seq, known_db.get(c_tgt, []))
            rec.update(nov_res)

            rec["ok"] = True
            rec["screen_verdict"] = "PASS" if (rec["novelty_verdict"] == "PASS" and rec["ipsae_binder"] >= 0.50) else "REVIEW"
            if rec["novelty_verdict"] == "FAIL":
                rec["screen_verdict"] = "FAIL_NOVELTY"

            print(f"  Fold: {fold_s:.1f}s | VRAM: {peak_vram:.2f} GB | pLDDT: {rec.get('plddt_mean')}% | ipTM: {rec['iptm']}")
            print(f"  ipSAE (binder): {rec['ipsae_binder']:.4f} | Pairs: {rec['ipsae_pairs']}")
            print(f"  Novelty: {rec['novelty_verdict']} (MaxID: {rec['max_identity_pct']}%, Cov: {rec['coverage_ours_pct']}%) -> {rec['screen_verdict']}")
            del result
        except torch.cuda.OutOfMemoryError as e:
            rec["error"] = "CUDA_OOM"
            print(f"  [ERROR] OOM on {total_res} residues: {e}")
        except Exception as e:
            rec["error"] = f"{type(e).__name__}: {e}"
            print(f"  [ERROR] {rec['error']}")
        finally:
            gc.collect(); torch.cuda.empty_cache()

        screened_results.append(rec)
        flush_checkpoint()

    # 9. Export Results & Rankings
    print("\n" + "=" * 80)
    print("Exporting Final Screening Reports & Rankings ...")
    print("=" * 80)

    def sort_key(r):
        verdict_rank = 0 if r.get("novelty_verdict") == "PASS" else 1
        ipsae = r.get("ipsae_binder") or 0.0
        iptm  = r.get("iptm") or 0.0
        return (verdict_rank, -ipsae, -iptm)

    sorted_candidates = sorted(screened_results, key=sort_key)
    csv_fields = [
        "rank", "name", "target", "length", "screen_verdict",
        "ipsae_binder", "ipsae_min", "iptm", "ptm", "plddt_mean",
        "novelty_verdict", "max_identity_pct", "coverage_ours_pct", "worst_hit_id",
        "fold_seconds", "vram_peak_gb", "error"
    ]

    rows_to_write = []
    for rk, c in enumerate(sorted_candidates, 1):
        rows_to_write.append({
            "rank": rk,
            "name": c.get("name"),
            "target": c.get("target"),
            "length": c.get("length"),
            "screen_verdict": c.get("screen_verdict"),
            "ipsae_binder": c.get("ipsae_binder"),
            "ipsae_min": c.get("ipsae_min"),
            "iptm": c.get("iptm"),
            "ptm": c.get("ptm"),
            "plddt_mean": c.get("plddt_mean"),
            "novelty_verdict": c.get("novelty_verdict"),
            "max_identity_pct": c.get("max_identity_pct"),
            "coverage_ours_pct": c.get("coverage_ours_pct"),
            "worst_hit_id": c.get("worst_hit_id"),
            "fold_seconds": c.get("fold_seconds"),
            "vram_peak_gb": c.get("vram_peak_gb"),
            "error": c.get("error", ""),
        })

    with open(csv_all_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=csv_fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows_to_write)
    print(f"  All Results CSV   : {csv_all_path} ({len(rows_to_write)} records)")

    top_pass = [r for r in rows_to_write if r.get("novelty_verdict") == "PASS"]
    csv_top_path = out_dir / "03_screening_top_candidates.csv"
    with open(csv_top_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=csv_fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(top_pass[:20])
    print(f"  Top Overall CSV   : {csv_top_path} ({min(len(top_pass), 20)} records)")

    for tgt_key in ("TL1A", "IL7RA"):
        tgt_cands = [r for r in top_pass if r.get("target", "").upper() == tgt_key][:20]
        if tgt_cands:
            tgt_csv = out_dir / f"03_screening_top_{tgt_key}.csv"
            with open(tgt_csv, "w", newline="", encoding="utf-8") as fh:
                w = csv.DictWriter(fh, fieldnames=csv_fields, extrasaction="ignore")
                w.writeheader()
                w.writerows(tgt_cands)
            print(f"  Top 20 {tgt_key} CSV  : {tgt_csv} ({len(tgt_cands)} records)")

    # JSON report
    report_data["finished_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    report_data["elapsed_min"]  = round(elapsed_min(), 2)
    report_data["cu_estimate"]  = round(cu_burned(), 3)
    report_data["summary_table"] = rows_to_write
    json_rep_path = out_dir / "03_screening_report.json"
    with open(json_rep_path, "w", encoding="utf-8") as fh:
        json.dump(report_data, fh, indent=2, default=str)
    print(f"  Report JSON       : {json_rep_path}")

    # ASCII Summary
    print("\n" + "=" * 96)
    print(f"{'Rk':<3} {'Candidate':<26} {'Tgt':<6} {'Len':>4} {'ipSAE':>8} {'ipTM':>6} {'pLDDT':>6} {'Novelty':<8} {'MaxID%':>7} {'Verdict'}")
    print("=" * 96)
    for r in rows_to_write[:25]:
        ipsae_str = f"{r['ipsae_binder']:.4f}" if r.get('ipsae_binder') is not None else "N/A"
        iptm_str  = f"{r['iptm']:.3f}" if r.get('iptm') is not None else "N/A"
        plddt_str = f"{r['plddt_mean']:.1f}" if r.get('plddt_mean') is not None else "N/A"
        id_str    = f"{r['max_identity_pct']:.1f}%" if r.get('max_identity_pct') is not None else "N/A"
        print(f"{r['rank']:<3} {r['name'][:25]:<26} {r['target']:<6} {r['length']:>4} "
              f"{ipsae_str:>8} {iptm_str:>6} {plddt_str:>6} {str(r['novelty_verdict']):<8} {id_str:>7} {r['screen_verdict']}")
    print("=" * 96)

    # Paste block for Claude
    print("\n========================================")
    print("PASTE CHO CLAUDE — BẮT ĐẦU")
    print("========================================")
    print(json.dumps({
        "script": "scripts/run_screening_colab.py",
        "target": args.target,
        "elapsed_min": report_data["elapsed_min"],
        "cu_estimate": report_data["cu_estimate"],
        "candidates_screened": len(screened_results),
        "top_candidates_count": len(top_pass),
        "top_5": rows_to_write[:5],
        "files": {
            "csv_all": str(csv_all_path),
            "csv_top": str(csv_top_path),
            "report_json": str(json_rep_path),
        }
    }, indent=2))
    print("========================================")
    print("PASTE CHO CLAUDE — KẾT THÚC")
    print("========================================")


if __name__ == "__main__":
    main()
