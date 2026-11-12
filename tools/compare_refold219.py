"""Compare the 219 aa official-construct re-fold against the 196 aa Mavenus scoring.

Input A: data/04_ipsae_IL7RA_recomputed.csv — ipSAE recomputed from Mavenus PAE matrices,
         receptor chain A = 196 aa (7OPB 36-231).
Input B: the Colab run's 03_screening_results.csv — same designs folded against the
         official construct, receptor chain A = 219 aa (P16871 Glu21-Asp239), MSA depth 921.

Prints the per-design delta, the control-anchored sanity check, and the final top 20 by the
219 aa numbers. Writes the merged table and a ready-to-submit CSV.

Usage: python tools/compare_refold219.py <03_screening_results.csv> [out_merged.csv] [out_submission.csv]
"""
from __future__ import annotations

import csv
import os
import sys

OLD = "data/04_ipsae_IL7RA_recomputed.csv"
CTRL_196 = 0.763802   # ctrl_7S5B scored from Mavenus PAE, 196 aa receptor
CTRL_219 = 0.748840   # ctrl_7S5B measured in Notebook 03 with the real ipsae.py, 219 aa


def norm(name: str) -> str:
    """c29_A17_5 / rank01_c29_... / c29 -> c29"""
    for part in str(name).split("_"):
        p = part.strip().lower()
        if p.startswith("c") and p[1:].isdigit():
            return f"c{int(p[1:]):02d}"
    return str(name).strip()


def main(argv: list[str]) -> int:
    if not argv:
        print(__doc__)
        return 2
    new_csv = argv[0]
    out_merged = argv[1] if len(argv) > 1 else "data/06_ipsae_IL7RA_219aa_merged.csv"
    out_sub = argv[2] if len(argv) > 2 else "data/submission_IL7RA_v3_219aa.csv"

    old = {norm(r["tag"]): r for r in csv.DictReader(open(OLD, encoding="utf-8-sig"))}
    new = {}
    for r in csv.DictReader(open(new_csv, encoding="utf-8-sig")):
        if (r.get("error") or "").strip():
            print(f"  [warn] {r.get('name')} errored: {r['error'][:80]}")
            continue
        new[norm(r["name"])] = r

    rows = []
    for tag, n in new.items():
        o = old.get(tag)
        if not o:
            print(f"  [warn] {tag} folded but absent from {OLD}")
            continue
        a = float(o["ipsae_min"])
        b = float(n["ipsae_min"])
        rows.append({
            "tag": tag, "backbone": o["backbone"], "length": o["length"],
            "ipsae_min_196": round(a, 6), "ipsae_min_219": round(b, 6),
            "delta": round(b - a, 6),
            "ipsae_binder_219": n.get("ipsae_binder", ""),
            "iptm_219": n.get("iptm", ""), "plddt_219": n.get("plddt_mean", ""),
            "verdict_219": n.get("screen_verdict", ""),
            "dist_N29": o["dist_N29"], "steric_atoms_lt2A": o["steric_atoms_lt2A"],
            "mavenus_pass_all": o["mavenus_pass_all"],
            "in_submission_v1": o["in_submission_v1"],
            "sequence": o["sequence"],
        })

    if not rows:
        print("no designs matched; nothing to do")
        return 1

    rows.sort(key=lambda r: -r["ipsae_min_219"])
    for i, r in enumerate(rows, 1):
        r["rank_219"] = i

    cols = ["rank_219", "tag", "backbone", "length", "ipsae_min_219", "ipsae_min_196", "delta",
            "ipsae_binder_219", "iptm_219", "plddt_219", "verdict_219", "dist_N29",
            "steric_atoms_lt2A", "mavenus_pass_all", "in_submission_v1", "sequence"]
    with open(out_merged, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader(); w.writerows(rows)

    deltas = [r["delta"] for r in rows]
    mean_d = sum(deltas) / len(deltas)
    print(f"\n{len(rows)} design folded on the 219 aa construct\n")
    hdr = f"{'rank':<6}{'tag':<6}{'bb':<8}{'219aa':<9}{'196aa':<9}{'delta':<9}{'v1':<4}"
    print(hdr); print("-" * len(hdr))
    for r in rows:
        print(f"{r['rank_219']:<6}{r['tag']:<6}{r['backbone']:<8}{r['ipsae_min_219']:<9.4f}"
              f"{r['ipsae_min_196']:<9.4f}{r['delta']:<+9.4f}{r['in_submission_v1'] or '-':<4}")

    print(f"\ndelta 219aa - 196aa : mean {mean_d:+.4f}  min {min(deltas):+.4f}  max {max(deltas):+.4f}")
    print(f"control anchor      : 7S5B was {CTRL_196:.4f} at 196 aa and {CTRL_219:.4f} at 219 aa, "
          f"i.e. {CTRL_219 - CTRL_196:+.4f}")
    print(f"designs >= 0.75     : {sum(1 for r in rows if r['ipsae_min_219'] >= 0.75)}")
    print(f"designs >= 0.30     : {sum(1 for r in rows if r['ipsae_min_219'] >= 0.30)}")
    print(f"designs == 0.00     : {sum(1 for r in rows if r['ipsae_min_219'] == 0.0)}")

    top = rows[:20]
    with open(out_sub, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\r\n")
        w.writerow(["name", "sequence"])
        for i, r in enumerate(top, 1):
            nm = f"r{i:02d}_{r['tag']}_{r['backbone']}_L{r['length']}_ipSAE{r['ipsae_min_219']:.3f}"
            w.writerow([nm[:64], r["sequence"]])

    gly = [r["tag"] for r in top if float(r["dist_N29"]) < 8.0]
    print(f"\nwrote {out_merged}")
    print(f"wrote {out_sub}  (top 20 by the 219 aa ipSAE)")
    print(f"of which glycan N29 < 8 A: {gly or 'none'}")
    print("\nNEXT: run tools/novelty_check.py on the new file before treating it as final.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
