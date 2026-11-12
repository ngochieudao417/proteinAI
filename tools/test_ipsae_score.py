"""Tests for tools/ipsae_score.py and tools/mavenus_to_ipsae.py.

Self-test answers are the ones fixed in docs/briefs/AGY-03-ipsae-adapter.md plus the
7S5B positive control, whose ipSAE our Notebook 03 measured with the real ipsae.py.
"""
import json
import math
import os
import subprocess
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ipsae_score as S
import mavenus_to_ipsae as A

ROOT = A.DEFAULT_ROOT
FAIL = []


def check(cond, msg):
    print(("  PASS  " if cond else "  FAIL  ") + msg)
    if not cond:
        FAIL.append(msg)


print("1. d0 formulas match Yang & Skolnick / ipsae.py")
check(math.isclose(S.calc_d0(26), 1.0), "calc_d0(26) floors at 1.0")
check(math.isclose(S.calc_d0(28), 1.24 * 13 ** (1 / 3) - 1.8),
      "calc_d0(28) uses the cube-root branch")
check(math.isclose(float(S.calc_d0_array(0)), 1.0), "calc_d0_array(0) clamps to 1.0, not 1.04")
check(math.isclose(float(S.calc_d0_array(55)), 1.24 * 40 ** (1 / 3) - 1.8), "calc_d0_array(55)")
check(math.isclose(S.ptm_func(0.0, 2.0), 1.0), "ptm_func(0) == 1")
check(math.isclose(S.ptm_func(2.0, 2.0), 0.5), "ptm_func(d0) == 0.5")

print("\n2. a perfect interface scores ~1, a hopeless one scores 0")
ch = np.array(["A"] * 5 + ["B"] * 5)
check(S.ipsae_from_pae(np.zeros((10, 10)), ch)[("A", "B")]["ipsae"] == 1.0,
      "all-zero PAE -> ipSAE 1.0")
check(S.ipsae_from_pae(np.full((10, 10), 30.0), ch)[("A", "B")]["ipsae"] == 0.0,
      "all-PAE-30 (above cutoff) -> ipSAE 0.0")

print("\n3. token counts fixed in the AGY-03 brief")
for name, n_res, blen in (("ctrl_7S5B", 251, 55), ("c29_A17_5", 294, 98)):
    d = A.load_mavenus(name, ROOT)
    check(d["n_res"] == n_res, f"{name}: {d['n_res']} residues == {n_res}")
    check(len(d["pae"]) == n_res and len(d["pae"][0]) == n_res, f"{name}: pae is {n_res}x{n_res}")
    check(d["chains"][0]["len"] == 196, f"{name}: receptor chain A is 196 aa")
    check(d["chains"][1]["len"] == blen, f"{name}: binder chain B is {blen} aa")

print("\n4. adapter output is self-consistent")
with tempfile.TemporaryDirectory() as td:
    r = A.write_ipsae_inputs("c29_A17_5", td, ROOT)
    j = json.load(open(r["json_path"], encoding="utf-8"))
    check(set(j) == {"pae", "plddt", "iptm", "ptm"}, "json has exactly the 4 AF2-path keys")
    check(len(j["plddt"]) == r["n_res"], "plddt is per-residue, length N")
    check(len(j["pae"]) == r["n_res"], "pae is N x N")
    check(abs(r["plddt_mean"] - 91.4768) < 0.01,
          f"plddt_mean {r['plddt_mean']:.4f} reproduces metrics plddt_0_100 91.4768")
    nat = sum(1 for ln in open(r["pdb_path"], encoding="utf-8") if ln.startswith("ATOM"))
    check(nat == r["n_atoms"], f"pdb holds all {r['n_atoms']} atoms")
    tail = [ln for ln in open(r["pdb_path"], encoding="utf-8") if ln.startswith("ATOM")][-1]
    check(tail[21] == "B" and int(tail[22:26]) == 98, "last ATOM is chain B residue 98")

print("\n5. 7S5B positive control reproduces the ipSAE measured with the real ipsae.py")
d, res = S.score_metrics_json(os.path.join(ROOT, "metrics_json", "ctrl_7S5B.metrics.json"))
ab, ba = res[("A", "B")]["ipsae"], res[("B", "A")]["ipsae"]
print(f"        A->B {ab:.6f}   B->A {ba:.6f}   (Notebook 03, 219 aa receptor: 0.74884 / 0.772791)")
check(0.70 < ab < 0.82, f"A->B {ab:.4f} within 0.07 of the Notebook 03 value")
check(0.74 < ba < 0.86, f"B->A {ba:.4f} within 0.07 of the Notebook 03 value")
check(S.summarise(res)["ipsae_min"] == min(ab, ba), "summarise picks the min asym pair")

print("\n6. the chain-swap is not symmetric, and the CLI runs")
check(abs(ab - ba) > 1e-6, "A->B differs from B->A, as ipSAE is asymmetric")
p = subprocess.run([sys.executable, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                 "ipsae_score.py"),
                    os.path.join(ROOT, "metrics_json", "c29_A17_5.metrics.json")],
                   capture_output=True, text=True)
check(p.returncode == 0 and "ipSAE_min=0.851" in p.stdout, "CLI prints c29 ipSAE_min 0.851")

print("\n" + "=" * 60)
print(f"{'ALL TESTS PASS' if not FAIL else str(len(FAIL)) + ' FAILURE(S): ' + '; '.join(FAIL)}")
raise SystemExit(1 if FAIL else 0)
