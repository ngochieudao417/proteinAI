"""ipSAE scorer — faithful reimplementation of DunbrackLab/IPSAE `ipsae.py` v4.

Why this exists: the official `ipsae.py` needs a `.pdb` + AF2-style `.json` pair and
shells out. The ipSAE score itself (`ipsae_d0res`, column `ipSAE` in ipsae.py's .txt
output) depends only on the PAE matrix and the per-token chain labels — inter-residue
distances enter pDockQ / LIS / the reported nres counters, never ipSAE. So a Mavenus
`*.metrics.json` (which carries `pae_matrix` + `pae_tokens`) can be scored directly.

Transcribed from ipsae.py v4 (MIT, Roland Dunbrack): ptm_func / calc_d0 / calc_d0_array
(~lines 110-135) and the per-residue + asym/max loops (~lines 726-850).

Calibration (see tools/test_ipsae_score.py): Mavenus `ctrl_7S5B` (7S5B binder 55 aa +
IL-7Ra 196 aa) scores ipSAE A->B 0.7638 / B->A 0.8138. Our own Notebook 03 run of the
same control on the full 219 aa construct gave 0.74884 / 0.772791 via the real ipsae.py.

CLI:  python tools/ipsae_score.py <metrics.json> [...]
"""
from __future__ import annotations

import json
import sys

import numpy as np

PAE_CUTOFF = 10.0
DIST_CUTOFF = 10.0  # recorded for provenance; unused by ipSAE_d0res


def ptm_func(x, d0):
    return 1.0 / (1.0 + (x / d0) ** 2.0)


def calc_d0(L, pair_type: str = "protein") -> float:
    """Yang & Skolnick d0, Proteins 57:702-710 (2004). Scalar form."""
    L = float(L)
    min_value = 2.0 if pair_type == "nucleic_acid" else 1.0
    d0 = 1.24 * (L - 15) ** (1.0 / 3.0) - 1.8 if L > 27 else 1.0
    return max(min_value, d0)


def calc_d0_array(L, pair_type: str = "protein"):
    """Vector form. Clamps L to >=26 then floors at 1.0, per ipsae.py 2026-01-03 fix."""
    L = np.maximum(np.array(L, dtype=float), 26.0)
    min_value = 2.0 if pair_type == "nucleic_acid" else 1.0
    return np.maximum(min_value, 1.24 * (L - 15) ** (1.0 / 3.0) - 1.8)


def ipsae_from_pae(pae, chains, pae_cutoff: float = PAE_CUTOFF) -> dict:
    """Score every ordered chain pair from a PAE matrix and per-token chain labels.

    Returns {(chain1, chain2): {ipsae, ipsae_d0chn, ipsae_d0dom, iptm_d0chn, n0res,
    d0res, n0chn, n0dom, nres1, nres2, asym_res_index}}, where `ipsae` is the value
    ipsae.py prints in its `ipSAE` column for the `asym` row of that ordered pair.
    """
    pae = np.asarray(pae, dtype=float)
    chains = np.asarray(chains)
    numres = len(chains)
    if pae.shape != (numres, numres):
        raise ValueError(f"pae shape {pae.shape} does not match {numres} chain labels")

    seen, uniq = set(), []
    for c in chains:
        if c not in seen:
            seen.add(c)
            uniq.append(c)

    out = {}
    for c1 in uniq:
        for c2 in uniq:
            if c1 == c2:
                continue
            m1, m2 = chains == c1, chains == c2
            n0chn = int(m1.sum() + m2.sum())
            d0chn = calc_d0(n0chn)
            ptm_d0chn = ptm_func(pae, d0chn)

            valid = np.outer(m1, m2) & (pae < pae_cutoff)
            rows1 = np.where(m1)[0]

            iptm_d0chn_byres = np.zeros(numres)
            ipsae_d0chn_byres = np.zeros(numres)
            set1, set2 = set(), set()
            for i in rows1:
                v = valid[i]
                if m2.any():
                    iptm_d0chn_byres[i] = ptm_d0chn[i, m2].mean()
                if v.any():
                    ipsae_d0chn_byres[i] = ptm_d0chn[i, v].mean()
                    set1.add(int(i))
                    set2.update(int(j) for j in np.where(v)[0])

            n0dom = len(set1) + len(set2)
            d0dom = calc_d0(n0dom)
            ptm_d0dom = ptm_func(pae, d0dom)

            n0res_byres = valid.sum(axis=1)
            d0res_byres = calc_d0_array(n0res_byres)

            ipsae_d0dom_byres = np.zeros(numres)
            ipsae_d0res_byres = np.zeros(numres)
            for i in rows1:
                v = valid[i]
                if not v.any():
                    continue
                ipsae_d0dom_byres[i] = ptm_d0dom[i, v].mean()
                ipsae_d0res_byres[i] = ptm_func(pae[i], d0res_byres[i])[v].mean()

            imax = int(np.argmax(ipsae_d0res_byres))
            out[(str(c1), str(c2))] = {
                "ipsae": float(ipsae_d0res_byres[imax]),
                "ipsae_d0chn": float(ipsae_d0chn_byres.max()),
                "ipsae_d0dom": float(ipsae_d0dom_byres.max()),
                "iptm_d0chn": float(iptm_d0chn_byres.max()),
                "n0res": int(n0res_byres[imax]),
                "d0res": float(d0res_byres[imax]),
                "n0chn": n0chn,
                "n0dom": n0dom,
                "nres1": len(set1),
                "nres2": len(set2),
                "asym_res_index": imax,
            }
    return out


def score_metrics_json(path: str, pae_cutoff: float = PAE_CUTOFF) -> tuple[dict, dict]:
    """Score a Mavenus `*.metrics.json` directly. Returns (metrics_dict, pair_scores)."""
    with open(path, encoding="utf-8") as fh:
        d = json.load(fh)
    chains = [t["chain_id"] for t in d["pae_tokens"]]
    return d, ipsae_from_pae(d["pae_matrix"], chains, pae_cutoff)


def summarise(pair_scores: dict) -> dict:
    """Collapse per-pair scores to the two numbers the screening gate uses."""
    vals = [v["ipsae"] for v in pair_scores.values()]
    return {"ipsae_max": max(vals) if vals else 0.0,
            "ipsae_min": min(vals) if vals else 0.0}


def main(argv: list[str]) -> int:
    if not argv:
        print(__doc__)
        return 2
    for p in argv:
        d, res = score_metrics_json(p)
        print(f"\n{d.get('name', p)}  iptm={d.get('iptm')}  ptm={d.get('ptm')}  "
              f"plddt={d.get('plddt_0_100')}")
        for (c1, c2), v in res.items():
            print(f"  {c1}->{c2} asym  ipSAE={v['ipsae']:.6f}  "
                  f"d0chn={v['ipsae_d0chn']:.6f}  d0dom={v['ipsae_d0dom']:.6f}  "
                  f"n0res={v['n0res']}  d0res={v['d0res']:.2f}")
        s = summarise(res)
        print(f"  ipSAE_max={s['ipsae_max']:.6f}  ipSAE_min={s['ipsae_min']:.6f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
