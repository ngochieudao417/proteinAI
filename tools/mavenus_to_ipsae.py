"""Adapter: Mavenus ESMFold2 output -> the (.pdb, .json) pair `ipsae.py` reads (AF2 path).

Satisfies docs/briefs/AGY-03-ipsae-adapter.md. This module only writes inputs; it never
runs ipsae.py. For scoring without the external script, see tools/ipsae_score.py, which
reimplements the ipSAE formula and is calibrated against the 7S5B control.

Mavenus delivers, per design `<name>`:
  <root>/metrics_json/<name>.metrics.json   pae_matrix, pae_tokens, iptm, ptm, plddt_0_100
  <root>/structures/<name>.cif              minimal mmCIF, one _atom_site loop_,
                                            B_iso_or_equiv column carries per-atom pLDDT (0-100)

`pae_tokens` order IS the row/column order of `pae_matrix`; the .cif residue order is
asserted against it rather than assumed.

CLI:  python tools/mavenus_to_ipsae.py <outdir> [name ...]     (no name = all designs)
"""
from __future__ import annotations

import glob
import json
import os
import sys

DEFAULT_ROOT = "Mavenus_output/IL7Ra_minibinder_FINAL"
REQUIRED_CIF_COLS = ("group_PDB", "type_symbol", "label_atom_id", "label_comp_id",
                     "auth_asym_id", "auth_seq_id", "B_iso_or_equiv",
                     "Cartn_x", "Cartn_y", "Cartn_z")


def _parse_cif_atom_site(path: str) -> tuple[dict, list[list[str]]]:
    """Read the single `_atom_site` loop_. Column indices come from the header, never hardcoded."""
    cols: list[str] = []
    rows: list[list[str]] = []
    in_loop = False
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            s = line.strip()
            if s.startswith("_atom_site."):
                cols.append(s.split(".", 1)[1])
                in_loop = True
                continue
            if in_loop and (s.startswith("ATOM") or s.startswith("HETATM")):
                rows.append(s.split())
            elif in_loop and s.startswith("#") and rows:
                break
    idx = {c: i for i, c in enumerate(cols)}
    missing = [c for c in REQUIRED_CIF_COLS if c not in idx]
    if missing:
        raise ValueError(f"{path}: _atom_site loop is missing columns {missing}")
    return idx, rows


def load_mavenus(name: str, root: str = DEFAULT_ROOT) -> dict:
    """Load one design. Returns the parsed metrics plus residue-grouped atoms in PAE order.

    Keys: name, metrics, chains, n_res, pae, iptm, ptm, residues (list of
    {chain, resnum, resname, plddt, atoms}), where `residues` is in pae_matrix order.
    """
    mpath = os.path.join(root, "metrics_json", f"{name}.metrics.json")
    cpath = os.path.join(root, "structures", f"{name}.cif")
    with open(mpath, encoding="utf-8") as fh:
        met = json.load(fh)
    idx, rows = _parse_cif_atom_site(cpath)

    grouped: dict[tuple[str, int], list[list[str]]] = {}
    order: list[tuple[str, int]] = []
    for r in rows:
        key = (r[idx["auth_asym_id"]], int(r[idx["auth_seq_id"]]))
        if key not in grouped:
            grouped[key] = []
            order.append(key)
        grouped[key].append(r)

    tok_order = [(t["chain_id"], int(t["residue_index"])) for t in met["pae_tokens"]]
    if order != tok_order:
        raise ValueError(
            f"{name}: .cif residue order does not match pae_tokens order "
            f"(cif n={len(order)} head={order[:3]}, tokens n={len(tok_order)} head={tok_order[:3]})"
        )

    residues = []
    for key in order:
        ats = grouped[key]
        ca = next((a for a in ats if a[idx["label_atom_id"]].strip().upper() == "CA"), None)
        if ca is not None:
            plddt = float(ca[idx["B_iso_or_equiv"]])
        else:
            plddt = sum(float(a[idx["B_iso_or_equiv"]]) for a in ats) / len(ats)
        residues.append({
            "chain": key[0], "resnum": key[1],
            "resname": ats[0][idx["label_comp_id"]][:3],
            "plddt": plddt,
            "atoms": [{
                "name": a[idx["label_atom_id"]].strip().upper(),
                "element": a[idx["type_symbol"]].strip().upper()[:2],
                "x": float(a[idx["Cartn_x"]]), "y": float(a[idx["Cartn_y"]]),
                "z": float(a[idx["Cartn_z"]]), "b": float(a[idx["B_iso_or_equiv"]]),
            } for a in ats],
        })

    pae = met["pae_matrix"]
    if len(pae) != len(residues):
        raise ValueError(f"{name}: pae is {len(pae)}x{len(pae)} but .cif has {len(residues)} residues")

    return {"name": name, "metrics": met, "chains": met["chains"], "n_res": len(residues),
            "pae": pae, "iptm": met.get("iptm", -1.0), "ptm": met.get("ptm", -1.0),
            "residues": residues}


def write_ipsae_inputs(name: str, outdir: str, root: str = DEFAULT_ROOT) -> dict:
    """Write <outdir>/<name>.pdb and <outdir>/<name>.json. Returns paths and counts."""
    d = load_mavenus(name, root)
    os.makedirs(outdir, exist_ok=True)
    pdb_path = os.path.join(outdir, f"{name}.pdb")
    json_path = os.path.join(outdir, f"{name}.json")

    lines, serial = [], 1
    for res in d["residues"]:
        ch = res["chain"][:1] or "A"
        for a in res["atoms"]:
            an = a["name"]
            an4 = f" {an:<3}" if len(an) < 4 else an[:4]
            lines.append(
                f"ATOM  {serial:5d} {an4} {res['resname']:>3s} {ch:1s}{res['resnum']:4d}    "
                f"{a['x']:8.3f}{a['y']:8.3f}{a['z']:8.3f}{1.00:6.2f}{a['b']:6.2f}"
                f"          {a['element']:>2s}"
            )
            serial += 1
    with open(pdb_path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines + ["TER", "END"]) + "\n")

    plddt = [round(r["plddt"], 2) for r in d["residues"]]
    with open(json_path, "w", encoding="utf-8") as fh:
        json.dump({"pae": d["pae"], "plddt": plddt,
                   "iptm": float(d["iptm"]), "ptm": float(d["ptm"])}, fh)

    return {"pdb_path": pdb_path, "json_path": json_path, "n_res": d["n_res"],
            "n_atoms": serial - 1, "chains": [(c["id"], c["len"]) for c in d["chains"]],
            "plddt_mean": sum(plddt) / len(plddt)}


def all_names(root: str = DEFAULT_ROOT) -> list[str]:
    return sorted(os.path.basename(p)[: -len(".metrics.json")]
                  for p in glob.glob(os.path.join(root, "metrics_json", "*.metrics.json")))


def main(argv: list[str]) -> int:
    if not argv:
        print(__doc__)
        return 2
    outdir, names = argv[0], argv[1:] or all_names()
    for n in names:
        r = write_ipsae_inputs(n, outdir)
        print(f"{n:<14} {r['n_res']:>4} res  {r['n_atoms']:>6} atoms  "
              f"plddt_mean={r['plddt_mean']:6.2f}  chains={r['chains']}")
    print(f"\n{len(names)} design(s) -> {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
