# IL-7Ra / CD127 minibinder — FINAL delivery

## Pipeline
RFdiffusion (40 backbones, 85-100 aa, 7OPB chain A + hotspots A97-A103)
-> ProteinMPNN (6 seq/backbone, top-12 backbones) -> hard sequence filters
-> ESMFold2 complex fold (1 seed / 1 sample / no template) -> structural hard-constraint validation.

## Files
- FINAL_ranked_37.csv    all 37 designs: pass_all first, then by ipTM
- SUBMISSION_top20.csv   the 20 designs that pass ALL hard constraints, ranked by ipTM
- SUBMISSION_top20.fasta same 20 as FASTA
- structures/*.cif       38 predicted complexes (chain A = receptor 196 aa, chain B = binder)
- metrics_json/*.json    ESMFold2 metrics + pae_matrix + pae_tokens (for local ipSAE)
- tables/*.csv           all intermediate QC tables

## Hard constraints checked on the folded all-atom models
1. Epitope engagement: >=3 of residues A97-A103 within 4.5 A  -> all 37 pass
2. Glycosylation clearance: >=8.0 A from N29 and N45 side chain -> 5 designs FAIL
3. Steric sanity: no inter-chain heavy-atom pair < 2.0 A -> 13 designs FAIL (ESMFold2 strain)
4. Free cysteine: 0 -> all 37 pass
=> 20 / 37 pass all four.

## Reference points
- 7S5B control, same pipeline: pLDDT 91.45 / pTM 0.929 / ipTM 0.937 / iface PAE 2.17 A
- Best passing design c29 (A17-5, 98 aa): pLDDT 91.48 / pTM 0.934 / ipTM 0.933 / iface PAE 2.17 A

## Important notes
1. NUMBERING: ESMFold2 renumbers both chains from 1. Receptor chain A 1-196 = 7OPB 36-231
   (offset +35). 7OPB numbering = construct/3DI3 numbering + 20.
2. ipTM / PAE are ESMFold2 values, not the official Dunbrack ipSAE. A local d0res proxy did
   not reproduce the reported 7S5B ipSAE of 0.7483, so ipSAE is left to your local ipsae.py.
3. c20 (ipTM 0.884) and c27 (ipTM 0.877) are top-2 by ipTM but FAIL the N29 glycan clearance
   (4.9 A and 5.8 A). They are kept in FINAL_ranked_37.csv with pass_all=0, not in the top-20.
4. All models are computational predictions; no experimental validation.
