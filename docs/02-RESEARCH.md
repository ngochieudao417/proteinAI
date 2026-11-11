# 02-RESEARCH — 长江杯 2026, minibinder tracks (TL1A + IL-7Rα)

> Research agent output. Date: **2026-11-11**. Deadline **2026-11-30** → **19 days**.
> Scope: Q1 ESMFold2 · Q2 ipSAE · Q3 TL1A · Q4 IL-7Rα · Q5 pipeline · Q6 gate exploitation.
> Every factual claim carries a URL. Anything I could not verify is marked **UNCONFIRMED**.

---

## ╔══ BOTTOM LINE — 5 most actionable conclusions ══╗

1. **"ESMFold2" is SOLVED. It is Biohub's ESMFold2 (released 2026-05-27), built on a frozen ESMC-6B language model + 48 looped folding layers + a diffusion decoder. 6.58 B params, MIT licence, weights public on Hugging Face (`biohub/ESMFold2`), code `pip install esm` (`Biohub/esm`).** It is NOT Lin et al. 2023 ESMFold and NOT ESM3. It takes optional per-chain MSAs (`ProteinInput(..., msa=MSA.from_a3m(...))`), does all-atom multimer, and returns `result.pae`, `result.plddt`, `result.ptm`, `result.iptm`. **The `-Fast` checkpoint silently ignores MSAs** → because the gate uses an MSA for the target, the gate must be running the **full** `biohub/ESMFold2`.
   - Sources: [HF model card](https://huggingface.co/biohub/ESMFold2) · [Biohub/esm README](https://github.com/Biohub/esm) · [preprint, Candido et al. 2026](https://www.biorxiv.org/content/10.64898/2026.06.03.729735v1.full) · [release announcement](https://biohub.org/news/world-model-of-protein-biology/)

2. **The L4 24 GB is a hard constraint, and TL1A-as-a-trimer probably does not fit.** Biohub maintainer `ebetica`: 24 GB VRAM handles "**complexes with ~500 total amino acids**"; H100-80 GB ≈ 2000 residues; **no CPU offload planned**. IL-7Rα(219) + binder(40–100) = **259–319 aa → fits**. TL1A trimer (3×180=540) + binder = **580–640 aa → does NOT fit** at stock settings. Mitigations: `dtype=bfloat16` (weights are stored F32, 26.3 GB — a naïve fp32 load OOMs *before inference*), fp8 ESMC backbone added in `esm` v3.4.2, and PR #385 `low_memory_mode=True` (ESMC CPU offload) which is **still open, not merged**.
   - Sources: [issue #320](https://github.com/Biohub/esm/issues/320) · [PR #385](https://github.com/Biohub/esm/issues/385) · [v3.4.2 release notes](https://github.com/Biohub/esm/releases) · [config.json `"dtype":"float32"`](https://huggingface.co/biohub/ESMFold2/raw/main/config.json)

3. **Target ipSAE ≥ 0.61, and treat 0.65+ as the competitive zone.** Three independent lines converge: (a) a meta-analysis of **3,766 experimentally characterised binders** recommends **AF3 `ipSAE_min` > 0.61** as the single best filter (1.4× average-precision gain over AF2 `pae_interaction`); (b) a competitor's public ESMFold2 run on 1,127 VHH designs has `ipsae_min` **p99 = 0.639, max = 0.698, median = 0.000**; (c) AsterFire's published TL1A numbers are 0.637 / 0.628. Minibinders score materially higher than VHH/scFv on ESMFold2 (paper hit rates 53.8–70 % vs 12.1–21.0 %), so for our two tracks aim **≥ 0.70** to be safe.
   - Sources: [meta-analysis](https://www.biorxiv.org/content/10.1101/2025.08.14.670059v1.full) · [competitor raw data](https://github.com/RFjason/binder-design/blob/HEAD/campaigns/tnfr1-vhh-20261103/results/all_candidates.csv) · [ESMFold2 preprint](https://www.biorxiv.org/content/10.64898/2026.06.03.729735v1.full)

4. **Pipeline recommendation: BoltzGen as primary generator → ESMFold2+ipSAE as the selector, with BindCraft as the fallback. Do NOT attempt the official ESMFold2-inversion binder design on Colab** — the official `binder_design.py` reports **51 GB VRAM, 27 GB even with `REUSE_ESMC=True`**, and runs `@app.cls(gpu="H100")`. BoltzGen gives 66–80 % nanomolar hit rates **at ≤20 designs per target** — exactly our submission budget — MIT licence, weights public. But BoltzGen itself wants ≥24 GB and "**a 24 GB card will OOM on anything but small systems**", so budget for a rented L40S/A100 for the generation step, or use the sponsor platforms (MatwingsVenus / AsterFire) which host BindCraft + BoltzGen + ESMFold2 + ipSAE for free during the contest.
   - Sources: [binder_design.py](https://github.com/Biohub/esm/blob/main/cookbook/tutorials/binder_design.py) · [BoltzGen paper](https://www.biorxiv.org/content/10.1101/2025.11.20.689494v2.full) · [BoltzGen repo](https://github.com/HannesStark/boltzgen) · [BindCraft paper](https://www.nature.com/articles/s41586-025-09429-6)

5. **Two epitope decisions, already resolved from structure.**
   - **TL1A**: the receptor (DcR3/DR3) site is a 25-residue patch that **straddles the groove between two trimer protomers** (I computed this from 3K51 assembly + 2RE9). In **2RE9 / submitted-construct numbering** the hotspots are **Y117, E119, P120, T121, D104, S105** (upper patch, groove-facing) and **N136, W137, F138, K169, E170** (lower patch, adjacent-protomer-facing). UniProt-187 and UniProt-190 are annotated by UniProt as "important for binding TNFRSF6B" → construct **S116, E119**.
   - **IL-7Rα**: the IL-7 site is in **domain D1** (FN3/D2 = UniProt 131–231). In **3DI3-chain-B / submitted-construct numbering**: **K77, K78, F79, L80, L81, I82, G83** (dominant), plus **S31, H33, L57, V58, E59**, plus D2 fringe **K138, Y139, Y192, F193**. This overlaps the clinical anti-CD127 "site-1/site-2b" epitope (reported residues H33, E75, F79, I82, K84, M144, R186, H191, Y192).
   - **Novelty risk is asymmetric**: there is a published de novo IL-7Rα minibinder (PDB **7S5B**, 55 aa, sequence below) — we must stay ≤70 % identical to it. I found **no published TL1A minibinder**, so TL1A novelty is low-risk.
   - Sources: [2RE9](https://www.rcsb.org/structure/2RE9) · [3K51](https://www.rcsb.org/structure/3K51) · [3DI3](https://www.rcsb.org/structure/3DI3) · [UniProt O95150](https://rest.uniprot.org/uniprotkb/O95150.json) · [UniProt P16871](https://rest.uniprot.org/uniprotkb/P16871.json) · [7S5B](https://www.rcsb.org/structure/7S5B)

## ╚═══════════════════════════════════════════════╝

---

# Q1 — What exactly is "ESMFold2"?

## 1.1 Identity (CONFIRMED)

| Field | Value | Source |
|---|---|---|
| Name | **ESMFold2** (and `ESMFold2-Fast`, `ESMFold2-Experimental`, `*-Cutoff2025`) | [HF biohub org](https://huggingface.co/api/models?author=biohub) |
| Releaser | **Biohub** (501(c)(3); the former EvolutionaryScale ESM team — repo still links `evolutionaryscale/wheels`, API "migrated from forge.evolutionaryscale.ai to biohub.ai") | [biohub.org news](https://biohub.org/news/world-model-of-protein-biology/), [README](https://github.com/Biohub/esm) |
| Public release | **2026-05-27** (HF repo created 2026-05-13) | [news](https://biohub.org/news/world-model-of-protein-biology/), HF API `createdAt` |
| Preprint | *Language Modeling Materializes a World Model of Protein Biology*, Candido, Hayes, Derry, Rao, Lin, Verkuil, … Rives (2026). DOI **10.64898/2026.06.03.729735** | [bioRxiv](https://www.biorxiv.org/content/10.64898/2026.06.03.729735v1.full) |
| Params | **6,579,755,843** (stored F32; 26.3 GB on disk) | [HF API](https://huggingface.co/api/models/biohub/ESMFold2) |
| Licence | **MIT** | [README §Licenses](https://github.com/Biohub/esm) |
| Weights public? | **Yes, ungated** (`gated: false`), 352 k downloads | HF API |
| Repo | `https://github.com/Biohub/esm` — PyPI `pip install esm` (`esm` v3.4.3, 2026-11-09) | [releases](https://github.com/Biohub/esm/releases) |
| Also in HF Transformers | ≥ **v5.16.0** (`EsmFold2Model`) | README |

**It is definitively not plain ESMFold.** Architecture per the preprint: representations pooled from **all layers of a frozen ESMC-6B**, compressed into a pair representation by a 4-layer encoder, iteratively updated by **48 folding layers** (24 in `-Fast`) with a stabilised recurrent update, then a **diffusion transformer with sliding-window atom attention** denoising to all-atom coordinates. Handles protein, DNA, RNA, ligands (SMILES or CCD), modified residues and covalent bonds. `config.json` confirms `folding_trunk_num_hidden_layers: 48`, `num_loops: 20` (default), `sliding_window: 128`.

The participant's 2026-11-05 warning is correct and the distinction matters enormously: ESMFold (2023) is MSA-free single-chain-only; ESMFold2 is an AF3-class all-atom cofolding model that **accepts MSAs** and beats AlphaFold3 on FoldBench protein–protein and antibody–antigen.

## 1.2 MSA support — the exact mechanism the gate is using (CONFIRMED)

MSA is an **optional per-chain** input. This maps 1:1 onto the competition's "MSA only for the target sequence, no MSA for the designed protein":

```python
from esm.models.esmfold2 import ESMFold2InputBuilder, EsmFold2Model, ProteinInput, StructurePredictionInput
from esm.utils.msa import MSA

model   = EsmFold2Model.from_pretrained("biohub/ESMFold2", device="cuda").eval()
model.set_kernel_backend("fused")          # faster above ~300 residues
builder = ESMFold2InputBuilder()

spi = StructurePredictionInput(sequences=[
    ProteinInput(id="A", sequence=TARGET, msa=MSA.from_a3m("target.a3m")),  # target: MSA
    ProteinInput(id="B", sequence=BINDER, msa=None),                        # design: no MSA
])
result = builder.fold(model, spi, num_loops=20, num_sampling_steps=200,
                      num_diffusion_samples=1, seed=0)
pae, plddt = result.pae, result.plddt   # ← ipSAE inputs
```

- CLI equivalent: `scripts/predict.py --msa CHAIN=path.a3m` (repeatable).
- `msa_max_depth` subsamples that many random rows per loop, always keeping the query.
- **Critical gotcha**: *"the `-Fast` checkpoint has no alignment reader, and will ignore an MSA without complaining. For alignments, load the full model."* → the gate must use `biohub/ESMFold2`, not `-Fast`.
- Sources: [esmfold2_local_gpu.ipynb](https://github.com/Biohub/esm/blob/main/cookbook/tutorials/esmfold2_local_gpu.ipynb), [cookbook a3m examples `g3l5_chainA.a3m` / `g3l5_chainB.a3m`](https://github.com/Biohub/esm/tree/main/cookbook/tutorials)

MSA gain on FoldBench (DockQ ≥ 0.23 pass rate): protein–protein **70 % → 76 %**; antibody–antigen **50 % → 53 %** (55 % at 20 loops). AlphaFold3 reference: 73 % / 47 %.

## 1.3 Confidence outputs (CONFIRMED)

`result` exposes `plddt` (per-token), `ptm`, `iptm`, **`pae`** (full N×N matrix), optional distogram. The Batch API additionally has `FoldMaxAccuracyConfig(include_pae=True, include_pair_chains_iptm=True, include_embeddings=True)`. **PAE is available locally** — verified in a third-party driver that asserts `raw_pae.shape == (count, count)`.

⚠️ **Token-order gotcha** (will silently corrupt your ipSAE if ignored): the PAE matrix is in *token* order from `builder.prepare_input()`, which is **not** the mmCIF atom/residue order. You must build the permutation from `chain.tokens[*].token_index` and re-index, then assert `ca_chains == token_chain_ids`. Working reference implementation: [`predict_esmfold2.py`](https://github.com/RFjason/binder-design/blob/HEAD/campaigns/tnfr1-vhh-20261103/scripts/hpc_snapshot/predict_esmfold2.py).

## 1.4 VRAM / runtime (partly CONFIRMED, partly inferred)

| Fact | Value | Source |
|---|---|---|
| Official 24 GB guidance | "**24 GB VRAM … can likely handle complexes with ~500 total amino acids**"; H100-80 GB ≈ 2000 residues; "**no plans to make CPU offloading possible**" | [issue #320](https://github.com/Biohub/esm/issues/320) |
| Stored dtype | `"dtype": "float32"`, index `total_size` = **26,319,023,372 B (26.3 GB)** | [config.json](https://huggingface.co/biohub/ESMFold2/raw/main/config.json) |
| fp8 ESMC backbone | Added in **esm v3.4.2** (2026-10-30): "fp8 inference for the ESMC backbone in ESMFold2, plus faster decode and caching" | [releases](https://github.com/Biohub/esm/releases) |
| `low_memory_mode=True` | Proposed in **PR #385** (ESMC CPU offload after LM feature extraction + chunked confidence). **Still open as of 2026-11-11** | [#385](https://github.com/Biohub/esm/issues/385) |
| Runtime | H100, single-sequence, 10 loops / 200 diffusion steps, 1024 residues: **ESMFold2 15.8 s**, `-Fast` 9.4 s (1.3× faster than AF3) | preprint |

**Our two complexes:**

| Complex | Total aa | 24 GB L4 verdict |
|---|---|---|
| IL-7Rα (219) + binder (40–100) | **259–319** | **Fits** (well under the ~500 guidance) |
| TL1A **monomer** (180) + binder | **220–280** | **Fits** |
| TL1A **trimer** (540) + binder | **580–640** | **Exceeds the stated 24 GB ceiling** |

Necessary L4 measures, in order: (1) load with `dtype=torch.bfloat16` — a naïve fp32 load needs 26.3 GB and OOMs before any inference; (2) `pip install -U esm>=3.4.2` and use the fp8 ESMC path; (3) `model.set_kernel_backend("fused")`; (4) if the trimer still OOMs, score TL1A as a monomer locally and move trimer scoring to a rented L40S/A100 or to the sponsor platforms. **UNCONFIRMED**: exact measured peak VRAM in bf16/fp8 at 280 aa and at 640 aa — nobody has published numbers; we must measure it ourselves on first run (log `torch.cuda.max_memory_allocated()`).

Colab is explicitly a supported path — the official local-GPU notebook has a `# Colab: uncomment` cell. Install footprint is ~5.6 GB of wheels (4.6 GB of that is the CUDA torch build) plus 26.3 GB of weights → **cache the weights on Google Drive**, you will exceed Colab's ephemeral disk otherwise.

## 1.5 Inference API (CONFIRMED, cost UNCONFIRMED)

- Biohub Platform: `https://biohub.ai`, keys at `/developer-console/api-keys`.
  - Sync API: `SequenceStructureForgeInferenceClient(model="esmfold2-fast-2026-05")`, **768-token input limit**.
  - **Batch API**: `ForgeBatchClient`, **4096-token limit**, queued tasks expire 24 h, results kept 31 days.
  - `fold_max_accuracy` = **ESMFold2 Pro**, checkpoint `esmfold2-2026-05-cutoff-2025`; the service generates its own MSA (UniRef30 + ColabFold; `max_prefilter=100000`, `max_accept=10000`, `diff=3000`, `expand_diff=3000`), runs 5 seeds and picks the best. **Not the competition config** (5 seeds, service-side MSA for every chain) — useful for *design validation*, useless for *gate reproduction*.
- Models are "freely available"; **UNCONFIRMED** whether the free tier has enough quota for thousands of folds. Guardrails block controlled pathogen/toxin sequences — irrelevant for our targets.
- Sources: [README §Batch API](https://github.com/Biohub/esm), [fold_max_accuracy.md](https://github.com/Biohub/esm/blob/main/cookbook/snippets/fold_max_accuracy.md)

## 1.6 Which checkpoint does the gate use? (UNCONFIRMED — ranked)

Published HF checkpoints: `ESMFold2`, `ESMFold2-Fast`, `ESMFold2-Experimental`, `ESMFold2-Experimental-Cutoff2025`, `ESMFold2-Experimental-Fast`, `ESMFold2-Experimental-Fast-Cutoff2025`, plus `*-base{300M,600M,6B}-step*` ablations.

| Rank | Candidate | Evidence |
|---|---|---|
| **1 (likely)** | `biohub/ESMFold2` (cutoff 2021) | Only non-Experimental checkpoint that reads MSAs; the canonical/most-downloaded release; `config.json` default `num_loops: 20` matches what independent reproductions assume |
| 2 | ESMFold2 Pro `esmfold2-2026-05-cutoff-2025` (platform-only) | Used by `fold_max_accuracy`; but that service forces 5 seeds + its own MSA for all chains, contradicting "1 seed … MSA only for the target" |
| 3 | `ESMFold2-Experimental*` | Has `MSAEncoder`; used by the official design code; unlikely for a neutral public gate |
| 4 | `ESMFold2-Fast` | **Ruled out** — silently ignores MSAs, so it cannot implement "MSA for the target sequence" |

**Action**: ask the organiser in the group chat for (a) checkpoint name, (b) `num_loops` / `num_sampling_steps`, (c) the a3m they use for each target, (d) how many target chains are folded, (e) which ipSAE column/direction. Until answered, reproduce with the defaults below.

## 1.7 MatwingsVenus / 晓鹜 announcement — BLOCKED

`https://mp.weixin.qq.com/s/5ef89061e1247218f0d0432148385a89` returns **参数错误** ("parameter error") to WebFetch — the WeChat article is not retrievable without a WeChat session. I could not find the text mirrored anywhere else. **UNCONFIRMED**: whether MatwingsVenus names the exact checkpoint/version. Their platform is `https://matvenus.com/zh` and per our own fact base it hosts ESMFold2 + ipSAE + RFdiffusion + RFantibody + BindCraft + BoltzGen + Boltz2 — **this is the cheapest way to see the organiser-blessed defaults**, because a sponsor platform advertising "ESMFold2 + ipSAE for the competition" will have hard-coded the gate config in its job form. **Recommend: register, open the ESMFold2 job form, and screenshot every default parameter.**

---

# Q2 — ipSAE: exact definition, implementation, calibration

## 2.1 Papers and repos (CONFIRMED)

| | |
|---|---|
| v1 paper | Dunbrack, *Rēs ipSAE loquunt: What's wrong with AlphaFold's ipTM score and how to fix it*, [bioRxiv 2025.02.10.637595](https://www.biorxiv.org/content/10.1101/2025.02.10.637595v1), [PMC11844409](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11844409/) |
| v1 repo | **https://github.com/DunbrackLab/IPSAE** (single file `ipsae.py`, numpy only) |
| v2 paper | *Rēs ipSAE adhūc loquuntur: What the directionality of predicted aligned error … can tell us*, [bioRxiv 2026.11.02.763449](https://www.biorxiv.org/content/10.64898/2026.11.02.763449v1) (posted 2026-11-02) |
| v2 repo | **https://github.com/DunbrackLab/IPSAE2** — `pip install .` → `ipsae` CLI, **alpha v0.1.10**, adds nucleic acid / ligand / ion support and `--group` |

**Which one does the gate use? UNCONFIRMED.** v2 is 9 days old and self-described as "not yet validated across every prediction format" — the gate almost certainly uses **v1 `ipsae.py`**. The independent reproduction I found pins v1 by SHA-256. Use v1 as primary; compute v2 as a sanity cross-check only.

## 2.2 Exact formula (read from `ipsae.py` source — CONFIRMED)

```python
ptm_func(x, d0) = 1.0 / (1 + (x/d0)**2)                       # line 111-112
calc_d0(L)      = max(1.0, 1.24*(L-15)**(1/3) - 1.8)          # line 116-124
                                                               # d0 = 2.0 if either chain is nucleic acid
```

For an ordered chain pair (chn1 = **aligned**, chn2 = **scored**), with a PAE cutoff `C`:

1. `valid[i][j] = (chain(i)==chn1) & (chain(j)==chn2) & (PAE[i][j] < C)`
2. For each residue `i` of chn1: `n0res_i = Σ_j valid[i][j]`, `d0res_i = calc_d0(n0res_i)`
3. `ipsae_byres[i] = mean_{j: valid[i][j]} ptm_func(PAE[i][j], d0res_i)` (0 if no valid j)
4. **`ipSAE(chn1→chn2) = max_i ipsae_byres[i]`**  ← a **max over aligned residues**, not a mean
5. Reported `asym` rows give both directions; the `max` row gives `max(A→B, B→A)`

Also emitted per pair: `ipSAE_d0chn` (d0 from `len(chn1)+len(chn2)`), `ipSAE_d0dom` (d0 from the union of residues with any interchain PAE<C), `ipTM_af` (predictor-native), `ipTM_d0chn`, `pDockQ`, `pDockQ2`, `LIS`, and the counts `n0res / n0chn / n0dom`, `d0res / d0chn / d0dom`, `nres1 / nres2`, `dist1 / dist2`. The variant the competition means by bare "ipSAE" is the **`d0res`** column (the 6th field, header `ipSAE`). **UNCONFIRMED** but strongly implied by the README's own column ordering.

## 2.3 Inputs and exact CLI (CONFIRMED)

```bash
pip install numpy
# AlphaFold2 / ColabFold
python ipsae.py scores_rank_001_..._seed_000.json  unrelaxed_rank_001_....pdb  15 15
# AlphaFold3 style (what you build from ESMFold2)
python ipsae.py fold_xxx_full_data_0.json  fold_xxx_model_0.cif  10 10
# Boltz
python ipsae.py pae_XXX_model_0.npz  XXX_model_0.cif  10 10
```
Positional args: `<pae/scores file> <pdb|cif> <pae_cutoff> <dist_cutoff>`. Outputs `<model>_<PAE>_<DIST>.txt` (chain-pair table), a by-residue file, and a PyMOL script.

**ESMFold2 → ipSAE adapter** (the one piece you must write; working third-party version linked in §1.3). Write an AF3-style JSON:
```json
{"pae": [[...]], "atom_plddts": [...], "token_chain_ids": [...], "token_res_ids": [...]}
```
then `python ipsae.py af3like_pae.json model.cif 10 10` and parse rows where field 4 == `asym`, taking field 5 as the directional ipSAE. Note `plddt` may come back on a 0–1 scale; rescale ×100 before writing `atom_plddts` (only affects pDockQ, not ipSAE).

IPSAE2 CLI (cross-check): `ipsae model.cif pae.json [summary.json] [--group B+C] [--output-file x.tsv]`, defaults PAE ≤10 Å, contact ≤5 Å; `ipsae12` aligns on chn1, `ipsae21` on chn2, reported `ipsae` = larger direction for same-type pairs.

## 2.4 Cutoffs (partly UNCONFIRMED)

- v1 paper: "cutoffs of **10 or 15 Å** may be most suitable".
- v2 default: **PAE ≤ 10 Å**, heavy-atom contact ≤ 5 Å.
- Independent gate reproduction uses **`10 10`**, variant `d0res`, score = `min(ipSAE_A→B, ipSAE_B→A)`.
- **Assume `10 10` + `d0res`. UNCONFIRMED whether the organiser takes the `max` row, the two `asym` rows' min, or the `A→B` direction only.** Hedge by requiring **both directions ≥ target** (optimising the min is strictly harder, so it also wins if they use max).
- **Which direction is "principled"?** The ipSAE2 paper ([bioRxiv 2026.11.02.763449](https://www.biorxiv.org/content/10.64898/2026.11.02.763449v1)) shows the PAE matrix is intrinsically asymmetric — "the predicted position of tokens in one molecule is evaluated relative to coordinate frames defined by tokens in another" — and recommends **aligning on the more complex molecule and scoring the simpler one**. For our case (target 180–219 aa vs binder 40–100 aa) that is **align on the target, score the binder** = `ipSAE(target→binder)`. v1 for same-type protein pairs instead reports the **larger** of the two directions. Note the new numeric thresholds in that paper (0.5 for AF3 at 93.8 % balanced accuracy, 0.8 for Boltz-2 at 86.2 %) are **protein–ligand** cutoffs validated against a 5 Å ligand RMSD — do **not** transfer them to protein–protein minibinder scoring.

## 2.5 How it differs from ipTM / pDockQ

| Metric | Basis | Failure mode ipSAE fixes |
|---|---|---|
| **ipTM** (AF2/AF3/ESMFold2 native) | TM-like sum over **all** interchain residue pairs, `d0` from **total** chain lengths | Disordered or non-interface regions dilute the score; a tiny real interface on two big chains scores low |
| **pDockQ** | logistic fit on interface pLDDT × log(#contacts) | pLDDT-driven; insensitive to interchain *positioning* |
| **pDockQ2** | per-pair PAE-based | no PAE-cutoff restriction, no adaptive `d0` |
| **LIS** | transform of interchain PAEs | no `d0` normalisation |
| **actifpTM** | also restricts to low-PAE pairs | does **not** rescale `d0` → weaker suppression of false interfaces |
| **ipSAE** | restricts to interchain pairs with PAE < C **and** sets `d0` from the count of those pairs | isolates the interface signal; adaptive `d0` |

Published separation example: RAF1–KSR1 (true) **ipSAE 0.73** vs ipTM 0.41; RAF1–RIPK1 (false) **ipSAE 0.04** vs ipTM 0.28. In the README's worked 3-chain example, real pairs give ipSAE 0.65–0.67 and the non-contacting B–C pair gives 0.32.

## 2.6 Calibration — what score do we actually need

**(a) Literature threshold.** Meta-analysis over **3,766 experimentally characterised de novo binders** ranks **AF3 `ipSAE_min` as the single best predictor of experimental success**, +1.4× average precision over AF2 `pae_interaction`, with the recommended filter **`ipSAE_min` > 0.61**. Companion filters: `ipSAE_min × interfaceΔG/ΔSASA < −1.5`; `LIS × input interface shape-complementarity > 0.42`; pre-filters `shape complementarity > 0.62`, `RMSD_binder(design vs AF3) < 3.73 Å`. Caveat the authors flag: precision at max-F1 for `ipSAE_min` ranges **0.1 → 1.0 across targets**, so the threshold is not target-transferable.

**(b) Real ESMFold2+ipSAE distribution.** A competitor published 1,127 scored VHH designs against TNFR1, produced with exactly the gate recipe. Distribution of `ipsae_min` (PAE/dist 10/10, `d0res`):

| p0 | p10 | p25 | **p50** | p75 | p90 | p95 | **p99** | max |
|---|---|---|---|---|---|---|---|---|
| 0.000 | 0.000 | 0.000 | **0.000** | 0.151 | 0.364 | 0.496 | **0.639** | **0.698** |

mean 0.097; **70 % of designs score < 0.05**; only 1.7 % reach ≥ 0.60. By generator: BindCraft2 n=250 mean 0.203 max 0.673 (6.0 % ≥0.6); RFD4 n=520 mean 0.070 max 0.698 (0.8 %); RFantibody n=357 mean 0.061 max 0.491 (**0 %**). Correlations in that dataset: `ipsae_min` vs native `ipTM` **r = 0.85**; vs `pose_dockq` **r = 0.56**; vs binder length **r = 0.17**; vs monomer pLDDT **r = −0.16**.

**(c) AsterFire's TL1A numbers** (0.637, 0.628, binder pLDDT 91.8, "credible interfaces") sit at **p99 of that VHH distribution** — i.e. genuinely good, and the public bar.

**Conclusion.** ESMFold2 minibinder hit rates are ~3–4× the scFv/VHH rates (53.8–70.0 % vs 12.1–21.0 %), so the minibinder ipSAE distribution will be shifted up. **Target ≥ 0.70, treat ≥ 0.61 as the floor, discard anything < 0.5.** With 1000 slots per track and (per the fact base) a large public field, assume the cut lands somewhere in **0.55–0.70** — **UNCONFIRMED**, nobody knows the field size.

---

# Q3 — TL1A / TNFSF15 target biology

## 3.1 Numbering — get this right or everything downstream is wrong (CONFIRMED)

UniProt **O95150**, 251 aa. Secreted/THD form = **72–251**; THD domain 95–251. Competition construct = **Leu72–Leu251 = 180 aa**, so:

```
construct index i  ⇄  UniProt (i + 71)      [construct 1 = Leu72]
2RE9 auth number   ==  construct index      (verified: 2RE9:A spans auth 17–180 → UniProt 88–251)
3K51 auth number   ==  construct index + 4  (verified: 3K51:A spans auth 25–184 → UniProt 92–251)
patent numbering   ==  UniProt + 2          (BMS SEQ ID 19; motif-matched, see §3.4)
```
**2RE9's author numbering is already construct-relative — hand RFdiffusion 2RE9 hotspots directly, no offset.** 3K51 needs −4.

Post-translational caveats for wet-lab: **N-glycosylation at UniProt N133 and N229** → construct **N62, N158**; **disulfide UniProt C162–C202** → construct **C91–C131**. Avoid designing over the glycan sites.

## 3.2 2RE9 — what it actually contains (CONFIRMED)

**"Crystal structure of TL1A at 2.1 Å"**, X-ray, **one polymer entity, 3 deposited chains A/B/C** = the physiological TL1A **homotrimer**, **apo — no receptor**. ([RCSB 2RE9](https://www.rcsb.org/structure/2RE9))

## 3.3 Is there a TL1A/DR3 complex? (answer: no — use DcR3)

- **No TL1A:DR3 (TNFRSF25) experimental structure exists** that I could find. Literature on DR3 selectivity works from *models* of TL1A:DR3 built on the DcR3 complex. **UNCONFIRMED** that none exists, but an RCSB full-text search and the 2022 engineered-TL1A paper both rely on models.
- **Mutagenesis pins down which parts of the footprint are DcR3-specific.** The engineered "Fc-scTL1A-**AKALFF**" variant carries **K111A, L123K, Q167A, S187L, E190F, N207F** (UniProt numbering) and **retains wild-type-like DR3 affinity while losing all measurable DcR3 binding**. ([Sci Rep 12:20538, 2022](https://www.nature.com/articles/s41598-022-24984-y)) Four of those six (**L123, S187, E190, N207**) are inside my computed 3K51 footprint; the other two (**K111, Q167**) fall inside the BMS antibody epitope regions (§3.5). In construct numbering: **K40, L52, Q96, S116, E119, N136**.
  → Interpretation: these six are **DcR3-weighted**, so the shared DR3 epitope is roughly the 3K51 footprint *minus* these positions. For our purposes this matters little (ranking is by KD to TL1A, not by DR3 selectivity), but it is useful if we want the +5k RMB cell-activity bonus — a binder that blocks DR3 signalling should engage the DcR3/DR3-**shared** residues (R32, E49, E51, G53, R85, G86, M87, D104, S105, D115, Y117, P118, P120, T121, W137, F138, Q150, Y167, T168, K169, E170 in construct numbering) rather than the six DcR3-specific ones.
- **TL1A:DcR3 = PDB 3K51**, 2.45 Å, chain A = TL1A secreted form, chain B = decoy receptor 3. ([RCSB 3K51](https://www.rcsb.org/structure/3K51)) DcR3 and DR3 compete for the same ligand surface (DcR3 is the decoy), so **3K51 is the best available proxy for the receptor-binding epitope**.

## 3.4 The epitope — computed, not guessed

I downloaded the 3K51 biological assembly and computed all TL1A residues within **4.5 Å (heavy atom)** of DcR3. 25 residues, and **the set reproduces the published critical-residue list exactly** (Arg36, Leu56, Gly57, Met91, Asp108, Ser120, Tyr121, 122, Glu123, 124, Asn140, Phe142, Thr172, Glu174 in 3K51 numbering) → high confidence.

| 3K51 auth | **construct idx (= 2RE9 auth)** | UniProt |
|---|---|---|
| R36 | **R32** | R103 |
| E53 E55 L56 G57 | **E49 E51 L52 G53** | E120 E122 L123 G124 |
| R89 G90 M91 | **R85 G86 M87** | R156 G157 M158 |
| D108 S109 | **D104 S105** | D175 S176 |
| D119 S120 Y121 P122 E123 P124 T125 | **D115 S116 Y117 P118 E119 P120 T121** | D186 S187 Y188 P189 E190 P191 T192 |
| N140 W141 F142 | **N136 W137 F138** | N207 W208 F209 |
| Q154 | **Q150** | Q221 |
| Y171 T172 K173 E174 | **Y167 T168 K169 E170** | Y238 T239 K240 E241 |

Independent corroboration: **UniProt O95150 annotates residues 187 and 190 as "Important for binding TNFRSF6B [DcR3]"** → construct **S116** and **E119**, both inside my computed patch. ([UniProt JSON](https://rest.uniprot.org/uniprotkb/O95150.json))

**The epitope straddles the inter-protomer groove.** I checked each epitope residue's distance to the neighbouring protomers in the 2RE9 trimer:
- Residues at the inter-protomer interface (≤4.5 Å to B or C): **S105, T121, N136, W137, F138**
- Very close to a neighbour (≤5 Å): **D104 (4.7 Å to C), T121 (4.1 Å to C), N136/W137/F138 (2.9–3.4 Å to B), E170 (4.6 Å to B), K169 (5.9 Å to B)**
- Clearly on one protomer's own face (≥9 Å from both neighbours): **R32, E49, E51, L52, Y117, E119, P120**

So the receptor site is the canonical TNF-superfamily **cleft between two adjacent protomers**, as expected. This has a direct consequence for the gate — see Q6 §6.4.

**Recommended RFdiffusion/BoltzGen hotspots (construct = 2RE9 numbering):**
- **Primary (upper patch, groove, therapeutically validated region):** `A104,A105,A115,A116,A117,A119,A121` — i.e. `D104 S105 D115 S116 Y117 E119 T121`
- **Secondary (lower patch, near-neighbour):** `A136,A137,A138,A169,A170`
- **Tertiary (single-protomer face, safest for a monomer-folded gate):** `A32,A49,A51,A52,A85,A86,A87,A117,A119,A120`
- **Avoid:** N62, N158 (glycans); C91/C131 (disulfide); the extreme C-terminus (170–180) which is partly buried in the trimer core.

A hydrophobic anchor is available: **W137 / F138 / M87 / L52** — the TNFα minibinder work succeeded by targeting exactly such "conserved hydrophobic surface" on a TNF-superfamily trimer. ([Nat Commun Biol 2025](https://www.nature.com/articles/s42003-025-09030-7))

## 3.5 Known anti-TL1A biologics (for epitope validation + novelty)

| Agent | Modality | Epitope / notes |
|---|---|---|
| **Tulisokibart** (PRA023 / MK-7240, Merck/Prometheus) | fully human anti-hTL1A IgG1κ | Neutralises TL1A; in Ph2/3 for UC & Crohn's. **Epitope not public — UNCONFIRMED.** ([DrugBank DB18719](https://go.drugbank.com/drugs/DB18719), [NCT04996797](https://cdn.clinicaltrials.gov/large-docs/97/NCT04996797/Prot_000.pdf)) |
| **Duvakitug** (TEV-48574, Teva/Sanofi) | anti-TL1A IgG | Ph2b in IBD; **epitope UNCONFIRMED** ([NCT04545385](https://cdn.clinicaltrials.gov/large-docs/85/NCT04545385/Prot_000.pdf)) |
| **10A4.F7 / 10A4.F7.2E8** (Bristol Myers Squibb) | mouse-derived mAb, patented | **Discontinuous epitope**, HDX-MS: peptides **102–116 `TVVRQTPTQHFKNQF`** and **166–180 `EIRQAGRPNKPDSIT`**; refined to **169–172 `QAGR`** and **113–116 `KNQF`** (patent numbering). Motif-matching to O95150 gives **UniProt 100–114 and 164–178**, i.e. **construct 29–43 and 93–107**; `QAGR` = construct **96–99**, `KNQF` = construct **40–43**. ([US10968279B2](https://patents.google.com/patent/US10968279B2)) |
| Br3s95, Tr3s90 + 4 others | **nanobodies**, GNN-designed | bind TL1A, KD ≈ **2.0 µM** — weak; useful as the affinity bar to beat ([bioRxiv 2025.06.08.658414](https://www.biorxiv.org/content/10.1101/2025.06.08.658414v1.full.pdf)) |

Note the BMS epitope (construct 29–43 + 93–107) **overlaps the receptor site only partially** (it shares R32 and the 93–107 stretch containing D104/S105) — it is a distinct, also-neutralising surface. That gives us a legitimate **second epitope option** if the groove proves hard to hit.

**Novelty check for TL1A: low risk.** I found **no published de novo minibinder for TL1A** (searched bioRxiv/Nature/PDB). The known binders are antibodies and nanobodies, which are a different modality and far above 100 aa — an alignment covering ≥80 % of our 40–100 aa design against a 130 aa VHH is implausible. **UNCONFIRMED** that none exists in a non-indexed source; still run a BLAST of every submission against the PDB + a nanobody database before submitting.

---

# Q4 — IL-7Rα / CD127 target biology

## 4.1 Numbering (CONFIRMED)

UniProt **P16871-1**, 459 aa. Signal peptide **1–20**; extracellular **21–239**; TM 240–264. Competition construct = **Glu21–Asp239 = 219 aa** = the **complete ectodomain**, and construct index 1 = UniProt 21, so:

```
construct index i  ⇄  UniProt (i + 20)
3DI3 chain-B auth  ==  construct index      (verified: 3DI3:B auth 17–209 → UniProt 37–229)
```
**3DI3's author numbering is already construct-relative** — and it is also the numbering used in the anti-CD127 antibody literature. One frame for everything.

**Domain split (CONFIRMED, UniProt):** Fibronectin type-III = **UniProt 131–231** → **construct 111–211 = D2**. Therefore **D1 ≈ construct 1–110**.
**Glycans:** UniProt N49, N65, N151, N182, N232, N233 → construct **N29, N45, N131, N162, N212, N213**.
**Disulfides:** UniProt 42–57, 74–82, 108–118 → construct **22–37, 54–62, 88–98**.

## 4.2 3DI3 — verified

**"Crystal structure of the complex of human interleukin-7 with glycosylated human interleukin-7 receptor alpha ectodomain"**, 2.9 Å. **Chain A = IL-7, chain B = IL-7Rα.** Confirmed — the competition's "3DI3 chain B" is the receptor. ([RCSB 3DI3](https://www.rcsb.org/structure/3DI3)) Related: [PNAS, structural reorganisation of the IL-7 signalling complex](https://www.pnas.org/doi/10.1073/pnas.1116582109); apo IL-7Rα ectodomain forms an "X" homodimer with C-termini 110 Å apart.

## 4.3 The IL-7 epitope on IL-7Rα — computed

4.5 Å heavy-atom interface, 3DI3 chain B vs chain A. **Construct index == 3DI3 auth:**

| Region | Residues (construct = 3DI3 auth) | UniProt | Domain |
|---|---|---|---|
| D1 loop 1 | **S31, H33** | S51, H53 | D1 |
| D1 loop 2 | **L57, V58, E59** | L77, V78, E79 | D1 |
| **D1 core strand (dominant)** | **K77, K78, F79, L80, L81, I82, G83** | K97–G103 | D1 |
| D2 fringe | **K138, Y139** | K158, Y159 | D2 |
| D2 fringe | **Y192, F193** | Y212, F213 | D2 |

IL-7 side of the same interface (3DI3 chain A auth): K10 Q11 S14 V15 L16 V18 S19 Q22 T72 D74 L77 H78 L80 K81 E84 G85 I88.

**Answer to the flagged question: the IL-7 site is carried overwhelmingly by D1** (12 of 16 contacts, including the whole hydrophobic `K77-KFLLIG-83` strand), with four peripheral D2 contacts. A 40–100 aa minibinder can cover the entire D1 patch.

**Independent corroboration from the clinical antibody literature:** anti-CD127 antibodies are claimed against epitopes `24SQLEVNGSQHSLTCAF39` and `73FIETKKFLLIGKSNIC88`, and a specific epitope **"H33, E75, F79, I82, K84, M144, R186, H191, Y192"** — note **H33, F79, I82, Y192 are in my computed IL-7 interface**, and `FLLIG` in the second peptide is exactly construct 79–83. ([US11440964](https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/11440964), [US10428152](https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/10428152), [WO2020154293A1](https://patents.google.com/patent/WO2020154293A1/en)) **Lusvertikimab (OSE-127)** is described as binding **"site-1 and site-2b of CD127"**, blocking heterodimerisation with γc without affecting TSLPR signalling — site-1/2b maps onto this same D1+D2-fringe surface. ([OSE Immunotherapeutics](https://www.ose-immuno.com/en/our-products/ose-127-modular/), [Blood 2024](https://ashpublications.org/blood/article/143/26/2735/515436/The-IL-7R-antagonist-lusvertikimab-reduces)) Another structure exists: **6P67 = IL-7Rα + anti-IL-7Rα 2B8 Fab**. ([PDBj 6P67](https://pdbj.org/mine/summary/6p67))

**Recommended hotspots (construct = 3DI3 auth):**
- **Primary (IL-7 site, therapeutically validated):** `B77,B79,B80,B81,B82,B83` (`K77 F79 L80 L81 I82 G83`) — hydrophobic, compact, exactly the kind of patch minibinders are good at
- **Extend for a larger interface:** add `B31,B33,B57,B58,B59`
- **Do not include** the D2 fringe (138/139/192/193) in the primary hotspot set — they are 40+ Å away along the receptor and will produce an implausibly elongated binder
- **Avoid:** N29, N45 (D1 glycans, close to the S31/H33 loop!); N131, N162; the 22–37 / 54–62 / 88–98 disulfides

⚠️ **Glycan hazard specific to IL-7Rα**: construct **N29** sits two residues from hotspot **S31**, and **N45** is in D1. The gate folds protein-only (no glycans) so ipSAE will not see them, but the wet-lab construct (expressed with C-terminal His, presumably in mammalian cells) **will** be glycosylated. A design that packs against N29/N45 will score well in silico and fail by SPR. Explicitly penalise binder contacts within ~8 Å of N29/N45 Cβ.

## 4.4 Novelty — **this is a real constraint for IL-7Rα**

De novo IL-7Rα minibinders are **published and in the PDB**. Cao et al., *Design of protein-binding proteins from the target structure alone*, Nature 605:551–560 (2022) designed IL-7Rα minibinders that "use two helices to bind IL-7Rα" with a crystal structure at 1.63 Å RMSD to the design model. ([Nature](https://www.nature.com/articles/s41586-022-04654-9), [Baker Lab PDF](https://www.bakerlab.org/wp-content/uploads/2022/03/Cao_etal_Nature2022_Design_of_binders_from_target_structure_alone.pdf))

**PDB 7S5B — "Unbound State of a De novo designed Protein Binder to the Human Interleukin-7 Receptor", 55 aa:**
```
SVIEKLRKLEKQARKQGDEVLVMLARMVLEYLEKGWVSEEDADESADRIEEVLKK
```
([RCSB 7S5B](https://www.rcsb.org/structure/7S5B), [FASTA](https://www.rcsb.org/fasta/entry/7S5B))

Also: a 2025 complementarity/dock-and-design paper designed binders to **two IL-7Rα epitopes** — one overlapping the native IL-7 site, one at a previously untargeted site — getting nanomolar binders from **8 sequences per site**. ([PMC12412469](https://pmc.ncbi.nlm.nih.gov/articles/PMC12412469/)) IL-7Rα also appears in RFdiffusion VHH work and in AlphaProteo / Chai-2 benchmarks.

**Action (mandatory before submission):** pairwise-align all 20 IL-7Rα candidates against (i) 7S5B's 55-mer above, (ii) every IL-7Rα minibinder sequence in the Cao 2022 supplementary tables, (iii) the 2025 complementarity-paper designs. Rule: if the aligned region covers ≥80 % of our design's full length, identity must be **≤70 %**. The danger is not literal copying — it is that de novo helical minibinders against the same hydrophobic patch converge on similar 3-helix sequences with high `E/K/L/A` content. 7S5B is a 3-helix `EKLRKLEKQARKQ…` pattern; if our generator emits an idealised ferredoxin/3-helix bundle with that composition, **spot identity can plausibly cross 70 % over a 44-residue window**. Enforce a hard filter, and bias toward 4-helix, helix-loop-helix-strand, or 70–100 aa topologies to break the alignment coverage condition.

---

# Q5 — Best pipeline that actually runs in 19 days

## 5.1 The candidates, with numbers

### (A) Official ESMFold2 inversion (Biohub `binder_design.py`) — **the gate's own generator, but H100-only**
- Hallucination by gradient descent on soft sequence logits through `EsmFold2ExperimentalModel`, with an **ESMC language-model loss** (masked-LM regularisation, 15 % mask, 4 passes), losses `{intra_contact 0.5, inter_contact 0.5, glob 0.2, epitope 0.5}`, `STEPS=150`, `lr=0.1`, cysteines hard-masked out (`logits[:,:,CYS_IDX] = -1e6`).
- **Explicit epitope conditioning**: `target_hotspot_ids=["Y117","E119",...]` (1-indexed within the provided target sequence, residue letter must match), `epitope_contact_distance=12.0 Å`, loss = mean over hotspots of the best-5 binder contacts.
- Built-in minibinder prompt: `PromptFactory(name="minibinder", template="{seq}", length_ranges={"seq": (60, 200)})` — note **their minibinder prior is 60–200 aa**; our cap is 100.
- Ranks samples by **ipTM** and a cheap `distogram_iptm_proxy` ("Algorithm 15").
- **VRAM, from the source comment:** target `cd45` (185 aa) + a VH-VL binder at `batch_size=1` → **51 GB, reduced to 27 GB with `REUSE_ESMC=True`**. Deployment is `@app.cls(gpu="H100", timeout=3600, cpu=16, memory=80GB)`.
- Published search cost: minibinders **~500 H100-h (low-compute) to ~2,400 H100-h (high-compute)** for 15 k–67 k candidates.
- Published experimental hit rates: **minibinders 53.8 % → 70.0 %**, scFvs 12.1 % → 21.0 %; validated affinities **68 pM – 70 nM**; cryo-EM-confirmed EGFR minibinder at 1.204 Å RMSD; PD-L1 minibinder IC50 1.6 nM vs atezolizumab 2.6 nM.
- **Verdict: perfectly aligned with the gate (it literally optimises the scorer), completely out of reach on an L4.** 27 GB > 24 GB even with the memory flag, and that was for a smaller complex than ours.
- Sources: [binder_design.py](https://github.com/Biohub/esm/blob/main/cookbook/tutorials/binder_design.py), [binder_design.ipynb](https://github.com/Biohub/esm/blob/main/cookbook/tutorials/binder_design.ipynb), [preprint](https://www.biorxiv.org/content/10.64898/2026.06.03.729735v1.full)

### (B) RFdiffusion family + ProteinMPNN + AF2 initial-guess — the classic
- **RFdiffusion3** open-sourced **2025-12-03**, weights + training code at [`RosettaCommons/foundry`](https://github.com/RosettaCommons/foundry) (production branch), BSD-3 covering code *and* weights, 10× faster inference than RFdiffusion2. ([IPD announcement](https://www.ipd.uw.edu/2025/12/rfdiffusion3-now-available/))
- **RFD4-Proteina** (IPD × NVIDIA) released **2026-10-23**, open and public. ([Rosetta Commons](https://rosettacommons.org/2026/10/23/rfd4-proteina-improves-computational-protein-design-through-open-science-and-collaboration/))
- Standard filters (Bennett/Cao convention): **`pae_interaction` < 10** (the single best classic predictor), **`plddt_binder` > 80**, **`rmsd_binder` (AF2 vs the backbone MPNN was given) < 2 Å**.
- Reported pass rates at `pAE≤10 & pLDDT≥80 & RMSD≤2`: **RFD4** 3/24 and 2/24 with its own sequences, **9/96** with ProteinMPNN sequences; **RFD3+ProteinMPNN 0/96**. Realistic in-silico yield ≈ **2–10 %**, and the *experimental* hit rate for classic RFdiffusion minibinders is low-single-digit % of designs ordered.
- **On the competitor's own ESMFold2+ipSAE numbers, RFD4 produced only 0.8 % of designs ≥ 0.6** (4/520 passing the full screen) and RFantibody produced **zero**. That is the strongest practical argument against making this the primary route.
- Colab: `ColabDesign` / `RFdiffusion.ipynb` (`sokrypton/ColabDesign`) runs on a T4/L4 for small targets, but a 180–219 aa target plus AF2 initial-guess filtering is slow and the maintained path is local.

### (C) BindCraft — AF2 backprop hallucination
- [Nature 2025](https://www.nature.com/articles/s41586-025-09429-6), [repo](https://github.com/martinpacesa/BindCraft) (MIT, actively maintained), [Colab notebook](https://github.com/martinpacesa/BindCraft/blob/main/notebooks/BindCraft.ipynb), [wiki](https://github.com/martinpacesa/BindCraft/wiki/De-novo-binder-design-with-BindCraft).
- **VRAM scales semi-linearly with (target + binder) length**: ~**16 GB** footprint typical, **~30 GB for a 500-residue complex**. Ideal GPUs L40/V100/A100/H100.
- **~138 s GPU per full hallucination** on benchmark targets; **50–80 backbones per H100-hour**; **24–72 GPU-hours for the recommended 100 accepted designs**.
- **Experimental success: mean 46.3 % across 12 targets (range 10–100 %)** — far better than classic RFdiffusion.
- **On ESMFold2+ipSAE, BindCraft2 was the best generator in the competitor's data**: mean 0.203, p90 0.564, max 0.673, **6.0 % ≥ 0.6**, and 15 of 19 screen-passers came from it, at **12.29 GPU-h for 549 scored designs (2.46 designs/parent)**.
- **L4 verdict: our complexes are 259–319 aa (IL-7Rα) and 220–280 aa (TL1A monomer) → comfortably inside 16–20 GB. BindCraft fits on an L4.** TL1A trimer target (540 aa) will not.

### (D) BoltzDesign1 — inverts Boltz-1 (AF3 reproduction)
- [Repo](https://github.com/yehlincho/BoltzDesign1) (public), **[official Colab](https://colab.research.google.com/github/yehlincho/BoltzDesign1/blob/main/Boltzdesign1.ipynb)**. Uses only the Pairformer + Confidence modules → much cheaper than full backprop. Strong in-silico success + diversity, generalises to non-protein targets.
- Weaker *experimental* validation record than BoltzGen. Good cheap diversity generator.

### (E) BoltzGen (2025-11) — **best published hit rate per design tested**
- [Paper](https://www.biorxiv.org/content/10.1101/2025.11.20.689494v2.full), [repo `HannesStark/boltzgen`](https://github.com/HannesStark/boltzgen) — **training code, inference code, weights and all designs under MIT**.
- **66 % success for both nanobody and miniprotein modalities against novel targets, nanomolar affinities, testing ≤15 designs per target**; **80 % on well-characterised targets (PD-L1, TNFα) with ≤20 designs per target**; SPR/BLI-confirmed nM binders for 4/5 targets.
- Ships a filtered, ranked, diversity-optimised output set; ranking weights were fit to maximise enrichment on **11,000 experimentally validated binder designs**.
- **VRAM: minimum 24 GB, 48 GB (L40S) recommended; "a 24 GB card will OOM on anything but small systems"; A10G-24 GB officially undersized.** ~**30–60 s per design on an L40S**. Model downloads ~6 GB (Boltz) + **~33 GB (it bundles ESMFold2 weights for its own filtering)**.
- **≤20 designs per target is exactly our submission budget** — this is the decisive argument.

### (F) 2026 newcomers worth knowing (did not fully verify)
- **BoltzProt-1** — "efficient de novo binder design with good developability", [tech report](https://boltz.bio/boltzprot1-technical-report.pdf), [bioRxiv 2026.06.23.733997](https://www.biorxiv.org/content/10.64898/2026.06.23.733997v1).
- **TorchCraft** — "unified binder design by inverting an all-atom structure predictor", [arXiv 2609.19770](https://arxiv.org/pdf/2609.19770).
- **EasyNano** — epitope-targeted design by **differentiable distogram optimisation with ESMFold2**, [arXiv 2606.12772](https://arxiv.org/pdf/2606.12772). Nanobody-CDR focused, so not directly a minibinder tool, but it is the cheapest known way to steer **ESMFold2 specifically** (distogram optimisation instead of full backprop) and `distogram_conditioning` is exposed in the official API. **UNCONFIRMED**: its VRAM/runtime and whether code is public — the PDF text layer did not yield these.
- **ProteinDJ** — modular high-performance RFdiffusion→MPNN→AF2 pipeline, [bioRxiv 2025.09.24.678028](https://www.biorxiv.org/content/10.1101/2025.09.24.678028v2.full.pdf). Useful if you want the classic route industrialised.

## 5.2 Recommendation

**Decouple the generator from the selector. The selector is non-negotiable: ESMFold2 + ipSAE, reproduced exactly.** The generator should be whatever produces the most diverse pool of plausible epitope-targeted minibinders per GPU-hour available to us.

**Plan (19 days, 2 targets, 20 sequences each):**

| Days | Step | Where |
|---|---|---|
| 1–2 | **Build the gate reproducer first.** Prepare target a3m (ColabFold MSA server / UniRef30) for TL1A-72-251 and IL-7Rα-21-239. Write the ESMFold2→AF3-JSON→`ipsae.py 10 10` adapter with the token-permutation assertion. Validate it on a **positive control** (3DI3 IL-7 + IL-7Rα should score high) and a **negative control** (IL-7Rα + an unrelated 60-mer should score ~0). Measure peak VRAM at 280 aa and 640 aa. | Colab L4 + Drive for the 26 GB weights |
| 2–3 | **Register on MatwingsVenus and AsterFire** and read off their ESMFold2/ipSAE job defaults. These are free during the contest and host BindCraft/BoltzGen/ESMFold2/ipSAE — they are the cheapest large-scale compute we have. Simultaneously post the five gate questions (§1.6) in the group chat. | sponsor platforms |
| 3–10 | **Generate.** Primary: **BoltzGen** with the hotspots from Q3/Q4 (on a sponsor platform or a rented L40S — it will OOM on L4). Parallel: **BindCraft on Colab L4** (fits: 259–319 aa complexes, ~138 s/hallucination, so ~25 designs/L4-hour → ~600 designs in a 24-hour-equivalent budget per target). Add **BoltzDesign1 Colab** purely for topology diversity. Push binder length to **80–100 aa** (see Q6 §6.2). | mixed |
| 10–15 | **Score everything through the gate reproducer**, rank by `min(ipSAE_A→B, ipSAE_B→A)`, keep anything ≥ 0.55. Secondary orthogonal filters (because the gate alone is exploitable): AF2/AF3 or Boltz-2 `pae_interaction < 10`, `plddt_binder > 80`, `rmsd_binder < 2 Å`, Rosetta/`freeSASA` interface ΔG/ΔSASA, shape complementarity > 0.62. | Colab L4 |
| 15–18 | **Select the 20 per track for a Pareto front, not a single metric** (Q6 §6.5): ipSAE high, interface ΔG/ΔSASA negative, ΔSAP (hydrophobic patch) moderate not extreme, ≥4 distinct topologies, ≥3 distinct epitope sub-patches, no unpaired cysteine, no N/Q deamidation or DP/DG motifs, pI away from 7, and **novelty check vs 7S5B + Cao 2022** for IL-7Rα. | CPU |
| 18–19 | Write the design summary PDF (≤4000 words), upload, confirm. Captain only; confirm is irreversible per track. | — |

**Why BoltzGen primary rather than the official ESMFold2 inversion:** the ESMFold2 inversion would be ideal (it directly optimises the scorer and has 53.8–70 % hit rates) but needs ≥27 GB and an H100, and we have neither locally. BoltzGen's "≤20 designs for 66–80 % success" is the only published protocol matched to a 20-sequence budget, and its weights are free. **Why BindCraft as the L4 workhorse:** it is the only hallucination method that demonstrably fits in 16–20 GB at our complex sizes *and* it is the generator that actually produced the highest ESMFold2-ipSAE scores in the one public dataset we have.

**Do not** rely on RFdiffusion3/RFD4 + ProteinMPNN as the primary route: 0.8 % of its designs cleared ipSAE 0.6 in the public dataset versus 6.0 % for BindCraft2, at 56.4 GPU-h vs 12.3 GPU-h. Keep it only as a diversity source if GPU time is free.

## 5.3 Colab notebooks, with 2026 maintenance status

| Tool | Notebook | Status |
|---|---|---|
| ESMFold2 local GPU | [`esmfold2_local_gpu.ipynb`](https://github.com/Biohub/esm/blob/main/cookbook/tutorials/esmfold2_local_gpu.ipynb) — has a `# Colab: uncomment` cell | **Current** (dated 2026-08-25 in the notebook) |
| ESMFold2 full tutorial | [`esmfold2.ipynb`](https://github.com/Biohub/esm/blob/main/cookbook/tutorials/esmfold2.ipynb) | Current |
| ESMFold2 binder design | [`binder_design.ipynb`](https://github.com/Biohub/esm/blob/main/cookbook/tutorials/binder_design.ipynb) | Current but **Modal/H100-targeted** |
| Apple Silicon | [`esmfold2_local_applesilicon.ipynb`](https://github.com/Biohub/esm/blob/main/cookbook/tutorials/esmfold2_local_applesilicon.ipynb) | Current (MLX) |
| BindCraft | [`notebooks/BindCraft.ipynb`](https://github.com/martinpacesa/BindCraft/blob/main/notebooks/BindCraft.ipynb) | **Maintained** (Nature 2025; repo active) |
| BoltzDesign1 | [Colab](https://colab.research.google.com/github/yehlincho/BoltzDesign1/blob/main/Boltzdesign1.ipynb) | Public, maintenance **UNCONFIRMED** |
| BoltzGen | no official Colab; CLI only, [repo](https://github.com/HannesStark/boltzgen) | Active (Nov 2025) |
| RFdiffusion / ColabDesign | `sokrypton/ColabDesign` RFdiffusion notebook | Still works; the maintained RFdiffusion3/RFD4 path is **local, via `RosettaCommons/foundry`** — no official Colab |
| ipSAE | none needed — pure numpy, run in any cell | v1 stable, v2 alpha |

---

# Q6 — Winning the gate *and* the KD ranking

The gate is deterministic (1 seed, 1 sample, no template, fixed MSA policy) so it is fully reproducible locally, and selecting our 20 by it is legitimate optimisation. But ipSAE has structural biases, and the final prize is measured KD. Here is where they agree and where they fight.

## 6.1 The decisive bias: ipSAE rewards a LARGE low-PAE neighbourhood, via `d0`

`d0res_i = max(1.0, 1.24·(n0res_i − 15)^(1/3) − 1.8)` where `n0res_i` = the number of residues in the scored chain with `PAE_i,j < 10` as seen from aligned residue `i`. Larger `n0res` ⇒ larger `d0` ⇒ `1/(1+(PAE/d0)²)` closer to 1 **for the same PAE**. Computed from the source:

| `n0res` | `d0` | ptm @PAE=3 | @5 | @8 | @10 |
|---|---|---|---|---|---|
| ≤25 | 1.00 (floor) | 0.100 | 0.038 | 0.015 | 0.010 |
| 30 | 1.26 | 0.150 | 0.060 | 0.024 | 0.016 |
| 40 | 1.83 | 0.270 | 0.118 | 0.050 | 0.032 |
| 60 | 2.61 | 0.431 | 0.214 | 0.096 | 0.064 |
| 100 | 3.65 | **0.597** | 0.348 | 0.172 | 0.118 |
| 180 | 5.00 | **0.735** | 0.500 | 0.281 | 0.200 |

**Two hard consequences:**
1. **Below `n0res`≈25 the `d0` floor of 1.0 Å makes a high ipSAE essentially unreachable** — even a perfect PAE of 3 Å caps the per-pair term at 0.10. **A small, tight interface is structurally penalised by this metric.** A 40 aa binder with a ~20-residue footprint is in the dead zone.
2. ipSAE is a **max over aligned residues**, so the score is set by the *single* binder (or target) residue with the broadest, most confident set of cross-chain partners — not by the mean quality of the interface. One residue sitting in the middle of a large, confidently-placed patch is worth more than a uniformly decent interface.

**Therefore: use the length allowance. Design 80–100 aa binders, not 40–60 aa, and aim for a wide, flat, confidently-packed epitope footprint.** This is the clearest legitimate lever we have, and it is cheap. Caveat: `n0res` counting admits residues up to PAE 10 Å, so blindly growing the halo also drags the *mean* down — the optimum is a compact high-confidence core (PAE 2–4 Å) embedded in a broad 60–120-residue PAE<10 neighbourhood. Empirically the achievable ceiling is ~0.70 (the competitor's 1,127-design max was 0.698).

## 6.2 Where ipSAE and real KD agree (good news)

- The 3,766-binder meta-analysis puts **`ipSAE_min` first among all tested in-silico features** for predicting experimental success — ahead of AF2 `pae_interaction`, ipTM, pLDDT and Rosetta energies. So optimising ipSAE is *not* orthogonal to real binding; it is the best single proxy anyone has published.
- ESMFold2's own designed minibinders (ranked by iPTM, a close cousin) achieved **68 pM – 70 nM** measured affinities at 53.8–70 % hit rates. The scorer's inverse really does produce binders.
- In the competitor's ESMFold2 data, `ipsae_min` correlates with `pose_dockq` at **r = 0.56** — it is tracking real pose correctness, not just confidence.

## 6.3 Where they fight (the hedges)

| Bias / failure mode | Evidence | Hedge |
|---|---|---|
| **Small tight interfaces are under-scored** (`d0` floor, §6.1) | derived from `ipsae.py` source | Go long (80–100 aa), broad footprint. Costs nothing. |
| **Interface hydrophobicity is where the metric is blind.** The meta-analysis found `ipSAE_min` average precision **"did not correlate with interface hydrophobicity (ΔSAP), despite its known association with binder success"** | [meta-analysis](https://www.biorxiv.org/content/10.1101/2025.08.14.670059v1.full) | ipSAE will not reward or punish hydrophobicity → **we must impose it ourselves**. Compute ΔSAP / interface ΔG/ΔSASA independently and require real buried hydrophobic surface. Both our epitopes offer one (TL1A W137/F138/M87/L52; IL-7Rα F79/L80/L81/I82). |
| **Exposed-hydrophobic false positives**: "membrane proteins can form apparently confident interfaces between exposed hydrophobic surfaces" and "small interfaces may be unstable or nonspecific" | [AF2/AF3 scoring-metric assessment](https://www.biorxiv.org/content/10.1101/2025.04.16.648930v1.full.pdf) | Cap binder surface hydrophobicity; reject designs whose own monomer has a large exposed apolar patch (these aggregate and fail developability, which is also worth +5k RMB). |
| **Idealised helical bundles are the known over-represented, over-confident class** in AF2-based design; the metric does not penalise them, and the IL-7Rα novelty rule bites hardest there (7S5B is a 3-helix `EKLRKLEK…` bundle) | §4.4; [mapping the evolution of designed binders](https://www.biorxiv.org/content/10.1101/2025.10.04.680454v1.full) | Deliberately diversify topology: mix 3-helix, 4-helix, helix-loop-helix-β, and 70–100 aa folds. Kills two birds (novelty + hedge). |
| **Confidence is not affinity.** "Although most sequences show relatively low ipTM … binders tend to have higher average ipTM than weak binders and non-binders" — i.e. the metric separates *binder vs non-binder*, not *10 nM vs 100 pM*. Final ranking is by KD alone. | [GFRAL/AF3 analysis](https://arxiv.org/pdf/2506.01995) | Do **not** submit the 20 highest-ipSAE designs. Submit a Pareto set (§6.5): clear the gate with margin, then maximise *predicted affinity* proxies (interface ΔG/ΔSASA, buried apolar area, H-bond count) within the surviving pool. The meta-analysis's best *combined* feature was exactly `ipSAE_min × interfaceΔG/ΔSASA < −1.5`. |
| **Target-dependent calibration.** `ipSAE_min` precision ranged **0.1 → 1.0 across targets**; predictability was best for "targets with high energy density … higher H-bond percentage or more negative ΔG/ΔSASA" | meta-analysis | Our TL1A groove is a classic high-energy-density TNF-family cleft (good). IL-7Rα's `FLLIG` strand is also high density (good). Both should be in the favourable regime. |
| **Determinism ≠ robustness.** 1 seed, 1 sample means the gate is a single lottery draw. In the competitor's re-evaluation, within-parent σ was 0.0046 vs overall σ 0.1463 (**32× ratio**) — so the ESMFold2 score is highly reproducible *for a fixed sequence* (good: our local reproduction will match) but varies enormously *across sequences from the same parent trajectory* (so sequence choice, not backbone choice, is where the score is won). | [competitor PR #5](https://github.com/RFjason/binder-design/pull/5) | Do many MPNN/sequence variants per accepted backbone and score every one. This is the single cheapest way to raise the max. |
| **Glycans invisible to the gate** | §4.3 | Penalise contacts near TL1A N62/N158 and IL-7Rα N29/N45 even though ipSAE will not. |
| **The His-tag is invisible to the gate but present in the assay** (TL1A N-term His, IL-7Rα C-term His) | competition spec | Avoid epitopes immediately adjacent to the tagged terminus. For IL-7Rα the C-terminus is the D2 end (construct ~210–219) → the D2-fringe residues 192/193/212/213 are near the C-term His; another reason to keep the primary hotspot in D1. For TL1A the N-term His is at construct position 1, which is the flexible stalk region — low risk. |

## 6.4 ⚠️ The biggest unresolved gate question: how many target chains does the gate fold for TL1A?

The rule says MSA is used for "仅靶点序列" — *the target sequence*, singular. For ipSAE you need ≥2 chains. For the five monomeric/VHH tracks the obvious implementation is **1 target chain + 1 design chain**. If the organiser does the same for TL1A, **the gate folds a TL1A MONOMER + binder**, and:
- The groove epitope (which requires two protomers, §3.4) **cannot exist** in the prediction. A groove-targeted design will fold against a half-site and score badly.
- The "right" design for the gate would target a **single-protomer face**, which may be a worse real epitope.
- Conversely, if the gate folds the trimer (3× target chain + binder = 580–640 aa), ipSAE is computed *per chain pair*, so a groove binder splits its contacts across two protomers and **each pairwise ipSAE is depressed** — unless the organiser uses IPSAE2's `--group A+B+C`, which v1 cannot do at all.
- Either way, **pairwise ipSAE structurally disfavours the therapeutically validated inter-protomer epitope on TL1A.** This is the sharpest gate-vs-biology conflict in our two tracks.

**Hedge for TL1A (do this):** split the 20 submissions — roughly **12 designs against a single-protomer face** (hotspots `A32,A49,A51,A52,A85,A86,A87,A117,A119,A120`, which are ≥9 Å from any neighbouring protomer and therefore score identically whether the gate folds a monomer or a trimer) and **8 designs against the groove** (`A104,A105,A115,A116,A117,A119,A121` + `A136,A137,A138`) which are biologically better and will win the KD ranking if they survive. Score every candidate **both ways** locally (monomer+binder and trimer+binder) and require that it clears the bar in **both** configurations — that makes us config-agnostic. Ask the organiser directly; this single answer is worth more than a week of compute.

## 6.5 Selection rule for the final 20 (per track)

```
HARD GATE   min(ipSAE_A→B, ipSAE_B→A) ≥ 0.61   [PAE cutoff 10, dist 10, d0res variant]
            scored in BOTH monomer and trimer configs for TL1A
HARD GATE   novelty: ≤70% identity to 7S5B / Cao-2022 designs over any ≥80%-coverage alignment
HARD GATE   length 40–100 aa (target 80–100), 20 standard AAs, no Cys unless deliberately paired,
            no signal peptide / tag / fusion
ORTHOGONAL  AF2-initial-guess or Boltz-2: pae_interaction < 10, plddt_binder > 80, rmsd_binder < 2 Å
ORTHOGONAL  interface ΔG/ΔSASA < −1.5 (combined feature ipSAE_min × ΔG/ΔSASA < −1.5)
ORTHOGONAL  shape complementarity > 0.62
DEVELOP.    moderate ΔSAP (no large exposed apolar patch), pI ∉ [6.5, 7.5], no NG/DG/DP motifs,
            no free Cys, no contacts within 8 Å of glycan Asn (TL1A N62/N158; IL-7Rα N29/N45)
DIVERSITY   ≥4 distinct topologies, ≥3 epitope sub-patches, pairwise identity among our own 20 < 60%
```
Rationale for diversity: ranking uses only each team's **single best** sequence per track, and KD is measured not predicted. A portfolio of 20 uncorrelated bets maximises P(at least one picomolar hit) far better than 20 near-copies of the top-ipSAE design.

---

# UNCONFIRMED / NEEDS FOLLOW-UP

**Ask the organiser (highest value first):**
1. **Exact ESMFold2 checkpoint** — `biohub/ESMFold2`, `ESMFold2-Experimental`, `ESMFold2-Experimental-Cutoff2025`, or platform-only "ESMFold2 Pro / `esmfold2-2026-05-cutoff-2025`"? Different cutoffs → different scores.
2. **How many target chains** are folded for the TL1A track — 1 (monomer) or 3 (trimer)? §6.4 — this changes which epitope we should target.
3. **ipSAE version and column**: `DunbrackLab/IPSAE` v1 or `IPSAE2` v0.1.10? Which variant — `ipSAE` (`d0res`), `ipSAE_d0chn`, or `ipSAE_d0dom`? Which row — `asym` (and then min or a fixed direction) or `max`? For a trimer, is `--group` used?
4. **`pae_cutoff` / `dist_cutoff`** — 10/10 or 15/15?
5. **`num_loops` / `num_sampling_steps`** — the `config.json` default is `num_loops: 20`; `num_sampling_steps` is not defaulted in the config.
6. **The exact a3m** used for each target (database, depth, `msa_max_depth`, `msa_column_mask_rate`). MSA depth materially changes the prediction, and we cannot reproduce the gate without matching it.
7. Will organisers publish the per-sequence scores or just the pass list?

**Could not verify:**
- **MatwingsVenus WeChat announcement** `mp.weixin.qq.com/s/5ef89061e1247218f0d0432148385a89` → returns **参数错误**; not mirrored anywhere I could find. Need a WeChat session or an in-platform look at `https://matvenus.com/zh`.
- **AsterFire's own page** `https://sidereus-ai.com/exhibition/kits` → **HTTP 502**. The "ipSAE 0.637 / 0.628, binder pLDDT 91.8" TL1A numbers reached me only via search-result summaries and a third-party GitHub PR, **not from AsterFire directly**. Treat as indicative, not authoritative.
- **Measured peak VRAM for ESMFold2 in bf16 / fp8** at 280 aa and 640 aa. Only the maintainer's rule of thumb ("~500 aa on 24 GB") exists. **Measure on first run.**
- **Whether `low_memory_mode` (PR #385) is merged** — it was still open on 2026-11-11. Re-check before relying on it.
- **Biohub platform free-tier quota** for ESMFold2 folds. If generous, the Batch API (4096-token limit) could do all our scoring without any local GPU.
- **Tulisokibart and duvakitug epitopes** — not public. The only structurally defined anti-TL1A epitope I found is BMS 10A4.F7 (US10968279B2). The 2-residue offset between the patent numbering and UniProt O95150 is resolved by motif matching, but the patent's SEQ ID 19 itself was not inspected.
- **Existence of a TL1A:DR3 experimental structure.** I found none; all DR3-interface claims in the literature derive from models built on 3K51. Also **not verified**: which single TL1A residue is "critical for DR3 but not DcR3" — the Sci Rep 2022 engineered-TL1A paper has it, but `ncbi.nlm.nih.gov` was blocked from this environment. Retry via [nature.com/articles/s41598-022-24984-y](https://www.nature.com/articles/s41598-022-24984-y).
- **ipSAE2's protein–protein thresholds.** Retrieved on retry: the direction convention (align on the more complex molecule) is confirmed, but the published cutoffs (0.5 AF3 / 0.8 Boltz-2) are **protein–ligand** only. No protein–protein threshold is given.
- **DR3-specific TL1A residues.** Resolved via the AKALFF mutant (§3.3) — but that identifies *DcR3-weighted* positions by subtraction, not a directly measured DR3 footprint. A true TL1A:DR3 structure still does not exist, and `ncbi.nlm.nih.gov` / `europepmc.org` were blocked from this environment (403), so the per-mutation KD table in Sci Rep 12:20538 was not read. Retry for the exact DR3 and DcR3 KD values.
- **Full Cao 2022 IL-7Rα minibinder sequence list** — only 7S5B's 55-mer was retrieved. Pull the Nature supplementary tables for the complete novelty-check set.
- **EasyNano** VRAM/runtime/code availability (arXiv 2606.12772) — PDF text layer insufficient.
- **Competitor data caveat**: the 1,127-design ipSAE distribution in §2.6 and the GPU-hour table in §5.1 come from a **third-party public repo** ([RFjason/binder-design](https://github.com/RFjason/binder-design)) whose author's gate assumptions (num_loops=20, sampling=200, msa_max_depth=1024, min-of-directions, 10/10) are **their reconstruction, not the organiser's published config**. Their data is VHH-against-TNFR1, not minibinder — the minibinder distribution will be shifted upward.
- **Field size / where the top-1000 cut lands.** Unknowable until results.
