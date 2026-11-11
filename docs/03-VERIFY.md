# 03 — VERIFY: independent verification of the fact base

> Agent: VERIFY. Run date **2026-11-11** (site clock confirmed: `/api/public/meta` → `"now":"2026-11-11T12:31:18+08:00"`).
> Method: public, unauthenticated reads only. No login attempted, no credentials entered, no state-mutating endpoint touched.
> Artifacts produced: `data/structures/`, `data/known_binders/`, `data/submission_template.csv`, `data/validate_submission_csv.py`.

---

## 0. Corrections to `docs/00-COMPETITION-FACTS.md` — read this first

| # | Fact base says | Reality | Severity |
|---|---|---|---|
| C1 | "`11QC` có vẻ là typo của PDB … khả năng là `1IQC`?" | **WRONG. `11QC` is a real PDB entry.** "Human transferrin receptor ectodomain", cryo-EM **2.4 Å**, chains **A, B** (homodimer), UniProt P02786, released 2026-05-27, from *Nat Commun* 2026 "Transferrin receptor 1 binds human parvovirus B19 VP1u to facilitate entry". Matches the competition's "11QC chain A、B / 二聚体" exactly. `1IQC` is a bacterial di-heme peroxidase — unrelated. Remove the typo note. | Medium (credibility) |
| C2 | "Công bố danh sách vào vòng wet-lab: **cuối 11/2026**" | **WRONG.** Official rules §三: 入选名单锁定 = **2026-12-01 → mid-December 2026**. | Medium |
| C3 | "Công bố kết quả **1/2027**" | **WRONG.** Official rules §三: 实验验证 mid-Dec-2026 → **est. 2027-02**; 实验结果发布 **est. 2027-02**. The site's own schedule page also now says **2027年2月**. The "2027-01" figure comes from the stale `/api/public/meta` timeline. | Medium |
| C4 | Material 1 "hỗ trợ CSV/XLSX/XLS" — correct, but "header phải có name và sequence; tối đa 20 dòng" | Mostly right; see §2 for the precise semantics (max 20 **sequences** = 21 lines incl. header; headers are case-insensitive and order-free; extra columns ignored). | Low |
| C5 | Design summary "≤4000 từ khi confirm; ≤8000 ở chỗ khác" | Confirmed both limits exist in the bundle, and they are **characters (字), not words**. Server rejects >4000 in one code path and >8000 in another. **Treat 4000 characters as the binding limit.** | Low |
| C6 | Tổng prize pool ~1.4M RMB | Confirmed, and it reconciles exactly: (6 + 2.5×2 + 1×3) = 14万/track × 7 = 98万, plus max 1万 bonus × 42 winning designs = 42万 → **140万**. | — |
| C7 | Range TL1A "Leu72–Leu251" / IL-7Rα "Glu21–Asp239" | **Both verified correct and both are biologically meaningful boundaries** (see §4) — not arbitrary. | — |
| C8 | "⚠️ Cần verify trực tiếp khi login" re: deadline | Resolved without logging in — see §1. | — |

Everything else I checked in the fact base (7 tracks, UniProt IDs, ranges, oligomer states, tags, 20 seq/track, ESMFold2+ipSAE config, seed-team rule, 7000/1000 wet-lab numbers, prize structure, duplicate rule >90% coverage & 100% identity) **matched the official sources**.

---

## 1. T1 — The deadline. **Conclusion: 2026-11-30 23:59:59 CST is the real deadline, but you MUST be confirmed by 2026-11-15.**

**Confidence: high (~90%) that the governing deadline is Nov 30; high (~80%) that the submit button nevertheless breaks at midnight on Nov 15 unless the organizers patch a config row.**

### Evidence

**(a) No new JS bundle.** Re-fetched `https://aiprotein.sjtu.edu.cn/` — asset hash **unchanged**: still `assets/index-BZIjDIfc.js`. Response headers show `last-modified: Tue, 10 Nov 2026 02:36:07 GMT`, i.e. the bundle was *redeployed yesterday* with byte-identical content (Vite content-hashes filenames). So as of yesterday's deploy the frontend was not corrected.

**(b) The authoritative rules document says Nov 30.** `GET https://aiprotein.sjtu.edu.cn/api/pages/rules` (public, HTTP 200, `updated_at: 2026-11-06T15:36:33+08:00` — i.e. edited *after* the 11-01 target release). Its §三 timetable table:

| 阶段 | 时间 |
|---|---|
| 报名 | 2026-09-24 至 **2026-11-30 23:59:59** |
| 赛题发布、序列设计 | 2026-11-01 至 **2026-11-30 23:59:59**（…完成序列提交） |

and §二.4: "报名与成员锁定：…截止（锁定）时间为 **2026 年 11 月 30 日 23:59:59**（北京时间）". §三.1: "本规则全部时间为北京时间（UTC+8），**截止时刻以官方服务器记录为准**". This page is explicitly part of the participation contract ("本规则与《隐私政策与免责声明》共同构成参赛协议").

**(c) The website's own rendered schedule page says Nov 30.** Hardcoded in the current bundle's 赛事须知/赛程安排 component: `2026年9月24日 报名启动 | 2026年11月1日 靶点发布 | 2026年11月30日 报名截止 / 序列提交截止 | 2026年12月中 入围名单公布 | 2026年12月下旬 湿实验验证 | 2027年2月 赛事结果公布`, with prize pool 140万. So the *display* layer is already updated to Nov 30.

**(d) But the server config that gates the submit button still says Nov 15.** `GET https://aiprotein.sjtu.edu.cn/api/public/meta` returns, right now:

```json
"phase":"submission", "submission_open":true, "now":"2026-11-11T12:31:18+08:00",
"timeline":[ … {"key":"submission","start":"2026-11-01","end":"2026-11-15","state":"current"},
             {"key":"review","start":"2026-11-16","end":"2026-11-30","state":"upcoming"},
             {"key":"wetlab","start":"2026-12-01","end":"2026-12-31"},
             {"key":"results","start":"2027-01-01","end":"2027-01-31"} ]
```

**This is the mechanism behind the 2026-11-08 group-chat complaint ("15号就变灰了").** In the bundle, every control on the submission page is gated on `meta.submission_open` (`disabled: !n.submission_open`), and the i18n table contains the *server* error string `"已过提交截止时间"` → "The submission deadline has passed." So the lockout is enforced server-side off this config, not just cosmetically.

**(e) Supporting evidence that `/api/public/meta` and `/api/home-sections` are stale, not authoritative.** `GET /api/home-sections` → `home_rules` still reads "序列提交截止时间为 **11 月 15 日**（北京时间）" — but the *same* text also says "奖金总额约 **100 万元**" and "不少于 **5000 条** 序列" and "本次赛事计划设置 **5 个以上**靶点，涵盖 Anthropic 公布的…靶点", all of which are superseded (140万 / 7000条 / 7 targets). That block is pre-announcement CMS text nobody refreshed. `/api/public/meta` is *partly* refreshed (its `stats` correctly says 140 万元) and partly not (wet-lab still "5000+ 条", timeline still 11-15, results still 2027-01). Also note the bundle contains a hardcoded *fallback* phase table (`submission` until `2026-11-15T23:59:59+08:00`) used only when meta is unavailable — same stale numbers.

**(f) The WeChat announcement.** `instruction.md` 2026-10-23, from the official account 长江杯AI蛋白质设计挑战赛: "⏰ 利好延期：报名&序列提交截止至**11月30日**". Consistent with (b) and (c). I could **not** independently retrieve the mirrored WeChat article — `mp.weixin.qq.com` returns an anti-bot page ("环境异常，完成验证后即可继续访问") to non-browser clients. **UNCONFIRMED by independent fetch**, but corroborated by (b) and (c), which are stronger sources anyway.

### Practical implication — act on this

1. **Treat 2026-11-15 23:59:59 CST as your hard deadline.** It is 4 days away. The Nov-30 date is almost certainly the real contractual one, but the submit path is gated by a config row that currently flips on 2026-11-16 00:00, and the organizers have had this reported since 11-08 and have not fixed it in the 11-10 redeploy.
2. **A saved draft probably does not count.** Rules §三.3: "以截止前最后一次**经系统成功确认**的完整版本为准" — the *confirmed* version is what counts. The UI distinguishes 保存 (draft, overwritable) from 确认提交 (locks the track). So you must reach the **confirmed** state, not merely saved.
3. **Confirming is reversible but awkward.** UI/server say "该赛题已确认提交，不能再修改", yet the admin-side string reads "已确认的序列不能退回草稿。如需更正，请重新确认并填写说明。" → re-confirmation with a mandatory note (field `maxLength: 2000`, `required`) is supported. So confirming early is recoverable, but likely needs organizer involvement.
4. **Recommended plan:** freeze candidate sets by 2026-11-13, confirm both tracks on 2026-11-14 (leaving one buffer day for site problems). If the window genuinely runs to Nov 30, request a re-confirm to upgrade. Do not gamble on Nov 30.

---

## 2. T2 — Submission format, extracted from the bundle's validation code

All of the below is read directly out of `assets/index-BZIjDIfc.js` (validator function `Bv`, CSV field splitter `tc`, length-spec builder `ec`, extension checks `Uv`/`qv`, spreadsheet reader `Vv`, constants `Lv`/`Xu`/`Rv`/`Iv`/`Zu`).

### The authoritative numbers come from the server, not the page text

`/api/public/meta` carries the live limits, and the page reads them:

```
max_sequences: 20
vhh_min_length: 100,       vhh_max_length: 150
minibinder_min_length: 40, minibinder_max_length: 100
```

**Resolution of the "40–180 / 100 条" confusion:** the string `每条序列长度 40–180，每队最多 100 条。` is a **hardcoded static hint string** in the submission-page UI (it has an English i18n twin: "Each sequence must be 40–180 residues long, with at most 100 sequences per team."). It is **not** derived from `meta` and matches nothing — not minibinder (40–100), not VHH (100–150), not the sequence cap (20). It is dead placeholder text. **The enforced values are the `meta` values, and they agree with the official rules §四(一) table (minibinder 40—100 aa, VHH 100—150 aa) and §一(三).3 (≤20 sequences per team per track). Ignore 40–180/100 entirely.**

### Exact validation semantics

| Item | Rule |
|---|---|
| File types | `.csv`, `.xlsx`, `.xls` (extension check, case-insensitive). XLSX/XLS are read with SheetJS, **first sheet only**, converted to CSV internally. |
| Max file size | **10 MB**, checked **twice**: on the raw file, and again on the UTF-8 byte length of the converted CSV text. |
| Encoding | **UTF-8.** A leading BOM (`﻿`) is stripped. Line endings `\r\n` or `\n` both fine. (GB18030 is accepted only on an unrelated *admin* roster import, not here.) |
| Required headers | Row 1 must contain **`name`** and **`sequence`**. Each header cell is `.trim().toLowerCase()`-ed → **case-insensitive** (`Name`, `SEQUENCE` fine). **Order does not matter** (looked up by index). |
| Extra columns | **Ignored.** Confirmed by code (only the `name`/`sequence` column indices are read) and by the UI hint: "UTF-8。第一行表头至少包含 name 和 sequence，**多余列会忽略**。" |
| Max rows | `最多提交 20 条序列（含表头共 21 行）`. Two checks: non-blank line count > 21 → reject; and a running data-row counter > 20 → reject. |
| Blank rows | A row where **both** name and sequence are empty is silently skipped. Fully blank lines skipped. |
| `name` | **1 to 64 characters inclusive** (`!x || x.length > 64`). Trimmed. Measured in JS UTF-16 code units — keep names ASCII to avoid surprises. |
| Duplicate names | Rejected within a file (`名称重复`). Rules also require design IDs unique. |
| `sequence` preprocessing | **All whitespace removed** (`/\s+/g` → ``) and **uppercased** — so lowercase input and internal spaces/newlines are tolerated. |
| Allowed characters | `new Set("ACDEFGHIKLMNPQRSTVWY")` — the **20 standard amino acids only**. No `X`, `B`, `J`, `O`, `U`, `Z`, no `*`, no `-`. Anything else → `含有非法氨基酸`. **Yes, there is client-side sequence character validation.** |
| Length bounds | `S.length < min || S.length > max` → fail. So **bounds are INCLUSIVE**: for minibinder, 40 aa and 100 aa both PASS. Independently corroborated by the admin UI hint for the same settings: "含最短和最长，**两端都算合格**。" |
| Atomicity | Any error → nothing is written: "本地核对未发现错误。**有错误时整份不会写入**。" |
| What is uploaded | For CSV the file itself is **not stored** — only the parsed text is posted (`POST /api/submissions/csv`): "文件不会保存到服务器，只提交其中的文本。" |
| One file per target | "按靶点上传 CSV 或 Excel。**每个靶点一份文件**…同一靶点重新上传需要填写说明。" |

### Design summary (Material 2)

- Text area `maxLength` = **8000** characters; hint "选填，最多 8000 字"; a second code path says "确认提交时必填，最多 **4000** 字". Server rejects both ">4000" and ">8000" in different handlers. **Stay ≤4000 characters to be safe.** Required at confirm time.
- Or a **PDF ≤20 MB**, one per target, re-uploadable (latest wins), posted base64.
- Re-upload note: required, `maxLength: 2000`.

### Template + validator (deliverables)

- **`data/submission_template.csv`** — header-only, exactly `name,sequence` + newline, plain UTF-8 (no BOM). Append up to 20 rows. Deliberately contains **no example row**, so there is nothing to forget to delete.
- **`data/validate_submission_csv.py`** — a faithful line-by-line reimplementation of the site's `Bv`/`tc`/`ec` validator (same BOM stripping, same quoted-field CSV parser with `""` escaping, same whitespace-strip/uppercase, same AA set, same inclusive bounds, same error strings). Run it before every upload:
  ```
  python3 data/validate_submission_csv.py your_file.csv
  ```
  Defaults are the minibinder settings (max 20, 40–100 aa). Output is either the parsed sequence count + `NONE -> would upload`, or the same `(line, name, message)` error tuples the site would show. Verified against the shipped template (→ `没有序列`, correct for an empty template) and against a 1-row 55-aa minibinder file (→ clean).

---

## 3. T3 — Reference structures (all in `data/structures/`)

### 2RE9 — TL1A reference. **Apo homotrimer. Directly usable; trim the N-terminal tag only.**

Source: `https://data.rcsb.org/rest/v1/core/entry/2RE9`, coordinates `https://files.rcsb.org/download/2RE9.pdb`.

- Title "Crystal structure of TL1A at 2.1 A", **X-ray, 2.1 Å**, 1 polymer entity, **3 chains A/B/C** = three copies of TL1A. *Confirmed TL1A / TNFSF15* ("TNF superfamily ligand TL1A", synonym "tumor necrosis factor (ligand) superfamily, member 15", keywords "VEGI, homotrimer"). Citation: *Biochem. Biophys. Res. Commun.* 2007, "X-ray crystal structure of TNF ligand family member TL1A at 2.1 A".
- **It is an apo trimer, NOT a complex with a binder.** The only heteroatoms are 1 Mg²⁺, 12 glycerols and waters. (RCSB's `pdbx_keywords` reads "HORMONE/GROWTH FACTOR COMPLEX" — that is a classification artefact, there is no second molecule.) This matches the competition's "多聚体形态: 三聚体".
- **Numbering.** `COMPND … FRAGMENT: SOLUBLE PART: RESIDUES 72-251`. `DBREF 2RE9 A 1 180 UNP Q8NFE9 … 72 251` — the file's DBREF points at the TrEMBL accession **Q8NFE9**, but RCSB/SIFTS maps the entity to **O95150 residues 72–251, length 180, entity coverage 0.994**. So: **PDB residue *i* = UniProt O95150 residue *i* + 71.** There is also a `GLY 0` expression tag (`SEQADV … EXPRESSION TAG`).
- **Sequence identity to our construct: exact.** I diffed the 2RE9 SEQRES (minus the Gly0 tag) against `data/TL1A_O95150_L72-L251.fasta`: **180/180 identical, zero mismatches.** Wild type — both native cysteines present.
- **Observed coordinates: PDB 17–180 = UniProt 88–251**, i.e. **164 of 180 residues (91.1%), with NO internal gaps.** Missing: the Gly0 tag and PDB 1–16 = **UniProt 72–87** (disordered N-terminal stalk), identically in all three chains.
- **Why the missing stretch doesn't matter:** UniProt annotates the **THD (TNF homology domain) as residues 95–251**. Residues 72–94 are the pre-THD stalk, which is genuinely disordered — that is why it's absent. The folded, bindable target surface is fully present.
- **Verdict: directly usable as RFdiffusion target input.** Keep all three chains (the TNF-superfamily receptor site lies across adjacent protomers), delete HETATMs (Mg, GOL, HOH) and the Gly0 tag, and remember your residue offset (+71) when mapping hotspots back to competition numbering.

### 3DI3 — IL-7Rα reference. **Complex. Chain B is the receptor. Usable, but fix one residue.**

Source: `https://data.rcsb.org/rest/v1/core/entry/3DI3`, coordinates `https://files.rcsb.org/download/3DI3.pdb`.

- Title "Crystal structure of the complex of human interleukin-7 with **glycosylated** human interleukin-7 receptor alpha ectodomain", **X-ray, 2.9 Å**. Citation: "Structural and Biophysical Studies of the Human IL-7/IL-7Ralpha Complex", *Structure* 2009. **Confirmed: it IS the IL-7/IL-7Rα complex.**
- **Chain assignment — the competition's "3DI3 chain B" is correct, chain B is the RECEPTOR:**
  - **Chain A = Interleukin-7** (the cytokine), UniProt **P13232**, `FRAGMENT: UNP RESIDUES 26 TO 177`, 152 aa mapped, engineered mutation **E131A** (`SEQADV ALA A 106 ← GLU 131`).
  - **Chain B = Interleukin-7 receptor subunit alpha** (CD127), UniProt **P16871**, `FRAGMENT: UNP RESIDUES 21 TO 239 (LIGAND BINDING ECTODOMAIN)` — **exactly the competition's Glu21–Asp239 range.**
- **Numbering.** `DBREF 3DI3 B 1 219 UNP P16871 … 21 239`, plus a **GSHM** expression tag at PDB −3…0. So **PDB residue *i* (chain B) = UniProt P16871 residue *i* + 20.**
- **⚠ One sequence difference from the competition construct.** Diffing 3DI3 entity-2 SEQRES (minus GSHM) against `data/IL7Ra_P16871_E21-D239.fasta`: **218/219 identical, one mismatch — UniProt residue 66: our construct has Ile, 3DI3 has Thr (I66T).** The PDB's `SEQADV` records do *not* flag it (they flag only `VAL B 118 ← ILE 138`, which no longer differs because P16871 was revised to SV=3 in 2018). The same T66 appears in 3UP1, so it is a property of that lab's construct, not a typo. Functionally near-silent (N65-I66-T67 and N65-T66-T67 are both N-glycosylation sequons) — but **if you use 3DI3 chain B as a design target, mutate T66→I to match the competition construct**, or at minimum do not place a hotspot there.
- **Observed coordinates, chain B: PDB 17–209 = UniProt 37–229**, i.e. **193 of 219 residues (88.1%), with NO internal gaps.** Missing: GSHM tag, PDB 1–16 = **UniProt 21–36** (N-terminal), and PDB 210–219 = **UniProt 230–239** (C-terminal, which includes glycosylation sites N232/N233).
- Chain A (IL-7) is much worse: observed 7–98 and 122–147 only — a **23-residue internal gap** (and 36 missing residues total). Irrelevant if you only use chain B.
- **Glycans are modelled**, as separate NAG chains C/D/E. I verified attachment by distance: **NAG C → Asn131(B) = UniProt N151; NAG D → Asn45(B) = UniProt N65; NAG E → Asn29(B) = UniProt N49.** None of the three sits in the IL-7 interface, but N151 is 7 residues from interface residue K158 and will partially shade the D2 side of the epitope in reality.
- **Verdict: usable but needs cleanup.** Strip chain A (unless you want it as a competition reference for the blocking epitope), strip NAG/HOH, strip the GSHM tag, mutate T66→I, and remember the +20 offset. **Consider 3UP1 or 7OPB instead (below) — both are better inputs.**

### Why I also fetched 2RJL, 3K51, 3UP1, 7OPB

I ran a systematic RCSB search for every structure mapping to each target accession (`search.rcsb.org/rcsbsearch/v2/query`, attribute `rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession`):

- **O95150 (TL1A) → 7 entries total: 2O0O, 2QE3, 2RE9, 2RJK, 2RJL, 3K51, 3MI8.**
- **P16871 (IL-7Rα) → 8 entries total: 3DI2, 3DI3, 3UP1, 5J11, 6P50, 6P67, 7OPB, 9XEA.**

| PDB | What it is | Why it matters |
|---|---|---|
| **2RJL** | TL1A ectodomain **monomer**, X-ray **2.05 Å** — the highest-resolution TL1A. Construct is UniProt 68–251 (184 aa incl. GSHM tag). **Double mutant C95S/C135S in their local numbering = UniProt C162S/C202S** — i.e. the native disulfide is engineered away. (2QE3 = WT monomer 2.5 Å; 2RJK = single C95S.) | Higher resolution than 2RE9 but **monomeric and disulfide-free** → wrong oligomeric state for a trimer target. Use 2RE9 for design; 2RJL is a useful cross-check / loop-modelling donor. Also tells you the C162–C202 disulfide is dispensable for folding. |
| **3K51** | **DcR3–TL1A complex**, 2.45 Å. Chain A = TL1A (UniProt 72–251, PDB offset −67), chain B = DcR3/TNFRSF6B (O95407, UniProt 30–195). (3MI8 is the same complex at 2.95 Å.) | **The only TL1A-with-a-binder structure that exists.** I computed the interface at 4.5 Å heavy-atom cutoff: TL1A residues contacting DcR3 = **UniProt 175, 176, 186–192, 221**. Caveat: 3K51 is a 1:1 slice of a 3:3 complex, so this is one protomer's share; the true receptor site in the trimer spans two adjacent protomers. This is your best proxy for the DR3/TNFRSF25 functional epitope. |
| **3UP1** | **Unliganded (apo) IL-7Rα ectodomain**, X-ray **2.15 Å**, 2 copies (chains A, B). Same construct as 3DI3 chain B: GSHM + UniProt 21–239, and it carries the same **I66T**. | Better resolution than 3DI3 (2.15 vs 2.9 Å) and apo. Good alternative target input, same caveats (T66, +20 offset). |
| **7OPB** | **"IL7R in complex with an antagonist"**, X-ray **2.144 Å**. Entity 1 = IL-7Rα, chains A/B/C, **UniProt 37–231, wild-type (has I66!)**. Entity 2 = **"IL7R binder", 55 residues, chains D/E/F**, expressed in *E. coli*, no UniProt. Citation: **Cao, Coventry, Goreshnik, Huang et al., "Design of protein-binding proteins from the target structure alone", Nature 605:551–560 (2022)**, https://doi.org/10.1038/s41586-022-04654-9 | **The single most important structure for us.** It is (a) the best IL-7Rα target model — highest resolution, WT sequence, no engineered T66, apo-equivalent receptor conformation; and (b) **it shows a published de novo minibinder bound to the exact epitope we would target.** See §4. |

Also noted and not downloaded: 3DI2 (unglycosylated IL-7/IL-7Rα, 2.7 Å, two copies), 5J11 (TSLP/TSLPR/IL-7Rα ternary, 2.56 Å — a different, competing interface), 6P50/6P67 (murine anti-IL-7Rα Fabs 4A10 / 2B8 — sequences saved, see §4), 9XEA (NMR of the IL-7R **transmembrane** dimer, out of range).

### 11QC (TfR1, not our track) — see correction C1. Real entry, not a typo.

One footnote for completeness: 11QC's entity aligns to **P02786 starting at residue 122**, while the competition range is Cys89–Phe760 — so 11QC does not cover residues 89–121. Not our problem.

---

## 4. T4 — Target sequences: independently verified, both exact

I re-derived both slices from the UniProt REST API (`https://rest.uniprot.org/uniprotkb/O95150.json`, `.../P16871.json`) and compared byte-for-byte.

### TL1A — `data/TL1A_O95150_L72-L251.fasta` ✅ **CORRECT**

- O95150 canonical length **251**; slice `[71:251]` → **180 aa**, starts **L** (res 72 = Leu), ends **L** (res 251 = Leu, the protein's C-terminus). File content matches the slice exactly. Independently matches 2RE9 SEQRES exactly (§3).
- **Isoform ambiguity: none in practice.** O95150 is reviewed Swiss-Prot, SV=2 (2008-05-20). It *does* declare alternative splicing with isoforms O95150-1/-2/-3, but **-1 is the Displayed (canonical) sequence** and -2/-3 are only "Described" (they differ in residues 1–85 / 1–100, i.e. upstream of or overlapping the start of our range, and have no separate FASTA). The competition wrote "O95150" without a suffix; the canonical -1 is the right and only sensible choice.
- **Boundary semantics — important.** TL1A is a **type II membrane protein**, so there is **no classical signal peptide**. UniProt: Topological domain 1–35 cytoplasmic; **Transmembrane 36–56** ("Helical; Signal-anchor for type II membrane protein"); Topological domain 57–251 extracellular. **Chain 72–251 is annotated by UniProt as the "secreted form".** So the competition range is *not* "right after the signal peptide" — it is the proteolytically released soluble ectodomain, starting 15 residues inside the extracellular region. Worth stating correctly in any write-up.
- **THD (TNF homology domain) = 95–251.** Residues **72–94 are a flexible stalk outside the folded domain** and are exactly the region unresolved in 2RE9 (72–87). Do not expect a designed binder to engage 72–94 productively; do not place hotspots there.
- **N-glycosylation in range: N133, N229** (2 sites). My independent N-X-[ST] sequon scan over the 180-mer found exactly {133, 229} — perfect agreement with UniProt.
- **Disulfides in range: one, C162–C202.** My Cys scan found exactly two cysteines, at 162 and 202 → one intrachain disulfide, **no free cysteine**. Note 2RJL/3K51 engineer both to Ser, so it is not required for folding. Rules §四(一) permit retaining native disulfides **but they must be disclosed** — this applies to the *target*, not your design; your minibinder should ideally be cysteine-free.
- **Accessibility read:** glycans at N133 and N229 are on the folded THD and will be present in the expressed target; the DcR3/DR3 functional epitope (UniProt ~175–176, 186–192, 221, §3) is distinct from both. Epitope looks accessible.

### IL-7Rα — `data/IL7Ra_P16871_E21-D239.fasta` ✅ **CORRECT**

- P16871 canonical length **459**; slice `[20:239]` → **219 aa**, starts **E** (res 21 = Glu), ends **D** (res 239 = Asp). File content matches exactly.
- **P16871-1 IS the canonical sequence — CONFIRMED.** Reviewed Swiss-Prot, SV=3 (last sequence update 2018-09-12). Isoforms P16871-1/-2/-3/-4 exist; **-1 is "Displayed"**, the others are "Described" and all differ only from residue **237 onward** (alternative sequences at 237–459, 237–252/253–459, 293–459) — i.e. **essentially outside our 21–239 window**, so even isoform confusion would barely perturb the construct. The competition's explicit "-1" suffix is correct.
- **Boundary semantics — both ends are exactly right.** UniProt: **Signal peptide 1–20**; Chain 21–459; **Topological domain 21–239 extracellular**; **Transmembrane 240–264**. So **Glu21 is the first residue after the signal peptide, and Asp239 is the last residue before the TM helix.** The competition range is precisely the complete ectodomain. (Contrast with GLP-1R, where the organizers deliberately extended the range into the TM region per the 2026-11-06 changelog.)
- **N-glycosylation in range: N49, N65, N151, N182, N232, N233 — six sites.** My sequon scan found exactly {49, 65, 151, 182, 232, 233}. Agreement is exact. Three are experimentally resolved in 3DI3 (N49, N65, N151 — verified by NAG–ND2 distance, §3). **This is a heavily glycosylated target** — six sequons over 219 residues — and that materially constrains where a designed epitope can live.
- **Disulfides in range: three — C42–C57, C74–C82, C108–C118.** My Cys scan found exactly six cysteines at {42, 57, 74, 82, 108, 118} → three disulfides, **no free cysteine**. All three are in D1.
- **Domain boundaries D1 / D2 within 21–239:**
  - **D2 = residues 131–231**, annotated by UniProt as a **Fibronectin type-III domain**. CONFIRMED from UniProt.
  - **D1 ≈ residues 21–130** — an FnIII-like domain completing the cytokine-receptor homology module. UniProt does **not** annotate D1 explicitly, so the exact N/C boundaries of D1 are **UNCONFIRMED**; I infer 21–130 as "everything before the annotated FnIII domain". A short linker 232–239 follows D2 before the TM.
  - All three disulfides and 2 of 6 glycans (N49, N65) are in D1; N151/N182 are in D2; N232/N233 are in the post-D2 linker.
- **The IL-7 binding epitope (= the functional blocking site).** I computed the IL-7/IL-7Rα interface from 3DI3 at a 4.5 Å heavy-atom cutoff. **16 IL-7Rα residues contact IL-7 (UniProt numbering): S51, H53, L77, V78, E79, K97, K98, F99, L100, L101, I102, G103 (12 in D1) and K158, Y159, Y212, F213 (4 in D2).** The site is dominated by a hydrophobic patch in D1 (F99-L100-L101-I102) — a classic designable minibinder epitope — and **none of the six glycosylation sequons falls inside it**, so the epitope is glycan-free and accessible, with the caveat that the N151 glycan sits near the D2 contacts (K158/Y159).

---

## 5. T5 — **THE BIG ONE: novelty.** IL-7Rα is a minefield; TL1A is clean.

### 5.1 IL-7Rα: there is a large, public, published de novo minibinder corpus for this exact target

**IL-7Rα was one of the 12 targets in Cao et al., Nature 605:551–560 (2022), "Design of protein-binding proteins from the target structure alone"** (Baker lab) — https://doi.org/10.1038/s41586-022-04654-9. This is the single most consequential fact in this report.

**Provenance of `data/known_binders/IL7Ra_denovo_minibinder_designs_IPD_Cao2022.fasta` (30 MB, 122,273 records).** I downloaded the paper's own supplementary archive — `https://files.ipd.uw.edu/pub/robust_de_novo_design_minibinders_2021/supplemental_files/scripts_and_main_pdbs.tar.gz` (65 MB; the download instructions are printed in Supplementary Information, `41586_2022_4654_MOESM1_ESM.pdf`) — and extracted `supplemental_files/design_models_sequence/`, which contains per-target design sequence files. The three IL-7Rα files are:

| file | rows |
|---|---|
| `IL7Ra.seq` (de novo scaffold library) | 14,999 |
| `IL7Ra_graft.seq` (motif-grafted designs) | 100,000 |
| `IL7Ra_ssm.seq` (site-saturation variants of hits) | 9,130 |
| **total rows / unique binder sequences** | 124,129 / **122,273** |

Each line is `binder_seq  target_seq  design_name`; I kept the binder column, de-duplicated, and tagged each record with its source set and original design name. **Length range 36–67 aa, median 63 aa, and 99.99% fall inside the competition's 40–100 aa minibinder window.** The target sequence embedded in every line is IL-7Rα starting at UniProt 37 and carrying **T66** — i.e. they designed against the 3DI3/3UP1 construct.

**Which ones were experimentally validated** (`data/known_binders/IL7Ra_denovo_minibinders_VALIDATED.fasta`, 2 records). From Supplementary Table 1 ("Topologies, initial amino acid sequences, final optimized amino acid sequences and physicochemical properties of the de novo miniprotein binders for all 12 targets"), the IL-7Rα row reads: binder **IL-7Rα_mb**, topology **HHH** (three-helix bundle), **SC50 (yeast) = 3 nM**, **13 mutations** between initial and optimized:

```
IL7Ra_mb_original   SVKKKVRKVEKKARKAGDELAVLLARRVLEALEKGLVSEEDADESADRIEEALKK   (55 aa)
IL7Ra_mb_optimized  SVIEKLRKLEKQARKQGDEVLVMLARMVLEYLEKGWVSEEDADESADRIEEVLKK   (55 aa)
```

I cross-validated these three independent ways, and all three agree, so I am confident they are transcribed correctly: (1) the optimized sequence is **byte-identical** to 7OPB entity 2 SEQRES from the RCSB API; (2) the Hamming distance between the two is **exactly 13**, matching the paper's stated mutation count; (3) the original sequence is **present verbatim** in `IL7Ra.seq` from the IPD archive.

**Affinity and epitope.** The paper reports a **picomolar-affinity** IL-7Rα binder that **targets the IL-7 binding site and blocks IL-7-induced STAT5 signalling**; the co-crystal structure (7OPB) matched the design model. Supplementary Table 1 gives the yeast-display **SC50 = 3 nM** for the optimized design (SC50 is a display-based proxy, not the solution KD; the "picomolar" figure is the paper's own characterization of the matured binder). So: **the epitope the known minibinder occupies is the IL-7 site I mapped in §4 — S51/H53/L77-E79/K97-G103 in D1 plus K158/Y159/Y212/F213 in D2.** That is precisely the epitope any sensible team would choose, because it is the functional blocking site and the one designable hydrophobic patch.

### 5.2 How the novelty rule actually reads

Official rules §四(二).2, verbatim from `/api/pages/rules`:

- 覆盖总长度 > 90% 且一致性 = 100% → duplicate (antibodies compared on CDR only).
- **参考库：minibinder 比对已知蛋白数据库，VHH 比对已知抗体数据库；两类库均包含该靶点已发表论文与公开专利中的序列。**
- **minibinder：与每个靶点的已知 minibinder 进行比对，比对区域覆盖总长度的至少 80% 时，序列一致性须 ≤70%。**
- Thresholds are applied against a **frozen** reference library; post-hoc discoveries don't change the threshold but can be used in plagiarism checks; teams get a right of reply.

So the rule is a **conjunction**: you fail only if (aligned region ≥80% of 总长度) **AND** (identity >70%). Note that **"总长度" is not disambiguated** — it does not say whose total length. The parent's reading (our full length) is the natural one and is consistent with the duplicate clause, but **this ambiguity is worth an explicit question to the organizers**, because it determines whether the mitigation in §5.4 works.

Crucially, "该靶点已发表论文…中的序列" plausibly **includes the IPD supplementary archive**, since it is the published paper's own supplementary data. Whether the organizers actually ingested 122k untested design models, or only the one binder pair reported in the main text and Supplementary Table 1, is **UNCONFIRMED and is the main uncertainty in the risk estimate below.**

### 5.3 Quantified risk — I measured it

Because there was no alignment tool available I wrote a 4-mer-indexed prefilter plus an exact Needleman–Wunsch (match +1 / mismatch −1 / gap −4) and ran three experiments against the IPD corpus. Identity is reported over the shorter sequence.

**Experiment A — how self-similar is the published corpus?** 120 random designs from `IL7Ra.seq` vs 20,000 other designs from the same set:

> mean **78.0%**, median **80.0%**, p90 90.7%, p99 98.4%, max 98.4%. **85.0% of designs exceed 70% identity** to at least one other design in the corpus.

**Experiment B — the realistic baseline.** 120 designs that Cao et al. made for *different* targets (TrkA, FGFR2, Tie2, H3, InsulinR, IGF1R, PDGFR, EGFR — 494,696 sequences available) vs 12,000 IL-7Rα designs. These are *independent* designs (different target, different interface) that nonetheless came from the **same scaffold library and the same method** — the closest available proxy for "we run a standard pipeline and compare to the known corpus":

> mean 48.3%, median **43.0%**, p90 **69.9%**, p99 77.8%, max 79.7%. **9.2% exceed 70%**; 29.2% exceed 60%.

**Experiment C — composition-matched negative control.** The same 120 sequences, shuffled, vs the same 12,000:

> mean 33.0%, median 32.7%, p90 39.2%, **max 43.1%, 0.0% exceed 70%**.

**Interpretation.** The gap between B (median 43%, 9% over threshold) and C (median 33%, 0% over threshold) is real signal: it is **shared scaffold provenance**, not just the E/K/A/L-rich compositional bias of idealized helical bundles. Experiment A's 85% figure is inflated by near-sibling designs (grafted variants and point mutants of the same scaffold) and is *not* the number that applies to us — but it is the number that would apply if we reused a public scaffold library.

### 5.4 My read: **material, not negligible — roughly 1-in-10 per sequence — and cheaply mitigable**

- **If the frozen library contains only the 2 validated binders:** risk is **near zero**. Two specific 55-mers; a fresh RFdiffusion backbone + ProteinMPNN sequence will land ~30–45% identical. Don't worry.
- **If the frozen library contains the full 122k IPD design set:** Experiment B says **~9% of designs from a comparable pipeline would exceed 70% identity** over ≥80% coverage. On 20 submitted IL-7Rα sequences that is **~2 expected rejections** — and those would be silent losses of wet-lab slots. **That is material.** The risk concentrates exactly where we want to be: everyone aiming at the IL-7 blocking epitope with a 3-helix bundle of ~55 aa is walking into the densest part of that corpus.
- **The risk is much worse if you reuse public scaffolds.** If any part of the pipeline uses the Baker-lab miniprotein scaffold set (or BindCraft/RFdiffusion runs that happen to regenerate those idealized HHH topologies at the same length), you move from Experiment B (9%) toward Experiment A (85%). **Do not warm-start, graft, or scaffold-match against anything from `robust_de_novo_design_minibinders_2021`.**

**What we have to do to be safe — concrete, in priority order:**

1. **Self-check every candidate against `data/known_binders/IL7Ra_denovo_minibinder_designs_IPD_Cao2022.fasta` before submitting.** This file exists precisely so we can do this. Install `blastp`/`mmseqs2`/`diamond` (none was available in this environment) and reject anything with ≥80% coverage and >70% identity. Use a safety margin — reject at **>60%**, since we don't know the organizers' exact alignment parameters and the threshold is a cliff.
2. **Design longer than 70 aa.** This is the strongest lever, and it is nearly free. All known IL-7Rα minibinders are **36–67 aa (median 63)**. If our design is **≥80 aa**, an alignment to a 67-aa known binder can cover at most 67/80 = **84%** of our length, and to the 55-aa validated binder at most **69%** — which is **below the 80% coverage trigger, so the rule cannot fire at all** (under the "总长度 = our length" reading). At 85–100 aa we are structurally immune to the entire corpus. Combined with the fact that 100 aa is the legal maximum, **designing 80–100 aa minibinders for IL-7Rα is the single highest-value decision available**, and it costs nothing in ipSAE terms (larger interfaces generally score better).
3. **Ask the organizers to disambiguate "总长度"** (ours, theirs, or the alignment's) and to state whether the frozen minibinder library includes the Cao et al. supplementary design models or only the reported binders. This is a legitimate clarification question and the answer changes our strategy. Also ask whether they will report a novelty failure or silently drop the sequence.
4. **Diversify topology deliberately.** Submit a mix — not 20 variants of one HHH bundle. The corpus is overwhelmingly HHH/HHHH and HEEH; mixed α/β or 4-helix topologies at 90+ aa are far from anything published.
5. **Consider targeting a secondary epitope on one or two designs** as orthogonal insurance, accepting lower functional relevance.

### 5.5 Also saved for IL-7Rα (antibodies — low novelty risk, useful context)

- `data/known_binders/IL7Ra_therapeutic_antibodies.fasta` (6 records): **Lusvertikimab / OSE-127** (OSE Immunotherapeutics, IgG4κ, Phase III), **Bempikibart** (Q32 Bio, IgG1κ, Phase II), **Crebankitug** (Pfizer/Zura Bio, IgG1λ, Phase II) — full VH and VL for each.
- `data/known_binders/IL7Ra_PDB_antibody_Fabs.fasta` (4 records): murine Fabs **4A10** (PDB **6P50**) and **2B8** (PDB **6P67**), heavy and light chains, both from IL-7Rα co-crystals.

**Source for the therapeutics: Thera-SAbDab** — bulk file `https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/static/downloads/TheraSAbDab_SeqStruc_OnlineDownload.csv` (1,133 therapeutics; Thera-SAbDab curates sequences from WHO INN proposed/recommended lists). Cross-validation: Bempikibart's VH CDR1 `DHAMH` and VL CDR1 `RASQGISSALA` / CDR2 `DASSLES` match the CDRs quoted in the anti-IL-7R patent literature (US 11008395 / WO2020154293A1), so the table is keyed correctly.

**Do therapeutic IgGs threaten a minibinder's novelty? Essentially no — the parent's instinct is right.** A VH domain is ~120 aa; for a 55–100 aa minibinder the coverage condition *can* in principle be satisfied by a sub-region of the VH. But identity would have to exceed 70% between an immunoglobulin β-sandwich framework and a de novo helical bundle, and my Experiment C control shows that unrelated sequences of *similar* fold/composition top out at 43%. Antibody-vs-minibinder identity will sit in the 20–35% band. **Ignore them for the minibinder tracks** (they matter enormously for the five VHH tracks, which we are not entering).

### 5.6 The AsterFire 16 designs: **NOT RETRIEVED — I did not fabricate them**

The WeChat post (`instruction.md` 2026-11-11 08:42, 折木) claims 16 TL1A minibinder designs with `candidate_0002` ipSAE 0.637, `candidate_0015` ipSAE 0.628 (binder pLDDT 91.8), `candidate_0010` with a minimal, non-overlapping epitope, and "对 84 条专利来源的抗 TL1A VH 做新颖性比对，16 条全部 PASS".

**I could not retrieve the actual sequences.** `https://mp.weixin.qq.com/s?__biz=MzcwMzM5OTYyNA==&mid=2247483976&...` returns an anti-bot interstitial ("环境异常，完成验证后即可继续访问") to non-browser clients. **The 16 sequences are UNCONFIRMED and I have not reproduced, guessed, or reconstructed any of them.** If they matter, the article needs to be opened in a real browser session (the `claude-in-chrome` tooling could do it).

Two things are still worth taking from the announcement, treating the numbers as the vendor's own claims:
- **An ipSAE of ~0.63 is being publicly marketed as a "credible interface" for TL1A.** Expect a large mass of submissions clustered at ipSAE 0.6–0.65, because many teams will run exactly this advertised one-click workflow. To clear the gate comfortably we should target **ipSAE ≥ 0.75** and treat 0.65 as the median of the pack, not as good.
- **Their novelty check was against 84 patent-derived anti-TL1A *VH* sequences** — i.e. antibodies, which per §5.5 are the wrong comparator for a minibinder and are nearly impossible to fail. **Their "16/16 PASS" is close to a tautology and provides no evidence of real novelty.** It does not test what the competition's minibinder rule tests.

### 5.7 TL1A: **no published de novo minibinder exists. This is the cleaner track.**

Stated plainly, as requested:

- **The PDB contains zero TL1A–binder complexes other than the natural decoy receptor.** The exhaustive RCSB search over O95150 returned 7 entries: **2O0O, 2QE3, 2RE9, 2RJK, 2RJL** are all apo TL1A (trimer or monomer, WT or Cys→Ser mutants), and **3K51, 3MI8** are TL1A–DcR3 complexes. **No antibody, no nanobody, no minibinder, no de novo design.**
- **Cao et al. 2022 did not include TL1A** among its 12 targets, so there is no IPD design corpus for it. A web search for published de novo TL1A/TNFSF15 minibinders (RFdiffusion/BindCraft, 2025–2026) returned nothing on target — only TL1A target-biology reviews and the therapeutic antibody pipeline. **I found no published de novo minibinder for TL1A. (Marked UNCONFIRMED in the strict sense that absence of evidence from a web search is not proof, but the PDB result is a strong independent negative.)**
- **I found no published anti-TL1A nanobody/VHH sequence either** — searches returned only anti-TNF-α nanobodies (a different target) and general nanobody reviews.
- What exists for TL1A, and is saved:
  - `data/known_binders/TL1A_therapeutic_antibodies.fasta` (6 records): **Tulisokibart / PRA023 / MK-7240** (Prometheus Biosciences→Merck, IgG1κ, Phase III), **Duvakitug / TEV-48574** (Sanofi/Teva, IgG1λ, Phase II), **Afimkibart** (Pfizer/Roche, IgG1κ, Phase III) — full VH and VL for each, from Thera-SAbDab. All three requested therapeutics were obtained; Afimkibart is a bonus third.
  - `data/known_binders/TL1A_natural_binder_DcR3.fasta` (2 records): DcR3 / TNFRSF6B (**O95407**) full-length plus the approximate ectodomain. **This is a natural 176-residue TNFR-family receptor, not a minibinder** — it is a legitimate reference for the functional epitope (see the 3K51 interface in §3) but it is implausible as a novelty comparator for a de novo design.

**Strategic consequence.** For TL1A the plausible frozen minibinder reference library is **empty or nearly so**, so the ≤70% identity rule is **effectively inert** — the only realistic novelty exposure is colliding with *other teams'* submissions under the separate duplicate rule (>90% coverage **and** 100% identity, first-come-first-served), which requires exact identity and is therefore only a risk if many teams converge on the identical AsterFire-advertised output. **TL1A is the lower-risk, cleaner track and should get priority**, with the caveat that it is also the track the sponsors have most heavily advertised (so expect more competitors and more duplicate-rule collisions). IL-7Rα carries a real novelty tax that must be actively managed per §5.4.

---

## 6. T6 — Oversubscription of the 1000-slot gate

### The arithmetic checks out

Rules §一(三).6: "赛事计划对不少于 **7000** 条候选序列开展湿实验验证：亲和力检测（SPR/BLI）每赛道约 **1000** 条；成药性与细胞活性检测每赛道取亲和力排名靠前的 **50** 个蛋白." → 7 × 1000 = 7000. Confirmed. (Note `/api/public/meta` still advertises the old "5000+ 条" — another stale field.)

**The seed teams eat into the 1000.** §五: "每赛道的湿实验入选序列，由种子团队序列与公开报名团队序列两部分构成". §五(一): 20 seed teams total, their compliant sequences are **all** admitted, exempt from ipSAE ranking. So **public slots per track = 1000 − (seed sequences in that track)**.

### Model and scenarios

Assumptions, stated: 1,200 registered teams (from the 2026-10-23 announcement "1200+队伍集结"; registration runs to 11-30 so this is a floor); teams enter ~2 tracks on average, so a uniform track would draw 2/7 ≈ 28.6% of submitting teams; a **minibinder bias multiplier** because minibinder design is far more mature than de novo VHH design (RFdiffusion/BindCraft vs the much harder nanobody problem) and because the sponsors published TL1A/minibinder tutorials specifically; serious teams max out near the 20-sequence cap; ~10–16 of the 20 seed teams enter a given minibinder track at 20 sequences each.

| scenario | submit rate | mb bias | teams in track | avg seq/team | public seqs | seed seqs | public slots | accept rate | required ipSAE percentile |
|---|---|---|---|---|---|---|---|---|---|
| optimistic | 25% | 1.0× | 86 | 15 | 1,287 | 240 | 760 | **59%** | ~p41 |
| **base** | 35% | 1.6× | 192 | 15 | 2,883 | 200 | 800 | **28%** | **~p72** |
| pessimistic | 50% | 2.0× | 343 | 18 | 6,178 | 280 | 720 | **12%** | ~p88 |
| worst | 65% | 2.2× | 491 | 18 | 8,834 | 320 | 680 | **8%** | ~p92 |

Sensitivity at a fixed 800 public slots: 2,000 public sequences → 40% accept (need ≥p60); 3,000 → 27% (≥p73); 4,000 → 20% (≥p80); 5,000 → 16% (≥p84); 6,000 → 13% (≥p87).

### Read

- **Central estimate: the gate is ~3–4× oversubscribed, so we need each sequence to sit around the 75th percentile of ipSAE or better, and realistically we should target p85–p90 to be robust to the pessimistic branch.** The spread between scenarios is dominated by one unknown — what fraction of the 1,200 registrants actually submit — so plan for the pessimistic branch.
- **Scoring is per sequence, not per team** (§五(二): "按单条序列统一打分排序"). Two consequences: (a) getting *some* sequences through is much easier than getting all 20 through, so there is no reason to submit fewer than 20; (b) because only each team's single best KD counts for the prize (§七(一).2), the optimal strategy is **20 sequences chosen for diversity of epitope and topology at high ipSAE**, not 20 near-copies of one top-scoring design — we want many independent shots at a low KD, not a hedge on one.
- **The screen is exactly reproducible, which is a large edge.** The config is fully specified — ESMFold2, ipSAE, 1 seed, 1 sample, no template, MSA for the target sequence only, no MSA for the designed protein. There is no seed noise and no sampling variance to beat. We can compute the organizers' exact score locally and submit only our top 20 by that metric, provided we pin down which model "ESMFold2" is and which ipSAE implementation they use (both still open in the fact base, §8).
- Calibration anchor: the only public ipSAE numbers for this competition are AsterFire's TL1A claims of **0.637 / 0.628** marketed as credible interfaces (§5.6). Treat ~0.63 as the expected centre of mass of the submission pool and **aim for ≥0.75**.

---

## 7. Residual open questions

| Question | Status |
|---|---|
| Will the organizers patch `/api/public/meta` before 2026-11-16 00:00 CST? | **Unknown — assume NO. Confirm by 2026-11-15.** |
| Does the frozen minibinder reference library include the 122k Cao et al. design models, or only the 2 reported binders? | **UNCONFIRMED — drives §5.4. Ask the organizers.** |
| Does "比对区域覆盖总长度的至少80%" mean our length, theirs, or the alignment's? | **Ambiguous in the rules. Ask.** Design ≥80 aa to be safe under any reading. |
| The 16 AsterFire TL1A sequences | **NOT RETRIEVED** (WeChat anti-bot). Needs a real browser session. Not fabricated. |
| Exact D1 N/C boundary of IL-7Rα | **UNCONFIRMED** (UniProt annotates only D2 = FnIII 131–231). D1 ≈ 21–130 by inference. |
| Which model is "ESMFold2"; which ipSAE implementation/version | Out of scope for VERIFY — still open per fact base §8. **Blocking for §6 self-selection.** |
| Whether a novelty failure is reported to the team or silently dropped | Unknown. Ask. |
