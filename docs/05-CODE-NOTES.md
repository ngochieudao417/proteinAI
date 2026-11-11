# 05-CODE-NOTES — API surface đã verify bằng source thật

> Code agent, 2026-11-11. Mọi dòng trong file này đọc từ **source thật**, không phải từ research doc:
> - `esm==3.4.3` wheel (`pip download esm==3.4.3` → giải nén → đọc `esm/models/esmfold2/*.py`, `esm/utils/msa/msa.py`, `esm/utils/structure/{input_builder,molecular_complex}.py`, `esm/models/hub.py`)
> - `https://huggingface.co/biohub/ESMFold2/raw/main/config.json` (commit `69869f737beffec5294845ede23db5fc0b4f509e`)
> - `https://huggingface.co/api/models/biohub/ESMFold2?blobs=true` (danh sách file + byte size)
> - `https://raw.githubusercontent.com/DunbrackLab/IPSAE/main/ipsae.py` (sha256 ghi trong report)
> - `https://raw.githubusercontent.com/Biohub/esm/main/cookbook/tutorials/esmfold2_local_gpu.ipynb`
> - `https://raw.githubusercontent.com/sokrypton/ColabFold/main/colabfold/colabfold.py` (protocol MSA API)
>
> Chỗ nào chưa chạy được trên GPU thật thì ghi rõ **CHƯA ĐO**.

---

## 0. CẬP NHẬT sau khi notebook 01 chạy thật trên L4 (2026-11-11)

### 0.1 Đã ĐO ĐƯỢC (số thật, giữ lại)

| Hạng mục | Giá trị |
|---|---|
| `esm==3.4.3` | cài sạch, torch **không đổi** (2.14.1+cu130) |
| Weights | tải 24.901 GB trong 230.8 s · stage local 127.9 s · 9 file verify byte-exact · cache tại `MyDrive/proteinAI/weights/ESMFold2` → **KHÔNG BAO GIỜ TẢI LẠI** |
| Load `dtype=bfloat16, esmc_precision="fp8"` | **105.7 s**, `vram_weights_gb = 6.413` |
| `n_params` | 266 656 029 · `n_buffer_elements` 6 346 017 830 (ESMC fp8 nằm trong buffer) |
| VRAM trống sau load | 15.396 GB · `ram_rss_gb` 9.37 |
| `kernel_backend_effective` | **`"fused"`** — Triton path CHẠY trên torch/CUDA này |
| `msa_encoder_present` | **true** · `offload_esmc_flag_present` true · `ccd_load_seconds` 8.3 |
| Thứ tự token | `is_permutation: true, is_identity: true` cả 5 lần → **rủi ro permutation ĐÃ ĐÓNG** |
| Chi phí | 14.63 phút ≈ 1.17 CU (còn ~88.8 CU) |

→ **§5 của file này (bảng VRAM maths) đã lạc hậu ở một chỗ quan trọng**: ước "bf16+fp8 ≈
6.8–7.0 GB" là đúng (đo 6.413), nhưng con số **26.3 GB fp32** gây hiểu nhầm — `esmc_precision`
quantize ESMC *bất kể* `dtype`, nên load fp32 cũng chỉ ≈ **6.95 GB**, không phải 26.3 GB.
Xem §0.2.

### 0.2 LỖI CHẶN & NGUYÊN NHÂN GỐC — `dtype=torch.bfloat16` là SAI

**Cả 5 lần fold đều chết cùng một lỗi**, không phụ thuộc MSA, không phụ thuộc kích thước:

```
RuntimeError: mat1 and mat2 must have the same dtype, but got Float and BFloat16
```

Đọc `esm/models/esmfold2/model.py` → `EsmFold2Model.forward`:

| dòng | code | |
|---|---|---|
| 973 | `with torch.amp.autocast("cuda", enabled=use_amp, dtype=torch.bfloat16):` | **mở** autocast |
| 974–1157 | inputs_embedder → pair init → lm_encoder → msa_encoder → recycle loop → parcae | trunk, trong autocast |
| 1159 | `z = z.float()` | pair rep **upcast fp32** |
| **1160** | *hết block* | **autocast ĐÓNG** |
| 1181 | `self.distogram_head(z + z.transpose(-2, -3))` | fp32 activation, **ngoài** autocast |
| 1185 | `self.structure_head.sample(z_trunk=z, s_inputs=x_inputs, ...)` | diffusion decoder |
| 1239–1241 | `self.confidence_head(..., z=z.detach().float(), ...)` | head sinh PAE/pLDDT |

**Ba head sau trunk được thiết kế để chạy fp32**, và source nói rõ ở nhiều chỗ:

- `DiffusionConditioning.forward` (`layers.py:1520-1534`) cưỡng chế
  `z_trunk.to(torch.float32)`, `relative_position_encoding.to(torch.float32)`,
  `s_inputs_eff.to(torch.float32)`, `t` fp32.
- `ConfidenceHead.forward` (`model.py:221`): `pair.add_(pair_delta.float())` rồi
  `self.pae_head(self.pae_ln(pair))` trên fp32.
- `InputsEmbedder.forward` (`layers.py:2136`): `torch.cat([a, aatype, profile, deletion_mean])`
  trong đó `a` bf16 nhưng `aatype`/`profile`/`deletion_mean` bị `.float()` ở
  `model.py:936/976/977` → **type promotion làm `x_inputs` = fp32**. Và `x_inputs` chính là
  `s_inputs` của structure_head (1187) + confidence_head (1240).

⇒ **Cả hai đầu vào của ba head (`z` và `x_inputs`) đều fp32 một cách có chủ ý.**

Nhưng `from_pretrained(..., dtype=torch.bfloat16)` gọi `model.to(torch.bfloat16)`
(`hub.py:194`) — cast **toàn bộ** model, kể cả ba head đó. Linear đầu tiên gặp phải là
`self.distogram_head` (`nn.Linear(256, 64)` trần) ở dòng 1181 → **fp32 × bf16 = đúng lỗi trên**.

**Ba hệ quả:**

1. **Diffusion decoder CHƯA TỪNG CHẠY.** Lỗi ở `distogram_head`, ngay sau trunk. Mọi số
   VRAM/thời gian của notebook 01 là **sàn dưới của sàn dưới** (7.611 GB @ 274 aa,
   activation 1.198 GB — thiếu cả diffusion lẫn confidence head).
2. **`esmc_precision` KHÔNG liên quan.** `esmc_precision="bf16"` (± `_offload_esmc`) sẽ chết
   **y nguyên**. Đừng tiêu compute thử.
3. **Cách sửa: `dtype=None` (fp32).**

### 0.3 Cấu hình ĐÚNG

```python
model = EsmFold2Model.from_pretrained(
    MODEL_DIR, device="cuda",
    dtype=None,                  # <-- fp32. KHÔNG truyền torch.bfloat16.
    esmc_precision="fp8", load_esmc=True,
).eval()
model.set_kernel_backend("fused")
```

`hub.py:194` → `if dtype is not None: model.to(dtype)` ⇒ `dtype=None` giữ nguyên fp32.
Rồi `set_esmc_precision("fp8")` vẫn cast ESMC bf16→fp8, nên ESMC **không** tốn 25 GB.

**Trunk vẫn chạy bf16** — không mất gì về tốc độ hay activation:

- autocast `model.py:973` bọc toàn bộ trunk → mọi `nn.Linear` ra bf16;
- `FoldingTrunk.forward` (`layers.py:2834-2840`) tự cast `pair` sang bf16 khi bật fused backend;
- fused kernel cast weight **tại call site** (`_to_fused_kernel_dtype`, `layers.py:2625-2630`),
  nên weight fp32 không phá Triton path. `fused_lnlin_swiglu` còn có
  `@torch.amp.custom_fwd(cast_inputs=torch.bfloat16)` (`kernels/fused_lnlin_swiglu.py:356`).

| Cấu hình | ESMC | Trunk+head | Tổng weights | Kết quả |
|---|---|---|---|---|
| `dtype=bfloat16, fp8` | 5.9 GB | 0.53 GB | 6.413 GB **(đo)** | ❌ chết ở `distogram_head` |
| **`dtype=None, fp8`** | 5.9 GB | 1.07 GB | **≈ 6.95 GB** | ✅ đường tác giả thiết kế |
| `dtype=bfloat16, fp8` + cast 3 head về fp32 | 5.9 GB | mixed | ≈ 6.5 GB | ✅ (fallback, nhẹ hơn 0.45 GB) |
| `dtype=None, bf16` | 12.7 GB | 1.07 GB | ≈ 13.8 GB | ✅ nhưng tốn gấp đôi, vô ích |

**Giá của fp32 chỉ là +0.53 GB weights.** Activation gần như không đổi (trunk vẫn bf16).

### 0.4 MSA — ĐÃ CÓ a3m THẬT, không còn phụ thuộc mạng

Notebook 01 rơi về `query_only_FALLBACK` (depth 1, `fetch_info` rỗng). **Đã điều tra: API
hoàn toàn bình thường** — gọi lại cùng protocol (`POST /ticket/msa` với
`q=">101\n<seq>\n"`, `mode=env` → poll `/ticket/<id>` → `GET /result/download/<id>` →
`uniref.a3m`) trả về ngay:

| target | depth | match columns | query == target | thời gian | file |
|---|---|---|---|---|---|
| IL-7Rα `P16871` E21-D239 (219 aa) | **921** | 219 ✓ | ✓ | 7 s | `data/msa/IL7Ra_P16871_E21-D239.uniref.a3m` |
| TL1A `O95150` L72-L251 (180 aa) | **2037** | 180 ✓ | ✓ | 29 s | `data/msa/TL1A_O95150_L72-L251.uniref.a3m` |

sha256: IL-7Rα `08e2a05b5b4787bb28c03302b1d56c3856980423a5582f2096bc6f5238b0be03` ·
TL1A `4464ad15470e2ecda72ea879beb5972059062085b04ff515cef26f642e092c27`

Cả hai **đã nhúng gzip+base64 vào notebook 02** (kèm verify sha256) → đường quan trọng
nhất không còn cần mạng. Thứ tự ưu tiên trong notebook: cache Drive (verify sha) →
payload nhúng → API (ghi **nguyên** request/response) → query-only (ghi rõ
`msa_depth_realistic: false`).

> Lưu ý: `mode=env` trả 3 file; ta dùng `uniref.a3m`. Depth 921/2037 đã quá
> `msa_max_depth=1024` nên `fold()` còn subsample lại mỗi loop.

Một chi tiết source quan trọng cho TL1A trimer — `prepare_input.py:1182-1197`:
> *"A single ProteinInput with id=['A','B','C',...] yields one item but multiple chains
> (one per id); broadcast the MSA across all of them."*

⇒ `ProteinInput(id=["A","B","C"], sequence=TL1A, msa=MSA_TL1A)` tự broadcast MSA cho cả 3
protomer. Không cần làm gì thêm. (Nhưng pairing theo taxonomy vẫn là -1 vì a3m ColabFold
không có token `key=` — xem §4.)

### 0.5 Notebook 02 đo gì (`notebooks/02_dtype_fix_and_measure.ipynb`)

| Cell | Việc |
|---|---|
| 1–3 | tiền kiểm · cài esm · **verify cache, tải 0 byte** (dừng cứng nếu thiếu) |
| 4 | sequence (thêm **DcR3** O95407 res30-205, 176 aa, receptor TL1A thật, PDB 3K51 — PROBE) |
| 5 | MSA thật nhúng sẵn (§0.4) |
| 6 | load fp32 + in dtype từng khối + **micro-repro** (clone `distogram_head` sang bf16, bắn fp32 vào → tái hiện lỗi trong vài ms, 0 giây GPU) |
| 7 | adapter ipSAE (giữ nguyên đường AF2 đã test) |
| 8 | `fold_once`: đo VRAM/thời gian **thật**, lưu `pae.npy`, bắt OOM riêng, fail-soft |
| 9 | **chốt dtype**: fold tối thiểu 80 aa / 1 loop / 4 step qua ladder 4 cấu hình, rẻ trước; chỉ reload khi buộc phải |
| 10 | **đối chứng dương/âm (quan trọng nhất)** + **tính tiền định** (fold y hệt 2 lần seed=0, so PAE từng bit) |
| 11 | ipSAE có MSA vs không MSA |
| 12 | trần residue 280→420→540→640→775, OOM thì thử lại cùng mốc với `_offload_esmc=True`+`chunk_size=32` |
| 13 | **monomer vs trimer** (cùng DcR3, tham số y hệt, cả hai bật offload) |
| 14 | REPORT → `out/02_report.json` |

Trần cứng 75 phút (≈ 6 CU). Mỗi nhóm có `have_budget()` gác, bỏ qua sạch nếu hết giờ.

### 0.6 Vì sao chọn DcR3 cho câu hỏi monomer-vs-trimer

Chưa có binder TL1A do ta design. Dùng binder bất kỳ sẽ cho ipSAE ≈ 0 ở **cả hai** cấu
hình → delta ≈ 0 → kết luận "hai cấu hình như nhau" là **giả**. DcR3 là receptor tự nhiên
đã có tinh thể phức (3K51) và — như cả họ receptor TNF — bám vào **rãnh giữa hai protomer
kề nhau**. Nếu monomer khác trimer thì nó lộ ra ở đây. Notebook vì vậy đọc **cả giá trị
tuyệt đối**, và nếu cả hai < 0.1 thì tự kết luận `KHONG_KET_LUAN_CA_HAI_GAN_0` thay vì
báo "giống nhau". DcR3 176 aa, 20+ cysteine, không khai báo disulfide → đây là điểm yếu
đã biết của probe này.

### 0.7 Hai discrepancy vs `config.json` (đã ghi vào report máy đọc được)

| Tham số | `builder.fold()` default | `config.json` | Ghi chú |
|---|---|---|---|
| `num_sampling_steps` | **200** | `structure_head.inference_num_steps` = **14** | docstring của `forward` nói default là `config.structure_head.inference_num_steps`, nhưng `fold()` truyền 200 **tường minh** → `fold()` thắng |
| `lm_dropout` | **0.3** | `lm_encoder.lm_dropout` = **0.25** | docstring `fold()`: *"paper folding-eval value; 0/None disables"*. **Ngẫu nhiên**, fresh mask mỗi loop |

`fold()` bọc forward trong `_seed_context(seed)` + `_lm_dropout_context(model, lm_dropout)`
(`processor.py:418-421`) ⇒ `seed=0` *nên* cho tiền định. **CHƯA ĐO** — cell 10 của notebook
02 trả lời bằng cách fold y hệt 2 lần rồi so PAE từng bit. Nếu KHÔNG tiền định thì câu
"1 seed / 1 sample" của cửa sàng là một lần **rút thăm**, và xếp hạng ổn định sẽ cần fold
nhiều lần → đắt gấp n lần.

### 0.8 CHƯA ĐO (notebook 02 sẽ trả lời)

- VRAM/thời gian **thật** khi diffusion decoder chạy (số cũ là sàn dưới)
- ipSAE binder 7S5B thật vs xáo trộn → thước đo có hoạt động không
- tính tiền định ở `seed=0`
- trần residue thật → TL1A trimer khả thi trên L4 không
- delta monomer vs trimer
- delta ipSAE có MSA vs không MSA

---

## 1. API surface chính xác (CONFIRMED từ source)

### Import

```python
from esm.models.esmfold2 import (
    EsmFold2Model, ESMFold2InputBuilder,
    ProteinInput, StructurePredictionInput,
    ESMFOLD2_HF_REPO,            # == "biohub/ESMFold2"
)
from esm.utils.msa import MSA                       # cũng re-export tại esm.models.esmfold2
from esm.utils.msa.msa import is_a3m_insertion
```

`esm/models/esmfold2/__init__.py` export đủ: `EsmFold2Model`, `EsmFold2ExperimentalModel`,
`ESMFold2InputBuilder`, `ProteinInput`, `RNAInput`, `DNAInput`, `LigandInput`, `MSA`,
`Modification`, `CovalentBond`, `DistogramConditioning`, `PocketConditioning`,
`StructurePredictionInput`, `ChainInfo`, `prepare_esmfold2_input`, `clean_esmfold2_input`,
`MolecularComplex`, `MolecularComplexResult`, `ESMFOLD2_HF_REPO`,
`ESMFOLD2_EXPERIMENTAL_HF_REPO`, `EsmFold2Config`, `load_ccd`.

### Load — `from_pretrained` CÓ tham số `dtype`

```python
# esm/models/esmfold2/model.py:589
@classmethod
def from_pretrained(cls,
    pretrained_model_name_or_path,
    *,
    load_esmc: bool = True,
    esmc_precision: str = "fp8",       # "fp8" | "bf16" | "fp32"
    config: EsmFold2Config | None = None,
    device: torch.device | str = "cpu",
    dtype: torch.dtype | None = None,
    **kwargs,                          # -> resolve_model_dir: revision, cache_dir, token,
)                                      #    local_files_only, force_download
```

Thứ tự thực tế bên trong (quan trọng — đây là lý do bf16 cứu được ta trên L4):

1. `resolve_model_dir()` — **local dir có `config.json` thì trả về nguyên path, không gọi mạng**.
2. `_load_pretrained(local_dir, config, device="cpu", dtype=dtype)` → `init_empty_weights()` (accelerate),
   đọc safetensors, `load_state_dict(strict=False, assign=True)`, **`del raw, adapted` TRƯỚC khi
   `model.to(dtype)`**, rồi `model.to("cpu")`, `.eval()`.
   → Cast sang bf16 xảy ra **trên CPU**, fp32 26.3 GB **không bao giờ chạm GPU**.
3. `set_esmc_precision(esmc_precision)` — vẫn trên CPU. `"fp8"` → cast bf16 rồi `convert_esmc_to_fp8_()`
   thay mọi projection trong `esmc.transformer.blocks` bằng `Fp8Linear` (weight lưu dạng `uint8` view
   nên `.to(dtype)` về sau **không** cast lại được).
4. `return model.to(device)`.

Load đúng cho L4:

```python
# !!! LẠC HẬU — dtype=torch.bfloat16 LÀM CHẾT fold. Xem §0.2 / §0.3.
# Cấu hình đúng:
model = EsmFold2Model.from_pretrained(
    MODEL_DIR, device="cuda", dtype=None, esmc_precision="fp8", load_esmc=True,
).eval()
model.set_kernel_backend("fused")     # None | "fused" | "fast"
```

### Memory / perf knob (research KHÔNG biết hai cái này — đều quan trọng)

| Knob | Source | Tác dụng |
|---|---|---|
| `model.set_chunk_size(int \| None)` | `model.py:700`, default `_DEFAULT_CHUNK_SIZE = 64` | chặn peak của các op L² (triangle / OPM / pair-transition). `None` = nhanh hơn ở L<600 nhưng dễ OOM. Hạ xuống 32 = đòn bẩy VRAM. |
| `model._offload_esmc = True` | `model.py:490` (khai báo), `:1063` (dùng trong **forward thường**, không riêng đường distributed), `:721` (restore) | sau khi tính LM hidden states một lần, đẩy ESMC (**≈6.4 GB fp8 / ≈12.7 GB bf16**) sang CPU + `empty_cache()`, giải phóng cho trunk + diffusion. Comment trong source: "Validated with CP: 2.14x end-to-end + 18.5 GB lower peak vs fp32". |

→ **Research nói PR #385 (`low_memory_mode`) còn open ⇒ "không có CPU offload". Sai.** Cơ chế offload
ESMC **đã có trong `esm` 3.4.3** dưới dạng `model._offload_esmc`, và nó chạy trong forward path thường.
Đây là đòn bẩy chính để đẩy trần residue lên.

### Builder + fold

```python
builder = ESMFold2InputBuilder(ccd_cache: Path | None = None)
```

⚠️ `__init__` gọi `load_ccd(ccd_cache)`:
- `ccd_cache=None` → `hf_hub_download(repo_id="biohub/ESMFold2", filename="ccd.pkl")` vào HF cache mặc định
- `ccd_cache=<dir>` → đọc **`<dir>/ccd.pkl`**, và **KHÔNG tải nếu thiếu** (raise `FileNotFoundError`)
- env `ESMCFOLD_CCD_PATH` đọc ở **module import time** (`conformers.py:18`), nên muốn dùng thì phải set
  **trước** `import esm...`

`ccd.pkl` = **417 MB** và **KHÔNG nằm trong `snapshot_download`** của `resolve_model_dir`
(`allow_patterns=["*.json", "*.safetensors"]`) → phải tải tay nếu muốn cache offline.

```python
# esm/models/esmfold2/processor.py:331 — mọi tham số sau `input` là keyword-only
builder.fold(model, input, *,
    num_loops: int = 20,
    num_sampling_steps: int = 200,
    num_diffusion_samples: int = 1,
    seed: int | None = None,
    noise_scale=None, step_scale=None, max_inference_sigma=None,
    lm_mask_pct=None,
    lm_dropout: float | None = 0.3,          # <-- NGẪU NHIÊN, mặc định 0.3
    early_exit=None,                         # deprecated, bị warn + bỏ qua
    msa_max_depth: int | None = 1024,
    msa_column_mask_rate: float = 0.1,
    msa_subsample_at_inference=None,         # deprecated
    include_embeddings: bool = False,
    complex_id: str = "pred",
) -> MolecularComplexResult | list[MolecularComplexResult]
```

`num_diffusion_samples == 1` → trả về **một** `MolecularComplexResult`, không phải list.

Hai nửa của `fold()` (dùng khi cần tensor thô):
`builder.prepare_input(input, seed=None, device=None) -> (features: dict, chain_infos: list[ChainInfo])`
và `builder.decode(output, features, chain_infos, *, num_diffusion_samples=1, complex_id="pred")`.

### Input dataclass

```python
# esm/utils/structure/input_builder.py
@dataclass
class ProteinInput:
    id: str | list[str]                 # list[str] = "một sequence, nhiều copy"
    sequence: str
    modifications: list[Modification] | None = None
    msa: MSA | None = None

@dataclass
class StructurePredictionInput:
    sequences: Sequence[ProteinInput | RNAInput | DNAInput | LigandInput]
    pocket: PocketConditioning | None = None        # NotImplementedError nếu set
    distogram_conditioning: list[DistogramConditioning] | None = None
    covalent_bonds: list[CovalentBond] | None = None
```

`MSA.from_a3m(path, remove_insertions: bool = False, max_sequences: int | None = None)` — cookbook
gọi bằng keyword: `MSA.from_a3m(path="chainA.a3m")`. `parse_fasta` **bỏ qua dòng bắt đầu bằng `#`**
→ a3m của ColabFold (có dòng `#219\t1` ở đầu) dùng được trực tiếp.

### Output

```python
# esm/utils/structure/molecular_complex.py:34
@dataclass
class MolecularComplexResult:
    complex: MolecularComplex
    plddt: torch.Tensor | None        # [L], THANG 0–1 (_categorical_mean(start=0.0, end=1.0))
    ptm: float | None
    iptm: float | None
    pae: torch.Tensor | None          # [L, L], đơn vị Å, thang 0–32
    pde: torch.Tensor | None
    distogram: torch.Tensor | None
    pair_chains_iptm: torch.Tensor | None     # [n_chains, n_chains]
    residue_index / entity_id / output_embedding_* ...
```

`result.complex.to_mmcif(*, name_null_chains: bool = False) -> str` (qua biotite).
`result.complex` fields: `sequence` (3-letter, per token), `chain_id` (int asym_id per token),
`token_to_atoms` ([L,2] start/end), `atom_positions`, `atom_names`, `atom_elements`, `plddt`,
`metadata.chain_lookup` (`{asym_id: chain_id_str}`).

---

## 2. Thứ tự token vs thứ tự residue — cái "gotcha" thật ra thế nào

Research nói: *"PAE ở token order của `prepare_input()`, KHÔNG phải mmCIF residue order"*.
Đọc source thì chính xác hơn là:

- `build_chains_from_input()` (`prepare_input.py:676`) đi theo đúng thứ tự `input.sequences`, gán
  `token_index` **tăng dần liên tục** và `asym_id` tăng dần.
- `build_molecular_complex_from_features()` (`output.py:18`) lặp `for ci in chain_infos` rồi
  `groupby(ci.tokens, key=residue_index)` — **cùng thứ tự đó**.

⇒ Với input **chỉ protein, không `Modification`, không ligand**, **token order == residue order ==
identity**. Permutation chỉ khác identity khi có residue bị atom-tokenized (PTM) hoặc ligand.

**Nhưng không được giả định.** Notebook làm thế này:

1. Gọi `builder.prepare_input(spi, seed=0)` riêng (CPU, rẻ) để lấy `chain_infos`, dựng
   `perm[k] = chain_infos[*].tokens[*].token_index` theo đúng thứ tự xuất file, rồi assert
   `perm` là permutation hợp lệ và ghi cờ `is_identity`.
2. **Tự ghi PDB theo thứ tự token** thay vì dùng `to_mmcif()`. Nhờ vậy permutation trở thành identity
   *by construction*, và nếu `perm` không identity thì PAE được reindex `pae[np.ix_(perm, perm)]`.

---

## 3. ipSAE v1 — đường AF2 là đường đúng cho ta (research chỉ nói đường AF3)

Research đề nghị viết JSON kiểu AF3:
`{"pae":..., "atom_plddts":..., "token_chain_ids":..., "token_res_ids":...}` + `.cif`.
Đọc `ipsae.py` thì:

- **`token_chain_ids` / `token_res_ids` KHÔNG được dùng.** ipsae lấy chain/resnum từ file structure.
- Đường **AF3** (`.cif` + `.json`) dựng `token_mask` bằng cách **quét từng dòng ATOM/HETATM của cif**,
  và index `atom_plddts` theo **`atom_num - 1`** (serial của ATOM). Hai chỗ này phụ thuộc vào cột
  `_atom_site.*` chính xác mà biotite sinh ra → thêm hai chỗ sai âm thầm.
- Đường **AF2** (`.pdb` + `.json`) đơn giản hơn nhiều (`ipsae.py:428`):
  ```json
  {"pae": [[...]],       // [numres, numres], theo đúng thứ tự residue trong PDB, KHÔNG subset
   "plddt": [...],       // per-RESIDUE, thang 0–100
   "iptm": 0.42, "ptm": 0.55}
  ```
  Không có bước `token_mask`. `pae_matrix = np.array(data['pae'])` dùng nguyên.

→ **Dùng đường AF2.** Ta tự ghi PDB cột-cố-định nên kiểm soát được 100%.
`parse_pdb_atom_line` đọc: `[6:11]` serial · `[12:16]` atom name · `[17:20]` resname · `[21]` chain ·
`[22:26]` resnum · `[30:38][38:46][46:54]` x/y/z.

**Bẫy shape phải assert trước khi gọi:** ipsae dựng `residues` từ atom `CA`, và `cb_residues` từ
`CB` (hoặc `CA` nếu `GLY`). Nếu hai danh sách lệch độ dài thì `np.ix_` broadcast sai **mà không báo lỗi**.
Notebook assert `n_CA == n_CB == n_residues == pae.shape[0]`.

**plddt thang:** ESMFold2 trả 0–1, ipsae/pDockQ cần 0–100 → `× 100`. Chỉ ảnh hưởng pDockQ, không
ảnh hưởng ipSAE, nhưng vẫn làm cho đúng.

**CLI + output:**
```bash
python ipsae.py <pae.json> <model.pdb> 10 10
# -> model_10_10.txt, model_10_10_byres.txt, model_10_10.pml
```
Tên file: `pae_string = str(int(cutoff))`, pad `"0"` nếu `<10` → cutoff 10/10 cho hậu tố `_10_10`.
Bảng (whitespace-split): `[0]=Chn1 [1]=Chn2 [2]=PAE [3]=Dist [4]=Type [5]=ipSAE(d0res)`
`[6]=ipSAE_d0chn [7]=ipSAE_d0dom [8]=ipTM_af [9]=ipTM_d0chn [10]=pDockQ [11]=pDockQ2 [12]=LIS` …
`Type` ∈ {`asym`, `max`}. **Research đúng chỗ này** (field 4 == asym, field 5 == ipSAE).

**ĐÃ TEST OFFLINE:** adapter (`write_pdb_token_order` + `run_ipsae`) chạy với một
`MolecularComplex` mock dựng theo đúng field layout thật → `ipsae.py` rc=0, parse ra
`{'A->B': 0.296, 'B->A': 0.314}` + dòng `max`. Nên phần adapter là phần ít rủi ro nhất của notebook.

---

## 4. MSA — cơ chế thật

`biohub/ESMFold2` `config.json`: **`msa_encoder.enabled = true`**, `overwrite=true`,
`hidden_size=128`, `num_hidden_layers=4`, `max_depth=1024`, `column_mask_rate=0.1`.
→ MSA path là thật trên checkpoint này.

Nhưng: `EsmFold2Model.__init__` chỉ dựng `self.msa_encoder` khi `msa_cfg.enabled`
(`model.py:531`), và forward chỉ dùng MSA khi `self.msa_encoder is not None and msa is not None`
(`model.py:1090`) — **không có warning nào nếu encoder là None**. Đúng y bệnh của `-Fast`.
→ Notebook **bắt buộc** in `model.msa_encoder is not None`.

`msa=None` cho một chain: `compute_msa_features` (`prepare_input.py:1192`) tự thay bằng
`MSA.from_sequences([item.sequence])` (depth 1), và `build_chains_from_input` phát
`UserWarning: No MSA provided for B, using single sequence mode`. **Warning này là đúng** —
cửa sàng cố tình không cho binder MSA.

Pairing theo taxonomy: `paired_msa.py` đọc token `key=N` trong FASTA header. a3m của ColabFold
**không có** `key=`, nên mọi row là unpaired (-1). Với ta (chỉ 1 chain có MSA) thì không ảnh hưởng.
Nếu sau này fold trimer mà **cả 3 protomer** đều có MSA, pairing sẽ thành vấn đề thật.

Chiều dài: `L = số ký tự không-insertion của row 0` (`is_a3m_insertion` = lowercase hoặc `.`).
Phải khớp `len(sequence)`, nếu không `prepare_input` sẽ lỗi index → notebook assert trước.

**ColabFold MSA API** (đọc từ `sokrypton/ColabFold/colabfold/colabfold.py:68`):
```
POST https://api.colabfold.com/ticket/msa   form: q="> 101\n<seq>\n", mode="env"
  -> {"id": ..., "status": "PENDING|RUNNING|COMPLETE|ERROR|MAINTENANCE|UNKNOWN|RATELIMIT"}
GET  https://api.colabfold.com/ticket/<id>          -> {"status": ...}
GET  https://api.colabfold.com/result/download/<id> -> tar.gz: uniref.a3m,
                                                      bfd.mgnify30.metaeuk30.smag30.a3m, pdb70.m8
```
Nên set `User-Agent` (server yêu cầu, sắp thành hard error). Notebook dùng `uniref.a3m`,
có timeout cứng, fallback sang a3m query-only và **ghi rõ `msa_depth_realistic: false`**.

---

## 5. VRAM maths (tính tay — con số thật phải chờ notebook)

Từ `config.json`: ESMC backbone `hidden_size=2560`, `num_hidden_layers=80`,
`intermediate_size=6912`, `num_attention_heads` theo cấu hình → mỗi layer ≈
`4·2560² (attn) + 3·2560·6912 (SwiGLU) ≈ 26.2 M + 53.1 M = 79.3 M` → **ESMC ≈ 6.35 B**.
Tổng model 6.58 B ⇒ trunk + heads ≈ **0.23 B**.

| Cấu hình | ESMC | Trunk | Tổng weights trên GPU |
|---|---|---|---|
| fp32 (naive, như cookbook) | 25.4 GB | 0.9 GB | **26.3 GB → OOM ngay trên L4 21.5 GB** |
| bf16 toàn bộ | 12.7 GB | 0.46 GB | **13.2 GB** → còn ~8 GB cho activation |
| bf16 trunk + **ESMC fp8** | 6.35 GB | 0.46 GB | **≈ 6.8–7.0 GB** → còn ~14.5 GB cho activation |
| bf16 + fp8 + `_offload_esmc` | 0 GB *trong pha trunk/diffusion* | 0.46 GB | ≈ 0.5 GB sau pha LM |

Activation (bf16, pair rep `d_pair=256`, `pae/distogram` 64 bin):

| L | pair rep L²·256·2 B | pae_logits L²·64·4 B |
|---|---|---|
| 280 | 40 MB | 20 MB |
| 320 | 52 MB | 26 MB |
| 595 | 181 MB | 91 MB |
| 640 | 210 MB | 105 MB |
| 820 | 345 MB | 172 MB |

Mỗi tensor nhỏ, nhưng trunk có 48 layer × nhiều temp, `set_chunk_size(64)` chặn phần L² tệ nhất,
và atom-attention là sliding window 128. Phép tính thô này cho thấy **"24 GB ≈ 500 aa" của maintainer
rất có thể là nói về load fp32/bf16 thuần, và bf16+fp8+offload có thể vượt khá xa 500 aa.**
**CHƯA ĐO.** Đó chính là lý do có cell 9.

Lưu ý RAM hệ thống: `_load_pretrained` đọc toàn bộ safetensors (mmap) rồi `model.to(bfloat16)` trên
CPU. Peak RSS thực ước **13–20 GB** (page cache thu hồi được). 52 GB là đủ, nhưng nếu runtime bị kill
lúc load thì nguyên nhân là đây.

---

## 6. Những chỗ research sai hoặc thiếu

| # | Research nói | Thực tế (source) |
|---|---|---|
| 1 | *"không có CPU offload, PR #385 còn open"* | **`model._offload_esmc = True` đã có trong `esm` 3.4.3**, chạy trong forward thường (`model.py:1063`), giải phóng ≈6.4 GB fp8. Đòn bẩy VRAM lớn nhất ta có. |
| 2 | Không nhắc `set_chunk_size` | Có, default 64, hạ được xuống 32 để bớt VRAM (`model.py:700`). |
| 3 | ~~Snippet `from_pretrained("biohub/ESMFold2", device="cuda")` thiếu `dtype` → OOM 26.3 GB. Phải `dtype=torch.bfloat16`.~~ | **DÒNG NÀY CỦA TÔI SAI — đã sửa, xem §0.2/§0.3.** Cookbook **đúng**: `dtype=None` là cấu hình tác giả thiết kế. Vì `esmc_precision` default `"fp8"` quantize ESMC *bất kể* `dtype`, load fp32 chỉ ≈ **6.95 GB**, không bao giờ là 26.3 GB. Ngược lại `dtype=torch.bfloat16` **làm chết fold** ở `distogram_head` (fp32 activation × bf16 weight). |
| 4 | Không nhắc `esmc_precision` | `from_pretrained(..., esmc_precision="fp8" \| "bf16" \| "fp32")`, default **đã là `"fp8"`**. fp8 GEMM cần compute capability ≥ 8.9 → L4 (8.9) có. Đây mới là knob tiết kiệm VRAM — **không phải `dtype`**. |
| 5 | Adapter ipSAE cần `token_chain_ids` / `token_res_ids` trong JSON | **Không dùng tới.** Và đường **AF2** (`.pdb`) đơn giản + ít rủi ro hơn đường AF3 (`.cif`). |
| 6 | Không nhắc `ccd.pkl` | `ESMFold2InputBuilder()` tải riêng **417 MB** `ccd.pkl`, **không** nằm trong `snapshot_download`. Phải cache riêng. |
| 7 | Không nhắc `lm_dropout` | `builder.fold()` default `lm_dropout=0.3` — **có ngẫu nhiên**. Deterministic chỉ khi truyền `seed`. `config.json` lại ghi `lm_encoder.lm_dropout=0.25`. → **Phải hỏi organizer.** |
| 8 | `num_sampling_steps` "không có default trong config" | `config.json` → `structure_head.inference_num_steps = 14`, nhưng `builder.fold()` default **200**. Hai số rất khác nhau → **phải hỏi organizer**. |
| 9 | — | `clean_esmfold2_input` tách chuỗi có `\|`/`:` thành nhiều chain và tự slice MSA theo từng chain. Có ích khi fold trimer. |

Research **đúng** ở: identity của model, repo HF, `msa_encoder` tồn tại, `-Fast` bỏ qua MSA không báo,
plddt thang 0–1, công thức ipSAE, layout cột của bảng ipsae, `set_kernel_backend("fused")`,
tên class `EsmFold2Model` / `ESMFold2InputBuilder`, `MSA.from_a3m`.

---

## 7. Manifest weights (dùng để verify cache)

`https://huggingface.co/api/models/biohub/ESMFold2?blobs=true`, commit `69869f737beffec5294845ede23db5fc0b4f509e`:

| file | byte |
|---|---|
| `config.json` | 8 214 |
| `model.safetensors.index.json` | 245 594 |
| `model-00001-of-00006.safetensors` | 4 955 672 788 |
| `model-00002-of-00006.safetensors` | 4 971 194 264 |
| `model-00003-of-00006.safetensors` | 4 997 460 576 |
| `model-00004-of-00006.safetensors` | 4 979 069 192 |
| `model-00005-of-00006.safetensors` | 4 934 556 472 |
| `model-00006-of-00006.safetensors` | 1 481 396 616 |
| `ccd.pkl` | 417 306 584 |
| **tổng** | **26 737 068 175 B = 24.90 GiB = 26.74 GB** |

Notebook verify **từng byte**; lệch một byte là tải lại đúng file đó. Thêm check: 8 byte đầu mỗi
safetensors = độ dài header JSON (bắt file cắt cụt), và `index.json` phải chỉ trỏ tới shard có trong manifest.

---

## 8. Notebook tiếp theo phụ thuộc kết quả smoke test ra sao

| Kết quả cell 9 (trần residue) | Notebook 02 làm gì |
|---|---|
| **≥ 640 aa chạy được** | Chấm cả hai track tại chỗ trên L4, TL1A fold **trimer**. 02 = gate reproducer đầy đủ + batch scorer. Đây là nhánh tốt nhất. |
| **~280–595 aa** | IL-7Rα chấm tại chỗ. TL1A chấm **monomer** tại chỗ (hedge §6.4 của research: 12 design nhắm mặt một-protomer, 8 nhắm groove) và đẩy phần trimer sang platform sponsor / L40S thuê. 02 = scorer + script export cho sponsor. |
| **< 280 aa, hoặc OOM cả ở 235 aa** | L4 không chấm nổi. 02 chuyển sang Biohub Batch API (giới hạn 4096 token) hoặc MatwingsVenus. Phải báo ngay. |

| Kết quả control check | Phản ứng |
|---|---|
| `ipsae_min(7S5B binder) ≥ 0.5` **và** gap vs scrambled `≥ 0.15` | Thước đo chạy đúng → đi tiếp sang generator (BindCraft trên L4). |
| binder thật điểm thấp (`< 0.5`) | **Setup của ta sai**, không phải design sai. Nghi can theo thứ tự: (a) MSA depth không thật, (b) chiều ipSAE / cutoff, (c) `lm_dropout`/`num_sampling_steps` lệch, (d) thứ tự token. Notebook 02 phải là notebook *debug thước đo*, không phải notebook sinh design. |
| gap < 0.15 | Thước đo không phân biệt nổi → không được dùng để chọn 20 sequence. |

| Kết quả MSA | Phản ứng |
|---|---|
| `msa_source = colabfold_uniref`, depth ≥ 50 | Lưu a3m lên Drive, dùng lại mãi. Coi là a3m chuẩn của ta cho đến khi organizer công bố a3m của họ. |
| `msa_source = query_only_FALLBACK` | **Chưa tái dựng được cửa sàng.** Việc đầu tiên của notebook 02: lấy MSA thật (ColabFold API lúc đỡ tải, hoặc `colabfold_batch`, hoặc platform sponsor). Mọi điểm trước đó chỉ là sàn dưới. |
| `msa_encoder_present = false` | Checkpoint sai → kiểm lại `MODEL_DIR` (có phải lỡ trỏ vào `-Fast`?). Chặn mọi việc khác. |

| Thời gian / fold đo được (cell 8 run A) | Ý nghĩa ngân sách |
|---|---|
| `t` giây cho ~274 aa | 1000 fold ≈ `t·1000/3600` GPU-giờ ≈ `t·1000/3600·4.8` CU. Với t=120 s → 33 giờ / **160 CU** → **vượt xa 90 CU hiện có**. Nghĩa là: không thể tự chấm hàng nghìn design trên L4; phải lọc rẻ trước (AF2/Boltz proxy, hoặc `num_loops` thấp hơn) rồi chỉ chấm top-N bằng cấu hình cửa sàng đầy đủ. **Con số này quyết định toàn bộ kiến trúc notebook 02.** |

---

## 9. Câu hỏi còn mở (CHƯA ĐO / phải hỏi organizer)

Đã đo được bằng notebook 01: peak VRAM, trần residue, thời gian/fold, MSA path, ipSAE path.

**Cập nhật**: notebook 01 chỉ đo được MSA path + ipSAE path + token order; VRAM/thời gian/trần
residue đều **sai vì diffusion decoder chưa chạy** (§0.2). Notebook 02 đo lại.

Không đo được bằng code, **phải hỏi organizer** (xếp theo giá trị):
1. `num_sampling_steps`: **14** (`config.json`) hay **200** (`fold()` default)? Khác nhau rất lớn về cả điểm và thời gian.
2. `lm_dropout`: 0.3 (`fold()`), 0.25 (`config.json`), hay 0 (tắt)? Nó là nguồn ngẫu nhiên.
3. Checkpoint nào: `biohub/ESMFold2` hay `*-Experimental*` / `*-Cutoff2025`?
4. a3m cụ thể cho từng target (database, depth, `msa_max_depth`, `msa_column_mask_rate`).
5. TL1A fold mấy chain target: 1 hay 3? (cấu hình ghi `仅靶点序列使用 MSA` — "sequence", số ít;
   cell 13 của notebook 02 đo delta ipSAE giữa hai cấu hình để biết câu này có chặn ta không)
6. ipSAE: v1 hay IPSAE2? cột `d0res`? hàng `asym` (min/direction nào) hay `max`? cutoff 10/10 hay 15/15?
