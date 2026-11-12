# AGY-03 — Adapter: Mavenus ESMFold2 output → input cho ipsae.py chính thức

Bạn là agent code. Làm đúng một việc, có self-test đáp án biết trước. Không suy đoán.

## Việc

Viết module `tools/mavenus_to_ipsae.py` chuyển output ESMFold2 của platform Mavenus
thành đúng cặp file mà `ipsae.py` (DunbrackLab/IPSAE) nhận ở **đường AF2**:
một file `.pdb` và một file `.json`.

Module KHÔNG chạy ipsae.py. Chỉ sinh input. Tôi chạy ipsae.py ở bước sau.

## ĐÃ BIẾT — đừng tìm lại

Những điều dưới đây đã verify, dùng luôn, đừng search lại, đừng hỏi lại:

1. **Input có sẵn tại** `Mavenus_output/IL7Ra_minibinder_FINAL/`:
   - `metrics_json/<name>.metrics.json` — 38 file
   - `structures/<name>.cif` — 38 file, cùng tên
   - Tên: `c01_A0_4` … `c37_A1_2`, cộng một chứng dương `ctrl_7S5B`

2. **Cấu trúc file metrics.json** (key ở tầng ngoài cùng):
   - `name` (str), `chains` (list, mỗi phần tử có `id`, `header`, `len`)
   - `pae_matrix` — list of list, float, đơn vị Ångström
   - `pae_tokens` — list, mỗi phần tử: `{"index": int, "chain_id": "A"|"B", "residue_index": int, ...}`
     Thứ tự phần tử trong `pae_tokens` **chính là** thứ tự hàng/cột của `pae_matrix`.
   - `ptm` (float), `iptm` (float), `plddt_0_100` (float — chỉ là **giá trị trung bình**, KHÔNG phải per-residue)

3. **File .cif** là định dạng minimal: chỉ có một block `loop_` `_atom_site.*`.
   Các cột gồm `group_PDB, type_symbol, label_atom_id, label_alt_id, label_comp_id,
   label_asym_id, label_seq_id, pdbx_PDB_ins_code, auth_seq_id, ...`
   Phải **đọc header của loop_ để lấy chỉ số cột**, không hardcode vị trí cột.
   Cột B-factor trong file này chứa **pLDDT theo từng atom** (thang 0–100).

4. **Format json mà ipsae.py đường AF2 cần** — đúng 4 key này:
   ```json
   {"pae": [[...]], "plddt": [...], "iptm": 0.0, "ptm": 0.0}
   ```
   - `pae`: ma trận vuông N×N, N = tổng số residue
   - `plddt`: list độ dài N, **per-residue**, thang 0–100
   - `plddt` phải lấy từ B-factor của .cif, gom theo residue (lấy giá trị của atom CA;
     nếu residue không có CA thì lấy trung bình các atom của residue đó)

5. **Format pdb mà ipsae.py cần**: ATOM record chuẩn cột cố định. Thứ tự residue trong pdb
   phải **khớp đúng** thứ tự hàng của ma trận `pae`. Chain id giữ nguyên `A` (receptor) và `B` (binder).

6. **Kích thước biết trước** (dùng làm self-test):
   | name | receptor len | binder len | N tổng | pae_shape |
   |---|---|---|---|---|
   | `ctrl_7S5B` | 196 | 55 | **251** | 251×251 |
   | `c29_A17_5` | 196 | 98 | **294** | 294×294 |

## API phải có

```python
def load_mavenus(name: str, root: str = "Mavenus_output/IL7Ra_minibinder_FINAL") -> dict
def write_ipsae_inputs(name: str, outdir: str, root: str = ...) -> dict
```

`write_ipsae_inputs` trả dict có: `pdb_path`, `json_path`, `n_res`, `n_atoms`,
`chain_lengths` (dict chain_id → số residue).

## Self-test bắt buộc — viết vào `tools/test_mavenus_to_ipsae.py`

Mỗi case có đáp án biết trước. In `PASS`/`FAIL` từng case và tổng cuối.

1. `ctrl_7S5B`: `n_res == 251`, `chain_lengths == {"A": 196, "B": 55}`
2. `c29_A17_5`: `n_res == 294`, `chain_lengths == {"A": 196, "B": 98}`
3. Ma trận `pae` trong json sinh ra phải **vuông** và cạnh `== n_res`
4. `len(json["plddt"]) == n_res`, và mọi giá trị trong khoảng `0 <= x <= 100`
5. PAE đối xứng gần đúng: `max(abs(pae[i][j] - pae[j][i])) < 0.5` trên toàn ma trận
6. Đường chéo PAE nhỏ: `max(pae[i][i]) < 2.0`
7. Thứ tự khớp: residue thứ k trong file pdb phải có `chain_id` bằng `pae_tokens[k]["chain_id"]`
   — kiểm trên cả 251 vị trí của `ctrl_7S5B`
8. File pdb parse lại được: số dòng bắt đầu bằng `ATOM` == `n_atoms`, và mọi dòng dài >= 54 ký tự
9. pLDDT của `ctrl_7S5B` phải có trung bình **trong khoảng 91.0–91.9**
   (giá trị đúng đã biết: `plddt_0_100 = 91.4492`). Case này bắt lỗi gom B-factor sai.

Chạy bằng: `.venv/bin/python3 tools/test_mavenus_to_ipsae.py`

## Ràng buộc

- Chỉ dùng standard library + numpy. **Không** dùng biopython, gemmi, mdanalysis.
- **TUYỆT ĐỐI KHÔNG SỬA/XOÁ**: `data/known_binders/`, `data/msa/`, `Mavenus_output/`,
  `data/submission_*.csv`, bất cứ file nào trong `scripts/`.
  Chỉ được tạo mới 2 file: `tools/mavenus_to_ipsae.py` và `tools/test_mavenus_to_ipsae.py`.
- Không bịa số. Chỗ nào không chắc thì ghi `UNCONFIRMED` trong comment, không đoán.
- Không chạy ipsae.py, không tải gì từ mạng, không chạy GPU.

## Xong là khi

9/9 self-test PASS, in ra bảng kết quả. Báo lại số case pass và đường dẫn 2 file đã tạo.
