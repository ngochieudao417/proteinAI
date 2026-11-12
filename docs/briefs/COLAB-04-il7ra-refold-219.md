# COLAB-04 — Re-fold IL-7Rα trên receptor 219 aa chính thức

Mục đích: đóng lỗ hổng cuối cùng của track IL-7Rα. MatvenUS fold trên receptor **196 aa**
(7OPB 36–231), còn cửa sàng của BTC dùng construct **219 aa** (P16871 Glu21–Asp239).
Cần số ipSAE đo trên đúng construct trước khi chốt 20 slot.

## Tại sao chưa thể bỏ qua

- `ipTM` trong file MatvenUS tính trên 196 aa, **không phải** thước đo cửa sàng.
- ipSAE đã tính lại được từ `pae_matrix` của MatvenUS (xem `data/04_ipsae_IL7RA_recomputed.csv`),
  nhưng vẫn trên 196 aa.
- Hiệu chuẩn bằng đối chứng dương: 7S5B trên 196 aa cho ipSAE_min **0,7638**; Notebook 03 của ta
  trên 219 aa cho **0,7488**. Lệch **−0,015**, tức receptor ngắn làm điểm **cao hơn ~0,015**.
  Nhỏ, nhưng chưa đo trực tiếp trên từng design.

## Việc

Chạy pipeline sàng lọc có sẵn trên Colab L4, target IL7RA, 37 design MatvenUS.

```bash
python scripts/run_screening_colab.py \
  --target IL7RA \
  --candidates-file /content/drive/MyDrive/proteinAI/data/05_il7ra_refold_candidates.csv \
  --out-dir /content/drive/MyDrive/proteinAI/out_il7ra_refold219 \
  --seed 0
```

Cần copy `data/05_il7ra_refold_candidates.csv` lên `MyDrive/proteinAI/data/` trước.
Weights đã cache ở `MyDrive/proteinAI/weights/ESMFold2` — **đừng tải lại**.

## Chi phí dự kiến

| | |
|---|---|
| Số design | 37 |
| Thời gian fold | ~40 s/design → ~25 phút |
| Load model | ~3 phút |
| CU tiêu tốn | **~2,2 CU** (còn 81,94) |

File đã xếp sẵn theo ipSAE giảm dần, nên nếu đứt mạng giữa đường thì các design tốt nhất
đã fold xong trước. Cell checkpoint-resume tự bỏ qua design đã có.

## ĐÃ BIẾT — đừng tìm lại

1. `scripts/run_screening_colab.py` đã nạp đúng construct 219 aa (`IL7RA_SEQ`, `len: 219`,
   `IL7Ra_P16871_E21-D239`) và MSA depth 921 từ `data/msa/`. Fold IL-7Rα dạng **monomer**
   chain A + binder chain B. Không cần sửa gì.
2. Cấu hình ESMFold2 đúng: `dtype=None, esmc_precision="fp8", kernel_backend="fused"`.
   VRAM phức 274–317 aa khoảng 8,3–8,7 GB, rất an toàn trên L4 22 GB.
3. Gate chuẩn: 20 loops, 200 diffusion steps, 1 sample, seed 0, no-template.
4. Đối chứng dương 7S5B trên 219 aa: ipSAE_binder 0,7728 / ipSAE_min 0,7488, ipTM 0,9374.
   Dùng con số này để xác nhận phiên chạy không lệch.
5. ipSAE cutoff 10/10, lấy **đường AF2** (`.pdb` + json có `pae/plddt/iptm/ptm`).
   **Không** dùng đường AF3 `.cif`.
6. Epitope đã xác minh đúng: MatvenUS dùng hotspot A97–A103 theo numbering 7OPB,
   bằng đúng `KKFLLIG` construct 77–83. Không có lỗi offset.

## ĐẦU RA

- `out_il7ra_refold219/03_screening_results.csv` — toàn bộ 37 design, có `ipsae_binder`, `ipsae_min`
- `out_il7ra_refold219/03_screening_top_IL7RA.csv` — top 20
- Dán khối `PASTE CHO CLAUDE` về đây.

## Tiêu chí "xong"

- 37/37 design có `ipsae_min` trên receptor 219 aa.
- Nếu `ipsae_min` của c29 vẫn ≥ 0,75 → chốt c29 làm leader, không cần làm gì thêm.
- Nếu top 20 theo số mới khác `data/submission_IL7RA_PROPOSED_v2.csv` → xếp lại theo số mới.
