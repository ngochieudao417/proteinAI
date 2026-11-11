# Quy ước làm việc: tôi viết notebook → bạn chạy → bạn gửi lại kết quả

Đã chốt: **không SSH vào Colab, không tailscale.** Lý do: Colab ToS xếp "remote control such as SSH shells"
vào mục bị cấm, và credit Colab là tài nguyên chính của dự án này — không đánh cược nó.
Tailnet của bạn là workspace tổ chức nên càng không chạm vào.

Vậy vòng lặp là: **tôi viết notebook → bạn chạy trên Colab → bạn gửi lại 1 block JSON hoặc 1 file.**
Hơi chậm hơn, nhưng an toàn và không mất gì về chất lượng kết quả.

## Mọi notebook tôi viết sẽ tuân đúng khuôn này

1. **Cell markdown đầu**: notebook này làm gì, cần GPU loại nào, chạy mất bao lâu, tốn khoảng bao nhiêu compute unit.
2. **Cell 1 — tiền kiểm**: kiểm GPU/VRAM/Drive trước khi cài gì. Thiếu điều kiện thì `raise SystemExit` kèm lời nhắc cụ thể,
   để bạn không ngồi chờ 20 phút rồi mới biết sai runtime.
3. **Các cell giữa**: mỗi cell làm **một** việc, chạy độc lập được, in tiến độ. Việc nào >10 phút thì in nhịp,
   để bạn biết nó đang chạy chứ không treo.
4. **Mọi output quan trọng ghi vào `/content/drive/MyDrive/proteinAI/`**, không bao giờ chỉ ghi `/content`.
   Runtime chết là mất sạch `/content`; coi như nó có thể chết bất cứ lúc nào.
5. **Cell cuối — BÁO CÁO**: luôn in ra một block có mốc rõ ràng:

   ```
   ========================================
   PASTE CHO CLAUDE — BẮT ĐẦU
   ========================================
   { ...json... }
   ========================================
   PASTE CHO CLAUDE — KẾT THÚC
   ========================================
   ```

   và đồng thời ghi cùng nội dung đó ra `out/<tên_notebook>_report.json` trên Drive.
   → Bạn chọn: copy block text, **hoặc** tải file JSON gửi tôi. Cái nào tiện hơn thì làm.

## Bạn gửi lại cho tôi thế nào

**Cách tiện nhất:** tải file `out/*_report.json` từ Drive rồi gửi thẳng vào chat.
**Cách nhanh:** copy nguyên block giữa hai dòng `PASTE CHO CLAUDE`.

Ba điều giúp tôi nhiều nhất:

- **Lỗi đỏ thì copy nguyên traceback.** Đừng tóm tắt thành "nó lỗi" — dòng cuối của traceback thường là
  toàn bộ thông tin tôi cần để sửa. `CUDA out of memory` và `ModuleNotFoundError` dẫn tới hai cách sửa hoàn toàn khác nhau.
- **Thời gian chạy thật** của cell nào lâu (Colab in sẵn ở góc dưới mỗi cell). Tôi dùng nó để tính budget cho các bước sau.
- Những gì **code không đọc được**, notebook sẽ đánh dấu `<<ĐIỀN TAY>>` — chủ yếu là số compute unit còn lại
  và screenshot từ web thi đấu.

Nếu notebook chạy nửa đường rồi runtime chết: cứ gửi phần output có được. Thiết kế notebook là
resume được từ Drive, không phải chạy lại từ đầu.

## File output cố định trên Drive

```
/content/drive/MyDrive/proteinAI/
  weights/   ← model weights, cache 1 lần dùng mãi (đừng xoá: tải lại rất lâu)
  data/      ← target sequence, structure, known-binder DB
  runs/      ← output thô từng lần chạy design (backbone PDB, v.v.)
  out/       ← báo cáo JSON + CSV submit cuối. Đây là thứ bạn gửi tôi.
  logs/      ← log dài, khi cần tôi sẽ nhờ bạn gửi file cụ thể
```

## Thứ tự notebook

| Notebook | Việc | Trạng thái |
|---|---|---|
| `00_probe_env.ipynb` | Đo GPU/VRAM/Drive/mạng → quyết định pipeline nào khả thi | ✅ sẵn sàng, **chạy cái này trước** |
| `01_*` | (chờ research agent chốt ESMFold2 + pipeline) sinh backbone binder | ⏳ đang chờ |
| `02_*` | Sinh sequence (ProteinMPNN) + fold lại + chấm điểm | ⏳ |
| `03_*` | Tái dựng đúng cửa sàng (ESMFold2 + ipSAE, 1 seed 1 sample no template) → chọn 20 sequence | ⏳ |
| `04_*` | Novelty check + xuất CSV đúng format submit | ⏳ |
