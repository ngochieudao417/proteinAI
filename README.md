# proteinAI — 长江杯 2026 全球AI蛋白质设计挑战赛

> 🔸 **TẠM DỪNG 2026-11-11.** Đọc **`docs/07-STATE-AND-HANDOFF.md`** trước khi làm gì tiếp —
> nó chứa trạng thái, bug đang chặn, 3 việc kế tiếp, và danh sách cạm bẫy đã trả giá để biết.

**Trạng thái:** ngày 11/11/2026. Mốc làm việc thật: **chốt 13/11, confirm 14/11** (không phải 30/11).
**Tư cách:** public team (tự đăng ký) → sequence phải qua **cửa sàng ESMFold2 + ipSAE** mới được đem đo KD thật.
**Track đã chọn:** 2 track minibinder — **TL1A** (O95150, 180 aa, trimer) và **IL-7Rα** (P16871, 219 aa).
**Compute:** Colab Pro (L4 24GB / T4), Google Drive >50GB.

## Hai mục tiêu, phải thắng cả hai

1. **Cửa sàng (in-silico).** 1000 slot/track, chấm bằng ipSAE trên structure do ESMFold2 dự đoán, cấu hình cố định: 1 seed, 1 sample, no template, MSA chỉ cho target. Vì nó *tiền định*, ta tái dựng được tại chỗ và chọn 20 sequence theo đúng thước đo đó.
2. **Xếp hạng cuối (wet-lab).** Thuần theo KD đo thật. Score in-silico cao mà binder không thật thì vô nghĩa.

Hai mục tiêu này **không đồng nhất** — chỗ chúng xung đột là nơi quyết định thắng thua.

## Cấu trúc repo

```
docs/
  00-COMPETITION-FACTS.md   ← fact base, trích từ JS bundle của site chính thức
  01-SPONSOR-SIGNUP.md      ← 👉 VIỆC CỦA BẠN: checklist đăng ký platform sponsor
  02-RESEARCH.md            ← research agent: ESMFold2 là gì, ipSAE, epitope, pipeline
  03-VERIFY.md              ← verify agent: deadline thật, structure, sequence, novelty DB
data/
  TL1A_O95150_L72-L251.fasta    ← target 1, đã cắt theo range, verify khớp đầu/cuối
  IL7Ra_P16871_E21-D239.fasta   ← target 2
  targets_raw_from_bundle.txt   ← raw data 7 target trích từ bundle
  structures/                   ← 2RE9, 3DI3
  known_binders/                ← DB để tự check novelty (≤70% identity)
notebooks/
  00_probe_env.ipynb        ← 👉 VIỆC CỦA BẠN: chạy trên Colab, paste kết quả lại
```

## Phân việc

Vòng lặp làm việc: **tôi viết notebook → bạn chạy trên Colab → bạn gửi lại JSON/file.**
Không SSH vào Colab (ToS cấm remote shell, rủi ro mất credit), không tailscale (tailnet là workspace tổ chức).
Chi tiết quy ước: `docs/04-WORKFLOW.md`.

| Ai | Việc |
|---|---|
| **Bạn** | Đăng ký platform sponsor (`docs/01-SPONSOR-SIGNUP.md`), chạy `notebooks/00_probe_env.ipynb` rồi gửi JSON, screenshot trang submit, bất cứ gì cần +86/WeChat/login |
| **research agent** | ESMFold2 là model nào · ipSAE formula + repo · epitope TL1A/IL-7Rα · chọn pipeline chạy được trên L4 · bias của ipSAE |
| **verify agent** | Deadline thật (15/11 vs 30/11) · format CSV chính xác · fetch + mô tả 2RE9/3DI3 · verify sequence · thu DB known binder · tính độ cạnh tranh của cửa sàng |
| **code agent** | (chờ research chốt pipeline) notebook sinh backbone → sequence → fold → score ipSAE → filter → CSV submit |
| **tôi (điều phối)** | fact base, cắt target sequence, review output của agent, chốt quyết định |

## Rủi ro đã nhận diện

- **ESMFold2 chưa rõ là model nào.** Config sàng ghi "MSA cho target" nhưng ESMFold gốc *không bao giờ* dùng MSA → "ESMFold2" không phải ESMFold v1. Nếu không tái dựng được đúng model, ta chấm bằng proxy và mất độ chính xác khi chọn 20 sequence. Đang truy.
- **Deadline mâu thuẫn.** Group chat nói 30/11, frontend còn code 15/11, có người báo nút submit xám từ 15/11. → **lập kế hoạch theo mốc 15/11 cho an toàn**, coi 30/11 là phần đệm.
- **Cửa sàng có thể quá tải nặng.** 1200+ team đăng ký; nếu nhiều team dồn vào cùng track thì 1000 slot rất chật. Verify agent đang ước lượng percentile cần đạt.
- **Colab Pro có hạn.** L4 ~4.8 compute unit/giờ. Pipeline 2 target cần cỡ 15–40 giờ GPU → phải biết số CU còn lại trước khi cam kết.
