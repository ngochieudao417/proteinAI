# Thuê GPU có SSH — hướng dẫn cho bạn

**Trả lời câu hỏi của bạn: ĐƯỢC, và đây là phương án tốt nhất.** Khác hẳn SSH vào Colab:
với máy thuê thì SSH là tính năng chính thức của dịch vụ, không vi phạm điều khoản nào.

## Vì sao nên làm, theo thứ tự quan trọng

**1. Thời gian — lý do mạnh nhất.** Mốc thực tế là **13–14/11**, còn 2–3 ngày (xem `00-COMPETITION-FACTS.md` §1:
API `/api/public/meta` vẫn ghi hết hạn 11-15 và chưa được sửa). Vòng lặp notebook-thủ-công cần một lượt
qua lại cho mỗi lần sửa. Qua SSH tôi chạy trực tiếp, sai thì sửa trong vài giây. Với 3 ngày, đây là khác biệt
giữa làm xong và không làm xong.

**2. Mở khoá TL1A — thứ L4 không thể làm.** Epitope TL1A nằm ở **groove giữa hai protomer** của trimer
(research agent tính được W137/F138 cách protomer kế bên 2,9–3,4 Å). Muốn đánh giá đúng một binder nhắm vào
groove đó thì **phải fold cả trimer**: 3×180 + binder = **580–640 residue**. L4 24GB kham ~500 residue → không đủ.
A100 80GB kham ~2000 → thoải mái. Nếu chỉ có L4, ta buộc phải fold TL1A dạng monomer và tự lừa mình rằng
điểm đó phản ánh thực tế.

**3. Hết lo compute unit.** 90 CU (~18,7 giờ L4) là ràng buộc chặt nhất hiện tại. Thuê GPU xoá nó.
Và không có giới hạn 12 giờ/session, không mất runtime giữa đường, đĩa persistent cho 26GB weights ESMFold2.

## Chi phí thực tế

Giá tháng 11/2026:

| Provider | A100 80GB | H100 | Ghi chú |
|---|---|---|---|
| **RunPod** Community | **$1,19/h** | $1,99/h | dễ dùng, SSH sẵn, ổn định hơn |
| **RunPod** Secure | $1,79/h | $2,89/h | |
| **Vast.ai** | $1,68/h | $0,90–1,55/h | rẻ nhất nhưng chất lượng host thay đổi |

Nhu cầu ước tính: **20–30 giờ** (scoring ESMFold2 + sinh design cho 2 target)
→ **khoảng $25–35 trên RunPod A100 80GB.**

Giải nhất một track là **60.000 RMB** (~$8.400). Tỉ lệ $30 so với $8.400 là 1:280.

Nguồn giá: [RunPod pricing](https://www.runpod.io/pricing) · [so sánh provider 11/2026](https://altstreet.investments/tools/gpu/prices) · [RunPod vs Vast.ai 2026](https://www.spheron.network/blog/runpod-vastai-2026/)

## Việc bạn làm (tôi không thể tự làm: cần tạo account + nhập thẻ)

Tôi **không** tạo account hay nhập thông tin thanh toán thay bạn được. Phần đó bạn làm.

### Đề xuất: RunPod (dễ hơn Vast.ai cho lần đầu)

1. Tạo account: https://www.runpod.io/ → nạp **$30** (đủ dư; pay-as-you-go, không dùng hết thì còn lại)
2. **Trước khi tạo pod**, thêm SSH public key: `Settings → SSH Public Keys`
   → tôi sẽ tạo keypair và đưa bạn **public key** để dán vào đó (chỉ public key, an toàn hoàn toàn)
3. `Deploy` → chọn:
   - GPU: **A100 80GB PCIe** (Community Cloud cho rẻ). Nếu hết hàng: **L40S 48GB** hoặc **A6000 48GB** cũng đủ.
     ⚠️ Đừng chọn dưới 40GB — sẽ lặp lại đúng vấn đề của L4.
   - Template: **RunPod PyTorch 2.x** (có CUDA sẵn)
   - Disk: **Container 50GB + Volume 100GB** (weights ESMFold2 F32 là 26GB, cần chỗ)
   - Bật **SSH over exposed TCP** (không chỉ web terminal)
4. Pod chạy → copy dòng `ssh root@<host> -p <port> -i ~/.ssh/...` mà RunPod hiển thị → **gửi dòng đó cho tôi**
5. Xong việc: **nhớ bấm Terminate**, không chỉ Stop — pod stopped vẫn tính tiền lưu trữ

### Lưu ý tiền

- Tính theo giây, chỉ khi pod chạy. Tôi sẽ nói rõ khi nào nên tắt.
- Đặt **spending limit** trong account để không bao giờ cháy quá dự kiến.
- Community Cloud có thể bị thu hồi máy (rẻ hơn vì vậy). Mọi thứ quan trọng tôi ghi ra volume persistent,
  nên bị thu hồi thì mất tối đa một bước, không mất cả tiến trình.

## Phân vai sau khi có máy

| Nơi | Việc |
|---|---|
| **GPU thuê** (A100) | sinh design + chấm ESMFold2/ipSAE cho **cả hai** target, gồm TL1A trimer full |
| **Colab Pro** (90 CU) | để dành làm dự phòng, hoặc chạy song song phần nhẹ. Không còn là đường tới hạn |
| **Platform sponsor** | đối chiếu chéo: nếu ipSAE của MatwingsVenus khớp số của ta thì ta đã tái dựng đúng cửa sàng. Rất đáng để verify |
| **Laptop của bạn** | novelty check (alignment trên CPU, không cần GPU) — chống 122.273 design Cao 2022 |

## Nếu bạn không muốn tiêu tiền

Vẫn làm được, kế hoạch B:
- IL-7Rα: chấm trên Colab L4 (259–319 aa fit) — track này không cần GPU lớn
- TL1A: đẩy sang platform sponsor (AsterFire nói compute không giới hạn, MatwingsVenus có sẵn ESMFold2+ipSAE),
  hoặc chấp nhận fold monomer và ghi rõ đó là xấp xỉ
- Thiệt hại thật: chậm hơn (mỗi vòng lặp mất một lượt qua lại), và TL1A bị đánh giá bằng thước đo gần đúng

Không có gì sai với kế hoạch B. Chỉ là với 3 ngày còn lại thì nó mỏng.
