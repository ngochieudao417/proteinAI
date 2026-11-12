# Checklist đăng ký platform sponsor — **bạn tự làm phần này**

Tôi không đăng ký account / nhập credential thay bạn được. Dưới đây là thứ tự ưu tiên + đúng những gì cần bấm.
Mỗi mục có ô tick — tick xong thì nói tôi biết, tôi sẽ điều chỉnh pipeline theo đúng cái gì bạn có.

**Lưu ý chung:** các platform này đều là của TQ, phần lớn yêu cầu **số điện thoại TQ (+86)** hoặc **WeChat** để đăng ký. Nếu bạn không có, hãy thử email trước; chỗ nào bắt buộc +86 thì ghi lại và bỏ qua — ta vẫn chạy được pipeline trên Colab mà không cần họ (xem mục "Nếu không đăng ký được gì" ở cuối).

---

## 🔴 ƯU TIÊN 0 — Trước hết: xác nhận đăng ký THI

Cái này quan trọng hơn mọi platform. Nếu chưa đăng ký thi thì mọi việc khác vô nghĩa.

- [ ] **Đã đăng ký đội chưa?** Trang: https://aiprotein.sjtu.edu.cn
  - Hạn đăng ký: **2026-11-30** (đã gia hạn từ 10-31)
  - Login: username = **email captain**, password ban đầu = **số điện thoại captain**
  - Sau khi đăng ký phải **chờ 1 ngày làm việc** mới login được
- [ ] Login được → **đổi password ngay** (không đổi thì không tạo team / submit được)
- [ ] Vào trang submit, **chụp screenshot** cho tôi xem:
  - Nút upload có bấm được không (có người báo nút xám từ 15/11)
  - Deadline hiển thị trên trang là ngày nào
  - Form yêu cầu đúng những field gì
- [ ] Nếu nút submit xám / deadline sai → nhắn official support trong group WeChat (có người đã báo, organizer trả lời "您私信我")

---

## 🟠 ƯU TIÊN 1 — MatwingsVenus® 晓鹜 (quan trọng nhất)

**Lý do ưu tiên:** đây là platform duy nhất đã công bố có **ESMFold2 + tính ipSAE** online — tức là đúng cái thước đo mà organizer dùng để lọc vào wet-lab. Có nó, ta tự check điểm trước khi submit thay

Tutorial: VHH https://mp.weixin.qq.com/s/wl-h4Av1CjpGimIh2k7McA · minibinder https://mp.weixin.qq.com/s/I1LNXzNwrEN0heBoLDMJOg · ESMFold2 https://mp.weixin.qq.com/s/5ef89061e1247218f0d0432148385a89

---

## 🟠 ƯU TIÊN 2 — AsterFire™ / 星使智算 (compute không giới hạn)

**Lý do:** họ nói rõ *"比赛期间算力不限额，免费用"* (compute không giới hạn trong thời gian thi). Và đã đăng case study thật cho **TL1A minibinder** — đúng 1 trong 2 track minibinder, có số ipSAE cụ thể (candidate_0002 ipSAE 0.637, candidate_0015 0.628, binder pLDDT 91.8) + đã novelty-check với 84 patent anti-TL1A VH. Đây là baseline tham chiếu rất giá trị.

- [ ] Đăng ký: https://sidereus-ai.com/ (workbench: https://sidereus-ai.com/workbench)
- [ ] **Điền form mở quyền miễn phí**: https://sidereus-ai.feishu.cn/share/base/form/shrcncNebxQyUFdRsrZZKSuhtah
- [ ] Báo tôi: có **Protenix**, **Mol-1**, post-processing filter, MD không? Có API/agent hook không (họ nói support "codex và dsh" tức gọi được từ agent bên ngoài → **nếu có, tôi dùng được trực tiếp**)
 
Tutorial: https://mp.weixin.qq.com/s/_Aemdfh_e8gDjyEc9GuBug · case TL1A: https://mp.weixin.qq.com/s/6629f3f0206b10b65c5f7ca3c31de966

---

## 🟡 ƯU TIÊN 3 — SPONGE Mokda™ (MD / đánh giá ổn định)

**Lý do:** dùng để lọc cuối (MD, binding stability). Không cần sớm, nhưng form duyệt có thể mất thời gian nên đăng ký trước.

- [ ] Đăng ký + tải: https://sponge-mokda.cn/ hoặc https://sponge-mokda.com/download/
- [ ] **Điền form mở quyền**: https://sidereus-ai.feishu.cn/share/base/form/shrcngswAROGVTiVVwfe0pqXP1b
- [ ] Lưu ý cách dùng: design nằm ở **toolbox → RFdiffusion**; MD thì phải bấm **建模 (modeling)** ở thanh trên để sinh input "raw" trước

Tutorial: https://mp.weixin.qq.com/s/e7W2nTAJeXV8TadSweGvaw · có GLP-1R cyclic-peptide scaffold tutorial

---

## 🟡 ƯU TIÊN 4 — 书生·端砚 InternAI Discovery

- [ ] Đăng ký: https://discovery.intern-ai.org.cn/chat
- [ ] Vào WeChat group "书生·端砚交流群" **số 2** (group 1 đã full) — cách nhận compute ("墨点值" / TokenPlan) nằm ở **group announcement**
- [ ] Giá trị thật: **protein dataset (sequence + structure)** + lit search. Nếu họ có list known VHH/minibinder theo target thì cực hữu ích cho novelty check → hỏi thử.

Tutorial: https://mp.weixin.qq.com/s/aOBqOq6ZWo9Q3KE5nOlCvg

---

## 🟢 ƯU TIÊN 5 — OriBench 元生·开物 (optional)

- [ ] Đăng ký / tải: https://oribench.lglab.ac.cn/ · docs: https://oribench.lglab.ac.cn/docs.html
- [ ] Họ nói *"比赛期间免费提供Token与平台服务，开箱即用，无需额外报名兑换"* → không cần form, login là dùng
- Thiên về LLM agent / drug-research skills hơn là structure modeling. Ưu tiên thấp.

---

## Những gì tôi cần từ bạn (theo mức hữu ích)

1. **Screenshot trang submit** sau khi login aiprotein.sjtu.edu.cn → chốt được deadline + format CSV thật
2. **Colab tier** của bạn: Free / Pro / Pro+? GPU nào bấm được (T4 / L4 / A100)? Còn bao nhiêu compute unit?
   → quyết định ta chạy RFdiffusion/BindCraft local-Colab hay dựa vào platform sponsor
3. **Google Drive** còn trống bao nhiêu GB (weights RFdiffusion + ESMFold ~5–20GB, cần cache để không tải lại mỗi session)
4. Account nào đăng ký **được** / **không được** (vướng +86)
5. Nếu MatwingsVenus / AsterFire có **API key hoặc endpoint** → paste cho tôi *cách gọi* (đừng paste chính key vào chat; để vào file `.env` trong repo rồi nói tôi tên biến)

## Nếu không đăng ký được platform nào

Vẫn làm được. Plan B hoàn toàn trên Colab GPU:
- RFdiffusion / RFantibody (weights public) → sinh backbone
- ProteinMPNN / LigandMPNN → sinh sequence
- ESMFold hoặc Boltz / Chai-1 (public weights) → predict phức hợp
- ipSAE (script public, từ repo của Dunbrack lab) → score
Thiếu duy nhất là **ESMFold2 chính xác như organizer dùng**; research agent đang truy xem nó là model nào và có public weights không. Nếu không có, ta dùng proxy score (ipSAE từ model gần nhất) + nhiều filter khác để bù.
