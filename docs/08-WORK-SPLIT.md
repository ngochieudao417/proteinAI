# Phân việc Claude ↔ Antigravity

## Thực tế kỹ thuật (đã kiểm 2026-11-11)

**Antigravity GUI: không có CLI.** `antigravity 2.16.0` là app Electron fork VS Code.
Package AUR không có `bin/` wrapper (chỗ VS Code để `cli.js`), và `app.asar` không chứa `cli.js`.
`antigravity --help` chỉ mở GUI. → không gọi agent của GUI bằng lệnh được.

**NHƯNG `gemini` CLI đã cài sẵn** ở `~/.local/bin/gemini` (+ node v26.7.0 qua mise), và nó **headless, script được**.
Hai thứ dùng chung thư mục config `~/.gemini/` (có cả `~/.gemini/antigravity/`), nên **dùng chung auth**.

→ **Claude gọi trực tiếp agent Gemini được.** Không cần copy-paste, không cần bạn làm trung gian.

### Lệnh chuẩn

```bash
export PATH="$HOME/.local/bin:$PATH"
cd /home/ngochieu/repo/proteinAI
gemini -p "$(cat docs/briefs/AGY-01-novelty-checker.md)" \
       --approval-mode auto_edit \
       --skip-trust \
       -o json \
       --include-directories /home/ngochieu/repo/proteinAI
```

Các flag quan trọng:

| Flag | Vì sao |
|---|---|
| `-p` | chế độ headless (không có nó thì vào interactive, treo) |
| `--approval-mode` | `plan` = read-only (phân tích) · `auto_edit` = tự duyệt sửa file · `yolo` = tự duyệt mọi thứ kể cả chạy lệnh |
| `--skip-trust` | tin workspace trong session này. Không có nó → bị hạ về `default` rồi treo chờ người duyệt |
| `-o json` | output máy đọc được, để Claude parse kết quả |
| `-m` | chọn model, nếu cần hạ xuống model rẻ cho việc cơ học |

**Mức approval theo loại việc:** phân tích/đọc → `plan`. Viết code + self-test → `auto_edit`.
Cần `pip install`/chạy lệnh hệ thống → `yolo`, và chỉ khi git đã sạch.

### An toàn: git là lưới đỡ

Repo **trước đây không có VCS** → agent auto-approve mà xoá nhầm là mất vĩnh viễn.
Đã `git init` + commit snapshot (47 file, commit `11f06e5`).
`data/known_binders/` (122k sequence) và `data/msa/` (MSA thật depth 921/2037) đều đắt và không lấy lại từ cache được.

**Luật: trước mỗi lần chạy agent ở `yolo`, `git status` phải sạch.** Hỏng thì `git checkout -- .`.

## Nguyên tắc phân việc: theo KIỂU THẤT BẠI, không theo độ khó

Đây là điểm quan trọng nhất của cả file này.

Sai lầm thường gặp là "việc dễ giao model yếu". Không đúng — vì có những việc **dễ làm nhưng sai âm thầm**,
và đó mới là loại giết dự án. Cách phân đúng:

| | Giao **agy** | Giữ cho **Claude** |
|---|---|---|
| Thất bại biểu hiện thế nào | **Ồn ào** — crash, test đỏ, output sai rõ ràng | **Âm thầm** — ra số hợp lý nhưng sai |
| Kiểm lại tốn gì | **Rẻ** — chạy self-test là biết | **Đắt** — phải đọc source mới biết |
| Ví dụ | alignment, parse file, build CSV, scaffold notebook, vòng lặp debug có traceback rõ | chọn epitope, đọc luật thi, verify claim bằng source, thiết kế thước đo |

**Hệ quả bắt buộc: mọi brief cho agy phải có self-test với đáp án biết trước.**
Cách dùng model yếu an toàn không phải là tin nó, mà là **làm cho output của nó kiểm được**.
Nếu tôi không nghĩ ra được self-test cho một việc → việc đó **không** giao agy.

### Ví dụ cụ thể trong dự án này

- Bug `dtype` vừa rồi: lỗi **ồn ào** (có traceback `mat1 and mat2 must have the same dtype`)
  → đáng lẽ giao được agy. Nhưng **nguyên nhân** nằm trong source `esm` (autocast đóng ở `model.py:1160`,
  ba head cố tình giữ fp32) → tìm nguyên nhân là việc Claude, *áp dụng* bản sửa là việc agy.
- Cạm bẫy `ESMFold2-Fast` bỏ qua MSA im lặng: thất bại **âm thầm**, không crash, số vẫn đẹp
  → tuyệt đối không giao agy phần chọn checkpoint / cấu hình chấm điểm.
- Novelty check: 122.273 sequence, so identity, ra một con số. Sai thì self-test bắt được ngay
  (nạp chính một known binder vào → phải ra 100%) → **giao agy, lý tưởng.**

## Hàng đợi việc

### Giao agy ngay

| Brief | Việc | Tại sao an toàn |
|---|---|---|
| `AGY-01` | **Novelty checker** — so 20 design của ta với 122.273 known binder, báo max identity + coverage theo đúng luật thi | Self-test: known binder → 100%; sequence random → thấp. Deterministic. |
| `AGY-02` | **Submission builder** — sinh + validate CSV đúng format site | Đã có `data/validate_submission_csv.py` làm bản tham chiếu để đối chiếu |

### Giao agy sau (khi có đầu vào)

- Scaffold notebook theo khuôn `04-WORKFLOW.md` (copy pattern từ notebook 01/02 đã có)
- Trích hotspot residue từ structure ra file spec máy đọc được cho RFdiffusion
- Soạn **design summary** (material 2 khi submit) từ dữ kiện tôi cung cấp
- Vòng lặp debug pipeline sinh design, khi lỗi có traceback rõ

### Claude giữ

- Đọc/diễn giải luật thi (đã làm, nhưng còn phát sinh)
- **Thiết kế và verify thước đo chấm điểm** (ESMFold2 + ipSAE) — nơi sai là âm thầm và đắt nhất
- Chiến lược epitope, hedge monomer-vs-trimer
- Review mọi output của agy
- Verify claim của bất kỳ agent nào bằng source thật
  *(research agent từng khẳng định "không có CPU offload" — sai; tôi từng khẳng định phải dùng bf16 — cũng sai.
  Cả hai đều bị bắt bằng cách đọc source. Đừng xây lên claim chưa kiểm.)*

## Luật viết brief cho agy

1. **Tự chứa.** agy không đọc được hội thoại này. Mọi số, mọi đường dẫn phải nằm trong brief.
2. **Mục "ĐÃ BIẾT — đừng tìm lại".** Nếu thiếu, nó sẽ đốt quota tìm lại thứ ta đã trả giá để biết
   (ESMFold2 là model nào, công thức ipSAE, epitope ở đâu).
3. **Một việc, một tiêu chí xong đo được.** Không giao việc kiểu "tối ưu pipeline".
4. **Self-test bắt buộc**, có đáp án biết trước ghi thẳng trong brief.
5. **Cấm phần suy đoán.** Ghi rõ: không chắc thì ghi `UNCONFIRMED`, **tuyệt đối không bịa sequence**.
6. **Nói rõ cái gì KHÔNG được sửa** — `data/known_binders/`, `data/msa/`, weights trên Drive.
