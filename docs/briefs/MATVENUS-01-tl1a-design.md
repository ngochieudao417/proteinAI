# MATVENUS-01 — Sinh design minibinder cho TL1A

> Cho bạn chạy trên https://matvenus.com/chat/agent
> **Prompt viết bằng tiếng Trung** (ngôn ngữ gốc của platform → giảm nguy cơ agent hiểu sai thuật ngữ).
> Phần giải thích tiếng Việt ở dưới mỗi khối, đọc để biết đang yêu cầu gì.

## Vì sao TL1A trước, không phải IL-7Rα

| | TL1A | IL-7Rα |
|---|---|---|
| De novo minibinder đã công bố | **0** | **122.273** (Cao 2022) |
| Rủi ro novelty | sạch | đã đo: ~9,2% design vượt ngưỡng |
| Đối thủ | ai cũng phải làm từ đầu | mọi team đều tìm được cùng bộ design đó |

TL1A là track **sạch** — 20 sequence của ta không phải cạnh tranh với 122k design có sẵn.
Đó là lợi thế thật, nên dồn vào đây trước.

---

## PROMPT 1 — Khởi động, kiểm năng lực platform

Dán nguyên khối này vào chat trước. Nó rẻ và trả lời câu quan trọng: platform có chạy được
**ESMFold2 + ipSAE đúng cấu hình cửa sàng** hay không.

```
我在参加"长江杯"全球AI蛋白质设计挑战赛，TL1A赛道（minibinder，40–100 aa）。
在开始设计之前，请先确认平台当前可用的工具与参数能力，逐项回答：

1. 可用的 binder 从头设计工具有哪些？是否包含 RFdiffusion、BindCraft、BoltzGen？
   各自能否指定 hotspot 残基、能否指定 binder 长度范围（我需要 85–100 aa）？
2. 序列设计工具：是否有 ProteinMPNN / AntiBMPNN？能否对同一 backbone 生成多条序列？
3. 结构预测与打分：平台的 ESMFold2 是哪一个 checkpoint？
   具体问题：是完整版 biohub/ESMFold2，还是 ESMFold2-Fast？
   这一点对我很关键，因为 -Fast 版本会静默忽略 MSA 输入。
4. ESMFold2 能否按以下精确配置运行：1 seed、1 sample、无模板（no template）、
   仅靶点链使用 MSA、设计的 binder 链不使用 MSA？
5. 能否计算 ipSAE 并导出数值？能否同时导出复合物 PDB 文件和 PAE 矩阵？
6. 批量能力：一次最多能跑多少个设计？是否有 API 或批处理，还是只能对话式逐个调用？
7. 我的账户当前可用积分是多少？上面每一步大致消耗多少积分？

请直接回答，暂时不要开始任何计算任务。
```

**Giải thích:** câu 3 và 4 là câu đáng tiền nhất. Nếu MatvenUS chạy được ESMFold2 **full** đúng
cấu hình cửa sàng, ta có bản tái dựng cửa sàng **miễn phí**, không đốt compute unit Colab.
Nếu nó dùng `-Fast` thì mọi điểm ipSAE nó đưa ra đều **không đáng tin** cho mục đích sàng lọc,
vì `-Fast` bỏ qua MSA mà không báo lỗi — ta sẽ tưởng đã tái dựng cửa sàng nhưng không.
Câu 6 quyết định ta làm 1 campaign lớn hay nhiều lượt nhỏ. Câu 7 để tôi chia ngân sách.

**Gửi lại tôi toàn bộ câu trả lời** trước khi chạy PROMPT 2.

---

## PROMPT 2 — Campaign thiết kế chính

Chỉ chạy sau khi có câu trả lời PROMPT 1 (vì có thể phải đổi tên tool cho khớp platform).

```
任务：为 TL1A 从头设计 minibinder。

== 靶点 ==
名称：TL1A / TNFSF15
UniProt：O95150，残基范围 Leu72–Leu251（共 180 aa）
参考结构：PDB 2RE9（apo 同源三聚体，2.1 Å，链 A/B/C；可观测残基 88–251）
重要：2RE9 的作者编号与下面给出的构建体编号一致，hotspot 可直接使用，无需偏移。
生物学形态：天然为三聚体。受体 DcR3 / DR3 的结合位点横跨两个相邻 protomer 之间的沟槽。

靶点序列（180 aa）：
LKGQEFAPSHQQVYAPLRADGDKPRAHLTVVRQTPTQHFKNQFPALHWEHELGLAFTKNRMNYTNKFLLIPESGDYFIYSQVTFRGMTSECSEIRQAGRPNKPDSITVVITKVTDSYPEPTQLLMGTKSVCEVGSNWFQPIYLGAMFSLQEGDKLMVNVSDISLVDYTKEDKTFFGAFLL

== 设计约束（比赛规则，必须严格遵守）==
- 长度：85–100 aa。请注意：比赛允许 40–100 aa，但我明确要求 85–100 aa，不要给我更短的设计。
- 单条多肽链，单个结合模块。
- 仅使用 20 种标准天然氨基酸。
- 不得包含信号肽、His-tag、Flag-tag、Fc 或任何融合片段 —— 只要成熟蛋白本体。
- 必须是从头设计，不得以已知的抗 TL1A 抗体或 DcR3 序列为起点做亲和力优化。

== 两组表位策略，请分别设计 ==

A 组（主力，约 70% 的设计量）—— 沟槽表位，靶向两个 protomer 之间：
   hotspot：D104、S105、Y117、E119、P120、T121、N136、W137、F138、K169、E170
   这是 DcR3/DR3 的结合位点，有临床验证价值。
   请以三聚体（至少两条相邻链）作为设计输入，否则无法正确建模这个界面。

B 组（对冲，约 30% 的设计量）—— 单 protomer 表面表位：
   仅靶向完全落在单一 protomer 上的表面，不跨越 protomer 边界。
   请避开上面 A 组中位于链间界面的残基：S105、T121、N136、W137、F138。
   目的：这组设计在"仅用单链靶点"的评分配置下也能被正确评估。

== 产出要求 ==
1. 请先生成大量 backbone（在积分允许范围内尽可能多，目标 500 个以上），
   再用 ProteinMPNN 类工具为每个 backbone 生成多条序列。
2. 用平台的筛选指标做初筛，然后给我排名最高的约 100 条序列。
3. 对这 100 条，请用 ESMFold2 + ipSAE 评分，配置严格为：
   1 seed、1 sample、无模板、仅靶点链使用 MSA、binder 链不使用 MSA。
4. 导出给我：
   - CSV：序列、长度、所属组别（A/B）、ipSAE、binder 的 pLDDT、iptm、以及平台的其他筛选分数
   - 每条设计的复合物 PDB 文件
   - 如果可以，PAE 矩阵（我需要在本地独立复算 ipSAE 做交叉验证）
5. 请明确告诉我：这次运行中 TL1A 是按单体还是三聚体折叠评分的。

== 参考基准 ==
另一个平台（AsterFire）公布过 TL1A minibinder 的结果：ipSAE 0.637 和 0.628，
binder pLDDT 91.8，被他们描述为"可信界面"。我的目标是 ipSAE ≥ 0.75，请以此为筛选门槛。

请先告诉我你的执行计划和预估积分消耗，等我确认后再开始运行。
```

**Giải thích các lựa chọn:**

- **Chia A/B 70-30**: epitope "đúng" về sinh học của TL1A nằm ở groove giữa hai protomer. Nhưng nếu
  cửa sàng chỉ fold **một** chuỗi target, binder nhắm groove sẽ bị chấm ipSAE thấp vì một nửa
  interface của nó không tồn tại trong model. Nhóm B là bảo hiểm cho khả năng đó. Đây là hedge, không
  phải do dự — ta chưa xác nhận được cửa sàng fold monomer hay trimer, và notebook 02 đang đo việc này.
- **Đòi ≥500 backbone**: tỉ lệ thành công của de novo binder design thấp. Muốn chọn ra 20 cái tốt thì
  phải sinh hàng trăm rồi lọc. Đừng để nó đưa 20 cái rồi dừng.
- **Đòi PDB + PAE matrix**: để tôi **tự tính lại ipSAE** bằng `DunbrackLab/IPSAE` ở local. Không tin
  con số platform đưa ra mà không đối chiếu — nếu ipSAE của họ khớp của tôi, ta biết đã tái dựng đúng
  cửa sàng; nếu lệch, ta biết có vấn đề *trước* khi submit.
- **Câu cuối "cho tôi kế hoạch trước khi chạy"**: chống việc nó đốt hết 100k điểm trong một lượt.

---

## Sau khi có kết quả

Gửi tôi CSV. Tôi sẽ:
1. Chạy `tools/novelty_check.py --target tl1a` trên toàn bộ (dù TL1A sạch, vẫn phải kiểm — xem
   `AGY-01-RESULT.md`, lập luận độ dài không đủ đảm bảo)
2. Tự tính lại ipSAE từ PDB + PAE, đối chiếu với số của platform
3. Chọn 20 sequence cuối, cân giữa nhóm A và B theo kết quả monomer-vs-trimer của notebook 02
4. Xuất CSV submit bằng `tools/build_submission.py` (AGY-02)
