# MATVENUS-03 — TL1A design campaign (bản v2, viết lại sau khi có bài học IL-7Rα)

> Dán vào https://matvenus.com/chat/agent
> Thay thế PROMPT 2 của `MATVENUS-01`. Lý do thay: `resp_01.md` + `rep02.md` đã lộ ra
> 6 giới hạn thật của platform và 1 sai số đo. Prompt cũ yêu cầu nhiều thứ platform không làm được.

## Những gì đã thay đổi so với brief cũ, và vì sao

| Chỗ sửa | Cũ | Mới | Lý do |
|---|---|---|---|
| Độ dài | 85–100 aa | **60–95 aa** | Sàn 85 aa được dựng ra chỉ để né luật novelty của IL-7Rα. TL1A có **0** de novo minibinder đã công bố nên không có gì để né. Mavenus đã chỉ ra `d0` của ipSAE phụ thuộc **số residue ở interface**, không phụ thuộc tổng chiều dài; khối lượng không tiếp xúc còn làm loãng pTM. Bộ minibinder de novo thành công đã công bố phần lớn nằm ở 50–70 aa. |
| Tỉ lệ A/B | 70 / 30 | **85 / 15** | Nhóm B là bảo hiểm cho khả năng không fold được trimer. Đã đo được: trimer + binder fold ở **13.7–14.0 GB**, 121–136 s trên L4 21.5 GB. Bảo hiểm đó không còn cần nhiều. |
| Đòi ipSAE từ platform | có | **bỏ** | Platform không có ipSAE native. Ta tự tính trên Colab. |
| Đòi cấu hình MSA | có | **bỏ** | Platform **không hề** có tham số MSA/template, càng không có MSA per-chain. Không đòi thứ không tồn tại. |
| Đòi PAE matrix | "nếu được" | **bắt buộc** | Đã xác nhận: 38/38 file trong `metrics_json/` của bộ IL-7Rα **đều có** `pae_matrix` + `pae_tokens`. Platform làm được việc này. |
| Chứng dương | không có | **DcR3 + TL1A** | Lặp lại cách `ctrl_7S5B` đã dùng cho IL-7Rα. Không có vạch tham chiếu cùng pipeline thì ngưỡng 0.70 vô nghĩa. |
| Topology nhóm B | đòi α/β pha | **bỏ** | Platform không ràng buộc được fold class. Bộ IL-7Rα chỉ ra 3/40. Đòi là đốt sampling. |

---

## PROMPT — dán nguyên khối

```
任务：为 TL1A 从头设计 minibinder。这是第二个赛道；IL-7Rα 赛道我们已经用你的平台跑完一轮，
下面我会直接写明从那一轮学到的东西，请不要重复那些弯路。

== 从 IL-7Rα 那一轮已经确认的事实，请直接采用，不要重新探索 ==
1. RFdiffusion 输出的复合物 PDB 只有主链原子（N/CA/C/O），没有侧链。
   因此几何 QC 必须用 CA–CA 判据，不能用侧链重原子 5 Å 阈值 —— 上一轮用错判据
   导致 12 个骨架被误判为"全部没有界面"，后来修正后结论完全反转。
2. 平台无法约束二级结构或 fold class。不要尝试指定 4 螺旋束或 α/β 混合拓扑，
   产出率只有 3/40。接受螺旋束为主的结果。
3. 平台没有 MSA / template 参数，也没有原生 ipSAE。这两点我已经接受，不要为此设计替代方案。
   我需要的是原始数据（PDB + PAE 矩阵），ipSAE 由我在本地用 DunbrackLab/IPSAE 官方脚本计算。
4. 你的 ESMFold2 导出格式是可用的：IL-7Rα 那一轮每个设计的 metrics JSON 里都含
   pae_matrix（N×N）和 pae_tokens（含 chain_id、residue_index，顺序与矩阵行列一致），
   以及 ptm、iptm、plddt。本轮请用完全相同的格式导出。
5. 积分：你的任务返回体里没有计费字段，所以不要估算积分。只报告 estimated_seconds，
   实际扣费我自己去账户页核对。

== 靶点 ==
名称：TL1A / TNFSF15
UniProt：O95150，残基范围 Leu72–Leu251（共 180 aa）
参考结构：PDB 2RE9（apo 同源三聚体，2.1 Å，链 A/B/C）
生物学形态：天然三聚体。受体 DcR3 / DR3 的结合位点横跨两个相邻 protomer 之间的沟槽。

靶点序列（180 aa，这是我的构建体编号 1–180 对应 Leu72–Leu251）：
LKGQEFAPSHQQVYAPLRADGDKPRAHLTVVRQTPTQHFKNQFPALHWEHELGLAFTKNRMNYTNKFLLIPESGDYFIYSQVTFRGMTSECSEIRQAGRPNKPDSITVVITKVTDSYPEPTQLLMGTKSVCEVGSNWFQPIYLGAMFSLQEGDKLMVNVSDISLVDYTKEDKTFFGAFLL

== 第一步，先做编号核验（这一步必须先做，不要跳过）==
上一轮在 IL-7Rα 上，7OPB 的编号比我的构建体编号 +20，直接套用 hotspot 会静默落到错误残基上。
这次请先逐位比对 2RE9 链 A 的可观测序列与上面这条 180 aa 序列，然后明确告诉我：
  a) 2RE9 的作者编号与我的构建体编号 1–180 的换算关系（偏移多少，或是否一致）
  b) 2RE9 中哪些残基是缺失/不可观测的
  c) 下面我给的 hotspot 列表，在 2RE9 的实际编号里应该写成什么
在你给出这三个答案并经我确认之前，不要提交任何生成任务。

== 表位策略 ==
A 组（主力，约 85% 的设计量）—— 沟槽表位，横跨两个相邻 protomer：
   hotspot（我的构建体编号）：D104、S105、Y117、E119、P120、T121、N136、W137、F138、K169、E170
   这是 DcR3/DR3 的结合位点，有临床验证价值。
   设计输入必须包含至少两条相邻链（建议整个三聚体 A/B/C），否则这个界面不存在。

B 组（对冲，约 15%）—— 单 protomer 表面表位：
   只靶向完全落在单一 protomer 上的表面。请避开上面那些位于链间界面的残基。

== 设计约束（比赛硬规则 + 我的选择）==
- 长度：60–95 aa。比赛允许 40–100，我要 60–95。请把这个区间作为硬约束。
- 单条多肽链，单个结合模块。
- 仅 20 种标准天然氨基酸。
- 不含半胱氨酸（Cys = 0）。
- 不得含信号肽、His-tag、Flag-tag、Fc 或任何融合片段。
- 必须从头设计，不得以已知抗 TL1A 抗体或 DcR3 序列为起点做优化。

== 需要避开的位点（平台只有正向 hotspot，所以这些靠事后筛）==
- 糖基化位点（我的构建体编号）：N62、N158 —— 设计需距其侧链 ≥8 Å
- 链内二硫键：C91、C131 —— 不要把界面压在上面
请用 CA–CA 口径先筛，并在报告里注明这是 CA 口径（侧链实际更远）。

== 分阶段执行，每阶段有 gate ==
S1 试点：RFdiffusion 150 个 backbone（A 组约 128 / B 组约 22），长度 60–95。
        然后做逐设计几何 QC（CA–CA 判据）：CA_min ≥4.2 Å，且 ≥3 个 hotspot 的 CA ≤10 Å，
        并统计受体侧接触残基数、hotspot 命中数、与 N62/N158 的距离。
        Gate：给我 QC 表，我确认后再进入 S2。不要自动放量。
S2 序列：对通过 QC 且按「界面规模 × hotspot 覆盖」排序的 top 30 骨架，
        用 ProteinMPNN 生成 6 条/骨架，去冗余，剔除含 Cys 的序列。
S3 折叠：对全部存活序列跑 ESMFold2 复合物预测，seed=1、samples=1。
        明确告诉我：TL1A 是按三聚体还是单体折叠的。我要的是三聚体。

== 必须一起跑的阳性对照（很重要）==
请用与 S3 完全相同的 pipeline 和参数，额外跑一个对照：
   TL1A 三聚体 + DcR3 的 TL1A 结合域（DcR3 = UniProt O95407，TNFRSF6B）
DcR3 是 TL1A 的天然高亲和力受体，是真实存在的结合物。
我需要它在同一条 pipeline 下的 ipTM / PAE，作为判断我的设计好坏的参照线。
没有这条对照，任何绝对数值都没有意义。请把它命名为 ctrl_DcR3 一起导出。

== 交付物 ==
1. CSV：序列、长度、组别(A/B)、pLDDT、pTM、ipTM、界面残基数、hotspot 命中列表、
   与 N62/N158 的距离、是否含 Cys
2. 每个设计的复合物结构文件（PDB 或 CIF 均可）
3. 每个设计的 metrics JSON，格式与 IL-7Rα 那一轮一致：必须含 pae_matrix 与 pae_tokens
4. ctrl_DcR3 的同样三件套
5. 明确写出 S3 用的是三聚体还是单体

请先只回答编号核验那三个问题（a/b/c）和你的 S1 执行计划，等我确认后再启动任何生成任务。
```

## Sau khi có kết quả, phần của tôi

1. Tự tính ipSAE từ `pae_matrix` + structure bằng `ipsae.py` chính thức, cấu hình `pae_cutoff=10, dist_cutoff=10`
2. Fold lại trên Colab bằng **đúng cấu hình cửa sàng** (MSA chỉ cho target, binder không MSA),
   vì platform không làm được bước này — đây mới là con số quyết định suất wet-lab
3. So với vạch `ctrl_DcR3` để biết thiết kế của ta đứng đâu so với binder thật
4. `tools/novelty_check.py --target tl1a`
5. Chọn 20, xuất bằng `tools/build_submission.py`
