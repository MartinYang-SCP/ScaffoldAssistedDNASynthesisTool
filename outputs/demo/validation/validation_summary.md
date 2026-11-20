# 验证摘要

- 结论：**PASS**
- 验证目录：`D:\THU\大四上\核酸纳米材料\小课题\outputs\demo`
- 拓扑：circular
- oligo 6 条（其中跨界 2 条）、staple 28 条、nick 6 条
- 错误 0、警告 0、提示 0

## 目标覆盖

| 项 | 值 |
|---|---|
| target_length | 728 |
| covered_bases | 728 |
| uncovered_bases | 0 |
| overlapping_bases | 0 |
| coverage_fraction | 1.0 |
| topology | circular |

## 各项检查

| 检查 | 状态 | 问题数 |
|---|---|---|
| sequence_basics | pass | 0 |
| target_coverage | pass | 0 |
| group_order | pass | 0 |
| boundary_spanning | pass | 0 |
| staple_length | pass | 0 |
| scaffold_mapping | pass | 0 |
| target_binding | pass | 0 |
| nicks | pass | 0 |

## 问题清单

没有发现问题。

## 热力学工具（未运行）

**idt_oligoanalyzer**：`not_run`，输入文件 `idt_oligoanalyzer_input.csv`

- 能查：GC content, Tm, hairpin, self-dimer, hetero-dimer
- 注意：这些工具只算二级结构与杂交倾向，不能验证连接化学，也不能预测支架组装效率

**nupack**：`not_run`，输入文件 `nupack_strands.txt`

- 能查：多链竞争性杂交, 替代配对, 复合物平衡浓度
- 注意：NUPACK 不模拟 DNA 连接酶，也不模拟 scaffold 的拓扑约束

> 本模块没有运行上述任一工具，报告里没有它们的结果
