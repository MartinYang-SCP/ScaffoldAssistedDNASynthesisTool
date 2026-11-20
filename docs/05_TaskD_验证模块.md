# Task D｜独立验证模块与热力学接口

## 0. 先修正一条规则：目标长度**不需要**是 26 的倍数

上一轮我把一个**设计选择**当成了物理定律，你的纠正是对的。真正必须成立的是这三条：

| 真正的约束 | 说明 |
|---|---|
| 每条 staple 的 26 nt 臂必须落在目标上，且**互不重叠** | 两条 staple 抢同一段目标才算冲突 |
| nick 不能压在 staple 臂的接缝上 | 接缝 = 26 nt 臂的端点，不是"26 的倍数"这个数本身 |
| 每条 oligo 至少跨过一个接缝 | 否则它只被一条 staple 锚定，容易松 |

于是：

- **目标总长 L 自由**，可以是 720、723、715 这种数；
- 只有"被 26 nt 臂覆盖的那部分"是 26 的整数倍，剩下的尾巴（≤ 25 nt）不锚定，**留作交界 linker 正好**；
- 接缝只决定"不能切在哪"，不决定"必须切几段"。

代码上做了三处改动：

1. `design_staples` 里 `group_length` 从 **error 降为 note**，并明确写出"末尾 N nt 不被臂覆盖"；
2. 新增 `suggest_oligo_cuts(length, seams=..., min_length=40, max_length=120, clearance=3)`，
   按上面三条约束给出切点建议（只建议，不改目标序列）；
3. 演示里实测：728 / 720 / 723 / 715 / 600 / 900 nt 都能切出合法分段，长度 100–120 nt。

```
728 → [108, 104, 104, 104, 104, 104, 100]   7 条
720 → [100, 104, 104, 103, 103, 103, 103]   7 条   ← 不是 26 的倍数，照样切
600 → [120, 120, 120, 120, 120]             5 条
```

顺带说明：目标 oligo 的**序列**本来就是从目标序列切片得到的（`target[interval]`，环状自动回绕），
没有任何重新生成——切点只是决定"切在哪儿"。

## 1. 独立性

`ssdna/validate.py` **不 import** `ssdna.design`，只读目录里的文件：

```
target.fasta / scaffold_working.fasta / target_oligos.csv / staples.csv / nicks.csv / design_report.json
```

换一套别人产出的同结构目录也能验。有测试断言源码里没有 `from .design` / `import design`，
以及删掉 `design_report.json` 后仍能跑（只是交界信息缺失，会报边界问题）。

## 2. 八个检查函数

```python
validate_sequence(sequence)
validate_length(sequence, expected_length)
validate_reverse_complement(seq1, seq2)
validate_target_coverage(target_sequence, oligos)      # → (issues, stats)
validate_group_order(oligos)
validate_boundary_spanning_oligos(oligos, boundaries)
validate_scaffold_mapping(scaffold, staples)
validate_nicks(oligos, nick_records)
```

每条 `Issue` 都带：`level`（error/warning/info）、`code`、`subject`（oligo id）、
`coordinates`、`message`、`suggestion`（修复建议）。测试里逐项注入错误验证它真的能抓到。

### 抓得到的错误（都有对应测试）

| code | 触发方式 |
|---|---|
| `oligo_sequence_mismatch` | 改一个碱基 |
| `uncovered_bases` / `overlapping_bases` | 删一条 oligo / 把起点往前挪 |
| `group_order` | A/B 标签互换 |
| `boundary_count` / `junction_mislabeled` | 把跨界 oligo 改成 ordinary |
| `not_reverse_complement` / `same_orientation_complement` | 篡改 staple 或写成同向互补 |
| `scaffold_site_collision` / `scaffold_site_overlap` | 两条 staple 抢同一段 scaffold |
| `target_site_overlap` | 两条 26 nt 臂重叠 |
| `nick_count` / `nick_not_recorded` / `nick_wrong_end` | 删一条 nick / 把 3′ 写成 5′ |
| `ligation_without_phosphate` | 打算连接但 5′-OH |

## 3. Nick 检查的两个要点

- **逐条列**相邻关系：上游必须是 3′、下游必须是 5′；坐标要与上游 `target_end` 对得上；
- **相邻 ≠ 已连接**：`intended_ligation` 只是意图，报告里从不把它当成"已经连上"；
  5′-P 未知时给 `info`，不是 warning——因为合成 oligo 默认就是 5′-OH，这是常态而非错误。

## 4. 热力学接口：只生成输入，绝不编结果

IDT 与 NUPACK 在本环境都不可调用，所以：

| 工具 | 产物 | 状态 |
|---|---|---|
| IDT OligoAnalyzer | `idt_oligoanalyzer_input.csv` | `not_run` |
| NUPACK | `nupack_strands.txt`（`name sequence` 两列） | `not_run` |

报告里写清"能查什么"（GC/Tm/hairpin/self-dimer/hetero-dimer；多链竞争性杂交）和
"不能查什么"：**两者都不能验证连接化学**，NUPACK 也不模拟连接酶与 scaffold 的拓扑约束。
有测试断言报告里不出现任何伪造的 `ΔG` / `hairpin_score` 字段。

## 5. 产物

- `validation_report.json`：状态、输入文件 sha256、计数、覆盖统计、8 项检查、全部 issue、thermo 状态
- `validation_summary.md`：人读版（表格 + 问题清单 + 工具运行说明）
- `potential_conflicts.csv`：所有需要人工关注的冲突（error + warning）

演示跑出来：**PASS**，0 error / 0 warning，覆盖 728/728 nt，6 条 oligo（2 条跨界）、28 条 staple、6 条 nick。

## 6. 测试

`python -m unittest discover -s tests` → **166 个测试全部通过**（Task D 新增 30 个，切点建议新增 3 个）。

## 7. 尚未解决

1. `scaffold.strand` 仍未定，所有 13 nt 臂的方向都悬着。
2. 热力学工具没装：要走 IDT / NUPACK 得手工跑，报告里已给出步骤。
3. 验证器目前只做"结构/序列/化学标注"层面的检查，不做二级结构预测（那是工具的事）。
