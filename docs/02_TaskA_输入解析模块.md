# Task A｜Scaffold 文件解析与基础数据结构

状态：**已完成**。不包含任何 oligo 设计逻辑。

## 已实现

| 能力 | 入口 |
|---|---|
| 读 SnapGene `.dna`（Biopython 优先，内置解析器兜底，两者交叉校验） | `ssdna.io.load_snapgene` |
| 读 FASTA（多记录必须显式指定，不自动取第一条） | `ssdna.io.load_fasta` |
| 数据模型 | `ssdna.models.Scaffold / Feature / Primer / Segment` |
| 序列工具：反向互补、规范化、环状切片、检查 | `ssdna.sequence` |
| 导出 FASTA / 注释 CSV / JSON | `ssdna.report.export_all` |
| 命令行 | `python -m ssdna.cli "minimalplasmid(1).dna" --outdir outputs` |

## 关键数据结构

```
Scaffold(id, sequence, length, topology, features, primers,
         molecule_type, source, sha256, parser_backend, notes,
         sequence_class, cross_check, warnings)
Feature(feature_type, label, strand, segments[], qualifiers, kind, origin)
Segment(start, end, strand, wraps, raw, origin)     # 0-based 半开；跨原点 wraps=True
Primer(name, sequence, site[], site_origin)          # site 恒为 derived
```

三条不可越界的约定：

1. `length != len(sequence)` 直接抛错，不做修正；
2. `topology` 只接受 `circular` / `linear`，其它一律拒绝；
3. `origin` 字段区分 **file（文件里有）** 与 **derived（推断出来的）**，两者绝不互相覆盖。

## 算法要点

- **容器解析**：`[type:1][len:4 BE][data]` 逐块读取，未知块保留不丢弃；序列块首字节是标志位（bit0 = 环状），长度由块长度决定。
- **坐标换算**：SnapGene 的 `1-based 闭区间` → `0-based 半开` 只在 `parse_feature_ranges` 里换算一次，原始字符串留在 `Segment.raw`。
- **跨原点**：`get_circular_subsequence` 用 `seq[start:]+seq[:end-n]` 回绕，起点越界/超长直接报错；不静默取模。
- **反向互补**：IUPAC 全表（A↔T, C↔G, U→A, R↔Y, S↔S, W↔W, K↔M, B↔V, D↔H, N↔N）。
- **双解析器交叉校验**：Biopython 给序列/拓扑/特征，内置解析器给 `Notes` 块与引物原始序列；两者序列必须完全一致（否则报错），特征集合不一致则显式警告并原样列在 `cross_check` 里。

## 测试结果

`python -m unittest discover -s tests` → **65 个测试全部通过**（Task A 新增 39 个，既有 26 个）。

覆盖：

- 线性区间提取、环状末端跨越（含 0..len 全组合穷举比对 `seq*2`）、零长度、越界拒绝；
- 反向互补的碱基、回文、对合性（两次还原）、IUPAC 简并码；
- 非法序列（`X`、`?`、`-`、空串、纯数字）必须抛错，简并碱基合法但要在检查里提示；
- SnapGene 读取：长度 1200、circular、5 条注释坐标、跨原点 MCS、5 条引物、Notes；
- 两个后端结果一致（序列/拓扑/特征坐标），R10 跨原点引物长度 120 nt；
- 三个导出文件的字段、编码（CSV 用 `utf-8-sig`）、坐标约定、sha256。

## 本文件读取到的事实

| 项 | 值 |
|---|---|
| 长度 / 拓扑 | 1200 nt / circular |
| GC | 57.83% |
| 注释 | ori `[78,667)`−、BleoR `[691,1066)`−、MCS `[0,37)` 与 `[1181,1200)`（跨原点拆两段）、cexuf1 `[1146,1166)` |
| 引物 | 测序-f2、F1、F2、R1、R10（后四条 120 nt；R10 与 F1 跨原点） |
| Notes | `Type=Synthetic`、`ConfirmedExperimentally=0`、`SequenceClass=UNA` |
| 交叉校验 | sequence OK / topology OK / features OK |

## 尚未解决

1. **Biopython 在本机装不上**：pip 的解包临时目录被 DSH 沙箱的 ACL 挡住（`Errno 13`）。已把 wheel 直接解到 `build/pylibs`（gitignored），运行时加 `PYTHONPATH` 即可；代码在缺 Biopython 时会自动回退并在输出里说明。用户自己机器上 `pip install -r requirements.txt` 即可。
2. **scaffold 的有效区域与链别仍未定**（`scaffold_region` / `scaffold_strand` 都是 `null`）：本模块只报告文件里写了什么，不做判断。
3. **引物结合位点是推出来的**：Biopython 按引物序列回贴到模板算出位置，`site_origin="derived"`；内置解析器不回溯，因此后端的引物只有序列没有位点。

## 附：输入文件被重新保存过一次（2026-11-19）

`minimalplasmid(1).dna` 在 **2026-11-19 12:20** 被外部重新保存过：

| 项 | 变化 |
|---|---|
| 文件字节数 | 22 716 → 22 586（−130） |
| 文件 sha256 | `5f09e2f4…` → `eef30e74…` |
| 序列长度 / 拓扑 | 1200 nt / circular（**未变**） |
| 5 条注释坐标 | ori `[78,667)`、BleoR `[691,1066)`、MCS `[0,37)`+`[1181,1200)`、cexuf1 `[1146,1166)`（**未变**） |
| 5 条引物结合位点 | 由序列回贴推得：[534:554]、[66:186]、join{[1086:1200],[0:6]}、[6:126]、join{[1146:1200],[0:66]}（**未变**） |
| 容器块 | type 2/3/5 被重写（+214 / +66），type 27 历史块与 type 14 消失 |

结论：**序列没有变，变的是 SnapGene 重新保存时重写的容器块**。

因此代码与参数做了两处加固：

- `Scaffold.sequence_sha256`（序列本身的 hash）成为身份标识，`Scaffold.sha256`（文件字节 hash）只作记录；
- `params/project.params.json` 里同时记 `scaffold_file_sha256` 与 `scaffold_sequence_sha256`，
  校验器与测试都用**序列** hash —— 以后文件被重新保存不会再误报，而序列真被改动时一定会报。

