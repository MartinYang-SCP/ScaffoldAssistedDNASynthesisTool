# 环状 scaffold + staple 的目标片段组装设计工具

用一条约 1.2 kb 的**环状单链** scaffold，靠统一设计的 39 nt staple（13 nt 贴 scaffold + 26 nt 贴目标），
把由多条 oligo 拼成的**目标片段**定位并闭合成环。

**这个程序验证的是序列设计与映射关系。** 它不能证明实际实验中所有 nick 都会被连接，
不能证明 DNA 在酵母中一定稳定存在或能表达，也不能替代任何湿实验。
最终能否成立，需要结合实验验证结果进一步判断。

---

## ⚠️ 先说清楚它做什么、不做什么

| 程序能做的 | 程序**不能**做的 |
|---|---|
| 解析 SnapGene `.dna`，检查序列/坐标/拓扑 | 证明连接酶会把所有 nick 连上 |
| 生成 target oligo 与 staple 的 5′→3′ 序列，并逐条验证互补关系 | 证明这个环在酵母里能稳定复制、能表达 |
| 检查覆盖、分组、交界、nick、scaffold 映射是否自洽 | 预测二级结构、组装产率、脱靶杂交（需 IDT/NUPACK 等外部工具） |
| 给出 Tm 分布与 scaffold 相位建议 | 用 Tm 反推"应该合成哪条序列"（序列没有自由度） |

**必须知道的四件事**：

1. 程序验证的是**序列设计与映射关系**；
2. 程序**不能证明**实际实验中所有 nick 都会被连接；
3. 程序**不能证明** DNA 在酵母中一定稳定存在或能够表达；
4. 需要结合**实验验证结果**进一步判断。

---

## 当前状态

| 阶段 | 内容 | 交付 |
|---|---|---|
| 0 | 结构理解与确认 | [01_结构确认_v2.md](docs/01_结构确认_v2.md) |
| A | 输入解析与数据模型 | [02_TaskA_输入解析模块.md](docs/02_TaskA_输入解析模块.md) |
| B | scaffold–staple–target 结构模型 | [03_TaskB_结构模型与映射.md](docs/03_TaskB_结构模型与映射.md) |
| C | 序列生成引擎 | [04_TaskC_序列生成引擎.md](docs/04_TaskC_序列生成引擎.md) |
| D | 独立验证模块 + 热力学接口 | [05_TaskD_验证模块.md](docs/05_TaskD_验证模块.md) |
| E | 整合 / 测试 / 交付 | 本 README + `main.py` |
| — | 外部评审意见处理 | [07_评审意见处理.md](docs/07_评审意见处理.md) |

**尚未生成任何真实实验序列**：默认配置里的关键参数都是 `null`，程序会在第 1 步停下并列出缺失项（现在是 8 项）。

---

## 安装

```bash
python -m pip install -r requirements.txt   # biopython>=1.81, pyyaml>=6.0
python --version                            # 需要 3.10+（用了 `X | None` 语法）
```

- `biopython` 用于 SnapGene 解析（交叉校验）与 Tm 计算；
- `pyyaml` 用于读配置文件（没有它时可用等价的 `.json` 配置）。

两个依赖都是**可降级**的：缺 Biopython 时自动退回内置解析器（`parser_backend=raw`，Tm 显示"未计算"），
缺 PyYAML 时报错并提示，不会静默用一个猜测的配置。

---

## 快速开始

```bash
# 1) 先看"缺参数就停"的行为（这是默认路径）
python main.py --config config/design_config.yaml
#   → [1/8] 读取设计参数并校验完整性 …
#     缺少 7 个必需参数，流程在这里停止 …（退出码 2，不生成任何文件）

# 2) 跑通完整流程（目标序列是合成序列，不是真实基因）
python main.py --config config/design_config.demo.yaml
#   → [1/8]…[8/8] 逐步打印，最后输出设计状态 OK / 验证状态 PASS

# 3) 只跑某个模块
python -m ssdna.design --demo --scan-phase            # 相位扫描，不生成序列
python -m ssdna.validate outputs/demo                  # 单独验证一个输出目录
python -m unittest discover -s tests -v                # 196 个测试
```

演示输出在 `outputs/demo/`，验证结果在 `outputs/demo/validation/`。

---

## 项目结构

```
.
├── main.py                      ★ 唯一入口：8 步工作流
├── data/
│   ├── scaffold.dna             输入：SnapGene 文件（1200 bp 环状）
│   └── target.fasta             演示用合成目标（728 nt，非真实基因）
├── config/
│   ├── design_config.yaml       模板（关键项为 null → 会停）
│   └── design_config.demo.yaml  演示配置（合成序列）
├── ssdna/                       ★ 核心包
│   ├── models.py                Scaffold / Feature / TargetRegion / TargetOligo / Staple / Nick / Junction
│   ├── sequence.py              reverse_complement / normalize / 环状切片 / 检查
│   ├── io.py                    load_snapgene（Biopython 优先 + 内置兜底 + 交叉校验）/ load_fasta
│   ├── snapgene_raw.py          内置 SnapGene 字节解析器（零依赖）
│   ├── mapping.py               装配层：13 nt 平铺 / oligo 划分 / nick / junction
│   ├── design.py                生成引擎：target oligo + staple + 9 项自检 + 切点建议 + 相位扫描
│   ├── thermo.py                Tm 体检（只做检查，不参与选序列）
│   ├── validate.py              独立验证（只读文件，不 import 生成器）
│   ├── report.py                FASTA / CSV / JSON 导出
│   └── cli.py                   `python -m ssdna.cli` 单文件解析入口
├── tests/                       196 个单元测试
├── outputs/                     运行产物（demo / validation）
├── docs/                        各阶段报告 + 5 张图
├── params/                      参数清单（JSON）与 design_config 旧版模板
└── requirements.txt
```

任务书建议的目录名与本项目实际文件的对应关系（功能一一对应，只是包名用了 `ssdna/` 而不是 `src/`）：

| 建议 | 实际 |
|---|---|
| `src/models.py` | `ssdna/models.py` |
| `src/sequence_utils.py` | `ssdna/sequence.py` |
| `src/scaffold_parser.py` | `ssdna/io.py` + `ssdna/snapgene_raw.py` |
| `src/mapping.py` | `ssdna/mapping.py` |
| `src/target_design.py` / `staple_design.py` | `ssdna/design.py`（同一流程内的两步） |
| `src/validation.py` | `ssdna/validate.py` |
| `src/export.py` | `ssdna/report.py` + `design.write_outputs` |
| `main.py` | `main.py` |
| `tests/test_*.py` | `tests/test_{snapgene_reader,task_a,mapping,design,validate,workflow}.py` |

---

## 工作流（`main.py`，8 步）

```
[1/8] 读取设计参数并校验完整性   ← 缺参数在这里停，退出码 2
[2/8] 读取 scaffold              ← SnapGene 解析 + 交叉校验
[3/8] 读取目标序列               ← FASTA
[4/8] 建立结构映射               ← A/B 分区 + 两个交界
[5/8] 生成 target oligos
[6/8] 生成 staples
[7/8] 导出结果                   ← CSV / FASTA / JSON
[8/8] 独立验证                   ← 只读上一步落盘的文件
```

**有两处顺序与任务书不同，原因写在这里**：

- 任务书把"读参数"放在第 3 步，但第 2、3 步的文件路径都来自配置，所以这里提前到第 1 步；
- 任务书把"验证"放在导出前，但验证模块被要求**独立**（只读文件、不碰生成器内存），
  所以它必须排在导出之后。其余步骤完全按任务书顺序。

任何一步失败都会给出具体原因（哪个 oligo、哪个坐标、什么类型、怎么修），不会静默继续。

---

## 配置文件

`config/design_config.yaml` 是模板，值为 `null` 的项全部是**阻塞项**：

```yaml
scaffold:
  file: data/scaffold.dna
  strand: null            # 实验室的单链对应文件正链(+)还是负链(-)
  region: null            # "whole" 或 [a,b)：文件坐标下的有效区域，必须显式给
staple:
  scaffold_binding_length: 13
  target_binding_length: 26
  total_length: 39        # 必须等于两者之和
  group_pattern: "AB"
  domain_order: "5p_target_26"
  covered_start: 0
target:
  source: null            # FASTA 路径
  group_boundaries: null  # [364] 表示 A=[0,364)、B=[364,728)
  oligo_cuts: null        # nick 位置，不能压 staple 接缝
  topology: circular
  five_prime_phosphate: null
  intended_ligation: null
analysis:
  tm: true
  na_mM: 50.0
  mg_mM: 0.0
  temperature_c: null
output:
  outdir: outputs/design
```

几个约束（程序会检查，不会自动折中）：

- `total_length == scaffold_binding_length + target_binding_length`；
- 每条 target oligo 至少跨过一个 staple 臂接缝，nick 不能压在接缝上；
- 目标长度**不必**是 26 的倍数：26 nt 臂铺不满的部分就是未被锚定的尾巴；
- `group_pattern: "AB"` 时两组 staple 条数只差 1，因此 `L(A) ≈ L(B)`；
- `scaffold.region` 不填就停：程序不替你假设"整段 1200 nt 都是 scaffold"。

**连接状态**：`nick` 的 `intended_ligation` 只是设计意图，`circle_closed` 只有在
`ligation_state == "ligated"`（实验确认）时才为真。程序不会把"计划连接"当成"已经连上"。

---

## 输出

| 文件 | 内容 |
|---|---|
| `target_oligos.csv` | id / 组 / 类型 / 区间 / 长度 / 方向 / 是否跨界 / 端基修饰 / **序列** |
| `staples.csv` | id / 组 / **39 nt 序列** / scaffold 区间 / target 区间 / domain_order / 两个 domain 的方向 |
| `nicks.csv` | 每条相邻关系的 3′/5′ 端、坐标、5′-P 状态、连接意图 |
| `design_config.yaml` | 实际生效的配置 |
| `design_report.json` | 状态、输入 sha256、计数、9 项自检、Tm 汇总 |
| `target.fasta` / `scaffold_working.fasta` | 实际使用的序列（按 strand 定向后的 scaffold） |
| `validation/validation_report.json` | 8 项独立检查、全部 issue、热力学工具状态 |
| `validation/validation_summary.md` | 人读版摘要 |
| `validation/potential_conflicts.csv` | 需要人工关注的冲突 |
| `validation/idt_oligoanalyzer_input.csv` | IDT 批量输入（**未运行**） |
| `validation/nupack_strands.txt` | NUPACK 链列表（**未运行**） |

所有坐标都是 **0-based 左闭右开**；跨原点用 `start > end` 或 `end > length` 表示，两份都写清楚。

---

## 测试

```bash
python -m unittest discover -s tests -v      # 196 个
```

任务书要求的十类测试与对应位置：

| 要求 | 测试 |
|---|---|
| 普通线性序列 | `test_mapping.TestLinearTarget`、`test_design.TestTargetOligos` |
| 环状 scaffold 跨边界 | `test_snapgene_reader.TestCircularAccess`、`test_task_a.TestCircularSubsequence`、`test_mapping.TestCircularScaffoldMapping` |
| reverse complement | `test_task_a.TestReverseComplement`、`test_validate.TestSequenceLevelChecks.test_reverse_complement` |
| A/B 分组 | `test_mapping.TestToyStructure.test_groups_alternate`、`test_validate.TestCoverage.test_group_swapped` |
| 两个 boundary-spanning oligo | `test_mapping.TestToyStructure.test_two_boundary_oligos`、`test_design.TestTargetOligos.test_boundary_oligos`、`test_validate.TestBoundaryOligos` |
| nick 关系 | `test_mapping.TestToyStructure.test_nicks_are_explicit_objects`、`test_validate.TestNicks` |
| staple 长度与 domain | `test_mapping.TestToyStructure.test_domain_lengths`、`test_design.TestStaples` |
| 错误方向 | `test_validate.TestSequenceLevelChecks.test_reverse_complement`（同向互补）、`test_validate.TestTamperedOligo.test_wrong_strand` |
| 不完整覆盖 | `test_validate.TestCoverage.test_dropped_oligo_leaves_gap`、`test_design.TestChecks` |
| 重复覆盖 | `test_validate.TestCoverage.test_overlap_is_detected`、`test_validate.TestStapleChecks.test_duplicate_scaffold_site` |

另外：`test_workflow.py` 端到端跑 `main.py`（含 `subprocess` 真实调用），验证缺参数时退出码 2 且不落盘。

---

## 设计模型摘要

- scaffold 被 13 nt 一块地铺满，相邻两块交替属于 A 组与 B 组 staple（13+13 = 26）；
- 每条 staple = `5′-[revcomp(目标 26)]-[revcomp(scaffold 13)]-3′`，两个 domain 反向平行，
  折返只影响空间走向，**不改变链的真实 5′→3′ 序列**；
- 目标是 A、B 两段，由两个交界接成环；交界**不是 nick**，它必须被一条 boundary-spanning oligo
  连续覆盖——"序列连续"不等于"已经共价闭环"；
- 唯一能用 Tm 做决策的地方是 scaffold 相位（`--scan-phase`），不是改序列。

---

## 假设与限制（请务必读）

1. **scaffold 的链别与有效区域仍未确定**（`scaffold.strand` 为 null）。选错会让所有 13 nt 臂反向作废。
2. **目标序列未定**。演示用的是合成序列；真实设计要满足"A/B 交界各被一条 oligo 覆盖"和"nick 避开接缝"。
3. **A/B 长度必须接近**（`|n_A − n_B| ≤ 1`），因为 13 nt 块严格交替。若真实需求是 2:1，需要改分组模式，
   那时"13+13 成对"只在局部成立，要重新定义铺法。
4. **折返处是一个 180° 的骨架转角**，其能量代价本模块不评估。
5. **13 nt 臂单独结合的 Tm 偏低**（演示里 29–56 °C），结构稳定性靠几十条 staple 的协同；
   本模块只报 Tm 分布，不判断"会不会散"。
6. **热力学工具未接**：IDT / NUPACK 只生成了输入文件，状态是 `not_run`，报告里没有任何伪造结果。
   这两个工具也**不能**验证连接化学。
7. **不模拟**：连接酶效率、酵母体内修复、复制与表达、脱靶杂交、二级结构。
8. **Biopython / PyYAML 在本机是通过 `build/pylibs` 提供的**（沙箱 pip 装不上）；
   代码里两条路径都测过，但你自己的环境请正常 `pip install -r requirements.txt`。

---

## 常见问题

**Q：为什么默认配置一跑就退出？**
A：设计参数没填全。这是有意的——见 README 顶部第 1 条，程序不会替你猜。

**Q：目标序列长度必须是 26 的倍数吗？**
A：不必。26 nt 臂按 26 的节距铺，铺不满的尾巴就是未被锚定的部分（可以留作交界 linker）。

**Q：能不能用 Tm 来挑 staple 序列？**
A：不能。staple 序列由 scaffold 与目标唯一确定，没有自由度。Tm 只能用来体检、选相位、选 nick 位置。
