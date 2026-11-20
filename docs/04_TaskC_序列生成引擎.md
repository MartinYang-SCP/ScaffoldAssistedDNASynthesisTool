# Task C｜Target oligo 与 staple 序列生成引擎

## 0. 先回答：要不要算 Tm、能不能据此设计

**要算，但不能"据此设计"——因为这套结构里序列没有自由度。**

staple 的 13 nt 必须是 scaffold 那 13 nt 的反向互补，26 nt 必须是目标那 26 nt 的反向互补。
这两个约束是硬的，没有任何"换一段序列把 Tm 调到 60 °C"的余地。所以 Tm 在这里只有两个用途：

| 用途 | 是否合法 | 说明 |
|---|---|---|
| 体检：找出 13 nt 臂的 Tm 离群者 | ✅ | 同一批 staple 里 Tm 跨 27 °C，说明锚定强度不均，组装会不同步 |
| 决策：**扫描 scaffold 相位** | ✅ | 唯一真正能"用 Tm 做设计"的地方：把目标放在 scaffold 的哪一段 |
| 决策：**选 nick 位置** | ✅ | 避免把 nick 放在弱臂旁边 |
| 优化：改 staple 序列让 Tm 更均匀 | ❌ | 不可能。改了就与 scaffold/目标不互补 |
| 优化：改目标序列 | ❌ | 违反"不为了满足长度/性质改动目标序列" |

演示（728 nt 合成目标、28 条 staple、50 mM Na⁺）：

```
13 nt 臂：n=28  mean=44.1  范围 [29.3, 56.0]  stdev=6.3
26 nt 臂：n=28  mean=61.2  范围 [51.0, 69.0]  stdev=4.8
```

13-mer 的跨度接近 27 °C，这不是"要不要优化"的问题，而是"这批 staple 只能靠协同作用"的定量证据。
而它**无法通过改序列解决**——只能靠换位置。

所以我把"用 Tm 做设计"做成了一个独立的、只读的操作：`scan_scaffold_phase()`。

```
$ python -m ssdna.design --demo --scan-phase
相位扫描（93 个位置，按 13 nt 臂 Tm 标准差升序）：
  covered_start= 689 stdev= 3.81 mean=49.33 range=[44.64, 58.63]
  covered_start= 702 stdev= 4.04 mean=49.18 range=[41.24, 58.63]
  covered_start= 715 stdev= 4.35 mean=48.85 range=[40.42, 58.63]
  ...
```

同一个设计，只把覆盖区从 0 挪到 689，13 nt 臂的 Tm 标准差就从 6.29 降到 3.81，
**一条序列都没改**。这才是 Tm 在本项目里的正确用法。

Tm 用 SantaLucia 最近邻（Biopython `Tm_NN`），未做 Mg²⁺ 完整修正，因此只用于**同批横向比较**，
不要拿绝对值当结论。没有 Biopython 时返回 `None` 并明说"没算"，不是 0。

## 1. 已实现的接口

```python
design_target_oligos(target_sequence, group_definitions, oligo_parameters,
                     boundary_definitions) -> (regions, oligos, checks)

design_staples(scaffold, target_oligos, staple_parameters,
               scaffold_strand="+", target_sequence=...) -> (staples, checks)

check_design(...) -> list[Check]      # 9 项生成后检查
scan_scaffold_phase(...) -> list[dict]  # 相位扫描（只读）
run_design(inputs, outdir=...) -> dict   # 全流程 + 四个产物
find_missing(config) -> list[str]        # 缺参数清单
```

## 2. Target oligo 生成

- 切点由用户给（`oligo_cuts`），**不自动等分**；
- 每条 oligo 的序列 = 目标序列在该区间的**原样切片**（环状自动回绕），不补不裁；
- `oligo_type` 由"区间是否严格包含某个交界坐标"推出，不是手工标；
- `target_segments` 记录它落在哪些 TargetRegion（跨界的会有两个）；
- 跨界位置记录在该区间 + `junctions` 表里，两条交界分别处理（环状 2 个、线性 1 个）。

演示里 `T03 = [300,420)` 跨 J1(364)，`T06 = [660,828)` 跨 J2(0)，各自 120/168 nt，
自身序列连续，没有"互补桥接"这种简化。

## 3. Staple 生成

```
scaffold 第 j 个 13 nt 块： [covered_start + C - 13j - 13, covered_start + C - 13j)
staple 序列：              5′-[revcomp(目标 26)]-[revcomp(scaffold 13)]-3′
```

- 长度从配置读，**没有写死**；`13/26/39` 是默认值，改配置即可（`total_length ≠ 13+26` 会报错）；
- 两个 domain 的方向都由反向互补决定（都是 `-`）；折返只影响空间走向，**不改变链的真实 5′→3′ 序列**；
- 每条 staple 都能追溯到 scaffold 坐标 + 目标坐标 + 两个 domain 的真实序列。

## 4. 配置

`params/design_config.yaml` 是模板，空值全部是 `null`（阻塞）。

```yaml
staple: {scaffold_binding_length: 13, target_binding_length: 26, total_length: 39, ...}
scaffold: {file: ..., strand: null, topology: circular}
target: {source: null, groups: [A, B], boundary_count: 2, oligo_cuts: null, ...}
```

`null` 的项会出现在缺失清单里：**填完之前不生成任何序列，也不创建输出目录**（已有测试断言）。

## 5. 生成后检查（9 项）

| code | 检查 |
|---|---|
| `staple_scaffold_length` / `staple_target_length` / `staple_total_length` | 13 / 26 / 39 |
| `not_complementary` | 逐条**重算** revcomp，不信构造过程 |
| `target_coverage` | oligo 总长 = 目标长 |
| `group_mismatch` | A/B 归属与起点所在区域一致 |
| `junction_coverage` | 每个交界被恰好一条 oligo 覆盖 |
| `scaffold_bounds` | 13 nt 块不越界、回绕正常 |
| `unassigned_target` | 没有未被覆盖的碱基 |
| `duplicate_scaffold_site` / `duplicate_target_site` | 没有两条 staple 抢同一段 |

## 6. 输出（`outputs/demo/` 为演示；正式输出在你自己指定的目录）

- `target_oligos.csv`：id / group / type / 区间 / 长度 / 方向 / 是否跨界 / 端基修饰 / **序列**
- `staples.csv`：id / group / **39 nt 序列** / scaffold 区间 / target 区间 / 长度 / domain_order / 两个 domain 的方向
- `design_config.yaml`：实际生效的配置（含派生值）
- `design_report.json`：状态、输入序列的 sha256、计数、9 项检查、Tm 汇总

## 7. 演示与"不许假装最终设计"

默认配置全是 `null`，`python -m ssdna.design --config params/design_config.yaml` 直接退出码 2。

要跑通流程必须显式加 `--demo`：目标序列是**确定性合成序列**（LCG 生成，不是任何基因），
配置里写入 `is_demo: true` 与说明，产物落在 `outputs/demo/`。这样"能跑"和"是真实设计"不会混淆。

## 8. 测试

`python -m unittest discover -s tests` → **133 个测试全部通过**（Task C 新增 26 个）。

覆盖：缺参数不落盘、oligo 区间连续与环状回绕、边界识别、序列等于原切片、
staple 计数/长度、两个臂的 revcomp 独立复算、scaffold 块不重叠且按 3′→5′ 反向推进、
翻转 `domain_order` 只交换两段、篡改一个碱基会被 `not_complementary` 抓到、
四个输出文件的字段与内容、YAML 往返、Tm 与相位扫描。

## 9. 尚未解决

1. `scaffold.strand` 仍未定 → 所有 13 nt 臂的方向都悬着（选错则整批作废）。
2. 目标序列未定；长度需满足 `L = 2C = 26n` 且 `|n_A − n_B| ≤ 1`。
3. `five_prime_phosphate` / `intended_ligation` 未定 → nick 表的化学列仍是 `unknown`。
4. 退火程序（温度曲线）没有建模：本模块只给 Tm 分布，不给程序建议。
