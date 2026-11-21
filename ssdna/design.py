"""Target oligo 与 staple 的序列生成引擎（Task C）。

设计立场（重要）：

1. **序列没有自由度**：staple 的 13 nt 与 26 nt 都由 scaffold / 目标序列唯一决定
   （反向互补），不存在"换一段序列把 Tm 调到 60 °C"这种事。所以 Tm 只能当**体检指标**
   与**少数自由度的取舍依据**（scaffold 相位、nick 位置、退火程序），不能当设计目标。
2. **不为了凑长度改动目标序列**：覆盖不满、长度不是 26 的倍数，只报告，不自动补/裁。
3. **参数不全就拒绝生成**：``find_missing`` 会列出全部缺失项，一条序列都不输出。
"""

from __future__ import annotations

import csv
import hashlib
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from .mapping import Check, _contains, _overlaps
from .models import BindingSite, Scaffold, Staple, TargetOligo, TargetRegion
from .sequence import get_circular_subsequence, normalize_sequence, reverse_complement

__all__ = [
    "DesignInputs",
    "DEFAULT_CONFIG",
    "load_config",
    "find_missing",
    "find_invalid",
    "design_target_oligos",
    "design_staples",
    "build_nicks",
    "check_design",
    "run_design",
    "to_report",
    "scan_scaffold_phase",
    "suggest_oligo_cuts",
]

DEFAULT_CONFIG: dict[str, Any] = {
    "schema": "scaffold_staple_design/design_config/v1",
    "scaffold": {
        "file": "minimalplasmid(1).dna",
        "strand": None,  # "+" | "-"：实验室提供的单链对应文件正链还是负链
        # "whole" 表示整段；[a, b) 表示文件坐标下的一段（b < a 表示跨原点）。
        # 必须显式给出——"整段"也要你写 "whole"，程序不替你假设。
        "region": None,
        "topology": "circular",
    },
    "staple": {
        "scaffold_binding_length": 13,
        "target_binding_length": 26,
        "total_length": 39,
        "group_pattern": "AB",
        "domain_order": "5p_target_26",
        "covered_start": 0,
    },
    "target": {
        "source": None,  # FASTA 路径；null 表示待提供
        "groups": ["A", "B"],
        "group_boundaries": None,  # [a]：A 为 [0,a)，B 为 [a,L)
        "oligo_cuts": None,  # 目标坐标上的切点（不含 0）
        "topology": "circular",
        "boundary_count": 2,
        "five_prime_phosphate": None,
        "intended_ligation": None,
    },
    "analysis": {
        "tm": True,
        "na_mM": 50.0,
        "mg_mM": 0.0,
        "dnac_nM": 250.0,
        "temperature_c": None,
    },
    "output": {"outdir": "outputs/design"},
}


# --------------------------------------------------------------------------- #
# 配置
# --------------------------------------------------------------------------- #
def _deep_merge(base: dict, override: dict) -> dict:
    out = dict(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(out.get(key), dict):
            out[key] = _deep_merge(out[key], value)
        else:
            out[key] = value
    return out


def load_config(path: str | Path) -> dict[str, Any]:
    """读 YAML（需要 PyYAML）；没有 PyYAML 时可以改用等价的 JSON。"""
    path = Path(path)
    text = path.read_text(encoding="utf-8")
    try:
        import yaml

        data = yaml.safe_load(text)
    except ImportError:
        try:
            data = json.loads(text)
        except json.JSONDecodeError as exc:  # pragma: no cover
            raise ImportError(
                "未安装 PyYAML，且该文件不是 JSON。"
                "请 pip install pyyaml，或改用 design_config.json"
            ) from exc
    if not isinstance(data, dict):
        raise ValueError(f"{path} 的顶层必须是映射（dict）")
    return _deep_merge(DEFAULT_CONFIG, data)


def find_missing(config: dict[str, Any]) -> list[str]:
    """列出阻塞参数（值为 None 的叶子）。只报告，不给默认值。"""
    missing: list[str] = []

    def walk(node: Any, prefix: str) -> None:
        if isinstance(node, dict):
            for key, value in node.items():
                if key.startswith("_"):
                    continue
                walk(value, f"{prefix}.{key}" if prefix else key)
        elif node is None:
            missing.append(prefix)

    walk(config, "")
    return sorted(missing)


def find_invalid(config: dict[str, Any]) -> dict[str, str]:
    """返回**类型/取值域**非法的参数映射（key -> 原因）。

    补充 find_missing 只查"缺没缺"、不查"对不对"的盲区：配置里把数字写成
    字符串、枚举值拼错、切点越界等，要在这里一步报清，而不是散落在生成流程深处
    抛难懂的 int()/float() ValueError。
    """
    invalid: dict[str, str] = {}

    def fail(key: str, why: str) -> None:
        invalid[key] = why

    sc = config.get("scaffold", {})
    if sc.get("strand") not in ("+", "-", None):
        fail("scaffold.strand", "必须是 '+' 或 '-'")
    region = sc.get("region")
    if region not in ("whole", None):
        if not (isinstance(region, list) and len(region) == 2
                and all(isinstance(v, int) and not isinstance(v, bool) for v in region)):
            fail("scaffold.region", "必须是 'whole' 或 [start, end)")

    st = config.get("staple", {})
    for key in ("scaffold_binding_length", "target_binding_length", "total_length",
                "covered_start"):
        value = st.get(key)
        if value is not None and (isinstance(value, bool) or not isinstance(value, int)):
            fail(f"staple.{key}", "必须是整数")
    if st.get("domain_order") not in ("5p_target_26", "5p_scaffold_13", None):
        fail("staple.domain_order", "必须是 5p_target_26 或 5p_scaffold_13")
    if st.get("group_pattern") is not None and any(
        ch not in "AB" for ch in str(st.get("group_pattern"))
    ):
        fail("staple.group_pattern", "只能含 A/B")
    if (st.get("total_length") is not None
            and isinstance(st.get("scaffold_binding_length"), int)
            and isinstance(st.get("target_binding_length"), int)
            and st["total_length"] != st["scaffold_binding_length"] + st["target_binding_length"]):
        fail("staple.total_length", "必须等于 scaffold_binding_length + target_binding_length")

    tg = config.get("target", {})
    if tg.get("topology") not in ("circular", "linear", None):
        fail("target.topology", "必须是 circular 或 linear")
    cuts = tg.get("oligo_cuts")
    if cuts is not None:
        if not isinstance(cuts, list) or any(
            isinstance(c, bool) or not isinstance(c, int) or c <= 0 for c in cuts
        ):
            fail("target.oligo_cuts", "必须是正整数列表")
    gb = tg.get("group_boundaries")
    if gb is not None and (not isinstance(gb, list) or len(gb) != 1
                           or not all(isinstance(v, int) and not isinstance(v, bool) for v in gb)):
        fail("target.group_boundaries", "必须恰好一个整数（A/B 分界）")

    an = config.get("analysis", {})
    for key in ("na_mM", "mg_mM", "dnac_nM", "temperature_c"):
        value = an.get(key)
        if value is not None and (isinstance(value, bool) or not isinstance(value, (int, float))):
            fail(f"analysis.{key}", "必须是数字")

    return invalid


# --------------------------------------------------------------------------- #
# 输入装配
# --------------------------------------------------------------------------- #
@dataclass
class DesignInputs:
    scaffold: Scaffold
    scaffold_strand: str
    working_sequence: str  # 已按 strand 定向的 5′→3′ scaffold 序列
    target_sequence: str
    target_length: int
    config: dict[str, Any] = field(default_factory=dict)
    notes: list[str] = field(default_factory=list)
    scaffold_region: str = "whole"
    scaffold_region_bounds: tuple[int, int] | None = None


def _orientation(seq: str, strand: str) -> str:
    if strand == "+":
        return seq
    if strand == "-":
        return reverse_complement(seq)
    raise ValueError(f"scaffold.strand 必须是 '+' 或 '-'，收到 {strand!r}")


# --------------------------------------------------------------------------- #
# 二、Target oligo 生成
# --------------------------------------------------------------------------- #
def design_target_oligos(
    target_sequence: str,
    group_definitions: list[dict[str, Any]],
    oligo_parameters: dict[str, Any],
    boundary_definitions: list[dict[str, Any]],
) -> tuple[list[TargetRegion], list[TargetOligo], list[Check]]:
    """把目标序列切成 oligo，并标出跨 A/B 边界的那些。

    Parameters
    ----------
    target_sequence:
        目标链的 5′→3′ 序列（**原样使用**，不补不裁）。
    group_definitions:
        ``[{"region_id": "target_A", "group": "A", "start": 0, "end": 364}, ...]``
    oligo_parameters:
        ``{"cuts": [100, 200, ...], "topology": "circular",
        "five_prime_modification": ..., "three_prime_modification": ...}``
    boundary_definitions:
        ``[{"junction_id": "J1", "position": 364, "upstream_group": "A",
        "downstream_group": "B"}, ...]``
    """
    seq = normalize_sequence(target_sequence)
    length = len(seq)
    topology = oligo_parameters.get("topology", "circular")
    checks: list[Check] = []

    regions = [
        TargetRegion(
            region_id=d["region_id"],
            group=d["group"],
            start=int(d["start"]),
            end=int(d["end"]),
            strand="+",
            sequence=get_circular_subsequence(
                seq, int(d["start"]) % length, int(d["end"]) - int(d["start"])
            ),
        )
        for d in group_definitions
    ]

    raw_cuts = oligo_parameters.get("cuts") or []
    if not raw_cuts:
        raise ValueError(
            "oligo_parameters['cuts'] 为空：不知道要切成几条（不自动等分）"
        )
    cuts: list[int] = []
    for value in raw_cuts:
        if isinstance(value, bool) or not isinstance(value, int):
            raise ValueError(f"切点必须是整数，收到 {value!r}")
        if not (0 < value < length):
            raise ValueError(
                f"切点 {value} 越界：必须满足 0 < cut < {length}。"
                "不接受负数或越界值——那会被取模悄悄变成另一个位置"
            )
        cuts.append(value)
    duplicates = sorted({c for c in cuts if cuts.count(c) > 1})
    if duplicates:
        raise ValueError(f"切点重复：{duplicates}（重复项会被静默合并，因此直接报错）")
    cuts.sort()

    bounds: list[tuple[int, int]] = []
    for i, start in enumerate(cuts):
        end = cuts[(i + 1) % len(cuts)]
        bounds.append((start, end if end > start else end + length))

    oligos: list[TargetOligo] = []
    for index, (start, end) in enumerate(bounds, start=1):
        piece = get_circular_subsequence(seq, start % length, end - start)
        spanning = [
            j
            for j in boundary_definitions
            if _contains(start, end, int(j["position"]) % length, length)
        ]
        segments = [
            r for r in regions if _overlaps(start, end, r.start, r.end, length)
        ]
        start_in_a = any(
            r.start <= (start % length) < r.end for r in regions if r.group == "A"
        )
        oligos.append(
            TargetOligo(
                oligo_id=f"T{index:02d}",
                group="A" if start_in_a else "B",
                oligo_type="boundary_spanning" if spanning else "ordinary",
                sequence=piece,
                target_start=start % length,
                target_end=end,
                target_segments=segments,
                strand="+",
                five_prime_modification=oligo_parameters.get(
                    "five_prime_modification", "unknown"
                ),
                three_prime_modification=oligo_parameters.get(
                    "three_prime_modification", "unknown"
                ),
            )
        )

    # 覆盖性检查
    covered = sum(o.length for o in oligos)
    if topology == "circular" and covered != length:
        checks.append(
            Check("error", "target_coverage",
                  f"oligo 总长 {covered} ≠ 目标长度 {length}")
        )
    for junction in boundary_definitions:
        hit = [
            o.oligo_id
            for o in oligos
            if _contains(o.target_start, o.target_end, int(junction["position"]) % length, length)
        ]
        if len(hit) != 1:
            checks.append(
                Check("error", "junction_coverage",
                      f"{junction['junction_id']} 被 {len(hit)} 条 oligo 覆盖")
            )

    return regions, oligos, checks



# --------------------------------------------------------------------------- #
# 三、Staple 生成
# --------------------------------------------------------------------------- #
def design_staples(
    scaffold: Scaffold,
    target_oligos: list[TargetOligo],
    staple_parameters: dict[str, Any],
    *,
    scaffold_strand: str = "+",
    target_sequence: str | None = None,
    working_sequence: str | None = None,
) -> tuple[list[Staple], list[Check]]:
    """按 13+26 平铺生成 staple。

    staple 序列 = ``5′-[revcomp(目标 26)]-[revcomp(scaffold 13)]-3′``；
    折返只影响空间走向，不改变这条链的真实 5′→3′ 序列。
    """
    arm_s = int(staple_parameters["scaffold_binding_length"])
    arm_t = int(staple_parameters["target_binding_length"])
    total = int(staple_parameters.get("total_length", arm_s + arm_t))
    pattern = staple_parameters.get("group_pattern", "AB")
    if not pattern or any(ch not in "AB" for ch in pattern):
        raise ValueError(f"group_pattern 只能含 A/B，收到 {pattern!r}")
    covered_start = int(staple_parameters.get("covered_start", 0))
    order = staple_parameters.get("domain_order", "5p_target_26")
    if order not in ("5p_target_26", "5p_scaffold_13"):
        raise ValueError(f"domain_order 必须是 5p_target_26 或 5p_scaffold_13，收到 {order!r}")
    checks: list[Check] = []

    if total != arm_s + arm_t:
        checks.append(
            Check("error", "total_length",
                  f"total_length={total} 与 {arm_s}+{arm_t}={arm_s + arm_t} 不一致")
        )
    if scaffold_strand not in ("+", "-"):
        raise ValueError(f"scaffold_strand 必须是 '+'/'-'，收到 {scaffold_strand!r}")

    # 工作序列：region 子区段时由调用方传入（以工作序列为 0 点），否则用整段定向序列。
    # 不能用 scaffold.sequence 重新算——否则 region 子区段的坐标会错位。
    working = (
        working_sequence
        if working_sequence is not None
        else _orientation(scaffold.sequence, scaffold_strand)
    )
    L = len(working)

    # 每组的目标跨度：优先用调用方给的 group_spans（环状目标时从 oligo 反推会被
    # 跨原点的区间带偏，所以不默认反推）
    span: dict[str, tuple[int, int]] = {
        str(g): (int(v[0]), int(v[1]))
        for g, v in (staple_parameters.get("group_spans") or {}).items()
    }
    if not span:
        for group in ("A", "B"):
            items = [o for o in target_oligos if o.group == group]
            if not items:
                continue
            span[group] = (
                min(o.target_start for o in items),
                max(o.target_end for o in items),
            )
    # 26 nt 臂按 26 的节距铺，铺到哪儿算哪儿；剩下的尾巴不锚定，**不是错误**
    leftovers: dict[str, int] = {}
    for group in ("A", "B"):
        if group not in span:
            continue
        s, e = span[group]
        leftovers[group] = (e - s) % arm_t
        if leftovers[group]:
            checks.append(
                Check(
                    "note",
                    "uncovered_target_tail",
                    f"{group} 组 {e - s} nt，{arm_t} nt 臂铺满 "
                    f"{(e - s) // arm_t} 条后，末尾 {leftovers[group]} nt 不被臂覆盖"
                    "（可以留作交界 linker，也可以调整分组边界）",
                )
            )

    n_by_group = {g: (e - s) // arm_t for g, (s, e) in span.items()}
    n_total = sum(n_by_group.values())

    # 严格校验：分组模式实际给每组分配的条数，必须等于目标跨度要求的条数。
    # 否则（例如 group_boundaries=[26]，要求 n_A=1/n_B=27，而 AB 交替给 14/14）
    # 会出现 A 臂铺进 B 区、目标漏铺，产生结构上不可能的序列而不报错。
    pattern_counts: dict[str, int] = {}
    for index in range(n_total):
        group = pattern[index % len(pattern)]
        pattern_counts[group] = pattern_counts.get(group, 0) + 1
    for group in ("A", "B"):
        required = n_by_group.get(group, 0)
        actual = pattern_counts.get(group, 0)
        if required != actual:
            checks.append(Check(
                "error", "group_pattern_mismatch",
                f"分组模式 {pattern!r} 给 {group} 组 {actual} 条 staple，"
                f"但目标跨度要求 {required} 条（n_A={n_by_group.get('A', 0)}, "
                f"n_B={n_by_group.get('B', 0)}）。请调整 group_boundaries 使两组条数"
                "只差 1，或改用匹配的分组模式。",
            ))
    if any(c.code == "group_pattern_mismatch" for c in checks):
        return [], checks

    covered = arm_s * n_total
    if covered > L:
        checks.append(
            Check("error", "scaffold_overflow",
                  f"需要覆盖 {covered} nt，超过 scaffold 长度 {L}")
        )

    target_seq = normalize_sequence(target_sequence) if target_sequence else None
    if target_seq is not None:
        for group in ("A", "B"):
            if group in span and span[group][1] > len(target_seq):
                checks.append(
                    Check("error", "target_out_of_range",
                          f"{group} 组结束于 {span[group][1]}，超过目标长度 {len(target_seq)}")
                )

    staples: list[Staple] = []
    counters = {"A": 0, "B": 0}
    for index in range(n_total):
        group = pattern[index % len(pattern)]
        k = counters[group]
        counters[group] += 1

        p_start = index * arm_s
        raw_start = covered_start + covered - p_start - arm_s
        s_start = raw_start % L
        scaffold_site = BindingSite(
            molecule="scaffold", start=s_start, end=s_start + arm_s,
            strand="-", molecule_length=L,
        )
        t_start = span[group][0] + k * arm_t
        target_site = BindingSite(
            molecule="target", start=t_start, end=t_start + arm_t,
            strand="-", molecule_length=len(target_seq) if target_seq else span["B"][1],
        )

        scaffold_piece = get_circular_subsequence(working, s_start, arm_s)
        if target_seq is None:
            raise ValueError("design_staples 需要 target_sequence 才能拼出真实序列")
        target_piece = get_circular_subsequence(target_seq, t_start, arm_t)

        target_arm = reverse_complement(target_piece)
        scaffold_arm = reverse_complement(scaffold_piece)
        sequence = (
            target_arm + scaffold_arm if order == "5p_target_26" else scaffold_arm + target_arm
        )

        staples.append(
            Staple(
                staple_id=f"S{index + 1:02d}",
                group=group,
                scaffold_binding_region=scaffold_site,
                target_binding_region=target_site,
                sequence=sequence,
                domain_order=order,
            )
        )

    return staples, checks


# --------------------------------------------------------------------------- #
# 五、生成后检查
# --------------------------------------------------------------------------- #
def check_design(
    scaffold: Scaffold,
    working_sequence: str,
    target_sequence: str,
    regions: list[TargetRegion],
    oligos: list[TargetOligo],
    staples: list[Staple],
    *,
    junctions: list[dict[str, Any]],
    arm_s: int,
    arm_t: int,
    nicks: list[Any] | None = None,
) -> list[Check]:
    checks: list[Check] = []
    L = len(working_sequence)
    T = len(target_sequence)

    # 1) 长度
    for staple in staples:
        if staple.scaffold_domain_length != arm_s:
            checks.append(Check("error", "staple_scaffold_length",
                                f"{staple.staple_id} 的 13 nt 臂长度不对"))
        if staple.target_domain_length != arm_t:
            checks.append(Check("error", "staple_target_length",
                                f"{staple.staple_id} 的 26 nt 臂长度不对"))
        if len(staple.sequence) != arm_s + arm_t:
            checks.append(Check("error", "staple_total_length",
                                f"{staple.staple_id} 序列长 {len(staple.sequence)}"))

    # 2) 反向互补关系（逐条实算，不信构造过程）
    for staple in staples:
        s_site, t_site = staple.scaffold_binding_region, staple.target_binding_region
        scaffold_piece = get_circular_subsequence(working_sequence, s_site.start, arm_s)
        target_piece = get_circular_subsequence(target_sequence, t_site.start, arm_t)
        expected = (
            reverse_complement(target_piece) + reverse_complement(scaffold_piece)
            if staple.domain_order == "5p_target_26"
            else reverse_complement(scaffold_piece) + reverse_complement(target_piece)
        )
        if staple.sequence != expected:
            checks.append(Check("error", "not_complementary",
                                f"{staple.staple_id} 与 scaffold/目标不互补"))

    # 3) 目标覆盖
    covered = sum(o.length for o in oligos)
    if covered != T:
        checks.append(Check("error", "target_coverage",
                            f"oligo 覆盖 {covered} nt ≠ 目标 {T} nt"))

    # 4) A/B 分组正确性
    a_region = [r for r in regions if r.group == "A"]
    for oligo in oligos:
        inside = [r for r in a_region if r.start <= oligo.target_start < r.end]
        expect = "A" if inside else "B"
        if oligo.group != expect:
            checks.append(Check("error", "group_mismatch",
                                f"{oligo.oligo_id} 起点在 {expect} 区，却标成 {oligo.group}"))

    # 5) 两个边界都被覆盖
    for junction in junctions:
        hit = [
            o.oligo_id
            for o in oligos
            if _contains(o.target_start, o.target_end, int(junction["position"]) % T, T)
        ]
        if len(hit) != 1:
            checks.append(Check("error", "junction_coverage",
                                f"{junction['junction_id']} 被 {len(hit)} 条 oligo 覆盖"))

    # 6) scaffold 映射越界
    for staple in staples:
        site = staple.scaffold_binding_region
        if site.length != arm_s or site.start < 0:
            checks.append(Check("error", "scaffold_bounds",
                                f"{staple.staple_id} 的 scaffold 区间越界"))
        occupied = {(site.start + i) % L for i in range(site.length)}
        if len(occupied) != arm_s:
            checks.append(Check("error", "scaffold_bounds",
                                f"{staple.staple_id} 的 scaffold 区间回绕异常"))

    # 7) 未分配的目标区间
    assigned = set()
    for oligo in oligos:
        for offset in range(oligo.length):
            assigned.add((oligo.target_start + offset) % T)
    if len(assigned) != T:
        missing = sorted(set(range(T)) - assigned)[:10]
        checks.append(Check("error", "unassigned_target",
                            f"{T - len(assigned)} nt 未被任何 oligo 覆盖，例如 {missing}"))

    # 8) 重复/冲突映射
    seen_scaffold: dict[int, str] = {}
    for staple in staples:
        start = staple.scaffold_binding_region.start
        if start in seen_scaffold:
            checks.append(Check("error", "duplicate_scaffold_site",
                                f"{staple.staple_id} 与 {seen_scaffold[start]} 抢同一段 scaffold"))
        seen_scaffold[start] = staple.staple_id
    seen_target: dict[int, str] = {}
    for staple in staples:
        start = staple.target_binding_region.start
        if start in seen_target:
            checks.append(Check("error", "duplicate_target_site",
                                f"{staple.staple_id} 与 {seen_target[start]} 抢同一段目标"))
        seen_target[start] = staple.staple_id

    # 9) nick 不得压在 staple 臂接缝上；每条 oligo 至少跨过一个接缝
    seams = set()
    for staple in staples:
        seams.add(staple.target_binding_region.start % T)
        seams.add(staple.target_binding_region.end % T)
    for nick in nicks or []:
        if nick.genomic_or_target_position % T in seams:
            checks.append(Check(
                "error", "nick_on_seam",
                f"{nick.nick_id} 落在 staple 臂接缝 "
                f"{nick.genomic_or_target_position % T} 处",
            ))
    for oligo in oligos:
        span = (oligo.target_end - oligo.target_start) % T or T
        if not any(
            0 < ((seam - oligo.target_start) % T) < span for seam in seams
        ):
            checks.append(Check(
                "error", "oligo_spans_no_seam",
                f"{oligo.oligo_id} 没有跨过任何 staple 臂接缝",
            ))

    # 10) 提示
    checks.append(Check("note", "engine",
                        "序列由 scaffold/目标唯一决定；Tm 仅作体检，未参与选序列"))
    return checks


# --------------------------------------------------------------------------- #
# 主流程 + 输出
# --------------------------------------------------------------------------- #
def run_design(
    inputs: DesignInputs,
    *,
    outdir: str | Path,
    progress: Any = None,
) -> dict[str, Any]:
    """完整流程。``progress`` 是可选回调，用于命令行逐步打印。"""

    def step(name: str) -> None:
        if progress is not None:
            progress(name)

    cfg = inputs.config
    staple_cfg = cfg["staple"]
    target_cfg = cfg["target"]
    arm_s = int(staple_cfg["scaffold_binding_length"])
    arm_t = int(staple_cfg["target_binding_length"])

    boundaries = target_cfg.get("group_boundaries") or []
    if len(boundaries) != 1:
        raise ValueError(
            f"target.group_boundaries 必须恰好一个元素（A/B 分界），收到 {boundaries!r}"
        )
    declared_boundary_count = int(target_cfg.get("boundary_count", 2))
    if cfg["target"]["topology"] == "circular" and declared_boundary_count != 2:
        raise ValueError(
            f"环状拓扑下 boundary_count 应为 2，配置里是 {declared_boundary_count}"
        )
    if cfg["target"]["topology"] == "linear" and declared_boundary_count != 1:
        raise ValueError(
            f"线性拓扑下 boundary_count 应为 1，配置里是 {declared_boundary_count}"
        )
    group_defs = []
    prev = 0
    for i, edge in enumerate(list(boundaries) + [inputs.target_length]):
        group = "A" if i == 0 else "B"
        group_defs.append(
            {"region_id": f"target_{group}", "group": group, "start": prev, "end": edge}
        )
        prev = edge

    junction_defs = [
        {"junction_id": f"J{i + 1}", "position": int(edge),
         "upstream_group": "A", "downstream_group": "B"}
        for i, edge in enumerate(boundaries)
    ]
    if cfg["target"]["topology"] == "circular":
        junction_defs.append(
            {"junction_id": f"J{len(junction_defs) + 1}", "position": 0,
             "upstream_group": "B", "downstream_group": "A"}
        )

    step("build_mapping")
    regions, oligos, checks_1 = design_target_oligos(
        inputs.target_sequence,
        group_defs,
        {
            "cuts": target_cfg["oligo_cuts"],
            "topology": cfg["target"]["topology"],
            "five_prime_modification": target_cfg.get("five_prime_phosphate") or "unknown",
            "three_prime_modification": str(
                target_cfg.get("three_prime_modification") or "OH"
            ),
        },
        junction_defs,
    )
    step("design_oligos")
    step("design_staples")
    staples, checks_2 = design_staples(
        inputs.scaffold,
        oligos,
        {
            "scaffold_binding_length": arm_s,
            "target_binding_length": arm_t,
            "total_length": int(staple_cfg.get("total_length", arm_s + arm_t)),
            "group_pattern": staple_cfg.get("group_pattern", "AB"),
            "covered_start": int(staple_cfg.get("covered_start", 0)),
            "domain_order": staple_cfg.get("domain_order", "5p_target_26"),
            "group_spans": {r.group: [r.start, r.end] for r in regions},
        },
        scaffold_strand=inputs.scaffold_strand,
        target_sequence=inputs.target_sequence,
        working_sequence=inputs.working_sequence,
    )
    step("check_design")
    nicks = build_nicks(
        oligos,
        topology=cfg["target"]["topology"],
        five_prime_phosphate=cfg["target"].get("five_prime_phosphate") or "unknown",
        intended_ligation=bool(cfg["target"].get("intended_ligation")),
        three_prime_hydroxyl=str(cfg["target"].get("three_prime_hydroxyl") or "unknown"),
        nick_overrides=cfg["target"].get("nick_overrides") or {},
    )
    checks_3 = check_design(
        inputs.scaffold,
        inputs.working_sequence,
        inputs.target_sequence,
        regions,
        oligos,
        staples,
        junctions=junction_defs,
        arm_s=arm_s,
        arm_t=arm_t,
        nicks=nicks,
    )

    payload = to_report(inputs, regions, oligos, staples, junction_defs,
                        checks_1 + checks_2 + checks_3, nicks)
    payload["_sequences"] = {
        "target": inputs.target_sequence,
        "scaffold_working": inputs.working_sequence,
        "oligos": {o.oligo_id: o.sequence for o in oligos},
        "staples": {
            s.staple_id: {
                "sequence": s.sequence,
                "scaffold_domain": (
                    s.sequence[arm_t : arm_t + arm_s]
                    if s.domain_order == "5p_target_26"
                    else s.sequence[:arm_s]
                ),
                "target_domain": (
                    s.sequence[:arm_t]
                    if s.domain_order == "5p_target_26"
                    else s.sequence[arm_s : arm_s + arm_t]
                ),
            }
            for s in staples
        },
    }
    analysis = inputs.config.get("analysis", {})
    if analysis.get("tm"):
        from .thermo import annotate_staples, summarize

        annotate_staples(
            payload["_sequences"]["staples"],
            na_mM=float(analysis.get("na_mM") or 50.0),
            mg_mM=float(analysis.get("mg_mM") or 0.0),
            dnac_nM=float(analysis.get("dnac_nM") or 250.0),
        )
        payload["_tm_summary"] = summarize(payload["_sequences"]["staples"])
        if payload["_tm_summary"].get("available"):
            payload["tm_summary"] = payload["_tm_summary"]
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    write_outputs(payload, out)
    return payload


def build_nicks(
    oligos,
    *,
    topology,
    five_prime_phosphate,
    intended_ligation,
    three_prime_hydroxyl="unknown",
    nick_overrides: dict[str, dict] | None = None,
):
    """相邻 target oligo 的 3′→5′ 相邻关系，逐条显式记录。

    默认值只用于"没单独配置"的 nick；``nick_overrides`` 可以按 ``N1`` / ``N2``
    逐条覆盖端基与连接状态，从而表达混合修饰与部分连接。

    注意：``intended_ligation`` 是**意图**，只翻译成 ``ligation_state="planned"``；
    要表达"实验已确认连上"，必须显式写 ``ligation_state="ligated"``。
    """
    from .models import Nick

    overrides = nick_overrides or {}
    ordered = sorted(oligos, key=lambda o: o.target_start)
    nicks = []
    for i, upstream in enumerate(ordered):
        if topology == "linear" and i == len(ordered) - 1:
            continue
        downstream = ordered[(i + 1) % len(ordered)]
        nick_id = f"N{i + 1}"
        override = dict(overrides.get(nick_id, {}))
        intended = bool(override.get("intended_ligation", intended_ligation))
        state = override.get("ligation_state") or ("planned" if intended else "none")
        nicks.append(
            Nick(
                nick_id=nick_id,
                upstream_oligo=upstream.oligo_id,
                downstream_oligo=downstream.oligo_id,
                # 用下游起点（已归一化到 [0, L)），避免跨原点 nick 存成未回绕值
                genomic_or_target_position=downstream.target_start,
                five_prime_phosphate=str(
                    override.get("five_prime_phosphate", five_prime_phosphate)
                ),
                three_prime_hydroxyl=str(
                    override.get("three_prime_hydroxyl", three_prime_hydroxyl)
                ),
                intended_ligation=intended,
                ligation_state=state,
            )
        )
    return nicks


def to_report(inputs, regions, oligos, staples, junctions, checks, nicks=None) -> dict[str, Any]:
    errors = [c for c in checks if c.level == "error"]
    return {
        "schema": "scaffold_staple_design/design_report/v1",
        "status": "BLOCKED" if errors else "OK",
        "config": inputs.config,
        "inputs": {
            "scaffold_source": inputs.scaffold.source,
            "scaffold_file_sha256": inputs.scaffold.sha256,
            "scaffold_sequence_sha256": inputs.scaffold.sequence_sha256,
            "scaffold_strand": inputs.scaffold_strand,
            "working_sequence_sha256": hashlib.sha256(
                inputs.working_sequence.encode()
            ).hexdigest(),
            "target_length": inputs.target_length,
            "target_sequence_sha256": hashlib.sha256(
                inputs.target_sequence.encode()
            ).hexdigest(),
        },
        "scaffold_region": {
            "mode": getattr(inputs, "scaffold_region", None) or "whole",
            "file_coordinates": getattr(inputs, "scaffold_region_bounds", None),
            "working_length": len(inputs.working_sequence),
            "coordinate_map": "staple 的 scaffold 坐标以 working_sequence 为 0 点；"
                              "换算回文件坐标需加上 file_coordinates[0]（整段时为 0）",
        },
        "regions": [
            {"region_id": r.region_id, "group": r.group, "start": r.start,
             "end": r.end, "length": r.length}
            for r in regions
        ],
        "target_oligos": [
            {
                "oligo_id": o.oligo_id, "group": o.group, "oligo_type": o.oligo_type,
                "target_start": o.target_start, "target_end": o.target_end,
                "length": o.length, "strand": o.strand,
                "crosses_junction": o.oligo_type == "boundary_spanning",
                "five_prime_modification": o.five_prime_modification,
                "three_prime_modification": o.three_prime_modification,
            }
            for o in oligos
        ],
        "staples": [
            {
                "staple_id": s.staple_id, "group": s.group,
                "scaffold_start": s.scaffold_binding_region.start,
                "scaffold_end": s.scaffold_binding_region.end,
                "scaffold_length": s.scaffold_domain_length,
                "target_start": s.target_binding_region.start,
                "target_end": s.target_binding_region.end,
                "target_length": s.target_domain_length,
                "total_length": s.total_length,
                "domain_order": s.domain_order,
                "scaffold_strand": s.scaffold_binding_region.strand,
                "target_strand": s.target_binding_region.strand,
            }
            for s in staples
        ],
        "junctions": junctions,
        "nicks": [
            {
                "nick_id": n.nick_id,
                "upstream_oligo": n.upstream_oligo,
                "upstream_end": n.upstream_end,
                "downstream_oligo": n.downstream_oligo,
                "downstream_end": n.downstream_end,
                "genomic_or_target_position": n.genomic_or_target_position,
                "five_prime_phosphate": n.five_prime_phosphate,
                "three_prime_hydroxyl": n.three_prime_hydroxyl,
                "intended_ligation": n.intended_ligation,
                "ligation_state": n.ligation_state,
            }
            for n in (nicks or [])
        ],
        "checks": [{"level": c.level, "code": c.code, "message": c.message} for c in checks],
        "error_count": len(errors),
        "counts": {
            "target_oligos": len(oligos),
            "boundary_oligos": sum(1 for o in oligos if o.oligo_type == "boundary_spanning"),
            "staples": len(staples),
            "staples_A": sum(1 for s in staples if s.group == "A"),
            "staples_B": sum(1 for s in staples if s.group == "B"),
            "covered_scaffold_nt": sum(s.scaffold_domain_length for s in staples),
        },
    }


def write_outputs(payload: dict[str, Any], outdir: Path) -> dict[str, Path]:
    paths: dict[str, Path] = {}

    oligo_csv = outdir / "target_oligos.csv"
    with oligo_csv.open("w", encoding="utf-8-sig", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow([
            "oligo_id", "group", "oligo_type", "target_start", "target_end",
            "length", "strand", "crosses_junction", "five_prime_modification",
            "three_prime_modification", "sequence",
        ])
        for o in payload["target_oligos"]:
            writer.writerow([
                o["oligo_id"], o["group"], o["oligo_type"], o["target_start"],
                o["target_end"], o["length"], o["strand"], o["crosses_junction"],
                o["five_prime_modification"], o["three_prime_modification"],
                payload["_sequences"]["oligos"][o["oligo_id"]],
            ])
    paths["target_oligos_csv"] = oligo_csv

    staple_csv = outdir / "staples.csv"
    with staple_csv.open("w", encoding="utf-8-sig", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow([
            "staple_id", "group", "sequence", "scaffold_start", "scaffold_end",
            "scaffold_length", "target_start", "target_end", "target_length",
            "total_length", "domain_order", "scaffold_strand", "target_strand",
            "tm_scaffold_domain", "tm_target_domain",
        ])
        for s in payload["staples"]:
            row = payload["_sequences"]["staples"][s["staple_id"]]
            writer.writerow([
                s["staple_id"], s["group"], row["sequence"], s["scaffold_start"],
                s["scaffold_end"], s["scaffold_length"], s["target_start"],
                s["target_end"], s["target_length"], s["total_length"],
                s["domain_order"], s["scaffold_strand"], s["target_strand"],
                row.get("tm_scaffold_domain", ""), row.get("tm_target_domain", ""),
            ])
    paths["staples_csv"] = staple_csv

    nicks_csv = outdir / "nicks.csv"
    with nicks_csv.open("w", encoding="utf-8-sig", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow([
            "nick_id", "upstream_oligo", "upstream_end", "downstream_oligo",
            "downstream_end", "genomic_or_target_position", "five_prime_phosphate",
            "three_prime_hydroxyl", "intended_ligation", "ligation_state",
        ])
        for n in payload["nicks"]:
            writer.writerow([
                n["nick_id"], n["upstream_oligo"], n["upstream_end"],
                n["downstream_oligo"], n["downstream_end"],
                n["genomic_or_target_position"], n["five_prime_phosphate"],
                n["three_prime_hydroxyl"], n["intended_ligation"],
                n["ligation_state"],
            ])
    paths["nicks_csv"] = nicks_csv

    yaml_path = outdir / "design_config.yaml"
    try:
        import yaml

        yaml_path.write_text(
            yaml.safe_dump(payload["config"], allow_unicode=True, sort_keys=False),
            encoding="utf-8",
        )
    except ImportError:  # pragma: no cover
        yaml_path.write_text(
            json.dumps(payload["config"], ensure_ascii=False, indent=2), encoding="utf-8"
        )
    paths["design_config_yaml"] = yaml_path

    target_fa = outdir / "target.fasta"
    with target_fa.open("w", encoding="ascii") as fh:
        fh.write(">target length=%d\n" % len(payload["_sequences"]["target"]))
        seq = payload["_sequences"]["target"]
        for i in range(0, len(seq), 70):
            fh.write(seq[i : i + 70] + "\n")
    paths["target_fasta"] = target_fa

    scaffold_fa = outdir / "scaffold_working.fasta"
    with scaffold_fa.open("w", encoding="ascii") as fh:
        seq = payload["_sequences"]["scaffold_working"]
        fh.write(">scaffold_working strand=%s length=%d\n"
                 % (payload["inputs"]["scaffold_strand"], len(seq)))
        for i in range(0, len(seq), 70):
            fh.write(seq[i : i + 70] + "\n")
    paths["scaffold_fasta"] = scaffold_fa

    report_path = outdir / "design_report.json"
    public = {k: v for k, v in payload.items() if not k.startswith("_")}
    report_path.write_text(
        json.dumps(public, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    paths["design_report_json"] = report_path
    return paths

# --------------------------------------------------------------------------- #
# 命令行
# --------------------------------------------------------------------------- #
def _synthetic_target(length: int, seed: int = 2026) -> str:
    """确定性的合成序列，**只用于跑通流程**，不是任何真实基因。"""
    state = seed
    out = []
    for _ in range(length):
        state = (1103515245 * state + 12345) & 0x7FFFFFFF
        out.append("ACGT"[(state >> 16) & 3])
    return "".join(out)


def demo_config(base: dict[str, Any]) -> dict[str, Any]:
    """把 null 填成合成值，并显式标记 is_demo。"""
    cfg = json.loads(json.dumps(base))
    cfg["is_demo"] = True
    cfg["_demo_note"] = "target 为合成序列，仅用于验证引擎，不是真实设计"
    cfg["scaffold"].update({"strand": "+", "region": "whole"})
    cfg["analysis"].update({"temperature_c": 37.0})
    cfg["target"].update({
        "source": {"inline": _synthetic_target(728)},
        "group_boundaries": [364],
        "oligo_cuts": [108, 212, 316, 420, 524, 628],
        "topology": "circular",
        "five_prime_phosphate": "no",
        "intended_ligation": False,
    })
    return cfg


def _load_target(cfg: dict[str, Any]) -> str:
    source = cfg["target"]["source"]
    if isinstance(source, dict) and "inline" in source:
        return normalize_sequence(source["inline"])
    from .io import load_fasta

    return load_fasta(source).sequence


def main(argv: list[str] | None = None) -> int:
    import argparse

    parser = argparse.ArgumentParser(description="生成 target oligo 与 staple 序列")
    parser.add_argument("--config", help="design_config.yaml（或 .json）")
    parser.add_argument("--outdir", default=None)
    parser.add_argument("--demo", action="store_true",
                        help="用合成目标序列跑通流程（输出会标记 is_demo）")
    parser.add_argument("--scan-phase", action="store_true",
                        help="扫描 scaffold 相位，比较 13 nt 臂的 Tm 分布（不生成序列）")
    args = parser.parse_args(argv)

    cfg = load_config(args.config) if args.config else json.loads(json.dumps(DEFAULT_CONFIG))
    if args.demo:
        cfg = demo_config(cfg)
    if args.outdir:
        cfg["output"]["outdir"] = args.outdir

    missing = find_missing(cfg)
    if missing:
        print(f"缺少 {len(missing)} 个必需参数，未生成任何序列：")
        for key in missing:
            print(f"  - {key}")
        print("\n补齐后重跑；工具不会代入默认值。")
        return 2

    from .io import load_snapgene

    scaffold = load_snapgene(cfg["scaffold"]["file"])
    target_sequence = _load_target(cfg)
    strand = cfg["scaffold"]["strand"]
    inputs = DesignInputs(
        scaffold=scaffold,
        scaffold_strand=strand,
        working_sequence=_orientation(scaffold.sequence, strand),
        target_sequence=target_sequence,
        target_length=len(target_sequence),
        config=cfg,
    )
    if args.scan_phase:
        rows = scan_scaffold_phase(
            scaffold,
            {
                "scaffold_binding_length": int(cfg["staple"]["scaffold_binding_length"]),
                "target_binding_length": int(cfg["staple"]["target_binding_length"]),
                "group_spans": {
                    "A": [0, int(cfg["target"]["group_boundaries"][0])],
                    "B": [int(cfg["target"]["group_boundaries"][0]), inputs.target_length],
                },
            },
            scaffold_strand=strand,
            na_mM=float(cfg["analysis"].get("na_mM") or 50.0),
            mg_mM=float(cfg["analysis"].get("mg_mM") or 0.0),
            dnac_nM=float(cfg["analysis"].get("dnac_nM") or 250.0),
        )
        print(f"相位扫描（{len(rows)} 个位置，按 13 nt 臂 Tm 标准差升序）：")
        for row in rows[:5]:
            print(f"  covered_start={row['covered_start']:4d} "
                  f"stdev={row['tm_stdev']:5.2f} mean={row['tm_mean']:5.2f} "
                  f"range=[{row['tm_min']}, {row['tm_max']}]")
        return 0

    payload = run_design(inputs, outdir=cfg["output"]["outdir"])

    print(f"[design] status={payload['status']} "
          f"oligos={payload['counts']['target_oligos']} "
          f"staples={payload['counts']['staples']} "
          f"covered={payload['counts']['covered_scaffold_nt']} nt")
    for check in payload["checks"]:
        if check["level"] == "error":
            print(f"  ERROR [{check['code']}] {check['message']}")
    tm = payload.get("tm_summary", {})
    if tm.get("available"):
        for key in ("tm_scaffold_domain", "tm_target_domain"):
            if key in tm:
                s = tm[key]
                print(f"  {key}: n={s['n']} mean={s['mean']} "
                      f"range=[{s['min']}, {s['max']}] stdev={s['stdev']}")
    print(f"  -> {cfg['output']['outdir']}")
    return 0 if payload["status"] == "OK" else 1


# --------------------------------------------------------------------------- #
# 唯一可以用 Tm 做决策的地方：scaffold 平铺相位
# --------------------------------------------------------------------------- #
def scan_scaffold_phase(
    scaffold: Scaffold,
    staple_parameters: dict[str, Any],
    *,
    scaffold_strand: str = "+",
    step: int = 13,
    na_mM: float = 50.0,
    mg_mM: float = 0.0,
    dnac_nM: float = 250.0,
) -> list[dict[str, Any]]:
    """扫描覆盖区起点，比较每批 13 nt 臂的 Tm 分布。

    这是"用 Tm 做设计"的**唯一合法方式**：序列不能改（必须互补），
    能改的只是把目标放在 scaffold 的哪一段、平铺从哪里起。
    判据：13 nt 臂的 Tm 标准差越小，整批 staple 的锚定强度越一致。

    注意：只用 Tm 排序，不改变任何一条 staple 的序列。
    """
    from statistics import fmean, pstdev

    from .thermo import tm_nn

    L = len(_orientation(scaffold.sequence, scaffold_strand))
    arm_s = int(staple_parameters["scaffold_binding_length"])
    arm_t = int(staple_parameters["target_binding_length"])

    spans = staple_parameters.get("group_spans") or {}
    n_total = sum((int(e) - int(s)) // arm_t for s, e in spans.values())
    covered = arm_s * n_total
    if covered > L:
        raise ValueError(f"覆盖 {covered} nt 超过 scaffold {L} nt")

    working = _orientation(scaffold.sequence, scaffold_strand)
    rows: list[dict[str, Any]] = []
    for start in range(0, L, max(1, step)):
        values = []
        for index in range(n_total):
            raw = start + covered - index * arm_s - arm_s
            s_start = raw % L
            piece = get_circular_subsequence(working, s_start, arm_s)
            value = tm_nn(piece, na_mM=na_mM, mg_mM=mg_mM, dnac_nM=dnac_nM)
            if value is not None:
                values.append(value)
        if not values:
            continue
        rows.append({
            "covered_start": start,
            "n": len(values),
            "tm_min": round(min(values), 2),
            "tm_max": round(max(values), 2),
            "tm_mean": round(fmean(values), 2),
            "tm_stdev": round(pstdev(values), 2) if len(values) > 1 else 0.0,
        })
    rows.sort(key=lambda r: r["tm_stdev"])
    return rows

# --------------------------------------------------------------------------- #
# 切点建议：目标是"够长 + nick 不压 staple 接缝"，不是"长度必须是 26 的倍数"
# --------------------------------------------------------------------------- #
def suggest_oligo_cuts(
    target_length: int,
    *,
    seams: list[int],
    min_length: int = 40,
    max_length: int = 120,
    clearance: int = 3,
    topology: str = "circular",
) -> list[int]:
    """给出满足三条约束的切点建议（只建议，不改目标序列）。

    约束：
    1. 每条 oligo 长度 ∈ [min_length, max_length]；
    2. nick 不能落在 staple 臂接缝上，且至少离开 ``clearance`` nt；
    3. 每条 oligo 内部至少包含一个接缝（否则它只被一条 staple 锚定）。

    目标总长**不需要**是 26 的整数倍：接缝只决定"不能切在哪"，不决定"必须切几个"。
    """
    if target_length <= 0:
        raise ValueError("target_length 必须为正")
    if min_length >= max_length:
        raise ValueError("min_length 必须小于 max_length")

    L = target_length
    seam_set = {s % L for s in seams}

    def near_seam(position: int) -> bool:
        return any(
            min(abs(position - seam), L - abs(position - seam)) <= clearance
            for seam in seam_set
        )

    def spans_seam(start: int, end: int) -> bool:
        if start < end:
            return any(start < seam < end for seam in seam_set)
        return any(seam > start for seam in seam_set) or any(seam < end for seam in seam_set)

    def build(n: int) -> list[int] | None:
        cuts: list[int] = []
        prev = 0
        for i in range(1, n):
            ideal = i * L / n
            best = None
            for delta in range(0, L):
                for candidate in (round(ideal) + delta, round(ideal) - delta):
                    if not (0 < candidate < L):
                        continue
                    if near_seam(candidate):
                        continue
                    if not (min_length <= candidate - prev <= max_length):
                        continue
                    if i == n - 1 and not (min_length <= L - candidate <= max_length):
                        continue
                    if not spans_seam(prev, candidate):
                        continue
                    best = candidate
                    break
                if best is not None:
                    break
            if best is None:
                return None
            cuts.append(best)
            prev = best
        if cuts:
            last_len = L - cuts[-1]
            if not (min_length <= last_len <= max_length) or not spans_seam(cuts[-1], L):
                return None
        return sorted(cuts)

    n_min = max(1, -(-L // max_length))
    for n in range(n_min, n_min + 9):
        cuts = build(n)
        if cuts is not None:
            return cuts
    raise ValueError(
        f"在 [{min_length}, {max_length}] nt 与接缝避让 {clearance} nt 的约束下切不开 {L} nt 的目标；"
        "请放宽长度区间或减少 staple（接缝太密）"
    )


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
