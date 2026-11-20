"""scaffold–staple–target 的映射层（Task B）。

本模块只做**结构表示与一致性检查**，不优化序列、不模拟反应、不生成真实碱基
（占位模式下所有序列都是 ``N``）。

核心映射（已确认的模型，见 docs/01_结构确认_v2.md）：

* scaffold 被 13 nt 一块地铺满，相邻两块交替属于 A 组与 B 组 staple；
* 每条 staple 的 13 nt 与 scaffold 上那 13 nt 反向互补，26 nt 与目标上 26 nt 反向互补；
* 装配体中 scaffold 物理方向 3′→5′，目标 5′→3′，所以 staple 的两个 domain 反向平行，
  折返在左端 → staple 序列 = ``revcomp(目标26) + revcomp(scaffold13)``。

一条硬规矩：**A/B 交界不是 nick**。它应当被某条 boundary-spanning oligo 连续覆盖，
此时交界处的键是"同一条合成链内部"的磷酸二酯键，而不是连接反应的产物。
是否真的连成环，看 ``Junction.bond_state`` 与全体 ``Nick.intended_ligation``，不看坐标是否相邻。
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from .models import BindingSite, Junction, Nick, Staple, TargetOligo, TargetRegion
from .sequence import get_circular_subsequence, reverse_complement

__all__ = [
    "PlanConfig",
    "AssemblyPlan",
    "plan_assembly",
    "validate_plan",
    "to_structure_json",
    "toy_plan",
]

PLACEHOLDER = "N"


@dataclass
class PlanConfig:
    """一份装配计划的输入参数。全部由用户给出，模块不做假设。"""

    scaffold_length: int
    covered_start: int  # 覆盖区在 scaffold 存储坐标里的起点
    covered_length: int  # C = 13 × n_staples
    n_staples: int
    target_A_length: int
    target_B_length: int
    target_topology: str  # "circular" | "linear"
    # 目标坐标下的 oligo 划分；end > length 表示跨原点。必须无缝、无重叠。
    oligo_intervals: list[tuple[int, int]] = field(default_factory=list)
    group_pattern: str = "AB"  # 按 13 nt 块依次循环的分组模式
    staple_domain_order: str = "5p_target_26"
    scaffold_arm: int = 13
    target_arm: int = 26
    placeholder: bool = True
    five_prime_phosphate: str = "unknown"
    intended_ligation: bool = False
    # 非占位模式下由调用方提供（本模块只做映射，不设计序列）
    scaffold_sequence: str | None = None
    target_sequence: str | None = None

    @property
    def target_length(self) -> int:
        return self.target_A_length + self.target_B_length

    @property
    def block(self) -> int:
        return self.scaffold_arm


@dataclass(frozen=True)
class Check:
    """一条检查结论。level="error" 表示结构不自洽，必须由用户决断。"""

    level: str  # "error" | "note"
    code: str
    message: str


@dataclass
class AssemblyPlan:
    config: PlanConfig
    regions: list[TargetRegion] = field(default_factory=list)
    staples: list[Staple] = field(default_factory=list)
    oligos: list[TargetOligo] = field(default_factory=list)
    nicks: list[Nick] = field(default_factory=list)
    junctions: list[Junction] = field(default_factory=list)
    checks: list[Check] = field(default_factory=list)

    # -- 便捷视图 ---------------------------------------------------------- #
    @property
    def target_length(self) -> int:
        return self.config.target_length

    @property
    def covered_length(self) -> int:
        return self.config.covered_length

    @property
    def is_circular(self) -> bool:
        return self.config.target_topology == "circular"

    @property
    def errors(self) -> list[Check]:
        return [c for c in self.checks if c.level == "error"]

    def staples_of(self, group: str) -> list[Staple]:
        return [s for s in self.staples if s.group == group]

    def boundary_oligos(self) -> list[TargetOligo]:
        return [o for o in self.oligos if o.oligo_type == "boundary_spanning"]

    @property
    def circle_closed(self) -> bool:
        """目标链是否已经**共价**闭合成环。

        需要同时满足：拓扑是环、且每一个 nick 都被显式标为已连接。
        序列覆盖连续（含 boundary oligo 跨过交界）**不算**共价闭环。
        """
        if not self.is_circular or not self.nicks:
            return False
        return all(nick.intended_ligation for nick in self.nicks)


# --------------------------------------------------------------------------- #
# 映射
# --------------------------------------------------------------------------- #
def _scaffold_block(
    index: int, cfg: PlanConfig
) -> BindingSite:
    """第 index 块 13 nt 在 scaffold 存储坐标里的位置。

    绘图坐标 p（0 = 覆盖区最左端，向右递增）对应存储坐标 ``covered_start + C - 1 - p``，
    因为装配体中 scaffold 的物理方向是 3′→5′。
    """
    p_start = index * cfg.block
    raw_start = cfg.covered_start + cfg.covered_length - p_start - cfg.block
    start = raw_start % cfg.scaffold_length
    return BindingSite(
        molecule="scaffold",
        start=start,
        end=start + cfg.block,
        strand="-",  # 与 staple 的 13 nt 反向平行
        molecule_length=cfg.scaffold_length,
    )


def _target_arm(
    group: str, k: int, cfg: PlanConfig
) -> BindingSite:
    offset = 0 if group == "A" else cfg.target_A_length
    start = offset + k * cfg.target_arm
    return BindingSite(
        molecule="target",
        start=start,
        end=start + cfg.target_arm,
        strand="-",  # 与 staple 的 26 nt 反向平行
        molecule_length=cfg.target_length,
    )



def _compose_staple(
    cfg: PlanConfig, scaffold_site: BindingSite, target_site: BindingSite
) -> str:
    """按 5′→3′ 拼出 staple：5′-[26 nt 目标臂]-[13 nt scaffold 臂]-3′。

    两个臂都是各自伙伴序列的反向互补；两段之间是**一个**磷酸二酯键，
    不是两条 oligo，也不是 overlap。
    """
    if not cfg.scaffold_sequence or not cfg.target_sequence:
        raise ValueError("非占位模式必须同时提供 scaffold_sequence 与 target_sequence")

    scaffold_piece = get_circular_subsequence(
        cfg.scaffold_sequence, scaffold_site.start, scaffold_site.length
    )
    target_piece = get_circular_subsequence(
        cfg.target_sequence, target_site.start, target_site.length
    )
    target_arm = reverse_complement(target_piece)
    scaffold_arm = reverse_complement(scaffold_piece)
    if cfg.staple_domain_order == "5p_scaffold_13":
        return scaffold_arm + target_arm
    return target_arm + scaffold_arm


def _placeholder(total: int) -> str:
    return PLACEHOLDER * total


def plan_assembly(cfg: PlanConfig) -> AssemblyPlan:
    """按配置生成 staple / oligo / nick / junction 的结构（不含真实序列）。"""
    if cfg.target_topology not in ("circular", "linear"):
        raise ValueError(f"未知目标拓扑：{cfg.target_topology!r}")
    if not cfg.group_pattern:
        raise ValueError("group_pattern 不能为空")

    plan = AssemblyPlan(config=cfg)

    # --- 目标区域 ---
    plan.regions = [
        TargetRegion(
            region_id="target_A",
            group="A",
            start=0,
            end=cfg.target_A_length,
            strand="+",
            sequence=_placeholder(cfg.target_A_length),
            placeholder=cfg.placeholder,
        ),
        TargetRegion(
            region_id="target_B",
            group="B",
            start=cfg.target_A_length,
            end=cfg.target_length,
            strand="+",
            sequence=_placeholder(cfg.target_B_length),
            placeholder=cfg.placeholder,
        ),
    ]

    # --- staple ---
    counters = {"A": 0, "B": 0}
    for index in range(cfg.n_staples):
        group = cfg.group_pattern[index % len(cfg.group_pattern)]
        if group not in ("A", "B"):
            raise ValueError(f"group_pattern 只能含 A/B，遇到 {group!r}")
        k = counters[group]
        counters[group] += 1

        scaffold_site = _scaffold_block(index, cfg)
        target_site = _target_arm(group, k, cfg)
        if cfg.placeholder:
            sequence = _placeholder(cfg.scaffold_arm + cfg.target_arm)
        else:
            sequence = _compose_staple(
                cfg,
                scaffold_site,
                target_site,
            )

        plan.staples.append(
            Staple(
                staple_id=f"S{index + 1:02d}",
                group=group,
                scaffold_binding_region=scaffold_site,
                target_binding_region=target_site,
                sequence=sequence,
                domain_order="5p_target_26",
                placeholder=cfg.placeholder,
            )
        )

    # --- target oligo 划分 ---
    junction_positions = [cfg.target_A_length] if cfg.target_A_length else []
    if cfg.target_topology == "circular":
        junction_positions.append(0)

    for i, (start, end) in enumerate(cfg.oligo_intervals):
        if end <= start:
            raise ValueError(f"oligo 区间 [{start}, {end}) 非法")
        segments = [
            region
            for region in plan.regions
            if _overlaps(start, end, region.start, region.end, cfg.target_length)
        ]
        start_group = "A" if start % cfg.target_length < cfg.target_A_length else "B"
        kind = "ordinary"
        for position in junction_positions:
            if _contains(start, end, position, cfg.target_length):
                kind = "boundary_spanning"
        plan.oligos.append(
            TargetOligo(
                oligo_id=f"T{i + 1:02d}",
                group=start_group,
                oligo_type=kind,
                sequence=_placeholder(end - start),
                target_start=start,
                target_end=end,
                target_segments=segments,
                strand="+",
                five_prime_modification=cfg.five_prime_phosphate,
                three_prime_modification="OH",
                placeholder=cfg.placeholder,
            )
        )

    # --- junction ---
    for ordinal, position in enumerate(junction_positions, start=1):
        covered_by = next(
            (
                oligo.oligo_id
                for oligo in plan.oligos
                if _contains(oligo.start, oligo.end, position, cfg.target_length)
            ),
            None,
        )
        plan.junctions.append(
            Junction(
                junction_id=f"J{ordinal}",
                position=position % cfg.target_length,
                upstream_group="A",
                downstream_group="B",
                upstream_end="3'",
                downstream_end="5'",
                covered_by=covered_by,
                bond_state="within_oligo" if covered_by else "nick",
            )
        )

    # --- nick：按目标坐标顺序找相邻 oligo ---
    ordered = sorted(plan.oligos, key=lambda o: o.start)
    for i, upstream in enumerate(ordered):
        nxt = ordered[(i + 1) % len(ordered)]
        if not plan.is_circular and i == len(ordered) - 1:
            continue
        position = upstream.end % cfg.target_length
        plan.nicks.append(
            Nick(
                nick_id=f"N{i + 1}",
                upstream_oligo=upstream.oligo_id,
                downstream_oligo=nxt.oligo_id,
                genomic_or_target_position=position,
                five_prime_phosphate=cfg.five_prime_phosphate,
                three_prime_hydroxyl="yes",
                intended_ligation=cfg.intended_ligation,
            )
        )

    plan.checks = validate_plan(plan)
    return plan


# --------------------------------------------------------------------------- #
# 区间工具（目标坐标是环状的）
# --------------------------------------------------------------------------- #
def _normalize_interval(start: int, end: int, length: int) -> list[tuple[int, int]]:
    """把可能跨原点的区间拆成 [0, length) 内的若干段。"""
    if length <= 0:
        return [(start, end)]
    pieces: list[tuple[int, int]] = []
    start %= length
    end %= length
    if start < end and end - start == 0:
        return [(0, length)]
    if start < end:
        return [(start, end)]
    # 跨原点
    pieces.append((start, length))
    pieces.append((0, end))
    return pieces


def _overlaps(a_start: int, a_end: int, b_start: int, b_end: int, length: int) -> bool:
    for s1, e1 in _normalize_interval(a_start, a_end, length):
        for s2, e2 in _normalize_interval(b_start, b_end, length):
            if max(s1, s2) < min(e1, e2):
                return True
    return False


def _contains(start: int, end: int, position: int, length: int) -> bool:
    """position 是否严格落在区间内部（环状语义，端点不算）。

    用"相对偏移"判断，避免把跨原点区间在 0 点处误判成端点：
    区间 [350, 420) 在长度 390 的环上覆盖 350..390 与 0..30，
    位置 0 的偏移是 40，严格在 0..70 之间。
    """
    span = (end - start) % length or length
    offset = (position - start) % length
    return 0 < offset < span


# --------------------------------------------------------------------------- #
# 一致性检查
# --------------------------------------------------------------------------- #
def validate_plan(plan: AssemblyPlan) -> list[Check]:
    """返回检查结论清单。只报告，不修复；error 级别必须由用户决断。"""
    cfg = plan.config
    out: list[Check] = []

    def err(code: str, message: str) -> None:
        out.append(Check("error", code, message))

    def note(code: str, message: str) -> None:
        out.append(Check("note", code, message))

    # 1) domain 长度
    for staple in plan.staples:
        if staple.scaffold_domain_length != cfg.scaffold_arm:
            err(
                "scaffold_domain_length",
                f"{staple.staple_id} 的 scaffold domain 是 "
                f"{staple.scaffold_domain_length} nt，预期 {cfg.scaffold_arm}",
            )
        if staple.target_domain_length != cfg.target_arm:
            err(
                "target_domain_length",
                f"{staple.staple_id} 的 target domain 是 "
                f"{staple.target_domain_length} nt，预期 {cfg.target_arm}",
            )

    # 2) 覆盖长度与 staple 数
    expected_covered = cfg.n_staples * cfg.block
    if cfg.covered_length != expected_covered:
        err(
            "covered_length",
            f"covered_length={cfg.covered_length} 与 13×n={expected_covered} 不一致",
        )
    if cfg.covered_length > cfg.scaffold_length:
        err(
            "covered_too_long",
            f"覆盖长度 {cfg.covered_length} 超过 scaffold 长度 {cfg.scaffold_length}",
        )

    # 3) 每组目标长度必须是 26 的整数倍，且与 staple 条数匹配
    for group, region_length in (("A", cfg.target_A_length), ("B", cfg.target_B_length)):
        count = len(plan.staples_of(group))
        expected = count * cfg.target_arm
        if region_length != expected:
            err(
                "group_length",
                f"{group} 组长度 {region_length} nt ≠ 26 × {count} = {expected} nt"
                "（26 nt 臂若首尾相接平铺，两者必须相等）",
            )

    # 4) 严格交替的分组模式 ⇒ 两组条数只差 1
    if set(cfg.group_pattern) == {"A", "B"} and len(cfg.group_pattern) == 2:
        n_a = len(plan.staples_of("A"))
        n_b = len(plan.staples_of("B"))
        if abs(n_a - n_b) > 1:
            err(
                "group_pattern",
                f"严格交替模式下 n_A={n_a}、n_B={n_b} 相差超过 1，"
                "两组长度无法同时满足；请确认分组模式或两组长度",
            )

    # 5) oligo 必须无缝、无重叠地铺满目标链
    ordered = sorted(plan.oligos, key=lambda o: o.start)
    length = cfg.target_length
    for i, oligo in enumerate(ordered):
        if not plan.is_circular and i == len(ordered) - 1:
            continue
        nxt = ordered[(i + 1) % len(ordered)]
        expected_next = oligo.end % length
        if nxt.start != expected_next:
            err(
                "oligo_gap",
                f"{oligo.oligo_id} 结束于 {oligo.end}，但下一条 "
                f"{nxt.oligo_id} 起于 {nxt.start}：存在间隙或重叠",
            )
    if plan.is_circular:
        covered = sum(o.length for o in plan.oligos)
        if covered != length:
            err(
                "oligo_coverage",
                f"环状目标链上 oligo 总长 {covered} ≠ 目标长度 {length}",
            )

    # 6) nick 不得落在 staple 臂接缝上；每条 oligo 至少跨过一个接缝
    seams = _staple_seams(plan)
    for nick in plan.nicks:
        if nick.position in seams:
            err(
                "nick_on_seam",
                f"{nick.nick_id} 落在 staple 臂接缝 {nick.position} 处"
                "（目标片段的 nick 不能与 staple 间的空隙重合）",
            )
    for oligo in plan.oligos:
        if not any(_contains(oligo.start, oligo.end, seam, length) for seam in seams):
            err(
                "oligo_spans_no_seam",
                f"{oligo.oligo_id} 没有跨过任何 staple 臂接缝"
                "（目标片段应 cover 至少一个两 staple 间的缝隙）",
            )

    # 7) 每个 junction 必须被恰好一条 oligo 连续覆盖
    for junction in plan.junctions:
        covering = [
            o.oligo_id
            for o in plan.oligos
            if _contains(o.start, o.end, junction.position, length)
        ]
        if len(covering) != 1:
            err(
                "junction_coverage",
                f"{junction.junction_id} 被 {len(covering)} 条 oligo 覆盖"
                "（交界需要且只需要一条 boundary-spanning oligo）",
            )
        elif junction.covered_by != covering[0]:
            err(
                "junction_covered_by",
                f"{junction.junction_id} 记录的覆盖者 {junction.covered_by} "
                f"与实际 {covering[0]} 不一致",
            )

    # 8) 提示级
    if cfg.placeholder:
        note(
            "placeholder_sequences",
            "占位模式：所有序列都是 N，互补性与热力学未评估（本任务只验证结构表示）",
        )
    if plan.is_circular and not plan.circle_closed:
        note(
            "not_covalently_closed",
            "目标链是环状但尚未共价闭合：所有 nick 的 intended_ligation 都为 False"
            "（序列被连续覆盖 ≠ 已经连成环）",
        )

    return out


def _staple_seams(plan: AssemblyPlan) -> set[int]:
    """所有 staple 臂在目标坐标上的接缝位置。"""
    seams: set[int] = set()
    for group in ("A", "B"):
        for staple in plan.staples_of(group):
            seams.add(staple.target_binding.start % plan.target_length)
            seams.add(staple.target_binding.end % plan.target_length)
    return seams


# --------------------------------------------------------------------------- #
# 结构关系 JSON
# --------------------------------------------------------------------------- #
def to_structure_json(plan: AssemblyPlan) -> dict[str, Any]:
    cfg = plan.config
    return {
        "schema": "scaffold_staple_design/assembly_structure/v1",
        "placeholder_sequences": cfg.placeholder,
        "assumptions": {
            "staple_domain_order": "5p_target_26（26 nt 在 5′ 侧，13 nt 在 3′ 侧）",
            "scaffold_physical_orientation": "装配体中 scaffold 为 3′→5′，目标为 5′→3′",
            "junction_is_not_a_nick": "A/B 交界由 boundary-spanning oligo 连续覆盖",
        },
        "target": {
            "topology": cfg.target_topology,
            "length": cfg.target_length,
            "regions": [
                {
                    "region_id": r.region_id,
                    "group": r.group,
                    "start": r.start,
                    "end": r.end,
                    "length": r.length,
                    "strand": r.strand,
                    "placeholder": r.placeholder,
                }
                for r in plan.regions
            ],
            "oligo_order": [o.oligo_id for o in sorted(plan.oligos, key=lambda o: o.start)],
        },
        "scaffold": {
            "length": cfg.scaffold_length,
            "topology": "circular",
            "covered_start": cfg.covered_start,
            "covered_length": cfg.covered_length,
            "uncovered_length": cfg.scaffold_length - cfg.covered_length,
            "covered_arc": [
                cfg.covered_start,
                (cfg.covered_start + cfg.covered_length) % cfg.scaffold_length,
            ],
        },
        "staples": [
            {
                "staple_id": s.staple_id,
                "group": s.group,
                "domain_order": s.domain_order,
                "total_length": s.total_length,
                "scaffold_start": s.scaffold_start,
                "scaffold_end": s.scaffold_binding.end,
                "scaffold_length": s.scaffold_domain_length,
                "scaffold_wraps": s.scaffold_binding.wraps,
                "scaffold_strand": s.scaffold_binding.strand,
                "target_start": s.target_start,
                "target_end": s.target_end,
                "target_length": s.target_domain_length,
                "target_strand": s.target_binding.strand,
                "placeholder": s.placeholder,
            }
            for s in plan.staples
        ],
        "staples_by_group": {
            "A": [s.staple_id for s in plan.staples_of("A")],
            "B": [s.staple_id for s in plan.staples_of("B")],
        },
        "oligos": [
            {
                "oligo_id": o.oligo_id,
                "group": o.group,
                "oligo_type": o.oligo_type,
                "target_start": o.start,
                "target_end": o.end,
                "length": o.length,
                "strand": o.strand,
                "regions": [r.region_id for r in o.target_segments],
                "five_prime_modification": o.five_prime_modification,
                "three_prime_modification": o.three_prime_modification,
            }
            for o in sorted(plan.oligos, key=lambda o: o.start)
        ],
        "boundary_spanning_oligos": [o.oligo_id for o in plan.boundary_oligos()],
        "junctions": [
            {
                "junction_id": j.junction_id,
                "position": j.position,
                "upstream_group": j.upstream_group,
                "upstream_end": j.upstream_end,
                "downstream_group": j.downstream_group,
                "downstream_end": j.downstream_end,
                "covered_by": j.covered_by,
                "bond_state": j.bond_state,
            }
            for j in plan.junctions
        ],
        "nicks": [
            {
                "nick_id": n.nick_id,
                "upstream_oligo": n.upstream_oligo,
                "upstream_end": n.upstream_end,
                "downstream_oligo": n.downstream_oligo,
                "downstream_end": n.downstream_end,
                "genomic_or_target_position": n.position,
                "five_prime_phosphate": n.five_prime_phosphate,
                "three_prime_hydroxyl": n.three_prime_hydroxyl,
                "intended_ligation": n.intended_ligation,
                "bonded": n.bonded,
            }
            for n in plan.nicks
        ],
        "checks": [
            {"level": c.level, "code": c.code, "message": c.message}
            for c in plan.checks
        ],
        "error_count": len(plan.errors),
        "circle_closed": plan.circle_closed,
    }


# --------------------------------------------------------------------------- #
# Toy example（不含真实序列）
# --------------------------------------------------------------------------- #
def toy_plan() -> AssemblyPlan:
    """15 条 staple 的最小完整例子：8 条 A + 7 条 B，环状目标 390 nt。

    选 15 条是因为它同时暴露三件事：
    严格交替 ⇒ n_A 与 n_B 只差 1（8/7），于是 L(A)=208、L(B)=182；
    交界 208 与 0 各被一条 boundary-spanning oligo 覆盖；
    所有序列都是占位符，互补性不参与判定。
    """
    n = 15
    cfg = PlanConfig(
        scaffold_length=1200,
        covered_start=0,
        covered_length=13 * n,
        n_staples=n,
        target_A_length=26 * 8,
        target_B_length=26 * 7,
        target_topology="circular",
        oligo_intervals=[
            (30, 70),
            (70, 140),
            (140, 190),
            (190, 226),  # 跨过 J1（位置 208）
            (226, 280),
            (280, 350),
            (350, 420),  # 跨过 J2（位置 0，即 390）
        ],
        group_pattern="AB",
        five_prime_phosphate="unknown",
    )
    return plan_assembly(cfg)


def main(argv: list[str] | None = None) -> int:
    import argparse

    parser = argparse.ArgumentParser(description="生成结构关系 JSON（toy example）")
    parser.add_argument("-o", "--out", default="outputs/structure_example.json")
    args = parser.parse_args(argv)

    plan = toy_plan()
    payload = to_structure_json(plan)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")

    print(f"[toy] staples={len(plan.staples)} (A={len(plan.staples_of('A'))}, "
          f"B={len(plan.staples_of('B'))}) oligos={len(plan.oligos)} "
          f"nicks={len(plan.nicks)} junctions={len(plan.junctions)}")
    print(f"      target {plan.target_length} nt, covered scaffold "
          f"{plan.covered_length}/{plan.config.scaffold_length} nt")
    for check in plan.checks:
        prefix = "ERROR" if check.level == "error" else "note "
        print(f"  {prefix} [{check.code}] {check.message}")
    print(f"  -> {out}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
