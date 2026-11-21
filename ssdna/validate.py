"""独立的验证模块（Task D）。

**独立性**：本模块只读文件（``target.fasta`` / ``scaffold_working.fasta`` / 三张 CSV /
``design_report.json``），不 import ``ssdna.design``，不依赖生成器的任何内部状态。
换一套别人的输出，只要目录结构一样，也能验。

**不做的事**：不修数据、不补序列、不把"坐标相邻"当成"共价连接"、不声称能验证连接化学。
热力学部分只生成标准输入文件与手动运行说明；没跑就明确写 ``not_run``，绝不编结果。
"""

from __future__ import annotations

import csv
import hashlib
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable

from .sequence import normalize_sequence, reverse_complement

__all__ = [
    "Issue",
    "Bundle",
    "load_bundle",
    "validate_sequence",
    "validate_length",
    "validate_reverse_complement",
    "validate_target_coverage",
    "validate_group_order",
    "validate_boundary_spanning_oligos",
    "validate_scaffold_mapping",
    "validate_nicks",
    "run_validation",
]


@dataclass(frozen=True)
class Issue:
    """一条检查结论。错误必须能指到 oligo、坐标、类型和修法。"""

    level: str  # "error" | "warning" | "info"
    code: str
    subject: str
    coordinates: str
    message: str
    suggestion: str = ""

    def to_dict(self) -> dict[str, str]:
        return {
            "level": self.level,
            "code": self.code,
            "subject": self.subject,
            "coordinates": self.coordinates,
            "message": self.message,
            "suggestion": self.suggestion,
        }


@dataclass
class Bundle:
    """验证所需的全部输入（全部来自磁盘）。"""

    directory: Path
    target_sequence: str
    scaffold_sequence: str
    oligos: list[dict[str, Any]] = field(default_factory=list)
    staples: list[dict[str, Any]] = field(default_factory=list)
    nicks: list[dict[str, Any]] = field(default_factory=list)
    reported_nicks: list[dict[str, Any]] = field(default_factory=list)
    junctions: list[dict[str, Any]] = field(default_factory=list)
    report: dict[str, Any] = field(default_factory=dict)
    files: dict[str, str] = field(default_factory=dict)


# --------------------------------------------------------------------------- #
# 读取（只读文件）
# --------------------------------------------------------------------------- #
def _read_fasta(path: Path) -> str:
    chunks: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith(">"):
            continue
        chunks.append(line)
    if not chunks:
        raise ValueError(f"{path} 里没有序列")
    return normalize_sequence("".join(chunks))


def _read_csv(path: Path) -> list[dict[str, Any]]:
    with path.open(encoding="utf-8-sig") as fh:
        return list(csv.DictReader(fh))


def _require_columns(path: Path, rows: list[dict[str, Any]], required: list[str]) -> None:
    """CSV 缺列时给出清晰错误，而不是让下游抛 KeyError。"""
    if not rows:
        return  # 空表由各自的检查逻辑处理
    missing = [c for c in required if c not in rows[0]]
    if missing:
        raise ValueError(
            f"{path} 缺少必需列 {missing}；期望列：{required}"
        )


def load_bundle(directory: str | Path) -> Bundle:
    directory = Path(directory)
    required = {
        "target": directory / "target.fasta",
        "scaffold": directory / "scaffold_working.fasta",
        "oligos": directory / "target_oligos.csv",
        "staples": directory / "staples.csv",
    }
    missing = [name for name, path in required.items() if not path.exists()]
    if missing:
        raise FileNotFoundError(
            f"{directory} 缺少验证所需的输入：{missing}；"
            "验证模块只读文件，不会去问生成器"
        )

    target_sequence = _read_fasta(required["target"])
    scaffold_sequence = _read_fasta(required["scaffold"])
    oligos = _read_csv(required["oligos"])
    staples = _read_csv(required["staples"])

    _require_columns(required["oligos"], oligos,
                     ["oligo_id", "sequence", "target_start", "target_end",
                      "oligo_type", "group", "strand", "length"])
    _require_columns(required["staples"], staples,
                     ["staple_id", "sequence", "scaffold_start", "target_start"])

    nicks_path = directory / "nicks.csv"
    nicks = _read_csv(nicks_path) if nicks_path.exists() else []
    _require_columns(nicks_path, nicks,
                     ["nick_id", "upstream_oligo", "downstream_oligo",
                      "genomic_or_target_position"])

    report: dict[str, Any] = {}
    report_path = directory / "design_report.json"
    if report_path.exists():
        report = json.loads(report_path.read_text(encoding="utf-8"))

    files = {
        name: hashlib.sha256(path.read_bytes()).hexdigest()
        for name, path in {
            "target": required["target"],
            "scaffold": required["scaffold"],
            "oligos": required["oligos"],
            "staples": required["staples"],
            "nicks": nicks_path,
        }.items()
        if path.exists()
    }

    return Bundle(
        directory=directory,
        target_sequence=target_sequence,
        scaffold_sequence=scaffold_sequence,
        oligos=oligos,
        staples=staples,
        nicks=nicks,
        reported_nicks=report.get("nicks", []),
        junctions=report.get("junctions", []),
        report=report,
        files=files,
    )


# --------------------------------------------------------------------------- #
# 一、基础序列检查
# --------------------------------------------------------------------------- #
def validate_sequence(sequence: str, *, subject: str = "", coordinates: str = "") -> list[Issue]:
    issues: list[Issue] = []
    if not sequence:
        issues.append(Issue("error", "empty_sequence", subject, coordinates,
                            "序列为空", "检查该 oligo 的源区间是否为空"))
        return issues
    if sequence != sequence.upper():
        issues.append(Issue("warning", "lowercase", subject, coordinates,
                            "序列含小写字符", "统一写成大写 5′→3′"))
    bad = sorted({ch for ch in sequence.upper() if ch not in "ACGT"})
    if bad:
        issues.append(Issue("error", "non_acgt", subject, coordinates,
                            f"含非 ACGT 字符 {bad}", "确认是否混入了简并码或格式字符"))
    return issues


def validate_length(sequence: str, expected_length: int, *,
                    subject: str = "", coordinates: str = "") -> list[Issue]:
    if len(sequence) == expected_length:
        return []
    return [Issue("error", "length_mismatch", subject, coordinates,
                  f"长度 {len(sequence)} nt，期望 {expected_length} nt"
                  f"（差 {len(sequence) - expected_length:+d} nt）",
                  "检查源区间与拼接，不要靠补碱基凑长度")]


def validate_reverse_complement(seq1: str, seq2: str, *,
                                subject: str = "", coordinates: str = "") -> list[Issue]:
    """检查 seq2 是否为 seq1 的反向互补。"""
    if reverse_complement(seq1.upper()) == seq2.upper():
        return []
    # 常见错误：写成了同向互补（忘记反向）
    if seq1.upper().translate(str.maketrans("ACGT", "TGCA")) == seq2.upper():
        return [Issue("error", "same_orientation_complement", subject, coordinates,
                      "互补但方向没反过来（写成了同向互补）",
                      "按 5′→3′ 写成 reverse_complement，而不是直接替换碱基")]
    if seq1.upper() == seq2.upper():
        return [Issue("error", "identical_not_complementary", subject, coordinates,
                      "两段序列完全相同，不是互补关系", "检查是否把模板链当成了配对链")]
    return [Issue("error", "not_reverse_complement", subject, coordinates,
                  "两段序列不构成反向互补", "逐位比对，确认模板区间没有选错")]


# --------------------------------------------------------------------------- #
# 二、目标序列覆盖
# --------------------------------------------------------------------------- #
def validate_target_coverage(
    target_sequence: str,
    oligos: list[dict[str, Any]],
    *,
    topology: str = "circular",
) -> tuple[list[Issue], dict[str, Any]]:
    issues: list[Issue] = []
    length = len(target_sequence)
    if not oligos:
        return [Issue("error", "no_oligos", "-", "-", "没有任何 oligo", "先生成 oligo")], {}

    counts = [0] * length
    owner: list[str] = [""] * length
    for oligo in oligos:
        start = int(oligo["target_start"])
        end = int(oligo["target_end"])
        span = (end - start) % length or length
        for offset in range(span):
            pos = (start + offset) % length
            counts[pos] += 1
            if not owner[pos]:
                owner[pos] = oligo["oligo_id"]

    gaps = [i for i, c in enumerate(counts) if c == 0]
    overlaps = [i for i, c in enumerate(counts) if c > 1]

    def _runs(indices: list[int]) -> str:
        if not indices:
            return "-"
        parts: list[str] = []
        run_start = prev = indices[0]
        for idx in indices[1:]:
            if idx == prev + 1:
                prev = idx
                continue
            parts.append(f"[{run_start},{prev + 1})")
            run_start = prev = idx
        parts.append(f"[{run_start},{prev + 1})")
        return ", ".join(parts[:8]) + (" …" if len(parts) > 8 else "")

    if gaps:
        issues.append(Issue("error", "uncovered_bases", "-", _runs(gaps),
                            f"{len(gaps)} nt 没有被任何 oligo 覆盖",
                            "补一条 oligo，或调整切点；不要自动补序列"))
    if overlaps:
        issues.append(Issue("error", "overlapping_bases", "-", _runs(overlaps),
                            f"{len(overlaps)} nt 被多条 oligo 重复覆盖",
                            "调整切点，让相邻 oligo 只共享端点不共享碱基"))

    # 序列一致性：每条 oligo 的序列必须等于目标区间切片
    doubled = target_sequence * 2
    for oligo in oligos:
        start = int(oligo["target_start"])
        span = int(oligo["oligo_length"]) if "oligo_length" in oligo else int(oligo["length"])
        expected = doubled[start : start + span]
        if oligo["sequence"] != expected:
            issues.append(Issue("error", "oligo_sequence_mismatch", oligo["oligo_id"],
                                f"[{start},{start + span})",
                                "oligo 序列与目标序列对应区间不一致",
                                "检查是否被反向互补过，或区间错位"))
        if oligo.get("strand", "+") != "+":
            issues.append(Issue("error", "oligo_wrong_strand", oligo["oligo_id"],
                                f"[{start},{start + span})",
                                f"oligo 方向标记为 {oligo['strand']}，与目标链 5′→3′ 不一致",
                                "target oligo 必须与目标链同向"))

    stats = {
        "target_length": length,
        "covered_bases": length - len(gaps),
        "uncovered_bases": len(gaps),
        "overlapping_bases": len(overlaps),
        "coverage_fraction": round((length - len(gaps)) / length, 6) if length else 0.0,
        "topology": topology,
    }
    return issues, stats


def validate_group_order(
    oligos: list[dict[str, Any]],
    *,
    groups: Iterable[str] = ("A", "B"),
) -> list[Issue]:
    """A 组必须整体在 B 组之前，且各自连续。"""
    issues: list[Issue] = []
    ordered = sorted(oligos, key=lambda o: int(o["target_start"]))
    seen = [o["group"] for o in ordered]
    blocks: list[str] = []
    for group in seen:
        if not blocks or blocks[-1] != group:
            blocks.append(group)
    expected = [g for g in groups]
    if blocks != expected:
        issues.append(Issue("error", "group_order", "-", "-",
                            f"分组块顺序是 {blocks}，期望 {expected}",
                            "按 A 段在前、B 段在后重排切点"))
    # 每组内部必须连续
    for group in expected:
        members = [o for o in ordered if o["group"] == group]
        for a, b in zip(members, members[1:]):
            if int(b["target_start"]) != int(a["target_end"]):
                issues.append(Issue("error", "group_gap", b["oligo_id"],
                                    f"[{a['target_end']},{b['target_start']})",
                                    f"{group} 组内部不连续", "让同组 oligo 首尾相接"))
    return issues


def validate_boundary_spanning_oligos(
    oligos: list[dict[str, Any]],
    boundaries: list[dict[str, Any]],
    *,
    target_length: int | None = None,
) -> list[Issue]:
    issues: list[Issue] = []
    boundary = [o for o in oligos if o.get("oligo_type") == "boundary_spanning"]
    if len(boundary) != len(boundaries):
        issues.append(Issue("error", "boundary_count", "-", "-",
                            f"跨界的 oligo 有 {len(boundary)} 条，交界有 {len(boundaries)} 个",
                            "每条交界需要且只需要一条 boundary-spanning oligo"))
    for junction in boundaries:
        position = int(junction["position"])
        hit = [
            o for o in oligos
            if _strictly_contains(int(o["target_start"]), int(o["target_end"]),
                                  position, target_length)
        ]
        if len(hit) != 1:
            issues.append(Issue("error", "junction_coverage", junction["junction_id"],
                                f"位置 {position}",
                                f"该交界被 {len(hit)} 条 oligo 覆盖",
                                "调整切点使交界严格落在一条 oligo 内部"))
            continue
        oligo = hit[0]
        if oligo.get("oligo_type") != "boundary_spanning":
            issues.append(Issue("error", "junction_mislabeled", oligo["oligo_id"],
                                f"[{oligo['target_start']},{oligo['target_end']})",
                                "该 oligo 跨过交界却没标成 boundary_spanning",
                                "标成 boundary_spanning 并记录跨界位置"))
        # 序列必须连续（长度等于区间长度，中间没有被拼断）
        if len(oligo["sequence"]) != int(oligo["length"]):
            issues.append(Issue("error", "boundary_not_contiguous", oligo["oligo_id"],
                                f"[{oligo['target_start']},{oligo['target_end']})",
                                "跨界 oligo 的序列长度与区间长度不符",
                                "它必须是一条真实的连续链，不能是两条拼起来的"))
    return issues


def _phosphate_state(value, *, hydroxyl: bool = False) -> str:
    """把各种写法归一到 yes / no / unknown。"""
    text = str(value).strip().lower()
    if text in ("yes", "true", "phos", "phosphate", "5'-p", "5p"):
        return "yes"
    if text in ("no", "false", "oh", "hydroxyl", "none", "-"):
        return "no"
    if hydroxyl and text in ("3'-oh", "3oh"):
        return "yes"
    return "unknown"


def _strictly_contains(start: int, end: int, position: int,
                       length: int | None) -> bool:
    if length:
        span = (end - start) % length or length
        offset = (position - start) % length
        return 0 < offset < span
    return start < position < end


# --------------------------------------------------------------------------- #
# 四、Staple 检查
# --------------------------------------------------------------------------- #
def validate_scaffold_mapping(
    scaffold_sequence: str,
    staples: list[dict[str, Any]],
    *,
    scaffold_arm: int = 13,
    target_arm: int = 26,
    domain_order: str = "5p_target_26",
) -> list[Issue]:
    issues: list[Issue] = []
    L = len(scaffold_sequence)
    doubled = scaffold_sequence * 2
    windows: dict[int, str] = {}

    for staple in staples:
        sid = staple["staple_id"]
        sequence = staple["sequence"]
        s_start = int(staple["scaffold_start"])
        t_start = int(staple["target_start"])

        issues.extend(
            validate_length(sequence, scaffold_arm + target_arm, subject=sid,
                            coordinates=f"scaffold {s_start}")
        )
        if domain_order == "5p_target_26":
            s_dom, t_dom = sequence[target_arm:], sequence[:target_arm]
        else:
            s_dom, t_dom = sequence[:scaffold_arm], sequence[scaffold_arm:]

        expected_s = doubled[s_start : s_start + scaffold_arm]
        if len(expected_s) < scaffold_arm:
            issues.append(Issue("error", "scaffold_out_of_range", sid,
                                f"scaffold [{s_start},{s_start + scaffold_arm})",
                                "结合位点越界", "把覆盖区起点往前挪或减少 staple 条数"))
        else:
            for issue in validate_reverse_complement(expected_s, s_dom,
                                                     subject=sid,
                                                     coordinates=f"scaffold [{s_start},{s_start + scaffold_arm})"):
                issues.append(issue)

        if s_start in windows:
            issues.append(Issue("error", "scaffold_site_collision", sid,
                                f"scaffold {s_start}",
                                f"与 {windows[s_start]} 抢同一个结合位点",
                                "调整平铺相位或删掉重复的一条"))
        windows[s_start] = sid

        _ = t_start  # 目标侧由 validate_target_binding 单独核对

    # 结合位点是否互相重叠：按位置占用表判断，**跨原点也算**
    # （只比较"排序后相邻"的区间会漏掉环上末段与首段的重叠）
    occupied: dict[int, str] = {}
    collisions: dict[tuple[str, str], list[int]] = {}
    for staple in sorted(staples, key=lambda s: int(s["scaffold_start"])):
        sid = staple["staple_id"]
        start = int(staple["scaffold_start"])
        for offset in range(scaffold_arm):
            position = (start + offset) % L
            owner = occupied.get(position)
            if owner is None:
                occupied[position] = sid
                continue
            if owner != sid:
                collisions.setdefault((owner, sid), []).append(position)
    for (first_id, second_id), positions in collisions.items():
        issues.append(Issue(
            "error", "scaffold_site_overlap", second_id,
            "scaffold " + ", ".join(str(p) for p in sorted(positions)[:6]),
            f"与 {first_id} 的结合位点重叠 {len(positions)} nt",
            "相邻 staple 的 13 nt 区间不能相交（环状序列的跨原点处也一样）",
        ))
    return issues


def validate_target_binding(
    target_sequence: str,
    staples: list[dict[str, Any]],
    *,
    target_arm: int = 26,
    domain_order: str = "5p_target_26",
) -> list[Issue]:
    """26 nt 臂是否与目标对应区间互补，以及是否有两条 staple 抢同一段目标。"""
    if domain_order not in ("5p_target_26", "5p_scaffold_13"):
        raise ValueError(f"domain_order 非法：{domain_order!r}")
    issues: list[Issue] = []
    T = len(target_sequence)
    doubled = target_sequence * 2
    windows: dict[int, str] = {}

    for staple in staples:
        sid = staple["staple_id"]
        sequence = staple["sequence"]
        t_start = int(staple["target_start"])
        t_dom = sequence[:target_arm] if domain_order == "5p_target_26" else sequence[-target_arm:]
        expected = doubled[t_start : t_start + target_arm]
        if len(expected) < target_arm:
            issues.append(Issue("error", "target_out_of_range", sid,
                                f"target [{t_start},{t_start + target_arm})",
                                "目标结合位点越界", "检查分组边界与臂数"))
            continue
        for issue in validate_reverse_complement(expected, t_dom, subject=sid,
                                                 coordinates=f"target [{t_start},{t_start + target_arm})"):
            issues.append(issue)
        if t_start in windows:
            issues.append(Issue("error", "target_site_collision", sid,
                                f"target {t_start}",
                                f"与 {windows[t_start]} 抢同一个目标区段",
                                "每条 26 nt 臂必须落在互不重叠的目标区段"))
        windows[t_start] = sid

    # 目标臂重叠：位置占用表，跨原点也算（与 scaffold 侧一致）
    occupied: dict[int, str] = {}
    collisions: dict[tuple[str, str], list[int]] = {}
    for staple in sorted(staples, key=lambda s: int(s["target_start"])):
        sid = staple["staple_id"]
        start = int(staple["target_start"])
        for offset in range(target_arm):
            position = (start + offset) % T
            owner = occupied.get(position)
            if owner is None:
                occupied[position] = sid
                continue
            if owner != sid:
                collisions.setdefault((owner, sid), []).append(position)
    for (first_id, second_id), positions in collisions.items():
        issues.append(Issue(
            "error", "target_site_overlap", second_id,
            "target " + ", ".join(str(p) for p in sorted(positions)[:6]),
            f"与 {first_id} 的目标臂重叠 {len(positions)} nt",
            "26 nt 臂首尾相接可以，重叠不行（环状跨原点处也一样）",
        ))
    return issues


# --------------------------------------------------------------------------- #
# 三、Nick 检查
# --------------------------------------------------------------------------- #
def validate_nicks(
    oligos: list[dict[str, Any]],
    nick_records: list[dict[str, Any]],
    *,
    topology: str = "circular",
    target_length: int | None = None,
) -> list[Issue]:
    issues: list[Issue] = []
    ordered = sorted(oligos, key=lambda o: int(o["target_start"]))
    expected = len(ordered) if topology == "circular" else len(ordered) - 1

    if len(nick_records) != expected:
        issues.append(Issue("error", "nick_count", "-", "-",
                            f"有 {len(nick_records)} 条 nick 记录，{topology} 拓扑下应为 {expected}",
                            "每条相邻关系都要单独记录"))

    by_pair = {(r["upstream_oligo"], r["downstream_oligo"]): r for r in nick_records}
    for i, upstream in enumerate(ordered):
        if topology == "linear" and i == len(ordered) - 1:
            continue
        downstream = ordered[(i + 1) % len(ordered)]
        key = (upstream["oligo_id"], downstream["oligo_id"])
        record = by_pair.get(key)
        if record is None:
            issues.append(Issue("warning", "nick_not_recorded", downstream["oligo_id"],
                                f"位置 {upstream['target_end']}",
                                f"{upstream['oligo_id']} 与 {downstream['oligo_id']} 相邻但没记录",
                                "补一条 nick 记录；相邻不等于已连接"))
            continue
        if record.get("upstream_end") != "3'":
            issues.append(Issue("error", "nick_wrong_end", record["nick_id"],
                                f"位置 {record.get('genomic_or_target_position')}",
                                f"上游端写成 {record.get('upstream_end')}，应为 3′",
                                "方向写错会导致连接化学不成立"))
        if record.get("downstream_end") != "5'":
            issues.append(Issue("error", "nick_wrong_end", record["nick_id"],
                                f"位置 {record.get('genomic_or_target_position')}",
                                f"下游端写成 {record.get('downstream_end')}，应为 5′",
                                "方向写错会导致连接化学不成立"))
        position = record.get("genomic_or_target_position")
        if position is not None and target_length:
            expected_pos = int(upstream["target_end"]) % target_length
            if int(position) % target_length != expected_pos:
                issues.append(Issue("error", "nick_position", record["nick_id"],
                                    f"记录 {position}，实际 {expected_pos}",
                                    "nick 坐标与相邻关系不一致", "按上游 3′ 端坐标重记"))
        # 与对应 oligo 的端基修饰交叉核对（两边都明确时才判断）
        upstream_oligo = next(
            (o for o in oligos if o["oligo_id"] == record.get("upstream_oligo")), None
        )
        downstream_oligo = next(
            (o for o in oligos if o["oligo_id"] == record.get("downstream_oligo")), None
        )
        claimed = _phosphate_state(record.get("five_prime_phosphate"))
        if downstream_oligo is not None:
            declared = _phosphate_state(downstream_oligo.get("five_prime_modification"))
            if "unknown" not in (claimed, declared) and claimed != declared:
                issues.append(Issue(
                    "warning", "nick_oligo_end_mismatch", record["nick_id"],
                    f"位置 {record.get('genomic_or_target_position')}",
                    f"nick 记录 5′-P={claimed}，但下游 oligo "
                    f"{downstream_oligo['oligo_id']} 声明为 {declared}",
                    "确认到底做没做磷酸化；两处必须一致",
                ))
        if upstream_oligo is not None:
            declared_3 = _phosphate_state(upstream_oligo.get("three_prime_modification"),
                                          hydroxyl=True)
            claimed_3 = _phosphate_state(record.get("three_prime_hydroxyl"), hydroxyl=True)
            if "unknown" not in (claimed_3, declared_3) and claimed_3 != declared_3:
                issues.append(Issue(
                    "warning", "nick_oligo_end_mismatch", record["nick_id"],
                    f"位置 {record.get('genomic_or_target_position')}",
                    f"nick 记录 3′-OH={claimed_3}，但上游 oligo "
                    f"{upstream_oligo['oligo_id']} 声明为 {declared_3}",
                    "3′ 端状态必须与合成记录一致",
                ))

        phosphate = str(record.get("five_prime_phosphate", "unknown")).lower()
        wants_ligation = str(record.get("intended_ligation")).lower() == "true"
        if phosphate in ("unknown", ""):
            issues.append(Issue("info", "phosphate_unknown", record["nick_id"],
                                f"位置 {position}",
                                "5′ 磷酸化状态未记录（合成 oligo 默认 5′-OH）",
                                "确认是否做过 5′ 磷酸化；不磷酸化就只能是 nick"))
        elif wants_ligation and phosphate in ("no", "false"):
            issues.append(Issue("warning", "ligation_without_phosphate", record["nick_id"],
                                f"位置 {position}",
                                "打算连接，但 5′ 端没有磷酸化",
                                "连接酶需要 5′-P；要么磷酸化，要么接受这条 nick"))
    return issues


# --------------------------------------------------------------------------- #
# 五、热力学接口（只生成输入，不伪造结果）
# --------------------------------------------------------------------------- #
def write_idt_input(bundle: Bundle, path: Path) -> Path:
    """IDT OligoAnalyzer 的批量输入：名字 + 序列，一行一条。"""
    rows = [("oligo_id", "sequence", "note")]
    for oligo in bundle.oligos:
        rows.append((oligo["oligo_id"], oligo["sequence"], "target oligo"))
    for staple in bundle.staples:
        rows.append((staple["staple_id"], staple["sequence"], "staple 39 nt"))
    with path.open("w", encoding="utf-8-sig", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerows(rows)
    return path


def write_nupack_input(bundle: Bundle, path: Path, *, temperature_c: float = 25.0) -> Path:
    """NUPACK 的链列表（``name sequence`` 两列，NUPACK 3/4 及网页版都接受）。"""
    lines = []
    for oligo in bundle.oligos:
        lines.append(f"{oligo['oligo_id']} {oligo['sequence']}")
    for staple in bundle.staples:
        lines.append(f"{staple['staple_id']} {staple['sequence']}")
    path.write_text("\n".join(lines) + "\n", encoding="ascii")
    return path


def _thermo_section(outdir: Path, bundle: Bundle) -> dict[str, Any]:
    idt = write_idt_input(bundle, outdir / "idt_oligoanalyzer_input.csv")
    nupack = write_nupack_input(bundle, outdir / "nupack_strands.txt")
    return {
        "idt_oligoanalyzer": {
            "status": "not_run",
            "input_file": idt.name,
            "what_it_would_check": ["GC content", "Tm", "hairpin", "self-dimer", "hetero-dimer"],
            "how_to_run": [
                "打开 IDT OligoAnalyzer (https://www.idtdna.com/calc/analyzer)",
                "批量模式粘贴 idt_oligoanalyzer_input.csv 的 sequence 列（或逐条粘贴）",
                "参数：Oligo conc 0.25 µM，Na+ 50 mM，Mg++ 0 mM（与本项目配置一致）",
                "导出 Analyze 结果，按 oligo_id 回填到本报告",
            ],
            "caveat": "这些工具只算二级结构与杂交倾向，不能验证连接化学，也不能预测支架组装效率",
        },
        "nupack": {
            "status": "not_run",
            "input_file": nupack.name,
            "what_it_would_check": ["多链竞争性杂交", "替代配对", "复合物平衡浓度"],
            "how_to_run": [
                "nupack 需单独安装（pip install nupack 或官网学术版）",
                "把 nupack_strands.txt 作为链列表，跑 complexes / mfe",
                "建议只用同一组的 staple + 目标片段做小体系，全体系规模太大",
            ],
            "caveat": "NUPACK 不模拟 DNA 连接酶，也不模拟 scaffold 的拓扑约束",
        },
        "not_claimed": "本模块没有运行上述任一工具，报告里没有它们的结果",
    }


# --------------------------------------------------------------------------- #
# 汇总
# --------------------------------------------------------------------------- #
def run_validation(
    bundle: Bundle,
    *,
    topology: str | None = None,
    scaffold_arm: int = 13,
    target_arm: int = 26,
    domain_order: str = "5p_target_26",
) -> dict[str, Any]:
    config = bundle.report.get("config", {})
    topo = topology or config.get("target", {}).get("topology", "circular")
    staple_cfg = config.get("staple", {})
    # 以报告里的实际配置为准，避免"生成器用了 5p_scaffold_13、验证器还按 5p_target_26 拆"
    scaffold_arm = int(staple_cfg.get("scaffold_binding_length") or scaffold_arm)
    target_arm = int(staple_cfg.get("target_binding_length") or target_arm)
    domain_order = str(staple_cfg.get("domain_order") or domain_order)
    if domain_order not in ("5p_target_26", "5p_scaffold_13"):
        raise ValueError(f"配置里 domain_order 非法：{domain_order!r}")
    T = len(bundle.target_sequence)

    issues: list[Issue] = []
    checks: list[dict[str, Any]] = []

    def record(name: str, found: list[Issue]) -> None:
        issues.extend(found)
        checks.append({
            "name": name,
            "status": "fail" if any(i.level == "error" for i in found)
            else ("warn" if found else "pass"),
            "issues": len(found),
        })

    if not bundle.staples:
        issues.append(Issue("error", "no_staples", "-", "-",
                            "staples.csv 为空：没有任何 staple",
                            "重新生成 staple；空设计不应判定为通过"))
    record("sequence_basics", [
        issue
        for item in list(bundle.oligos) + list(bundle.staples)
        for issue in validate_sequence(item["sequence"],
                                       subject=item.get("oligo_id") or item.get("staple_id", ""),
                                       coordinates="-")
    ])

    coverage_issues, coverage_stats = validate_target_coverage(
        bundle.target_sequence, bundle.oligos, topology=topo
    )
    record("target_coverage", coverage_issues)
    record("group_order", validate_group_order(bundle.oligos))
    record("boundary_spanning", validate_boundary_spanning_oligos(
        bundle.oligos, bundle.junctions, target_length=T))
    record("staple_length", [
        issue
        for staple in bundle.staples
        for issue in validate_length(staple["sequence"], scaffold_arm + target_arm,
                                     subject=staple["staple_id"],
                                     coordinates=f"scaffold {staple['scaffold_start']}")
    ])
    record("scaffold_mapping", validate_scaffold_mapping(
        bundle.scaffold_sequence, bundle.staples,
        scaffold_arm=scaffold_arm, target_arm=target_arm, domain_order=domain_order))
    record("target_binding", validate_target_binding(
        bundle.target_sequence, bundle.staples,
        target_arm=target_arm, domain_order=domain_order))
    record("nicks", validate_nicks(
        bundle.oligos, bundle.nicks or bundle.reported_nicks,
        topology=topo, target_length=T))

    counts = {"error": 0, "warning": 0, "info": 0}
    for issue in issues:
        counts[issue.level] = counts.get(issue.level, 0) + 1

    outdir = bundle.directory / "validation"
    outdir.mkdir(parents=True, exist_ok=True)
    thermo = _thermo_section(outdir, bundle)

    payload: dict[str, Any] = {
        "schema": "scaffold_staple_design/validation_report/v1",
        "status": "PASS" if counts["error"] == 0 else "FAIL",
        "validated_directory": str(bundle.directory),
        "input_files": bundle.files,
        "topology": topo,
        "counts": {
            "oligos": len(bundle.oligos),
            "boundary_oligos": sum(
                1 for o in bundle.oligos if o.get("oligo_type") == "boundary_spanning"),
            "staples": len(bundle.staples),
            "nicks": len(bundle.nicks or bundle.reported_nicks),
            "issues": counts,
        },
        "coverage": coverage_stats,
        "checks": checks,
        "issues": [issue.to_dict() for issue in issues],
        "thermo": thermo,
        "independence": (
            "本模块只读 target.fasta / scaffold_working.fasta / 三张 CSV / design_report.json，"
            "不 import 生成器的内部对象"
        ),
    }

    (outdir / "validation_report.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")

    _write_conflicts(issues, outdir / "potential_conflicts.csv")
    _write_summary(payload, outdir / "validation_summary.md", bundle)
    return payload


def _write_conflicts(issues: list[Issue], path: Path) -> Path:
    with path.open("w", encoding="utf-8-sig", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["level", "code", "subject", "coordinates", "message", "suggestion"])
        for issue in issues:
            if issue.level == "info" and issue.code in ("phosphate_unknown",):
                continue
            writer.writerow([issue.level, issue.code, issue.subject,
                             issue.coordinates, issue.message, issue.suggestion])
    return path


def _write_summary(payload: dict[str, Any], path: Path, bundle: Bundle) -> Path:
    lines = [
        "# 验证摘要",
        "",
        f"- 结论：**{payload['status']}**",
        f"- 验证目录：`{payload['validated_directory']}`",
        f"- 拓扑：{payload['topology']}",
        f"- oligo {payload['counts']['oligos']} 条（其中跨界 "
        f"{payload['counts']['boundary_oligos']} 条）、staple {payload['counts']['staples']} 条、"
        f"nick {payload['counts']['nicks']} 条",
        f"- 错误 {payload['counts']['issues']['error']}、警告 "
        f"{payload['counts']['issues']['warning']}、提示 {payload['counts']['issues']['info']}",
        "",
        "## 目标覆盖",
        "",
        "| 项 | 值 |",
        "|---|---|",
    ]
    for key, value in payload["coverage"].items():
        lines.append(f"| {key} | {value} |")

    lines += ["", "## 各项检查", "", "| 检查 | 状态 | 问题数 |", "|---|---|---|"]
    for check in payload["checks"]:
        lines.append(f"| {check['name']} | {check['status']} | {check['issues']} |")

    lines += ["", "## 问题清单", ""]
    if not payload["issues"]:
        lines.append("没有发现问题。")
    else:
        lines += ["| 级别 | 代码 | 对象 | 坐标 | 说明 | 建议 |", "|---|---|---|---|---|---|"]
        for issue in payload["issues"]:
            lines.append(
                f"| {issue['level']} | {issue['code']} | {issue['subject']} | "
                f"{issue['coordinates']} | {issue['message']} | {issue['suggestion']} |"
            )

    lines += ["", "## 热力学工具（未运行）", ""]
    for tool, info in payload["thermo"].items():
        if not isinstance(info, dict):
            continue
        lines.append(f"**{tool}**：`{info['status']}`，输入文件 `{info['input_file']}`")
        lines.append("")
        lines.append(f"- 能查：{', '.join(info['what_it_would_check'])}")
        lines.append(f"- 注意：{info['caveat']}")
        lines.append("")
    lines.append(f"> {payload['thermo']['not_claimed']}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


# --------------------------------------------------------------------------- #
# 命令行
# --------------------------------------------------------------------------- #
def main(argv: list[str] | None = None) -> int:
    import argparse

    parser = argparse.ArgumentParser(description="独立验证 target oligo / staple 设计")
    parser.add_argument("directory", help="设计输出目录（含 target.fasta / *.csv）")
    args = parser.parse_args(argv)

    bundle = load_bundle(args.directory)
    payload = run_validation(bundle)
    print(f"[validate] {payload['status']}  "
          f"errors={payload['counts']['issues']['error']} "
          f"warnings={payload['counts']['issues']['warning']}")
    for issue in payload["issues"]:
        if issue["level"] == "error":
            print(f"  ERROR [{issue['code']}] {issue['subject']} {issue['coordinates']}: "
                  f"{issue['message']}")
    print(f"  -> {Path(args.directory) / 'validation'}")
    return 0 if payload["status"] == "PASS" else 1


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
