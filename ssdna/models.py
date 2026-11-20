"""公共数据模型：scaffold、特征、片段。

三条硬约定（贯穿全项目）：

* 序列一律按 5′→3′ 存储，不存"反向链"的副本；
* 坐标一律 0-based 左闭右开 ``[start, end)``，环状序列允许 ``end > length`` 表示跨原点；
* 区分**文件里的原始数据**（``origin="file"``）与**推断出来的数据**（``origin="derived"``），
  两者绝不混在一起、也不互相覆盖。
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass, field
from typing import Any

__all__ = ["Segment", "Feature", "Primer", "Scaffold"]


@dataclass(frozen=True)
class Segment:
    """一个位置区间。跨原点时 ``wraps=True`` 且 ``end < start``。"""

    start: int
    end: int
    strand: str = "."  # "+" / "-" / "."
    wraps: bool = False
    raw: str = ""  # 文件里的原始写法（SnapGene 是 1-based 闭区间）
    origin: str = "file"  # "file" | "derived"

    def length_on(self, circular_length: int) -> int:
        """按给定环长计算长度（跨原点也成立）。"""
        if self.wraps:
            return (self.end - self.start) % circular_length
        return self.end - self.start

    @property
    def length(self) -> int:
        if self.wraps:
            raise ValueError(
                "跨原点片段的长度依赖环长，请用 length_on(circular_length)"
            )
        return self.end - self.start


@dataclass
class Feature:
    """一条注释（不含 primer）。"""

    feature_type: str
    label: str
    strand: str
    segments: list[Segment] = field(default_factory=list)
    qualifiers: dict[str, str] = field(default_factory=dict)
    kind: str = "annotation"  # "annotation" | "primer"
    origin: str = "file"

    @property
    def note(self) -> str:
        return self.qualifiers.get("note", "")

    @property
    def start(self) -> int | None:
        return self.segments[0].start if self.segments else None

    @property
    def end(self) -> int | None:
        return self.segments[0].end if self.segments else None

    def total_length(self, circular_length: int) -> int:
        return sum(seg.length_on(circular_length) for seg in self.segments)


@dataclass
class Primer:
    """SnapGene 里记录的引物。

    ``sequence`` 来自文件的 primer 块（事实）；``site`` 是**按序列回贴到模板**推出来的
    结合位点（推断），因此 ``site_origin`` 固定为 ``"derived"``。
    """

    name: str
    sequence: str
    site: list[Segment] = field(default_factory=list)
    site_origin: str = "derived"


@dataclass
class Scaffold:
    """scaffold 的输入数据模型。

    ``length`` 必须等于 ``len(sequence)``，构造时校验；不做任何"看起来对"的修正。
    """

    id: str
    sequence: str
    length: int
    topology: str  # "circular" | "linear"
    features: list[Feature] = field(default_factory=list)
    primers: list[Primer] = field(default_factory=list)
    molecule_type: str = "DNA"
    molecule_type_source: str = "default"
    source: str = ""
    sha256: str = ""
    parser_backend: str = ""
    notes: dict[str, str] = field(default_factory=dict)
    sequence_class: str = ""
    cross_check: dict[str, Any] = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)

    def __post_init__(self) -> None:
        if self.length != len(self.sequence):
            raise ValueError(
                f"length={self.length} 与序列实际长度 {len(self.sequence)} 不一致；"
                "拒绝静默修正"
            )
        if self.topology not in ("circular", "linear"):
            raise ValueError(f"未知拓扑：{self.topology!r}（只接受 circular / linear）")

    # -- 便捷视图 ---------------------------------------------------------- #
    @property
    def is_circular(self) -> bool:
        return self.topology == "circular"

    @property
    def gc_content(self) -> float:
        if not self.sequence:
            return 0.0
        gc = sum(base in "GC" for base in self.sequence.upper())
        return gc / len(self.sequence)

    def subsequence(self, start: int, length: int) -> str:
        """取 ``length`` 个碱基；环状序列自动跨原点。"""
        from .sequence import get_circular_subsequence

        return get_circular_subsequence(
            self.sequence, start, length, circular=self.is_circular
        )

    @property
    def sequence_sha256(self) -> str:
        """序列本身的 sha256（判定"还是不是同一条序列"用它，不用文件字节）。"""
        return hashlib.sha256(self.sequence.encode("ascii")).hexdigest()

    def digest(self) -> str:
        return self.sequence_sha256

    def summary(self) -> str:
        return (
            f"{self.id}: {self.length} nt, {self.topology}, "
            f"GC {self.gc_content:.1%}, {len(self.features)} 注释, "
            f"{len(self.primers)} 引物, backend={self.parser_backend}"
        )


# --------------------------------------------------------------------------- #
# 装配层模型（Task B）
#
# 坐标约定与文件特征不同，这里统一为：
#     start ∈ [0, L)，length = end - start，end 可以 >= L（跨原点），wraps = end > L。
# 这样一个区间只需要 (start, length) 两个数就能完整表达，且能直接喂给
# Scaffold.subsequence()。文件特征那套（start > end）是 SnapGene 自己的写法，
# 只在 Segment 里保留。
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class BindingSite:
    """某个分子上的一段连续结合位点。"""

    molecule: str  # "scaffold" | "target"
    start: int  # 0-based 半开
    end: int
    strand: str  # 相对该分子存储的 5′→3′ 方向："+"/"-"
    molecule_length: int

    @property
    def length(self) -> int:
        return self.end - self.start

    @property
    def wraps(self) -> bool:
        return self.end > self.molecule_length

    @property
    def normalized_start(self) -> int:
        return self.start % self.molecule_length

    def __post_init__(self) -> None:
        if self.start < 0 or self.end < self.start:
            raise ValueError(f"非法区间 [{self.start}, {self.end})")


@dataclass
class TargetRegion:
    """目标链上的一段：A 组或 B 组。坐标在**目标分子自身**的 5′→3′ 系里。"""

    region_id: str
    group: str  # "A" | "B"
    start: int
    end: int
    strand: str  # 相对目标分子存储方向，恒为 "+"
    sequence: str = ""
    placeholder: bool = False

    @property
    def length(self) -> int:
        return self.end - self.start


@dataclass
class TargetOligo:
    """实际要合成的一条 target oligo（5′→3′）。"""

    oligo_id: str
    group: str
    oligo_type: str  # "ordinary" | "boundary_spanning"
    sequence: str
    target_start: int
    target_end: int
    target_segments: list[TargetRegion] = field(default_factory=list)
    strand: str = "+"
    five_prime_modification: str = "unknown"
    three_prime_modification: str = "unknown"
    placeholder: bool = False

    # 旧名（结构模型内部计算用），保持可用
    @property
    def start(self) -> int:
        return self.target_start

    @property
    def end(self) -> int:
        return self.target_end

    @property
    def length(self) -> int:
        return self.target_end - self.target_start

    @property
    def spans_junctions(self) -> bool:
        return self.oligo_type == "boundary_spanning"


@dataclass
class Staple:
    """普通 staple：13 nt 贴 scaffold + 26 nt 贴目标，同一条链上连续。"""

    staple_id: str
    group: str  # 服务于哪一组目标："A" | "B"
    scaffold_binding_region: BindingSite
    target_binding_region: BindingSite
    sequence: str
    domain_order: str = "5p_target_26"  # "5p_target_26" | "5p_scaffold_13"
    placeholder: bool = False

    # 旧名（内部计算与既有测试用），保持可用
    @property
    def scaffold_binding(self) -> BindingSite:
        return self.scaffold_binding_region

    @property
    def target_binding(self) -> BindingSite:
        return self.target_binding_region

    @property
    def scaffold_start(self) -> int:
        return self.scaffold_binding_region.start

    @property
    def target_start(self) -> int:
        return self.target_binding_region.start

    @property
    def target_end(self) -> int:
        return self.target_binding_region.end

    @property
    def scaffold_domain_length(self) -> int:
        return self.scaffold_binding_region.length

    @property
    def target_domain_length(self) -> int:
        return self.target_binding_region.length

    @property
    def total_length(self) -> int:
        return self.scaffold_binding_region.length + self.target_binding.length


@dataclass
class Nick:
    """两条相邻 target oligo 之间的相邻关系。默认**不是**共价连接。"""

    nick_id: str
    upstream_oligo: str
    downstream_oligo: str
    genomic_or_target_position: int  # 本项目中固定为目标坐标系
    upstream_end: str = "3'"
    downstream_end: str = "5'"
    five_prime_phosphate: str = "unknown"  # "yes" / "no" / "unknown"
    three_prime_hydroxyl: str = "yes"
    intended_ligation: bool = False

    @property
    def position(self) -> int:
        return self.genomic_or_target_position

    @property
    def bonded(self) -> bool:
        """是否已经共价连接：只有显式 intended_ligation 才算。"""
        return self.intended_ligation


@dataclass
class Junction:
    """A/B 交界。它**不是** nick，而是一个"应当被某条 oligo 连续覆盖"的位置。"""

    junction_id: str
    position: int
    upstream_group: str
    downstream_group: str
    upstream_end: str  # "3'"
    downstream_end: str  # "5'"
    covered_by: str | None = None  # 覆盖它的 boundary-spanning oligo id
    bond_state: str = "nick"  # "within_oligo" | "nick" | "ligated"

    def __post_init__(self) -> None:
        if self.bond_state not in ("within_oligo", "nick", "ligated"):
            raise ValueError(f"未知 bond_state：{self.bond_state!r}")

    @property
    def covered(self) -> bool:
        return self.covered_by is not None
