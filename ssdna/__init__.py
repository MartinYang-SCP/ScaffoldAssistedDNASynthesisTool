"""scaffold–staple 目标片段组装设计工具。

当前阶段：**输入层（Task A）**。只读文件、建数据模型、出报告，不生成任何设计序列。

公共 API
--------
>>> from ssdna import load_snapgene, export_all
>>> scaffold = load_snapgene("minimalplasmid(1).dna")
>>> scaffold.length, scaffold.topology
(1200, 'circular')
>>> paths = export_all(scaffold, "outputs")
"""

from __future__ import annotations

from .io import BackendUnavailable, load_fasta, load_snapgene
from .models import (
    BindingSite,
    Feature,
    Junction,
    Nick,
    Primer,
    Scaffold,
    Segment,
    Staple,
    TargetOligo,
    TargetRegion,
)
from .report import (
    build_info,
    export_all,
    write_fasta,
    write_features_csv,
    write_info_json,
)
from .sequence import (
    check_sequence,
    gc_content,
    get_circular_subsequence,
    normalize_sequence,
    reverse_complement,
)

__version__ = "0.1.0"

__all__ = [
    "__version__",
    "Scaffold",
    "Feature",
    "Primer",
    "Segment",
    "BindingSite",
    "TargetRegion",
    "TargetOligo",
    "Staple",
    "Nick",
    "Junction",
    "PlanConfig",
    "AssemblyPlan",
    "plan_assembly",
    "toy_plan",
    "to_structure_json",
    "load_snapgene",
    "load_fasta",
    "BackendUnavailable",
    "reverse_complement",
    "normalize_sequence",
    "check_sequence",
    "get_circular_subsequence",
    "gc_content",
    "build_info",
    "export_all",
    "write_fasta",
    "write_features_csv",
    "write_info_json",
]

_MAPPING_EXPORTS = {
    "AssemblyPlan",
    "PlanConfig",
    "plan_assembly",
    "toy_plan",
    "to_structure_json",
}


def __getattr__(name: str):
    """装配层按需导入：避免 ``python -m ssdna.mapping`` 时的重复导入。"""
    if name in _MAPPING_EXPORTS:
        from . import mapping

        return getattr(mapping, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
