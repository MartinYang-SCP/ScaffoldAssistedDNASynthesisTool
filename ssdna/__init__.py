"""scaffold–staple 目标片段组装设计工具。

一条约 1.2 kb 的环状单链 scaffold，用统一设计的 39 nt staple
（13 nt 贴 scaffold + 26 nt 贴目标）把由多条 oligo 组成的目标片段定位并闭合成环。

模块分四层：

* 输入层：``io`` / ``snapgene_raw``（SnapGene 双解析器交叉校验）、``sequence``；
* 数据模型：``models``（Scaffold / Feature / TargetRegion / TargetOligo / Staple /
  Nick / Junction）；
* 装配与生成：``mapping``（结构模型）、``design``（序列生成引擎 + 切点建议 + 相位扫描）；
* 验证与报告：``validate``（独立验证）、``thermo``（Tm 体检）、``report``（导出）。

唯一入口：``python main.py --config config/design_config.yaml``。
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

__version__ = "0.2.0"

__all__ = [
    "__version__",
    # 数据模型
    "Scaffold", "Feature", "Primer", "Segment", "BindingSite",
    "TargetRegion", "TargetOligo", "Staple", "Nick", "Junction",
    # 结构模型（mapping）
    "Check", "PlanConfig", "AssemblyPlan", "plan_assembly", "validate_plan",
    "toy_plan", "to_structure_json",
    # 输入
    "load_snapgene", "load_fasta", "BackendUnavailable",
    # 序列
    "reverse_complement", "normalize_sequence", "check_sequence",
    "get_circular_subsequence", "gc_content",
    # 生成
    "design_target_oligos", "design_staples", "build_nicks", "check_design",
    "run_design", "scan_scaffold_phase", "suggest_oligo_cuts",
    "load_config", "find_missing", "find_invalid",
    # 验证
    "validate_plan", "run_validation",
    # 报告
    "build_info", "export_all", "write_fasta", "write_features_csv", "write_info_json",
]

# 按需导入：mapping / design / validate 里都有主流程，避免 `python -m ssdna.x` 时重复导入
_LAZY_EXPORTS = {
    "Check": ("ssdna.mapping", "Check"),
    "PlanConfig": ("ssdna.mapping", "PlanConfig"),
    "AssemblyPlan": ("ssdna.mapping", "AssemblyPlan"),
    "plan_assembly": ("ssdna.mapping", "plan_assembly"),
    "validate_plan": ("ssdna.mapping", "validate_plan"),
    "toy_plan": ("ssdna.mapping", "toy_plan"),
    "to_structure_json": ("ssdna.mapping", "to_structure_json"),
    "design_target_oligos": ("ssdna.design", "design_target_oligos"),
    "design_staples": ("ssdna.design", "design_staples"),
    "build_nicks": ("ssdna.design", "build_nicks"),
    "check_design": ("ssdna.design", "check_design"),
    "run_design": ("ssdna.design", "run_design"),
    "scan_scaffold_phase": ("ssdna.design", "scan_scaffold_phase"),
    "suggest_oligo_cuts": ("ssdna.design", "suggest_oligo_cuts"),
    "load_config": ("ssdna.design", "load_config"),
    "find_missing": ("ssdna.design", "find_missing"),
    "find_invalid": ("ssdna.design", "find_invalid"),
    "run_validation": ("ssdna.validate", "run_validation"),
}


def __getattr__(name: str):
    """按需导入，避免 ``python -m ssdna.mapping`` / ``ssdna.design`` 的重复导入。"""
    if name in _LAZY_EXPORTS:
        import importlib

        module_name, attr = _LAZY_EXPORTS[name]
        module = importlib.import_module(module_name)
        return getattr(module, attr)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
