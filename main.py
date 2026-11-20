#!/usr/bin/env python
"""项目唯一入口：按 8 个步骤跑完整工作流。

    python main.py --config config/design_config.yaml

步骤：
    1 读取 scaffold        2 读取目标序列        3 读取设计参数
    4 建立结构映射          5 生成 target oligos   6 生成 staples
    7 执行验证              8 导出结果

任何关键参数缺失都会**在第 3 步停下**并列出缺失项，不生成、不落盘、不猜。
"""

from __future__ import annotations

import argparse
import json
import pathlib
import sys
from typing import Any

ROOT = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from ssdna import design, validate  # noqa: E402
from ssdna.design import DesignInputs  # noqa: E402
from ssdna.io import load_fasta, load_snapgene  # noqa: E402
from ssdna.sequence import normalize_sequence, reverse_complement  # noqa: E402

BANNER = "=" * 72


def _resolve(path: str, *, base: pathlib.Path) -> pathlib.Path:
    """配置文件里的相对路径：先按当前目录，再按项目根目录。"""
    candidate = pathlib.Path(path)
    if candidate.is_absolute():
        return candidate
    for anchor in (pathlib.Path.cwd(), base, ROOT):
        resolved = (anchor / candidate).resolve()
        if resolved.exists():
            return resolved
    return (base / candidate).resolve()


def _load_target(cfg: dict[str, Any], *, base: pathlib.Path) -> str:
    source = cfg["target"]["source"]
    if isinstance(source, dict) and "inline" in source:
        return normalize_sequence(source["inline"])
    return load_fasta(_resolve(str(source), base=base)).sequence


def _orient(sequence: str, strand: str) -> str:
    if strand == "+":
        return sequence
    if strand == "-":
        return reverse_complement(sequence)
    raise ValueError(f"scaffold.strand 必须是 '+' 或 '-'，收到 {strand!r}")


def run(config_path: pathlib.Path, *, verbose: bool = True) -> int:
    base = config_path.resolve().parent

    def say(text: str) -> None:
        if verbose:
            print(text, flush=True)

    say(BANNER)
    say(f"配置：{config_path}")
    say(BANNER)

    # ---- 步骤 1：读参数并校验（先做，因为 2/3 的路径都来自它） ----
    say("[1/8] 读取设计参数并校验完整性 …")
    cfg = design.load_config(config_path)
    missing = design.find_missing(cfg)
    if missing:
        say(f"  缺少 {len(missing)} 个必需参数，流程在这里停止：")
        for key in missing:
            say(f"    - {key}")
        say("  未读取 scaffold、未读取目标、未生成任何序列、未创建输出目录。")
        return 2
    say("  参数齐备")

    # ---- 步骤 1：scaffold ----
    say("[2/8] 读取 scaffold …")
    scaffold_path = _resolve(str(cfg["scaffold"]["file"]), base=base)
    scaffold = load_snapgene(scaffold_path)
    say(f"  {scaffold.id}: {scaffold.length} nt, {scaffold.topology}, "
        f"backend={scaffold.parser_backend}")

    # ---- 步骤 2：目标序列 ----
    say("[3/8] 读取目标序列 …")
    target_sequence = _load_target(cfg, base=base)
    say(f"  {len(target_sequence)} nt")

    # ---- 步骤 4-6：映射 / oligo / staple ----
    strand = cfg["scaffold"]["strand"]
    inputs = DesignInputs(
        scaffold=scaffold,
        scaffold_strand=strand,
        working_sequence=_orient(scaffold.sequence, strand),
        target_sequence=target_sequence,
        target_length=len(target_sequence),
        config=cfg,
    )

    labels = {
        "build_mapping": "[4/8] 建立结构映射（A/B 分区 + 两个交界）…",
        "design_oligos": "[5/8] 生成 target oligos …",
        "design_staples": "[6/8] 生成 staples …",
        "check_design": "      生成后自检 …",
        "export": "[7/8] 导出结果（CSV / FASTA / JSON）…",
    }

    def progress(name: str) -> None:
        say(labels.get(name, f"  … {name}"))

    outdir = _resolve(str(cfg["output"]["outdir"]), base=base)
    payload = design.run_design(inputs, outdir=outdir, progress=progress)

    for issue in payload["checks"]:
        if issue["level"] == "error":
            say(f"  ERROR [{issue['code']}] {issue['message']}")

    # ---- 步骤 8：独立验证（只读上一步落盘的文件） ----
    say("[8/8] 独立验证（只读上一步的文件）…")
    bundle = validate.load_bundle(outdir)
    validation = validate.run_validation(bundle)
    say(f"  {validation['status']}  errors={validation['counts']['issues']['error']} "
        f"warnings={validation['counts']['issues']['warning']}")
    for issue in validation["issues"]:
        if issue["level"] == "error":
            say(f"  ERROR [{issue['code']}] {issue['subject']} "
                f"{issue['coordinates']}: {issue['message']}")

    # ---- 汇总 ----
    say(BANNER)
    say(f"设计状态 : {payload['status']}")
    say(f"验证状态 : {validation['status']}")
    say(f"oligo    : {payload['counts']['target_oligos']} 条"
        f"（跨界 {payload['counts']['boundary_oligos']} 条）")
    say(f"staple   : {payload['counts']['staples']} 条"
        f"（A {payload['counts']['staples_A']} / B {payload['counts']['staples_B']}）")
    say(f"覆盖     : {payload['counts']['covered_scaffold_nt']} nt scaffold")
    tm = payload.get("tm_summary", {})
    if tm.get("available"):
        for key in ("tm_scaffold_domain", "tm_target_domain"):
            if key in tm:
                s = tm[key]
                say(f"{key:18s}: mean {s['mean']} °C, range [{s['min']}, {s['max']}]")
    say(f"输出     : {outdir}")
    say(f"验证报告 : {outdir / 'validation'}")
    say(BANNER)

    failed = payload["status"] != "OK" or validation["status"] != "PASS"
    return 1 if failed else 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="main.py",
        description="scaffold–staple–target 目标片段组装设计工具（完整工作流）",
    )
    parser.add_argument("--config", default="config/design_config.yaml",
                        help="设计配置文件（默认 config/design_config.yaml）")
    parser.add_argument("--quiet", action="store_true", help="只输出最终状态")
    args = parser.parse_args(argv)
    return run(pathlib.Path(args.config), verbose=not args.quiet)


if __name__ == "__main__":
    raise SystemExit(main())
