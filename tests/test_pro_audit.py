"""pro 模型审核后新增的回归测试（2026-11-21）。

每条对应一处本次审核发现并修复的问题；修复前失败、修复后通过。
"""

from __future__ import annotations

import pathlib
import sys
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
_VENDORED = ROOT / "build" / "pylibs"
if _VENDORED.exists():
    sys.path.insert(0, str(_VENDORED))

import main as cli  # noqa: E402
from ssdna import design, validate  # noqa: E402
from ssdna.io import load_snapgene  # noqa: E402
from ssdna.mapping import PlanConfig, plan_assembly  # noqa: E402
from ssdna.sequence import (  # noqa: E402
    normalize_sequence,
    reverse_complement,
)

TMP = ROOT / "build" / "tmp" / "pro_audit"


def _demo_inputs():
    scaffold = load_snapgene(ROOT / "data" / "scaffold.dna")
    cfg = design.demo_config(design.DEFAULT_CONFIG)
    target = cfg["target"]["source"]["inline"]
    working, _, _ = cli._resolve_region(scaffold, cfg, "+")
    return design.DesignInputs(
        scaffold=scaffold, scaffold_strand="+", working_sequence=working,
        target_sequence=target, target_length=len(target), config=cfg,
    )


class TestRegionSubregion(unittest.TestCase):
    """C1：region 子区段必须真正进入 staple 生成。"""

    def test_subregion_design_is_consistent(self) -> None:
        scaffold = load_snapgene(ROOT / "data" / "scaffold.dna")
        cfg = design.demo_config(design.DEFAULT_CONFIG)
        cfg["scaffold"]["region"] = [100, 500]
        target = cfg["target"]["source"]["inline"]
        working, mode, bounds = cli._resolve_region(scaffold, cfg, "+")
        self.assertEqual(mode, "subregion")
        self.assertEqual(bounds, (100, 500))
        self.assertEqual(len(working), 400)

        inputs = design.DesignInputs(
            scaffold=scaffold, scaffold_strand="+", working_sequence=working,
            target_sequence=target, target_length=len(target), config=cfg,
            scaffold_region="subregion", scaffold_region_bounds=(100, 500),
        )
        payload = design.run_design(inputs, outdir=TMP / "subregion")
        self.assertEqual(payload["status"], "OK", [c["code"] for c in payload["checks"]])
        # 所有 staple 的 scaffold 坐标都应落在工作序列 [0, 400) 内
        self.assertEqual(len(payload["staples"]), 28)
        for s in payload["staples"]:
            self.assertLess(s["scaffold_start"] + s["scaffold_length"], 400)

    def test_empty_region_is_rejected(self) -> None:
        scaffold = load_snapgene(ROOT / "data" / "scaffold.dna")
        with self.assertRaises(ValueError):
            cli._resolve_region(scaffold, {"scaffold": {"region": [5, 5]}}, "+")

    def test_linear_scaffold_does_not_wrap(self) -> None:
        scaffold = load_snapgene(ROOT / "data" / "scaffold.dna")
        scaffold.topology = "linear"
        with self.assertRaises(IndexError):
            cli._resolve_region(scaffold, {"scaffold": {"region": [1190, 10]}}, "+")


class TestPatternCountValidation(unittest.TestCase):
    """agent1-C1：分组模式与目标跨度不一致时必须报错。"""

    def test_uneven_groups_are_rejected(self) -> None:
        scaffold = load_snapgene(ROOT / "data" / "scaffold.dna")
        cfg = design.demo_config(design.DEFAULT_CONFIG)
        cfg["target"]["group_boundaries"] = [26]  # A 26 nt，B 702 nt，AB 交替无法满足
        target = cfg["target"]["source"]["inline"]
        inputs = design.DesignInputs(
            scaffold=scaffold, scaffold_strand="+", working_sequence=scaffold.sequence,
            target_sequence=target, target_length=len(target), config=cfg,
        )
        payload = design.run_design(inputs, outdir=TMP / "uneven")
        self.assertEqual(payload["status"], "BLOCKED")
        self.assertTrue(any(c["code"] == "group_pattern_mismatch" for c in payload["checks"]))


class TestSequenceHardening(unittest.TestCase):
    """sequence.py：数字不再静默删除；表外字符不再静默透传。"""

    def test_digits_are_rejected(self) -> None:
        with self.assertRaises(ValueError):
            normalize_sequence("ACGT12ACGT")

    def test_reverse_complement_rejects_illegal_chars(self) -> None:
        for seq in ("AT-GC", "AT?GC", "AT1GC"):
            with self.assertRaises(ValueError, msg=seq):
                reverse_complement(seq)

    def test_reverse_complement_keeps_iupac(self) -> None:
        self.assertEqual(reverse_complement("ACGTRYN"), "NRYACGT")


class TestFindInvalid(unittest.TestCase):
    def test_bad_types_are_caught(self) -> None:
        cfg = design.demo_config(design.DEFAULT_CONFIG)
        cfg["staple"]["scaffold_binding_length"] = "abc"
        cfg["staple"]["domain_order"] = "wrong"
        invalid = design.find_invalid(cfg)
        self.assertIn("staple.scaffold_binding_length", invalid)
        self.assertIn("staple.domain_order", invalid)

    def test_clean_config_has_no_invalid(self) -> None:
        cfg = design.demo_config(design.DEFAULT_CONFIG)
        self.assertEqual(design.find_invalid(cfg), {})


class TestMappingAlignment(unittest.TestCase):
    """mapping 与 design 的 domain_order / junction / 尾巴策略对齐。"""

    def test_domain_order_is_honored(self) -> None:
        cfg = PlanConfig(
            scaffold_length=1200, covered_start=0, covered_length=195, n_staples=15,
            target_A_length=208, target_B_length=182, target_topology="circular",
            oligo_intervals=[(30, 70), (70, 140), (140, 190), (190, 226),
                             (226, 280), (280, 350), (350, 420)],
            staple_domain_order="5p_scaffold_13",
        )
        plan = plan_assembly(cfg)
        self.assertEqual(plan.staples[0].domain_order, "5p_scaffold_13")

    def test_junction_groups_differ(self) -> None:
        plan = plan_assembly(PlanConfig(
            scaffold_length=1200, covered_start=0, covered_length=195, n_staples=15,
            target_A_length=208, target_B_length=182, target_topology="circular",
            oligo_intervals=[(30, 70), (70, 140), (140, 190), (190, 226),
                             (226, 280), (280, 350), (350, 420)],
        ))
        groups = {(j.position, j.upstream_group, j.downstream_group) for j in plan.junctions}
        self.assertIn((208, "A", "B"), groups)
        self.assertIn((0, "B", "A"), groups)

    def test_leftover_is_note_not_error(self) -> None:
        # 组长度 210（26*8=208，余 2）应是 note，不是 error
        plan = plan_assembly(PlanConfig(
            scaffold_length=1200, covered_start=0, covered_length=195, n_staples=15,
            target_A_length=210, target_B_length=182, target_topology="circular",
            oligo_intervals=[(30, 70), (70, 140), (140, 190), (190, 226),
                             (226, 280), (280, 350), (350, 420)],
        ))
        self.assertNotIn("group_length", {c.code for c in plan.errors})


class TestValidateGuards(unittest.TestCase):
    def test_empty_staples_is_not_pass(self) -> None:
        from ssdna.validate import Bundle

        outdir = TMP / "empty_staples"
        payload = design.run_design(_demo_inputs(), outdir=outdir)
        # 清空 staples.csv
        path = outdir / "staples.csv"
        header = path.read_text(encoding="utf-8").splitlines()[0]
        path.write_text(header + "\n", encoding="utf-8")
        bundle = Bundle(
            directory=outdir,
            target_sequence=payload["_sequences"]["target"],
            scaffold_sequence=payload["_sequences"]["scaffold_working"],
            oligos=validate._read_csv(outdir / "target_oligos.csv"),
            staples=validate._read_csv(path),
            nicks=validate._read_csv(outdir / "nicks.csv"),
            report=payload,
        )
        result = validate.run_validation(bundle)
        self.assertEqual(result["status"], "FAIL")
        self.assertTrue(any(i["code"] == "no_staples" for i in result["issues"]))

    def test_target_wrap_overlap_detected(self) -> None:
        # 目标 728 nt，两条臂 [720,746) 与 [0,26) 跨原点重叠 18 nt
        seq = "ACGT" * 182
        staples = [
            {"staple_id": "S1", "sequence": "A" * 39, "target_start": 720},
            {"staple_id": "S2", "sequence": "A" * 39, "target_start": 0},
        ]
        issues = validate.validate_target_binding(seq, staples, target_arm=26)
        self.assertIn("target_site_overlap", {i.code for i in issues})


if __name__ == "__main__":
    unittest.main()
