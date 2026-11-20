"""针对外部评审提出的 6 个缺陷的回归测试。

每个用例都对应一条具体指控，修复前会失败、修复后必须通过。
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
from ssdna.sequence import reverse_complement  # noqa: E402

DNA_FILE = ROOT / "data" / "scaffold.dna"
TARGET = ("ACGT" * 182)[:728]


def _oligos_and_staples(domain_order: str = "5p_target_26"):
    regions, oligos, _ = design.design_target_oligos(
        TARGET,
        [
            {"region_id": "target_A", "group": "A", "start": 0, "end": 364},
            {"region_id": "target_B", "group": "B", "start": 364, "end": 728},
        ],
        {"cuts": [108, 212, 316, 420, 524, 628], "topology": "circular"},
        [{"junction_id": "J1", "position": 364}, {"junction_id": "J2", "position": 0}],
    )
    scaffold = load_snapgene(DNA_FILE)
    staples, _ = design.design_staples(
        scaffold,
        oligos,
        {
            "scaffold_binding_length": 13,
            "target_binding_length": 26,
            "domain_order": domain_order,
            "group_spans": {"A": [0, 364], "B": [364, 728]},
        },
        target_sequence=TARGET,
    )
    return scaffold, regions, oligos, staples


class TestLigationIsNotClosure(unittest.TestCase):
    """① 不能把"计划连接"当成"已经连接"。"""

    def _plan(self, **kwargs):
        base = dict(
            scaffold_length=1200, covered_start=0, covered_length=195, n_staples=15,
            target_A_length=208, target_B_length=182, target_topology="circular",
            oligo_intervals=[(30, 70), (70, 140), (140, 190), (190, 226),
                             (226, 280), (280, 350), (350, 420)],
        )
        base.update(kwargs)
        return plan_assembly(PlanConfig(**base))

    def test_intent_without_phosphate_is_not_closed(self) -> None:
        plan = self._plan(intended_ligation=True, five_prime_phosphate="no")
        self.assertFalse(plan.circle_closed)
        self.assertFalse(plan.nicks[0].bonded)
        self.assertEqual(plan.nicks[0].ligation_state, "planned")

    def test_explicit_ligated_is_closed(self) -> None:
        plan = self._plan()
        for nick in plan.nicks:
            nick.ligation_state = "ligated"
        self.assertTrue(plan.circle_closed)

    def test_nick_state_comes_from_override(self) -> None:
        scaffold, _, oligos, _ = _oligos_and_staples()
        nicks = design.build_nicks(
            oligos, topology="circular", five_prime_phosphate="no",
            intended_ligation=True,
            nick_overrides={"N2": {"ligation_state": "ligated",
                                   "five_prime_phosphate": "yes"}},
        )
        self.assertEqual(nicks[0].ligation_state, "planned")
        self.assertEqual(nicks[1].ligation_state, "ligated")
        self.assertEqual(nicks[1].five_prime_phosphate, "yes")
        self.assertFalse(nicks[0].bonded)


class TestRegionMustBeExplicit(unittest.TestCase):
    """② scaffold 有效区域不能默认取整段。"""

    def test_region_is_a_blocking_parameter(self) -> None:
        cfg = design.load_config(ROOT / "config" / "design_config.yaml")
        self.assertIn("scaffold.region", design.find_missing(cfg))

    def test_whole_and_subregion(self) -> None:
        scaffold = load_snapgene(DNA_FILE)
        whole, mode, bounds = cli._resolve_region(
            scaffold, {"scaffold": {"region": "whole"}}, "+")
        self.assertEqual(mode, "whole")
        self.assertIsNone(bounds)
        self.assertEqual(len(whole), 1200)

        piece, mode, bounds = cli._resolve_region(
            scaffold, {"scaffold": {"region": [100, 400]}}, "+")
        self.assertEqual(mode, "subregion")
        self.assertEqual(bounds, (100, 400))
        self.assertEqual(piece, scaffold.sequence[100:400])

    def test_reversed_strand_flips_the_subregion(self) -> None:
        scaffold = load_snapgene(DNA_FILE)
        piece, _, _ = cli._resolve_region(
            scaffold, {"scaffold": {"region": [100, 400]}}, "-")
        self.assertEqual(piece, reverse_complement(scaffold.sequence[100:400]))

    def test_bad_region_is_rejected(self) -> None:
        scaffold = load_snapgene(DNA_FILE)
        with self.assertRaises(ValueError):
            cli._resolve_region(scaffold, {"scaffold": {"region": "most"}}, "+")


class TestCutsAreStrict(unittest.TestCase):
    """③ 切点不允许静默取模、去重或越界。"""

    def _run(self, cuts):
        return design.design_target_oligos(
            TARGET,
            [{"region_id": "target_A", "group": "A", "start": 0, "end": 364},
             {"region_id": "target_B", "group": "B", "start": 364, "end": 728}],
            {"cuts": cuts, "topology": "circular"},
            [{"junction_id": "J1", "position": 364}],
        )

    def test_negative_cut_rejected(self) -> None:
        with self.assertRaises(ValueError) as ctx:
            self._run([-628, 100])
        self.assertIn("越界", str(ctx.exception))

    def test_out_of_range_cut_rejected(self) -> None:
        with self.assertRaises(ValueError):
            self._run([828])

    def test_duplicate_cut_rejected(self) -> None:
        with self.assertRaises(ValueError) as ctx:
            self._run([100, 100, 300])
        self.assertIn("重复", str(ctx.exception))

    def test_non_integer_cut_rejected(self) -> None:
        with self.assertRaises(ValueError):
            self._run([100.5])  # type: ignore[list-item]


class TestDomainOrderConsistency(unittest.TestCase):
    """④ 切到 5p_scaffold_13 后，验证器必须按实际顺序拆 domain。"""

    def test_reversed_order_passes_validation(self) -> None:
        _, _, _, staples = _oligos_and_staples("5p_scaffold_13")
        rows = [
            {"staple_id": s.staple_id, "sequence": s.sequence,
             "target_start": s.target_binding_region.start}
            for s in staples
        ]
        issues = validate.validate_target_binding(
            TARGET, rows, target_arm=26, domain_order="5p_scaffold_13")
        self.assertEqual(issues, [], [i.message for i in issues][:3])

    def test_wrong_order_still_detected(self) -> None:
        """反向构造的设计用默认顺序去验证，必须报错——防止验证器"总是通过"。"""
        _, _, _, staples = _oligos_and_staples("5p_scaffold_13")
        rows = [
            {"staple_id": s.staple_id, "sequence": s.sequence,
             "target_start": s.target_binding_region.start}
            for s in staples
        ]
        issues = validate.validate_target_binding(
            TARGET, rows, target_arm=26, domain_order="5p_target_26")
        self.assertTrue(issues)

    def test_target_domain_slicing(self) -> None:
        _, _, _, staples = _oligos_and_staples("5p_scaffold_13")
        staple = staples[0]
        start = staple.target_binding_region.start
        self.assertEqual(staple.sequence[13:39], reverse_complement(TARGET[start:start + 26]))


class TestCircularOverlap(unittest.TestCase):
    """⑤ 环状 scaffold 上跨原点的重叠必须检出。"""

    def test_wrap_around_overlap(self) -> None:
        staples = [
            {"staple_id": "S1", "sequence": "A" * 39, "scaffold_start": 95,
             "target_start": 0, "target_length": 50},
            {"staple_id": "S2", "sequence": "A" * 39, "scaffold_start": 0,
             "target_start": 26, "target_length": 50},
        ]
        issues = validate.validate_scaffold_mapping(
            "A" * 100, staples, scaffold_arm=13, target_arm=26)
        self.assertIn("scaffold_site_overlap", {i.code for i in issues})

    def test_near_miss_is_not_reported(self) -> None:
        staples = [
            {"staple_id": "S1", "sequence": "A" * 39, "scaffold_start": 95,
             "target_start": 0, "target_length": 50},
            {"staple_id": "S2", "sequence": "A" * 39, "scaffold_start": 8,
             "target_start": 26, "target_length": 50},
        ]
        issues = validate.validate_scaffold_mapping(
            "A" * 100, staples, scaffold_arm=13, target_arm=26)
        self.assertNotIn("scaffold_site_overlap", {i.code for i in issues})


class TestPerNickChemistry(unittest.TestCase):
    """⑥ 每条 nick 的端基状态要能独立配置。"""

    def test_defaults_are_unknown(self) -> None:
        _, _, oligos, _ = _oligos_and_staples()
        nicks = design.build_nicks(
            oligos, topology="circular", five_prime_phosphate="unknown",
            intended_ligation=False)
        self.assertEqual(nicks[0].three_prime_hydroxyl, "unknown")

    def test_mixed_states(self) -> None:
        _, _, oligos, _ = _oligos_and_staples()
        nicks = design.build_nicks(
            oligos, topology="circular", five_prime_phosphate="no",
            intended_ligation=False,
            nick_overrides={"N3": {"five_prime_phosphate": "yes",
                                   "three_prime_hydroxyl": "yes"}})
        states = {(n.nick_id, n.five_prime_phosphate) for n in nicks}
        self.assertIn(("N1", "no"), states)
        self.assertIn(("N3", "yes"), states)


class TestSeamRuleRunsInMainPath(unittest.TestCase):
    """架构问题：主流程此前没跑"nick 避开接缝"这条规则。"""

    def test_nick_on_seam_is_reported_by_check_design(self) -> None:
        scaffold, regions, oligos, staples = _oligos_and_staples()
        seams = sorted({s.target_binding_region.start % 728 for s in staples})
        bad_nick = design.build_nicks(
            oligos, topology="circular", five_prime_phosphate="no",
            intended_ligation=False)[0]
        object.__setattr__(bad_nick, "genomic_or_target_position", seams[1])
        checks = design.check_design(
            scaffold, scaffold.sequence, TARGET, regions, oligos, staples,
            junctions=[{"junction_id": "J1", "position": 364},
                       {"junction_id": "J2", "position": 0}],
            arm_s=13, arm_t=26, nicks=[bad_nick],
        )
        self.assertIn("nick_on_seam", {c.code for c in checks})


if __name__ == "__main__":
    unittest.main()
