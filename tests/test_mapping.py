"""Task B 的验收测试：结构模型与映射。

两部分：
1. 占位模式（``placeholder=True``）——只验证**结构表示**是否正确；
2. 真实序列模式——用一条**人工合成的小随机序列**（不是任何真实基因）验证
   "13 nt 与 scaffold 互补 / 26 nt 与 target 互补 / 拼接顺序为 5′-[26]-[13]-3′"。

运行：``python -m unittest discover -s tests -v``
"""

from __future__ import annotations

import json
import pathlib
import random
import sys
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from ssdna.mapping import (  # noqa: E402
    PlanConfig,
    plan_assembly,
    to_structure_json,
    toy_plan,
    validate_plan,
)
from ssdna.sequence import get_circular_subsequence, reverse_complement  # noqa: E402

TMP = ROOT / "build" / "tmp"
TMP.mkdir(parents=True, exist_ok=True)


def synthetic(n: int, seed: int) -> str:
    """确定性伪随机序列；只用于测试映射，不代表任何真实序列。"""
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


class TestToyStructure(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.plan = toy_plan()

    def test_staple_counts(self) -> None:
        self.assertEqual(len(self.plan.staples), 15)
        self.assertEqual(len(self.plan.staples_of("A")), 8)
        self.assertEqual(len(self.plan.staples_of("B")), 7)

    def test_domain_lengths(self) -> None:
        for staple in self.plan.staples:
            self.assertEqual(staple.scaffold_domain_length, 13)
            self.assertEqual(staple.target_domain_length, 26)
            self.assertEqual(staple.total_length, 39)

    def test_scaffold_tiling_is_exact(self) -> None:
        """13 nt 一块、无缝、无重叠地铺满覆盖区。"""
        covered: set[int] = set()
        for staple in self.plan.staples:
            site = staple.scaffold_binding
            for offset in range(site.length):
                position = (site.start + offset) % self.plan.config.scaffold_length
                self.assertNotIn(position, covered, "scaffold 上出现重叠")
                covered.add(position)
        self.assertEqual(len(covered), self.plan.covered_length)
        self.assertEqual(self.plan.covered_length, 13 * 15)

    def test_groups_alternate(self) -> None:
        groups = [s.group for s in self.plan.staples]
        self.assertEqual(groups[:4], ["A", "B", "A", "B"])

    def test_target_regions(self) -> None:
        regions = {r.group: r for r in self.plan.regions}
        self.assertEqual((regions["A"].start, regions["A"].end), (0, 208))
        self.assertEqual((regions["B"].start, regions["B"].end), (208, 390))
        self.assertEqual(self.plan.target_length, 390)

    def test_staples_point_at_their_own_region(self) -> None:
        for staple in self.plan.staples_of("A"):
            self.assertLessEqual(staple.target_end, 208)
        for staple in self.plan.staples_of("B"):
            self.assertGreaterEqual(staple.target_start, 208)
            self.assertLessEqual(staple.target_end, 390)

    def test_oligos_tile_the_circle(self) -> None:
        ordered = sorted(self.plan.oligos, key=lambda o: o.start)
        self.assertEqual(sum(o.length for o in ordered), 390)
        for a, b in zip(ordered, ordered[1:]):
            self.assertEqual(a.end, b.start)
        # 最后一条跨原点
        self.assertEqual((ordered[-1].start, ordered[-1].end), (350, 420))

    def test_two_boundary_oligos(self) -> None:
        boundary = self.plan.boundary_oligos()
        self.assertEqual(len(boundary), 2)
        self.assertEqual(len(self.plan.junctions), 2)

    def test_junction_is_covered_not_nicked(self) -> None:
        for junction in self.plan.junctions:
            self.assertTrue(junction.covered)
            self.assertEqual(junction.bond_state, "within_oligo")
        # 交界处不应该再出现一条 nick
        nick_positions = {n.position for n in self.plan.nicks}
        for junction in self.plan.junctions:
            self.assertNotIn(junction.position, nick_positions)

    def test_nicks_are_explicit_objects(self) -> None:
        self.assertEqual(len(self.plan.nicks), 7)
        for nick in self.plan.nicks:
            self.assertEqual(nick.upstream_end, "3'")
            self.assertEqual(nick.downstream_end, "5'")
            self.assertFalse(nick.bonded)  # 默认未连接

    def test_circle_is_not_closed_by_coverage(self) -> None:
        self.assertTrue(self.plan.is_circular)
        self.assertFalse(self.plan.circle_closed)

    def test_no_errors(self) -> None:
        errors = [c.message for c in self.plan.errors]
        self.assertEqual(errors, [], errors)

    def test_notes_are_informational(self) -> None:
        codes = {c.code for c in self.plan.checks if c.level == "note"}
        self.assertIn("placeholder_sequences", codes)
        self.assertIn("not_covalently_closed", codes)


class TestStructureJson(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = to_structure_json(toy_plan())

    def test_required_sections(self) -> None:
        for key in (
            "target",
            "scaffold",
            "staples",
            "oligos",
            "junctions",
            "nicks",
            "checks",
            "boundary_spanning_oligos",
        ):
            self.assertIn(key, self.payload)

    def test_staple_pairing_fields(self) -> None:
        first = self.payload["staples"][0]
        for key in (
            "staple_id",
            "group",
            "scaffold_start",
            "scaffold_length",
            "target_start",
            "target_end",
            "total_length",
            "domain_order",
        ):
            self.assertIn(key, first)
        self.assertEqual(first["total_length"], 39)
        self.assertEqual(first["scaffold_strand"], "-")
        self.assertEqual(first["target_strand"], "-")

    def test_boundary_list_has_two(self) -> None:
        self.assertEqual(len(self.payload["boundary_spanning_oligos"]), 2)

    def test_nick_and_junction_records(self) -> None:
        nick = self.payload["nicks"][0]
        self.assertEqual(
            set(nick),
            {
                "nick_id",
                "upstream_oligo",
                "upstream_end",
                "downstream_oligo",
                "downstream_end",
                "genomic_or_target_position",
                "five_prime_phosphate",
                "three_prime_hydroxyl",
                "intended_ligation",
                "bonded",
            },
        )
        junction = self.payload["junctions"][0]
        self.assertEqual(junction["bond_state"], "within_oligo")
        self.assertIn(junction["covered_by"], self.payload["boundary_spanning_oligos"])

    def test_json_serialisable(self) -> None:
        text = json.dumps(self.payload, ensure_ascii=False)
        self.assertGreater(len(text), 1000)
        path = TMP / "structure_probe.json"
        path.write_text(text, encoding="utf-8")
        self.assertTrue(json.loads(path.read_text(encoding="utf-8")))


class TestCircularScaffoldMapping(unittest.TestCase):
    def test_block_wraps_across_origin(self) -> None:
        cfg = PlanConfig(
            scaffold_length=1200,
            covered_start=1190,
            covered_length=195,
            n_staples=15,
            target_A_length=208,
            target_B_length=182,
            target_topology="circular",
            oligo_intervals=[
                (30, 70), (70, 140), (140, 190), (190, 226),
                (226, 280), (280, 350), (350, 420),
            ],
        )
        plan = plan_assembly(cfg)
        # 最后一块从 1190 开始、到 1203 结束 —— 跨过原点
        last = plan.staples[-1].scaffold_binding
        self.assertEqual((last.start, last.end), (1190, 1203))
        self.assertTrue(last.wraps)
        self.assertEqual(last.length, 13)

    def test_all_blocks_within_single_turn(self) -> None:
        plan = toy_plan()
        for staple in plan.staples:
            site = staple.scaffold_binding
            self.assertGreaterEqual(site.start, 0)
            self.assertLess(site.start, plan.config.scaffold_length)


class TestSequenceMapping(unittest.TestCase):
    """真实序列模式：验证互补关系与 5′-[26]-[13]-3′ 的拼接顺序。"""

    @classmethod
    def setUpClass(cls) -> None:
        cls.scaffold_seq = synthetic(1200, seed=7)
        cls.target_seq = synthetic(390, seed=11)
        cfg = PlanConfig(
            scaffold_length=1200,
            covered_start=0,
            covered_length=195,
            n_staples=15,
            target_A_length=208,
            target_B_length=182,
            target_topology="circular",
            oligo_intervals=[
                (30, 70), (70, 140), (140, 190), (190, 226),
                (226, 280), (280, 350), (350, 420),
            ],
            placeholder=False,
            scaffold_sequence=cls.scaffold_seq,
            target_sequence=cls.target_seq,
        )
        cls.plan = plan_assembly(cfg)

    def test_staple_length(self) -> None:
        for staple in self.plan.staples:
            self.assertEqual(len(staple.sequence), 39)
            self.assertEqual(set(staple.sequence), set("ACGT"))

    def test_26nt_arm_is_revcomp_of_target(self) -> None:
        for staple in self.plan.staples:
            arm = staple.sequence[:26]
            site = staple.target_binding
            expected = reverse_complement(
                get_circular_subsequence(self.target_seq, site.start, site.length)
            )
            self.assertEqual(arm, expected)

    def test_13nt_arm_is_revcomp_of_scaffold(self) -> None:
        for staple in self.plan.staples:
            arm = staple.sequence[26:]
            site = staple.scaffold_binding
            expected = reverse_complement(
                get_circular_subsequence(self.scaffold_seq, site.start, site.length)
            )
            self.assertEqual(arm, expected)

    def test_staple_arms_are_antiparallel(self) -> None:
        """13 nt 与 26 nt 分别与 scaffold / 目标反向互补。"""
        staple = self.plan.staples[0]
        self.assertNotEqual(staple.scaffold_binding.strand, "+")
        self.assertNotEqual(staple.target_binding.strand, "+")

    def test_group_pattern_swaps_still_consistent(self) -> None:
        cfg = PlanConfig(
            scaffold_length=1200,
            covered_start=0,
            covered_length=195,
            n_staples=15,
            target_A_length=208,
            target_B_length=182,
            target_topology="circular",
            oligo_intervals=[
                (30, 70), (70, 140), (140, 190), (190, 226),
                (226, 280), (280, 350), (350, 420),
            ],
            group_pattern="ABBB",
        )
        plan = plan_assembly(cfg)
        self.assertEqual([s.group for s in plan.staples[:4]], ["A", "B", "B", "B"])
        # ABBB 模式下 n_A=4，与 208 nt（需要 8 条）不再匹配，应报错
        self.assertTrue(any(c.code == "group_length" for c in plan.errors))


class TestLinearTarget(unittest.TestCase):
    def test_single_junction_and_no_wrap_nick(self) -> None:
        cfg = PlanConfig(
            scaffold_length=1200,
            covered_start=0,
            covered_length=195,
            n_staples=15,
            target_A_length=208,
            target_B_length=182,
            target_topology="linear",
            oligo_intervals=[
                (0, 60), (60, 140), (140, 190), (190, 226),
                (226, 280), (280, 350), (350, 390),
            ],
        )
        plan = plan_assembly(cfg)
        self.assertFalse(plan.is_circular)
        self.assertEqual(len(plan.junctions), 1)
        self.assertEqual(len(plan.nicks), len(plan.oligos) - 1)
        self.assertEqual([j.position for j in plan.junctions], [208])
        self.assertEqual(plan.junctions[0].covered_by, "T04")
        self.assertFalse(plan.circle_closed)


class TestValidatorCatchesErrors(unittest.TestCase):
    def _plan(self, **overrides):
        base = dict(
            scaffold_length=1200,
            covered_start=0,
            covered_length=195,
            n_staples=15,
            target_A_length=208,
            target_B_length=182,
            target_topology="circular",
            oligo_intervals=[
                (30, 70), (70, 140), (140, 190), (190, 226),
                (226, 280), (280, 350), (350, 420),
            ],
        )
        base.update(overrides)
        return plan_assembly(PlanConfig(**base))

    def codes(self, plan) -> set[str]:
        return {c.code for c in plan.errors}

    def test_group_length_mismatch(self) -> None:
        plan = self._plan(target_A_length=260, target_B_length=130)
        self.assertIn("group_length", self.codes(plan))

    def test_covered_length_mismatch(self) -> None:
        plan = self._plan(covered_length=200)
        self.assertIn("covered_length", self.codes(plan))

    def test_nick_lands_on_staple_seam(self) -> None:
        # 182 是一条 staple 臂的接缝
        plan = self._plan(
            oligo_intervals=[
                (30, 70), (70, 140), (140, 182), (182, 226), (226, 280),
                (280, 350), (350, 420),
            ]
        )
        self.assertIn("nick_on_seam", self.codes(plan))

    def test_junction_not_covered(self) -> None:
        # 去掉跨 208 的那条 oligo 的跨越性：把切点改到 208 本身
        plan = self._plan(
            oligo_intervals=[
                (30, 70), (70, 140), (140, 208), (208, 280),
                (280, 350), (350, 420),
            ]
        )
        self.assertIn("junction_coverage", self.codes(plan))

    def test_oligo_gap(self) -> None:
        plan = self._plan(
            oligo_intervals=[
                (30, 70), (70, 140), (140, 190), (190, 226),
                (226, 280), (280, 350), (352, 420),
            ]
        )
        self.assertIn("oligo_gap", self.codes(plan))

    def test_oligo_spans_no_seam(self) -> None:
        # 156 是接缝，切一条 150..170 的短 oligo 会跨过它；改成 160..175 不跨接缝
        plan = self._plan(
            oligo_intervals=[
                (30, 70), (70, 140), (140, 160), (160, 175),
                (175, 226), (226, 280), (280, 350), (350, 420),
            ]
        )
        self.assertIn("oligo_spans_no_seam", self.codes(plan))

    def test_ligation_makes_circle_closed(self) -> None:
        plan = self._plan(intended_ligation=True)
        self.assertTrue(plan.circle_closed)
        self.assertNotIn(
            "not_covalently_closed", {c.code for c in plan.checks if c.level == "note"}
        )


class TestSpecFieldNames(unittest.TestCase):
    """任务书点名的字段必须原样存在（不是同义改写）。"""

    def test_target_region_fields(self) -> None:
        region = toy_plan().regions[0]
        for name in ("region_id", "sequence", "start", "end", "strand", "group"):
            self.assertTrue(hasattr(region, name), name)

    def test_target_oligo_fields(self) -> None:
        oligo = toy_plan().oligos[0]
        for name in (
            "oligo_id", "sequence", "group", "oligo_type", "target_start",
            "target_end", "target_segments", "strand",
            "five_prime_modification", "three_prime_modification",
        ):
            self.assertTrue(hasattr(oligo, name), name)

    def test_staple_fields(self) -> None:
        staple = toy_plan().staples[0]
        for name in (
            "staple_id", "sequence", "scaffold_binding_region",
            "target_binding_region", "scaffold_start", "target_start",
            "target_end", "total_length",
        ):
            self.assertTrue(hasattr(staple, name), name)

    def test_nick_fields(self) -> None:
        nick = toy_plan().nicks[0]
        for name in (
            "nick_id", "upstream_oligo", "downstream_oligo", "upstream_end",
            "downstream_end", "genomic_or_target_position",
            "five_prime_phosphate", "three_prime_hydroxyl", "intended_ligation",
        ):
            self.assertTrue(hasattr(nick, name), name)

    def test_boundary_oligo_spans_two_regions_with_contiguous_sequence(self) -> None:
        plan = toy_plan()
        boundary = plan.boundary_oligos()
        self.assertTrue(boundary)
        for oligo in boundary:
            self.assertEqual(len(oligo.target_segments), 2, oligo.oligo_id)
            # 自身序列连续：长度就等于区间长度，中间没有被拆开
            self.assertEqual(len(oligo.sequence), oligo.length)


class TestNoImplicitAssumptions(unittest.TestCase):
    def test_no_overlap_assumed_between_oligos(self) -> None:
        """相邻 oligo 之间只有坐标相接，没有 overlap 字段，也不会自动合并。"""
        plan = toy_plan()
        ordered = sorted(plan.oligos, key=lambda o: o.start)
        for a, b in zip(ordered, ordered[1:]):
            self.assertEqual(a.end, b.start)
        ids = {o.oligo_id for o in plan.oligos}
        self.assertEqual(len(ids), len(plan.oligos))

    def test_staple_type_is_uniform(self) -> None:
        """不允许出现 crossover staple：所有 staple 都是 13+26。"""
        plan = toy_plan()
        for staple in plan.staples:
            self.assertEqual(staple.domain_order, "5p_target_26")
            self.assertEqual(staple.total_length, 39)

    def test_validate_is_pure(self) -> None:
        plan = toy_plan()
        before = len(plan.checks)
        validate_plan(plan)
        self.assertEqual(len(plan.checks), before)


if __name__ == "__main__":
    unittest.main()
