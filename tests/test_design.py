"""Task C 的验收测试：target oligo 与 staple 序列生成。"""

from __future__ import annotations

import csv
import json
import pathlib
import sys
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
_VENDORED = ROOT / "build" / "pylibs"
if _VENDORED.exists():
    sys.path.insert(0, str(_VENDORED))

from ssdna import design  # noqa: E402
from ssdna.io import load_snapgene  # noqa: E402
from ssdna.sequence import get_circular_subsequence, reverse_complement  # noqa: E402

TMP = ROOT / "build" / "tmp" / "design"
DNA_FILE = ROOT / "minimalplasmid(1).dna"


def demo_inputs():
    cfg = design.demo_config(design.DEFAULT_CONFIG)
    scaffold = load_snapgene(DNA_FILE)
    target = cfg["target"]["source"]["inline"]
    return cfg, design.DesignInputs(
        scaffold=scaffold,
        scaffold_strand="+",
        working_sequence=scaffold.sequence,
        target_sequence=target,
        target_length=len(target),
        config=cfg,
    )


class TestMissingParamsGate(unittest.TestCase):
    def test_default_config_is_blocked(self) -> None:
        missing = design.find_missing(json.loads(json.dumps(design.DEFAULT_CONFIG)))
        self.assertIn("scaffold.strand", missing)
        self.assertIn("target.source", missing)
        self.assertIn("target.group_boundaries", missing)
        self.assertIn("target.oligo_cuts", missing)

    def test_cli_refuses_and_writes_nothing(self) -> None:
        outdir = TMP / "blocked"
        if outdir.exists():
            for p in outdir.rglob("*"):
                p.unlink()
            outdir.rmdir()
        code = design.main(["--outdir", str(outdir)])
        self.assertEqual(code, 2)
        self.assertFalse(outdir.exists(), "缺参数时不应该创建输出")

    def test_demo_config_is_complete(self) -> None:
        cfg = design.demo_config(json.loads(json.dumps(design.DEFAULT_CONFIG)))
        self.assertEqual(design.find_missing(cfg), [])


class TestTargetOligos(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.cfg, cls.inputs = demo_inputs()
        cls.regions, cls.oligos, cls.checks = design.design_target_oligos(
            cls.inputs.target_sequence,
            [
                {"region_id": "target_A", "group": "A", "start": 0, "end": 364},
                {"region_id": "target_B", "group": "B", "start": 364, "end": 728},
            ],
            {
                "cuts": [100, 200, 300, 420, 540, 660],
                "topology": "circular",
                "five_prime_modification": "no",
                "three_prime_modification": "OH",
            },
            [
                {"junction_id": "J1", "position": 364},
                {"junction_id": "J2", "position": 0},
            ],
        )

    def test_count_and_coverage(self) -> None:
        self.assertEqual(len(self.oligos), 6)
        self.assertEqual(sum(o.length for o in self.oligos), 728)

    def test_intervals_are_contiguous(self) -> None:
        ordered = sorted(self.oligos, key=lambda o: o.target_start)
        for a, b in zip(ordered, ordered[1:]):
            self.assertEqual(b.target_start, a.target_end)
        self.assertEqual(ordered[-1].target_end, 100 + 728)

    def test_sequence_is_the_original_slice(self) -> None:
        target = self.inputs.target_sequence
        for oligo in self.oligos:
            expected = get_circular_subsequence(
                target, oligo.target_start, oligo.length
            )
            self.assertEqual(oligo.sequence, expected)

    def test_boundary_oligos(self) -> None:
        boundary = [o for o in self.oligos if o.oligo_type == "boundary_spanning"]
        self.assertEqual(len(boundary), 2)
        for oligo in boundary:
            self.assertEqual(len(oligo.target_segments), 2)
            self.assertEqual(len(oligo.sequence), oligo.length)

    def test_group_assignment(self) -> None:
        for oligo in self.oligos:
            self.assertEqual(oligo.group, "A" if oligo.target_start < 364 else "B")


class TestStaples(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cfg, inputs = demo_inputs()
        cls.inputs = inputs
        cls.regions, cls.oligos, _ = design.design_target_oligos(
            inputs.target_sequence,
            [
                {"region_id": "target_A", "group": "A", "start": 0, "end": 364},
                {"region_id": "target_B", "group": "B", "start": 364, "end": 728},
            ],
            {"cuts": [100, 200, 300, 420, 540, 660], "topology": "circular"},
            [{"junction_id": "J1", "position": 364}, {"junction_id": "J2", "position": 0}],
        )
        cls.staples, cls.checks = design.design_staples(
            inputs.scaffold,
            cls.oligos,
            {
                "scaffold_binding_length": 13,
                "target_binding_length": 26,
                "total_length": 39,
                "group_pattern": "AB",
                "covered_start": 0,
                "domain_order": "5p_target_26",
                "group_spans": {"A": [0, 364], "B": [364, 728]},
            },
            scaffold_strand="+",
            target_sequence=inputs.target_sequence,
        )

    def test_counts(self) -> None:
        self.assertEqual(len(self.staples), 28)
        self.assertEqual(sum(1 for s in self.staples if s.group == "A"), 14)
        self.assertEqual(sum(1 for s in self.staples if s.group == "B"), 14)

    def test_lengths(self) -> None:
        for staple in self.staples:
            self.assertEqual(staple.scaffold_domain_length, 13)
            self.assertEqual(staple.target_domain_length, 26)
            self.assertEqual(len(staple.sequence), 39)

    def test_domains_are_reverse_complements(self) -> None:
        """独立复算：26 nt == revcomp(目标片段)，13 nt == revcomp(scaffold 片段)。"""
        scaffold_seq = self.inputs.working_sequence
        target_seq = self.inputs.target_sequence
        for staple in self.staples:
            t_site = staple.target_binding_region
            s_site = staple.scaffold_binding_region
            t_piece = get_circular_subsequence(target_seq, t_site.start, 26)
            s_piece = get_circular_subsequence(scaffold_seq, s_site.start, 13)
            self.assertEqual(staple.sequence[:26], reverse_complement(t_piece))
            self.assertEqual(staple.sequence[26:], reverse_complement(s_piece))

    def test_scaffold_blocks_do_not_overlap(self) -> None:
        occupied: set[int] = set()
        for staple in self.staples:
            for offset in range(13):
                pos = (staple.scaffold_binding_region.start + offset) % 1200
                self.assertNotIn(pos, occupied)
                occupied.add(pos)
        self.assertEqual(len(occupied), 364)

    def test_coverage_walks_backwards_on_scaffold(self) -> None:
        """第 0 条从覆盖区右端开始，向左推进（因为 scaffold 在装配体中 3′→5′）。"""
        starts = [s.scaffold_binding_region.start for s in self.staples[:3]]
        self.assertEqual(starts, [351, 338, 325])

    def test_flip_domain_order_changes_only_order(self) -> None:
        flipped, _ = design.design_staples(
            self.inputs.scaffold,
            self.oligos,
            {
                "scaffold_binding_length": 13,
                "target_binding_length": 26,
                "domain_order": "5p_scaffold_13",
                "group_spans": {"A": [0, 364], "B": [364, 728]},
            },
            scaffold_strand="+",
            target_sequence=self.inputs.target_sequence,
        )
        for a, b in zip(self.staples, flipped):
            self.assertEqual(a.sequence, b.sequence[13:] + b.sequence[:13])


class TestChecks(unittest.TestCase):
    def test_demo_has_no_errors(self) -> None:
        cfg, inputs = demo_inputs()
        payload = design.run_design(inputs, outdir=TMP / "demo")
        errors = [c for c in payload["checks"] if c["level"] == "error"]
        self.assertEqual(errors, [], errors)
        self.assertEqual(payload["status"], "OK")

    def test_length_mismatch_is_caught(self) -> None:
        cfg, inputs = demo_inputs()
        regions, oligos, _ = design.design_target_oligos(
            inputs.target_sequence,
            [
                {"region_id": "target_A", "group": "A", "start": 0, "end": 364},
                {"region_id": "target_B", "group": "B", "start": 364, "end": 728},
            ],
            {"cuts": [100, 200, 300, 420, 540, 660], "topology": "circular"},
            [{"junction_id": "J1", "position": 364}],
        )
        staples, staple_checks = design.design_staples(
            inputs.scaffold, oligos,
            {
                "scaffold_binding_length": 13,
                "target_binding_length": 26,
                "total_length": 40,  # 故意与 13+26 不符
                "group_spans": {"A": [0, 364], "B": [364, 728]},
            },
            target_sequence=inputs.target_sequence,
        )
        self.assertIn("total_length", {c.code for c in staple_checks})
        checks = design.check_design(
            inputs.scaffold, inputs.working_sequence, inputs.target_sequence,
            regions, oligos, staples,
            junctions=[{"junction_id": "J1", "position": 364}],
            arm_s=13, arm_t=26,
        )
        # 长度不一致由 design_staples 报出（上面已断言），这里只确认 check_design 不误报
        self.assertNotIn("not_complementary", {c.code for c in checks})

    def test_tampered_sequence_is_caught(self) -> None:
        cfg, inputs = demo_inputs()
        regions, oligos, _ = design.design_target_oligos(
            inputs.target_sequence,
            [
                {"region_id": "target_A", "group": "A", "start": 0, "end": 364},
                {"region_id": "target_B", "group": "B", "start": 364, "end": 728},
            ],
            {"cuts": [100, 200, 300, 420, 540, 660], "topology": "circular"},
            [{"junction_id": "J1", "position": 364}],
        )
        staples, _ = design.design_staples(
            inputs.scaffold, oligos,
            {
                "scaffold_binding_length": 13,
                "target_binding_length": 26,
                "group_spans": {"A": [0, 364], "B": [364, 728]},
            },
            target_sequence=inputs.target_sequence,
        )
        # 篡改一条 staple 的最后一个碱基
        bad = staples[0]
        bad.sequence = bad.sequence[:-1] + ("A" if bad.sequence[-1] != "A" else "C")
        checks = design.check_design(
            inputs.scaffold, inputs.working_sequence, inputs.target_sequence,
            regions, oligos, staples,
            junctions=[{"junction_id": "J1", "position": 364}],
            arm_s=13, arm_t=26,
        )
        self.assertIn("not_complementary", {c.code for c in checks})


class TestOutputs(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cfg, inputs = demo_inputs()
        cls.payload = design.run_design(inputs, outdir=TMP / "out")
        cls.outdir = TMP / "out"

    def test_files(self) -> None:
        for name in ("target_oligos.csv", "staples.csv",
                     "design_config.yaml", "design_report.json"):
            path = self.outdir / name
            self.assertTrue(path.exists(), name)
            self.assertGreater(path.stat().st_size, 0)

    def test_oligo_csv(self) -> None:
        with (self.outdir / "target_oligos.csv").open(encoding="utf-8-sig") as fh:
            rows = list(csv.DictReader(fh))
        self.assertEqual(len(rows), 6)
        self.assertEqual(sum(r["oligo_type"] == "boundary_spanning" for r in rows), 2)
        for row in rows:
            self.assertEqual(len(row["sequence"]), int(row["length"]))

    def test_staple_csv(self) -> None:
        with (self.outdir / "staples.csv").open(encoding="utf-8-sig") as fh:
            rows = list(csv.DictReader(fh))
        self.assertEqual(len(rows), 28)
        for row in rows:
            self.assertEqual(len(row["sequence"]), 39)
            self.assertEqual(row["total_length"], "39")

    def test_report(self) -> None:
        payload = json.loads(
            (self.outdir / "design_report.json").read_text(encoding="utf-8"))
        self.assertEqual(payload["status"], "OK")
        self.assertEqual(payload["counts"]["staples"], 28)
        self.assertEqual(payload["counts"]["boundary_oligos"], 2)
        self.assertIn("target_sequence_sha256", payload["inputs"])

    def test_yaml_roundtrip(self) -> None:
        path = self.outdir / "design_config.yaml"
        text = path.read_text(encoding="utf-8")
        self.assertIn("scaffold_binding_length: 13", text.replace('"', "").replace("'", ""))
        try:
            import yaml
        except ImportError:
            self.skipTest("无 PyYAML")
        data = yaml.safe_load(text)
        self.assertEqual(data["staple"]["target_binding_length"], 26)


class TestPhaseScan(unittest.TestCase):
    def test_scan_returns_sorted_rows(self) -> None:
        from ssdna import thermo

        if not thermo.available():
            self.skipTest("未安装 Biopython")
        cfg, inputs = demo_inputs()
        rows = design.scan_scaffold_phase(
            inputs.scaffold,
            {
                "scaffold_binding_length": 13,
                "target_binding_length": 26,
                "group_spans": {"A": [0, 364], "B": [364, 728]},
            },
            step=13,
        )
        self.assertGreater(len(rows), 10)
        stdevs = [r["tm_stdev"] for r in rows]
        self.assertEqual(stdevs, sorted(stdevs))
        for row in rows:
            self.assertLessEqual(row["tm_min"], row["tm_max"])

    def test_scan_does_not_generate_sequences(self) -> None:
        """扫描只读 scaffold，不写任何序列文件。"""
        from ssdna import thermo

        if not thermo.available():
            self.skipTest("未安装 Biopython")
        cfg, inputs = demo_inputs()
        before = set(p.name for p in (TMP / "out").glob("*")) if (TMP / "out").exists() else set()
        design.scan_scaffold_phase(
            inputs.scaffold,
            {"scaffold_binding_length": 13, "target_binding_length": 26,
             "group_spans": {"A": [0, 364], "B": [364, 728]}},
            step=13,
        )
        after = set(p.name for p in (TMP / "out").glob("*")) if (TMP / "out").exists() else set()
        self.assertEqual(before, after)


class TestCutSuggestion(unittest.TestCase):
    """目标长度不必是 26 的倍数：接缝只决定"不能切在哪"。"""

    SEAMS = [26 * k for k in range(1, 29)]

    def _check(self, length: int, cuts: list[int]) -> None:
        lens = [cuts[0]] + [b - a for a, b in zip(cuts, cuts[1:])] + [length - cuts[-1]]
        self.assertTrue(all(40 <= L <= 120 for L in lens), lens)
        for cut in cuts:
            self.assertGreater(min(abs(cut - s) for s in self.SEAMS), 3)
        bounds = [0] + cuts + [length]
        for a, b in zip(bounds, bounds[1:]):
            self.assertTrue(any(a < s < b for s in self.SEAMS), (a, b))

    def test_exact_multiple_of_26(self) -> None:
        cuts = design.suggest_oligo_cuts(728, seams=self.SEAMS)
        self._check(728, cuts)

    def test_not_a_multiple_of_26(self) -> None:
        """720 / 723 这类长度同样能切，不需要凑成 26 的倍数。"""
        for length in (720, 723, 715):
            cuts = design.suggest_oligo_cuts(length, seams=self.SEAMS)
            self._check(length, cuts)

    def test_impossible_constraints_raise(self) -> None:
        with self.assertRaises(ValueError):
            design.suggest_oligo_cuts(60, seams=list(range(0, 60, 2)),
                                      min_length=40, max_length=50)


class TestThermo(unittest.TestCase):
    def test_tm_available(self) -> None:
        from ssdna import thermo

        if not thermo.available():
            self.skipTest("未安装 Biopython")
        self.assertIsNotNone(thermo.tm_nn("ACGTACGTACGTACG"))
        self.assertGreater(thermo.tm_nn("GCGCGCGCGCGCG"), thermo.tm_nn("ATATATATATATA"))

    def test_summary_marks_availability(self) -> None:
        from ssdna import thermo

        summary = thermo.summarize({"S01": {"tm_scaffold_domain": 40.0,
                                            "tm_target_domain": 60.0}})
        self.assertIn("available", summary)


if __name__ == "__main__":
    unittest.main()
