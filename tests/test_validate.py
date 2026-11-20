"""Task D 的验收测试：独立验证模块。

思路：先复制一份 **合法** 的 demo 设计，再逐项注入一种错误，
确认验证器能独立抓到，且错误里带 oligo id / 坐标 / 类型 / 建议。
"""

from __future__ import annotations

import csv
import json
import pathlib
import shutil
import sys
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
_VENDORED = ROOT / "build" / "pylibs"
if _VENDORED.exists():
    sys.path.insert(0, str(_VENDORED))

from ssdna import validate  # noqa: E402

DEMO = ROOT / "outputs" / "demo"
WORK = ROOT / "build" / "tmp" / "validate"


def fresh_bundle(name: str) -> pathlib.Path:
    dest = WORK / name
    if dest.exists():
        shutil.rmtree(dest)
    dest.mkdir(parents=True)
    for path in DEMO.iterdir():
        if path.is_file():
            shutil.copy2(path, dest / path.name)
    return dest


def read_csv(path: pathlib.Path) -> list[dict]:
    with path.open(encoding="utf-8-sig") as fh:
        return list(csv.DictReader(fh))


def write_csv(path: pathlib.Path, rows: list[dict], columns: list[str]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def run(directory: pathlib.Path) -> dict:
    return validate.run_validation(validate.load_bundle(directory))


def codes(payload: dict) -> set[str]:
    return {i["code"] for i in payload["issues"]}


class TestCleanBundlePasses(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.directory = fresh_bundle("clean")
        cls.payload = run(cls.directory)

    def test_status(self) -> None:
        self.assertEqual(self.payload["status"], "PASS")
        self.assertEqual(self.payload["counts"]["issues"]["error"], 0)

    def test_all_checks_ran(self) -> None:
        names = {c["name"] for c in self.payload["checks"]}
        self.assertEqual(
            names,
            {
                "sequence_basics", "target_coverage", "group_order",
                "boundary_spanning", "staple_length", "scaffold_mapping",
                "target_binding", "nicks",
            },
        )

    def test_coverage_is_complete(self) -> None:
        cov = self.payload["coverage"]
        self.assertEqual(cov["uncovered_bases"], 0)
        self.assertEqual(cov["overlapping_bases"], 0)
        self.assertEqual(cov["coverage_fraction"], 1.0)

    def test_three_outputs(self) -> None:
        out = self.directory / "validation"
        for name in ("validation_report.json", "validation_summary.md",
                     "potential_conflicts.csv"):
            self.assertTrue((out / name).exists(), name)


class TestSequenceLevelChecks(unittest.TestCase):
    def test_empty_and_illegal(self) -> None:
        self.assertTrue(validate.validate_sequence(""))
        self.assertTrue(validate.validate_sequence("ACGTX"))
        self.assertTrue(validate.validate_sequence("acgt"))

    def test_length(self) -> None:
        self.assertEqual(validate.validate_length("ACGT", 4), [])
        issues = validate.validate_length("ACG", 4)
        self.assertEqual(issues[0].code, "length_mismatch")
        self.assertIn("1", issues[0].message)

    def test_reverse_complement(self) -> None:
        self.assertEqual(validate.validate_reverse_complement("ATGC", "GCAT"), [])
        # 同向互补（忘了反向）
        same = validate.validate_reverse_complement("ATGC", "TACG")
        self.assertEqual(same[0].code, "same_orientation_complement")
        # 完全一样
        ident = validate.validate_reverse_complement("ATGC", "ATGC")
        self.assertEqual(ident[0].code, "identical_not_complementary")

    def test_issue_has_all_fields(self) -> None:
        issue = validate.validate_reverse_complement("ATGC", "TTTT")[0]
        payload = issue.to_dict()
        for key in ("level", "code", "subject", "coordinates", "message", "suggestion"):
            self.assertIn(key, payload)
        self.assertTrue(payload["suggestion"])


class TestTamperedOligo(unittest.TestCase):
    def test_sequence_mismatch(self) -> None:
        directory = fresh_bundle("bad_oligo_seq")
        rows = read_csv(directory / "target_oligos.csv")
        rows[0]["sequence"] = rows[0]["sequence"][:-1] + (
            "A" if rows[0]["sequence"][-1] != "A" else "C")
        write_csv(directory / "target_oligos.csv", rows, list(rows[0]))
        payload = run(directory)
        self.assertIn("oligo_sequence_mismatch", codes(payload))
        self.assertEqual(payload["status"], "FAIL")
        bad = [i for i in payload["issues"] if i["code"] == "oligo_sequence_mismatch"][0]
        self.assertEqual(bad["subject"], "T01")
        self.assertTrue(bad["coordinates"])

    def test_wrong_strand(self) -> None:
        directory = fresh_bundle("bad_strand")
        rows = read_csv(directory / "target_oligos.csv")
        rows[0]["strand"] = "-"
        write_csv(directory / "target_oligos.csv", rows, list(rows[0]))
        payload = run(directory)
        self.assertIn("oligo_wrong_strand", codes(payload))


class TestCoverage(unittest.TestCase):
    def test_dropped_oligo_leaves_gap(self) -> None:
        directory = fresh_bundle("missing_oligo")
        rows = read_csv(directory / "target_oligos.csv")
        kept = [r for r in rows if r["oligo_id"] != "T02"]
        write_csv(directory / "target_oligos.csv", kept, list(rows[0]))
        payload = run(directory)
        self.assertIn("uncovered_bases", codes(payload))
        self.assertGreater(payload["coverage"]["uncovered_bases"], 0)

    def test_overlap_is_detected(self) -> None:
        directory = fresh_bundle("overlap")
        rows = read_csv(directory / "target_oligos.csv")
        rows[1]["target_start"] = str(int(rows[1]["target_start"]) - 10)
        rows[1]["length"] = str(int(rows[1]["length"]) + 10)
        write_csv(directory / "target_oligos.csv", rows, list(rows[0]))
        payload = run(directory)
        self.assertIn("overlapping_bases", codes(payload))

    def test_group_swapped(self) -> None:
        directory = fresh_bundle("group_order")
        rows = read_csv(directory / "target_oligos.csv")
        for row in rows:
            row["group"] = "B" if row["group"] == "A" else "A"
        write_csv(directory / "target_oligos.csv", rows, list(rows[0]))
        payload = run(directory)
        self.assertIn("group_order", codes(payload))


class TestBoundaryOligos(unittest.TestCase):
    def test_boundary_label_removed(self) -> None:
        directory = fresh_bundle("boundary_label")
        rows = read_csv(directory / "target_oligos.csv")
        for row in rows:
            if row["oligo_type"] == "boundary_spanning":
                row["oligo_type"] = "ordinary"
                row["crosses_junction"] = "False"
                break
        write_csv(directory / "target_oligos.csv", rows, list(rows[0]))
        payload = run(directory)
        self.assertIn("boundary_count", codes(payload))
        self.assertIn("junction_mislabeled", codes(payload))

    def test_junction_covered_twice(self) -> None:
        directory = fresh_bundle("junction_twice")
        report = json.loads((directory / "design_report.json").read_text(encoding="utf-8"))
        report["junctions"].append({"junction_id": "JX", "position": 100,
                                    "upstream_group": "A", "downstream_group": "B"})
        (directory / "design_report.json").write_text(
            json.dumps(report, ensure_ascii=False), encoding="utf-8")
        payload = run(directory)
        self.assertIn("boundary_count", codes(payload))


class TestStapleChecks(unittest.TestCase):
    def test_tampered_staple_arm(self) -> None:
        directory = fresh_bundle("bad_staple")
        rows = read_csv(directory / "staples.csv")
        rows[0]["sequence"] = rows[0]["sequence"][:-1] + (
            "A" if rows[0]["sequence"][-1] != "A" else "C")
        write_csv(directory / "staples.csv", rows, list(rows[0]))
        payload = run(directory)
        self.assertIn("not_reverse_complement", codes(payload))

    def test_duplicate_scaffold_site(self) -> None:
        directory = fresh_bundle("dup_site")
        rows = read_csv(directory / "staples.csv")
        rows[1]["scaffold_start"] = rows[0]["scaffold_start"]
        write_csv(directory / "staples.csv", rows, list(rows[0]))
        payload = run(directory)
        self.assertIn("scaffold_site_collision", codes(payload))

    def test_staple_length(self) -> None:
        directory = fresh_bundle("short_staple")
        rows = read_csv(directory / "staples.csv")
        rows[0]["sequence"] = rows[0]["sequence"][:38]
        write_csv(directory / "staples.csv", rows, list(rows[0]))
        payload = run(directory)
        self.assertIn("length_mismatch", codes(payload))

    def test_overlapping_target_arms(self) -> None:
        directory = fresh_bundle("target_overlap")
        rows = read_csv(directory / "staples.csv")
        rows[1]["target_start"] = str(int(rows[0]["target_start"]) + 13)
        write_csv(directory / "staples.csv", rows, list(rows[0]))
        payload = run(directory)
        self.assertIn("target_site_overlap", codes(payload))


class TestNicks(unittest.TestCase):
    def test_missing_record(self) -> None:
        directory = fresh_bundle("missing_nick")
        rows = read_csv(directory / "nicks.csv")
        write_csv(directory / "nicks.csv", rows[1:], list(rows[0]))
        payload = run(directory)
        self.assertIn("nick_count", codes(payload))
        self.assertIn("nick_not_recorded", codes(payload))

    def test_wrong_end_direction(self) -> None:
        directory = fresh_bundle("wrong_end")
        rows = read_csv(directory / "nicks.csv")
        rows[0]["upstream_end"] = "5'"
        write_csv(directory / "nicks.csv", rows, list(rows[0]))
        payload = run(directory)
        self.assertIn("nick_wrong_end", codes(payload))

    def test_ligation_without_phosphate(self) -> None:
        directory = fresh_bundle("ligation_no_p")
        rows = read_csv(directory / "nicks.csv")
        for row in rows:
            row["intended_ligation"] = "True"
        write_csv(directory / "nicks.csv", rows, list(rows[0]))
        payload = run(directory)
        self.assertIn("ligation_without_phosphate", codes(payload))

    def test_missing_phosphate_is_only_a_warning(self) -> None:
        directory = fresh_bundle("phos_unknown")
        rows = read_csv(directory / "nicks.csv")
        for row in rows:
            row["five_prime_phosphate"] = "unknown"
        write_csv(directory / "nicks.csv", rows, list(rows[0]))
        payload = run(directory)
        self.assertIn("phosphate_unknown", codes(payload))
        self.assertEqual(payload["status"], "PASS")  # 提示不算失败


class TestThermoInterfaces(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.directory = fresh_bundle("thermo")
        cls.payload = run(cls.directory)

    def test_status_is_not_run(self) -> None:
        for tool in ("idt_oligoanalyzer", "nupack"):
            self.assertEqual(self.payload["thermo"][tool]["status"], "not_run")

    def test_input_files_exist(self) -> None:
        out = self.directory / "validation"
        self.assertTrue((out / "idt_oligoanalyzer_input.csv").exists())
        self.assertTrue((out / "nupack_strands.txt").exists())

    def test_no_fabricated_numbers(self) -> None:
        text = json.dumps(self.payload["thermo"], ensure_ascii=False)
        self.assertNotIn("hairpin_score", text)
        self.assertNotIn("delta_g", text.lower())
        self.assertIn("不能验证连接化学", text)

    def test_nupack_lines_have_name_and_sequence(self) -> None:
        out = self.directory / "validation"
        for line in (out / "nupack_strands.txt").read_text(encoding="ascii").splitlines():
            name, seq = line.split()
            self.assertTrue(name)
            self.assertTrue(set(seq) <= set("ACGT"))


class TestIndependence(unittest.TestCase):
    def test_validator_does_not_import_design(self) -> None:
        source = (ROOT / "ssdna" / "validate.py").read_text(encoding="utf-8")
        self.assertNotIn("from .design", source)
        self.assertNotIn("import design", source)

    def test_works_without_design_report(self) -> None:
        """没有 design_report.json 也能验（交界信息缺失时只报边界问题）。"""
        directory = fresh_bundle("no_report")
        (directory / "design_report.json").unlink()
        payload = run(directory)
        self.assertIn("status", payload)
        self.assertGreater(payload["counts"]["oligos"], 0)

    def test_refuses_when_inputs_missing(self) -> None:
        directory = fresh_bundle("empty_bundle")
        for path in directory.iterdir():
            path.unlink()
        with self.assertRaises(FileNotFoundError):
            validate.load_bundle(directory)


if __name__ == "__main__":
    unittest.main()
