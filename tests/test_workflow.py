"""Task E 的验收测试：完整命令行工作流。

覆盖任务书要求的十类测试（其余几类在 test_mapping / test_design / test_validate 里，
README 有对照表），外加端到端的 ``python main.py`` 行为。
"""

from __future__ import annotations

import json
import pathlib
import shutil
import subprocess
import sys
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
_VENDORED = ROOT / "build" / "pylibs"
if _VENDORED.exists():
    sys.path.insert(0, str(_VENDORED))

import main as cli  # noqa: E402

DEMO_CONFIG = ROOT / "config" / "design_config.demo.yaml"
TEMPLATE_CONFIG = ROOT / "config" / "design_config.yaml"
WORK = ROOT / "build" / "tmp" / "workflow"


def _clean(path: pathlib.Path) -> None:
    if not path.exists():
        return
    for item in sorted(path.rglob("*"), reverse=True):
        if item.is_file():
            item.unlink()
        elif item.is_dir():
            try:
                item.rmdir()
            except OSError:
                pass


class TestEndToEnd(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        _clean(WORK)
        WORK.mkdir(parents=True, exist_ok=True)
        config = DEMO_CONFIG.read_text(encoding="utf-8")
        config = config.replace("outdir: outputs/demo", f"outdir: {WORK.as_posix()}")
        config = config.replace("data/scaffold.dna", str(ROOT / "data" / "scaffold.dna"))
        config = config.replace("data/target.fasta", str(ROOT / "data" / "target.fasta"))
        cls.config = WORK / "config.yaml"
        cls.config.write_text(config, encoding="utf-8")
        cls.code = cli.run(cls.config, verbose=False)

    def test_exit_code(self) -> None:
        self.assertEqual(self.code, 0)

    def test_all_outputs(self) -> None:
        for name in (
            "target_oligos.csv", "staples.csv", "nicks.csv",
            "design_config.yaml", "design_report.json",
            "target.fasta", "scaffold_working.fasta",
        ):
            self.assertTrue((WORK / name).exists(), name)

    def test_validation_outputs(self) -> None:
        out = WORK / "validation"
        for name in ("validation_report.json", "validation_summary.md",
                     "potential_conflicts.csv"):
            self.assertTrue((out / name).exists(), name)
        report = json.loads((out / "validation_report.json").read_text(encoding="utf-8"))
        self.assertEqual(report["status"], "PASS")

    def test_design_report(self) -> None:
        payload = json.loads((WORK / "design_report.json").read_text(encoding="utf-8"))
        self.assertEqual(payload["status"], "OK")
        self.assertEqual(payload["counts"]["boundary_oligos"], 2)
        self.assertEqual(payload["counts"]["staples_A"], payload["counts"]["staples_B"])


class TestMissingParamsStops(unittest.TestCase):
    def test_template_stops_before_reading_anything(self) -> None:
        outdir = ROOT / "outputs" / "workflow_should_not_exist"
        if outdir.exists():
            _clean(outdir)
        code = cli.run(TEMPLATE_CONFIG, verbose=False)
        self.assertEqual(code, 2)
        self.assertFalse(outdir.exists())

    def test_missing_list_is_complete(self) -> None:
        from ssdna import design

        cfg = design.load_config(TEMPLATE_CONFIG)
        missing = design.find_missing(cfg)
        for key in ("scaffold.strand", "target.source", "target.group_boundaries",
                    "target.oligo_cuts", "target.five_prime_phosphate"):
            self.assertIn(key, missing)


class TestSubprocessEntryPoint(unittest.TestCase):
    """``python main.py`` 必须能直接跑（不能只靠 import）。"""

    def test_cli_runs(self) -> None:
        result = subprocess.run(
            [sys.executable, "main.py", "--config", str(DEMO_CONFIG), "--quiet"],
            cwd=ROOT, capture_output=True, text=True, encoding="utf-8",
            errors="replace", timeout=180,
        )
        self.assertIn(result.returncode, (0, 1), result.stderr[-2000:])

    def test_cli_reports_missing(self) -> None:
        result = subprocess.run(
            [sys.executable, "main.py", "--config", str(TEMPLATE_CONFIG)],
            cwd=ROOT, capture_output=True, text=True, encoding="utf-8",
            errors="replace", timeout=180,
        )
        self.assertEqual(result.returncode, 2)
        self.assertIn("缺少", result.stdout)


class TestDataAssets(unittest.TestCase):
    def test_scaffold_asset_matches_the_original(self) -> None:
        import hashlib

        data = (ROOT / "data" / "scaffold.dna").read_bytes()
        original = (ROOT / "minimalplasmid(1).dna").read_bytes()
        self.assertEqual(hashlib.sha256(data).hexdigest(),
                         hashlib.sha256(original).hexdigest())

    def test_target_fasta_is_labelled_synthetic(self) -> None:
        text = (ROOT / "data" / "target.fasta").read_text(encoding="ascii")
        self.assertTrue(text.startswith(">"))
        self.assertIn("synthetic", text.splitlines()[0])
        body = "".join(line for line in text.splitlines() if not line.startswith(">"))
        self.assertEqual(len(body), 728)

    def test_config_templates_present(self) -> None:
        self.assertTrue(TEMPLATE_CONFIG.exists())
        self.assertTrue(DEMO_CONFIG.exists())


if __name__ == "__main__":
    unittest.main()
