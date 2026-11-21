"""``tools/snapgene_reader`` 的单元测试。

运行： ``python -m unittest discover -s tests -v``
"""

from __future__ import annotations

import contextlib
import hashlib
import io
import json
import pathlib
import shutil
import sys
import uuid
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from tools import check_params  # noqa: E402
from tools.snapgene_reader import (  # noqa: E402
    SnapGeneError,
    read_snapgene,
    validate,
)

DNA_FILE = ROOT / "minimalplasmid(1).dna"

# 临时文件放在工作区内：受 DSH 工作区写入策略约束，系统临时目录未必可写。
TMP_ROOT = ROOT / "build" / "tmp"
TMP_ROOT.mkdir(parents=True, exist_ok=True)

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L",
    "CTA": "L", "CTG": "L", "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V", "TCT": "S", "TCC": "S",
    "TCA": "S", "TCG": "S", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "GCT": "A", "GCC": "A",
    "GCA": "A", "GCG": "A", "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q", "AAT": "N", "AAC": "N",
    "AAA": "K", "AAG": "K", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W", "CGT": "R", "CGC": "R",
    "CGA": "R", "CGG": "R", "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def translate(seq: str) -> str:
    return "".join(CODON_TABLE.get(seq[i : i + 3], "?") for i in range(0, len(seq) - 2, 3))


class TestSnapGeneFacts(unittest.TestCase):
    """对上传文件本身的事实性断言。"""

    @classmethod
    def setUpClass(cls) -> None:
        cls.record = read_snapgene(DNA_FILE)
        cls.by_name = {}
        for feature in cls.record.features:
            cls.by_name.setdefault(feature.name, []).append(feature)

    def test_sha256_matches_attachment(self) -> None:
        self.assertEqual(
            self.record.sha256,
            hashlib.sha256(DNA_FILE.read_bytes()).hexdigest(),
        )

    def test_length_and_topology(self) -> None:
        self.assertEqual(self.record.length, 1200)
        self.assertEqual(self.record.topology, "circular")

    def test_alphabet_is_canonical_dna(self) -> None:
        self.assertEqual(sorted(set(self.record.sequence)), ["A", "C", "G", "T"])

    def test_first_and_last_bases(self) -> None:
        # 文件编号原点落在 MCS 内部，首碱基为 T、末碱基为 G
        self.assertTrue(self.record.sequence.startswith("TCGACTCTAGAGGATCCCGGG"))
        self.assertTrue(self.record.sequence.endswith("TGCATGCCTGCAGG"))

    def test_feature_names(self) -> None:
        # 两条 MCS 记录同名（SnapGene 里跨原点的 MCS 被拆成两条），按名称归并后为 4 组
        self.assertEqual(sorted(self.by_name), ["BleoR", "MCS", "cexuf1", "ori"])
        self.assertEqual(len(self.by_name["MCS"]), 2)
        self.assertEqual(len(self.record.features), 5)

    def test_ori_coordinates(self) -> None:
        (ori,) = self.by_name["ori"]
        self.assertEqual(ori.strand, "-")
        self.assertEqual((ori.segments[0].start, ori.segments[0].end), (78, 667))
        self.assertEqual(ori.segments[0].length, 589)

    def test_cexuf1_coordinates(self) -> None:
        (feat,) = self.by_name["cexuf1"]
        self.assertEqual((feat.segments[0].start, feat.segments[0].end), (1146, 1166))
        self.assertEqual(feat.segments[0].length, 20)

    def test_mcs_is_split_across_the_origin(self) -> None:
        # MCS 在文件里被写成两段：1-37 与 1182-1200，合起来跨越 0 点
        mcs = sorted(self.by_name["MCS"], key=lambda f: f.segments[0].start)
        coords = [(f.segments[0].start, f.segments[0].end) for f in mcs]
        self.assertEqual(coords, [(0, 37), (1181, 1200)])
        self.assertEqual(sum(f.segments[0].length for f in mcs), 56)
        # 跨原点拼接后应连续：环上 [1181, 1200) + [0, 37) 是同一段
        joined = self.record.subsequence(1181, 1200 + 37)
        self.assertEqual(len(joined), 56)
        self.assertEqual(joined, self.record.sequence[1181:] + self.record.sequence[:37])

    def test_bleor_translation_only_works_on_reverse_strand(self) -> None:
        """这条断言同时锁定了两件事：链别语义（directionality=2 -> "-"）
        与坐标约定（1-based 闭区间 -> 0-based 半开）。"""
        (bleor,) = self.by_name["BleoR"]
        seg = bleor.segments[0]
        self.assertEqual((seg.start, seg.end), (691, 1066))
        self.assertEqual(bleor.strand, "-")
        fragment = self.record.sequence[seg.start : seg.end]
        self.assertEqual(len(fragment), 375)
        expected = bleor.qualifiers["translation"]
        self.assertEqual(translate(revcomp(fragment)), expected[:125])
        self.assertNotEqual(translate(fragment), expected[:125])

    def test_primers_are_cloning_primers_not_staples(self) -> None:
        names = sorted(p.name for p in self.record.primers)
        self.assertEqual(names, sorted(["测序-f2", "F1", "F2", "R1", "R10"]))
        for primer in self.record.primers:
            self.assertGreaterEqual(len(primer.sequence), 20)
            self.assertLessEqual(len(primer.sequence), 120)

    def test_notes(self) -> None:
        self.assertEqual(self.record.notes["Type"], "Synthetic")
        self.assertEqual(self.record.notes["ConfirmedExperimentally"], "0")
        self.assertEqual(self.record.sequence_class, "UNA")

    def test_interpretation_is_not_guessed(self) -> None:
        payload = self.record.to_json(include_sequence=False)
        self.assertIsNone(payload["interpretation"]["scaffold_region"])
        self.assertIsNone(payload["interpretation"]["scaffold_strand"])

    def test_json_roundtrip(self) -> None:
        payload = json.loads(json.dumps(self.record.to_json(), ensure_ascii=False))
        self.assertEqual(payload["length"], 1200)
        self.assertEqual(len(payload["features"]), 5)


class TestCircularAccess(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.record = read_snapgene(DNA_FILE)

    def test_wrap_slicing_matches_naive_concatenation(self) -> None:
        doubled = self.record.sequence * 2
        for start in range(0, 1200, 137):
            for length in (1, 13, 26, 39, 300):
                self.assertEqual(
                    self.record.subsequence(start, start + length),
                    doubled[start : start + length],
                )

    def test_full_circle_slice(self) -> None:
        self.assertEqual(self.record.subsequence(0, 1200), self.record.sequence)
        # 起点 1100 加 300 nt，回绕后应等于 [1100,1200) + [0,200)
        self.assertEqual(
            self.record.subsequence(1100, 1400),
            self.record.sequence[1100:] + self.record.sequence[:200],
        )

    def test_out_of_range(self) -> None:
        with self.assertRaises(IndexError):
            self.record.subsequence(0, 1201)


class TestRobustness(unittest.TestCase):
    """坏文件必须报错，好文件必须原样读出。

    临时文件直接写在工作区的 build/tmp 下：本机 ``tempfile.TemporaryDirectory``
    建出的目录会被 DSH 的 Windows 沙箱判为不可写/不可列举，与其绕开权限不如
    换个位置，测试仍然验证同一件事。
    """

    def _scratch(self, name: str) -> pathlib.Path:
        path = TMP_ROOT / f"{uuid.uuid4().hex[:12]}_{name}"
        self.addCleanup(lambda: path.unlink(missing_ok=True))
        return path

    def test_rejects_non_snapgene_file(self) -> None:
        bogus = self._scratch("x.dna")
        bogus.write_bytes(b"not a snapgene file at all")
        with self.assertRaises(SnapGeneError):
            read_snapgene(bogus)

    def test_rejects_truncated_block(self) -> None:
        truncated = self._scratch("t.dna")
        truncated.write_bytes(DNA_FILE.read_bytes()[:200])
        with self.assertRaises(SnapGeneError):
            read_snapgene(truncated)

    def test_reader_does_not_modify_input(self) -> None:
        copy = self._scratch("copy.dna")
        shutil.copyfile(DNA_FILE, copy)
        before = hashlib.sha256(copy.read_bytes()).hexdigest()
        read_snapgene(copy)
        after = hashlib.sha256(copy.read_bytes()).hexdigest()
        self.assertEqual(before, after)

    def test_sequence_is_returned_verbatim(self) -> None:
        record = read_snapgene(DNA_FILE)
        raw = DNA_FILE.read_bytes()
        # 序列块前 24 字节后即为序列本体，逐字节比对，证明没有做任何改写
        self.assertIn(record.sequence.encode("ascii"), raw)


class TestValidate(unittest.TestCase):
    def test_uploaded_file_has_no_structural_warnings(self) -> None:
        record = read_snapgene(DNA_FILE)
        self.assertEqual(validate(record), [])
        self.assertEqual(record.warnings, [])

    def test_detects_translation_length_drift(self) -> None:
        record = read_snapgene(DNA_FILE)
        for feature in record.features:
            if feature.name == "BleoR":
                feature.qualifiers["translation"] = "M" * 20
        warnings = validate(record)
        self.assertEqual(len(warnings), 1, warnings)
        self.assertIn("翻译", warnings[0])


class TestCheckParams(unittest.TestCase):
    """权威配置只有 config/design_config.yaml；缺参/非法参都会被拦截。"""

    def setUp(self) -> None:
        from ssdna.design import load_config as _lc
        self.config = _lc(ROOT / "config" / "design_config.yaml")

    def test_template_is_blocked(self) -> None:
        blocked = check_params.blocking(self.config)
        self.assertIn("scaffold.strand", blocked)
        self.assertIn("scaffold.region", blocked)
        self.assertIn("target.source", blocked)
        self.assertIn("target.oligo_cuts", blocked)

    def test_demo_config_is_clean(self) -> None:
        from ssdna.design import load_config as _lc
        demo = _lc(ROOT / "config" / "design_config.demo.yaml")
        self.assertEqual(check_params.blocking(demo), [])

    def test_invalid_values_are_blocked(self) -> None:
        bad = dict(self.config)
        bad["staple"] = dict(self.config["staple"], scaffold_binding_length="abc")
        self.assertIn("staple.scaffold_binding_length", check_params.blocking(bad))

    def test_cli_exit_code_is_two_when_blocked(self) -> None:
        buffer = io.StringIO()
        with contextlib.redirect_stdout(buffer):
            code = check_params.main([str(ROOT / "config" / "design_config.yaml")])
        self.assertEqual(code, 2)
        self.assertIn("阻塞项", buffer.getvalue())


if __name__ == "__main__":
    unittest.main()
