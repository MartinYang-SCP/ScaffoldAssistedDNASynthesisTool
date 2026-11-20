"""Task A 的验收测试：输入解析模块。

覆盖需求里点明的五类用例（线性提取、环状跨越提取、反向互补、非法输入、SnapGene 读取），
外加数据模型约束、三个导出文件、以及 Biopython ↔ 内置解析器的交叉校验。

运行：``python -m unittest discover -s tests -v``
"""

from __future__ import annotations

import csv
import json
import pathlib
import struct
import sys
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

# 本机把 Biopython 解在工作区的 build/pylibs（pip 无写权限）；有就加上，没有就跳过 bio 用例
_VENDORED = ROOT / "build" / "pylibs"
if _VENDORED.exists():
    sys.path.insert(0, str(_VENDORED))

try:  # noqa: SIM105
    import Bio  # noqa: F401,E402

    HAS_BIOPYTHON = True
except ImportError:  # pragma: no cover
    HAS_BIOPYTHON = False

from ssdna import (  # noqa: E402
    Scaffold,
    export_all,
    get_circular_subsequence,
    load_fasta,
    load_snapgene,
    normalize_sequence,
    reverse_complement,
)
from ssdna.models import Segment  # noqa: E402

DNA_FILE = ROOT / "minimalplasmid(1).dna"
OUTPUTS = ROOT / "outputs"
TMP = ROOT / "build" / "tmp"
TMP.mkdir(parents=True, exist_ok=True)


class TestReverseComplement(unittest.TestCase):
    def test_basic(self) -> None:
        self.assertEqual(reverse_complement("ATGC"), "GCAT")

    def test_palindrome(self) -> None:
        self.assertEqual(reverse_complement("GAATTC"), "GAATTC")

    def test_involution(self) -> None:
        for seq in ("ACGTN", "AAAACCCCGGGGTTTT", "atgc"):
            self.assertEqual(reverse_complement(reverse_complement(seq)), seq)

    def test_iupac_codes(self) -> None:
        # N 互补成 N，R（A/G）互补成 Y（C/T）
        self.assertEqual(reverse_complement("ACGTRY"), "RYACGT")

    def test_length_preserved(self) -> None:
        seq = "ACGTTGCA" * 150
        self.assertEqual(len(reverse_complement(seq)), len(seq))


class TestCircularSubsequence(unittest.TestCase):
    SEQ = "ACGTACGTAC"  # 长度 10

    def test_linear_inside(self) -> None:
        self.assertEqual(get_circular_subsequence(self.SEQ, 2, 4), "GTAC")

    def test_linear_full(self) -> None:
        self.assertEqual(get_circular_subsequence(self.SEQ, 0, 10), "ACGTACGTAC")

    def test_wrap_across_origin(self) -> None:
        # 8..10 + 0..2
        self.assertEqual(get_circular_subsequence(self.SEQ, 8, 4), "ACAC")
        self.assertEqual(
            get_circular_subsequence(self.SEQ, 8, 4),
            self.SEQ[8:] + self.SEQ[:2],
        )

    def test_wrap_many_starts(self) -> None:
        doubled = self.SEQ * 2
        for start in range(len(self.SEQ)):
            for length in range(0, len(self.SEQ) + 1):
                self.assertEqual(
                    get_circular_subsequence(self.SEQ, start, length),
                    doubled[start : start + length],
                )

    def test_zero_length(self) -> None:
        self.assertEqual(get_circular_subsequence(self.SEQ, 3, 0), "")

    def test_rejects_too_long(self) -> None:
        with self.assertRaises(IndexError):
            get_circular_subsequence(self.SEQ, 0, 11)

    def test_rejects_bad_start(self) -> None:
        with self.assertRaises(IndexError):
            get_circular_subsequence(self.SEQ, 10, 1)
        with self.assertRaises(IndexError):
            get_circular_subsequence(self.SEQ, -1, 1)

    def test_linear_mode_refuses_wrap(self) -> None:
        with self.assertRaises(IndexError):
            get_circular_subsequence(self.SEQ, 8, 4, circular=False)


class TestNormalizeAndCheck(unittest.TestCase):
    def test_strips_whitespace_and_digits(self) -> None:
        self.assertEqual(normalize_sequence("acgt\n  12 acgt\n"), "ACGTACGT")

    def test_case_folding(self) -> None:
        self.assertEqual(normalize_sequence("acgtACGT"), "ACGTACGT")

    def test_illegal_characters_raise(self) -> None:
        for bad in ("ACGTX", "AC GT?", "ACGT-ACGT", "ACGTZ"):
            with self.assertRaises(ValueError, msg=bad):
                normalize_sequence(bad)

    def test_empty_raises(self) -> None:
        for empty in ("", "   ", "12345"):
            with self.assertRaises(ValueError):
                normalize_sequence(empty)

    def test_degenerate_bases_are_legal(self) -> None:
        self.assertEqual(normalize_sequence("acgtnrykmswbdhv"), "ACGTNRYKMSWBDHV")

    def test_check_sequence_reports_degenerate(self) -> None:
        from ssdna import check_sequence

        self.assertTrue(check_sequence("ACGTNN"))
        self.assertEqual(check_sequence("ACGT"), [])


class TestScaffoldModel(unittest.TestCase):
    def test_length_mismatch_is_rejected(self) -> None:
        with self.assertRaises(ValueError):
            Scaffold(id="x", sequence="ACGT", length=5, topology="circular")

    def test_unknown_topology_is_rejected(self) -> None:
        with self.assertRaises(ValueError):
            Scaffold(id="x", sequence="ACGT", length=4, topology="triangular")

    def test_segment_wrap_length(self) -> None:
        seg = Segment(start=1181, end=37, wraps=True)
        self.assertEqual(seg.length_on(1200), 56)
        with self.assertRaises(ValueError):
            _ = seg.length


class TestSnapGeneRead(unittest.TestCase):
    def test_reads_the_uploaded_file(self) -> None:
        scaffold = load_snapgene(DNA_FILE, backend="raw")
        self.assertEqual(scaffold.length, 1200)
        self.assertEqual(scaffold.topology, "circular")
        self.assertEqual(len(scaffold.sequence), 1200)
        self.assertEqual(set(scaffold.sequence), set("ACGT"))
        self.assertEqual(scaffold.id, "minimalplasmid_1")  # FASTA 合法 ID：括号被归一成下划线

    def test_features_and_coordinates(self) -> None:
        scaffold = load_snapgene(DNA_FILE, backend="raw")
        table = {
            f.label: (f.feature_type, f.strand, f.segments[0].start, f.segments[0].end)
            for f in scaffold.features
        }
        self.assertEqual(table["ori"], ("rep_origin", "-", 78, 667))
        self.assertEqual(table["BleoR"], ("CDS", "-", 691, 1066))
        self.assertEqual(table["cexuf1"], ("misc_feature", ".", 1146, 1166))
        # MCS 在文件里被拆成两段（跨原点）
        mcs = [f for f in scaffold.features if f.label == "MCS"]
        self.assertEqual(len(mcs), 2)
        self.assertEqual(
            sorted((f.segments[0].start, f.segments[0].end) for f in mcs),
            [(0, 37), (1181, 1200)],
        )

    def test_primers_and_notes(self) -> None:
        scaffold = load_snapgene(DNA_FILE, backend="raw")
        self.assertEqual(
            sorted(p.name for p in scaffold.primers), sorted(["测序-f2", "F1", "F2", "R1", "R10"])
        )
        self.assertEqual(scaffold.notes.get("Type"), "Synthetic")
        self.assertEqual(scaffold.sequence_class, "UNA")

    def test_circular_subsequence_on_scaffold(self) -> None:
        scaffold = load_snapgene(DNA_FILE, backend="raw")
        joined = scaffold.subsequence(1181, 56)
        self.assertEqual(len(joined), 56)
        self.assertEqual(joined, scaffold.sequence[1181:] + scaffold.sequence[:37])

    def test_does_not_guess_scaffold_region(self) -> None:
        scaffold = load_snapgene(DNA_FILE, backend="raw")
        info = export_all(scaffold, TMP / "probe")
        payload = json.loads(pathlib.Path(info["info_json"]).read_text(encoding="utf-8"))
        self.assertIsNone(payload["scaffold_region_interpretation"]["scaffold_region"])
        self.assertIsNone(payload["scaffold_region_interpretation"]["scaffold_strand"])

    def test_rejects_non_snapgene(self) -> None:
        bogus = TMP / "bogus.dna"
        bogus.write_bytes(b"not a snapgene file")
        from ssdna.snapgene_raw import SnapGeneError

        with self.assertRaises(SnapGeneError):
            load_snapgene(bogus, backend="raw")

    @unittest.skipUnless(HAS_BIOPYTHON, "未安装 Biopython")
    def test_biopython_backend_matches_raw_backend(self) -> None:
        bio = load_snapgene(DNA_FILE, backend="biopython")
        raw = load_snapgene(DNA_FILE, backend="raw")
        self.assertEqual(bio.sequence, raw.sequence)
        self.assertEqual(bio.topology, raw.topology)
        self.assertEqual(
            sorted((f.label, f.segments[0].start, f.segments[0].end) for f in bio.features),
            sorted((f.label, f.segments[0].start, f.segments[0].end) for f in raw.features),
        )
        self.assertTrue(bio.cross_check["sequence_identical"])
        self.assertTrue(bio.cross_check["topology_identical"])
        self.assertTrue(bio.cross_check["features_identical"], bio.cross_check)

    @unittest.skipUnless(HAS_BIOPYTHON, "未安装 Biopython")
    def test_biopython_recovers_primer_binding_sites(self) -> None:
        bio = load_snapgene(DNA_FILE, backend="biopython")
        by_name = {p.name: p for p in bio.primers}
        self.assertEqual(len(by_name), 5)
        # R10 是跨原点的引物：join{[1086:1200](-), [0:6](-)}，长 120 nt
        r10 = by_name["R10"].site[0]
        self.assertTrue(r10.wraps)
        self.assertEqual(r10.length_on(1200), 120)
        self.assertEqual(by_name["R10"].site_origin, "derived")

    @unittest.skipUnless(HAS_BIOPYTHON, "未安装 Biopython")
    def test_auto_backend_uses_biopython(self) -> None:
        scaffold = load_snapgene(DNA_FILE)
        self.assertEqual(scaffold.parser_backend, "biopython+raw")


class TestFasta(unittest.TestCase):
    def test_roundtrip(self) -> None:
        path = TMP / "t.fa"
        path.write_text(">seq1 desc\nACGT\nACGT\n", encoding="ascii")
        scaffold = load_fasta(path)
        self.assertEqual(scaffold.sequence, "ACGTACGT")
        self.assertEqual(scaffold.topology, "linear")
        self.assertEqual(scaffold.id, "seq1")

    def test_multi_record_requires_id(self) -> None:
        path = TMP / "multi.fa"
        path.write_text(">a\nACGT\n>b\nTTTT\n", encoding="ascii")
        with self.assertRaises(ValueError):
            load_fasta(path)
        self.assertEqual(load_fasta(path, record_id="b").sequence, "TTTT")

    def test_not_fasta(self) -> None:
        path = TMP / "nofasta.txt"
        path.write_text("ACGT\n", encoding="ascii")
        with self.assertRaises(ValueError):
            load_fasta(path)


class TestExports(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.scaffold = load_snapgene(DNA_FILE, backend="raw")
        cls.paths = export_all(cls.scaffold, OUTPUTS)

    def test_three_files_exist(self) -> None:
        for key in ("fasta", "features_csv", "info_json"):
            self.assertTrue(pathlib.Path(self.paths[key]).exists(), key)
            self.assertGreater(pathlib.Path(self.paths[key]).stat().st_size, 0)

    def test_fasta_content(self) -> None:
        text = pathlib.Path(self.paths["fasta"]).read_text(encoding="ascii")
        lines = text.strip().splitlines()
        self.assertTrue(lines[0].startswith(">"))
        self.assertIn("minimalplasmid", lines[0])
        self.assertIn("length=1200", lines[0])
        self.assertIn("topology=circular", lines[0])
        body = "".join(lines[1:])
        self.assertEqual(body, self.scaffold.sequence)

    def test_csv_header_and_rows(self) -> None:
        with pathlib.Path(self.paths["features_csv"]).open(encoding="utf-8-sig") as fh:
            rows = list(csv.DictReader(fh))
        self.assertEqual(
            rows[0].keys(),
            {
                "feature_type", "label", "start", "end", "strand", "note", "kind",
                "segment_index", "segment_count", "length", "wraps", "origin",
                "raw_range", "sequence",
            },
        )
        annotations = [r for r in rows if r["kind"] == "annotation"]
        primers = [r for r in rows if r["kind"] == "primer"]
        self.assertEqual(len(annotations), 5)
        self.assertEqual(len(primers), 5)
        ori = next(r for r in annotations if r["label"] == "ori")
        self.assertEqual((ori["start"], ori["end"]), ("78", "667"))
        self.assertEqual(ori["raw_range"], "79-667")
        # primer 行在 raw backend 下没有回贴位点，但序列必须在
        self.assertTrue(all(r["sequence"] for r in primers))

    def test_json_fields(self) -> None:
        payload = json.loads(
            pathlib.Path(self.paths["info_json"]).read_text(encoding="utf-8")
        )
        self.assertEqual(payload["length"], 1200)
        self.assertEqual(payload["topology"], "circular")
        self.assertEqual(payload["molecule_type"], "DNA")
        self.assertEqual(len(payload["sequence"]), 1200)
        self.assertEqual(len(payload["features"]), 5)
        self.assertIn("0-based", payload["coordinate_convention"])
        self.assertEqual(payload["file_sha256"], self.scaffold.sha256)
        self.assertEqual(payload["sequence_sha256"], self.scaffold.sequence_sha256)


class TestSnapGeneContainer(unittest.TestCase):
    def test_sequence_block_length_matches_declared_length(self) -> None:
        raw = DNA_FILE.read_bytes()
        header_len = struct.unpack(">I", raw[1:5])[0]
        offset = 5 + header_len
        self.assertEqual(raw[offset], 0)  # type = 0
        block_len = struct.unpack(">I", raw[offset + 1 : offset + 5])[0]
        self.assertEqual(block_len, 1201)  # 1 字节标志 + 1200 nt


if __name__ == "__main__":
    unittest.main()
