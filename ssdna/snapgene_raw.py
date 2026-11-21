"""SnapGene ``.dna`` 文件的字节级解析器（不依赖第三方库）。

这是 SnappyGene 容器格式的**独立实现**，与 Biopython 的 ``SnapGeneIO`` 互为交叉校验：
两者对同一文件给出的序列、特征坐标、链别必须一致，否则 ``ssdna.io`` 会报警。

保留它的第二个原因：Biopython 不暴露 SnapGene 的 ``Notes`` 块（``Type=Synthetic``、
``ConfirmedExperimentally=0``、``SequenceClass=UNA``）与引物原始序列，这些信息对判断
"这条序列是不是实测过的"很关键。
"""

from __future__ import annotations

import gzip
import hashlib
import json
import re
import struct
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Any, Iterable

__all__ = [
    "SnapGeneError",
    "Segment",
    "Feature",
    "Primer",
    "SnapGeneFile",
    "read_snapgene",
    "parse_feature_ranges",
]

MAGIC = b"\x09"
SIGNATURE = b"SnapGene"
VALID_BASES = frozenset("ACGTUNRYKMSWBDHV")

BLOCK_NAMES = {
    0: "sequence",
    2: "features_legacy_binary",
    3: "features_legacy_binary_2",
    5: "primers_xml",
    6: "notes_xml",
    8: "additional_sequence_properties_xml",
    10: "features_xml",
    13: "enzyme_set",
    14: "custom_enzyme_set_xml",
    17: "alignable_sequences_xml",
    27: "history_gzip",
}


class SnapGeneError(ValueError):
    """文件不符合 SnapGene 容器格式。"""


# --------------------------------------------------------------------------- #
# 数据模型
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class Segment:
    """一个特征片段。``start``/``end`` 为 0-based 左闭右开。

    ``wraps`` 为 True 表示该片段跨过环状序列的原点（0 点），即 end <= start
    时按环长回绕解释；此时 ``raw`` 保留文件里的原始写法。
    """

    start: int
    end: int
    strand: str  # "+" / "-" / "."
    wraps: bool = False
    raw: str = ""

    @property
    def length(self) -> int:
        if self.wraps:
            raise ValueError("跨原点片段的长度需结合环长计算，见 length_on(circular_length)")
        return self.end - self.start

    def length_on(self, circular_length: int) -> int:
        if self.wraps:
            return (self.end - self.start) % circular_length
        return self.end - self.start


@dataclass
class Feature:
    name: str
    type: str
    strand: str
    segments: list[Segment] = field(default_factory=list)
    qualifiers: dict[str, str] = field(default_factory=dict)
    identifier: str | None = None


@dataclass
class Primer:
    name: str
    sequence: str
    description: str = ""
    date_added: str = ""
    identifier: str | None = None


@dataclass
class SnapGeneFile:
    source: str
    sha256: str
    length: int
    topology: str  # "circular" / "linear"
    sequence: str
    features: list[Feature] = field(default_factory=list)
    primers: list[Primer] = field(default_factory=list)
    notes: dict[str, str] = field(default_factory=dict)
    sequence_class: str = ""
    header_block_length: int = 19
    block_inventory: list[dict[str, Any]] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    # -- 便捷视图 ---------------------------------------------------------- #
    @property
    def is_circular(self) -> bool:
        return self.topology == "circular"

    def subsequence(self, start: int, end: int) -> str:
        """取 ``[start, end)``，环状序列自动回绕。不做大小写或字符修正。"""
        n = self.length
        length = end - start
        if length < 0:
            raise IndexError(f"区间 [{start}, {end}) 的 end 小于 start")
        if not self.is_circular:
            if not (0 <= start <= end <= n):
                raise IndexError(f"区间 [{start}, {end}) 超出线性序列长度 {n}")
            return self.sequence[start:end]
        if length > n:
            raise IndexError(f"区间长度 {length} 超过环长 {n}")
        s = start % n
        if s + length <= n:
            return self.sequence[s : s + length]
        return self.sequence[s:] + self.sequence[: s + length - n]

    def to_json(self, include_sequence: bool = True) -> dict[str, Any]:
        data = {
            "schema": "scaffold_staple_design/input_parse/v1",
            "source": self.source,
            "sha256": self.sha256,
            "topology": self.topology,
            "length": self.length,
            "sequence_class": self.sequence_class,
            "header_block_length": self.header_block_length,
            "block_inventory": self.block_inventory,
            "notes": self.notes,
            "features": [
                {
                    "name": f.name,
                    "type": f.type,
                    "strand": f.strand,
                    "segments": [asdict(s) for s in f.segments],
                    "qualifiers": f.qualifiers,
                }
                for f in self.features
            ],
            "primers": [
                {
                    "name": p.name,
                    "sequence": p.sequence,
                    "description": p.description,
                    "date_added": p.date_added,
                }
                for p in self.primers
            ],
            "warnings": self.warnings,
            "coordinate_convention": (
                "所有 start/end 均为 0-based 左闭右开半开区间 [start, end)。"
                "SnapGene 原始 1-based 闭区间字符串保留在 segments[].raw。"
            ),
            # 显式留空：本解析器不猜测 scaffold 的真实区域与链别。
            "interpretation": {
                "scaffold_region": None,
                "scaffold_strand": None,
                "note": (
                    "文件只提供一条环状双链质粒序列，无法从中判定实验室实际提供的"
                    "环状单链 scaffold 是哪一条链、覆盖哪一段。必须由用户在 params 中指定。"
                ),
            },
        }
        if include_sequence:
            data["sequence"] = self.sequence
        return data


# --------------------------------------------------------------------------- #
# 底层读取
# --------------------------------------------------------------------------- #
def _read_blocks(raw: bytes) -> tuple[int, list[tuple[int, int, bytes]]]:
    if not raw[:1] == MAGIC or raw[5:13] != SIGNATURE:
        raise SnapGeneError("不是 SnapGene .dna 文件：魔数或 'SnapGene' 签名不匹配")
    header_len = struct.unpack(">I", raw[1:5])[0]
    start = 5 + header_len
    if start > len(raw):
        raise SnapGeneError("头部块长度超出文件大小")

    blocks: list[tuple[int, int, bytes]] = []
    i = start
    while i < len(raw):
        if i + 5 > len(raw):
            raise SnapGeneError(f"偏移 {i} 处的数据块头部被截断")
        btype = raw[i]
        blen = struct.unpack(">I", raw[i + 1 : i + 5])[0]
        if i + 5 + blen > len(raw):
            raise SnapGeneError(
                f"偏移 {i} 处的 type={btype} 块声明长度 {blen} 超出文件末尾"
            )
        blocks.append((btype, i, raw[i + 5 : i + 5 + blen]))
        i += 5 + blen
    if not blocks:
        raise SnapGeneError("文件不含任何数据块")
    return start, blocks


def parse_feature_ranges(segment_element: ET.Element) -> Segment:
    """把 SnapGene ``<Segment range="79-667">`` 转成 0-based 半开区间。

    SnapGene 写的是 1-based 闭区间；跨原点时写成 ``end < start`` 的形式。
    """
    raw = segment_element.attrib.get("range")
    try:
        if raw is None:
            raise ValueError("缺少 range 属性")
        a_str, b_str = raw.split("-")
        a, b = int(a_str), int(b_str)
    except (ValueError, AttributeError) as exc:
        raise SnapGeneError(f"特征区间 {raw!r} 无法解析：{exc}") from exc
    # 1-based 闭区间 [a, b] -> 0-based 半开 [a-1, b)
    start, end = a - 1, b
    return Segment(start=start, end=end, strand=".", wraps=end <= start, raw=raw)



def _html_to_text(html: str) -> str:
    return re.sub(r"<[^>]+>", "", html).replace("&nbsp;", " ").strip()


def _parse_features(xml_bytes: bytes) -> list[Feature]:
    try:
        root = ET.fromstring(xml_bytes)
    except ET.ParseError as exc:
        raise SnapGeneError(f"特征 XML 块解析失败：{exc}") from exc
    features: list[Feature] = []
    for el in root.iter("Feature"):
        directionality = el.attrib.get("directionality", "")
        # 实测校验：directionality="2" 对应反向链（BleoR CDS 的翻译只有反向读框正确）
        strand = {"1": "+", "2": "-"}.get(directionality, ".")
        quals: dict[str, str] = {}
        for q in el.findall("Q"):
            qname = q.attrib.get("name", "")
            value_el = q.find("V")
            if value_el is None:
                continue
            if value_el.attrib.get("text") is not None:
                quals[qname] = _html_to_text(value_el.attrib["text"])
            elif value_el.attrib.get("int") is not None:
                quals[qname] = value_el.attrib["int"]
        segments = []
        for seg_el in el.findall("Segment"):
            seg = parse_feature_ranges(seg_el)
            segments.append(Segment(**{**asdict(seg), "strand": strand}))
        features.append(
            Feature(
                name=el.attrib.get("name", ""),
                type=el.attrib.get("type", ""),
                strand=strand,
                segments=segments,
                qualifiers=quals,
                identifier=el.attrib.get("recentID"),
            )
        )
    return features


def _parse_primers(xml_bytes: bytes) -> list[Primer]:
    try:
        root = ET.fromstring(xml_bytes)
    except ET.ParseError as exc:
        raise SnapGeneError(f"引物 XML 块解析失败：{exc}") from exc
    primers: list[Primer] = []
    for el in root.iter("Primer"):
        primers.append(
            Primer(
                name=el.attrib.get("name", ""),
                sequence=el.attrib.get("sequence", ""),
                description=_html_to_text(el.attrib.get("description", "")),
                date_added=el.attrib.get("dateAdded", ""),
                identifier=el.attrib.get("recentID"),
            )
        )
    return primers


def _parse_notes(xml_bytes: bytes) -> tuple[dict[str, str], str]:
    try:
        root = ET.fromstring(xml_bytes)
    except ET.ParseError as exc:
        raise SnapGeneError(f"备注 XML 块解析失败：{exc}") from exc
    notes: dict[str, str] = {}
    sequence_class = ""
    for child in root:
        text = (child.text or "").strip()
        key = child.tag
        if key == "SequenceClass":
            sequence_class = text
        notes[key] = text
        for attr, value in child.attrib.items():
            notes[f"{key}@{attr}"] = value
    return notes, sequence_class


# --------------------------------------------------------------------------- #
# 主入口
# --------------------------------------------------------------------------- #
def read_snapgene(path: str | Path) -> SnapGeneFile:
    """读取一个 SnapGene 文件，返回未加解释的原始事实。"""
    path = Path(path)
    raw = path.read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    header_len, blocks = _read_blocks(raw)

    record = SnapGeneFile(
        source=str(path),
        sha256=digest,
        length=0,
        topology="unknown",
        sequence="",
        header_block_length=header_len,
    )

    for btype, offset, data in blocks:
        record.block_inventory.append(
            {
                "type": btype,
                "name": BLOCK_NAMES.get(btype, "unknown"),
                "offset": offset,
                "length": len(data),
            }
        )
        if btype == 0:
            # 序列块：data[0] 为标志位（bit0 = 环状），其后紧跟序列本身，
            # 不含额外的长度字段；长度由块长度决定。
            if not data:
                raise SnapGeneError("type=0 序列块为空")
            flags = data[0]
            record.topology = "circular" if flags & 0x01 else "linear"
            try:
                record.sequence = data[1:].decode("ascii")
            except UnicodeDecodeError as exc:
                raise SnapGeneError(f"序列块含非 ASCII 字节：{exc}") from exc
            record.length = len(record.sequence)
        elif btype == 10:
            record.features = _parse_features(data)
        elif btype == 5:
            record.primers = _parse_primers(data)
        elif btype == 6:
            record.notes, record.sequence_class = _parse_notes(data)
        elif btype == 27:
            try:
                record.block_inventory[-1]["decompressed_length"] = len(
                    gzip.decompress(data)
                )
            except OSError:
                record.warnings.append("type=27 历史记录块不是有效的 gzip 数据")

    if not record.sequence:
        raise SnapGeneError("未找到序列数据块（type=0）")

    record.warnings.extend(validate(record))
    return record


def validate(record: SnapGeneFile) -> list[str]:
    """结构性自检；只报告，不修改任何数据。"""
    warnings: list[str] = []
    bad = sorted(set(record.sequence.upper()) - VALID_BASES)
    if bad:
        warnings.append(f"序列含非标准碱基字符：{bad}")
    if record.topology == "unknown":
        warnings.append("未解析出拓扑（环状/线性）")
    for feature in record.features:
        for seg in feature.segments:
            if seg.wraps:
                continue
            if seg.start < 0 or seg.end > record.length:
                warnings.append(
                    f"特征 {feature.name!r} 的区间 {seg.raw} 超出序列长度 {record.length}"
                )
    # 文件中记录的翻译长度与按区间翻译长度不一致时提示（不修改注释）
    for feature in record.features:
        tr = feature.qualifiers.get("translation", "")
        if tr and feature.segments and not feature.segments[0].wraps:
            nt_len = sum(s.length for s in feature.segments)
            expected = nt_len // 3
            if abs(len(tr) - expected) > 1:
                warnings.append(
                    f"特征 {feature.name!r} 记录翻译长 {len(tr)} aa，"
                    f"注释区间 {nt_len} nt 可译 {expected} aa，两者不完全一致"
                )
    return warnings


def _cli(argv: Iterable[str] | None = None) -> int:
    import argparse

    parser = argparse.ArgumentParser(description="只读解析 SnapGene .dna 文件")
    parser.add_argument("path", help=".dna 文件路径")
    parser.add_argument("-o", "--out", help="输出 JSON 路径（默认打印到标准输出）")
    parser.add_argument(
        "--no-sequence", action="store_true", help="输出中省略完整序列"
    )
    args = parser.parse_args(list(argv) if argv is not None else None)

    record = read_snapgene(args.path)
    payload = json.dumps(
        record.to_json(include_sequence=not args.no_sequence),
        ensure_ascii=False,
        indent=2,
    )
    if args.out:
        Path(args.out).write_text(payload, encoding="utf-8")
    else:
        print(payload)
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_cli())
