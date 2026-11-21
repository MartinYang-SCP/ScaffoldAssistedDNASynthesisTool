"""输入层：把 ``.dna`` / ``.fasta`` 读成 :class:`ssdna.models.Scaffold`。

SnapGene 文件走双通道：

1. **Biopython**（``Bio.SeqIO.read(path, "snapgene")``）解析序列、拓扑、特征；
2. **内置字节解析器**（:mod:`ssdna.snapgene_raw`）解析 Biopython 不暴露的部分
   （``Notes`` 块、引物原始序列、数据块清单），并**交叉校验**两者的结果。

两者对序列或坐标不一致时：序列不一致直接报错；特征不一致给出显式警告。
Biopython 不可用时自动退回内置解析器，并在 ``parser_backend`` 与 ``warnings`` 里说清楚。
"""

from __future__ import annotations

import hashlib
import pathlib
import re
from typing import Any

from . import snapgene_raw
from .models import Feature, Primer, Scaffold, Segment
from .sequence import check_sequence, get_circular_subsequence, normalize_sequence

__all__ = ["BackendUnavailable", "load_snapgene", "load_fasta"]



def _sanitize_id(text: str) -> str:
    """把文件名转成合法的 FASTA ID（只保留字母数字与 _ . -）。"""
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", text).strip("_")
    return cleaned or "scaffold"


class BackendUnavailable(ImportError):
    """请求的解析后端不可用。"""


# --------------------------------------------------------------------------- #
# Biopython 通道
# --------------------------------------------------------------------------- #
def _require_biopython():
    try:
        from Bio import SeqIO  # noqa: F401
    except ImportError as exc:  # pragma: no cover - 取决于环境
        raise BackendUnavailable(
            "未安装 Biopython。安装方式：pip install biopython；"
            "或者把 --backend 设为 raw 使用内置解析器。"
        ) from exc
    from Bio import SeqIO

    return SeqIO


def _location_segments(location: Any, length: int) -> tuple[list[Segment], bool]:
    """把 Biopython 的 Location 转成我们的 Segment 列表。"""
    parts = list(getattr(location, "parts", []) or [location])
    strand_value = getattr(location, "strand", None)
    strand = {1: "+", -1: "-"}.get(strand_value, ".")

    segments: list[Segment] = []
    for part in parts:
        start = int(part.start)
        end = int(part.end)
        segments.append(
            Segment(start=start, end=end, strand=strand, wraps=False, raw="")
        )

    # 跨原点：join{[1086:1200], [0:6]} 这种形式，第一段起点大于后续段起点
    wraps = len(segments) > 1 and segments[0].start > segments[-1].start
    if wraps and segments:
        first = segments[0]
        segments[0] = Segment(
            start=first.start,
            end=segments[-1].end if segments[-1].end < first.start else first.end,
            strand=strand,
            wraps=True,
            raw="",
        )
        segments = segments[:1]
    return segments, wraps


def _from_biopython(path: pathlib.Path) -> dict[str, Any]:
    SeqIO = _require_biopython()
    record = SeqIO.read(str(path), "snapgene")

    features: list[Feature] = []
    primer_sites: dict[str, list[Segment]] = {}
    for feature in record.features:
        label = ""
        if feature.qualifiers.get("label"):
            label = str(feature.qualifiers["label"][0])
        qualifiers = {}
        for key, value in feature.qualifiers.items():
            if key == "label":
                continue
            if isinstance(value, (list, tuple)):
                qualifiers[key] = "; ".join(str(v) for v in value)
            else:
                qualifiers[key] = str(value)
        segments, _ = _location_segments(feature.location, len(record.seq))

        if feature.type == "primer_bind":
            primer_sites[label] = segments
            continue

        features.append(
            Feature(
                feature_type=feature.type,
                label=label,
                strand=segments[0].strand if segments else ".",
                segments=segments,
                qualifiers=qualifiers,
                kind="annotation",
                origin="file",
            )
        )

    return {
        "sequence": str(record.seq),
        "topology": record.annotations.get("topology", "linear"),
        "molecule_type": record.annotations.get("molecule_type", "DNA"),
        "features": features,
        "primer_sites": primer_sites,
        "annotations": dict(record.annotations),
    }


# --------------------------------------------------------------------------- #
# 内置通道
# --------------------------------------------------------------------------- #
def _from_raw(path: pathlib.Path) -> dict[str, Any]:
    record = snapgene_raw.read_snapgene(path)
    features = [
        Feature(
            feature_type=f.type,
            label=f.name,
            strand=f.strand,
            segments=[
                Segment(
                    start=s.start,
                    end=s.end,
                    strand=s.strand,
                    wraps=s.wraps,
                    raw=s.raw,
                    origin="file",
                )
                for s in f.segments
            ],
            qualifiers=dict(f.qualifiers),
            kind="annotation",
            origin="file",
        )
        for f in record.features
    ]
    return {
        "sequence": record.sequence,
        "topology": record.topology,
        "molecule_type": "DNA",
        "features": features,
        "primers": [Primer(name=p.name, sequence=p.sequence) for p in record.primers],
        "notes": dict(record.notes),
        "sequence_class": record.sequence_class,
        "block_inventory": list(record.block_inventory),
        "raw_warnings": list(record.warnings),
    }


# --------------------------------------------------------------------------- #
# 交叉校验
# --------------------------------------------------------------------------- #
def _index(features: list[Feature]) -> dict[tuple, Feature]:
    index: dict[tuple, Feature] = {}
    for feature in features:
        key = (
            feature.feature_type,
            feature.label,
            tuple((s.start, s.end) for s in feature.segments),
        )
        index[key] = feature
    return index


def _strands_relaxed(a: str, b: str) -> bool:
    """文件没写方向（"."）而 Biopython 默认成 "+"，算可解释差异。"""
    return {a, b} == {".", "+"}


def _cross_check(bio: dict[str, Any], raw: dict[str, Any]) -> dict[str, Any]:
    bio_index = _index(bio["features"])
    raw_index = _index(raw["features"])

    common = bio_index.keys() & raw_index.keys()
    only_bio = sorted(bio_index.keys() - raw_index.keys())
    only_raw = sorted(raw_index.keys() - bio_index.keys())
    relaxed = [k for k in common if bio_index[k].strand != raw_index[k].strand]
    hard = [k for k in relaxed if not _strands_relaxed(
        bio_index[k].strand, raw_index[k].strand
    )]

    return {
        "sequence_identical": bio["sequence"] == raw["sequence"],
        "topology_identical": bio["topology"] == raw["topology"],
        "features_identical": not only_bio and not only_raw and not hard,
        "strand_relaxed": [
            {"feature": k[1], "biopython": bio_index[k].strand, "raw": raw_index[k].strand}
            for k in relaxed
            if k not in hard
        ],
        "only_in_biopython": only_bio,
        "only_in_raw_parser": only_raw,
        "strand_conflicts": [k[1] for k in hard],
    }


def _enrich_from_raw(features: list[Feature], raw_features: list[Feature]) -> None:
    """把内置解析器读到的信息补进 Biopython 的结果。

    只补两类**元数据**，不改变任何坐标：

    * ``Segment.raw`` —— 文件里的原始 1-based 区间字符串，便于人工核对；
    * ``strand`` —— 当文件本来没写方向（``"."``）而 Biopython 默认成 ``"+"`` 时，
      以文件的实际情况为准。
    """
    index = _index(raw_features)
    for feature in features:
        key = (
            feature.feature_type,
            feature.label,
            tuple((s.start, s.end) for s in feature.segments),
        )
        source = index.get(key)
        if source is None:
            continue
        if _strands_relaxed(feature.strand, source.strand):
            feature.strand = source.strand
            feature.segments = [
                Segment(
                    start=seg.start,
                    end=seg.end,
                    strand=source.strand,
                    wraps=seg.wraps,
                    raw=seg.raw,
                    origin=seg.origin,
                )
                for seg in feature.segments
            ]
        if len(source.segments) != len(feature.segments):
            continue
        feature.segments = [
            Segment(
                start=seg.start,
                end=seg.end,
                strand=seg.strand,
                wraps=seg.wraps,
                raw=src.raw,
                origin=seg.origin,
            )
            for seg, src in zip(feature.segments, source.segments)
        ]


# --------------------------------------------------------------------------- #
# 公共入口
# --------------------------------------------------------------------------- #
def load_snapgene(
    path: str | pathlib.Path,
    *,
    backend: str = "auto",
    scaffold_id: str | None = None,
) -> Scaffold:
    """读取 SnapGene ``.dna`` 文件。

    Parameters
    ----------
    backend:
        ``"auto"``（默认，优先 Biopython）、``"biopython"``、``"raw"``。
    scaffold_id:
        默认用文件名（不含扩展名）。
    """
    if backend not in ("auto", "biopython", "raw"):
        raise ValueError(f"backend 必须是 auto/biopython/raw，收到 {backend!r}")
    path = pathlib.Path(path)
    if not path.exists():
        raise FileNotFoundError(f"找不到文件：{path}")

    raw_digest = hashlib.sha256(path.read_bytes()).hexdigest()
    raw = _from_raw(path)  # 始终跑一遍：Notes / 引物序列 / 交叉校验都靠它

    warnings: list[str] = list(raw.get("raw_warnings", []))
    cross_check: dict[str, Any] = {}
    backend_used = "raw"

    bio: dict[str, Any] | None = None
    if backend in ("auto", "biopython"):
        try:
            bio = _from_biopython(path)
        except BackendUnavailable:
            if backend == "biopython":
                raise
            warnings.append(
                "Biopython 不可用，已退回内置解析器（parser_backend=raw）；"
                "引物仅有序列、没有回贴位点。"
            )
        except Exception as exc:  # Biopython 解析失败时也退回
            if backend == "biopython":
                raise
            warnings.append(f"Biopython 解析失败（{exc}），已退回内置解析器。")

    if bio is not None:
        cross_check = _cross_check(bio, raw)
        if not cross_check["sequence_identical"]:
            raise snapgene_raw.SnapGeneError(
                "两个解析器给出的序列不一致：这不该发生，请保留文件并上报。"
            )
        if not cross_check["topology_identical"]:
            warnings.append(
                f"两个解析器对拓扑判断不一致："
                f"biopython={bio['topology']} raw={raw['topology']}"
            )
        if not cross_check["features_identical"]:
            warnings.append(
                "两个解析器对特征集合判断不一致，详见 cross_check 字段（未做任何自动合并）。"
            )

        _enrich_from_raw(bio["features"], raw["features"])
        features = bio["features"]
        sequence = bio["sequence"]
        topology = bio["topology"]
        molecule_type = bio["molecule_type"]
        molecule_type_source = "biopython annotation"
        primers = [
            Primer(
                name=p.name,
                sequence=p.sequence,
                site=bio["primer_sites"].get(p.name, []),
                site_origin="derived",
            )
            for p in raw.get("primers", [])
        ]
        backend_used = "biopython+raw"
    else:
        features = raw["features"]
        sequence = raw["sequence"]
        topology = raw["topology"]
        molecule_type = raw["molecule_type"]
        molecule_type_source = "default（内置解析器不解码 molecule_type）"
        primers = list(raw.get("primers", []))
        backend_used = "raw"

    if topology not in ("circular", "linear"):
        raise snapgene_raw.SnapGeneError(f"无法识别拓扑：{topology!r}")

    scaffold = Scaffold(
        id=scaffold_id or _sanitize_id(path.stem),
        sequence=sequence,
        length=len(sequence),
        topology=topology,
        features=features,
        primers=primers,
        molecule_type=molecule_type,
        molecule_type_source=molecule_type_source,
        source=str(path),
        sha256=raw_digest,
        parser_backend=backend_used,
        notes=dict(raw.get("notes", {})),
        sequence_class=raw.get("sequence_class", ""),
        cross_check=cross_check,
        warnings=warnings,
    )
    scaffold.warnings.extend(check_sequence(scaffold.sequence))
    return scaffold


def load_fasta(
    path: str | pathlib.Path,
    *,
    record_id: str | None = None,
    topology: str = "linear",
) -> Scaffold:
    """读取 FASTA。多记录文件必须用 ``record_id`` 指定，否则报错（不猜）。"""
    path = pathlib.Path(path)
    if not path.exists():
        raise FileNotFoundError(f"找不到文件：{path}")

    records: list[tuple[str, str]] = []
    header: str | None = None
    chunks: list[str] = []
    for line in path.read_text(encoding="utf-8-sig").splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(chunks)))
            header = line[1:].strip()
            chunks = []
        else:
            if header is None:
                raise ValueError(f"{path} 不是 FASTA：第一行不是 '>' 开头")
            chunks.append(line)
    if header is not None:
        records.append((header, "".join(chunks)))

    if not records:
        raise ValueError(f"{path} 里没有序列记录")

    if record_id is None:
        if len(records) > 1:
            names = [h.split()[0] for h, _ in records]
            raise ValueError(
                f"{path} 含 {len(records)} 条记录 {names}；"
                "请用 record_id 指定要读哪一条（不自动取第一条）"
            )
        header, text = records[0]
    else:
        matched = [(h, t) for h, t in records if h.split()[0] == record_id]
        if not matched:
            raise ValueError(f"记录 {record_id!r} 不存在")
        if len(matched) > 1:
            raise ValueError(f"记录 {record_id!r} 出现 {len(matched)} 次，无法确定要读哪一条")
        header, text = matched[0]

    sequence = normalize_sequence(text)
    scaffold = Scaffold(
        id=record_id or (header.split()[0] if header else _sanitize_id(path.stem)),
        sequence=sequence,
        length=len(sequence),
        topology=topology,
        source=str(path),
        sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        parser_backend="fasta",
        warnings=[],
    )
    scaffold.warnings.extend(check_sequence(sequence))
    return scaffold


def circular_slice(sequence: str, start: int, length: int) -> str:
    """便捷包装：环状切片（供外部脚本使用）。"""
    return get_circular_subsequence(sequence, start, length, circular=True)
