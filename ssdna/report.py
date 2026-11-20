"""导出：FASTA / CSV / JSON。

CSV 用 ``utf-8-sig`` 写，Excel 直接双击不会乱码（注释里有中文与特殊符号）。
坐标列一律 0-based 半开区间，并在每一份输出的头部/字段里重复声明。
"""

from __future__ import annotations

import csv
import json
import pathlib
from typing import Any

from dataclasses import replace

from .models import Feature, Scaffold, Segment

__all__ = [
    "build_info",
    "write_fasta",
    "write_features_csv",
    "write_info_json",
    "export_all",
    "FEATURE_COLUMNS",
]

COORDINATE_NOTE = "坐标均为 0-based 左闭右开 [start, end)；跨原点时 start > end 且 wraps=True"

FEATURE_COLUMNS = [
    "feature_type",
    "label",
    "start",
    "end",
    "strand",
    "note",
    "kind",
    "segment_index",
    "segment_count",
    "length",
    "wraps",
    "origin",
    "raw_range",
    "sequence",
]


def _segment_row(
    scaffold: Scaffold,
    feature: Feature,
    segment: Segment | None,
    index: int,
    count: int,
) -> list[Any]:
    if segment is None:
        start = end = length = ""
        wraps = ""
        raw = ""
        origin = feature.origin
    else:
        start = segment.start
        end = segment.end
        length = segment.length_on(scaffold.length)
        wraps = "yes" if segment.wraps else ""
        raw = segment.raw
        origin = segment.origin
    return [
        feature.feature_type,
        feature.label,
        start,
        end,
        feature.strand,
        feature.note,
        feature.kind,
        index,
        count,
        length,
        wraps,
        origin,
        raw,
        "",  # sequence 列只对 primer 行有值
    ]


def write_features_csv(scaffold: Scaffold, path: str | pathlib.Path) -> pathlib.Path:
    """导出注释表（含 primer 行，用 kind 列区分）。"""
    path = pathlib.Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    rows: list[list[Any]] = []
    for feature in scaffold.features:
        count = len(feature.segments)
        for index, segment in enumerate(feature.segments, start=1):
            rows.append(_segment_row(scaffold, feature, segment, index, count))

    for primer in scaffold.primers:
        site = primer.site
        if not site:
            fake = Feature(
                feature_type="primer_bind",
                label=primer.name,
                strand=".",
                segments=[],
                qualifiers={"note": "结合位点未回贴（backend=raw）"},
                kind="primer",
                origin="file",
            )
            row = _segment_row(scaffold, fake, None, 1, 1)
            row[-1] = primer.sequence
            rows.append(row)
            continue
        for index, segment in enumerate(site, start=1):
            fake = Feature(
                feature_type="primer_bind",
                label=primer.name,
                strand=segment.strand,
                segments=[segment],
                qualifiers={"note": "结合位点由引物序列回贴推出（derived）"},
                kind="primer",
                origin="derived",
            )
            row = _segment_row(
                scaffold, fake, replace(segment, origin="derived"), index, len(site)
            )
            row[-1] = primer.sequence
            rows.append(row)

    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(FEATURE_COLUMNS)
        writer.writerows(rows)
    return path


def write_fasta(
    scaffold: Scaffold,
    path: str | pathlib.Path,
    *,
    width: int = 70,
) -> pathlib.Path:
    path = pathlib.Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    header = (
        f">{scaffold.id} length={scaffold.length} topology={scaffold.topology} "
        f"source={pathlib.Path(scaffold.source).name} sha256={scaffold.sequence_sha256[:12]}"
    )
    lines = [header]
    for i in range(0, scaffold.length, width):
        lines.append(scaffold.sequence[i : i + width])
    path.write_text("\n".join(lines) + "\n", encoding="ascii")
    return path


def build_info(scaffold: Scaffold, *, include_sequence: bool = True) -> dict[str, Any]:
    def segment_dict(seg: Segment) -> dict[str, Any]:
        return {
            "start": seg.start,
            "end": seg.end,
            "length": seg.length_on(scaffold.length),
            "wraps": seg.wraps,
            "raw_range": seg.raw,
            "origin": seg.origin,
        }

    info: dict[str, Any] = {
        "schema": "scaffold_staple_design/scaffold_info/v1",
        "id": scaffold.id,
        "source": scaffold.source,
        "file_sha256": scaffold.sha256,
        "sequence_sha256": scaffold.sequence_sha256,
        "parser_backend": scaffold.parser_backend,
        "length": scaffold.length,
        "topology": scaffold.topology,
        "molecule_type": scaffold.molecule_type,
        "molecule_type_source": scaffold.molecule_type_source,
        "gc_content": round(scaffold.gc_content, 6),
        "sequence_class": scaffold.sequence_class,
        "features": [
            {
                "feature_type": f.feature_type,
                "label": f.label,
                "strand": f.strand,
                "kind": f.kind,
                "origin": f.origin,
                "start": f.start,
                "end": f.end,
                "total_length": f.total_length(scaffold.length),
                "segments": [segment_dict(s) for s in f.segments],
                "qualifiers": f.qualifiers,
            }
            for f in scaffold.features
        ],
        "primers": [
            {
                "name": p.name,
                "sequence": p.sequence,
                "site_origin": p.site_origin,
                "site": [segment_dict(s) for s in p.site],
            }
            for p in scaffold.primers
        ],
        "notes": scaffold.notes,
        "cross_check": scaffold.cross_check,
        "warnings": scaffold.warnings,
        "coordinate_convention": COORDINATE_NOTE,
        "scaffold_region_interpretation": {
            "scaffold_region": None,
            "scaffold_strand": None,
            "note": (
                "本模块只报告文件里写了什么，不判断哪一段才是实验室实际使用的"
                "环状单链 scaffold。该决定必须由用户在 params 中显式给出。"
            ),
        },
    }
    if include_sequence:
        info["sequence"] = scaffold.sequence
    return info


def write_info_json(
    scaffold: Scaffold,
    path: str | pathlib.Path,
    *,
    include_sequence: bool = True,
) -> pathlib.Path:
    path = pathlib.Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(
        build_info(scaffold, include_sequence=include_sequence),
        ensure_ascii=False,
        indent=2,
    )
    path.write_text(payload, encoding="utf-8")
    return path


def export_all(
    scaffold: Scaffold,
    outdir: str | pathlib.Path,
    *,
    stem: str = "scaffold",
) -> dict[str, pathlib.Path]:
    outdir = pathlib.Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    return {
        "fasta": write_fasta(scaffold, outdir / f"{stem}.fasta"),
        "features_csv": write_features_csv(scaffold, outdir / f"{stem}_features.csv"),
        "info_json": write_info_json(scaffold, outdir / f"{stem}_info.json"),
    }
