"""向后兼容入口：解析器已迁到 :mod:`ssdna.snapgene_raw`。

保留这个文件是为了让既有脚本（``python -m tools.snapgene_reader ...``）与
``tests/test_snapgene_reader.py`` 继续可用。
"""

from __future__ import annotations

from ssdna.snapgene_raw import (  # noqa: F401
    BLOCK_NAMES,
    Feature,
    Primer,
    Segment,
    SnapGeneError,
    SnapGeneFile,
    _cli,
    parse_feature_ranges,
    read_snapgene,
    validate,
)

__all__ = [
    "Segment",
    "Feature",
    "Primer",
    "SnapGeneFile",
    "SnapGeneError",
    "read_snapgene",
    "parse_feature_ranges",
    "validate",
]


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_cli())
