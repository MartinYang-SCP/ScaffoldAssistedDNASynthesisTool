"""命令行入口。

    python -m ssdna.cli "minimalplasmid(1).dna" --outdir outputs
    python -m ssdna.cli --fasta template.fa --topology circular --outdir outputs
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import __version__
from .io import load_fasta, load_snapgene
from .report import export_all


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="ssdna",
        description="读取 scaffold 文件并导出 FASTA / 注释 CSV / JSON（不生成任何设计序列）",
    )
    parser.add_argument("path", nargs="?", help="输入文件（.dna 或 .fasta）")
    parser.add_argument("--fasta", help="以 FASTA 方式读取（覆盖自动判断）")
    parser.add_argument("--topology", choices=["circular", "linear"], default="linear",
                        help="FASTA 输入的拓扑（默认 linear，不猜）")
    parser.add_argument("--backend", choices=["auto", "biopython", "raw"], default="auto",
                        help=".dna 解析后端（默认 auto：优先 Biopython）")
    parser.add_argument("--id", dest="scaffold_id", help="scaffold id（默认用文件名）")
    parser.add_argument("--outdir", default="outputs", help="输出目录（默认 outputs/）")
    parser.add_argument("--no-sequence", action="store_true",
                        help="JSON 中省略完整序列")
    parser.add_argument("--version", action="version", version=f"ssdna {__version__}")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    if args.fasta:
        scaffold = load_fasta(args.fasta, record_id=args.scaffold_id, topology=args.topology)
    else:
        if not args.path:
            raise SystemExit("请给出输入文件，或用 --fasta 指定")
        source = Path(args.path)
        if source.is_dir():
            raise SystemExit(f"{args.path} 是目录，请给出 .dna 或 .fasta 文件")
        if source.suffix.lower() in {".fa", ".fasta", ".fna"}:
            scaffold = load_fasta(
                args.path, record_id=args.scaffold_id, topology=args.topology
            )
        else:
            scaffold = load_snapgene(
                args.path, backend=args.backend, scaffold_id=args.scaffold_id
            )

    paths = export_all(scaffold, args.outdir)

    print(f"[ssdna {__version__}] {scaffold.summary()}")
    print(f"  molecule_type : {scaffold.molecule_type} ({scaffold.molecule_type_source})")
    print(f"  file sha256   : {scaffold.sha256}")
    print(f"  seq  sha256   : {scaffold.sequence_sha256}")
    if scaffold.notes:
        keys = ("Type", "ConfirmedExperimentally", "SequenceClass", "LastModified")
        brief = {k: scaffold.notes[k] for k in keys if k in scaffold.notes}
        print(f"  notes         : {brief}")
    if scaffold.cross_check:
        cc = scaffold.cross_check
        print(
            "  cross-check   : "
            f"sequence={'OK' if cc['sequence_identical'] else 'MISMATCH'}, "
            f"topology={'OK' if cc['topology_identical'] else 'MISMATCH'}, "
            f"features={'OK' if cc['features_identical'] else 'MISMATCH'}"
        )
    for warning in scaffold.warnings:
        print(f"  ! {warning}")
    for name, path in paths.items():
        print(f"  -> {name:12s} {path}")

    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
