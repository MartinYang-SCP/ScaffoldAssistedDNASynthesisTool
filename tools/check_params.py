"""参数完整性校验：把"不得自行猜测"变成一条可执行的规则。

    python tools/check_params.py [params/project.params.json]

规则：
  * 任何值为 ``null`` 的键都是"未确定"；
  * 未确定的键只要不在 ``output`` 段（输出格式允许有默认值），就视为阻塞；
  * 有阻塞项时退出码为 2，并逐条列出 —— 绝不代入默认值、绝不猜。

以 ``_`` 开头的键是注释/选项说明，不参与校验。
"""

from __future__ import annotations

import json
import pathlib
import sys
from typing import Any

ROOT = pathlib.Path(__file__).resolve().parents[1]
DEFAULT_PARAMS = ROOT / "params" / "project.params.json"

# 这些段允许留空（有合理默认值，不影响结构正确性）
NON_BLOCKING_SECTIONS = {"output"}


def iter_unresolved(node: Any, path: str = "") -> list[tuple[str, str]]:
    """返回 [(路径, 所在段)]，只收集值为 None 的键。"""
    found: list[tuple[str, str]] = []
    if isinstance(node, dict):
        for key, value in node.items():
            if key.startswith("_"):
                continue
            child = f"{path}.{key}" if path else key
            if value is None:
                found.append((child, path.split(".")[0] if path else key))
            else:
                found.extend(iter_unresolved(value, child))
    return found


def blocking(params: dict[str, Any]) -> list[tuple[str, str]]:
    return [
        item
        for item in iter_unresolved(params)
        if item[1] not in NON_BLOCKING_SECTIONS
    ]


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    path = pathlib.Path(argv[0]) if argv else DEFAULT_PARAMS
    params = json.loads(path.read_text(encoding="utf-8"))

    pending = blocking(params)
    if not pending:
        print("参数齐备，可以进入序列生成阶段。")
        return 0

    print(f"阻塞项 {len(pending)} 个 —— 不生成任何序列：")
    current = None
    for key, section in pending:
        if section != current:
            print(f"\n[{section}]")
            current = section
        print(f"  - {key}")
    print("\n这些值必须由用户给出；工具不会代入默认值。")
    return 2


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
