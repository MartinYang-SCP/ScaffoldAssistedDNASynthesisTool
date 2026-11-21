"""参数完整性校验：把"不得自行猜测"变成一条可执行的规则。

权威配置只有一份：``config/design_config.yaml``（经 ``ssdna.design.load_config`` 读取）。
历史上还有 ``params/project.params.json`` 与 ``params/design_config.yaml`` 两份，
已删除——多份来源只会"改一处忘另一处"。

    python tools/check_params.py [config/design_config.yaml]

规则：
  * ``find_missing``：值为 ``null`` 的键 = 未确定，是阻塞项；
  * ``find_invalid``：值非空但类型/取值域不对（数字写成字符串、枚举拼错、切点越界），
    同样是阻塞项；
  * 有任一阻塞项时退出码 2，逐条列出，不代入默认值。
"""

from __future__ import annotations

import pathlib
import sys
from typing import Any

ROOT = pathlib.Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
DEFAULT_CONFIG = ROOT / "config" / "design_config.yaml"


def blocking(config: dict[str, Any]) -> list[str]:
    """返回所有阻塞项路径（缺失 + 非法），已排序。"""
    from ssdna.design import find_invalid, find_missing

    return sorted(set(find_missing(config)) | set(find_invalid(config)))


def main(argv: list[str] | None = None) -> int:
    from ssdna.design import find_invalid, find_missing, load_config

    argv = list(sys.argv[1:] if argv is None else argv)
    path = pathlib.Path(argv[0]) if argv else DEFAULT_CONFIG
    config = load_config(path)

    missing = find_missing(config)
    invalid = find_invalid(config)

    if not missing and not invalid:
        print("参数齐备，可以进入序列生成阶段。")
        return 0

    print(f"阻塞项 {len(missing) + len(invalid)} 个 —— 不生成任何序列：")
    if missing:
        print("\n缺失（值为 null）：")
        for key in missing:
            print(f"  - {key}")
    if invalid:
        print("\n非法（类型/取值不对）：")
        for key, why in invalid.items():
            print(f"  - {key}: {why}")
    print("\n这些值必须由用户给出或修正；工具不会代入默认值。")
    return 2


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
