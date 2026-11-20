"""热力学体检（**只做检查，不参与选序列**）。

为什么不能"据此设计"：

* 这套结构里 staple 的序列是**被唯一确定**的——13 nt 必须是 scaffold 那段的互补，
  26 nt 必须是目标那段的互补。没有"换一段序列把 Tm 调到 60 °C"的自由度；
* 唯一的自由度只有三个，而且都不是 Tm 直接决定的：
  ① scaffold 平铺的**相位**（从覆盖区哪一端开始铺 13 nt 块）；
  ② **nick 位置**（目标分段）；
  ③ 退火**程序**（升降温曲线、保持温度）。

所以这里算 Tm 的用途是：

1. 找出 13 nt 臂里 Tm 的**离群者**（GC 极高/极低），提示人工看一眼；
2. 看两组 staple 的 Tm **分布是否均匀**（差异大 ⇒ 组装不同步）；
3. 给退火温度一个下限参考（通常取最低 13 nt 臂 Tm 附近再往上）。

计算用 SantaLucia 最近邻参数（Biopython ``Tm_NN``），未做 Mg²⁺ 的完整修正，
所以数值只用于**同一批序列内部横向比较**，不要拿绝对值当结论。
"""

from __future__ import annotations

import statistics

__all__ = ["tm_nn", "annotate_staples", "summarize", "available"]

try:  # pragma: no cover - 取决于环境
    from Bio.SeqUtils import MeltingTemp as _mt

    _HAS_BIO = True
except ImportError:  # pragma: no cover
    _mt = None
    _HAS_BIO = False


def available() -> bool:
    return _HAS_BIO


def tm_nn(
    sequence: str,
    *,
    na_mM: float = 50.0,
    mg_mM: float = 0.0,
    dnac_nM: float = 250.0,
) -> float | None:
    """最近邻 Tm（°C）。没有 Biopython 时返回 None，不猜。"""
    if not _HAS_BIO or not sequence:
        return None
    try:  # 短序列在某些盐条件下会抛错，交给调用方处理
        return float(
            _mt.Tm_NN(
                sequence,
                Na=float(na_mM),
                Mg=float(mg_mM),
                dnac1=float(dnac_nM),
                dnac2=float(dnac_nM),
            )
        )
    except Exception:  # noqa: BLE001 - 参数越界等
        return None


def annotate_staples(
    staples: dict[str, dict],
    *,
    na_mM: float,
    mg_mM: float,
    dnac_nM: float,
) -> None:
    """就地给每条 staple 加上 tm_scaffold_domain / tm_target_domain。"""
    for row in staples.values():
        row["tm_scaffold_domain"] = tm_nn(
            row["scaffold_domain"], na_mM=na_mM, mg_mM=mg_mM, dnac_nM=dnac_nM
        )
        row["tm_target_domain"] = tm_nn(
            row["target_domain"], na_mM=na_mM, mg_mM=mg_mM, dnac_nM=dnac_nM
        )


def summarize(staples: dict[str, dict]) -> dict:
    """给出两组 domain 的 Tm 分布。数值只用于同批横向比较。"""
    if not _HAS_BIO:
        return {
            "available": False,
            "note": "未安装 Biopython，未计算 Tm（不是 0，是没算）",
        }

    out: dict = {"available": True, "note": "SantaLucia 最近邻；只做同批横向比较"}
    for key in ("tm_scaffold_domain", "tm_target_domain"):
        values = [row[key] for row in staples.values() if row.get(key) is not None]
        if not values:
            continue
        out[key] = {
            "n": len(values),
            "min": round(min(values), 2),
            "max": round(max(values), 2),
            "mean": round(statistics.fmean(values), 2),
            "stdev": round(statistics.pstdev(values), 2) if len(values) > 1 else 0.0,
            "outliers": [
                row_id
                for row_id, row in staples.items()
                if row.get(key) is not None
                and abs(row[key] - statistics.fmean(values)) > 2 * (
                    statistics.pstdev(values) or 1e9
                )
            ],
        }
    return out
