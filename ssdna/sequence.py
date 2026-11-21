"""序列级工具：反向互补、规范化、环状切片、检查。

这里刻意把"规范化"和"检查"分开：

* ``normalize_sequence`` 只做**不改变碱基身份**的处理（去空白/换行/FASTA 行号数字、
  统一大写）。遇到任何无法归一的字符直接抛错，不猜、不替换。
* ``check_sequence`` 只报告，不修改。

反向互补保留 IUPAC 简并码（N/R/Y/...），因为 scaffold 的注释里可能出现它们；
输出仍是 5′→3′，不生成"反向链"副本。
"""

from __future__ import annotations

__all__ = [
    "reverse_complement",
    "normalize_sequence",
    "check_sequence",
    "get_circular_subsequence",
    "gc_content",
]

# IUPAC 互补表（含小写）：A↔T, C↔G, U→A, R↔Y, S↔S, W↔W, K↔M, B↔V, D↔H, N↔N
_COMPLEMENT = str.maketrans(
    "ACGTURYSWKMBDHVNacgturyswkmbdhvn",
    "TGCAAYRSWMKVHDBNtgcaayrswmkvhdbn",
)
#: 允许出现的字符（大写形式）
ALLOWED_BASES = frozenset("ACGTUNRYKMSWBDHV")
#: 可以用于常规合成的字符
CANONICAL_BASES = frozenset("ACGT")


def reverse_complement(sequence: str) -> str:
    """返回反向互补序列（仍按 5′→3′ 书写）。

    >>> reverse_complement("ATGC")
    'GCAT'
    """
    if not isinstance(sequence, str):
        raise TypeError("sequence 必须是 str")
    bad = sorted(set(sequence.upper()) - ALLOWED_BASES)
    if bad:
        raise ValueError(
            f"序列含无法互补的字符 {bad}；只接受 ACGT 及 IUPAC 简并码"
        )
    return sequence.translate(_COMPLEMENT)[::-1]


def normalize_sequence(sequence: str) -> str:
    """去空白、统一大写。

    只改变**书写形式**，不改变任何一个碱基。数字、非法字母都报错并给位置，
    不做去除、不替换成 N——"静默删掉数字"会掩盖用户输入里的错误。
    """
    if not isinstance(sequence, str):
        raise TypeError("sequence 必须是 str")

    kept: list[str] = []
    for char in sequence:
        if char.isspace():
            continue
        kept.append(char)

    if not kept:
        raise ValueError("序列为空（去除空白后没有剩余字符）")

    text = "".join(kept).upper()
    bad = sorted(set(text) - ALLOWED_BASES)
    if bad:
        positions = [i for i, ch in enumerate(text) if ch in bad][:10]
        hint = "（含数字？FASTA 序列不含行号，数字视为非法）" if any(
            ch.isdigit() for ch in bad
        ) else ""
        raise ValueError(
            f"序列含非法碱基字符 {bad}{hint}；示例位置（归一化后 0-based）：{positions}。"
            "不做替换、不丢弃，请先修正输入。"
        )
    return text


def check_sequence(sequence: str) -> list[str]:
    """非致命提示清单；不修改任何数据。"""
    if not isinstance(sequence, str):
        raise TypeError("sequence 必须是 str")
    warnings: list[str] = []
    if not sequence:
        warnings.append("序列为空")
        return warnings

    upper = sequence.upper()
    bad = sorted(set(upper) - ALLOWED_BASES)
    if bad:
        warnings.append(f"含非法碱基字符：{bad}")

    degenerate = sorted(set(upper) & (ALLOWED_BASES - CANONICAL_BASES))
    if degenerate:
        warnings.append(f"含简并碱基 {degenerate}（可合成性/配对需单独处理）")

    return warnings


def get_circular_subsequence(
    sequence: str,
    start: int,
    length: int,
    *,
    circular: bool = True,
) -> str:
    """取 ``start`` 起 ``length`` 个碱基；``circular=True`` 时跨原点回绕。

    例如长度 10 的环、``start=8, length=4`` 返回 ``seq[8:] + seq[:2]``。
    不静默修正：起点越界、长度超过环长都直接报错。
    """
    if not isinstance(sequence, str):
        raise TypeError("sequence 必须是 str")
    if not sequence:
        raise ValueError("序列为空")
    if isinstance(start, bool) or not isinstance(start, int):
        raise TypeError(f"start 必须是整数，收到 {start!r}")
    if isinstance(length, bool) or not isinstance(length, int):
        raise TypeError(f"length 必须是整数，收到 {length!r}")
    if length < 0:
        raise ValueError(f"length 不能为负：{length}")

    n = len(sequence)
    if not (0 <= start < n):
        raise IndexError(f"start={start} 越界（序列长度 {n}，要求 0 <= start < {n}）")

    if length > n:
        raise IndexError(f"length={length} 超过序列长度 {n}")

    end = start + length
    if end <= n:
        return sequence[start:end]
    if not circular:
        raise IndexError(
            f"区间 [{start}, {end}) 超出线性序列长度 {n}；如需跨原点请设 circular=True"
        )
    return sequence[start:] + sequence[: end - n]


def gc_content(sequence: str) -> float:
    if not sequence:
        return 0.0
    upper = sequence.upper()
    return sum(base in "GC" for base in upper) / len(upper)
