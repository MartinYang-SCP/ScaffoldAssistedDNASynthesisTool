"""生成阶段 0 的三张说明图（不涉及任何序列设计）。

    python tools/make_figures.py

输出：
    docs/figures/fig1_structure_redraw.png   示意图的结构重绘 + 方向标注 + 冲突点
    docs/figures/fig2_scaffold_circular_map.png  scaffold 文件的环状图谱
    docs/figures/fig3_staple_geometry.png    39-mer 跨两条链的两种几何
"""

from __future__ import annotations

import json
import math
import pathlib

from PIL import Image, ImageDraw, ImageFont

ROOT = pathlib.Path(__file__).resolve().parents[1]
OUT = ROOT / "docs" / "figures"
FONT = r"C:\Windows\Fonts\msyh.ttc"
FONT_BOLD = r"C:\Windows\Fonts\msyhbd.ttc"

BLACK = (26, 26, 26)
RED = (192, 57, 43)
RED_LIGHT = (231, 126, 113)
BLUE = (47, 111, 191)
BLUE_LIGHT = (150, 180, 220)
GRAY = (120, 120, 120)
LIGHT = (238, 238, 238)
GREEN = (39, 139, 86)
GOLD = (200, 150, 20)


def font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont:
    return ImageFont.truetype(FONT_BOLD if bold else FONT, size)


def text(d: ImageDraw.ImageDraw, xy, s, size=22, fill=BLACK, bold=False, anchor="la"):
    d.text(xy, s, font=font(size, bold), fill=fill, anchor=anchor)


def arrow(d: ImageDraw.ImageDraw, p0, p1, color, width=5, head=18, dashed=False):
    if dashed:
        _dashed_line(d, p0, p1, color, width, dash=14, gap=10)
        return
    d.line([p0, p1], fill=color, width=width)
    ang = math.atan2(p1[1] - p0[1], p1[0] - p0[0])
    a1 = (p1[0] - head * math.cos(ang - 0.42), p1[1] - head * math.sin(ang - 0.42))
    a2 = (p1[0] - head * math.cos(ang + 0.42), p1[1] - head * math.sin(ang + 0.42))
    d.polygon([p1, a1, a2], fill=color)


def _dashed_line(d, p0, p1, color, width, dash=14, gap=10):
    x0, y0 = p0
    x1, y1 = p1
    total = math.hypot(x1 - x0, y1 - y0)
    if total == 0:
        return
    ux, uy = (x1 - x0) / total, (y1 - y0) / total
    pos = 0.0
    while pos < total:
        end = min(pos + dash, total)
        d.line(
            [(x0 + ux * pos, y0 + uy * pos), (x0 + ux * end, y0 + uy * end)],
            fill=color,
            width=width,
        )
        pos = end + gap


def dashed_arc(d, box, start, end, color, width=4):
    """PIL 的 arc 不支持虚线，这里按角度分段画。"""
    for a in range(int(start), int(end), 12):
        d.arc(box, a, min(a + 7, end), fill=color, width=width)


def brace(d, x, y0, y1, color, width=4, tick=18):
    d.line([(x, y0), (x, y1)], fill=color, width=width)
    d.line([(x, y0), (x + tick, y0)], fill=color, width=width)
    d.line([(x, y1), (x + tick, y1)], fill=color, width=width)


# --------------------------------------------------------------------------- #
# 图 1：结构重绘
# --------------------------------------------------------------------------- #
def figure1() -> pathlib.Path:
    W, H = 2000, 1580
    img = Image.new("RGB", (W, H), "white")
    d = ImageDraw.Draw(img)

    text(d, (60, 36), "图 1｜示意图的结构重绘：四类链、5′→3′ 方向与四处冲突点", 40, BLACK, True)
    text(
        d,
        (62, 92),
        "按你上传的原图逐元素重画，仅把线条画直、把方向标注写清楚；带「？」的问题见文档 §8。",
        22,
        GRAY,
    )

    X0, X1 = 200, 1800
    yA, yA1, yA2 = 300, 400, 500
    yS = 620
    yB1, yB2, yB = 740, 840, 960

    # --- A 组目标（黑）---
    arrow(d, (430, yA), (1010, yA), BLACK, 6)
    text(d, (400, yA - 46), "5′", 26, BLACK, True)
    text(d, (1000, yA - 46), "3′", 26, BLACK, True)
    arrow(d, (1050, yA), (1560, yA), BLACK, 6)
    text(d, (1042, yA - 46), "5′", 26, BLACK, True)
    text(d, (1570, yA - 46), "3′", 26, BLACK, True)
    d.line([(1030, yA - 26), (1030, yA + 26)], fill=BLACK, width=4)
    text(d, (1030, yA + 34), "nick", 16, GRAY, anchor="ma")
    text(d, (430, yA - 90), "A 组目标 oligo（黑色，可多条，带 nick）", 24, BLACK, True)

    # --- staple A（红，U 形）---
    arrow(d, (1520, yA1), (360, yA1), RED, 6)
    d.line([(360, yA1), (360, yA2)], fill=RED, width=6)
    arrow(d, (360, yA2), (1500, yA2), RED, 6)
    text(d, (1532, yA1 - 20), "5′", 26, RED, True)
    text(d, (1512, yA2 - 20), "3′", 26, RED, True)
    text(d, (520, yA1 - 44), "26 nt：与 A 组目标配对", 22, RED, True)
    text(d, (520, yA2 + 22), "13 nt：与 scaffold 配对", 22, RED, True)
    brace(d, 300, yA1 - 6, yA2 + 6, RED)
    text(d, (150, (yA1 + yA2) // 2), "折返", 20, RED, anchor="rm")

    # --- scaffold（蓝）---
    arrow(d, (X0, yS), (X1, yS), BLUE, 7)
    text(d, (X0 - 34, yS - 22), "5′", 30, BLUE, True, anchor="rm")
    text(d, (X1 + 12, yS - 22), "3′", 30, BLUE, True)
    text(
        d,
        (X0, yS + 26),
        "scaffold（蓝色）：环状单链，约 1.2 kb；此处把它线性展开，左右两端其实是同一个环的断点",
        24,
        BLUE,
        True,
    )
    brace(d, 1560, yS - 60, yS - 12, GRAY)
    text(d, (1566, yS - 54), "冗余 / 未覆盖区（长度待定）", 20, GRAY)

    # --- staple B（红，U 形）---
    arrow(d, (1120, yB1), (1480, yB1), RED, 6)
    d.line([(1480, yB1), (1480, yB2)], fill=RED, width=6)
    arrow(d, (1480, yB2), (380, yB2), RED, 6)
    text(d, (1104, yB1 - 20), "5′", 26, RED, True)
    text(d, (490, yB2 + 22), "3′", 26, RED, True)
    text(d, (1100, yB1 - 44), "13 nt：与 scaffold 配对", 22, RED, True)
    text(d, (480, yB2 + 50), "26 nt：与 B 组目标配对", 22, RED, True)
    brace(d, 1500, yB1 - 6, yB2 + 6, RED)

    # --- B 组目标（黑）---
    arrow(d, (420, yB), (1010, yB), BLACK, 6)
    text(d, (390, yB + 22), "5′", 26, BLACK, True)
    text(d, (1000, yB + 22), "3′", 26, BLACK, True)
    arrow(d, (1050, yB), (1560, yB), BLACK, 6)
    text(d, (1042, yB + 22), "5′", 26, BLACK, True)
    text(d, (1570, yB + 22), "3′", 26, BLACK, True)
    d.line([(1030, yB - 26), (1030, yB + 26)], fill=BLACK, width=4)
    text(d, (1030, yB + 52), "nick", 16, GRAY, anchor="ma")
    text(d, (420, yB + 62), "B 组目标 oligo（黑色，可多条，带 nick）", 24, BLACK, True)

    # --- 两个交界（虚线）---
    dashed_arc(d, (60, 200, 470, 1060), 90, 270, GRAY, 4)
    text(d, (86, 176), "J1？", 26, GRAY, True)
    text(d, (86, 212), "首尾连接示意：", 18, GRAY)
    text(d, (86, 236), "虚线 = 无碱基、无长度", 18, GRAY)
    dashed_arc(d, (1530, 200, 1940, 1060), 270, 450, GRAY, 4)
    text(d, (1790, 176), "J2？", 26, GRAY, True)

    # --- 冲突标记 ---
    marks = [
        ((330, yA2 - 40), "①"),
        ((980, yB1 - 34), "②"),
        ((960, 530), "③"),
        ((330, 572), "④"),
    ]
    for (x, y), tag in marks:
        d.ellipse([x - 20, y - 20, x + 20, y + 20], fill="white", outline=GOLD, width=4)
        text(d, (x, y), tag, 22, GOLD, True, anchor="mm")

    text(d, (992, 502), "原图在此附近标了 266 bp / 126 bp；", 20, (150, 60, 40))
    text(d, (992, 530), "两者都不可能是 staple 长度（39 nt）。", 20, (150, 60, 40))

    # --- 图例 / 说明 ---
    top = 1180
    d.line([(60, top - 24), (W - 60, top - 24)], fill=LIGHT, width=3)
    text(d, (60, top), "图例与冲突点", 30, BLACK, True)

    legend = [
        (BLACK, "黑 = 目标 oligo（target oligo / boundary-spanning target oligo）"),
        (RED, "红 = 普通 staple，39 nt = 26 nt（目标臂）+ 13 nt（scaffold 臂），同一条链"),
        (BLUE, "蓝 = scaffold，环状单链，5′→3′ 按选定链别标注"),
        (GRAY, "虚线 = 首尾连接示意，只表示拓扑相邻，不代表碱基、不代表共价键"),
    ]
    for i, (color, label) in enumerate(legend):
        y = top + 52 + i * 40
        d.line([(70, y + 12), (130, y + 12)], fill=color, width=6)
        text(d, (144, y), label, 24, BLACK)

    conflicts = [
        "① 图上 staple 的 13 nt 臂（下臂，5′→3′ 由左向右）与 scaffold（5′→3′ 由左向右）同向，无法互补配对。",
        "② 下侧 staple 同理：13 nt 臂（5′→3′ 由左向右）与 scaffold 同向；且 26 nt 臂与 B 组目标也同向。",
        "③ 图中 266 bp / 126 bp 不可能是 staple 长度（39 nt），只能是目标区段的长度；其归属待确认。",
        "④ 左侧折返标为“无碱基”，真实 DNA 链不存在这种片段，说明它只是示意图记号。",
    ]
    for i, line in enumerate(conflicts):
        text(d, (60, top + 230 + i * 36), line, 22, (150, 60, 40))

    out = OUT / "fig1_structure_redraw.png"
    img.save(out)
    return out


# --------------------------------------------------------------------------- #
# 图 2：scaffold 环状图谱
# --------------------------------------------------------------------------- #
def figure2() -> pathlib.Path:
    data = json.loads(
        (ROOT / "outputs" / "scaffold_info.json").read_text(encoding="utf-8")
    )
    length = data["length"]
    W, H = 1400, 1400
    img = Image.new("RGB", (W, H), "white")
    d = ImageDraw.Draw(img)

    cx, cy, R = 700, 640, 400
    text(d, (60, 40), "图 2｜minimalplasmid(1).dna 环状图谱（0-based，半开区间）", 38, BLACK, True)
    text(
        d,
        (62, 96),
        f"长度 {length} bp ｜ 拓扑 circular ｜ GC "
        f"{100 * sum(c in 'GC' for c in data['sequence']) / length:.1f}% ｜ "
        f"seq sha256 {data['sequence_sha256'][:12]}…",
        22,
        GRAY,
    )

    d.ellipse([cx - R, cy - R, cx + R, cy + R], outline=LIGHT, width=3)

    # 刻度
    for pos in range(0, length, 100):
        ang = -90 + 360 * pos / length
        rad = math.radians(ang)
        outer = (cx + (R + 16) * math.cos(rad), cy + (R + 16) * math.sin(rad))
        inner = (cx + R * math.cos(rad), cy + R * math.sin(rad))
        d.line([inner, outer], fill=(200, 200, 200), width=2)
        if pos % 200 == 0:
            lp = (cx + (R + 44) * math.cos(rad), cy + (R + 44) * math.sin(rad))
            text(d, lp, f"{pos}", 18, GRAY, anchor="mm")

    # 特征弧
    def arc_for(start, end, color, width, offset=0):
        a0 = -90 + 360 * start / length
        a1 = -90 + 360 * end / length
        box = [cx - R - offset, cy - R - offset, cx + R + offset, cy + R + offset]
        d.arc(box, a0, a1, fill=color, width=width)

    palette = {"ori": GOLD, "BleoR": GREEN, "MCS": BLUE, "cexuf1": GRAY}
    arcs = []
    for feature in data["features"]:
        for seg in feature["segments"]:
            starts = [seg["start"]]
            ends = [seg["end"]]
            if seg["wraps"]:
                ends = [length]
                starts = [seg["start"]]
            arcs.append((feature["label"], starts[0], ends[0], feature["strand"]))

    for name, s, e, strand in arcs:
        color = palette.get(name, GRAY)
        offset = {"ori": -6, "BleoR": -14, "MCS": -22, "cexuf1": -30}.get(name, -6)
        arc_for(s, e, color, 10, offset)

    # 0 / length 断点
    d.line([(cx, cy - R - 34), (cx, cy - R + 40)], fill=RED, width=3)
    text(
        d,
        (cx + 14, cy - R - 54),
        "编号原点 0 ≡ 1200（MCS 在此被拆成两段）",
        20,
        RED,
        True,
    )

    # 图注
    items = [
        ("MCS  0–37 与 1181–1200（跨原点，共 56 bp）", BLUE),
        ("ori  78–667（589 bp，反向链）", GOLD),
        ("BleoR  691–1066（375 bp，反向链，CDS）", GREEN),
        ("cexuf1  1146–1166（20 bp）", GRAY),
    ]
    y = 1112
    text(d, (60, y - 44), "注释（坐标为 0-based 半开区间）", 24, BLACK, True)
    for label, color in items:
        d.line([(70, y + 12), (130, y + 12)], fill=color, width=8)
        text(d, (146, y), label, 22, BLACK)
        y += 38
    text(
        d,
        (60, y + 8),
        "注意：本图只陈述文件里有什么。哪一段才是实验室实际提供的环状单链 scaffold，文件里没有信息，需你确认。",
        22,
        (150, 60, 40),
    )

    out = OUT / "fig2_scaffold_circular_map.png"
    img.save(out)
    return out


# --------------------------------------------------------------------------- #
# 图 3：39-mer 的两种几何
# --------------------------------------------------------------------------- #
def figure3() -> pathlib.Path:
    W, H = 2000, 1000
    img = Image.new("RGB", (W, H), "white")
    d = ImageDraw.Draw(img)
    text(d, (60, 40), "图 3｜一条 39-mer 要同时配对两条链，只有一种几何是无 crossover 的", 38, BLACK, True)
    text(
        d,
        (62, 96),
        "staple = 5′-[26 nt 目标臂]-[13 nt scaffold 臂]-3′（端序待确认）；两个 domain 共用一条连续磷酸二酯骨架。",
        22,
        GRAY,
    )

    # ---- 左：反向平行，U 形折返 ----
    box_l = (60, 170, 980, 720)
    d.rounded_rectangle(box_l, 18, outline=GREEN, width=3)
    text(d, (90, 190), "情形 A：两条伙伴链反向平行 → U 形折返，可行", 26, GREEN, True)

    x0, x1 = 200, 820
    yT, yS = 320, 520

    arrow(d, (x0, yT), (x1, yT), BLACK, 5)
    text(d, (x0 - 10, yT - 40), "5′", 22, BLACK, True, anchor="rm")
    text(d, (x1 + 10, yT - 40), "3′", 22, BLACK, True)
    text(d, (x0, yT - 78), "目标链（与 scaffold 反向平行）", 19, BLACK)

    arrow(d, (x1, yS), (x0, yS), BLUE, 5)
    text(d, (x0 - 10, yS + 14), "3′", 22, BLUE, True, anchor="rm")
    text(d, (x1 + 10, yS + 14), "5′", 22, BLUE, True)
    text(d, (x0, yS + 46), "scaffold（5′→3′ 由右向左）", 19, BLUE)

    arrow(d, (x1 - 10, yT + 34), (x0 + 10, yT + 34), RED, 5)
    d.line([(x0 + 10, yT + 34), (x0 + 10, yS - 34)], fill=RED, width=5)
    arrow(d, (x0 + 10, yS - 34), (x1 - 10, yS - 34), RED, 5)
    text(d, (x0 + 22, yT + 8), "staple（红）", 20, RED, True)

    text(
        d,
        (90, 600),
        "两个 domain 在空间上反向平行，骨头在末端折返 180°；\n"
        "配对的连续性由同一条链提供，不存在 helix 间的切换。",
        21,
        BLACK,
    )

    # ---- 右：同向平行 ----
    box_r = (1020, 170, 1940, 720)
    d.rounded_rectangle(box_r, 18, outline=RED, width=3)
    text(d, (1050, 190), "情形 B：两条伙伴链同向平行 → 必须跨 helix 切换", 26, RED, True)

    x0, x1 = 1150, 1800
    yT, yS = 330, 520
    arrow(d, (x0, yT), (x1, yT), BLACK, 5)
    text(d, (x0 - 10, yT - 40), "5′", 22, BLACK, True, anchor="rm")
    text(d, (x1 + 10, yT - 40), "3′", 22, BLACK, True)
    text(d, (x0, yT - 78), "目标链（与 scaffold 同向平行）", 19, BLACK)

    arrow(d, (x0, yS), (x1, yS), BLUE, 5)
    text(d, (x0 - 10, yS + 14), "5′", 22, BLUE, True, anchor="rm")
    text(d, (x1 + 10, yS + 14), "3′", 22, BLUE, True)
    text(d, (x0, yS + 46), "scaffold（5′→3′ 由左向右）", 19, BLUE)

    arrow(d, (x1 - 10, yT + 34), (x0 + 10, yT + 34), RED, 5)
    d.line([(x0 + 10, yT + 34), (x0 + 10, yS - 34)], fill=RED, width=5)
    arrow(d, (x0 + 10, yS - 34), (x1 - 10, yS - 34), RED, 5)

    # 冲突叉号
    d.line([(1400, 470), (1480, 520)], fill=RED, width=10)
    d.line([(1480, 470), (1400, 520)], fill=RED, width=10)
    text(
        d,
        (1050, 600),
        "同向平行的两条链无法被一条只有两个 domain 的链同时配对：\n"
        "配对方向在空间上必须一致，而两条链方向相反 → 中途必须切换 helix。\n"
        "而 helix 间切换就是 crossover，与“不存在 crossover staple”冲突。",
        21,
        (150, 60, 40),
    )

    text(
        d,
        (60, 780),
        "结论：以 scaffold 的 5′→3′ 为基准，被 staple 覆盖的那段目标链必须与 scaffold 反向平行。\n"
        "示意图把三条链都标成「5′→3′ 左→右」，三者不能同时成立（详见文档 §8）。",
        24,
        BLACK,
        True,
    )
    text(
        d,
        (60, 850),
        "对工具的影响：几何只用一个参数记录——scaffold 段与目标段在装配体中是否反向平行；"
        "它决定 staple 序列是 revcomp(scaffold)+revcomp(target) 还是需要按反向重排，也决定坐标推进方向。",
        22,
        GRAY,
    )

    out = OUT / "fig3_staple_geometry.png"
    img.save(out)
    return out


def _on_white(d, xy, s, size=20, fill=BLACK, bold=True):
    """带白色底衬的文字，避免压在细线上看不清。"""
    f = font(size, bold)
    bbox = d.textbbox(xy, s, font=f)
    d.rectangle([bbox[0] - 6, bbox[1] - 4, bbox[2] + 6, bbox[3] + 4], fill="white")
    d.text(xy, s, font=f, fill=fill)


def figure4() -> pathlib.Path:
    """修正后的 staple 平铺模型（对应文档 01）。

    坐标轴单位 = nt（绘图坐标，左端为 0）。scaffold 物理方向 3′→5′ 向右，
    A/B 两条目标链为 5′→3′ 向右，因此 13 nt 域与 26 nt 臂反向平行，折返在左端。
    """
    W, H = 2000, 1120
    img = Image.new("RGB", (W, H), "white")
    d = ImageDraw.Draw(img)

    X0 = 130
    SCALE = 14.0  # px per nt
    N = 104

    def X(nt: float) -> float:
        return X0 + nt * SCALE

    yA, yS, yD, yB = 250, 400, 455, 610

    text(d, (60, 34), "图 4｜修正后的平铺模型：13+13 铺满 scaffold，26 nt 臂铺满目标段", 38, BLACK, True)
    text(
        d,
        (62, 88),
        "坐标轴单位 nt；这是“存储方向”的域图 —— 三条链的箭头都指向各自 5′→3′ 的书写方向，配对时由反向互补自动决定谁反向。",
        21,
        GRAY,
    )

    # 坐标轴
    d.line([(X(0), 178), (X(N), 178)], fill=(170, 170, 170), width=2)
    for nt in range(0, N + 1, 4):
        h = 10 if nt % 13 == 0 else 6
        d.line([(X(nt), 178 - h), (X(nt), 178)], fill=(190, 190, 190), width=1)
    for nt in range(0, N + 1, 8):
        text(d, (X(nt), 154), str(nt), 16, GRAY, anchor="md")

    # A 组 26 nt 臂（上）
    for k in range(4):
        s = 26 * k
        arrow(d, (X(s), yA), (X(s + 26), yA), BLUE, 7)
    text(d, (X(0), yA - 34), "A 组 staple 的 26 nt 目标臂（与 A 段互补，平铺，节距 26）", 20, BLUE, True)

    # B 组 26 nt 臂（下）
    for k in range(4):
        s = 26 * k + 13
        if s >= N:
            break
        arrow(d, (X(s), yB), (X(min(s + 26, N)), yB), RED, 7)
    text(d, (X(13), yB + 18), "B 组 staple 的 26 nt 目标臂（与 B 段互补，平铺，节距 26，相对 A 错开 13）", 20, RED, True)

    # scaffold
    arrow(d, (X(0), yS), (X(N), yS), GREEN, 4)
    text(d, (X(0), yS - 40), "scaffold（存储 5′→3′；装配时物理方向 3′→5′ 向右）", 20, GREEN, True)

    # 13 nt scaffold 结合域：A(蓝)/B(红) 交替
    h = 16
    for k in range(4):
        s = 26 * k
        d.rectangle([X(s), yD - h, X(s + 13), yD + h], fill=BLUE)
        d.rectangle([X(s + 13), yD - h, X(s + 26), yD + h], fill=RED)
        text(d, (X(s + 6.5), yD), "A", 15, "white", True, anchor="mm")
        text(d, (X(s + 19.5), yD), "B", 15, "white", True, anchor="mm")
    _on_white(d, (X(0), yD + 34), "13 nt scaffold 结合域：A 组 / B 组交替，13+13 = 26，正好铺满、无空隙", 20, BLACK, True)

    # 折返连接（细线，画长了只是为看清）
    for k in range(4):
        s = 26 * k
        d.line([(X(s), yD - h), (X(s), yA)], fill=BLUE_LIGHT, width=3)
        d.line([(X(s + 13), yD - h), (X(s + 13), yB)], fill=RED_LIGHT, width=3)
    text(
        d,
        (X(60), 300),
        "竖线 = 13 nt 与 26 nt 之间的那一个磷酸二酯键（画长了只为看清）",
        18,
        GRAY,
    )

    # 右侧推导
    bx0, bx1 = 1290, 1940
    d.rounded_rectangle([bx0, 660, bx1, 1060], 14, outline=(210, 210, 210), width=2)
    text(d, (bx0 + 24, 674), "由此推出的三个硬约束", 24, BLACK, True)
    lines = [
        "① 每条 staple 13 nt，相邻两块不留空隙",
        "   → 覆盖长度 C = 13 n（n = staple 条数）",
        "",
        "② 每组 26 nt 臂首尾相接平铺",
        "   → L(A) = 26 a，L(B) = 26 b",
        "   → 目标总长 = 26(a+b) = 2C",
        "   → 理论上限 2×1200 = 2400 nt（你说的“1.2k×2”）",
        "",
        "③ 交替铺满要求两组条数只差 1：",
        "   |n_A − n_B| ≤ 1  →  L(A) ≈ L(B)",
        "   例：208 + 182（8+7 条），或 364 + 364（14+14 条）",
        "   266 / 126（2:1）与这条不相容",
    ]
    y = 712
    for line in lines:
        color = (150, 60, 40) if line.startswith("③") or line.strip().startswith("392") or line.strip().startswith("最近") or line.strip().startswith("→ a") else BLACK
        text(d, (bx0 + 24, y), line, 19, color)
        y += 26

    out = OUT / "fig4_staple_tiling.png"
    img.save(out)
    return out


def figure5() -> pathlib.Path:
    """结构示意图：目标环上的 oligo / nick / 交界，以及 scaffold 上的 13 nt 平铺。"""
    import json as _json

    data = _json.loads(
        (ROOT / "outputs" / "structure_example.json").read_text(encoding="utf-8")
    )
    W, H = 1700, 1160
    img = Image.new("RGB", (W, H), "white")
    d = ImageDraw.Draw(img)

    text(d, (60, 34), "图 5｜结构示意图（toy example，不含真实序列）", 38, BLACK, True)
    text(
        d,
        (62, 90),
        "目标环 390 nt = A 208 + B 182；7 条 target oligo、7 个 nick、2 个 A/B 交界；"
        "scaffold 覆盖 195 / 1200 nt，15 条 staple。",
        21,
        GRAY,
    )

    def ang(pos, length):
        return -90 + 360 * pos / length

    # ---------- 左：目标环 ----------
    cx, cy, R = 430, 500, 250
    length = data["target"]["length"]
    d.ellipse([cx - R, cy - R, cx + R, cy + R], outline=(230, 230, 230), width=2)

    for oligo in data["oligos"]:
        color = BLUE if oligo["group"] == "A" else RED
        width = 14 if oligo["oligo_type"] == "boundary_spanning" else 9
        a0 = ang(oligo["target_start"], length)
        a1 = ang(oligo["target_end"] % length or 360, length)
        rr = R + (14 if oligo["oligo_type"] == "boundary_spanning" else 0)
        if oligo["target_end"] > length:
            d.arc([cx - rr, cy - rr, cx + rr, cy + rr], a0, 360, fill=color, width=width)
            d.arc([cx - rr, cy - rr, cx + rr, cy + rr], -90, a1, fill=color, width=width)
        else:
            d.arc([cx - rr, cy - rr, cx + rr, cy + rr], a0, a1, fill=color, width=width)

    for nick in data["nicks"]:
        a = math.radians(ang(nick["genomic_or_target_position"], length))
        p0 = (cx + (R - 22) * math.cos(a), cy + (R - 22) * math.sin(a))
        p1 = (cx + (R + 22) * math.cos(a), cy + (R + 22) * math.sin(a))
        d.line([p0, p1], fill=BLACK, width=4)

    for junction in data["junctions"]:
        a = math.radians(ang(junction["position"], length))
        px, py = cx + R * math.cos(a), cy + R * math.sin(a)
        d.ellipse([px - 13, py - 13, px + 13, py + 13], fill="white", outline=GOLD, width=5)
        text(d, (px, py), "J", 15, GOLD, True, anchor="mm")

    text(d, (cx, cy - 40), "目标链", 26, BLACK, True, anchor="mm")
    text(d, (cx, cy), "环状 390 nt", 22, GRAY, anchor="mm")
    text(d, (cx, cy + 34), "7 条 oligo", 22, GRAY, anchor="mm")
    text(d, (cx, cy + 68), "7 个 nick", 22, GRAY, anchor="mm")

    # ---------- 右：scaffold ----------
    cx2, cy2, R2 = 1240, 500, 250
    d.ellipse([cx2 - R2, cy2 - R2, cx2 + R2, cy2 + R2], outline=(225, 225, 225), width=3)
    total = data["scaffold"]["length"]
    for staple in data["staples"]:
        color = BLUE if staple["group"] == "A" else RED
        start = staple["scaffold_start"]
        a0 = ang(start, total)
        a1 = ang((start + 13) % total, total)
        if a1 < a0:
            a1 += 360
        d.arc([cx2 - R2, cy2 - R2, cx2 + R2, cy2 + R2], a0, a1, fill=color, width=14)

    text(d, (cx2, cy2 - 40), "scaffold", 26, BLACK, True, anchor="mm")
    text(d, (cx2, cy2), f"环状 {total} nt", 22, GRAY, anchor="mm")
    text(d, (cx2, cy2 + 34), "覆盖 195 nt", 22, GREEN, True, anchor="mm")
    text(d, (cx2, cy2 + 68), "15 条 staple × 13 nt", 22, GRAY, anchor="mm")
    text(d, (cx2, cy2 + 102), "其余 1005 nt 冗余", 22, GRAY, anchor="mm")

    # ---------- 底部图例 ----------
    top = 810
    d.line([(60, top - 20), (W - 60, top - 20)], fill=LIGHT, width=3)
    legend = [
        (BLUE, "A 组 target oligo / A 组 staple 的 13 nt 块"),
        (RED, "B 组 target oligo / B 组 staple 的 13 nt 块"),
        (GOLD, "J = A/B 交界（由 boundary-spanning oligo 连续覆盖，不是 nick）"),
        (BLACK, "短径向线 = nick（3′-OH / 5′-P 相邻但未共价连接）"),
    ]
    for i, (color, label) in enumerate(legend):
        y = top + 10 + i * 40
        d.line([(70, y + 12), (130, y + 12)], fill=color, width=8 if color != BLACK else 3)
        text(d, (146, y), label, 23, BLACK)

    notes = [
        "外移到 R+14 的那两段弧 = boundary-spanning oligo（跨过 J）；交界处序列连续，但这不等于环已经共价闭合。",
        "scaffold 上每条 staple 占 13 nt，相邻两块交替属于 A/B，13+13 = 26 正好铺满，互不侵占。",
    ]
    for i, line in enumerate(notes):
        text(d, (60, top + 180 + i * 34), line, 21, (150, 60, 40))

    out = OUT / "fig5_structure_example.png"
    img.save(out)
    return out


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    for fn in (figure1, figure2, figure3, figure4, figure5):
        print("written:", fn())


if __name__ == "__main__":
    main()
