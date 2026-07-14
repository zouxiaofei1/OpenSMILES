# -*- coding: utf-8 -*-
"""Parse workstate.md → accuracy series + chart + analysis markdown."""
from __future__ import annotations

import csv
import json
import os
import re
from collections import OrderedDict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.patches as mpatches  # noqa: E402
from matplotlib import font_manager  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
WS = ROOT / "workstate.md"
OUT_DIR = ROOT / "tools"


def _setup_font() -> None:
    for fp in (
        r"C:\Windows\Fonts\msyh.ttc",
        r"C:\Windows\Fonts\simhei.ttf",
        r"C:\Windows\Fonts\simsun.ttc",
    ):
        if os.path.exists(fp):
            font_manager.fontManager.addfont(fp)
            prop = font_manager.FontProperties(fname=fp)
            plt.rcParams["font.family"] = prop.get_name()
            break
    plt.rcParams["axes.unicode_minus"] = False


def parse_workstate(text: str) -> list[dict]:
    lines = [ln.strip() for ln in text.splitlines() if ln.strip().startswith("[#")]
    pat_arrow = re.compile(
        r"dual\s+(\d+(?:\.\d+)?)%\s*(?:\((\d+)\))?\s*→\s*"
        r"(\d+(?:\.\d+)?)%\s*(?:\((\d+)\))?",
        re.I,
    )
    pat_flat = re.compile(
        r"dual\s+(\d+(?:\.\d+)?)%\s*(?:\((\d+)\))?\s*(?:持平|held)",
        re.I,
    )
    pat_fail = re.compile(r"fails\s+(\d+)\s*→\s*(\d+)", re.I)
    pat_commit = re.compile(r"^\[#([0-9a-f]+)\]", re.I)
    pat_tag = re.compile(r"\[(IUPAC[^\]]*|架构[^\]]*)\]")

    rows: list[dict] = []
    prev_pct, prev_ok = 0.3, None
    rows.append(
        {
            "round": 0,
            "commit": "baseline",
            "title": "起点（日志前）",
            "tag": "baseline",
            "dual_pct": 0.3,
            "ok_dual": None,
            "fails": None,
            "fails_before": None,
            "delta_pct": 0.0,
            "delta_ok": None,
            "is_arch": False,
            "is_flat": False,
            "raw": "",
        }
    )

    for i, ln in enumerate(lines, 1):
        m_c = pat_commit.search(ln)
        commit = m_c.group(1) if m_c else ""
        m_i = pat_tag.search(ln)
        tag = m_i.group(1) if m_i else ""
        parts = re.split(r"\]\s*", ln, maxsplit=2)
        title = ""
        if len(parts) >= 3:
            title = re.split(r"\s*\[", parts[2], maxsplit=1)[0].strip()
        is_arch = any(
            k in (tag + title)
            for k in (
                "架构",
                "registry",
                "Kind",
                "ring_producers",
                "fused56",
                "LeafHandler",
                "ParentKind",
            )
        )

        m_a = pat_arrow.search(ln)
        m_f = pat_flat.search(ln)
        m_fail = pat_fail.search(ln)
        fails_after = int(m_fail.group(2)) if m_fail else None
        fails_before = int(m_fail.group(1)) if m_fail else None

        if m_a:
            pct_b, ok_b, pct_a, ok_a = (
                m_a.group(1),
                m_a.group(2),
                m_a.group(3),
                m_a.group(4),
            )
            dual_pct = float(pct_a)
            ok_dual = int(ok_a) if ok_a else None
            delta_pct = float(pct_a) - float(pct_b)
            delta_ok = (int(ok_a) - int(ok_b)) if (ok_a and ok_b) else None
            is_flat = abs(delta_pct) < 1e-9 and (delta_ok in (None, 0))
        elif m_f:
            dual_pct = float(m_f.group(1))
            ok_dual = int(m_f.group(2)) if m_f.group(2) else prev_ok
            delta_pct = 0.0
            delta_ok = 0
            is_flat = True
            if m_f.group(2) is None:
                ok_dual = prev_ok
        else:
            m_any = re.search(
                r"dual\s+(\d+(?:\.\d+)?)%\s*(?:\((\d+)\))?", ln
            )
            if not m_any:
                raise ValueError(f"unparsed line {i}: {ln[:120]}")
            dual_pct = float(m_any.group(1))
            ok_dual = int(m_any.group(2)) if m_any.group(2) else prev_ok
            delta_pct = dual_pct - prev_pct
            delta_ok = (
                (ok_dual - prev_ok)
                if (ok_dual is not None and prev_ok is not None)
                else None
            )
            is_flat = abs(delta_pct) < 1e-9 and (delta_ok in (None, 0))

        rows.append(
            {
                "round": i,
                "commit": commit,
                "title": title,
                "tag": tag,
                "dual_pct": dual_pct,
                "ok_dual": ok_dual,
                "fails": fails_after,
                "fails_before": fails_before,
                "delta_pct": round(delta_pct, 4),
                "delta_ok": delta_ok,
                "is_arch": is_arch,
                "is_flat": is_flat,
                "raw": ln,
            }
        )
        prev_pct = dual_pct
        if ok_dual is not None:
            prev_ok = ok_dual

    n_dual = 4062
    for r in rows:
        if r["ok_dual"] is None and r["fails"] is not None:
            r["ok_dual"] = n_dual - r["fails"]
        if r["ok_dual"] is None:
            r["ok_dual_est"] = int(round(r["dual_pct"] / 100.0 * n_dual))
        else:
            r["ok_dual_est"] = r["ok_dual"]

    for i in range(1, len(rows)):
        a, b = rows[i - 1], rows[i]
        if b["delta_ok"] is None:
            b["delta_ok"] = int(b["ok_dual_est"] - a["ok_dual_est"])
    return rows


def phase_of(r: dict) -> str:
    rd = r["round"]
    if rd == 0:
        return "0 起点"
    if rd <= 12:
        return "1 开链一元 FG"
    if rd <= 23:
        return "2 单环碳环"
    if rd <= 45:
        return "3 多官能/不饱和/N-杂"
    if rd <= 58:
        return "4 苯系与芳 FG"
    if rd <= 91:
        return "5 杂芳/稠环保留母体"
    if rd <= 110:
        return "6 递归取代基"
    if rd <= 121:
        return "7 架构加固"
    if rd <= 122:
        return "8 词干 C11–C35"
    return "9 阴离子/不饱和贯通"


PHASE_COLORS = {
    "0 起点": "#9e9e9e",
    "1 开链一元 FG": "#1f77b4",
    "2 单环碳环": "#ff7f0e",
    "3 多官能/不饱和/N-杂": "#2ca02c",
    "4 苯系与芳 FG": "#d62728",
    "5 杂芳/稠环保留母体": "#9467bd",
    "6 递归取代基": "#8c564b",
    "7 架构加固": "#e377c2",
    "8 词干 C11–C35": "#17becf",
    "9 阴离子/不饱和贯通": "#bcbd22",
}


# 高增益轮次的机制标签（顺序=优先级；更具体的规则在前）
THEME_RULES: list[tuple[str, tuple[str, ...]]] = [
    ("架构/registry", ("架构", "Kind 注册", "ring_producers", "fused56", "LeafHandler", "ParentKind", "候选收集", "多候选生产者")),
    ("词干/词表批量扩展", ("C11", "C20+", "C11–C35", "倍数词头", "半系统表", "zh_stem", "程序 compose")),
    ("母体形态变体", ("羧酸根", "酸根", "anion", "oyl chloride", "酰氯", "anhydride", "酸酐", "alkanoate/酸根")),
    ("正交特征组合", ("羟基烯酸", "烯酸前缀", "polyalkenol", "alkenol", "alkenoic", "alkenal", "alkenenitrile", "alkenoate", "alkenedioic", "不饱和酸", "不饱和醛", "不饱和腈", "不饱和醇", "不饱和酯", "不饱和二酸", "amino/oxo", "prefix_alkenoic", "E/Z", "(E)/(Z)")),
    ("酸母体前缀提取", ("hydroxyalkanoic", "aminoalkanoic", "oxoalkanoic", "羟基酸", "氨基酸", "酮酸", "提取 hydroxy", "提取 amino", "提取 oxo")),
    ("通用前缀贯通", ("卤素前缀", "直链单烷基", "fluoro/chloro", "nitro", "methoxy/ethoxy", "alkoxy", "trifluoromethyl", "三氟甲基", "isobutyl", "sec-butyl", "neopentyl", "tert-butyl")),
    ("保留母体+取代放宽", ("保留母体", "phenol", "aniline", "benzoic", "benzaldehyde", "acetophenone", "上限", "2→3", "1→2", "简单取代", "xylene", "多取代苯", "简单苯", "alkylbenzene", "aminopheno")),
    ("新高频 FG 类", ("一元羧酸", "一元酮", "一元醛", "一元酯", "一元伯胺", "一元醇", "一元烯烃", "一元炔烃", "二元醇", "二元羧酸", "二元酮", "二元伯胺", "三元醇", "硫醇", "二烷基醚", "二烷基硫醚", "伯酰胺", "一元腈")),
    ("单点保留环母体", ("pyridine", "furan", "thiophene", "pyrrole", "imidazole", "pyrazole", "indole", "naphthalene", "quinoline", "quinazoline", "quinoxaline", "anthracene", "diazine", "morpholine", "piperidine", "oxolane", "aziridine", "pyrimidin", "pyridazine", "pyrazine")),
]


def theme_of(r: dict) -> str:
    blob = f"{r.get('tag', '')} {r.get('title', '')}"
    for name, keys in THEME_RULES:
        if any(k.lower() in blob.lower() for k in keys):
            return name
    return "其他"


def phase_stats(rows: list[dict]) -> list[tuple]:
    order: list[str] = []
    for r in rows:
        if r["phase"] not in order:
            order.append(r["phase"])
    out = []
    for ph in order:
        rs = [r for r in rows if r["phase"] == ph]
        d_ok = sum((r.get("delta_ok") or 0) for r in rs if r["round"] > 0)
        d_pct = sum((r.get("delta_pct") or 0) for r in rs if r["round"] > 0)
        n = sum(1 for r in rs if r["round"] > 0)
        out.append((ph, d_pct, d_ok, n, rs[0]["dual_pct"], rs[-1]["dual_pct"], rs[0]["round"], rs[-1]["round"]))
    return out


def plot(rows: list[dict], out: Path) -> None:
    _setup_font()
    rounds = [r["round"] for r in rows]
    pcts = [r["dual_pct"] for r in rows]
    oks = [r.get("ok_dual_est") for r in rows]
    deltas = [r.get("delta_ok") or 0 for r in rows]
    order: list[str] = []
    for r in rows:
        if r["phase"] not in order:
            order.append(r["phase"])
    stats = phase_stats(rows)

    fig = plt.figure(figsize=(14, 9.2), dpi=140)
    gs = fig.add_gridspec(3, 1, height_ratios=[3.2, 1.35, 1.7], hspace=0.38)
    ax = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1], sharex=ax)
    ax3 = fig.add_subplot(gs[2])

    for ph in order:
        rs = [r["round"] for r in rows if r["phase"] == ph]
        if not rs:
            continue
        c = PHASE_COLORS.get(ph, "#ccc")
        ax.axvspan(min(rs) - 0.5, max(rs) + 0.5, color=c, alpha=0.12, lw=0)
        ax2.axvspan(min(rs) - 0.5, max(rs) + 0.5, color=c, alpha=0.10, lw=0)

    ax.plot(
        rounds,
        pcts,
        color="#1565c0",
        lw=2.0,
        marker="o",
        ms=3.2,
        zorder=3,
        label="dual 准确率 %",
    )
    for r in rows:
        if (r.get("delta_ok") or 0) >= 8 or (r.get("delta_pct") or 0) >= 0.25:
            ax.scatter(
                [r["round"]],
                [r["dual_pct"]],
                s=46,
                color="#c62828",
                zorder=4,
                edgecolors="white",
                linewidths=0.6,
            )
        if (r.get("delta_ok") or 0) >= 10 or (r.get("delta_pct") or 0) >= 0.35:
            ax.annotate(
                f"+{r['delta_ok']}  {r['dual_pct']}%\n{r['title'][:16]}",
                xy=(r["round"], r["dual_pct"]),
                xytext=(6, 12),
                textcoords="offset points",
                fontsize=7,
                color="#b71c1c",
                arrowprops=dict(arrowstyle="-", color="#ef9a9a", lw=0.6),
            )
        if r.get("is_arch"):
            ax.axvline(r["round"], color="#ad1457", ls="--", lw=0.65, alpha=0.5, zorder=1)

    ax.set_ylabel("dual 准确率 (%)")
    ax.set_title(
        "NamePredict benchmark dual 准确率 vs 迭代轮次（workstate.md）",
        fontsize=13,
        pad=10,
    )
    ax.set_ylim(0, max(pcts) * 1.18)
    ax.grid(True, axis="y", alpha=0.35)
    ax.annotate(
        f"终态 {pcts[-1]}%  ({oks[-1]} / ~4062)",
        xy=(rounds[-1], pcts[-1]),
        xytext=(-130, -28),
        textcoords="offset points",
        fontsize=9,
        color="#0d47a1",
        arrowprops=dict(arrowstyle="->", color="#0d47a1", lw=0.8),
    )
    handles = [
        mpatches.Patch(color=PHASE_COLORS[p], alpha=0.55, label=p) for p in order
    ]
    handles.append(ax.lines[0])
    ax.legend(handles=handles, loc="upper left", fontsize=7.5, framealpha=0.92)

    colors = [
        "#2e7d32" if d > 0 else ("#c62828" if d < 0 else "#9e9e9e") for d in deltas
    ]
    ax2.bar(rounds, deltas, color=colors, width=0.8, zorder=2)
    ax2.axhline(0, color="#424242", lw=0.8)
    ax2.set_ylabel("Δok_dual（条）")
    ax2.set_xlabel("迭代轮次（日志序号，0=起点）")
    ax2.grid(True, axis="y", alpha=0.3)
    ax2.set_title("每轮 dual 通过条数变化（绿=增益，红=回退，灰=持平）", fontsize=10)

    labels = [p[0] for p in stats]
    vals = [p[2] for p in stats]
    cols = [PHASE_COLORS.get(p[0], "#888") for p in stats]
    bars = ax3.barh(labels, vals, color=cols, alpha=0.88)
    ax3.set_xlabel("阶段累计 Δok_dual（条）")
    ax3.set_title("各阶段 dual 净增益对比", fontsize=10)
    ax3.grid(True, axis="x", alpha=0.3)
    for bar, item in zip(bars, stats):
        ph, d_pct, d_ok, n, sp, ep, sr, er = item
        ax3.text(
            bar.get_width() + 0.6,
            bar.get_y() + bar.get_height() / 2,
            f"+{d_ok} 条 / +{d_pct:.1f}pp  ({n}轮 → {ep}%)",
            va="center",
            fontsize=8,
        )

    fig.savefig(out, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def theme_stats(rows: list[dict]) -> list[dict]:
    by: dict[str, dict] = {}
    for r in rows:
        if r["round"] <= 0:
            continue
        th = r.get("theme") or theme_of(r)
        slot = by.setdefault(
            th,
            {
                "theme": th,
                "n": 0,
                "delta_ok": 0,
                "delta_pct": 0.0,
                "pos": 0,
                "flat": 0,
                "neg": 0,
                "max_ok": 0,
                "examples": [],
            },
        )
        dok = r.get("delta_ok") or 0
        slot["n"] += 1
        slot["delta_ok"] += dok
        slot["delta_pct"] += r.get("delta_pct") or 0.0
        if dok > 0:
            slot["pos"] += 1
        elif dok < 0:
            slot["neg"] += 1
        else:
            slot["flat"] += 1
        if dok > slot["max_ok"]:
            slot["max_ok"] = dok
        slot["examples"].append((dok, r["round"], r["title"], r["commit"]))
    out = list(by.values())
    for s in out:
        s["eff"] = s["delta_ok"] / s["n"] if s["n"] else 0.0
        s["examples"] = sorted(s["examples"], key=lambda x: -x[0])[:3]
    out.sort(key=lambda s: (-s["delta_ok"], -s["eff"]))
    return out


def write_analysis(
    rows: list[dict], stats: list[tuple], themes: list[dict], md_path: Path
) -> None:
    gains = sorted(
        [r for r in rows if r["round"] > 0],
        key=lambda r: -(r.get("delta_ok") or 0),
    )[:12]
    top10 = gains[:10]
    negs = [r for r in rows if (r.get("delta_ok") or 0) < 0]
    flats = [r for r in rows if r["round"] > 0 and (r.get("delta_ok") == 0)]
    start, end = rows[0], rows[-1]
    net_pp = end["dual_pct"] - start["dual_pct"]
    net_ok = (end.get("ok_dual_est") or 0) - (start.get("ok_dual_est") or 0)

    def win_avg(w: int) -> tuple[float, int]:
        recent = [r for r in rows if r["round"] > len(rows) - 1 - w]
        s = sum(r.get("delta_ok") or 0 for r in recent)
        return s / max(len(recent), 1), s

    # top10 theme concentration
    top_theme_count: dict[str, int] = {}
    top_theme_ok: dict[str, int] = {}
    for r in top10:
        th = r.get("theme") or theme_of(r)
        top_theme_count[th] = top_theme_count.get(th, 0) + 1
        top_theme_ok[th] = top_theme_ok.get(th, 0) + (r.get("delta_ok") or 0)

    lines = []
    lines.append("# workstate 准确率–迭代分析\n\n")
    lines.append("## 总览\n\n")
    lines.append(f"- 数据源: `workstate.md` 日志 **{len(rows)-1}** 轮（含起点共 {len(rows)} 点）\n")
    lines.append("- 指标: benchmark **dual**（中英文双过；n_dual ≈ 4062）\n")
    lines.append(
        f"- 轨迹: **{start['dual_pct']}% → {end['dual_pct']}%**"
        f"（约 {start.get('ok_dual_est')} → {end.get('ok_dual_est')} 条）\n"
    )
    lines.append(f"- 净增益: **+{net_pp:.1f} pp / +{net_ok} 条**\n")
    lines.append("- 图: [`tools/workstate_accuracy_chart.png`](workstate_accuracy_chart.png)\n")
    lines.append("- 序列: [`tools/workstate_accuracy_series.csv`](workstate_accuracy_series.csv)\n\n")

    lines.append("## 阶段汇总\n\n")
    lines.append("| 阶段 | 轮次 | dual% | Δpp | Δok | 轮数 | 效率(条/轮) |\n")
    lines.append("|---|---:|---:|---:|---:|---:|---:|\n")
    for ph, d_pct, d_ok, n, sp, ep, sr, er in stats:
        eff = (d_ok / n) if n else 0
        lines.append(
            f"| {ph} | {sr}–{er} | {sp}→{ep} | +{d_pct:.2f} | +{d_ok} | {n} | {eff:.2f} |\n"
        )

    lines.append("\n## Top 增益轮次（按 Δok）\n\n")
    lines.append("| 轮次 | commit | Δok | Δpp | dual% | 主题 | 选题 |\n")
    lines.append("|---:|---|---:|---:|---:|---|---|\n")
    for r in gains:
        th = r.get("theme") or theme_of(r)
        lines.append(
            f"| {r['round']} | `{r['commit']}` | +{r.get('delta_ok')} | "
            f"{r.get('delta_pct'):+.2f} | {r['dual_pct']} | {th} | {r['title'][:60]} |\n"
        )

    lines.append("\n## Top10 高增益提升与什么有关\n\n")
    lines.append(
        "按 **Δok** 取 top10（合计 +"
        f"{sum(r.get('delta_ok') or 0 for r in top10)} 条，"
        f"约占全程净增益 {100 * sum(r.get('delta_ok') or 0 for r in top10) / max(net_ok, 1):.0f}%）。"
        "主题分布：\n\n"
    )
    lines.append("| 主题 | top10 次数 | top10 合计 Δok | 机制解释 |\n|---|---:|---:|---|\n")
    theme_explain = {
        "词干/词表批量扩展": "一次改动覆盖多个碳数/倍数词头，命中 benchmark 长尾词干面",
        "母体形态变体": "已有母体的电荷/盐/酰卤等同族形态，复用既有链选择与词干",
        "正交特征组合": "已有 FG × 不饱和/立体/前缀 的笛卡尔积；修门控即解锁一批组合分子",
        "酸母体前缀提取": "酸保留母体 + L3 前缀抽取；hydroxy/amino/oxo 同类模式可复制",
        "通用前缀贯通": "卤/烷基等横切前缀一次贯通 L2/L3/L5，多母体共享",
        "保留母体+取代放宽": "已有 retained parent 的取代上限/允许集放宽，覆盖度跳变",
        "新高频 FG 类": "空白官能团类首次落地，覆盖 benchmark 高频一元母体",
        "单点保留环母体": "每轮一个环骨架；单点增益小、持平多",
        "架构/registry": "当轮 dual 常持平，价值在后续扩展斜率",
        "其他": "难归类或混合改动",
    }
    for th, cnt in sorted(top_theme_count.items(), key=lambda x: -top_theme_ok[x[0]]):
        lines.append(
            f"| {th} | {cnt} | +{top_theme_ok[th]} | {theme_explain.get(th, '')} |\n"
        )

    lines.append("\n### 共性规律（可操作）\n\n")
    lines.append(
        "1. **批量面 > 单点**：词干表、取代上限、前缀门控一次改动扫过多个分子；"
        "单点 quinazoline 类常 0～+1。\n"
    )
    lines.append(
        "2. **正交组合 > 新环骨架**：在已通的酸/醇母体上放开 C=C、立体、amino/oxo/OH 前缀，"
        "近期 r123–126 合计 +56，效率远高于继续堆未取代杂芳。\n"
    )
    lines.append(
        "3. **同族形态复用**：羧酸→羧酸根、酸→酰氯/酐，复用 L1/L2/L5 路径，"
        "边际实现成本低、Δok 高。\n"
    )
    lines.append(
        "4. **横切前缀仍值钱**：早期卤/烷基贯通仍在 top；"
        "后续同类是 nitro/alkoxy/haloalkyl/立体前缀在更多母体上的统一。\n"
    )
    lines.append(
        "5. **架构轮不进 top10 但必要**：无当轮分数时别用 dual 否决；"
        "用它降低下一次“正交组合/同族形态”的实现成本。\n"
    )

    lines.append("\n## 主题全量效率（预测输入）\n\n")
    lines.append("| 主题 | 轮数 | 合计Δok | 效率(条/轮) | 正/平/负 | 最大单轮 | 代表 |\n")
    lines.append("|---|---:|---:|---:|---:|---:|---|\n")
    for s in themes:
        ex = s["examples"][0][2][:28] if s["examples"] else ""
        lines.append(
            f"| {s['theme']} | {s['n']} | +{s['delta_ok']} | {s['eff']:.2f} | "
            f"{s['pos']}/{s['flat']}/{s['neg']} | +{s['max_ok']} | {ex} |\n"
        )

    lines.append("\n## 能否预测最佳改动项？\n\n")
    lines.append(
        "**可以做方向级预测，不能精确到“下一轮 +N 条”。** "
        "历史 top 提升高度集中在少数机制；用机制先验 + fail 桶可排序候选，"
        "但具体 Δok 仍取决于 fail 面有多大、门控是否误伤。\n\n"
    )
    lines.append("### 预测启发式（按期望 Δok 降序）\n\n")
    lines.append(
        "| 优先级 | 改动形态 | 历史依据 | 预期 |\n"
        "|---:|---|---|---|\n"
        "| 1 | **词干/词表批量**（更长链、支化倍数、zh/en 对称表） | r122 +42，效率极高 | 若 fail 仍有长链词干缺口，+10～40 |\n"
        "| 2 | **正交特征贯通**（已有母体 × 不饱和/立体/多前缀） | r125 +21, r126 +6；阶段9 效率 ~14 条/轮 | 每打通一条门控 +5～25 |\n"
        "| 3 | **同族形态**（酸根/盐/酰卤/酯对称臂） | r123 +24 羧酸根 | 一类形态 +8～25 |\n"
        "| 4 | **通用前缀/leaf 横切**（haloalkyl、nitro、alkoxy 全母体） | 早期卤/烷基 top；横切摊销 | +5～15，偶发回归风险 |\n"
        "| 5 | **保留母体取代放宽**（上限 2→3、允许集） | r66 +10, r46/47 +9 | +5～12，需防过度 claim |\n"
        "| 6 | **新高频 FG 空白类** | r6/7/11 酸/酮/酯 | 仅当 fail 桶显示空白类够大 |\n"
        "| 7 | **单点未取代环母体** | 大量 0～+1 持平轮 | 默认低优先级 |\n"
        "| 8 | **纯架构** | 当轮≈0 | 为 1–4 铺路时做，不指望 dual |\n\n"
    )
    lines.append("### 下一轮候选（由日志规律外推，非已实现）\n\n")
    lines.append(
        "1. **alkynoic / 炔醇 / 多炔**：复制 alkenoic/alkenol 正交路径（阶段9 同构）。\n"
        "2. **多前缀酸**（hydroxy+amino、halo+oxo、nitro 酸）与 **立体前缀** 在更多 FG 上统一。\n"
        "3. **羧酸盐/金属盐**、**多羧酸根**：r123 同族延伸。\n"
        "4. **酯/酰胺 N-侧与酰基侧不对称复杂化**（已有简单 N-alkyl 基础）。\n"
        "5. **词干缺口扫 fail**：C36+、iso/sec/tert 与中文数字词干是否仍成批失败。\n"
        "6. 避免：再堆一个未取代稠杂芳（历史效率最低区间）。\n\n"
    )
    lines.append(
        "### 可操作的选题公式\n\n"
        "```\n"
        "期望收益 ≈ fail桶大小 × 覆盖倍率 × 实现通用度\n"
        "覆盖倍率: 词表批量≈高, 正交门控≈高, 同族形态≈中高, 单点母体≈低\n"
        "实现通用度: 横切 L2/L3/L5 共享 > 单 kind 特判\n"
        "风险折扣: 触及母体优先级/开链vs环边界时 ×0.5（历史上唯一回退在此）\n"
        "```\n\n"
        "实操：每轮先从最新 fail 聚成桶（按 parent kind / 缺失前缀 / 词干长度），"
        "只做 **桶≥8 且属于优先级 1–5** 的项；架构轮穿插在连续 2 次特判之后。\n"
    )

    lines.append("\n## 回退 / 持平\n\n")
    if negs:
        lines.append("| 轮次 | commit | Δok | dual% | 选题 |\n|---:|---|---:|---:|---|\n")
        for r in negs:
            lines.append(
                f"| {r['round']} | `{r['commit']}` | {r.get('delta_ok')} | "
                f"{r['dual_pct']} | {r['title'][:70]} |\n"
            )
    else:
        lines.append("- 无回退轮\n")
    lines.append(
        f"\n- 持平轮: **{len(flats)} / {len(rows)-1}**"
        f"（{100 * len(flats) / (len(rows) - 1):.1f}%）\n"
    )

    a10, s10 = win_avg(10)
    a20, s20 = win_avg(20)
    a30, s30 = win_avg(30)
    lines.append("\n## 边际效率\n\n")
    lines.append("| 窗口 | 总 Δok | 平均条/轮 |\n|---|---:|---:|\n")
    lines.append(f"| 近 10 轮 | {s10} | {a10:.2f} |\n")
    lines.append(f"| 近 20 轮 | {s20} | {a20:.2f} |\n")
    lines.append(f"| 近 30 轮 | {s30} | {a30:.2f} |\n")
    lines.append(f"| 全程 | {net_ok} | {net_ok / max(len(rows)-1, 1):.2f} |\n")

    lines.append("\n## 解读要点\n\n")
    lines.append(
        "1. **早期开链 FG（r1–12）**：空白类覆盖，dual 0.3%→~2.3%。\n"
    )
    lines.append(
        "2. **苯系+芳 FG（r46–58）与杂芳保留（r59–91）**：第二段升浪；"
        "后期单点环母体效率衰减。\n"
    )
    lines.append(
        "3. **递归取代基（r92–110）**：正增益但持平增多，深度/嵌套覆盖变窄。\n"
    )
    lines.append(
        "4. **架构轮 dual 常持平**：可扩展性投资；词干 C11–C35（r122 +42）"
        "是日志最大单轮涨幅。\n"
    )
    lines.append(
        "5. **阶段9（r123–126）再次陡升**：羧酸根 +24、羟基烯酸 +21、"
        "amino/oxo 烯酸 +6——说明当前主战场是 **已有母体的形态/正交扩展**，"
        "不是新环。\n"
    )
    lines.append(
        "6. **回退极少**（registry 叶边界 -1）：风险集中在母体/侧链边界，"
        "正交门控放行时要写回归用例。\n"
    )
    lines.append(
        f"7. 当前 dual **{end['dual_pct']}%**（~{end.get('ok_dual_est')} 条）；"
        "长尾仍是复杂药物/糖肽/多环立体；优先 **fail 桶驱动的批量/正交/同族**，"
        "停掉低效单点保留母体。\n"
    )

    md_path.write_text("".join(lines), encoding="utf-8")


def main() -> None:
    rows = parse_workstate(WS.read_text(encoding="utf-8"))
    for r in rows:
        r["phase"] = phase_of(r)
        r["theme"] = theme_of(r)

    json_path = OUT_DIR / "workstate_accuracy_series.json"
    csv_path = OUT_DIR / "workstate_accuracy_series.csv"
    chart_path = OUT_DIR / "workstate_accuracy_chart.png"
    md_path = OUT_DIR / "workstate_accuracy_analysis.md"

    json_path.write_text(
        json.dumps(rows, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    fields = [k for k in rows[0] if k != "raw"]
    with csv_path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k) for k in fields})

    stats = phase_stats(rows)
    themes = theme_stats(rows)
    plot(rows, chart_path)
    write_analysis(rows, stats, themes, md_path)

    print(f"points={len(rows)} rounds={len(rows)-1}")
    print(
        f"start={rows[0]['dual_pct']}% end={rows[-1]['dual_pct']}% "
        f"ok={rows[-1].get('ok_dual')}"
    )
    print("wrote", chart_path)
    print("wrote", md_path)
    print("wrote", csv_path)
    for item in stats:
        ph, d_pct, d_ok, n, sp, ep, sr, er = item
        print(f"{ph}: r{sr}-{er} {sp}%->{ep}% dpp={d_pct:+.2f} dok={d_ok:+} n={n}")
    print("--- themes by delta_ok ---")
    for s in themes:
        print(
            f"{s['theme']}: n={s['n']} dok=+{s['delta_ok']} "
            f"eff={s['eff']:.2f} max=+{s['max_ok']}"
        )
    print("--- top10 ---")
    gains = sorted(
        [r for r in rows if r["round"] > 0],
        key=lambda r: -(r.get("delta_ok") or 0),
    )[:10]
    for r in gains:
        print(
            f"r{r['round']} +{r.get('delta_ok')} {r.get('theme')} | {r['title'][:50]}"
        )


if __name__ == "__main__":
    main()
