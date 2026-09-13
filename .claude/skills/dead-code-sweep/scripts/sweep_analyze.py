"""交叉验证 + 分档 + 生成 dead-code-report.md。

三档口径:
  A 从未执行      —— 全量跑过也没进过该函数（覆盖证据）
  B 执行但无效果  —— 进过函数，但注入 return None 后逐行预测不变
  C 有用          —— 注入后预测改变，或注入后模块导入失败（导入期必需）

"A/B 类"只是候选，**能不能删还要看有没有下游引用**：被 tests/server/benchmarks
引用的，删了会 fail 测试或打断服务端，需同步修改。
"""

from __future__ import annotations

import json
import re
import sys
from collections import defaultdict
from pathlib import Path

import mutlib
import sweepcfg

REPORT = sweepcfg.ROOT / "dead-code-report.md"


def load_downstream() -> dict[str, str]:
    """扫主仓库的下游文件（tests/benchmarks/server/tools 顶层 .py）。

    注意 tools/<pkg>-v2、-v3 是独立引擎副本，不是本包的下游，须排除。
    """
    out: dict[str, str] = {}
    for d in sweepcfg.DOWNSTREAM:
        base = sweepcfg.REPO / d
        if not base.is_dir():
            continue
        it = base.glob("*.py") if d == "tools" else base.rglob("*.py")
        for p in it:
            if "__pycache__" in p.parts:
                continue
            rel = p.relative_to(sweepcfg.REPO).as_posix()
            if any(rel.startswith(x.rstrip("/") + "/") for x in sweepcfg.DOWNSTREAM_EXCLUDE):
                continue
            out[rel] = p.read_text(encoding="utf-8", errors="replace")
    return out


def hits(name: str, downstream: dict[str, str]) -> list[str]:
    pat = re.compile(r"\b" + re.escape(name) + r"\b")
    return [f for f, t in downstream.items() if pat.search(t)]


def main() -> None:
    if not sweepcfg.RESULTS.is_file():
        print("还没有扫描结果，先跑 sweep_drive.py")
        return
    mut = json.loads(sweepcfg.RESULTS.read_text(encoding="utf-8"))
    cov: dict[str, int] = {}
    if sweepcfg.COVERAGE.is_file():
        cov = json.loads(sweepcfg.COVERAGE.read_text(encoding="utf-8"))
    executed = {tuple(k.split("::", 1)) for k in cov if not k.endswith("::<module>")}
    imported_mods = {k.split("::")[0] for k in cov if k.endswith("::<module>")}

    def is_executed(key: str) -> bool:
        f, q = key.split("::", 1)
        if (f, q) in executed:
            return True
        # 嵌套函数 co_qualname 形如 outer.<locals>.inner，按末段宽松匹配同文件同名
        return any(rf == f and rq.split(".")[-1] == q.split(".")[-1] for rf, rq in executed)

    downstream = load_downstream()
    recs = {r["key"]: r for r in mutlib.collect()}
    shown = {k: v for k, v in mut.items() if k not in sweepcfg.SHIMS}

    cand: list[dict] = []
    for key, v in sorted(shown.items()):
        if v.get("status") != "same":
            continue
        rec = recs.get(key, {})
        name = key.split("::", 1)[1].split(".")[-1]
        cand.append({
            "key": key,
            "file": key.split("::", 1)[0].replace(f"{sweepcfg.SRC.name}/{sweepcfg.PKG}/", ""),
            "qual": key.split("::", 1)[1],
            "cls": "B 执行但无效果" if is_executed(key) else "A 从未执行",
            "downstream": hits(name, downstream),
            "rec": rec,
        })

    clean = [c for c in cand if not c["downstream"]]
    linked = [c for c in cand if c["downstream"]]
    inconsistent = [k for k, v in shown.items() if v.get("status") == "diff" and not is_executed(k)]
    essential = [k for k, v in shown.items() if v.get("status") == "import_error"]
    never_mod = sorted({c["file"] for c in cand if not is_executed(c["key"])})

    def block(items: list[dict]) -> str:
        g: dict[str, list[dict]] = defaultdict(list)
        for c in items:
            g[c["file"]].append(c)
        out = []
        for f in sorted(g):
            out.append(f"\n**`{f}`**\n")
            out.append("| 函数 | 档 | 下游引用 |")
            out.append("|---|---|---|")
            for c in sorted(g[f], key=lambda x: x["qual"]):
                ref = ", ".join(f"`{x}`" for x in c["downstream"][:3]) or "—"
                if len(c["downstream"]) > 3:
                    ref += f" 等 {len(c['downstream'])} 处"
                out.append(f"| `{c['qual']}` | {c['cls']} | {ref} |")
            out.append("")
        return "\n".join(out)

    md = [
        "# 死代码扫描报告（变异法）",
        "",
        f"- 包: `{sweepcfg.SRC.name}/{sweepcfg.PKG}`  函数总数 **{len({r['key'] for r in mutlib.collect()})}**",
        f"- 判无效果（注入 `return None` 后逐行预测不变）: **{len(cand)}**",
        f"  - A 从未执行: {sum(1 for c in cand if c['cls'].startswith('A'))}",
        f"  - B 执行但无效果: {sum(1 for c in cand if c['cls'].startswith('B'))}",
        f"- 其中 **零下游引用（可直接清理）: {len(clean)}**；被下游引用（需同步改）: {len(linked)}",
        "",
        "## 一、可直接清理（零下游引用）",
        block(clean),
        "## 二、需同步修改下游才能删",
        "",
        "下列函数在 `tests/` / `server/` / `benchmarks/` 里被引用，删除会 fail 测试或打断服务端，"
        "需连同调用方一起改。",
        block(linked),
        "## 三、注意",
        "",
        f"- **导入期必需（别动）**: {len(essential)} 个 —— 注入后模块直接导入失败，"
        + ("`" + "`, `".join(essential) + "`" if essential else "（无）"),
        f"- **口径自洽性**: 从未执行却被判定为「有改变」的函数 = {len(inconsistent)} 个"
        + ("（应为 0；不为 0 说明覆盖率与变异口径打架，先查清）" if inconsistent else " ✓"),
        f"- **整模块未被导入**: {', '.join('`' + m + '`' for m in never_mod) or '（无）'}",
        "",
        "## 四、口径边界",
        "",
        f"- 这是「在 `{sweepcfg.BENCH_DATA.name}` 这份 benchmark 口径下没被用到」，"
        "**不等于绝对死代码**——benchmark 未覆盖的输入形态仍可能触达。",
        "- B 类不是「没跑到」，是「跑到了但结果没人用/被下游覆盖」，删之前建议确认那条"
        "分支真的不会让别的输入退化。",
        "",
        "## 五、复现",
        "",
        "```bash",
        f"cd {sweepcfg.ROOT}",
        "python sweep_base.py && python sweep_cov.py && python sweep_base.py   # 基线与探针序",
        "python sweep_drive.py                                                 # 全量扫描",
        "python sweep_verify.py                                                # 批量注入终局确认",
        "python sweep_analyze.py                                               # 本报告",
        "```",
        "",
    ]
    REPORT.write_text("\n".join(md), encoding="utf-8")
    sweepcfg.ANALYSIS.write_text(json.dumps(cand, ensure_ascii=False, indent=1), encoding="utf-8")
    print(f"无效果 {len(cand)} 个（A {sum(1 for c in cand if c['cls'].startswith('A'))} / "
          f"B {sum(1 for c in cand if c['cls'].startswith('B'))}）")
    print(f"零引用可删 {len(clean)} | 需同步改下游 {len(linked)} | 导入期必需 {len(essential)} | "
          f"口径矛盾 {len(inconsistent)}")
    print(f"报告 -> {REPORT}")


if __name__ == "__main__":
    main()
