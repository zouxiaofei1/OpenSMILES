"""Build tools/benchmark_pred_preview.html from predicted benchmark rows.

Columns: SMILES | predicted IUPAC | gold name | structure (SmilesDrawer).
ok/match = prediction matches gold (dual when zh gold exists), not merely non-empty return.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TOOLS = ROOT / "tools"
DATA = TOOLS / "benchmark_pred_preview_data.json"
DRAWER_JS = TOOLS / "smiles-drawer.min.js"
OUT = TOOLS / "benchmark_pred_preview.html"
_WS = re.compile(r"\s+")


def normalize_en(name: str) -> str:
    s = (name or "").strip().lower()
    s = s.replace("–", "-").replace("—", "-")
    return _WS.sub(" ", s).replace(" ,", ",")


def normalize_zh(name: str) -> str:
    return (name or "").strip()


def score_row(r: dict) -> dict:
    pred_en = r.get("pred_en") or ""
    pred_zh = r.get("pred_zh") or ""
    gold_en = r.get("gold_en") or ""
    gold_zh = r.get("gold_zh") or ""
    en_ok = bool(gold_en) and normalize_en(pred_en) == normalize_en(gold_en)
    zh_ok = bool(gold_zh) and normalize_zh(pred_zh) == normalize_zh(gold_zh)
    if gold_en and gold_zh:
        dual = en_ok and zh_ok
    elif gold_en:
        dual = en_ok
    elif gold_zh:
        dual = zh_ok
    else:
        dual = False
    return {
        "s": r.get("smiles") or "",
        "en": pred_en,
        "zh": pred_zh,
        "ge": gold_en,
        "gz": gold_zh,
        "en_ok": en_ok,
        "zh_ok": zh_ok,
        "ok": dual,
        "ret": bool(r.get("ok")),
    }


def main() -> None:
    rows = json.loads(DATA.read_text(encoding="utf-8"))
    payload = [score_row(r) for r in rows]
    n_ok = sum(1 for r in payload if r["ok"])
    n_en = sum(1 for r in payload if r["en_ok"])
    data_js = json.dumps(payload, ensure_ascii=False, separators=(",", ":"))
    drawer = DRAWER_JS.read_text(encoding="utf-8") if DRAWER_JS.exists() else ""
    if not drawer:
        raise SystemExit(f"missing {DRAWER_JS}")
    html = (
        _HTML_HEAD
        + data_js
        + ";\n</script>\n<script>\n"
        + drawer
        + "\n</script>\n"
        + _HTML_TAIL
    )
    OUT.write_text(html, encoding="utf-8")
    print(
        f"wrote {OUT} ({OUT.stat().st_size / 1024:.1f} KB, rows={len(payload)}, "
        f"dual_ok={n_ok}, en_ok={n_en})"
    )


_HTML_HEAD = r"""<!DOCTYPE html>
<html lang="zh-CN">
<head>
<meta charset="UTF-8" />
<meta name="viewport" content="width=device-width, initial-scale=1" />
<title>NamePredict 测试数据预览</title>
<style>
  :root {
    --bg: #0f1419;
    --panel: #1a2332;
    --border: #2d3a4f;
    --text: #e7ecf3;
    --muted: #8b9bb4;
    --accent: #5b9fd4;
    --row-hover: #243044;
    --good: #7dcea0;
    --bad: #e89a9a;
  }
  * { box-sizing: border-box; }
  body {
    margin: 0;
    font-family: "Segoe UI", "PingFang SC", "Microsoft YaHei", sans-serif;
    background: var(--bg);
    color: var(--text);
    line-height: 1.45;
  }
  header {
    position: sticky; top: 0; z-index: 20;
    background: rgba(15,20,25,.96);
    border-bottom: 1px solid var(--border);
    padding: 12px 20px;
    display: flex; flex-wrap: wrap; gap: 12px; align-items: center;
  }
  h1 { margin: 0; font-size: 1.05rem; font-weight: 600; }
  .meta { color: var(--muted); font-size: .85rem; }
  .controls {
    display: flex; flex-wrap: wrap; gap: 8px; align-items: center;
    margin-left: auto;
  }
  input[type="search"], select, button {
    background: var(--panel); color: var(--text);
    border: 1px solid var(--border); border-radius: 6px;
    padding: 6px 10px; font-size: .9rem;
  }
  input[type="search"] { min-width: 200px; }
  button { cursor: pointer; }
  button:disabled { opacity: .45; cursor: not-allowed; }
  .chip {
    font-size: .78rem; padding: 2px 8px; border-radius: 999px;
    border: 1px solid var(--border); color: var(--muted);
  }
  .chip b { color: var(--text); font-weight: 600; }
  main { padding: 0 12px 40px; }
  table {
    width: 100%; border-collapse: collapse; table-layout: fixed;
    font-size: .88rem;
  }
  thead th {
    position: sticky; top: 62px; z-index: 10;
    background: #152032; border-bottom: 1px solid var(--border);
    text-align: left; padding: 10px 12px; font-weight: 600;
    color: var(--muted); font-size: .75rem; letter-spacing: .04em;
    text-transform: uppercase;
  }
  tbody td {
    border-bottom: 1px solid var(--border);
    padding: 10px 12px; vertical-align: middle;
    word-break: break-word;
  }
  tbody tr:hover { background: var(--row-hover); }
  tbody tr.match { background: rgba(61,154,106,.06); }
  tbody tr.miss { background: rgba(196,92,92,.04); }
  col.smiles { width: 22%; }
  col.pred { width: 28%; }
  col.gold { width: 28%; }
  col.struct { width: 22%; }
  .smiles {
    font-family: ui-monospace, "Cascadia Code", Consolas, monospace;
    font-size: .76rem; color: #c5d4e8; word-break: break-all;
  }
  .idx {
    display: inline-block; min-width: 2.6em; color: var(--muted);
    font-size: .75rem; margin-right: 6px;
  }
  .name-en { font-weight: 500; }
  .name-zh { color: var(--muted); font-size: .86rem; margin-top: 2px; }
  .gold-en { font-weight: 500; color: #c8d8ef; }
  .gold-zh { color: #9aacc4; font-size: .86rem; margin-top: 2px; }
  .badge {
    display: inline-block; font-size: .7rem; padding: 1px 6px;
    border-radius: 4px; margin-left: 6px; vertical-align: middle;
  }
  .badge.ok { background: rgba(61,154,106,.22); color: var(--good); }
  .badge.fail { background: rgba(196,92,92,.22); color: var(--bad); }
  .badge.ret { background: rgba(91,159,212,.18); color: #9ec9ea; }
  .struct-wrap {
    display: flex; justify-content: center; align-items: center;
    min-height: 150px; background: #ffffff; border-radius: 8px;
    border: 1px solid var(--border); padding: 4px;
  }
  .struct-wrap svg {
    width: 220px; height: 140px; max-width: 100%; display: block;
    background: #ffffff;
  }
  .struct-err { color: var(--bad); font-size: .8rem; }
  .empty { text-align: center; color: var(--muted); padding: 48px 16px; }
  footer {
    color: var(--muted); font-size: .78rem; text-align: center;
    padding: 16px; border-top: 1px solid var(--border);
  }
</style>
</head>
<body>
<header>
  <div>
    <h1>NamePredict · 测试数据预览</h1>
    <div class="meta">
      来源 <code>data/merged_benchmark.json</code> ·
      四列：SMILES / 预测IUPAC / 正确命名 / 结构 ·
      <b>match</b>=与金标一致（不是“有返回值”）
    </div>
  </div>
  <div class="controls">
    <span class="chip">共 <b id="total">0</b></span>
    <span class="chip">显示 <b id="shown">0</b></span>
    <span class="chip">match <b id="okn">0</b></span>
    <span class="chip">en✓ <b id="enn">0</b></span>
    <input id="q" type="search" placeholder="筛选 SMILES / 预测 / 正确命名…" />
    <select id="filter">
      <option value="all">全部</option>
      <option value="ok">仅 match</option>
      <option value="fail">仅 miss</option>
      <option value="wrong">有返回但错</option>
    </select>
    <select id="pageSize">
      <option value="50">每页 50</option>
      <option value="100" selected>每页 100</option>
      <option value="200">每页 200</option>
      <option value="500">每页 500</option>
    </select>
    <button id="prev" type="button">上一页</button>
    <span class="chip">页 <b id="page">1</b> / <b id="pages">1</b></span>
    <button id="next" type="button">下一页</button>
  </div>
</header>
<main>
  <table>
    <colgroup>
      <col class="smiles" />
      <col class="pred" />
      <col class="gold" />
      <col class="struct" />
    </colgroup>
    <thead>
      <tr>
        <th>SMILES</th>
        <th>本程序预测 IUPAC</th>
        <th>正确命名</th>
        <th>结构</th>
      </tr>
    </thead>
    <tbody id="tbody"></tbody>
  </table>
  <div id="empty" class="empty" hidden>无匹配行</div>
</main>
<footer>
  正确命名 = benchmark 金标 english_name / chinese_name。
  match = 预测与金标一致。结构由内嵌 SmilesDrawer 渲染 canvas[data-smiles]。
</footer>
<script>
window.__ROWS__ = """

_HTML_TAIL = r"""
<script>
(function () {
  const ROWS = window.__ROWS__ || [];
  const tbody = document.getElementById("tbody");
  const empty = document.getElementById("empty");
  const qEl = document.getElementById("q");
  const filterEl = document.getElementById("filter");
  const pageSizeEl = document.getElementById("pageSize");
  const prevBtn = document.getElementById("prev");
  const nextBtn = document.getElementById("next");
  const totalEl = document.getElementById("total");
  const shownEl = document.getElementById("shown");
  const oknEl = document.getElementById("okn");
  const ennEl = document.getElementById("enn");
  const pageEl = document.getElementById("page");
  const pagesEl = document.getElementById("pages");
  let page = 1;

  function filteredFast() {
    const q = (qEl.value || "").trim().toLowerCase();
    const f = filterEl.value;
    const out = [];
    for (let i = 0; i < ROWS.length; i++) {
      const r = ROWS[i];
      if (f === "ok" && !r.ok) continue;
      if (f === "fail" && r.ok) continue;
      if (f === "wrong" && !(r.ret && !r.ok)) continue;
      if (q) {
        const hay = [r.s, r.en, r.zh, r.ge, r.gz].join(" ").toLowerCase();
        if (hay.indexOf(q) < 0) continue;
      }
      out.push({ r: r, abs: i });
    }
    return out;
  }

  function badge(r) {
    if (r.ok) return ' <span class="badge ok">match</span>';
    if (r.ret) return ' <span class="badge fail">miss</span><span class="badge ret">returned</span>';
    return ' <span class="badge fail">miss</span>';
  }

  // SvgDrawer draws into <svg> directly (canvas path uses SVG→canvas and often stays blank).
  function drawAll() {
    if (typeof SmilesDrawer === "undefined" || !SmilesDrawer.SvgDrawer) {
      tbody.querySelectorAll(".struct-wrap").forEach(function (w) {
        w.innerHTML = '<span class="struct-err">SmilesDrawer 未加载</span>';
      });
      return;
    }
    const opts = { width: 220, height: 140, bondThickness: 1.1, padding: 6 };
    const drawer = new SmilesDrawer.SvgDrawer(opts);
    const svgs = tbody.querySelectorAll("svg[data-smiles]");
    for (let i = 0; i < svgs.length; i++) {
      (function (svg) {
        const smi = svg.getAttribute("data-smiles") || "";
        if (!smi) {
          svg.parentNode.innerHTML = '<span class="struct-err">空 SMILES</span>';
          return;
        }
        SmilesDrawer.parse(
          smi,
          function (tree) {
            try {
              // SvgDrawer.draw(data, target, theme, weights=null, infoOnly=false)
              drawer.draw(tree, svg, "light");
            } catch (err) {
              svg.parentNode.innerHTML = '<span class="struct-err">draw err</span>';
              console.warn(smi, err);
            }
          },
          function (err) {
            svg.parentNode.innerHTML = '<span class="struct-err">parse fail</span>';
            console.warn(smi, err);
          }
        );
      })(svgs[i]);
    }
  }

  function render() {
    const items = filteredFast();
    const ps = Math.max(1, parseInt(pageSizeEl.value, 10) || 100);
    const pages = Math.max(1, Math.ceil(items.length / ps));
    if (page > pages) page = pages;
    if (page < 1) page = 1;
    const start = (page - 1) * ps;
    const slice = items.slice(start, start + ps);

    totalEl.textContent = String(ROWS.length);
    shownEl.textContent = String(items.length);
    oknEl.textContent = String(ROWS.filter(function (r) { return r.ok; }).length);
    ennEl.textContent = String(ROWS.filter(function (r) { return r.en_ok; }).length);
    pageEl.textContent = String(page);
    pagesEl.textContent = String(pages);
    prevBtn.disabled = page <= 1;
    nextBtn.disabled = page >= pages;

    tbody.innerHTML = "";
    empty.hidden = slice.length > 0;

    for (let k = 0; k < slice.length; k++) {
      const item = slice[k];
      const r = item.r;
      const abs = item.abs;
      const tr = document.createElement("tr");
      tr.className = r.ok ? "match" : "miss";

      // SMILES
      const tdS = document.createElement("td");
      tdS.innerHTML = '<span class="idx">#' + (abs + 1) + "</span>";
      const sm = document.createElement("span");
      sm.className = "smiles";
      sm.textContent = r.s || "";
      tdS.appendChild(sm);

      // predicted
      const tdP = document.createElement("td");
      const en = (r.en || "").trim() || "(空)";
      const zh = (r.zh || "").trim() || "—";
      tdP.innerHTML = '<div class="name-en"></div><div class="name-zh"></div>';
      const pe = tdP.querySelector(".name-en");
      pe.appendChild(document.createTextNode(en));
      pe.insertAdjacentHTML("beforeend", badge(r));
      tdP.querySelector(".name-zh").textContent = zh;

      // gold / correct name
      const tdG = document.createElement("td");
      const ge = (r.ge || "").trim() || "—";
      const gz = (r.gz || "").trim() || "—";
      tdG.innerHTML = '<div class="gold-en"></div><div class="gold-zh"></div>';
      tdG.querySelector(".gold-en").textContent = ge;
      tdG.querySelector(".gold-zh").textContent = gz;

      // structure: svg[data-smiles] for SmilesDrawer.SvgDrawer
      const tdM = document.createElement("td");
      const wrap = document.createElement("div");
      wrap.className = "struct-wrap";
      const svg = document.createElementNS("http://www.w3.org/2000/svg", "svg");
      svg.setAttribute("data-smiles", r.s || "");
      svg.setAttribute("width", "220");
      svg.setAttribute("height", "140");
      svg.setAttribute("xmlns", "http://www.w3.org/2000/svg");
      wrap.appendChild(svg);
      tdM.appendChild(wrap);

      tr.append(tdS, tdP, tdG, tdM);
      tbody.appendChild(tr);
    }

    requestAnimationFrame(function () { drawAll(); });
  }

  let timer = null;
  function schedule() {
    clearTimeout(timer);
    timer = setTimeout(function () { page = 1; render(); }, 120);
  }
  qEl.addEventListener("input", schedule);
  filterEl.addEventListener("change", function () { page = 1; render(); });
  pageSizeEl.addEventListener("change", function () { page = 1; render(); });
  prevBtn.addEventListener("click", function () { page -= 1; render(); });
  nextBtn.addEventListener("click", function () { page += 1; render(); });

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", render);
  } else {
    render();
  }
})();
</script>
</body>
</html>
"""


if __name__ == "__main__":
    main()
