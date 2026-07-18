"""Build tools/benchmark_pred_preview_claimable3000.html for first 3000 gold rows.

Uses current worktree SMILESNNamer (universal claimable-block path).
Columns: SMILES | predicted IUPAC | gold name | structure (SmilesDrawer).
"""
from __future__ import annotations

import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TOOLS = ROOT / "tools"
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

try:
    from rdkit import RDLogger

    RDLogger.DisableLog("rdApp.*")
except Exception:
    pass

from namepredict.constants import normalize_en, normalize_zh  # noqa: E402

DATA = ROOT / "data" / "merged_benchmark.json"
if not DATA.is_file():
    DATA = Path(r"E:/dev/chem/data/merged_benchmark.json")
DRAWER_JS = TOOLS / "smiles-drawer.min.js"
OUT = TOOLS / "benchmark_pred_preview_claimable3000.html"
OUT_JSON = TOOLS / "benchmark_pred_preview_claimable3000_data.json"
LIMIT = 3000
_PROGRESS_EVERY = 100
_WORKER_NAMER = None


def _default_workers() -> int:
    n = os.cpu_count() or 1
    return max(1, n - 1) if n > 1 else 1


def _init_worker() -> None:
    global _WORKER_NAMER
    try:
        from rdkit import RDLogger

        RDLogger.DisableLog("rdApp.*")
    except Exception:
        pass
    if str(SRC) not in sys.path:
        sys.path.insert(0, str(SRC))
    from namepredict.namer import SMILESNNamer

    _WORKER_NAMER = SMILESNNamer()


def _score_pred(pred_en: str, pred_zh: str, gold_en: str, gold_zh: str) -> dict:
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
        "en_ok": en_ok,
        "zh_ok": zh_ok,
        "ok": dual,
        "ret": bool(pred_en or pred_zh),
    }


def _score_one(item: tuple[int, dict]) -> tuple[int, dict]:
    """Worker: (index, gold_row) -> (index, preview_row)."""
    idx, row = item
    smi = str(row.get("smiles") or "")
    gold_en = row.get("english_name") or ""
    gold_zh = row.get("chinese_name") or ""
    try:
        r = _WORKER_NAMER.name(smi)
        pred_en, pred_zh = r.en or "", r.zh or ""
    except Exception:
        pred_en, pred_zh = "", ""
    sc = _score_pred(pred_en, pred_zh, gold_en, gold_zh)
    return idx, {
        "s": smi,
        "en": pred_en,
        "zh": pred_zh,
        "ge": gold_en,
        "gz": gold_zh,
        **sc,
    }


def _filled_payload(out: list[dict | None]) -> list[dict]:
    """Completed rows in original index order (skips unfinished slots)."""
    return [x for x in out if x is not None]


def _print_progress(done: int, total: int, t0: float, out: list[dict | None]) -> None:
    elapsed = time.perf_counter() - t0
    rate = done / elapsed if elapsed > 0 else 0.0
    filled = _filled_payload(out)
    ok_en = sum(1 for x in filled if x["en_ok"])
    ok_zh = sum(1 for x in filled if x["zh_ok"])
    ok_d = sum(1 for x in filled if x["ok"])
    n = max(1, len(filled))
    print(
        f"[{done}/{total}] {rate:.1f} row/s  "
        f"EN={ok_en}/{n} ({100 * ok_en / n:.1f}%)  "
        f"ZH={ok_zh}/{n} ({100 * ok_zh / n:.1f}%)  "
        f"dual={ok_d}/{n} ({100 * ok_d / n:.1f}%)",
        flush=True,
    )


def _write_preview(payload: list[dict], *, total_target: int, drawer: str) -> None:
    """Rewrite HTML + JSON from current completed rows (safe to open mid-run)."""
    n = len(payload)
    n_ok = sum(1 for r in payload if r["ok"])
    n_en = sum(1 for r in payload if r["en_ok"])
    n_zh = sum(1 for r in payload if r["zh_ok"])
    OUT_JSON.write_text(json.dumps(payload, ensure_ascii=False), encoding="utf-8")
    data_js = json.dumps(payload, ensure_ascii=False, separators=(",", ":"))
    # Inject live progress into the sticky header meta line.
    head = _HTML_HEAD.replace(
        "前 3000 条",
        f"已完成 {n}/{total_target} 条"
        f"（dual {n_ok}/{max(1, n)} = {100 * n_ok / max(1, n):.1f}%）",
    )
    html = (
        head
        + data_js
        + ";\n</script>\n<script>\n"
        + drawer
        + "\n</script>\n"
        + _HTML_TAIL
    )
    # Atomic-ish replace on Windows: write temp then replace.
    tmp = OUT.with_suffix(".html.tmp")
    tmp.write_text(html, encoding="utf-8")
    tmp.replace(OUT)
    print(
        f"  wrote {OUT.name} rows={n}/{total_target}  "
        f"dual={n_ok} ({100 * n_ok / max(1, n):.1f}%)  "
        f"en={n_en} ({100 * n_en / max(1, n):.1f}%)  "
        f"zh={n_zh} ({100 * n_zh / max(1, n):.1f}%)  "
        f"size={OUT.stat().st_size / 1024:.1f} KB",
        flush=True,
    )


def _checkpoint(done: int, total: int, t0: float, out: list[dict | None], drawer: str) -> None:
    _print_progress(done, total, t0, out)
    _write_preview(_filled_payload(out), total_target=total, drawer=drawer)


def _predict_rows(
    rows: list[dict],
    *,
    drawer: str,
    workers: int | None = None,
) -> list[dict]:
    n = len(rows)
    n_workers = max(1, workers if workers is not None else _default_workers())
    t0 = time.perf_counter()
    print(f"predict with workers={n_workers}", flush=True)
    out: list[dict | None] = [None] * n
    items = list(enumerate(rows))

    if n_workers == 1 or n <= 1:
        _init_worker()
        for i, item in enumerate(items, 1):
            idx, row_out = _score_one(item)
            out[idx] = row_out
            if i % _PROGRESS_EVERY == 0 or i == n:
                _checkpoint(i, n, t0, out, drawer)
        return _filled_payload(out)

    with ProcessPoolExecutor(max_workers=n_workers, initializer=_init_worker) as pool:
        futures = [pool.submit(_score_one, item) for item in items]
        done = 0
        for fut in as_completed(futures):
            idx, row_out = fut.result()
            out[idx] = row_out
            done += 1
            if done % _PROGRESS_EVERY == 0 or done == n:
                _checkpoint(done, n, t0, out, drawer)
    return _filled_payload(out)


def main() -> None:
    if not DATA.is_file():
        raise SystemExit(f"missing benchmark data: {DATA}")
    if not DRAWER_JS.is_file():
        raise SystemExit(f"missing {DRAWER_JS}")
    rows = json.loads(DATA.read_text(encoding="utf-8"))
    if not isinstance(rows, list):
        raise SystemExit("expected JSON list")
    rows = rows[:LIMIT]
    workers = int(os.environ.get("PREVIEW_WORKERS", "0")) or None
    drawer = DRAWER_JS.read_text(encoding="utf-8")
    print(f"predict first {len(rows)} rows from {DATA}", flush=True)
    print(f"HTML checkpoint every {_PROGRESS_EVERY} results -> {OUT}", flush=True)
    payload = _predict_rows(rows, drawer=drawer, workers=workers)
    # Final write already happened at done==n; rewrite once more for clean header.
    _write_preview(payload, total_target=len(rows), drawer=drawer)
    print(f"done. open {OUT}", flush=True)


_HTML_HEAD = r"""<!DOCTYPE html>
<html lang="zh-CN">
<head>
<meta charset="UTF-8" />
<meta name="viewport" content="width=device-width, initial-scale=1" />
<title>NamePredict · Claimable-Block 预览 (前3000)</title>
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
    <h1>NamePredict · Claimable-Block 预览（前 3000）</h1>
    <div class="meta">
      分支 universal-claimable-block · 前 3000 条
      <code>data/merged_benchmark.json</code> ·
      四列：SMILES / 预测IUPAC / 正确命名 / 结构 ·
      <b>match</b>=与金标一致（不是“有返回值”）
    </div>
  </div>
  <div class="controls">
    <span class="chip">共 <b id="total">0</b></span>
    <span class="chip">显示 <b id="shown">0</b></span>
    <span class="chip">match <b id="okn">0</b></span>
    <span class="chip">en✓ <b id="enn">0</b></span>
    <span class="chip">zh✓ <b id="zhn">0</b></span>
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
  Claimable-block 路径：owned_atoms + CoverageLedger 门控 + 候选重试。
  正确命名 = benchmark 金标 english_name / chinese_name。
  match = 预测与金标一致。结构由内嵌 SmilesDrawer 渲染。
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
  const zhnEl = document.getElementById("zhn");
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

  function esc(s) {
    return String(s || "")
      .replace(/&/g, "&amp;")
      .replace(/</g, "&lt;")
      .replace(/>/g, "&gt;")
      .replace(/"/g, "&quot;");
  }

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
              drawer.draw(tree, svg, "light");
            } catch (err) {
              svg.parentNode.innerHTML = '<span class="struct-err">draw err</span>';
            }
          },
          function () {
            svg.parentNode.innerHTML = '<span class="struct-err">parse fail</span>';
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
    zhnEl.textContent = String(ROWS.filter(function (r) { return r.zh_ok; }).length);
    pageEl.textContent = String(page);
    pagesEl.textContent = String(pages);
    prevBtn.disabled = page <= 1;
    nextBtn.disabled = page >= pages;

    if (!slice.length) {
      tbody.innerHTML = "";
      empty.hidden = false;
      return;
    }
    empty.hidden = true;
    let html = "";
    for (let j = 0; j < slice.length; j++) {
      const abs = slice[j].abs;
      const r = slice[j].r;
      const cls = r.ok ? "match" : "miss";
      const predEn = esc(r.en) || "<span class='name-zh'>(空)</span>";
      const predZh = esc(r.zh);
      const goldEn = esc(r.ge) || "<span class='name-zh'>(无英标)</span>";
      const goldZh = esc(r.gz);
      html +=
        '<tr class="' + cls + '">' +
        '<td><span class="idx">#' + (abs + 1) + '</span><span class="smiles">' + esc(r.s) + "</span></td>" +
        '<td><div class="name-en">' + predEn + badge(r) + "</div>" +
        (predZh ? '<div class="name-zh">' + predZh + "</div>" : "") + "</td>" +
        '<td><div class="gold-en">' + goldEn + "</div>" +
        (goldZh ? '<div class="gold-zh">' + goldZh + "</div>" : "") + "</td>" +
        '<td><div class="struct-wrap"><svg data-smiles="' + esc(r.s) + '"></svg></div></td>' +
        "</tr>";
    }
    tbody.innerHTML = html;
    drawAll();
  }

  qEl.addEventListener("input", function () { page = 1; render(); });
  filterEl.addEventListener("change", function () { page = 1; render(); });
  pageSizeEl.addEventListener("change", function () { page = 1; render(); });
  prevBtn.addEventListener("click", function () { page -= 1; render(); });
  nextBtn.addEventListener("click", function () { page += 1; render(); });
  render();
})();
</script>
</body>
</html>
"""


if __name__ == "__main__":
    main()
