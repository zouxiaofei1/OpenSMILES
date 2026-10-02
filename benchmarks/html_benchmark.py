"""单文件前端 index.min.html 跑分：用 merged_benchmark 给内嵌 JS 命名管线打分。

页面把 JS 版命名管线（模块表 + RDKit wasm）以 base85+brotli 内联在 <script> 里，
本脚本解出这些载荷、交给 Node 逐行命名，再用 benchmark.score_record 评分。
口径与 benchmarks.benchmark_parallel 一致：每行新建 namer（等价每行清空 cache）。

用法:
  python -m benchmarks.html_benchmark --data benchmarks/merged_benchmark.json --time
  python -m benchmarks.html_benchmark --html tmp/index.min.html --workers 8 --json
  python -m benchmarks.html_benchmark --baseline benchmarks/.last_run.json
"""

from __future__ import annotations
import argparse
import json
import shutil
import subprocess
import sys
import tempfile
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any

_ROOT = Path(__file__).resolve().parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from benchmarks.benchmark import (  # noqa: E402
    _empty_bucket,
    _finalize,
    _is_fail,
    _load_rows,
    _print_summary,
    _row_key,
    _tally,
    score_record,
)

_DEFAULT_DATA = _ROOT / "benchmarks" / "merged_benchmark.json"
_HTML_CANDIDATES = (
    _ROOT / "tmp" / "index.min.html",
    _ROOT / "tools" / "index.min.html",
)

# 页面内联 JS 管线的 Node 侧装载器：解载荷 → 起 RDKit → 逐行命名 → 写预测。
_RUNNER_JS = r"""
// 从 index.min.html 解出模块表与 RDKit wasm，命名一行数据集分片。
import fs from "node:fs";
import path from "node:path";
import zlib from "node:zlib";
import { createRequire } from "node:module";
import { fileURLToPath } from "node:url";

function opt(name, def) {
  const i = process.argv.indexOf("--" + name);
  return i >= 0 ? process.argv[i + 1] : def;
}
const HTML = opt("html"), DATA = opt("data"), OUT = opt("out");
const SHARD = parseInt(opt("shard", "0"), 10);
const SHARDS = parseInt(opt("shards", "1"), 10);
const LIMIT = opt("limit", "") === "" ? null : parseInt(opt("limit"), 10);

const X85A = "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ!#$%&()*+,-.:;=?@^_[]{}";
const X85R = new Int16Array(128);
for (let i = 0; i < 85; i++) X85R[X85A.charCodeAt(i)] = i;

// 自定义 base85 文本 -> 字节（末组按 ascii85 规则补 'u' 位）
function x85(s) {
  s = s.replace(/\s+/g, "");
  const n = s.length, full = (n / 5) | 0, rem = n - full * 5;
  const out = new Uint8Array(full * 4 + (rem > 0 ? rem - 1 : 0));
  let o = 0, i = 0;
  for (let g = 0; g < full; g++) {
    let v = 0;
    for (let k = 0; k < 5; k++) v = (v * 85 + X85R[s.charCodeAt(i++)]) >>> 0;
    out[o++] = v >>> 24; out[o++] = (v >>> 16) & 255;
    out[o++] = (v >>> 8) & 255; out[o++] = v & 255;
  }
  if (rem > 0) {
    let w = 0;
    for (let j = 0; j < 5; j++) {
      const d = j < rem ? X85R[s.charCodeAt(i++)] : 84;
      w = (w * 85 + d) >>> 0;
    }
    for (let m = 0; m < rem - 1; m++) out[o++] = (w >>> (24 - 8 * m)) & 255;
  }
  return out;
}

const html = fs.readFileSync(HTML, "utf8");
function grab(name) {
  const m = html.match(new RegExp(name + "\\s*=\\s*\"([^\"]*)\""));
  if (!m) throw new Error("页面缺少 " + name);
  return m[1];
}
const MODE = (html.match(/const MODE\s*=\s*"(\w+)"/) || [])[1] || "min";
// 两种载荷入口：min 版走内联 brotli，plain 版原样返回
function decode(payload) {
  const bytes = Buffer.from(x85(payload));
  return MODE === "min" ? zlib.brotliDecompressSync(bytes) : bytes;
}
const MODS = opt("mods", "");
// 模块来源：默认取页面内联载荷；--mods 给目录则直接读该目录，便于改源码时免重打包
const MODULES = MODS
  ? Object.fromEntries(fs.readdirSync(MODS).filter((f) => f.endsWith(".js"))
      .map((f) => [f, fs.readFileSync(path.join(MODS, f), "utf8")]))
  : JSON.parse(decode(grab("MODULES_PAYLOAD")).toString("utf8"));
const wasmBytes = decode(grab("WASM_PAYLOAD"));

// RDKit 胶水脚本内联在页面里，取出后用 new Function 求值成 initRDKitModule
const gs = html.indexOf("var initRDKitModule=");
if (gs < 0) throw new Error("页面缺少 RDKit 胶水脚本");
const glue = html.slice(gs, html.indexOf("</script>", gs));
const here = fileURLToPath(import.meta.url);
const initRDKitModule = new Function(
  "require", "module", "exports", "__dirname", "__filename",
  glue + "\n;return initRDKitModule;"
)(createRequire(import.meta.url), { exports: {} }, {}, path.dirname(here), here);

const RDKit = await initRDKitModule({
  instantiateWasm(imports, onSuccess) {
    WebAssembly.instantiate(wasmBytes, imports).then((r) => onSuccess(r.instance));
  },
});

// 与页面一致的 CommonJS 装载器：说明符按文件名平铺
const _cache = {};
function req(spec) {
  const base = spec.split("/").pop();
  const name = base.endsWith(".js") ? base : base + ".js";
  if (!Object.prototype.hasOwnProperty.call(MODULES, name)) throw new Error("模块缺失 " + name);
  if (_cache[name]) return _cache[name].exports;
  const m = { exports: {} };
  _cache[name] = m;
  new Function("module", "exports", "require", MODULES[name])(m, m.exports, req);
  return m.exports;
}
const entry = req("namer.js");

let rows = JSON.parse(fs.readFileSync(DATA, "utf8"));
if (LIMIT !== null) rows = rows.slice(0, LIMIT);

const preds = [];
let done = 0;
for (let i = SHARD; i < rows.length; i += SHARDS) {
  const row = rows[i];
  let en = "", zh = "", success = false, err = "";
  try {
    const namer = new entry.SMILESNNamer({ rdkit: RDKit, cache: null });
    const r = namer.name(String(row.smiles || ""));
    en = r.en || "";
    zh = r.zh || "";
    success = !!r.success;
  } catch (e) {
    err = (e && e.message) || String(e);
  }
  preds.push({ i, en, zh, success, err });
  done++;
  if (done % 100 === 0) console.log("PROGRESS " + done);
}
fs.writeFileSync(OUT, JSON.stringify(preds), "utf8");
console.log("PROGRESS " + done);
"""


def _find_html(explicit: str | None) -> Path:
    """定位单文件前端：显式路径优先，否则 tmp/、tools/ 依次查找。"""
    if explicit:
        p = Path(explicit)
        if not p.is_file():
            raise FileNotFoundError(f"找不到页面: {p}")
        return p
    for cand in _HTML_CANDIDATES:
        if cand.is_file():
            return cand
    raise FileNotFoundError("找不到 index.min.html（tmp/ 或 tools/）")


def _node_exe() -> str:
    """取 node 可执行文件路径。"""
    exe = shutil.which("node")
    if not exe:
        raise RuntimeError("跑 JS 管线需要 node，但 PATH 里没有")
    return exe


def _progress_bar(done: int, total: int, width: int = 28) -> str:
    """文本进度条。"""
    frac = 1.0 if total <= 0 else min(1.0, done / total)
    filled = int(width * frac)
    return f"[{'#' * filled}{'-' * (width - filled)}] {done}/{total} ({100.0 * frac:.1f}%)"


def _run_shard(
    node: str, runner: Path, html: Path, data: Path, out: Path,
    shard: int, shards: int, limit: int | None, on_progress, mods: str | None = None,
) -> float:
    """跑一个分片，边读 PROGRESS 行边回报；返回耗时秒。"""
    args = [node, str(runner), "--html", str(html), "--data", str(data),
            "--out", str(out), "--shard", str(shard), "--shards", str(shards)]
    if limit is not None:
        args += ["--limit", str(limit)]
    if mods:
        args += ["--mods", str(mods)]
    t0 = time.perf_counter()
    proc = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            text=True, encoding="utf-8", bufsize=1)
    for line in proc.stdout:
        if line.startswith("PROGRESS "):
            on_progress(int(line.split()[1]))
    err = proc.stderr.read()
    if proc.wait() != 0:
        raise RuntimeError(f"分片 {shard} 失败: {err.strip()[:400]}")
    return time.perf_counter() - t0


def run_benchmark_html(
    html_path: Path,
    data_path: Path,
    limit: int | None = None,
    workers: int = 1,
    quiet: bool = False,
    mods: str | None = None,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """跑全量行并评分；返回 (report, 合并后的逐行预测)。"""
    rows = _load_rows(data_path, limit)
    n_rows = len(rows)
    workers = max(1, min(workers, max(1, n_rows)))
    node = _node_exe()
    total = _empty_bucket()
    preds: dict[int, dict[str, Any]] = {}
    with tempfile.TemporaryDirectory(prefix="html_bm_") as tmp:
        tmpdir = Path(tmp)
        runner = tmpdir / "runner.mjs"
        runner.write_text(_RUNNER_JS, encoding="utf-8")
        outs = [tmpdir / f"preds_{i}.json" for i in range(workers)]
        lock = threading.Lock()
        t0 = time.perf_counter()
        counter = {"done": 0}

        def on_progress(k: int) -> None:
            if quiet:
                return
            with lock:
                counter["done"] += 1
                if counter["done"] % 100 == 0 or counter["done"] == n_rows:
                    el = time.perf_counter() - t0
                    rate = counter["done"] / el if el > 0 else 0.0
                    print(f"progress {_progress_bar(counter['done'], n_rows)}  "
                          f"{rate:.1f} row/s  elapsed={el:.1f}s", flush=True)

        with ThreadPoolExecutor(max_workers=workers) as pool:
            futs = [pool.submit(_run_shard, node, runner, html_path, data_path,
                                outs[i], i, workers, limit, on_progress, mods)
                    for i in range(workers)]
            for fut in as_completed(futs):
                fut.result()
        for out in outs:
            for p in json.loads(out.read_text(encoding="utf-8")):
                preds[p["i"]] = p

    ordered = [preds[i] for i in range(n_rows)]
    by_source: dict[str, dict[str, int]] = {}
    by_tier: dict[Any, dict[str, int]] = {}
    fails: list[dict] = []
    entries: list[dict] = []
    for row, p in zip(rows, ordered):
        score = score_record(p["en"], p["zh"], row)
        _tally(total, score)
        _tally(by_source.setdefault(str(row.get("source") or "unknown"), _empty_bucket()), score)
        _tally(by_tier.setdefault(row.get("tier", 0), _empty_bucket()), score)
        entry = {
            "key": _row_key(row), "id": row.get("id"), "smiles": row.get("smiles"),
            "source": row.get("source"), "tier": row.get("tier"),
            "english_name": row.get("english_name"), "chinese_name": row.get("chinese_name"),
            "pred_en": p["en"], "pred_zh": p["zh"], "js_err": p.get("err", ""),
            "en_ok": score["en_ok"], "zh_ok": score["zh_ok"], "dual_ok": score["dual_ok"],
        }
        entries.append(entry)
        if _is_fail(score):
            fails.append(entry)
    report = _finalize(total, by_source, by_tier, fails, entries)
    report["elapsed_sec"] = round(time.perf_counter() - t0, 3)
    report["workers"] = workers
    report["html"] = str(html_path)
    report["n_rows"] = n_rows
    return report, ordered


def _print_vs_baseline(report: dict[str, Any], baseline: Path, entries: list[dict]) -> None:
    """与 src 上次跑分 snapshot 对比总分并列出差异行。"""
    snap = json.loads(baseline.read_text(encoding="utf-8"))
    items = snap.get("items") or {}
    print(f"\n=== 对比基线 {baseline.name} ===")
    print(f"src  dual={snap.get('ok_dual')}/{snap.get('n_dual')}  "
          f"html dual={report['ok_dual']}/{report['n_dual']}")
    delta = int(report["ok_dual"]) - int(snap.get("ok_dual") or 0)
    print(f"差值 {delta:+d} 行（src 对而 html 错 {sum(1 for e in entries if e['dual_ok'] is False and (items.get(e['key']) or {}).get('dual_ok'))} 行）")


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Index.min.html JS 命名管线跑分")
    p.add_argument("--html", default=None, help="index.min.html 路径（默认 tmp/ 或 tools/）")
    p.add_argument("--data", type=Path, default=_DEFAULT_DATA, help="基准 JSON")
    p.add_argument("--limit", type=int, default=None, help="只跑前 N 行")
    p.add_argument("--workers", type=int, default=None, help="node 分片数（默认 min(cpu-1, 8)）")
    p.add_argument("--out-preds", type=Path, default=None, help="把逐行预测写到该文件")
    p.add_argument("--baseline", type=Path, default=None, help="src 跑分 snapshot，用于对比")
    p.add_argument("--mods", default=None, help="直接读该目录的 .js 模块（改源码时免重打包）")
    p.add_argument("--json", action="store_true", help="输出完整 report JSON")
    p.add_argument("--time", action="store_true", help="在 stderr 打印耗时")
    return p


def main(argv: list[str] | None = None) -> None:
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass
    args = _build_parser().parse_args(argv)
    cpu = __import__("os").cpu_count() or 1
    workers = args.workers if args.workers is not None else max(1, min(cpu - 1, 16))
    try:
        html_path = _find_html(args.html)
    except FileNotFoundError as exc:
        print(f"error: {exc}", file=sys.stderr)
        sys.exit(2)
    t0 = time.perf_counter()
    print(f"html benchmark start: html={html_path.name} data={args.data.name} workers={workers}", flush=True)
    report, preds = run_benchmark_html(html_path, args.data, args.limit, workers, mods=args.mods)
    elapsed = time.perf_counter() - t0
    n_rows = report["n_rows"]
    report["avg_ms_per_row"] = round(elapsed * 1000.0 / n_rows, 3) if n_rows else 0.0
    if args.out_preds:
        args.out_preds.write_text(
            json.dumps([{**p, "id": report["results"][p["i"]]["id"]} for p in preds],
                       ensure_ascii=False), encoding="utf-8")
    if args.json:
        slim = {k: v for k, v in report.items() if k not in ("results", "fails")}
        print(json.dumps(slim, ensure_ascii=False))
    else:
        _print_summary(report)
    if args.baseline:
        _print_vs_baseline(report, args.baseline, report["results"])
    if args.time:
        print(f"workers={workers} elapsed={elapsed:.2f}s "
              f"avg={report['avg_ms_per_row']:.2f}ms/row (n={n_rows})", file=sys.stderr)


if __name__ == "__main__":
    main()
