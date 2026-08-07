/* Benchmark Run page: launch a full accuracy run, poll progress, show result. */
import { state, api, API } from "./core.js";

function br$(id) {
  return document.getElementById("br-" + id);
}

export function stopBrPolling() {
  if (state.brPollTimer) {
    clearInterval(state.brPollTimer);
    state.brPollTimer = null;
  }
  state.brRunning = false;
}

function fmtPct(ok, n) {
  if (!n) return "—";
  return (100 * ok / n).toFixed(1) + "% (" + ok + "/" + n + ")";
}

function fmtSec(s) {
  if (s == null || isNaN(s)) return "";
  if (s < 60) return s.toFixed(0) + "s";
  var m = Math.floor(s / 60);
  var sec = Math.round(s % 60);
  return m + "m " + sec + "s";
}

async function startBenchmarkRun() {
  if (state.brRunning) return;

  var workersEl = document.getElementById("br-workers");
  var timeoutEl = document.getElementById("br-timeout");
  var limitEl = document.getElementById("br-limit");

  var workers = parseInt(workersEl && workersEl.value || "0", 10) || 0;
  var timeout = parseFloat(timeoutEl && timeoutEl.value || "1.0") || 1.0;
  var limitVal = (limitEl && limitEl.value || "").trim();
  var limit = limitVal ? parseInt(limitVal, 10) || null : null;

  var btn = br$("run");
  if (btn) btn.disabled = true;

  var statusEl = br$("status");
  if (statusEl) { statusEl.textContent = "启动中…"; statusEl.style.color = "#f59e0b"; }

  // Hide previous result
  var resultEl = br$("result");
  if (resultEl) resultEl.classList.add("hidden");

  try {
    var body = { workers: workers, timeout: timeout };
    if (limit != null) body.limit = limit;
    var data = await api(API.benchmarkRun, { method: "POST", body: JSON.stringify(body) });
    if (!data || !data.ok) {
      var msg = (data && data.error) || "启动失败";
      if (statusEl) { statusEl.textContent = "错误: " + msg; statusEl.style.color = "#ef4444"; }
      if (btn) btn.disabled = false;
      return;
    }
  } catch (err) {
    if (statusEl) { statusEl.textContent = "请求失败: " + (err.message || String(err)); statusEl.style.color = "#ef4444"; }
    if (btn) btn.disabled = false;
    return;
  }

  state.brRunning = true;
  if (statusEl) { statusEl.textContent = "运行中…"; statusEl.style.color = "#f59e0b"; }

  // Show progress bar
  var progEl = br$("progress");
  if (progEl) progEl.classList.remove("hidden");

  // Start polling
  state.brPollTimer = setInterval(pollBrStatus, 1500);
}

async function pollBrStatus() {
  try {
    var status = await api(API.benchmarkRunStatus);
  } catch (_) {
    return;
  }

  if (!status) return;

  // Update progress bar
  var done = status.done || 0;
  var total = status.total || 0;
  var pct = total > 0 ? Math.round(100 * done / total) : 0;

  var bar = br$("progress-bar");
  if (bar) bar.style.width = pct + "%";

  var text = br$("prog-text");
  if (text) text.textContent = done + " / " + total + " (" + pct + "%)";

  var rate = br$("prog-rate");
  if (rate && status.rate) rate.textContent = status.rate.toFixed(1) + " row/s";

  var eta = br$("prog-eta");
  if (eta && status.eta) eta.textContent = "eta " + fmtSec(status.eta);

  // Running accuracy
  var sofarEl = br$("prog-sofar");
  if (sofarEl && status.sofar) {
    var sf = status.sofar;
    sofarEl.textContent =
      "实时 EN=" + fmtPct(sf.en_ok, sf.en_n) +
      "  ZH=" + fmtPct(sf.zh_ok, sf.zh_n) +
      "  dual=" + fmtPct(sf.dual_ok, sf.dual_n);
  }

  // Check if done
  if (!status.running) {
    stopBrPolling();

    var statusEl = br$("status");
    var btn = br$("run");
    if (btn) btn.disabled = false;

    if (status.error) {
      if (statusEl) { statusEl.textContent = "错误: " + status.error; statusEl.style.color = "#ef4444"; }
      return;
    }

    // Fetch final result
    try {
      var res = await api(API.benchmarkRunResult);
      if (res && res.ok && res.result) {
        showBrResult(res.result);
        if (statusEl) { statusEl.textContent = "完成"; statusEl.style.color = "#86efac"; }
      } else {
        if (statusEl) { statusEl.textContent = "完成 (无结果)"; statusEl.style.color = "#f59e0b"; }
      }
    } catch (err) {
      if (statusEl) { statusEl.textContent = "获取结果失败: " + (err.message || String(err)); statusEl.style.color = "#ef4444"; }
    }
  }
}

function showBrResult(r) {
  var resultEl = br$("result");
  if (resultEl) resultEl.classList.remove("hidden");

  // Accuracy cards
  var accEn = br$("acc-en");
  if (accEn) accEn.textContent = (r.acc_en != null ? r.acc_en.toFixed(1) + "%" : "—");
  var accEnDet = br$("acc-en-detail");
  if (accEnDet) accEnDet.textContent = fmtPct(r.ok_en, r.n_en);

  var accZh = br$("acc-zh");
  if (accZh) accZh.textContent = (r.acc_zh != null ? r.acc_zh.toFixed(1) + "%" : "—");
  var accZhDet = br$("acc-zh-detail");
  if (accZhDet) accZhDet.textContent = fmtPct(r.ok_zh, r.n_zh);

  var accDual = br$("acc-dual");
  if (accDual) accDual.textContent = (r.acc_dual != null ? r.acc_dual.toFixed(1) + "%" : "—");
  var accDualDet = br$("acc-dual-detail");
  if (accDualDet) accDualDet.textContent = fmtPct(r.ok_dual, r.n_dual);

  var fails = br$("fails");
  if (fails) fails.textContent = String(r.fails != null ? r.fails : "—");

  // Timing
  var timing = br$("timing");
  if (timing) {
    var parts = [];
    if (r.workers) parts.push("workers=" + r.workers);
    if (r.elapsed_sec) parts.push("elapsed " + r.elapsed_sec.toFixed(1) + "s");
    if (r.avg_ms_per_row) parts.push(r.avg_ms_per_row.toFixed(1) + "ms/row");
    if (r.n) parts.push("n=" + r.n);
    timing.textContent = parts.join("  ·  ");
  }
}

export function bindBenchmarkRun() {
  // Benchmark Run controls
  var brRun = document.getElementById("br-run");
  if (brRun) {
    brRun.addEventListener("click", startBenchmarkRun);
  }
}
