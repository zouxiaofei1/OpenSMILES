/* Layer Benchmark page: per-layer LLM sample + score + single-molecule analyze. */
import { $, state, api, API, escapeHtml } from "./core.js";

function lb$(id) {
  return document.getElementById("lb-" + id);
}

export function stopLbPolling() {
  if (state.lbPollTimer) {
    clearInterval(state.lbPollTimer);
    state.lbPollTimer = null;
  }
  state.lbGenerating = false;
}

function showLbProgress(total) {
  var prog = lb$("progress");
  if (prog) prog.classList.remove("hidden");
  var bar = lb$("progress-bar");
  if (bar) bar.style.width = "0%";
  var text = lb$("progress-text");
  if (text && total) text.textContent = "0 / " + total + " (0%)";
}

function hideLbProgress() {
  var prog = lb$("progress");
  if (prog) prog.classList.add("hidden");
}

function setLbStatus(msg, color) {
  var statusEl = lb$("status");
  if (!statusEl) return;
  statusEl.textContent = msg;
  if (color) statusEl.style.color = color;
  else statusEl.style.color = "";
}

var LB_LAYER_META = {
  0: { title: "Layer 0 — Molecule Preprocessor", desc: "SMILES 解析 / 无机-有机分类 / 盐解离" },
  1: { title: "Layer 1 — Structural Analyzer", desc: "官能团 / 环系统 / 不饱和度分析" },
};

function setLbLayer(layer) {
  state.lbLayer = layer;
  document.querySelectorAll(".lb-item").forEach(function (b) {
    var l = parseInt(b.getAttribute("data-lb-layer"), 10);
    if (l === layer) b.classList.add("active");
    else b.classList.remove("active");
  });
  renderLbPanel(layer);
}

function lbPlaceholderHtml(layer) {
  return (
    '<div class="card"><div class="card-head row-between"><h3 class="card-title">Layer ' + layer +
    "</h3><span class=\"chip\">未实装</span></div>" +
    '<p class="muted-text">功能尚未实装。布局与 Layer 0/1 相同，后续将接入该层输出与 LLM 数据的对比。</p></div>'
  );
}

function lbPanelHtml(layer) {
  var m = LB_LAYER_META[layer];
  return (
    '<div class="card" style="margin-bottom: var(--space-2xl)">' +
    '<div class="card-head row-between"><h3 class="card-title">单物质 Layer ' + layer +
    ' 分析</h3><span class="chip">LLM</span></div>' +
    '<div class="lb-one-controls">' +
    '<input id="lb-one-smiles" class="input mono" type="text" placeholder="输入 SMILES，如 CCO" autocomplete="off" spellcheck="false" />' +
    '<button id="lb-one-run" class="btn btn-primary" type="button">分析</button></div>' +
    '<div id="lb-one-result" class="lb-one-result hidden" aria-live="polite"></div></div>' +
    '<div class="card" style="margin-bottom: var(--space-2xl)">' +
    '<div class="card-head row-between"><h3 class="card-title">' + m.title +
    '</h3><span id="lb-status" class="chip">就绪</span></div>' +
    '<p class="muted-text small">' + m.desc + "</p>" +
    '<div class="lb-controls">' +
    '<div class="lb-field">' +
    '<label class="field-label" for="lb-ratio">抽样比例</label>' +
    '<input id="lb-ratio" class="input" type="range" min="1" max="100" value="10" step="1" />' +
    '<span id="lb-ratio-val" class="muted-text small mono">10%</span>' +
    '<span class="muted-text small">按 tier 分层随机抽取，覆盖难易不均</span></div>' +
    '<div class="lb-actions">' +
    '<button id="lb-refresh" class="btn btn-secondary" type="button">刷新 benchmark 数据</button>' +
    '<button id="lb-score" class="btn btn-primary" type="button">跑分</button></div></div>' +
    '<div id="lb-progress" class="hidden" style="margin-top: var(--space-xl)">' +
    '<div class="br-progress-bar-wrap"><div id="lb-progress-bar" class="br-progress-bar" style="width:0%"></div></div>' +
    '<div class="lb-progress-stats"><span id="lb-progress-text" class="muted-text small">0 / 0</span>' +
    '<span class="muted-text small">正在调用 LLM 生成 Layer ' + layer +
    ' 基准数据…</span></div></div></div>' +
    '<div id="lb-result" class="hidden">' +
    '<div class="card" style="margin-bottom: var(--space-2xl)">' +
    '<div class="card-head row-between"><h3 class="card-title">跑分结果：真实 Layer ' + layer +
    ' vs LLM 基准</h3><span id="lb-score-info" class="muted-text small mono"></span></div>' +
    '<div id="lb-score-table-wrap"></div></div>' +
    '<div class="card"><div class="card-head row-between"><h3 class="card-title">不一致示例</h3>' +
    '<span id="lb-diff-count" class="muted-text small"></span></div>' +
    '<div class="benchmark-table-wrap" id="lb-diff-wrap"></div></div></div>'
  );
}

export function renderLbPanel(layer) {
  var content = document.getElementById("lb-content");
  if (!content) return;
  content.innerHTML = (layer === 0 || layer === 1) ? lbPanelHtml(layer) : lbPlaceholderHtml(layer);
  hideLbProgress();
  setLbStatus("就绪");
  if (layer === 0 || layer === 1) loadLayerData(layer);
}

async function loadLayerData(layer) {
  try {
    var data = await api(API.layerData + "?layer=" + layer);
    if (!data || !data.ok) return;
    var ratioEl = lb$("ratio");
    var ratioVal = lb$("ratio-val");
    if (data.ratio > 0 && ratioEl) {
      ratioEl.value = String(Math.round(data.ratio));
      if (ratioVal) ratioVal.textContent = Math.round(data.ratio) + "%";
    }
    var truth = data.truth || {};
    var n = data.sample_n || 0;
    if (truth.running || (n > 0 && truth.done < truth.total)) {
      state.lbGenerating = true;
      showLbProgress(truth.total);
      startLbPolling();
      setLbStatus("生成基准中… " + (truth.done || 0) + "/" + (truth.total || n), "#f59e0b");
    } else if (n > 0 && truth.ready) {
      setLbStatus("基准就绪 · " + n + " 条");
    } else if (n > 0) {
      setLbStatus(n + " 条样本");
    }
    if (data.score && data.score.ok) {
      state.lbScore = data.score;
      renderLbScore(data.score);
    }
  } catch (_) {
    /* ignore */
  }
}

async function startLayerSample() {
  if (state.lbGenerating) return;
  var btn = lb$("refresh");
  if (btn) btn.disabled = true;
  setLbStatus("抽样中…", "#f59e0b");
  var ratioEl = lb$("ratio");
  var ratio = ratioEl ? parseFloat(ratioEl.value) : 10;
  if (!ratio || isNaN(ratio)) ratio = 10;
  try {
    var data = await api(API.layerSample, {
      method: "POST",
      body: JSON.stringify({ ratio: ratio, layer: state.lbLayer }),
    });
    if (!data || !data.ok) {
      setLbStatus("错误: " + ((data && data.error) || "未知"), "#ef4444");
      return;
    }
    setLbStatus("生成基准中… " + data.n + " 条", "#f59e0b");
    showLbProgress(data.n);
    state.lbGenerating = true;
    startLbPolling();
  } catch (err) {
    setLbStatus("请求失败: " + (err.message || String(err)), "#ef4444");
  } finally {
    if (btn) btn.disabled = false;
  }
}

async function pollLayerStatus() {
  var status = null;
  try {
    status = await api(API.layerSampleStatus + "?layer=" + state.lbLayer);
  } catch (_) {
    return;
  }
  if (!status || !status.ok) return;
  var done = status.done || 0;
  var total = status.total || 0;
  var pct = total > 0 ? Math.round(100 * done / total) : 0;
  var bar = lb$("progress-bar");
  if (bar) bar.style.width = pct + "%";
  var text = lb$("progress-text");
  if (text) text.textContent = done + " / " + total + " (" + pct + "%)";
  if (status.ready) {
    stopLbPolling();
    hideLbProgress();
    setLbStatus("基准就绪 · " + done + " 条");
    return;
  }
  if (status.error) {
    stopLbPolling();
    hideLbProgress();
    setLbStatus("生成失败: " + status.error, "#ef4444");
    return;
  }
  // 仍在运行或等待首条进度落盘：继续轮询；成败由后端 ready/error 兜底
}

function startLbPolling() {
  if (state.lbPollTimer) return;
  state.lbPollTimer = setInterval(pollLayerStatus, 3000);
}

async function runLayerScore() {
  var btn = lb$("score");
  if (btn) btn.disabled = true;
  setLbStatus("跑分中…", "#f59e0b");
  try {
    var data = await api(API.layerScore, {
      method: "POST",
      body: JSON.stringify({ layer: state.lbLayer }),
    });
    if (!data || !data.ok) {
      setLbStatus("错误: " + ((data && data.error) || "未知"), "#ef4444");
      return;
    }
    setLbStatus("完成", "#86efac");
    state.lbScore = data;
    renderLbScore(data);
  } catch (err) {
    setLbStatus("跑分失败: " + (err.message || String(err)), "#ef4444");
  } finally {
    if (btn) btn.disabled = false;
  }
}

function lbScoreFields(layer) {
  if (layer === 1) {
    return [
      ["parseable", "可解析 (parseable)"],
      ["fg", "官能团集合一致"],
      ["n_rings", "环数"],
      ["n_ring_systems", "环系"],
      ["has_double_bond", "C=C 双键"],
      ["has_triple_bond", "C≡C 三键"],
      ["complete", "全部字段一致"],
    ];
  }
  return [
    ["parseable", "可解析 (parseable)"],
    ["classification", "有机 / 无机分类"],
    ["salt_present", "是否为盐"],
    ["salt_type", "盐类型（均为盐时）"],
    ["n_fragments", "片段数"],
    ["complete", "全部字段一致"],
  ];
}

function lbScoreRows(score) {
  var fields = lbScoreFields(state.lbLayer);
  return fields
    .map(function (pair) {
      var key = pair[0], label = pair[1];
      var m = score[key] || {};
      var pct = m.pct != null ? m.pct.toFixed(1) + "%" : "—";
      var det = m.n ? "(" + m.ok + "/" + m.n + ")" : "";
      var bw = m.n ? Math.round(100 * m.ok / m.n) : 0;
      return (
        '<div class="lb-score-row">' +
        '<span class="lb-score-label">' + label + "</span>" +
        '<span class="lb-score-bar-wrap"><span class="lb-score-bar" style="width:' + bw + '%"></span></span>' +
        '<span class="lb-score-value mono">' + pct + " " + det + "</span>" +
        "</div>"
      );
    })
    .join("");
}

function lbSideHtml(s, layer) {
  if (!s) return '<span class="muted-text">—</span>';
  if (layer === 1) {
    var parts = [];
    parts.push(s.parseable ? "解析✓" : "解析✗");
    parts.push("FG:[" + ((s.functional_groups || []).join(",")) + "]");
    parts.push("环:" + (s.n_rings != null ? s.n_rings : "?"));
    parts.push("环系:" + (s.n_ring_systems != null ? s.n_ring_systems : "?"));
    parts.push((s.has_double_bond ? "双键" : "") + (s.has_triple_bond ? "三键" : ""));
    return escapeHtml(parts.join(" · "));
  }
  var parts = [];
  parts.push(s.parseable ? "解析✓" : "解析✗");
  parts.push(s.classification || "?");
  parts.push(s.salt ? "盐:" + s.salt : "无盐");
  parts.push("frag=" + (s.n_fragments != null ? s.n_fragments : "?"));
  return escapeHtml(parts.join(" · "));
}

function renderLbScore(score) {
  var resultEl = lb$("result");
  if (!resultEl) return;
  resultEl.classList.remove("hidden");
  var info = lb$("score-info");
  if (info) info.textContent = "n=" + score.n + " · ratio=" + (score.ratio != null ? score.ratio + "%" : "?");
  var wrap = lb$("score-table-wrap");
  if (wrap) wrap.innerHTML = lbScoreRows(score);
  var diffWrap = lb$("diff-wrap");
  var diffCount = lb$("diff-count");
  var diffs = score.diffs || [];
  if (diffCount) diffCount.textContent = diffs.length + " 条";
  if (!diffWrap) return;
  if (!diffs.length) {
    diffWrap.innerHTML = '<p class="muted-text small" style="padding:1rem">无不一致示例</p>';
    return;
  }
  var rows = diffs
    .map(function (d) {
      var layer = state.lbLayer;
      return (
        "<tr>" +
        '<td class="mono small">' + escapeHtml(d.smiles || "") + "</td>" +
        '<td class="small">' + lbSideHtml(d.real, layer) + "</td>" +
        '<td class="small">' + lbSideHtml(d.llm, layer) + "</td>" +
        "</tr>"
      );
    })
    .join("");
  diffWrap.innerHTML =
    '<table class="benchmark-table lb-diff-table"><thead><tr>' +
    "<th>SMILES</th><th>真实 L0</th><th>LLM 基准</th></tr></thead><tbody>" +
    rows +
    "</tbody></table>";
}

function lbOneHtml(r) {
  var saltHtml = r.salt
    ? Object.keys(r.salt)
        .map(function (k) {
          return "<dt>" + escapeHtml(k) + "</dt><dd class=\"mono\">" + escapeHtml(String(r.salt[k])) + "</dd>";
        })
        .join("")
    : '<dt>salt</dt><dd class="mono muted-text">null</dd>';
  return (
    '<dl class="kv">' +
    "<div><dt>parseable</dt><dd class=\"mono\">" + (r.parseable ? "true" : "false") + "</dd></div>" +
    "<div><dt>classification</dt><dd class=\"mono\">" + escapeHtml(r.classification || "?") + "</dd></div>" +
    "<div><dt>n_fragments</dt><dd class=\"mono\">" + escapeHtml(String(r.n_fragments != null ? r.n_fragments : "?")) + "</dd></div>" +
    "<div><dt>organic_smiles</dt><dd class=\"mono\">" + escapeHtml(r.organic_smiles || "null") + "</dd></div>" +
    saltHtml +
    "<div><dt>notes</dt><dd>" + escapeHtml(r.notes || "") + "</dd></div>" +
    "</dl>"
  );
}

async function runLayerOne() {
  var input = document.getElementById("lb-one-smiles");
  var smiles = (input && input.value || "").trim();
  if (!smiles) return;
  var btn = document.getElementById("lb-one-run");
  var resultEl = document.getElementById("lb-one-result");
  if (btn) btn.disabled = true;
  if (resultEl) {
    resultEl.classList.remove("hidden");
    resultEl.innerHTML = '<p class="muted-text small">分析中…</p>';
  }
  try {
    var data = await api(API.layerAnalyzeOne, {
      method: "POST",
      body: JSON.stringify({ smiles: smiles, layer: state.lbLayer }),
    });
    if (!data || !data.ok) {
      if (resultEl) resultEl.textContent = "错误: " + ((data && data.error) || "未知");
      return;
    }
    if (resultEl) resultEl.innerHTML = lbOneHtml(data.result);
  } catch (err) {
    if (resultEl) resultEl.textContent = "请求失败: " + (err.message || String(err));
  } finally {
    if (btn) btn.disabled = false;
  }
}

export function bindLayerBenchmark() {
  // Layer Benchmark controls（面板为动态渲染，用事件委托）
  document.querySelectorAll(".lb-item").forEach(function (b) {
    b.addEventListener("click", function () {
      setLbLayer(parseInt(b.getAttribute("data-lb-layer"), 10));
    });
  });
  var lbContent = document.getElementById("lb-content");
  if (lbContent) {
    lbContent.addEventListener("click", function (ev) {
      var t = ev.target;
      if (t.closest && t.closest("#lb-refresh")) { startLayerSample(); return; }
      if (t.closest && t.closest("#lb-score")) { runLayerScore(); return; }
      if (t.closest && t.closest("#lb-one-run")) { runLayerOne(); }
    });
    lbContent.addEventListener("input", function (ev) {
      if (ev.target && ev.target.id === "lb-ratio") {
        var v = parseInt(ev.target.value, 10) || 0;
        var valEl = document.getElementById("lb-ratio-val");
        if (valEl) valEl.textContent = v + "%";
      }
    });
    lbContent.addEventListener("keydown", function (ev) {
      if (ev.key === "Enter" && ev.target && ev.target.id === "lb-one-smiles") runLayerOne();
    });
  }
}
