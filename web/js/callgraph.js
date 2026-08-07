/* Call Graph page: sample call-graph data, ranking, callee chain, layered SVG. */
import { $, state, escapeHtml, api, API } from "./core.js";

function renderCgMeta(meta) {
  var el = $("cg-meta");
  if (!el || !meta) return;
  var cache = meta.cached ? "缓存命中" : "本次采样";
  el.textContent =
    "分子 " + meta.n_actual + "/" + meta.n +
    " · 总耗时 " + (meta.total_s * 1000).toFixed(1) + "ms" +
    " · 函数 " + meta.n_nodes + "/" + meta.n_nodes_total +
    " · 调用边 " + meta.n_edges +
    " · " + cache +
    " · 接口耗时 " + meta.elapsed_ms + "ms";
}

function cgNFromSlider() {
  var pct = parseFloat(($("cg-n") && $("cg-n").value) || 5);
  var total = state.cgSourceTotal || 4070;
  var n = Math.max(1, Math.round(total * pct / 100));
  var v = $("cg-n-val");
  if (v) v.textContent = pct + "% ≈ " + n + " 分子";
  return n;
}

function cgRankRow(n, idx, val) {
  var mod = (n.module || "").split("/").pop().replace(".py", "");
  var name = n.label + (n.line != null ? " @ " + mod + ":" + n.line : " (" + mod + ")");
  return (
    '<div class="cg-rank-row cg-rank-clickable" data-id="' + n.id + '" title="点击查看下层调用链">' +
    '<span class="cg-rank-idx mono">' + idx + "</span>" +
    '<span class="cg-rank-name mono" title="' + escapeHtml(n.module || "") + '">' + escapeHtml(name) + "</span>" +
    '<span class="cg-rank-val mono">' + val + "</span>" +
    "</div>"
  );
}

/* ── 下层调用链：展开/折叠树 + 右键进入 + SVG 高亮 ─────────────────── */
var cgRootId = null, cgSelId = null, cgOpen = {};

function cgChainById() {
  var idx = {};
  ((state.cgData && state.cgData.nodes) || []).forEach(function (n) { idx[n.id] = n; });
  return idx;
}

function cgChainShort(nd) {
  var mod = (nd.module || "").split("/").pop().replace(".py", "");
  return nd.label + " @" + mod + ":" + nd.line;
}

// 每个节点的直接 callee 边，按累计耗时降序
function cgOutEdges(byId, data) {
  var out = {};
  (data && data.edges || []).forEach(function (e) {
    (out[e.from] = out[e.from] || []).push(e);
  });
  Object.keys(out).forEach(function (f) {
    out[f].sort(function (a, b) {
      var sa = byId[a.to] ? byId[a.to].cum_s : 0;
      var sb = byId[b.to] ? byId[b.to].cum_s : 0;
      return sb - sa;
    });
  });
  return out;
}

function cgChev(id, hasKids) {
  if (!hasKids) return '<span class="cg-chev cg-chev-none" aria-hidden="true"></span>';
  return '<span class="cg-chev" aria-hidden="true">' + (cgOpen[id] ? "▾" : "▸") + "</span>";
}

function cgTreeHtml(byId, out, root) {
  function fmt(nd) {
    return cgChainShort(nd) + "  " + nd.cum_pct + "% (" + nd.ncalls.toLocaleString() + "×)";
  }
  function row(nd, depth) {
    var pad = depth * 16 + 8;
    var kids = out[nd.id] || [];
    var sel = nd.id === cgSelId ? " cg-sel-row" : "";
    return (
      '<div class="cg-chain-row cg-chain-name' + sel + '" data-id="' + nd.id + '" style="padding-left:' + pad + 'px">' +
      cgChev(nd.id, kids.length > 0) +
      '<span class="mono">' + escapeHtml(fmt(nd)) + "</span>" +
      "</div>"
    );
  }
  function walk(nd, depth, path) {
    var kids = out[nd.id] || [];
    if (!kids.length) return "";
    var html = "";
    kids.forEach(function (e) {
      var nd2 = byId[e.to];
      if (!nd2) return;
      if (path[nd2.id]) {
        html += '<div class="cg-chain-row cg-chain-loop" style="padding-left:' + (depth * 16 + 24) + 'px">↺ 循环</div>';
        return;
      }
      html += row(nd2, depth);
      if (cgOpen[nd2.id]) {
        var p2 = Object.assign({}, path);
        p2[nd2.id] = true;
        html += walk(nd2, depth + 1, p2);
      }
    });
    return html;
  }
  var rootNd = byId[root];
  if (!rootNd) return '<p class="muted-text small">无数据</p>';
  var rootSel = rootNd.id === cgSelId ? " cg-sel-row" : "";
  var head =
    '<div class="cg-chain-row cg-chain-root cg-chain-name' + rootSel + '" data-id="' + rootNd.id + '">' +
    cgChev(rootNd.id, (out[rootNd.id] || []).length > 0) +
    "<b>" + escapeHtml(fmt(rootNd)) + "</b></div>";
  var body = cgOpen[rootNd.id] ? walk(rootNd, 1, { [rootNd.id]: true }) : "";
  return head + body;
}

// 进入某个函数：切换树根、默认展开第一层、选中并高亮
function renderCalleeChain(rootId) {
  var el = $("cg-chain");
  if (!el) return;
  cgRootId = rootId;
  cgOpen[rootId] = true; // 进入后默认展开第一层
  cgSelId = rootId;
  var byId = cgChainById();
  var out = cgOutEdges(byId, state.cgData);
  el.innerHTML = cgTreeHtml(byId, out, rootId);
  applyCgHighlight();
}

function renderCalleeChainByName(query) {
  var data = state.cgData;
  if (!data || !data.nodes) return;
  var q = String(query || "").trim().toLowerCase();
  if (!q) return;
  var hits = data.nodes.filter(function (n) {
    return n.label.toLowerCase().indexOf(q) !== -1 ||
      (n.module + ":" + n.line).toLowerCase().indexOf(q) !== -1;
  });
  if (!hits.length) {
    var el = $("cg-chain");
    if (el) el.innerHTML = '<p class="muted-text small">未找到匹配函数："' + escapeHtml(query) + '"</p>';
    return;
  }
  renderCalleeChain(hits[0].id);
  if (hits.length > 1) {
    var el = $("cg-chain");
    if (el) el.innerHTML += '<p class="muted-text small">匹配 ' + hits.length + " 个函数，已展开第一个；可用 文件:行号 精确定位</p>";
  }
}

/* ── SVG 图上高亮选中节点 ─────────────────────────────────────────── */
function applyCgHighlight() {
  var view = $("cg-svg-view");
  if (!view) return;
  var prev = view.querySelectorAll(".cg-sel");
  for (var i = 0; i < prev.length; i++) prev[i].classList.remove("cg-sel");
  if (cgSelId == null) return;
  var g = view.querySelector('[id="cg' + cgSelId + '"]');
  if (g) g.classList.add("cg-sel");
}

/* ── 右键菜单：进入该函数 ─────────────────────────────────────────── */
var _cgCtx = null;
function showCgCtxMenu(x, y, id) {
  if (!_cgCtx) {
    _cgCtx = document.createElement("div");
    _cgCtx.className = "cg-ctx-menu";
    _cgCtx.innerHTML =
      '<div class="cg-ctx-item" data-act="enter">进入此函数查看下层调用链</div>' +
      '<div class="cg-ctx-item" data-act="copy">复制函数名</div>';
    document.body.appendChild(_cgCtx);
    _cgCtx.addEventListener("click", function (ev) {
      var item = ev.target.closest(".cg-ctx-item");
      if (!item || _cgCtx._id == null) return;
      var act = item.getAttribute("data-act");
      var id = _cgCtx._id;
      hideCgCtxMenu();
      if (act === "enter") {
        renderCalleeChain(id);
      } else if (act === "copy") {
        var nd = cgChainById()[id];
        var name = nd ? (nd.label + " @" + nd.line) : String(id);
        try { navigator.clipboard.writeText(name); } catch (e) {}
      }
    });
    window.addEventListener("click", hideCgCtxMenu);
    window.addEventListener("blur", hideCgCtxMenu);
  }
  _cgCtx._id = id;
  _cgCtx.hidden = false;
  var r = _cgCtx.getBoundingClientRect();
  var nx = Math.min(x, window.innerWidth - r.width - 8);
  var ny = Math.min(y, window.innerHeight - r.height - 8);
  _cgCtx.style.left = Math.max(0, nx) + "px";
  _cgCtx.style.top = Math.max(0, ny) + "px";
}
function hideCgCtxMenu() {
  if (_cgCtx) { _cgCtx.hidden = true; _cgCtx._id = null; }
}

function bindCgChainControls() {
  var go = $("cg-chain-go");
  var input = $("cg-chain-input");
  function run() {
    if (input) renderCalleeChainByName(input.value);
  }
  if (go) go.addEventListener("click", run);
  if (input) {
    input.addEventListener("keydown", function (ev) {
      if (ev.key === "Enter") run();
    });
  }
  // rank 排行行：点击进入该函数的调用链
  ["cg-rank-calls", "cg-rank-time"].forEach(function (id) {
    var list = $(id);
    if (list) {
      list.addEventListener("click", function (ev) {
        var row = ev.target.closest(".cg-rank-clickable");
        if (row) renderCalleeChain(parseInt(row.getAttribute("data-id"), 10));
      });
    }
  });
  // 下层调用链树：容器级事件委托（点击展开/折叠 + 选中，右键进入）
  var chain = $("cg-chain");
  if (!chain) return;
  chain.addEventListener("click", function (ev) {
    var row = ev.target.closest(".cg-chain-name");
    if (!row) return;
    var id = parseInt(row.getAttribute("data-id"), 10);
    cgOpen[id] = !cgOpen[id];
    cgSelId = id;
    var byId = cgChainById();
    var out = cgOutEdges(byId, state.cgData);
    chain.innerHTML = cgTreeHtml(byId, out, cgRootId);
    applyCgHighlight();
  });
  chain.addEventListener("contextmenu", function (ev) {
    var row = ev.target.closest(".cg-chain-name");
    if (!row) return;
    ev.preventDefault();
    showCgCtxMenu(ev.clientX, ev.clientY, parseInt(row.getAttribute("data-id"), 10));
  });
}

function renderCgRanks() {
  var data = state.cgData;
  if (!data || !data.nodes) return;
  var t = state.cgThreshold;
  var nodes = data.nodes.filter(function (n) {
    return n.cum_pct >= t;
  });
  function build(key, field, fmt) {
    var st = state.cgRank[key] || { asc: false, all: false };
    var arr = nodes.slice().sort(function (a, b) {
      return st.asc ? a[field] - b[field] : b[field] - a[field];
    });
    var shown = st.all ? arr : arr.slice(0, 20);
    return shown.map(function (n, i) { return cgRankRow(n, i + 1, fmt(n)); }).join("");
  }
  var cEl = $("cg-rank-calls");
  var tEl = $("cg-rank-time");
  if (cEl) {
    var rc = build("calls", "ncalls", function (n) { return n.ncalls.toLocaleString() + "×"; });
    cEl.innerHTML = rc || '<p class="muted-text small">无节点</p>';
  }
  if (tEl) {
    var rt = build("time", "cum_s", function (n) { return n.cum_pct + "%"; });
    tEl.innerHTML = rt || '<p class="muted-text small">无节点</p>';
  }
  updateCgRankBtns();
}

function updateCgRankBtns() {
  ["calls", "time"].forEach(function (key) {
    var st = state.cgRank[key] || { asc: false, all: false };
    var dir = $("cg-rank-" + key + "-dir");
    if (dir) dir.textContent = st.asc ? "↑ 升序" : "↓ 降序";
    var more = $("cg-rank-" + key + "-more");
    if (more) more.textContent = st.all ? "全部" : "Top 20";
  });
}

function bindCgRankButtons() {
  ["calls", "time"].forEach(function (key) {
    var dir = $("cg-rank-" + key + "-dir");
    if (dir) {
      dir.addEventListener("click", function () {
        state.cgRank[key] = state.cgRank[key] || { asc: false, all: false };
        state.cgRank[key].asc = !state.cgRank[key].asc;
        renderCgRanks();
      });
    }
    var more = $("cg-rank-" + key + "-more");
    if (more) {
      more.addEventListener("click", function () {
        state.cgRank[key] = state.cgRank[key] || { asc: false, all: false };
        state.cgRank[key].all = !state.cgRank[key].all;
        renderCgRanks();
      });
    }
  });
}

function startCgProgress() {
  var wrap = $("cg-progress");
  var bar = $("cg-progress-bar");
  var text = $("cg-progress-text");
  if (wrap) wrap.classList.remove("hidden");
  if (bar) bar.style.width = "0%";
  if (text) text.textContent = "0 / 0 (0%)";
  stopCgProgress();
  state.cgPollTimer = setInterval(pollCgProgress, 400);
}

function stopCgProgress() {
  if (state.cgPollTimer) {
    clearInterval(state.cgPollTimer);
    state.cgPollTimer = null;
  }
}

async function pollCgProgress() {
  var p;
  try {
    p = await api(API.callGraphProgress);
  } catch (_) {
    return;
  }
  if (!p || !p.ok) return;
  var done = p.done || 0;
  var total = p.total || 0;
  var pct = total > 0 ? Math.round(100 * done / total) : 0;
  var bar = $("cg-progress-bar");
  if (bar) bar.style.width = pct + "%";
  var text = $("cg-progress-text");
  if (text) text.textContent = done + " / " + total + " (" + pct + "%)";
  if (!p.active) stopCgProgress();
}

export async function loadCallGraph(force) {
  if (state.cgLoading) return;
  var n = cgNFromSlider();
  var st = $("cg-status");
  var btn = $("cg-regenerate");
  var empty = $("cg-empty");
  state.cgLoading = true;
  if (btn) btn.disabled = true;
  if (st) { st.textContent = "采样中…（约 1–3s）"; st.style.color = "#f59e0b"; }
  if (empty) empty.hidden = true;
  startCgProgress();
  try {
    var q = "?n=" + n + "&floor_pct=" + state.cgThreshold + (force ? "&refresh=1" : "");
    var data = await api(API.callGraph + q);
    if (!data || !data.ok) throw new Error((data && data.error) || "无数据");
    state.cgData = data;
    if (data.meta && data.meta.source_total) {
      state.cgSourceTotal = data.meta.source_total;
    }
    cgNFromSlider();
    // 重采样后旧分层 SVG 失效
    state.cgSvg = null;
    var s = $("cg-svg-status");
    if (s) { s.textContent = "数据已更新，请重新生成"; s.style.color = "#f59e0b"; }
    renderCgMeta(data.meta);
    renderCgRanks();
    // 采样完成：进度条收满
    var bar = $("cg-progress-bar");
    if (bar) bar.style.width = "100%";
    var pt = $("cg-progress-text");
    var n_actual = (data.meta && data.meta.n_actual) || 0;
    if (pt) pt.textContent = n_actual ? n_actual + " / " + n_actual + " (100%)" : "完成";
    if (st) { st.textContent = "就绪"; st.style.color = "#86efac"; }
    return true;
  } catch (err) {
    var cgErrMsg = err.message || String(err);
    if (/busy|in progress/i.test(cgErrMsg)) {
      // A historical sampling for another commit is running; retry shortly.
      if (empty) { empty.hidden = false; empty.textContent = "上一历史采样进行中，2 秒后自动重试…"; }
      if (st) { st.textContent = "等待中…"; st.style.color = "#f59e0b"; }
      setTimeout(function () { loadCallGraph(); }, 2000);
      return false;
    }
    if (empty) { empty.hidden = false; empty.textContent = "加载失败: " + cgErrMsg; }
    if (st) { st.textContent = "错误"; st.style.color = "#ef4444"; }
    return false;
  } finally {
    state.cgLoading = false;
    if (btn) btn.disabled = false;
    stopCgProgress();
    var wrap = $("cg-progress");
    if (wrap) {
      setTimeout(function () {
        wrap.classList.add("hidden");
      }, 600);
    }
  }
}

function bindCgControls() {
  $("cg-regenerate") &&
    $("cg-regenerate").addEventListener("click", function () {
      regenerate();
    });
  $("cg-threshold") &&
    $("cg-threshold").addEventListener("input", function (ev) {
      var tv = parseFloat(ev.target.value);
      state.cgThreshold = isNaN(tv) ? 0 : tv;
      var v = $("cg-threshold-val");
      if (v) v.textContent = state.cgThreshold.toFixed(1) + "%";
      renderCgRanks();
    });
  $("cg-n") &&
    $("cg-n").addEventListener("input", function () {
      cgNFromSlider();
    });
  $("cg-n") &&
    $("cg-n").addEventListener("change", function () {
      loadCallGraph();
    });
}

async function regenerate() {
  // 强制重新采样（更新数据 + 排行），再用新缓存重新生成分层 SVG
  await loadCallGraph(true);
  await loadCallGraphSvg();
}

function setupSvgPanZoom(view) {
  var svg = view.querySelector("svg");
  if (!svg) return;
  var st = { tx: 0, ty: 0, s: 1, dragging: false, sx: 0, sy: 0, stx: 0, sty: 0 };
  function apply() {
    svg.style.transform = "translate(" + st.tx + "px," + st.ty + "px) scale(" + st.s + ")";
  }
  function clamp(v, a, b) {
    return Math.max(a, Math.min(b, v));
  }
  function fit() {
    var vb = svg.viewBox.baseVal;
    var pad = 24;
    var r = view.getBoundingClientRect();
    var s = Math.min((r.width - pad * 2) / vb.width, (r.height - pad * 2) / vb.height);
    st.s = clamp(s, 0.1, 8);
    st.tx = (r.width - vb.width * st.s) / 2;
    st.ty = (r.height - vb.height * st.s) / 2;
    apply();
  }
  view.addEventListener(
    "wheel",
    function (e) {
      e.preventDefault();
      var r = view.getBoundingClientRect();
      var mx = e.clientX - r.left;
      var my = e.clientY - r.top;
      var factor = Math.pow(1.0015, -e.deltaY);
      var ns = clamp(st.s * factor, 0.1, 8);
      var wx = (mx - st.tx) / st.s;
      var wy = (my - st.ty) / st.s;
      st.s = ns;
      st.tx = mx - wx * ns;
      st.ty = my - wy * ns;
      apply();
    },
    { passive: false }
  );
  view.addEventListener("mousedown", function (e) {
    if (e.button !== 0) return;
    st.dragging = true;
    st.sx = e.clientX;
    st.sy = e.clientY;
    st.stx = st.tx;
    st.sty = st.ty;
    view.classList.add("dragging");
    e.preventDefault();
  });
  window.addEventListener("mousemove", function (e) {
    if (!st.dragging) return;
    st.tx = st.stx + (e.clientX - st.sx);
    st.ty = st.sty + (e.clientY - st.sy);
    apply();
  });
  window.addEventListener("mouseup", function () {
    st.dragging = false;
    view.classList.remove("dragging");
  });
  view.addEventListener("dblclick", fit);
  fit();
}

export async function loadCallGraphSvg(force) {
  if (state.cgSvgLoading) return;
  var st = $("cg-svg-status");
  var btn = $("cg-regenerate");
  state.cgSvgLoading = true;
  if (btn) btn.disabled = true;
  if (st) { st.textContent = "生成中…"; st.style.color = "#f59e0b"; }
  var n = cgNFromSlider();
  try {
    var layerEl = $("cg-layer");
    var layer = layerEl ? layerEl.value : "";
    var q = "?n=" + n + "&floor_pct=" + state.cgThreshold;
    if (layer !== "") q += layer === "core" ? "&layer=-1" : "&layer=" + layer;
    q += force ? "&refresh=1" : "";
    var data = await api(API.callGraphSvg + q);
    if (!data || !data.ok) throw new Error((data && data.error) || "无数据");
    state.cgSvg = data;
    var view = $("cg-svg-view");
    if (view) {
      view.innerHTML = data.svg;
      view.classList.remove("hidden");
      setupSvgPanZoom(view);
      applyCgHighlight(); // 若已有选中函数，在 SVG 图上高亮
    }
    if (st) { st.textContent = "就绪 · " + data.meta.n_nodes + " 节点"; st.style.color = "#86efac"; }
  } catch (err) {
    if (st) { st.textContent = "失败: " + (err.message || String(err)); st.style.color = "#ef4444"; }
  } finally {
    state.cgSvgLoading = false;
    if (btn) btn.disabled = false;
  }
}

export function bindCallGraph() {
  bindCgControls();
  bindCgRankButtons();
  bindCgChainControls();
}
