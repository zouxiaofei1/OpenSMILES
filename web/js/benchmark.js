/* Benchmark Preview page: load/refresh, poll generation, filter + render table. */
import { state, api, API } from "./core.js";

function bm$(id) {
  return document.getElementById("bm-" + id);
}

function bmEsc(s) {
  return String(s || "")
    .replace(/&/g, "&amp;")
    .replace(/</g, "&lt;")
    .replace(/>/g, "&gt;")
    .replace(/"/g, "&quot;");
}

function bmDataFile() {
  var sel = bm$("data");
  return (sel && sel.value) || "merged_benchmark.json";
}

function bmPreviewUrl(apiUrl) {
  return apiUrl + "?data_file=" + encodeURIComponent(bmDataFile());
}

async function loadBmDatasets() {
  var sel = bm$("data");
  if (!sel) return;
  try {
    var data = await api(API.benchmarkDatasets);
    if (data && data.ok && Array.isArray(data.datasets)) {
      sel.innerHTML = data.datasets.map(function (name) {
        var selAttr = name === "merged_benchmark.json" ? " selected" : "";
        return '<option value="' + name + '"' + selAttr + ">" + name + "</option>";
      }).join("");
    } else {
      sel.innerHTML = '<option value="">无可用测试文件</option>';
    }
  } catch (_) {
    sel.innerHTML = '<option value="">加载失败</option>';
  }
}

function showBmSkeleton() {
  var tbody = bm$("tbody");
  if (!tbody) return;
  var html = "";
  for (var i = 0; i < 8; i++) {
    html +=
      '<tr class="skeleton-row">' +
      '<td><div class="skeleton skeleton-cell" style="width:80%"></div></td>' +
      '<td><div class="skeleton skeleton-cell" style="width:90%"></div><div class="skeleton skeleton-cell" style="width:60%;margin-top:6px"></div></td>' +
      '<td><div class="skeleton skeleton-cell" style="width:85%"></div><div class="skeleton skeleton-cell" style="width:55%;margin-top:6px"></div></td>' +
      '<td><div class="skeleton skeleton-cell" style="width:100%;height:126px"></div></td>' +
      "</tr>";
  }
  tbody.innerHTML = html;
}

export async function loadBenchmark() {
  var loading = bm$("loading");
  var empty = bm$("empty");
  if (loading) { loading.hidden = false; loading.style.display = ""; }
  if (empty) empty.hidden = true;
  showBmSkeleton();
  try {
    var data = await api(bmPreviewUrl(API.benchmarkPreview));
    state.bmRows = (data && data.rows) || [];
    state.bmLoaded = true;
    state.bmPage = 1;
    // Track generation state from response
    state.bmGenerating = !!(data && data.generating);
    state.bmGenDone = (data && data.gen_done) || 0;
    state.bmGenTotal = (data && data.gen_total) || 0;
    if (loading) { loading.hidden = true; loading.style.display = "none"; }
    renderBenchmark();
    // Auto-poll if generation is running or data is stale
    if (state.bmGenerating || (data && data.stale && state.bmRows.length === 0)) {
      startBmPolling();
    }
    // Auto-start generation if no cache at all
    if (!state.bmGenerating && state.bmRows.length === 0 && (data && data.source_total > 0)) {
      refreshBenchmark();
    }
  } catch (err) {
    if (loading) { loading.hidden = true; loading.style.display = "none"; }
    if (empty) {
      empty.hidden = false;
      empty.textContent = "加载失败: " + (err.message || String(err));
    }
    if (bm$("tbody")) bm$("tbody").innerHTML = "";
  }
}

function startBmPolling() {
  if (state.bmPollTimer) return;
  state.bmPollTimer = setInterval(async function () {
    try {
      var status = await api(bmPreviewUrl(API.benchmarkStatus));
      state.bmGenerating = !!(status && status.running);
      state.bmGenDone = (status && status.done) || 0;
      state.bmGenTotal = (status && status.total) || 0;
      if (!state.bmGenerating) {
        // Generation finished — reload data
        stopBmPolling();
        var data = await api(bmPreviewUrl(API.benchmarkPreview));
        state.bmRows = (data && data.rows) || [];
        state.bmGenDone = state.bmGenTotal;
        renderBenchmark();
        return;
      }
      // Still generating — re-render to update progress bar
      renderBenchmark();
    } catch (_) {
      /* ignore poll errors */
    }
  }, 2000);
}

export function stopBmPolling() {
  if (state.bmPollTimer) {
    clearInterval(state.bmPollTimer);
    state.bmPollTimer = null;
  }
  state.bmGenerating = false;
}

async function refreshBenchmark() {
  var btn = bm$("refresh");
  var name = bmDataFile();
  if (btn) btn.disabled = true;
  stopBmPolling();
  try {
    var data = await api(API.benchmarkRefresh + "?data_file=" + encodeURIComponent(name), { method: "POST" });
    if (data && data.ok) {
      state.bmGenerating = true;
      state.bmGenDone = 0;
      state.bmGenTotal = data.total || 0;
      // Reload to pick up any partial cache
      var preview = await api(API.benchmarkPreview + "?data_file=" + encodeURIComponent(name));
      state.bmRows = (preview && preview.rows) || [];
      renderBenchmark();
      startBmPolling();
    } else {
      var msg = (data && data.error) || "未知错误";
      if (/busy|in progress/i.test(msg)) {
        // A historical generation for another commit is running; retry shortly.
        setTimeout(refreshBenchmark, 2000);
        return;
      }
      alert("刷新失败: " + msg);
    }
  } catch (err) {
    alert("刷新请求失败: " + (err.message || String(err)));
  } finally {
    if (btn) btn.disabled = false;
  }
}

function bmFiltered() {
  var q = (bm$("q") && bm$("q").value || "").trim().toLowerCase();
  var f = bm$("filter") ? bm$("filter").value : "all";
  var out = [];
  for (var i = 0; i < state.bmRows.length; i++) {
    var r = state.bmRows[i];
    if (f === "ok" && !r.ok) continue;
    if (f === "fail" && r.ok) continue;
    if (f === "wrong" && !(r.ret && !r.ok)) continue;
    if (q) {
      var hay = [r.s, r.en, r.zh, r.ge, r.gz].join(" ").toLowerCase();
      if (hay.indexOf(q) < 0) continue;
    }
    out.push({ r: r, abs: i });
  }
  return out;
}

function bmBadge(r) {
  if (r.ok) return ' <span class="bm-badge bm-ok">match</span>';
  if (r.ret) return ' <span class="bm-badge bm-fail">miss</span><span class="bm-badge bm-ret">returned</span>';
  return ' <span class="bm-badge bm-fail">miss</span>';
}

function bmDrawAll() {
  if (typeof SmilesDrawer === "undefined" || !SmilesDrawer.SvgDrawer) {
    var wraps = document.querySelectorAll("#bm-tbody .bm-struct-wrap");
    for (var i = 0; i < wraps.length; i++) {
      wraps[i].innerHTML = '<span class="bm-struct-err">SmilesDrawer 未加载</span>';
    }
    return;
  }
  var opts = { width: 200, height: 126, bondThickness: 1.1, padding: 6 };
  var drawer = new SmilesDrawer.SvgDrawer(opts);
  var svgs = document.querySelectorAll("#bm-tbody svg[data-smiles]");
  for (var i = 0; i < svgs.length; i++) {
    (function (svg) {
      var smi = svg.getAttribute("data-smiles") || "";
      if (!smi) {
        svg.parentNode.innerHTML = '<span class="bm-struct-err">空 SMILES</span>';
        return;
      }
      SmilesDrawer.parse(
        smi,
        function (tree) {
          try {
            drawer.draw(tree, svg, "light");
          } catch (e) {
            svg.parentNode.innerHTML = '<span class="bm-struct-err">draw err</span>';
          }
        },
        function () {
          svg.parentNode.innerHTML = '<span class="bm-struct-err">parse fail</span>';
        }
      );
    })(svgs[i]);
  }
}

function renderBenchmark() {
  var items = bmFiltered();
  var ps = Math.max(1, parseInt(bm$("pageSize") ? bm$("pageSize").value : "100", 10) || 100);
  var pages = Math.max(1, Math.ceil(items.length / ps));
  if (state.bmPage > pages) state.bmPage = pages;
  if (state.bmPage < 1) state.bmPage = 1;
  var start = (state.bmPage - 1) * ps;
  var slice = items.slice(start, start + ps);

  // Stats
  if (bm$("total")) bm$("total").textContent = String(state.bmRows.length);
  if (bm$("shown")) bm$("shown").textContent = String(items.length);
  if (bm$("okn")) bm$("okn").textContent = String(state.bmRows.filter(function (r) { return r.ok; }).length);
  if (bm$("enn")) bm$("enn").textContent = String(state.bmRows.filter(function (r) { return r.en_ok; }).length);
  if (bm$("zhn")) bm$("zhn").textContent = String(state.bmRows.filter(function (r) { return r.zh_ok; }).length);
  if (bm$("page")) bm$("page").textContent = String(state.bmPage);
  if (bm$("pages")) bm$("pages").textContent = String(pages);
  if (bm$("prev")) bm$("prev").disabled = state.bmPage <= 1;
  if (bm$("next")) bm$("next").disabled = state.bmPage >= pages;

  // Generation progress
  var progEl = bm$("gen-progress");
  if (state.bmGenerating && state.bmGenTotal > 0) {
    if (!progEl) {
      progEl = document.createElement("span");
      progEl.id = "bm-gen-progress";
      progEl.className = "chip";
      var controls = document.querySelector(".benchmark-controls");
      if (controls) controls.appendChild(progEl);
    }
    var pct = Math.round(state.bmGenDone / state.bmGenTotal * 100);
    progEl.innerHTML = '生成中 <b>' + state.bmGenDone + '/' + state.bmGenTotal + '</b> (' + pct + '%)';
    progEl.style.color = '#f59e0b';
  } else if (progEl) {
    if (state.bmGenDone >= state.bmGenTotal && state.bmGenTotal > 0) {
      progEl.innerHTML = '生成完成 <b>' + state.bmGenTotal + '</b> 条';
      progEl.style.color = '#86efac';
    } else {
      progEl.remove();
    }
  }

  var tbody = bm$("tbody");
  var empty = bm$("empty");
  if (!tbody) return;

  if (!slice.length) {
    tbody.innerHTML = "";
    if (empty) empty.hidden = false;
    return;
  }
  if (empty) empty.hidden = true;

  var html = "";
  for (var j = 0; j < slice.length; j++) {
    var abs = slice[j].abs;
    var r = slice[j].r;
    var cls = r.ok ? "bm-match" : "bm-miss";
    var predEn = bmEsc(r.en) || "<span style='color:#94a3b8'>(空)</span>";
    var predZh = bmEsc(r.zh);
    var goldEn = bmEsc(r.ge) || "<span style='color:#94a3b8'>(无英标)</span>";
    var goldZh = bmEsc(r.gz);
    html +=
      '<tr class="' + cls + '">' +
      '<td><span class="bm-idx">#' + (abs + 1) + '</span><span class="bm-smiles">' + bmEsc(r.s) + "</span></td>" +
      '<td><div class="bm-name-en">' + predEn + bmBadge(r) + "</div>" +
      (predZh ? '<div class="bm-name-zh">' + predZh + "</div>" : "") + "</td>" +
      '<td><div class="bm-gold-en">' + goldEn + "</div>" +
      (goldZh ? '<div class="bm-gold-zh">' + goldZh + "</div>" : "") + "</td>" +
      '<td><div class="bm-struct-wrap"><svg data-smiles="' + bmEsc(r.s) + '"></svg></div></td>' +
      "</tr>";
  }
  tbody.innerHTML = html;
  requestAnimationFrame(function () { bmDrawAll(); });
}

export function bindBenchmark() {
  // 数据文件下拉：加载可用列表；切换时重载对应文件
  if (bm$("data")) {
    loadBmDatasets();
    bm$("data").addEventListener("change", function () {
      stopBmPolling();
      state.bmRows = [];
      state.bmPage = 1;
      state.bmLoaded = false;
      state.bmGenerating = false;
      loadBenchmark();
    });
  }
  if (bm$("q")) {
    var bmTimer = null;
    bm$("q").addEventListener("input", function () {
      clearTimeout(bmTimer);
      bmTimer = setTimeout(function () { state.bmPage = 1; renderBenchmark(); }, 150);
    });
  }
  if (bm$("filter")) {
    bm$("filter").addEventListener("change", function () { state.bmPage = 1; renderBenchmark(); });
  }
  if (bm$("pageSize")) {
    bm$("pageSize").addEventListener("change", function () { state.bmPage = 1; renderBenchmark(); });
  }
  if (bm$("prev")) {
    bm$("prev").addEventListener("click", function () { state.bmPage -= 1; renderBenchmark(); });
  }
  if (bm$("next")) {
    bm$("next").addEventListener("click", function () { state.bmPage += 1; renderBenchmark(); });
  }
  if (bm$("refresh")) {
    bm$("refresh").addEventListener("click", refreshBenchmark);
  }
}
