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

/* 预览响应的唯一落地点。有三条路径会拿到新的 rows(首载 / 生成完成 / 刷新后取回),
   必须都走这里: 少更新一次 bmFiles, 行内 f 就会解析成 undefined, 文件排除会静默失效
   (Set.has(undefined) 恒为 false, 不报错)。 */
function bmApplyData(data, resetPage) {
  state.bmRows = (data && data.rows) || [];
  state.bmFiles = (data && data.files) || [];
  state.bmLoaded = true;
  if (resetPage) state.bmPage = 1;
  state.bmGenerating = !!(data && data.generating);
  state.bmGenDone = (data && data.gen_done) || 0;
  state.bmGenTotal = (data && data.gen_total) || 0;
  bmRenderSide();
}

export async function loadBenchmark() {
  var loading = bm$("loading");
  var empty = bm$("empty");
  if (loading) { loading.hidden = false; loading.style.display = ""; }
  if (empty) empty.hidden = true;
  showBmSkeleton();
  try {
    var data = await api(bmPreviewUrl(API.benchmarkPreview));
    bmApplyData(data, true);
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
        bmApplyData(data, false);
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
      bmApplyData(preview, false);
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

/* 排序键取值: 相似度缺失(该语言无金标名)返回 NaN, 由排序器放到最后。 */
function bmSortValue(r, key) {
  var v = r[key];
  return typeof v === "number" ? v : NaN;
}

/* 按 bm-sort 排序 items({r, abs}); 空值=源序, 同值回落到源序保证稳定。 */
function bmSortItems(items) {
  var spec = bm$("sort") ? bm$("sort").value : "";
  if (!spec) return items;
  var cut = spec.lastIndexOf("-");
  var key = spec.slice(0, cut);
  var dir = spec.slice(cut + 1) === "desc" ? -1 : 1;
  return items.slice().sort(function (a, b) {
    var va = bmSortValue(a.r, key);
    var vb = bmSortValue(b.r, key);
    var na = isNaN(va);
    var nb = isNaN(vb);
    if (na || nb) {
      // 缺失值恒排最后, 不随升降序翻转(升序看最差样本时不该被"无金标"占满首页)
      if (na && nb) return a.abs - b.abs;
      return na ? 1 : -1;
    }
    if (va === vb) return a.abs - b.abs;
    return (va - vb) * dir;
  });
}

/* 重原子数筛选: 双滑块取闭区间 [lo, hi]; 无 hac 字段的行(旧缓存)不受限。 */
function bmHacRange() {
  var loEl = bm$("hac-lo");
  var hiEl = bm$("hac-hi");
  if (!loEl || !hiEl || loEl.disabled) return null;
  return { lo: parseInt(loEl.value, 10) || 0, hi: parseInt(hiEl.value, 10) || 0 };
}

/* 行的重原子数域: {lo, hi}; 无任何数值型 hac 时返回 null。 */
function bmHacBounds() {
  var lo = Infinity;
  var hi = -Infinity;
  var n = 0;
  for (var i = 0; i < state.bmRows.length; i++) {
    var v = state.bmRows[i].hac;
    if (typeof v !== "number") continue;
    n += 1;
    if (v < lo) lo = v;
    if (v > hi) hi = v;
  }
  return n ? { lo: lo, hi: hi } : null;
}

function bmHacLabel(lo, hi) {
  var loEl = bm$("hac-lo");
  var hiEl = bm$("hac-hi");
  if (lo == null) {
    lo = loEl ? parseInt(loEl.value, 10) || 0 : 0;
    hi = hiEl ? parseInt(hiEl.value, 10) || 0 : 0;
  }
  if (bm$("hac-lo-v")) bm$("hac-lo-v").textContent = String(lo);
  if (bm$("hac-hi-v")) bm$("hac-hi-v").textContent = String(hi);
}

/* 默认取用的重原子数区间; 滑块可拖范围仍是本批数据的 [最小, 最大]。 */
const BM_HAC_DEFAULT = { lo: 1, hi: 25 };

/* 数据换了(重新加载/切文件)后重置滑块: 域=本批 hac 的 [最小, 最大], 取值=默认区间夹进域内。 */
function bmHacReset() {
  var loEl = bm$("hac-lo");
  var hiEl = bm$("hac-hi");
  if (!loEl || !hiEl) return;
  var b = bmHacBounds();
  if (!b) {
    loEl.disabled = true;
    hiEl.disabled = true;
    loEl.min = hiEl.min = loEl.max = hiEl.max = "0";
    loEl.value = hiEl.value = "0";
    bmHacLabel(0, 0);
    return;
  }
  var lo = Math.min(Math.max(BM_HAC_DEFAULT.lo, b.lo), b.hi);
  var hi = Math.min(Math.max(BM_HAC_DEFAULT.hi, b.lo), b.hi);
  loEl.disabled = false;
  hiEl.disabled = false;
  loEl.min = hiEl.min = String(b.lo);
  loEl.max = hiEl.max = String(b.hi);
  loEl.value = String(lo);
  hiEl.value = String(hi);
  bmHacLabel(lo, hi);
}

/* 拖动一端越过另一端时把另一端顶开, 保持 lo <= hi。 */
function bmHacInput(which) {
  var loEl = bm$("hac-lo");
  var hiEl = bm$("hac-hi");
  if (!loEl || !hiEl) return;
  var lo = parseInt(loEl.value, 10) || 0;
  var hi = parseInt(hiEl.value, 10) || 0;
  if (which === "lo" && lo > hi) {
    hi = lo;
    hiEl.value = String(hi);
  } else if (which === "hi" && hi < lo) {
    lo = hi;
    loEl.value = String(lo);
  }
  bmHacLabel(lo, hi);
  state.bmPage = 1;
  renderBenchmark();
}

/* ── 左栏排除面板 ────────────────────────────────────────────────
   勾中的文件/特征不出现在表里。排除只活在 bmFiltered() 里, 绝不改写 state.bmRows ——
   一旦过滤掉源数组, 重原子滑块区间会重置、「共 N 条」失真、排除也无法撤销。 */

const BM_FEAT_TOP = 30; // 默认列出的特征条数, 其余靠「展开全部」

/* 特征 → 在 bmRows 里的出现次数。 */
function bmFeatCounts() {
  var counts = {};
  for (var i = 0; i < state.bmRows.length; i++) {
    var feats = state.bmRows[i].feat || [];
    for (var k = 0; k < feats.length; k++) {
      counts[feats[k]] = (counts[feats[k]] || 0) + 1;
    }
  }
  return counts;
}

/* 一条勾选项的 HTML; kind 决定写进哪个排除集合。 */
function bmSideItem(kind, name, count, extraHtml) {
  var checked = (kind === "file" ? state.bmExFiles : state.bmExFeats).has(name);
  return '<li><label class="bm-side-item" title="' + bmEsc(name) + '">' +
    '<input type="checkbox" data-kind="' + kind + '" data-name="' + bmEsc(name) + '"' +
    (checked ? " checked" : "") + " />" +
    '<span class="bm-side-name">' + bmEsc(name) + (extraHtml || "") + "</span>" +
    '<span class="bm-side-cnt">' + count + "</span></label></li>";
}

/* 文件列表: 计数在 bmRows 上数, 不能用过滤后的行 —— 否则计数被它自己控制的排除集合影响。
   多文件归属是常态(merged 的行同时属于若干子集), 所以各项计数之和会大于总行数。 */
function bmRenderFiles() {
  var ul = bm$("file-list");
  if (!ul) return;
  var counts = {};
  var i, k;
  for (i = 0; i < state.bmRows.length; i++) {
    var f = state.bmRows[i].f || [];
    for (k = 0; k < f.length; k++) counts[f[k]] = (counts[f[k]] || 0) + 1;
  }
  var cur = bmDataFile();
  var html = "";
  for (i = 0; i < state.bmFiles.length; i++) {
    var name = state.bmFiles[i];
    // 勾选当前数据文件必然清空列表(每行的 f 都含它), 标出来免得用户以为坏了
    html += bmSideItem("file", name, counts[i] || 0, name === cur ? "<em>当前</em>" : "");
  }
  ul.innerHTML = html || '<li class="bm-side-none">无数据集文件</li>';
  if (bm$("file-n")) bm$("file-n").textContent = String(state.bmFiles.length);
}

/* 特征列表: 按出现次数降序(同次数再按名, 否则每次重绘顺序都会跳), 默认只列前 BM_FEAT_TOP 个,
   展开或搜索时列全部。已勾选项去重后置顶恒可见 —— 否则勾了长尾再收起会「看不见但仍在生效」。 */
function bmRenderFeats() {
  var ul = bm$("feat-list");
  if (!ul) return;
  var counts = bmFeatCounts();
  // 勾了但当前数据里没有的特征也列出来(计数 0), 否则用户没地方取消它
  state.bmExFeats.forEach(function (n) {
    if (!(n in counts)) counts[n] = 0;
  });
  var names = Object.keys(counts);
  var q = (state.bmFeatQuery || "").trim().toLowerCase();
  if (q) {
    names = names.filter(function (n) { return n.toLowerCase().indexOf(q) >= 0; });
  }
  names.sort(function (a, b) {
    return counts[b] - counts[a] || (a < b ? -1 : a > b ? 1 : 0);
  });

  var pinned = [];
  var rest = [];
  for (var i = 0; i < names.length; i++) {
    (state.bmExFeats.has(names[i]) ? pinned : rest).push(names[i]);
  }
  var expanded = state.bmFeatExpanded || !!q;
  var html = "";
  pinned.concat(expanded ? rest : rest.slice(0, BM_FEAT_TOP)).forEach(function (n) {
    html += bmSideItem("feat", n, counts[n]);
  });
  ul.innerHTML = html || '<li class="bm-side-none">' +
    (state.bmRows.length ? "无匹配特征" : "暂无数据") + "</li>";
  if (bm$("feat-n")) bm$("feat-n").textContent = String(Object.keys(counts).length);
  var more = bm$("feat-more");
  if (more) {
    more.hidden = !!q || rest.length <= BM_FEAT_TOP;
    more.textContent = state.bmFeatExpanded ? "收起" : "展开全部";
  }
}

/* 只重绘两个列表的 innerHTML: 搜索框 #bm-feat-q 是 index.html 里的静态节点, 不能被这里碰到
   (否则每敲一个字就丢焦点)。也不要把它挂到 renderBenchmark() 上 —— 那个函数在输入、拖滑块、
   生成轮询(每 2s)时都会跑, 重绘会把侧栏滚动位置打回顶部。 */
function bmRenderSide() {
  var side = bm$("side");
  if (side) side.classList.toggle("collapsed", !state.bmSideOpen);
  var toggle = bm$("side-toggle");
  if (toggle) {
    toggle.classList.toggle("is-on", !!state.bmSideOpen);
    toggle.setAttribute("aria-expanded", state.bmSideOpen ? "true" : "false");
  }
  bmRenderFiles();
  bmRenderFeats();
}

function bmFiltered() {
  var q = (bm$("q") && bm$("q").value || "").trim().toLowerCase();
  var f = bm$("filter") ? bm$("filter").value : "all";
  var hac = bmHacRange();
  // 预取引用: 这个函数被输入/拖滑块/2s 轮询高频调用, 循环里别再取属性、也别造闭包。
  var exF = state.bmExFiles;
  var exT = state.bmExFeats;
  var files = state.bmFiles;
  var useF = exF.size > 0;
  var useT = exT.size > 0;
  var exN = 0;
  var out = [];
  for (var i = 0; i < state.bmRows.length; i++) {
    var r = state.bmRows[i];
    // 排除: 按文件归属(行内 f 是 state.bmFiles 的下标) / 按特征
    if (useF) {
      var rf = r.f || [];
      var hitF = false;
      for (var k = 0; k < rf.length; k++) {
        if (exF.has(files[rf[k]])) { hitF = true; break; }
      }
      if (hitF) { exN += 1; continue; }
    }
    if (useT) {
      var rt = r.feat || [];
      var hitT = false;
      for (var m = 0; m < rt.length; m++) {
        if (exT.has(rt[m])) { hitT = true; break; }
      }
      if (hitT) { exN += 1; continue; }
    }
    if (f === "ok" && !r.ok) continue;
    if (f === "fail" && r.ok) continue;
    if (hac && typeof r.hac === "number" && (r.hac < hac.lo || r.hac > hac.hi)) continue;
    if (q) {
      var hay = [r.s, r.en, r.zh, r.ge, r.gz].join(" ").toLowerCase();
      if (hay.indexOf(q) < 0) continue;
    }
    out.push({ r: r, abs: i });
  }
  state.bmExCount = exN;
  return bmSortItems(out);
}

/* 相似度/复杂度徽标: 预测名 ↔ 正确名的字符级相似度, 无金标名的语言不显示。
   分级配色与 Namer 页相似度条同色系(高绿 / 中琥珀 / 低红)。 */
function bmSimHtml(r) {
  var parts = [];
  [["EN", r.sim_en], ["ZH", r.sim_zh]].forEach(function (pair) {
    var v = pair[1];
    if (typeof v !== "number") return;
    var pct = Math.max(0, Math.min(100, Math.round(v * 100)));
    var cls = pct >= 90 ? "bm-sim-hi" : pct >= 70 ? "bm-sim-mid" : "bm-sim-lo";
    parts.push(
      '<span class="bm-sim ' + cls + '" title="' + pair[0] + " 相似度 " + pct + '%">' +
        pair[0] + " " + pct + "%</span>"
    );
  });
  if (typeof r.cx === "number") {
    parts.push('<span class="bm-cx" title="BertzCT 复杂度">cx ' + Math.round(r.cx) + "</span>");
  }
  return parts.length ? '<div class="bm-sims">' + parts.join("") + "</div>" : "";
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
  // 行数组换了(重载/切数据文件/生成完成后重取) → 滑块域跟着换; 单纯重渲染不动用户选择
  if (state.bmHacRows !== state.bmRows) {
    state.bmHacRows = state.bmRows;
    bmHacReset();
  }
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
  if (bm$("fusedn")) bm$("fusedn").textContent = String(state.bmRows.filter(function (r) { return r.fused; }).length);
  // 其余 chips 都是全量口径, 只有「显示」看排除 —— 用这一项说明差额从哪来
  if (bm$("exn")) bm$("exn").textContent = String(state.bmExCount || 0);
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
      '<td><span class="bm-idx">#' + (abs + 1) + '</span><span class="bm-smiles">' + bmEsc(r.s) + "</span>" +
      bmSimHtml(r) + "</td>" +
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
  if (bm$("sort")) {
    bm$("sort").addEventListener("change", function () { state.bmPage = 1; renderBenchmark(); });
  }
  if (bm$("pageSize")) {
    bm$("pageSize").addEventListener("change", function () { state.bmPage = 1; renderBenchmark(); });
  }
  if (bm$("hac-lo")) {
    bm$("hac-lo").addEventListener("input", function () { bmHacInput("lo"); });
  }
  if (bm$("hac-hi")) {
    bm$("hac-hi").addEventListener("input", function () { bmHacInput("hi"); });
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

  // ── 左栏排除面板 ──
  if (bm$("side-toggle")) {
    bm$("side-toggle").addEventListener("click", function () {
      state.bmSideOpen = !state.bmSideOpen;
      bmRenderSide();
    });
  }
  if (bm$("ex-clear")) {
    bm$("ex-clear").addEventListener("click", function () {
      state.bmExFiles.clear();
      state.bmExFeats.clear();
      state.bmPage = 1;
      bmRenderSide();
      renderBenchmark();
    });
  }
  if (bm$("feat-more")) {
    bm$("feat-more").addEventListener("click", function () {
      state.bmFeatExpanded = !state.bmFeatExpanded;
      bmRenderSide();
    });
  }
  if (bm$("feat-q")) {
    var featTimer = null;
    bm$("feat-q").addEventListener("input", function () {
      clearTimeout(featTimer);
      featTimer = setTimeout(function () {
        state.bmFeatQuery = bm$("feat-q").value || "";
        // 只重绘左栏: 特征搜索不改变哪些行会被显示, 表格不用重算
        bmRenderSide();
      }, 150);
    });
  }
  // 勾选/取消: 浏览器已经改好了复选框本身, 这里只更新数据与表格 —— 不重绘列表,
  // 否则用户勾第二项时侧栏会跳回顶部。
  bm$("file-list") && bm$("file-list").addEventListener("change", onExcludeChange);
  bm$("feat-list") && bm$("feat-list").addEventListener("change", onExcludeChange);
}

/* 左栏复选框变更(事件委托): 增删对应排除集合后重渲染表格。 */
function onExcludeChange(ev) {
  var el = ev.target;
  if (!el || el.tagName !== "INPUT") return;
  var name = el.getAttribute("data-name") || "";
  var set = el.getAttribute("data-kind") === "file" ? state.bmExFiles : state.bmExFeats;
  if (el.checked) set.add(name);
  else set.delete(name);
  state.bmPage = 1;
  renderBenchmark();
}
