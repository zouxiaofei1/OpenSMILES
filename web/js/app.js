/* ChemAgent Namer — client */
(function () {
  "use strict";

  const API = {
    name: "/api/v1/name",
    benchmarkPreview: "/api/v1/benchmark-preview",
    benchmarkRefresh: "/api/v1/benchmark-preview/refresh",
    benchmarkStatus: "/api/v1/benchmark-preview/status",
    benchmarkRun: "/api/v1/benchmark-run",
    benchmarkRunStatus: "/api/v1/benchmark-run/status",
    benchmarkRunResult: "/api/v1/benchmark-run/result",
    layerSample: "/api/v1/layer-benchmark/sample",
    layerSampleStatus: "/api/v1/layer-benchmark/sample-status",
    layerData: "/api/v1/layer-benchmark/data",
    layerScore: "/api/v1/layer-benchmark/score",
    layerScoreResult: "/api/v1/layer-benchmark/score-result",
    layerAnalyzeOne: "/api/v1/layer-benchmark/analyze-one",
    codeAnalysis: "/api/v1/code-analysis",
    callGraph: "/api/v1/call-graph",
    callGraphSvg: "/api/v1/call-graph/svg",
    debug: "/api/v1/debug",
    settings: "/api/v1/settings",
    wiki: "/api/v1/wiki",
    wikiDoc: "/api/v1/wiki/doc",
  };

  const LIVE_NAME_DEBOUNCE_MS = 20;

  const state = {
    namerHistory: [],
    liveNameEnabled: true,
    ketcherReady: false,
    lastNamedSmiles: "",
    nameReqSeq: 0,
    ketcherBridge: null,
    liveDebounceTimer: null,
    liveSmilesTimer: null,
    suppressSmilesLive: false,
    // benchmark
    currentPage: "namer",
    bmRows: [],
    bmLoaded: false,
    bmPage: 1,
    bmGenerating: false,
    bmGenDone: 0,
    bmGenTotal: 0,
    bmGenTotal: 0,
    bmPollTimer: null,
    // benchmark run
    brRunning: false,
    brPollTimer: null,
    // layer benchmark
    lbLayer: 0,
    lbGenerating: false,
    lbPollTimer: null,
    lbScore: null,
    // code analysis
    caData: null,
    // call graph
    cgData: null,
    cgThreshold: 0,
    cgLoading: false,
    cgSvg: null,
    cgSvgLoading: false,
    cgSourceTotal: 0,
    wikiCurrent: "",
    cgRank: {
      calls: { asc: false, all: false },
      time: { asc: false, all: false },
    },
  };

  const $ = (id) => document.getElementById(id);

  function escapeHtml(str) {
    return String(str)
      .replace(/&/g, "&amp;")
      .replace(/</g, "&lt;")
      .replace(/>/g, "&gt;")
      .replace(/"/g, "&quot;");
  }

  async function api(url, opts) {
    const res = await fetch(url, {
      headers: { Accept: "application/json", ...(opts && opts.body ? { "Content-Type": "application/json" } : {}) },
      ...opts,
    });
    if (!res.ok) {
      let detail = res.statusText;
      try {
        const j = await res.json();
        detail = j.detail || JSON.stringify(j);
      } catch (_) {
        /* ignore */
      }
      throw new Error(typeof detail === "string" ? detail : JSON.stringify(detail));
    }
    if (res.status === 204) return null;
    return res.json();
  }

  /* ---------- Namer ---------- */

  function setNamerError(msg) {
    const errEl = $("namer-error");
    if (!errEl) return;
    if (!msg) {
      errEl.hidden = true;
      errEl.textContent = "";
      return;
    }
    errEl.hidden = false;
    errEl.textContent = msg;
  }

  function setKetcherControlsEnabled(on) {
    ["btn-load-smiles", "btn-clear-ketcher"].forEach((id) => {
      const el = $(id);
      if (el) el.disabled = !on;
    });
  }

  function showKetcherFallback(show, msg) {
    const fb = $("ketcher-fallback");
    const frame = $("ketcher-frame");
    const msgEl = $("ketcher-fallback-msg");
    if (fb) {
      fb.hidden = !show;
      fb.style.display = show ? "flex" : "none";
    }
    if (frame) frame.style.visibility = show ? "hidden" : "visible";
    if (msgEl && msg) msgEl.textContent = msg;
  }

  function ensureKetcher(forceRetry) {
    if (!window.ChemNamerKetcher) {
      showKetcherFallback(
        true,
        "前端 bridge 未加载（/js/namer-ketcher.js）。请硬刷新页面。"
      );
      return;
    }
    if (forceRetry && state.ketcherBridge) {
      try {
        state.ketcherBridge.destroy();
      } catch (_) {}
      state.ketcherBridge = null;
      state.ketcherReady = false;
      setKetcherControlsEnabled(false);
    }
    if (state.ketcherBridge) {
      state.ketcherBridge.init();
      return;
    }
    const iframe = $("ketcher-frame");
    if (!iframe) return;
    showKetcherFallback(false);
    state.ketcherBridge = window.ChemNamerKetcher.createBridge({
      iframe,
      src: window.ChemNamerKetcher.DEFAULT_SRC,
      onReady: () => {
        state.ketcherReady = true;
        setKetcherControlsEnabled(true);
        showKetcherFallback(false);
      },
      onError: (err) => {
        state.ketcherReady = false;
        setKetcherControlsEnabled(false);
        const message =
          (err && err.message ? err.message : String(err)) || "Ketcher 加载失败";
        showKetcherFallback(true, message + "。文本 SMILES 命名仍可用。");
      },
      onChange: () => {
        scheduleLiveName();
      },
    });
    state.ketcherBridge.init();
  }

  async function resolveSmilesForName() {
    let fromEditor = "";
    if (state.ketcherBridge && state.ketcherBridge.isReady()) {
      fromEditor = await state.ketcherBridge.getSmiles();
    }
    if (fromEditor) {
      const input = $("smiles-input");
      if (input) input.value = fromEditor;
      return fromEditor;
    }
    return (($("smiles-input") && $("smiles-input").value) || "").trim();
  }

  async function runName(smiles, opts) {
    const fromLive = !!(opts && opts.fromLive);
    const CK = window.ChemNamerKetcher;
    const seq = CK ? CK.nextReqSeq(state.nameReqSeq) : state.nameReqSeq + 1;
    state.nameReqSeq = seq;
    setNamerError("");
    if (!smiles) {
      if (!fromLive) setNamerError("请绘制或输入 SMILES");
      return;
    }
    const btn = $("btn-name");
    if (!fromLive && btn) btn.disabled = true;
    try {
      const result = await api(API.name, {
        method: "POST",
        body: JSON.stringify({ smiles }),
      });
      if (seq !== state.nameReqSeq) return;
      showNamerResult(result);
      state.lastNamedSmiles = smiles;
      state.namerHistory.unshift({
        smiles,
        en: result.en,
        zh: result.zh,
        success: result.success,
        time_ms: result.time_ms,
      });
      state.namerHistory = state.namerHistory.slice(0, 20);
      renderNamerHistory();
    } catch (err) {
      if (seq !== state.nameReqSeq) return;
      setNamerError(err.message || String(err));
    } finally {
      if (!fromLive && btn) btn.disabled = false;
    }
  }

  function scheduleLiveName() {
    if (!state.liveNameEnabled) return;
    if (state.liveDebounceTimer) clearTimeout(state.liveDebounceTimer);
    state.liveDebounceTimer = setTimeout(async () => {
      state.liveDebounceTimer = null;
      if (!state.liveNameEnabled) return;
      let smiles = "";
      if (state.ketcherBridge && state.ketcherBridge.isReady()) {
        smiles = await state.ketcherBridge.getSmiles();
      }
      const CK = window.ChemNamerKetcher;
      if (CK && CK.shouldSkipLiveName(true, smiles, state.lastNamedSmiles)) return;
      if (smiles) {
        const input = $("smiles-input");
        if (input) {
          state.suppressSmilesLive = true;
          input.value = smiles;
        }
      }
      await runName(smiles, { fromLive: true });
    }, LIVE_NAME_DEBOUNCE_MS);
  }

  function scheduleLiveSmilesName() {
    if (!state.liveNameEnabled) return;
    if (state.suppressSmilesLive) {
      state.suppressSmilesLive = false;
      return;
    }
    if (state.liveSmilesTimer) clearTimeout(state.liveSmilesTimer);
    state.liveSmilesTimer = setTimeout(async () => {
      state.liveSmilesTimer = null;
      if (!state.liveNameEnabled) return;
      const smiles = (($("smiles-input") && $("smiles-input").value) || "").trim();
      if (!smiles) return;
      const CK = window.ChemNamerKetcher;
      if (CK && CK.shouldSkipLiveName(true, smiles, state.lastNamedSmiles)) return;
      // Load SMILES into Ketcher (onChange will debounce but shouldSkipLiveName will skip)
      if (state.ketcherBridge && state.ketcherBridge.isReady()) {
        try {
          await state.ketcherBridge.setMolecule(smiles);
        } catch (_) {
          /* ignore ketcher load errors during live input */
        }
      }
      // Run naming
      await runName(smiles, { fromLive: true });
    }, LIVE_NAME_DEBOUNCE_MS);
  }

  function renderNamerHistory() {
    const ul = $("namer-history");
    if (!ul) return;
    if (!state.namerHistory.length) {
      ul.innerHTML = '<li class="muted-text small">尚无命名记录</li>';
      return;
    }
    ul.innerHTML = state.namerHistory
      .map((h) => {
        const ok = h.success ? "ok" : "fail";
        return (
          `<li>` +
          `<span class="h-smiles">${escapeHtml(h.smiles)}</span>` +
          `<span class="badge ${ok}">${h.success ? "ok" : "fail"} · ${Math.round(h.time_ms || 0)}ms</span>` +
          `<span class="h-names">${escapeHtml(h.en || "—")} / ${escapeHtml(h.zh || "—")}</span>` +
          `</li>`
        );
      })
      .join("");
  }

  async function onName(ev) {
    ev.preventDefault();
    const smiles = await resolveSmilesForName();
    await runName(smiles, { fromLive: false });
  }

  function showNamerResult(result) {
    const box = $("namer-result");
    if (!box) return;
    // Force re-trigger reveal animation
    box.classList.remove("hidden");
    box.style.animation = "none";
    box.offsetHeight; // force reflow
    box.style.animation = "";
    const badge = $("namer-success");
    if (badge) {
      badge.textContent = result.success ? "success" : "failed";
      badge.className = "badge " + (result.success ? "ok" : "fail");
    }
    if ($("namer-time")) $("namer-time").textContent = Math.round(result.time_ms || 0) + " ms";
    if ($("namer-en")) $("namer-en").textContent = result.en || "—";
    if ($("namer-zh")) $("namer-zh").textContent = result.zh || "—";
    if ($("namer-source")) $("namer-source").textContent = result.source || "—";
  }

  /* ---------- Page switching ---------- */

  function switchPage(name) {
    state.currentPage = name;
    // Update tab states
    document.querySelectorAll(".tab").forEach(function (t) {
      var page = t.getAttribute("data-page");
      if (page === name) {
        t.setAttribute("aria-current", "page");
      } else {
        t.removeAttribute("aria-current");
      }
    });
    // Show/hide pages
    var namerLayout = document.getElementById("namer-layout");
    var bmPage = document.getElementById("benchmark-page");
    var brPage = document.getElementById("benchmark-run-page");
    var lbPage = document.getElementById("layer-benchmark-page");
    var caPage = document.getElementById("code-analysis-page");
    var cgPage = document.getElementById("call-graph-page");
    var debugPage = document.getElementById("debug-page");
    var settingsPage = document.getElementById("settings-page");
    var wikiPage = document.getElementById("wiki-page");
    // Hide all first
    if (namerLayout) namerLayout.classList.add("hidden");
    if (bmPage) bmPage.classList.add("hidden");
    if (brPage) brPage.classList.add("hidden");
    if (lbPage) lbPage.classList.add("hidden");
    if (caPage) caPage.classList.add("hidden");
    if (cgPage) cgPage.classList.add("hidden");
    if (debugPage) debugPage.classList.add("hidden");
    if (settingsPage) settingsPage.classList.add("hidden");
    if (wikiPage) wikiPage.classList.add("hidden");
    // Stop polling
    stopBmPolling();
    stopBrPolling();
    stopLbPolling();
    // Show active page
    if (name === "benchmark") {
      if (bmPage) bmPage.classList.remove("hidden");
      if (!state.bmLoaded) {
        loadBenchmark();
      }
    } else if (name === "benchmark-run") {
      if (brPage) brPage.classList.remove("hidden");
    } else if (name === "layer-benchmark") {
      if (lbPage) lbPage.classList.remove("hidden");
      renderLbPanel(state.lbLayer);
    } else if (name === "code-analysis") {
      if (caPage) caPage.classList.remove("hidden");
      loadCodeAnalysis();
    } else if (name === "call-graph") {
      if (cgPage) cgPage.classList.remove("hidden");
      loadCallGraph();
    } else if (name === "debug") {
      if (debugPage) debugPage.classList.remove("hidden");
    } else if (name === "settings") {
      if (settingsPage) settingsPage.classList.remove("hidden");
      loadSettings();
    } else if (name === "wiki") {
      if (wikiPage) wikiPage.classList.remove("hidden");
      loadWiki();
    } else {
      if (namerLayout) namerLayout.classList.remove("hidden");
    }
  }

  /* ---------- Benchmark Preview ---------- */

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

  async function loadBenchmark() {
    var loading = bm$("loading");
    var empty = bm$("empty");
    if (loading) { loading.hidden = false; loading.style.display = ""; }
    if (empty) empty.hidden = true;
    showBmSkeleton();
    try {
      var data = await api(API.benchmarkPreview);
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
        var status = await api(API.benchmarkStatus);
        state.bmGenerating = !!(status && status.running);
        state.bmGenDone = (status && status.done) || 0;
        state.bmGenTotal = (status && status.total) || 0;
        if (!state.bmGenerating) {
          // Generation finished — reload data
          stopBmPolling();
          var data = await api(API.benchmarkPreview);
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

  function stopBmPolling() {
    if (state.bmPollTimer) {
      clearInterval(state.bmPollTimer);
      state.bmPollTimer = null;
    }
    state.bmGenerating = false;
  }

  async function refreshBenchmark() {
    var btn = bm$("refresh");
    if (btn) btn.disabled = true;
    stopBmPolling();
    try {
      var data = await api(API.benchmarkRefresh, { method: "POST" });
      if (data && data.ok) {
        state.bmGenerating = true;
        state.bmGenDone = 0;
        state.bmGenTotal = data.total || 0;
        // Reload to pick up any partial cache
        var preview = await api(API.benchmarkPreview);
        state.bmRows = (preview && preview.rows) || [];
        renderBenchmark();
        startBmPolling();
      } else {
        var msg = (data && data.error) || "未知错误";
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

  /* ---------- Pipeline Debug ---------- */

  var debugDebounceTimer = null;
  var DEBUG_DEBOUNCE_MS = 400;

  function scheduleLiveDebug() {
    if (debugDebounceTimer) clearTimeout(debugDebounceTimer);
    debugDebounceTimer = setTimeout(function () {
      debugDebounceTimer = null;
      runDebug();
    }, DEBUG_DEBOUNCE_MS);
  }

  async function runDebug() {
    var input = document.getElementById("debug-smiles");
    if (!input) return;
    var smiles = input.value.trim();
    if (!smiles) return;

    var btn = document.getElementById("debug-run");
    if (btn) btn.disabled = true;

    try {
      var res = await api(API.debug, { method: "POST", body: JSON.stringify({ smiles: smiles }) });
      renderDebug(res);
    } catch (err) {
      var out = document.getElementById("debug-output");
      if (out) out.innerHTML = '<div class="debug-empty"><h2>Error</h2><p>' + escapeHtml(String(err)) + '</p></div>';
    } finally {
      if (btn) btn.disabled = false;
    }
  }

  function renderDebug(data) {
    var out = document.getElementById("debug-output");
    if (!out) return;
    if (!data.success) {
      out.innerHTML = '<div class="debug-empty"><h2>Error</h2><p>' + escapeHtml(data.error || "Unknown error") + '</p></div>';
      return;
    }

    var layers = data.layers;
    var layerDefs = [
      { key: "L0_preprocess",   icon: "l0", name: "Layer 0 — Preprocess",     desc: "SMILES → Mol (RDKit sanitization, kekulization)" },
      { key: "L1_analyze",      icon: "l1", name: "Layer 1 — Analyze",         desc: "Mol → Info dict (FGs, ring systems, unsaturation)" },
      { key: "L2_parent",       icon: "l2", name: "Layer 2 — Parent Select",   desc: "Info → Parent candidates → Selected parent hydride" },
      { key: "L3_substituents", icon: "l3", name: "Layer 3 — Substituents",    desc: "Parent + Info → Substituent extraction + Coverage ledger" },
      { key: "L4_numbering",    icon: "l4", name: "Layer 4 — Numbering",       desc: "Parent + Substituents → Numbered dict (locants, chain orient)" },
      { key: "L5_assemble",     icon: "l5", name: "Layer 5 — Assembly",        desc: "Numbered dict → NameResult (en + zh + stereo)" },
    ];

    var html = '<div class="summary-row">';
    html += '<span class="chip info">Total: <b>' + data.total_time_ms + ' ms</b></span>';
    for (var i = 0; i < layerDefs.length; i++) {
      var ld = layerDefs[i];
      var l = layers[ld.key];
      if (l) {
        html += '<span class="chip">' + ld.key.replace("_"," ").replace("_"," ") + ': <b>' + (l.time_ms || "?") + ' ms</b></span>';
      }
    }
    html += '</div>';

    html += '<div class="layer-list">';
    for (var j = 0; j < layerDefs.length; j++) {
      var ld2 = layerDefs[j];
      var l2 = layers[ld2.key];
      if (!l2) continue;
      var hasError = !!l2.error;
      html += '<div class="layer-card' + (hasError ? ' error' : '') + '" data-layer="' + ld2.key + '">';
      html += '<div class="layer-head" onclick="ChemNamerDebug.toggleLayer(this)">';
      html += '<div class="layer-icon ' + ld2.icon + '">' + ld2.key[1] + '</div>';
      html += '<div class="layer-info"><div class="layer-name">' + ld2.name + '</div>';
      html += '<div class="layer-desc">' + ld2.desc + '</div></div>';
      html += '<div class="layer-meta">';
      if (hasError) {
        html += '<span class="badge fail">ERROR</span>';
      }
      html += '<span class="layer-time">' + (l2.time_ms || "?") + ' ms</span>';
      html += '<svg class="layer-chevron" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><polyline points="6 9 12 15 18 9"/></svg>';
      html += '</div></div>';
      html += '<div class="layer-content">' + renderDebugLayer(ld2.key, l2) + '</div>';
      html += '</div>';
    }
    html += '</div>';

    out.innerHTML = html;
  }

  function renderDebugLayer(key, l) {
    if (l.error) {
      return '<div style="color:var(--color-destructive-soft);font-family:var(--mono);">' + escapeHtml(l.error) + '</div>';
    }
    switch (key) {
      case "L0_preprocess":   return renderDebugL0(l);
      case "L1_analyze":      return renderDebugL1(l);
      case "L2_parent":       return renderDebugL2(l);
      case "L3_substituents": return renderDebugL3(l);
      case "L4_numbering":    return renderDebugL4(l);
      case "L5_assemble":     return renderDebugL5(l);
      default:                return '<pre class="json-block">' + JSON.stringify(l, null, 2) + '</pre>';
    }
  }

  function renderDebugL0(l) {
    return '<dl class="debug-kv">' +
      '<dt>SMILES</dt><dd>' + escapeHtml(l.smiles || "") + '</dd>' +
      '<dt>Atoms</dt><dd>' + l.num_atoms + '</dd>' +
      '<dt>Bonds</dt><dd>' + l.num_bonds + '</dd>' +
      '<dt>Composition</dt><dd>' + escapeHtml(l.composition || "") + '</dd>' +
      '</dl>';
  }

  function renderDebugL1(l) {
    var html = '<dl class="debug-kv">';
    html += '<dt>n_carbons</dt><dd>' + l.n_carbons + '</dd>';
    html += '<dt>n_ring_systems</dt><dd>' + l.n_ring_systems + '</dd>';
    html += '<dt>n_rings</dt><dd>' + l.n_rings + '</dd>';
    html += '<dt>molecule</dt><dd>' + (l.mol ? (l.mol.num_atoms + ' atoms, ' + l.mol.num_bonds + ' bonds') : '—') + '</dd>';
    html += '</dl>';

    // Ring systems
    if (l.ring_systems && l.ring_systems.length) {
      html += '<div class="debug-sub"><div class="debug-sub-title">Ring Systems (' + l.ring_systems.length + ')</div>';
      for (var i = 0; i < l.ring_systems.length; i++) {
        var rs = l.ring_systems[i];
        var hets = (rs.hetero_atoms || []).map(function (h) { return 'idx=' + h.idx + '(Z=' + h.Z + ')'; }).join(', ');
        html += '<div class="candidate-card">';
        html += '<span class="badge">' + rs.topology + '</span> ';
        html += '<span class="chip">' + rs.n_rings + ' rings</span> ';
        html += '<span class="chip">' + rs.n_atoms + ' atoms</span> ';
        html += '<span class="chip ' + (rs.is_aromatic_mancude ? 'accent' : '') + '">mancude=' + rs.is_aromatic_mancude + '</span>';
        if (hets) html += '<span class="chip info">hetero: ' + escapeHtml(hets) + '</span>';
        html += '<div style="margin-top:4px"><span class="atom-list">atoms=[' + (rs.atom_ids || []).slice(0, 30).join(',') + (rs.atom_ids && rs.atom_ids.length > 30 ? '...' : '') + ']</span></div>';
        html += '</div>';
      }
      html += '</div>';
    }

    // Rings
    if (l.rings && l.rings.length) {
      html += '<div class="debug-sub"><div class="debug-sub-title">SSSR Rings (' + l.rings.length + ')</div>';
      html += '<table class="subst-table"><thead><tr><th>#</th><th>Size</th><th>Atoms</th></tr></thead><tbody>';
      for (var j = 0; j < l.rings.length; j++) {
        var r = l.rings[j];
        var ids = r.atom_ids || [];
        html += '<tr><td>' + j + '</td><td>' + ids.length + '-membered</td><td class="atom-list">[' + ids.join(',') + ']</td></tr>';
      }
      html += '</tbody></table></div>';
    }

    // FG flags
    html += '<div class="debug-sub"><div class="debug-sub-title">Functional Groups</div>';
    var fgMap = [
      ['has_acid','COOH'],['has_alcohol','OH'],['has_alkene','C=C'],['has_alkyne','C≡C'],
      ['has_amide','CON'],['has_amine','NH2/NH'],['has_anhydride','(CO)2O'],['has_boronic','B(OH)2'],
      ['has_carbamate','OCON'],['has_carbonate','OCOO'],['has_ester','COOR'],['has_ether','C-O-C'],
      ['has_guanidine','N-C(=N)N'],['has_hydrazine','N-N'],['has_isocyanate','NCO'],
      ['has_isothiocyanate','NCS'],['has_ketone','C=O'],['has_nitrile','C≡N'],['has_nitro','NO2'],
      ['has_phosphate','OPO3'],['has_sulfide','C-S-C'],['has_sulfonamide','SO2N'],
      ['has_sulfonate','SO3R'],['has_sulfone','SO2'],['has_sulfonic_acid','SO3H'],
      ['has_sulfonyl_chloride','SO2Cl'],['has_sulfoxide','SO'],['has_thiol','SH'],['has_urea','NCON'],
    ];
    html += '<div class="summary-row">';
    for (var k = 0; k < fgMap.length; k++) {
      var f = fgMap[k];
      if (l[f[0]] === true) html += '<span class="chip accent">' + f[1] + '</span>';
    }
    html += '</div>';
    html += '<span class="small muted-text">hydroxyls=' + (l.hydroxyls||[]).length +
      ' amines=' + (l.amines||[]).length +
      ' carboxyls=' + (l.carboxyls||[]).length +
      ' esters=' + (l.esters||[]).length +
      ' ethers=' + (l.ethers||[]).length +
      ' ketones=' + (l.ketones||[]).length +
      ' amides=' + (l.amides||[]).length +
      ' dbl_bonds=' + (l.double_bonds||[]).length + '</span>';
    html += '</div>';

    return html;
  }

  function renderDebugL2(l) {
    var html = '<dl class="debug-kv">';
    html += '<dt>n_candidates</dt><dd>' + l.n_candidates + '</dd>';
    html += '<dt>selected kind</dt><dd style="color:var(--color-accent-soft);font-weight:600;">' + escapeHtml(l.selected && l.selected.kind || '—') + '</dd>';
    html += '</dl>';

    if (l.candidates && l.candidates.length) {
      html += '<div class="debug-sub"><div class="debug-sub-title">All Candidates</div>';
      for (var i = 0; i < l.candidates.length; i++) {
        var c = l.candidates[i];
        var sel = c.is_selected;
        html += '<div class="candidate-card' + (sel ? ' selected' : '') + '">';
        html += '<div class="candidate-kind">' + escapeHtml(c.kind || '?');
        if (sel) html += '<span class="sel-tag">SELECTED</span>';
        html += '</div>';
        html += '<pre class="json-block inline">' + JSON.stringify(c, null, 2) + '</pre>';
        html += '</div>';
      }
      html += '</div>';
    }

    return html;
  }

  function renderDebugL3(l) {
    var html = '<dl class="debug-kv">';
    html += '<dt>substituents</dt><dd>' + l.n_substituents + '</dd>';
    html += '<dt>coverage complete</dt><dd><span class="status-dot ' + (l.coverage && l.coverage.complete ? 'ok' : 'fail') + '"></span>' + (l.coverage && l.coverage.complete ? 'YES' : 'NO') + '</dd>';
    if (l.coverage && l.coverage.gap && l.coverage.gap.length) {
      html += '<dt>gap atoms</dt><dd class="atom-list">[' + l.coverage.gap.join(',') + ']</dd>';
    }
    if (l.coverage && l.coverage.overlap && l.coverage.overlap.length) {
      html += '<dt>overlap atoms</dt><dd class="atom-list">[' + l.coverage.overlap.join(',') + ']</dd>';
    }
    html += '<dt>owned atoms</dt><dd class="atom-list">[' + (l.owned_atoms||[]).slice(0,40).join(',') + ((l.owned_atoms||[]).length > 40 ? '...' : '') + ']</dd>';
    html += '</dl>';

    if (l.substituents && l.substituents.length) {
      html += '<div class="debug-sub"><div class="debug-sub-title">Substituent List</div>';
      html += '<table class="subst-table"><thead><tr><th>#</th><th>EN</th><th>ZH</th><th>Kind</th><th>Attach</th><th>Paren</th><th>Atoms</th></tr></thead><tbody>';
      for (var i = 0; i < l.substituents.length; i++) {
        var s = l.substituents[i];
        html += '<tr>' +
          '<td>' + i + '</td>' +
          '<td style="font-weight:500;">' + escapeHtml(s.en || '?') + '</td>' +
          '<td>' + escapeHtml(s.zh || '') + '</td>' +
          '<td>' + escapeHtml(s.kind || '') + '</td>' +
          '<td>' + (s.attach_idx !== undefined ? s.attach_idx : '—') + '</td>' +
          '<td>' + (s.paren ? 'yes' : '') + '</td>' +
          '<td class="atom-list">[' + (s.atoms||[]).join(',') + ']</td>' +
          '</tr>';
      }
      html += '</tbody></table></div>';
    }

    return html;
  }

  function renderDebugL4(l) {
    var html = '<dl class="debug-kv">';
    html += '<dt>parent kind</dt><dd style="color:var(--color-accent-soft);font-weight:600;">' + escapeHtml(l.parent_kind || '—') + '</dd>';
    html += '<dt>chain</dt><dd class="atom-list">[' + (l.chain||[]).join(',') + ']</dd>';
    html += '</dl>';

    if (l.locants && Object.keys(l.locants).length) {
      html += '<div class="debug-sub"><div class="debug-sub-title">Locants</div>';
      html += '<dl class="debug-kv">';
      var locKeys = Object.keys(l.locants);
      for (var i = 0; i < locKeys.length; i++) {
        var k = locKeys[i];
        html += '<dt>' + escapeHtml(k) + '</dt><dd>' + escapeHtml(String(l.locants[k])) + '</dd>';
      }
      html += '</dl></div>';
    }

    if (l.substituents && l.substituents.length) {
      html += '<div class="debug-sub"><div class="debug-sub-title">Numbered Substituents (' + l.substituents.length + ')</div>';
      html += '<table class="subst-table"><thead><tr><th>#</th><th>EN</th><th>Locant</th><th>Kind</th></tr></thead><tbody>';
      for (var j = 0; j < l.substituents.length; j++) {
        var s = l.substituents[j];
        html += '<tr>' +
          '<td>' + j + '</td>' +
          '<td style="font-weight:500;">' + escapeHtml(s.en || '?') + '</td>' +
          '<td>' + (s.locant !== undefined ? s.locant : '—') + '</td>' +
          '<td>' + escapeHtml(s.kind || '') + '</td>' +
          '</tr>';
      }
      html += '</tbody></table></div>';
    }

    return html;
  }

  function renderDebugL5(l) {
    var html = '';
    if (l.en || l.zh) {
      html += '<div style="display:flex;flex-direction:column;gap:var(--space-md);">';

      html += '<div style="background:var(--color-surface-2);border:1px solid rgba(34,197,94,0.3);border-radius:var(--radius);padding:var(--space-lg);">';
      html += '<div class="field-label">English Name</div>';
      html += '<div style="font-family:var(--mono);font-size:15px;font-weight:500;color:var(--color-foreground);word-break:break-all;">' + escapeHtml(l.en) + '</div>';
      html += '</div>';

      html += '<div style="background:var(--color-surface-2);border:1px solid rgba(56,189,248,0.3);border-radius:var(--radius);padding:var(--space-lg);">';
      html += '<div class="field-label">Chinese Name</div>';
      html += '<div style="font-family:var(--mono);font-size:15px;font-weight:500;color:var(--color-foreground);word-break:break-all;">' + escapeHtml(l.zh) + '</div>';
      html += '</div>';

      html += '<dl class="debug-kv">';
      html += '<dt>success</dt><dd>' + l.success + '</dd>';
      html += '<dt>source</dt><dd>' + escapeHtml(l.source || '') + '</dd>';
      html += '</dl>';

      html += '</div>';
    }
    if (l.meta && Object.keys(l.meta).length) {
      html += '<div class="debug-sub"><div class="debug-sub-title">Meta</div>';
      html += '<pre class="json-block inline">' + JSON.stringify(l.meta, null, 2) + '</pre></div>';
    }
    return html || '<span class="muted-text">(no output)</span>';
  }

  /* ---------- Call Graph ---------- */

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

  function cgChainById() {
    var idx = {};
    ((state.cgData && state.cgData.nodes) || []).forEach(function (n) { idx[n.id] = n; });
    return idx;
  }

  function cgChainShort(nd) {
    var mod = (nd.module || "").split("/").pop().replace(".py", "");
    return nd.label + " @" + mod + ":" + nd.line;
  }

  function renderCalleeChain(rootId) {
    var data = state.cgData;
    if (!data || !data.nodes) return;
    var byId = cgChainById();
    var root = byId[rootId];
    if (!root) return;
    var outEdges = {};
    data.edges.forEach(function (e) {
      (outEdges[e.from] = outEdges[e.from] || []).push(e);
    });
    function fmt(nd) {
      return cgChainShort(nd) + "  " + nd.cum_pct + "% (" + nd.ncalls.toLocaleString() + "×)";
    }
    function walk(nid, depth) {
      var callees = (outEdges[nid] || []).slice().sort(function (a, b) {
        return byId[b.to].cum_s - byId[a.to].cum_s;
      });
      var html = "";
      callees.forEach(function (e) {
        var nd = byId[e.to];
        var pad = depth * 16 + 8;
        html +=
          '<div class="cg-chain-row cg-chain-name" style="padding-left:' + pad + 'px" data-id="' + nd.id + '">' +
          '<span class="mono">' + escapeHtml(fmt(nd)) + "</span>" +
          (e.cum_pct ? '<span class="muted-text small">  ←' + e.calls + "×</span>" : "") +
          "</div>";
        if (seen[nd.id]) {
          html += '<div class="cg-chain-row cg-chain-loop" style="padding-left:' + (pad + 16) + 'px">↺ 循环</div>';
        } else {
          seen[nd.id] = true;
          html += walk(nd.id, depth + 1);
        }
      });
      return html;
    }
    var seen = {};
    seen[rootId] = true;
    var el = $("cg-chain");
    var head = '<div class="cg-chain-row cg-chain-root"><b>' + escapeHtml(fmt(root)) + "</b></div>";
    el.innerHTML = head + walk(rootId, 0);
    el.querySelectorAll(".cg-chain-name").forEach(function (row) {
      row.addEventListener("click", function () {
        renderCalleeChain(parseInt(row.getAttribute("data-id"), 10));
      });
    });
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
    ["cg-rank-calls", "cg-rank-time"].forEach(function (id) {
      var list = $(id);
      if (list) {
        list.addEventListener("click", function (ev) {
          var row = ev.target.closest(".cg-rank-clickable");
          if (row) renderCalleeChain(parseInt(row.getAttribute("data-id"), 10));
        });
      }
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

  async function loadCallGraph(force) {
    if (state.cgLoading) return;
    var n = cgNFromSlider();
    var st = $("cg-status");
    var btn = $("cg-regenerate");
    var empty = $("cg-empty");
    state.cgLoading = true;
    if (btn) btn.disabled = true;
    if (st) { st.textContent = "采样中…（约 1–3s）"; st.style.color = "#f59e0b"; }
    if (empty) empty.hidden = true;
    try {
      var q = "?n=" + n + "&floor_pct=0" + (force ? "&refresh=1" : "");
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
      if (st) { st.textContent = "就绪"; st.style.color = "#86efac"; }
    } catch (err) {
      if (empty) { empty.hidden = false; empty.textContent = "加载失败: " + (err.message || String(err)); }
      if (st) { st.textContent = "错误"; st.style.color = "#ef4444"; }
    } finally {
      state.cgLoading = false;
      if (btn) btn.disabled = false;
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

  async function loadCallGraphSvg(force) {
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
      }
      if (st) { st.textContent = "就绪 · " + data.meta.n_nodes + " 节点"; st.style.color = "#86efac"; }
    } catch (err) {
      if (st) { st.textContent = "失败: " + (err.message || String(err)); st.style.color = "#ef4444"; }
    } finally {
      state.cgSvgLoading = false;
      if (btn) btn.disabled = false;
    }
  }

  /* ---------- Init ---------- */

  function bind() {
    $("namer-form") && $("namer-form").addEventListener("submit", onName);
    $("btn-clear-history") &&
      $("btn-clear-history").addEventListener("click", () => {
        state.namerHistory = [];
        renderNamerHistory();
      });
    $("btn-load-smiles") &&
      $("btn-load-smiles").addEventListener("click", async () => {
        const smiles = (($("smiles-input") && $("smiles-input").value) || "").trim();
        if (!smiles) {
          setNamerError("请输入 SMILES 再载入画板");
          return;
        }
        if (!state.ketcherBridge || !state.ketcherBridge.isReady()) {
          setNamerError("Ketcher 未就绪");
          return;
        }
        try {
          setNamerError("");
          await state.ketcherBridge.setMolecule(smiles);
        } catch (err) {
          setNamerError(err.message || "载入结构失败");
        }
      });
    $("btn-clear-ketcher") &&
      $("btn-clear-ketcher").addEventListener("click", async () => {
        if (!state.ketcherBridge || !state.ketcherBridge.isReady()) return;
        try {
          await state.ketcherBridge.clear();
          // Also clear the SMILES input
          const input = $("smiles-input");
          if (input) input.value = "";
          setNamerError("");
        } catch (err) {
          setNamerError(err.message || "清空失败");
        }
      });
    $("btn-retry-ketcher") &&
      $("btn-retry-ketcher").addEventListener("click", () => {
        ensureKetcher(true);
      });
    $("live-name-toggle") &&
      $("live-name-toggle").addEventListener("change", (ev) => {
        state.liveNameEnabled = !!(ev.target && ev.target.checked);
        if (!state.liveNameEnabled) {
          if (state.liveDebounceTimer) {
            clearTimeout(state.liveDebounceTimer);
            state.liveDebounceTimer = null;
          }
          if (state.liveSmilesTimer) {
            clearTimeout(state.liveSmilesTimer);
            state.liveSmilesTimer = null;
          }
        }
      });

    // Live SMILES input: when user edits the SMILES text box, auto-name + load into Ketcher
    $("smiles-input") &&
      $("smiles-input").addEventListener("input", () => {
        scheduleLiveSmilesName();
      });

    // Tab switching
    document.querySelectorAll(".tab").forEach(function (t) {
      t.addEventListener("click", function () {
        var page = t.getAttribute("data-page");
        if (page) switchPage(page);
      });
    });

    // Benchmark controls
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

    // Benchmark Run controls
    var brRun = document.getElementById("br-run");
    if (brRun) {
      brRun.addEventListener("click", startBenchmarkRun);
    }

    // Debug page controls
    var debugRun = document.getElementById("debug-run");
    if (debugRun) debugRun.addEventListener("click", runDebug);
    var debugSmiles = document.getElementById("debug-smiles");
    if (debugSmiles) {
      debugSmiles.addEventListener("input", scheduleLiveDebug);
    }

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

  /* ---------- Benchmark Run ---------- */

  function br$(id) {
    return document.getElementById("br-" + id);
  }

  function stopBrPolling() {
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

  /* ---------- Layer Benchmark ---------- */

  function lb$(id) {
    return document.getElementById("lb-" + id);
  }

  function stopLbPolling() {
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

  function renderLbPanel(layer) {
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

  /* ---------- Code Analysis ---------- */

  var CA_COLORS = ["#22c55e", "#38bdf8", "#f59e0b", "#a78bfa", "#f472b6", "#2dd4bf"];

  // 深浅变体：amt>0 向白混合（变浅），amt<0 向黑混合（变深）
  function caShade(hex, amt) {
    var c = parseInt(hex.slice(1), 16);
    var r = (c >> 16) & 255, g = (c >> 8) & 255, b = c & 255;
    if (amt >= 0) {
      r = Math.round(r + (255 - r) * amt);
      g = Math.round(g + (255 - g) * amt);
      b = Math.round(b + (255 - b) * amt);
    } else {
      amt = -amt;
      r = Math.round(r * (1 - amt));
      g = Math.round(g * (1 - amt));
      b = Math.round(b * (1 - amt));
    }
    return "#" + ((1 << 24) + (r << 16) + (g << 8) + b).toString(16).slice(1);
  }

  function ca$(id) {
    return document.getElementById(id);
  }

  function caSkeleton() {
    return (
      '<div class="ca-skeleton" aria-hidden="true">' +
      '<div class="skeleton" style="width:180px;height:180px;border-radius:50%;margin:0 auto"></div>' +
      "</div>" +
      '<div class="ca-list">' +
      Array.from({ length: 6 })
        .map(function () {
          return (
            '<div class="ca-row"><div class="skeleton" style="width:12px;height:12px;border-radius:50%"></div>' +
            '<div class="skeleton" style="width:80px;height:14px"></div>' +
            '<div class="skeleton" style="width:60px;height:14px"></div>' +
            '<div class="skeleton" style="flex:1;height:8px;border-radius:4px"></div></div>'
          );
        })
        .join("") +
      "</div>"
    );
  }

  async function loadCodeAnalysis() {
    var donut = ca$("ca-donut");
    var list = ca$("ca-list");
    var summary = ca$("ca-summary");
    if (donut) donut.innerHTML = caSkeleton();
    if (list) list.innerHTML = "";
    if (summary) { summary.classList.remove("hidden"); summary.textContent = "统计中…"; }
    try {
      var data = await api(API.codeAnalysis);
      if (!data || !data.ok || !data.total) throw new Error("无数据");
      state.caData = data;
      if (summary) renderCaSummary(summary, data.total);
      if (donut) donut.innerHTML = renderCaDonut(data.layers, data.total.code);
      if (list) list.innerHTML = renderCaList(data.layers);
    } catch (err) {
      if (summary) {
        summary.classList.remove("hidden");
        summary.textContent = "加载失败: " + (err.message || String(err));
      }
      if (donut) donut.innerHTML = '<p class="muted-text">暂无统计数据</p>';
    }
  }

  function renderCaSummary(el, total) {
    el.innerHTML =
      '<div class="ca-stat"><span class="ca-stat-label">总文件</span><b class="mono">' + (total.file_count != null ? total.file_count : total.files) + "</b></div>" +
      '<div class="ca-stat"><span class="ca-stat-label">代码行</span><b class="mono">' + total.code.toLocaleString() + "</b></div>" +
      '<div class="ca-stat"><span class="ca-stat-label">注释行</span><b class="mono">' + total.comment.toLocaleString() + "</b></div>" +
      '<div class="ca-stat"><span class="ca-stat-label">总行数</span><b class="mono">' + total.lines.toLocaleString() + "</b></div>";
  }

  function renderCaDonut(layers, totalCode) {
    var size = 220, cx = 110, cy = 110, r = 78, w = 26;
    var circ = 2 * Math.PI * r;
    // 最大扇区从 12 点起，顺时针按占比降序
    var sorted = layers.slice().sort(function (a, b) { return b.code - a.code; });
    var acc = 0;
    var segs = sorted.map(function (l) {
      var frac = totalCode > 0 ? l.code / totalCode : 0;
      var len = frac * circ;
      var seg =
        '<circle cx="' + cx + '" cy="' + cy + '" r="' + r + '" fill="none" ' +
        'stroke="' + CA_COLORS[l.layer] + '" stroke-width="' + w + '" ' +
        'stroke-dasharray="' + len.toFixed(2) + ' ' + (circ - len).toFixed(2) + '" ' +
        'stroke-dashoffset="' + (-acc).toFixed(2) + '" ' +
        'transform="rotate(-90 ' + cx + ' ' + cy + ')" class="ca-seg">' +
        "<title>Layer " + l.layer + ": " + l.code + " 行 (" + l.share + "%)</title></circle>";
      acc += len;
      return seg;
    });
    return (
      '<svg viewBox="0 0 ' + size + " " + size + '" class="ca-donut" role="img" aria-label="各 Layer 代码行占比">' +
      segs.join("") +
      '<text x="' + cx + '" y="' + (cy - 4) + '" class="ca-donut-total" text-anchor="middle">' + totalCode.toLocaleString() + "</text>" +
      '<text x="' + cx + '" y="' + (cy + 18) + '" class="ca-donut-cap" text-anchor="middle">code lines</text>' +
      "</svg>"
    );
  }

  function renderCaList(layers) {
    var maxCode = layers.reduce(function (m, l) { return Math.max(m, l.code); }, 0);
    return (
      '<div class="ca-list-head">' +
      '<span></span><span>Layer</span><span class="mono">文件</span><span class="mono">代码行</span><span>占比</span><span class="mono">%</span><span></span>' +
      "</div>" +
      layers
        .map(function (l) {
          var bw = maxCode > 0 ? Math.round(100 * l.code / maxCode) : 0;
          return (
            '<div class="ca-row ca-row-toggle" data-layer="' + l.layer + '" role="button" tabindex="0" aria-expanded="false">' +
            '<span class="ca-dot" style="background:' + CA_COLORS[l.layer] + '"></span>' +
            '<span class="ca-name">Layer ' + l.layer + "</span>" +
            '<span class="ca-files mono">' + l.file_count + "</span>" +
            '<span class="ca-lines mono">' + l.code.toLocaleString() + "</span>" +
            '<span class="ca-bar-wrap"><span class="ca-bar" style="width:' + bw + "%;background:" + CA_COLORS[l.layer] + '"></span></span>' +
            '<span class="ca-share mono">' + l.share + "%</span>" +
            '<span class="ca-chevron" aria-hidden="true"></span>' +
            "</div>" +
            '<div class="ca-file-sub hidden" data-file-layer="' + l.layer + '"></div>'
          );
        })
        .join("")
    );
  }

  function caFileRows(layer) {
    var l = (state.caData && state.caData.layers || []).find(function (x) { return x.layer === layer; });
    if (!l || !l.files || !l.files.length) return '<p class="muted-text small">无文件</p>';
    var color = CA_COLORS[layer];
    var shade = caShade(color, 0.6); // 浅色 = 总行数（含注释/空行）
    var maxTotal = l.files.reduce(function (m, f) { return Math.max(m, f.code + f.comment + f.blank); }, 0);
    return (
      '<div class="ca-file-row ca-file-head">' +
      '<span class="ca-file-name">文件名</span><span class="mono">代码行</span><span>分布</span><span class="mono">总行数</span>' +
      "</div>" +
      l.files
        .map(function (f) {
          var total = f.code + f.comment + f.blank;
          var codeW = maxTotal > 0 ? Math.round(100 * f.code / maxTotal) : 0;
          var otherW = maxTotal > 0 ? Math.round(100 * (total - f.code) / maxTotal) : 0;
          return (
            '<div class="ca-file-row">' +
            '<span class="ca-file-name mono">' + escapeHtml(f.name) + "</span>" +
            '<span class="ca-file-lines mono" style="color:' + color + '">' + f.code.toLocaleString() + "</span>" +
            '<span class="ca-file-bar-wrap">' +
            '<span class="ca-file-bar" style="width:' + codeW + "%;background:" + color + '"></span>' +
            '<span class="ca-file-bar ca-file-bar-shade" style="width:' + otherW + "%;background:" + shade + '"></span>' +
            "</span>" +
            '<span class="ca-file-total mono" style="color:' + shade + '">' + total.toLocaleString() + "</span>" +
            "</div>"
          );
        })
        .join("")
    );
  }

  function toggleCaRow(row) {
    var layer = parseInt(row.getAttribute("data-layer"), 10);
    var expanded = row.getAttribute("aria-expanded") === "true";
    row.setAttribute("aria-expanded", String(!expanded));
    row.classList.toggle("open", !expanded);
    var sub = document.querySelector('.ca-file-sub[data-file-layer="' + layer + '"]');
    if (sub) {
      if (!expanded && !sub.getAttribute("data-rendered")) {
        sub.innerHTML = caFileRows(layer);
        sub.setAttribute("data-rendered", "1");
      }
      sub.classList.toggle("hidden", expanded);
    }
  }

  function bindCaToggle() {
    var list = ca$("ca-list");
    if (!list) return;
    list.addEventListener("click", function (ev) {
      var row = ev.target.closest(".ca-row-toggle");
      if (row) toggleCaRow(row);
    });
    list.addEventListener("keydown", function (ev) {
      if (ev.key !== "Enter" && ev.key !== " ") return;
      var row = ev.target.closest(".ca-row-toggle");
      if (row) {
        ev.preventDefault();
        toggleCaRow(row);
      }
    });
  }

  function bindCaRefresh() {
    var btn = $("ca-refresh");
    if (!btn) return;
    btn.addEventListener("click", async function () {
      var st = $("ca-refresh-status");
      if (st) st.textContent = "刷新中…";
      await loadCodeAnalysis();
      if (st) st.textContent = "已刷新 " + new Date().toLocaleTimeString();
    });
  }

  /* ---------- Settings ---------- */

  function settings$(id) {
    return document.getElementById("settings-" + id);
  }

  function applyTheme(theme) {
    var t = theme === "light" ? "light" : "dark";
    document.documentElement.setAttribute("data-theme", t);
    try { localStorage.setItem("chem-theme", t); } catch (_) {}
  }

  function initTheme() {
    var t = null;
    try { t = localStorage.getItem("chem-theme"); } catch (_) {}
    if (t !== "light" && t !== "dark") t = "dark";
    applyTheme(t);
    var sel = settings$("theme");
    if (sel) sel.value = t;
  }

  function fillApiKeyHint(s) {
    var hint = settings$("api-key-hint");
    if (!hint) return;
    hint.textContent = s && s.api_key_configured
      ? "已配置（……" + (s.api_key_hint || "") + "），留空保持不变"
      : "未配置 API Key";
  }

  async function loadSettings() {
    try {
      var data = await api(API.settings);
      if (!data || !data.ok) throw new Error((data && data.error) || "加载失败");
      var s = data.settings || {};
      if (settings$("model")) settings$("model").value = s.model || "";
      if (settings$("base-url")) settings$("base-url").value = s.base_url || "";
      if (settings$("concurrency")) settings$("concurrency").value = s.concurrency || 1;
      if (settings$("theme")) settings$("theme").value = s.theme || "dark";
      fillApiKeyHint(s);
      applyTheme(s.theme || "dark");
    } catch (err) {
      var st = settings$("status");
      if (st) st.textContent = "加载失败: " + (err.message || err);
    }
  }

  async function saveSettings() {
    var body = {
      model: settings$("model") ? settings$("model").value.trim() : "",
      base_url: settings$("base-url") ? settings$("base-url").value.trim() : "",
      concurrency: parseInt((settings$("concurrency") && settings$("concurrency").value) || "1", 10) || 1,
      theme: settings$("theme") ? settings$("theme").value : "dark",
    };
    var ak = settings$("api-key");
    if (ak && ak.value && ak.value.trim()) body.api_key = ak.value.trim();

    var st = settings$("status");
    if (st) { st.textContent = "保存中…"; st.className = "chip"; }
    try {
      var data = await api(API.settings, { method: "POST", body: JSON.stringify(body) });
      if (!data || !data.ok) throw new Error((data && data.error) || "保存失败");
      if (ak) ak.value = "";
      var s = data.settings || {};
      fillApiKeyHint(s);
      applyTheme(s.theme || "dark");
      if (st) { st.textContent = "已保存 " + new Date().toLocaleTimeString(); st.className = "chip accent"; }
    } catch (err) {
      if (st) { st.textContent = "保存失败: " + (err.message || err); st.className = "chip warn"; }
    }
  }

  function bindSettingsControls() {
    var saveBtn = document.getElementById("settings-save");
    if (saveBtn) saveBtn.addEventListener("click", saveSettings);
    var themeSel = settings$("theme");
    if (themeSel) {
      themeSel.addEventListener("change", function () {
        applyTheme(themeSel.value);
      });
    }
  }

  /* ---------- Wiki ---------- */

  function renderWikiTree(entry, depth) {
    var html = "";
    if (entry.name !== "wiki") {
      html +=
        '<details class="wiki-dir" data-path="' + escapeHtml(entry.path) + '"' + (depth === 0 ? " open" : "") + ">" +
        '<summary style="padding-left:' + (depth * 14 + 4) + 'px">' + escapeHtml(entry.name) + "</summary>";
    }
    entry.files.forEach(function (f) {
      html += '<div class="wiki-file" data-path="' + escapeHtml(f.path) + '" style="padding-left:' + ((depth + 1) * 14 + 8) + 'px">' + escapeHtml(f.name) + "</div>";
    });
    entry.dirs.forEach(function (d) {
      html += renderWikiTree(d, depth + 1);
    });
    if (entry.name !== "wiki") html += "</details>";
    return html;
  }

  async function loadWiki() {
    try {
      var data = await api(API.wiki);
      if (!data || !data.ok) throw new Error((data && data.error) || "加载失败");
      var treeEl = $("wiki-tree");
      if (treeEl) treeEl.innerHTML = renderWikiTree(data.tree, 0);
      loadWikiDoc("index.md");
    } catch (err) {
      var doc = $("wiki-doc");
      if (doc) doc.innerHTML = '<p class="muted-text">加载失败: ' + escapeHtml(err.message || err) + "</p>";
    }
  }

  async function loadWikiDoc(path) {
    // 解析相对路径：含 / 按 wiki 根相对；纯文件名按「当前文档同目录」优先、根兜底
    var cur = state.wikiCurrent || "";
    var curDir = cur.indexOf("/") !== -1 ? cur.substring(0, cur.lastIndexOf("/")) : "";
    var candidates = path.indexOf("/") !== -1
      ? [path]
      : (curDir ? [curDir + "/" + path, path] : [path]);

    var doc = $("wiki-doc");
    for (var i = 0; i < candidates.length; i++) {
      var data = null;
      try {
        data = await api(API.wikiDoc + "?path=" + encodeURIComponent(candidates[i]));
      } catch (e) {
        data = null;
      }
      if (data && data.ok) {
        state.wikiCurrent = data.path;
        var body = typeof marked !== "undefined" && marked.parse
          ? marked.parse(data.content)
          : escapeHtml(data.content);
        if (doc) doc.innerHTML = '<div class="wiki-title">' + escapeHtml(data.name) + "</div>" + body;
        document.querySelectorAll(".wiki-file.active").forEach(function (el) {
          el.classList.remove("active");
        });
        var sel = document.querySelector('.wiki-file[data-path="' + CSS.escape(data.path) + '"]');
        if (sel) sel.classList.add("active");
        return;
      }
    }
    if (doc) doc.innerHTML = '<p class="muted-text">无法打开文档: ' + escapeHtml(path) + "</p>";
  }

  function bindWikiControls() {
    var tree = $("wiki-tree");
    if (tree) {
      tree.addEventListener("click", function (ev) {
        var f = ev.target.closest(".wiki-file");
        if (f) loadWikiDoc(f.getAttribute("data-path"));
      });
    }
    var doc = $("wiki-doc");
    if (doc) {
      doc.addEventListener("click", function (ev) {
        var a = ev.target.closest("a");
        if (!a) return;
        var href = (a.getAttribute("href") || "").trim();
        // 外部链接 / 锚点保留默认行为
        if (!href || /^[a-z]+:/i.test(href) || href.charAt(0) === "#") return;
        // 仅拦截指向 wiki 内 .md 的相对链接
        if (!/\.md(?:[#?]|$)/i.test(href)) return;
        ev.preventDefault();
        loadWikiDoc(href.split(/[#?]/)[0]);
      });
    }
  }

  async function init() {
    bind();
    initTheme();
    bindCaToggle();
    bindCaRefresh();
    bindCgControls();
    bindCgRankButtons();
    bindCgChainControls();
    bindSettingsControls();
    bindWikiControls();
    renderNamerHistory();
    ensureKetcher();
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", init);
  } else {
    init();
  }

  // Expose debug toggleLayer for inline onclick handlers
  window.ChemNamerDebug = {
    toggleLayer: function (headEl) {
      headEl.parentElement.classList.toggle("open");
    }
  };
})();
