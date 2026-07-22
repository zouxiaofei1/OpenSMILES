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
    // Hide all first
    if (namerLayout) namerLayout.classList.add("hidden");
    if (bmPage) bmPage.classList.add("hidden");
    if (brPage) brPage.classList.add("hidden");
    // Stop polling
    stopBmPolling();
    stopBrPolling();
    // Show active page
    if (name === "benchmark") {
      if (bmPage) bmPage.classList.remove("hidden");
      if (!state.bmLoaded) {
        loadBenchmark();
      }
    } else if (name === "benchmark-run") {
      if (brPage) brPage.classList.remove("hidden");
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

  async function init() {
    bind();
    renderNamerHistory();
    ensureKetcher();
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", init);
  } else {
    init();
  }
})();
