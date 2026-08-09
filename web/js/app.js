/* ChemAgent Namer — entry point: page switching + control wiring.
   Feature modules live in sibling files (namer/benchmark/...); this file only
   imports them and ties the global layout together. */
import { state } from "./core.js";
import { bindNamer, ensureKetcher } from "./namer.js";
import { loadBenchmark, stopBmPolling, bindBenchmark } from "./benchmark.js";
import { stopBrPolling, bindBenchmarkRun } from "./benchmark-run.js";
import { renderLbPanel, stopLbPolling, bindLayerBenchmark } from "./layer-benchmark.js";
import { bindDebug } from "./debug.js";
import { loadCodeAnalysis, bindCodeAnalysis } from "./code-analysis.js";
import { loadCallGraph, loadCallGraphSvg, bindCallGraph } from "./callgraph.js";
import { initTheme, loadSettings, bindSettings } from "./settings.js";
import { loadWiki, bindWiki } from "./wiki.js";
import { loadIupac, bindIupac } from "./iupac.js";
import { initHistoryPicker } from "./history.js";

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
  var iupacPage = document.getElementById("iupac-page");
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
  if (iupacPage) iupacPage.classList.add("hidden");
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
    // 采样完成后自动生成分层 SVG，避免进入页面后空白、需手动点「重新生成」
    loadCallGraph().then(function (ok) {
      if (ok && !state.cgSvg) loadCallGraphSvg();
    });
  } else if (name === "debug") {
    if (debugPage) debugPage.classList.remove("hidden");
  } else if (name === "settings") {
    if (settingsPage) settingsPage.classList.remove("hidden");
    loadSettings();
  } else if (name === "wiki") {
    if (wikiPage) wikiPage.classList.remove("hidden");
    loadWiki();
  } else if (name === "iupac") {
    if (iupacPage) iupacPage.classList.remove("hidden");
    loadIupac();
  } else {
    if (namerLayout) namerLayout.classList.remove("hidden");
  }
}

/* ---------- Init ---------- */

function bind() {
  // Tab switching
  document.querySelectorAll(".tab").forEach(function (t) {
    t.addEventListener("click", function () {
      var page = t.getAttribute("data-page");
      if (page) switchPage(page);
    });
  });

  // Per-page control wiring lives in each module.
  bindNamer();
  bindBenchmark();
  bindBenchmarkRun();
  bindLayerBenchmark();
  bindDebug();
  bindCodeAnalysis();
  bindCallGraph();
  bindSettings();
  bindWiki();
  bindIupac();
}

async function init() {
  bind();
  initTheme();
  ensureKetcher();
  initHistoryPicker();
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
