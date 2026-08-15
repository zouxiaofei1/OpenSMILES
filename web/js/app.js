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

/* ---------- Responsive nav drawer ---------- */

/* 抽屉从 .topbar-nav 克隆 tabs（单一事实来源），克隆带 .tab 类，
   因此下方的 switchPage 状态同步与点击绑定会自动覆盖抽屉项。 */
function buildNavDrawer() {
  var target = document.getElementById("nav-drawer-tabs");
  if (!target) return;
  target.innerHTML = "";
  document.querySelectorAll(".topbar-nav .tab").forEach(function (t) {
    var clone = t.cloneNode(true);
    clone.classList.add("drawer-tab");
    target.appendChild(clone);
  });
}

function openNav(open) {
  var btn = document.getElementById("nav-toggle");
  var drawer = document.getElementById("nav-drawer");
  if (!btn || !drawer) return;
  btn.setAttribute("aria-expanded", String(open));
  btn.setAttribute("aria-label", open ? "关闭菜单" : "打开菜单");
  btn.title = open ? "关闭菜单" : "打开菜单";
  drawer.setAttribute("aria-hidden", String(!open));
  drawer.classList.toggle("open", open);
  if (open) {
    // 锁定背景滚动；抽屉出现后把焦点移入首项
    document.body.style.overflow = "hidden";
    var first = drawer.querySelector(".tab");
    if (first) setTimeout(function () { first.focus(); }, 160);
  } else {
    document.body.style.overflow = "";
    btn.focus();
  }
}

function bindNavToggle() {
  var btn = document.getElementById("nav-toggle");
  var drawer = document.getElementById("nav-drawer");
  if (!btn || !drawer) return;

  btn.addEventListener("click", function (ev) {
    ev.stopPropagation();
    openNav(btn.getAttribute("aria-expanded") !== "true");
  });

  // 点击抽屉外关闭
  document.addEventListener("click", function (ev) {
    if (
      drawer.classList.contains("open") &&
      !drawer.contains(ev.target) &&
      ev.target !== btn &&
      !btn.contains(ev.target)
    ) {
      openNav(false);
    }
  });

  // Esc 关闭
  document.addEventListener("keydown", function (ev) {
    if (ev.key === "Escape" && drawer.classList.contains("open")) openNav(false);
  });

  // 点击抽屉内任意 tab 后关闭
  drawer.addEventListener("click", function (ev) {
    if (ev.target.closest(".tab")) openNav(false);
  });
}

/* 顶部栏滚动投影：内容区滚动超过阈值时给 topbar 加深阴影。 */
function bindScrollElevation() {
  var topbar = document.querySelector(".topbar");
  var main = document.querySelector(".main");
  if (!topbar) return;
  function onScroll() {
    var scrolled =
      (window.scrollY > 4) || (main && main.scrollTop > 4);
    topbar.classList.toggle("scrolled", scrolled);
  }
  document.addEventListener("scroll", onScroll, { passive: true });
  if (main) main.addEventListener("scroll", onScroll, { passive: true });
  onScroll();
}

/* 中屏导航条滚动溢出检测：仅在真正溢出时加边缘渐变遮罩。 */
function syncNavOverflow() {
  var nav = document.querySelector(".topbar-nav");
  if (!nav) return;
  var over = nav.scrollWidth > nav.clientWidth + 1;
  nav.classList.toggle("overflowing", over);
}

/* ---------- Init ---------- */

function bind() {
  buildNavDrawer(); // 先于 tab 点击绑定，保证抽屉项也能被绑定
  // Tab switching
  document.querySelectorAll(".tab").forEach(function (t) {
    t.addEventListener("click", function () {
      var page = t.getAttribute("data-page");
      if (page) switchPage(page);
    });
  });
  bindNavToggle();
  bindScrollElevation();
  syncNavOverflow();
  window.addEventListener("resize", function () {
    syncNavOverflow();
    // 视口变宽恢复内联导航时，自动收起抽屉避免两套导航共存
    var drawer = document.getElementById("nav-drawer");
    if (window.innerWidth > 768 && drawer && drawer.classList.contains("open")) {
      openNav(false);
    }
  });
  // 字体异步加载后 tab 实际宽度变化，需重算溢出遮罩
  if (document.fonts && document.fonts.ready) {
    document.fonts.ready.then(syncNavOverflow);
  }

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
