/* History picker: git commit dropdown (datalist) + current-version badge.
   Selecting a commit points the three data pages (code-analysis / benchmark /
   call-graph) at that past version. The UI itself stays on the current code. */
import { state, api, API } from "./core.js";
import { loadBenchmark, stopBmPolling } from "./benchmark.js";
import { loadCodeAnalysis } from "./code-analysis.js";
import { loadCallGraph, loadCallGraphSvg } from "./callgraph.js";

let commits = [];
const byLabel = new Map(); // option.label -> full hash

function displayOf(c) {
  return (c.short || "") + " · " + (c.date || "").slice(0, 10) + " · " + (c.subject || "");
}

function resolveInputValue(v) {
  v = (v || "").trim();
  if (!v) return null; // empty -> HEAD
  if (/^HEAD(\s|$)/i.test(v)) return null;
  if (byLabel.has(v)) return byLabel.get(v);
  const hit = commits.find(function (c) {
    return v === c.short || (c.hash && c.hash.startsWith(v)) || (c.subject || "").startsWith(v);
  });
  return hit ? hit.hash : undefined; // undefined = unrecognised
}

function updateBadge(hash) {
  const b = document.getElementById("history-badge");
  if (!b) return;
  if (!hash) {
    b.textContent = "HEAD";
    b.classList.remove("info");
    return;
  }
  const c = commits.find(function (x) { return x.hash === hash; });
  b.textContent = "查看 " + (c ? c.short + " · " + (c.subject || "") : hash.slice(0, 7));
  b.classList.add("info");
}

function applyCommit(hash) {
  state.currentCommit = hash; // null = HEAD
  // Invalidate the three data pages so they refetch with ?commit=
  state.bmLoaded = false;
  state.bmRows = [];
  state.bmGenerating = false;
  state.caData = null;
  state.cgData = null;
  state.cgSvg = null;
  state.cgSourceTotal = 0;
  stopBmPolling();
  updateBadge(hash);

  var page = state.currentPage;
  if (page === "benchmark") {
    loadBenchmark();
  } else if (page === "code-analysis") {
    loadCodeAnalysis();
  } else if (page === "call-graph") {
    loadCallGraph().then(function (ok) {
      if (ok && !state.cgSvg) loadCallGraphSvg();
    });
  }
}

function onPick() {
  var input = document.getElementById("history-commit-input");
  var h = resolveInputValue(input.value);
  if (h === undefined) { input.value = ""; return; }
  applyCommit(h);
}

function onEnter(ev) {
  if (ev.key !== "Enter") return;
  var input = document.getElementById("history-commit-input");
  var h = resolveInputValue(input.value);
  if (h === undefined) return;
  applyCommit(h);
}

export async function initHistoryPicker() {
  var input = document.getElementById("history-commit-input");
  var list = document.getElementById("history-commit-list");
  if (!input || !list) return;
  try {
    var data = await api(API.gitCommits);
    commits = (data && data.commits) || [];
    list.innerHTML = "";
    var headOpt = document.createElement("option");
    headOpt.value = "HEAD · 当前工作区代码";
    list.appendChild(headOpt);
    for (var i = 0; i < commits.length; i++) {
      var c = commits[i];
      var o = document.createElement("option");
      o.value = displayOf(c);
      byLabel.set(o.value, c.hash);
      list.appendChild(o);
    }
    input.disabled = false;
  } catch (err) {
    input.placeholder = "提交列表加载失败";
    input.disabled = true;
  }
  input.addEventListener("change", onPick);
  input.addEventListener("keydown", onEnter);
}
