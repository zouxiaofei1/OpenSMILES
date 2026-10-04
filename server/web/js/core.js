/* OpenSMILES Namer — shared core: API endpoints, global state, DOM/HTTP helpers.
   All feature modules import from here; nothing in core imports back. */
export const API = {
  health: "/health",
  name: "/api/v1/name",
  resolveName: "/api/v1/name/resolve",
  pubchemIupac: "/api/v1/name/pubchem-iupac",
  locantsSvg: "/api/v1/name/locants-svg",
  atomIdsSvg: "/api/v1/name/atom-ids-svg",
  benchmarkPreview: "/api/v1/benchmark-preview",
  benchmarkDatasets: "/api/v1/benchmark-preview/datasets",
  benchmarkDiff: "/api/v1/benchmark-preview/diff",
  benchmarkRefresh: "/api/v1/benchmark-preview/refresh",
  benchmarkStatus: "/api/v1/benchmark-preview/status",
  codeAnalysis: "/api/v1/code-analysis",
  callGraph: "/api/v1/call-graph",
  callGraphProgress: "/api/v1/call-graph/progress",
  callGraphSvg: "/api/v1/call-graph/svg",
  debug: "/api/v1/debug",
  debugPrint: "/api/v1/debug-print",
  debugPrintStream: "/api/v1/debug-print-stream",
  debugPrintBatchStream: "/api/v1/debug-print-batch-stream",
  gitCommits: "/api/v1/git/commits",
};

export const state = {
  currentCommit: null /* null = HEAD; set by the history picker */,
  namerHistory: [],
  ketcherReady: false,
  lastNamedSmiles: "",
  nameReqSeq: 0,
  ketcherBridge: null,
  liveDebounceTimer: null,
  liveSmilesTimer: null,
  suppressSmilesLive: false,
  ketcherMuted: false, // 载入键入 SMILES / 点「从 SMILES 载入」期间静默其 onChange, 避免重排串写回输入框
  ketcherLoadCount: 0, // 进行中的载入数(引用计数); 全部结束后才解除静默
  ketcherMuteTimer: null, // 解除静默的延时器(距最后一次被吞 change 的安静期)
  namerOrient: true, // L4 编号图: 稠环按 preferred_orientation 水平摆放
  engine: "src", // 命名引擎: "src" | "v2" | "v3" | "ml"; namer 页开关会写回并持久化到 localStorage
  // benchmark
  currentPage: "namer",
  bmRows: [],
  bmLoaded: false,
  bmPage: 1,
  bmGenerating: false,
  bmGenDone: 0,
  bmGenTotal: 0,
  bmHacRows: null, // 上次初始化重原子滑块域时的 bmRows 引用; 变了就重置区间
  bmPollTimer: null,
  bmFiles: [], // 与 bmRows 同一次响应里的数据集文件名数组; 行内 f 是它的下标
  bmExFiles: new Set(), // 左栏「按文件排除」勾选的文件名
  bmExFeats: new Set(), // 左栏「按特征排除」勾选的特征名
  bmFeatExpanded: false, // 特征列表是否已「展开全部」
  bmFeatQuery: "", // 特征搜索词(只影响左栏列表, 不影响表格)
  bmSideOpen: true, // 左栏是否展开
  bmExCount: 0, // 上次过滤被排除掉的行数, 供「已排除」chip 显示
  bmDiffs: new Map(), // 已高亮的行: abs(bmRows 下标) → gold_diff; 换数据文件时清空
  bmDiffAll: false, // 「全部差异」开关: 开着时本页每行都自动高亮
  bmDiffPending: new Set(), // 已发出、尚未返回的差异请求行下标, 防重复请求
  // code analysis
  caData: null,
  // call graph
  cgData: null,
  cgThreshold: 0,
  cgLoading: false,
  cgSvg: null,
  cgSvgLoading: false,
  cgSourceTotal: 0,
  cgRank: {
    calls: { asc: false, all: false },
    time: { asc: false, all: false },
  },
};

export const $ = (id) => document.getElementById(id);

export function escapeHtml(str) {
  return String(str)
    .replace(/&/g, "&amp;")
    .replace(/</g, "&lt;")
    .replace(/>/g, "&gt;")
    .replace(/"/g, "&quot;");
}

export async function api(url, opts) {
  // When a historical commit is selected, point the three data pages at it.
  // Endpoints that don't declare ?commit= simply ignore the extra query param.
  if (state.currentCommit && !/[?&]commit=/.test(url)) {
    url += (url.indexOf("?") === -1 ? "?" : "&") + "commit=" + encodeURIComponent(state.currentCommit);
  }
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

/* ---------- 后端健康横幅 ---------- */

const HEALTH_RECHECK_MS = 3000; // 异常期间的重试间隔；正常时不再轮询
let healthTimer = null;

/* 探测 /health 并反映到顶部横幅：连不上的服务、编译不过的 src 都直接说出来，
   而不是让用户对着一个没反应的页面猜。异常期间每 3s 重试，恢复后自动收起。 */
export async function checkBackendHealth() {
  const el = $("backend-banner");
  const text = $("backend-banner-text");
  if (!el || !text) return;
  let msg = "";
  let kind = "";
  try {
    const res = await fetch(API.health, { cache: "no-store" });
    const j = await res.json();
    if (j.src === "broken") {
      kind = "is-degraded";
      msg = "src/ 编译失败，命名接口不可用：" + (j.src_error || "导入 opensmiles 出错");
    }
  } catch (_) {
    kind = "is-down";
    msg = "后端未启动或已断开（" + API.health + "）：页面数据无法加载";
  }
  if (healthTimer) {
    clearTimeout(healthTimer);
    healthTimer = null;
  }
  el.hidden = !msg;
  el.classList.toggle("is-degraded", kind === "is-degraded");
  el.classList.toggle("is-down", kind === "is-down");
  text.textContent = msg;
  if (msg) healthTimer = setTimeout(checkBackendHealth, HEALTH_RECHECK_MS);
}
