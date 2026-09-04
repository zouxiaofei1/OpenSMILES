/* ChemAgent Namer — shared core: API endpoints, global state, DOM/HTTP helpers.
   All feature modules import from here; nothing in core imports back. */
export const API = {
  name: "/api/v1/name",
  pubchemIupac: "/api/v1/name/pubchem-iupac",
  locantsSvg: "/api/v1/name/locants-svg",
  atomIdsSvg: "/api/v1/name/atom-ids-svg",
  benchmarkPreview: "/api/v1/benchmark-preview",
  benchmarkDatasets: "/api/v1/benchmark-preview/datasets",
  benchmarkRefresh: "/api/v1/benchmark-preview/refresh",
  benchmarkStatus: "/api/v1/benchmark-preview/status",
  benchmarkRun: "/api/v1/benchmark-run",
  benchmarkRunDatasets: "/api/v1/benchmark-run/datasets",
  benchmarkRunStatus: "/api/v1/benchmark-run/status",
  benchmarkRunResult: "/api/v1/benchmark-run/result",
  codeAnalysis: "/api/v1/code-analysis",
  callGraph: "/api/v1/call-graph",
  callGraphProgress: "/api/v1/call-graph/progress",
  callGraphSvg: "/api/v1/call-graph/svg",
  debug: "/api/v1/debug",
  debugPrint: "/api/v1/debug-print",
  debugPrintStream: "/api/v1/debug-print-stream",
  settings: "/api/v1/settings",
  wiki: "/api/v1/wiki",
  wikiDoc: "/api/v1/wiki/doc",
  iupac: "/api/v1/iupac",
  iupacDoc: "/api/v1/iupac/doc",
  gitCommits: "/api/v1/git/commits",
};

export const state = {
  currentCommit: null /* null = HEAD; set by the history picker */,
  namerHistory: [],
  liveNameEnabled: true,
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
  brDatasetsLoaded: false,
  // code analysis
  caData: null,
  // call graph
  cgData: null,
  cgThreshold: 0,
  cgLoading: false,
  cgSvg: null,
  cgSvgLoading: false,
  cgSourceTotal: 0,
  docsSource: "wiki", // 文档页当前来源: "wiki" | "iupac"
  wikiCurrent: "", // 文档页 · 各来源记住上次打开的文档
  iupacCurrent: "",
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
