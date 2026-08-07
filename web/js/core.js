/* ChemAgent Namer — shared core: API endpoints, global state, DOM/HTTP helpers.
   All feature modules import from here; nothing in core imports back. */
export const API = {
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
  callGraphProgress: "/api/v1/call-graph/progress",
  callGraphSvg: "/api/v1/call-graph/svg",
  debug: "/api/v1/debug",
  settings: "/api/v1/settings",
  wiki: "/api/v1/wiki",
  wikiDoc: "/api/v1/wiki/doc",
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
