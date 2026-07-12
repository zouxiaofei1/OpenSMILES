/* ChemAgent Console — client */
(function () {
  "use strict";

  const API = {
    state: "/api/v1/loop/state",
    start: "/api/v1/loop/start",
    pause: "/api/v1/loop/pause",
    stop: "/api/v1/loop/stop",
    sessions: "/api/v1/sessions",
    session: (n) => `/api/v1/sessions/${n}`,
    name: "/api/v1/name",
    events: "/api/v1/events",
  };

  const STATUS_LABEL = {
    running: "运行中",
    paused: "已暂停",
    stopping: "停止中",
    stopped: "已停止",
    idle: "空闲",
    gate_pass: "已提交",
    gate_fail: "已回退",
    error: "错误",
  };

  const state = {
    loop: null,
    sessions: [],
    selectedIter: null,
    selectedSession: null,
    dualHistory: [],
    namerHistory: [],
    busy: false,
    sseConnected: false,
    pollTimer: null,
  };

  const $ = (id) => document.getElementById(id);

  function pct(v) {
    if (v == null || Number.isNaN(Number(v))) return "—";
    const n = Number(v);
    if (n <= 1 && n >= 0) return (n * 100).toFixed(2) + "%";
    return n.toFixed(4);
  }

  function fmtNum(v, digits = 4) {
    if (v == null || Number.isNaN(Number(v))) return "—";
    return Number(v).toFixed(digits);
  }

  function statusLabel(s) {
    return STATUS_LABEL[s] || s || "—";
  }

  function setBusy(on) {
    state.busy = on;
    ["btn-start", "btn-pause", "btn-stop", "btn-name"].forEach((id) => {
      const el = $(id);
      if (el) el.disabled = on && id !== "btn-name" ? on : el.disabled;
    });
    if ($("btn-start")) $("btn-start").disabled = on;
    if ($("btn-pause")) $("btn-pause").disabled = on;
    if ($("btn-stop")) $("btn-stop").disabled = on;
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

  /* ---------- Top bar / Overview state ---------- */

  function renderStatus(status) {
    const pill = $("status-pill");
    const text = $("status-text");
    if (!pill || !text) return;
    const s = status || "stopped";
    pill.className = "status-pill status-" + s;
    text.textContent = statusLabel(s);
  }

  function renderState(data) {
    state.loop = data || {};
    const d = state.loop;
    renderStatus(d.status);

    const dual = d.last_dual;
    const best = d.best_dual;
    if ($("kpi-dual")) $("kpi-dual").textContent = pct(dual);
    if ($("kpi-best")) $("kpi-best").textContent = pct(best);
    if ($("kpi-iter")) $("kpi-iter").textContent = String(d.iter ?? 0);
    if ($("kpi-k")) $("kpi-k").textContent = String(d.K ?? 5);

    if ($("ov-dual")) $("ov-dual").textContent = pct(dual);
    if ($("ov-best")) $("ov-best").textContent = pct(best);
    if ($("ov-target")) $("ov-target").textContent = pct(d.target);
    if ($("ov-no-improve")) $("ov-no-improve").textContent = String(d.no_improve ?? 0);
    if ($("ov-iter")) $("ov-iter").textContent = String(d.iter ?? 0);
    if ($("ov-cluster")) $("ov-cluster").textContent = d.last_cluster != null ? String(d.last_cluster) : "—";
    if ($("ov-state-json")) $("ov-state-json").textContent = JSON.stringify(d, null, 2);
  }

  async function fetchState() {
    try {
      const data = await api(API.state);
      renderState(data);
    } catch (err) {
      appendLog("error", { message: "fetchState: " + err.message });
    }
  }

  /* ---------- Sessions ---------- */

  function dualFromSession(s) {
    const after = s && s.bench_after;
    if (after && after.dual != null) return Number(after.dual);
    const before = s && s.bench_before;
    if (before && before.dual != null) return Number(before.dual);
    return null;
  }

  function renderSessionList() {
    const list = $("session-list");
    const empty = $("session-empty");
    const count = $("session-count");
    if (!list) return;

    const sessions = state.sessions.slice().sort((a, b) => Number(b.iter) - Number(a.iter));
    if (count) count.textContent = String(sessions.length);

    if (!sessions.length) {
      list.innerHTML = "";
      const el = document.createElement("div");
      el.className = "empty-state";
      el.id = "session-empty";
      el.innerHTML =
        '<p class="empty-title">暂无循环 · 点击 Start 开始</p>' +
        '<p class="empty-sub">循环启动后，每轮 session 会显示在此</p>';
      list.appendChild(el);
      return;
    }

    list.innerHTML = "";
    sessions.forEach((s) => {
      const btn = document.createElement("button");
      btn.type = "button";
      btn.className = "session-item" + (state.selectedIter === s.iter ? " selected" : "");
      btn.setAttribute("role", "option");
      btn.setAttribute("aria-selected", state.selectedIter === s.iter ? "true" : "false");
      btn.title = `iter ${s.iter} · ${s.status || s.decision || ""}`;

      const st = s.status || s.decision || "unknown";
      const dual = dualFromSession(s);
      btn.innerHTML =
        `<span class="dot ${escapeAttr(st)}" aria-hidden="true"></span>` +
        `<span class="session-meta">` +
        `<span class="session-title mono">#${s.iter}</span>` +
        `<span class="session-sub">${escapeHtml(statusLabel(st))} · ${escapeHtml(s.decision || "—")}</span>` +
        `</span>` +
        `<span class="session-delta">${dual != null ? pct(dual) : "—"}</span>`;

      btn.addEventListener("click", () => selectSession(s.iter));
      list.appendChild(btn);
    });
  }

  function escapeHtml(str) {
    return String(str)
      .replace(/&/g, "&amp;")
      .replace(/</g, "&lt;")
      .replace(/>/g, "&gt;")
      .replace(/"/g, "&quot;");
  }

  function escapeAttr(str) {
    return String(str).replace(/[^a-zA-Z0-9_-]/g, "");
  }

  async function fetchSessions() {
    try {
      const data = await api(API.sessions);
      const list = Array.isArray(data) ? data : [];
      // Enrich summaries with bench duals (sidebar + chart); cap work
      const enriched = await Promise.all(
        list.slice(0, 50).map(async (s) => {
          if (s.bench_after || s.bench_before) return s;
          try {
            const d = await api(API.session(s.iter));
            return {
              ...s,
              status: d.status || s.status,
              decision: d.decision || s.decision,
              bench_before: d.bench_before,
              bench_after: d.bench_after,
            };
          } catch (_) {
            return s;
          }
        })
      );
      state.sessions = enriched;
      state.dualHistory = enriched
        .map((s) => ({ iter: s.iter, dual: dualFromSession(s) }))
        .filter((p) => p.dual != null)
        .sort((a, b) => a.iter - b.iter);
      renderSessionList();
      renderDualChart();
    } catch (err) {
      appendLog("error", { message: "fetchSessions: " + err.message });
    }
  }

  async function selectSession(iterN) {
    state.selectedIter = iterN;
    renderSessionList();
    try {
      const detail = await api(API.session(iterN));
      state.selectedSession = detail;
      // track dual for chart
      const d = dualFromSession(detail);
      if (d != null) {
        const idx = state.dualHistory.findIndex((p) => p.iter === iterN);
        const point = { iter: iterN, dual: d };
        if (idx >= 0) state.dualHistory[idx] = point;
        else state.dualHistory.push(point);
        state.dualHistory.sort((a, b) => a.iter - b.iter);
      }
      renderCycle(detail);
      renderPrompt(detail);
      renderBench(detail);
      renderFailClusters(detail);
      renderDualChart();
    } catch (err) {
      appendLog("error", { message: "selectSession: " + err.message });
    }
  }

  function kvHtml(entries) {
    return entries
      .map(
        ([k, v]) =>
          `<div><dt>${escapeHtml(k)}</dt><dd class="mono">${escapeHtml(v == null ? "—" : String(v))}</dd></div>`
      )
      .join("");
  }

  function renderCycle(s) {
    const empty = $("cycle-empty");
    const content = $("cycle-content");
    if (!s) {
      if (empty) empty.classList.remove("hidden");
      if (content) content.classList.add("hidden");
      return;
    }
    if (empty) empty.classList.add("hidden");
    if (content) content.classList.remove("hidden");

    const st = s.status || s.decision || "—";
    const badge = $("cycle-status-badge");
    if (badge) {
      badge.textContent = statusLabel(st);
      badge.className = "badge " + (st === "gate_pass" || s.decision === "commit" ? "ok" : st === "gate_fail" || s.decision === "revert" ? "fail" : "warn");
    }
    if ($("cycle-iter")) $("cycle-iter").textContent = "#" + s.iter;
    if ($("cycle-decision")) $("cycle-decision").textContent = "decision: " + (s.decision || "—");

    const bb = s.bench_before || {};
    const ba = s.bench_after || {};
    if ($("cycle-bench")) {
      $("cycle-bench").innerHTML = kvHtml([
        ["dual", `${fmtNum(bb.dual)} → ${fmtNum(ba.dual)}`],
        ["en", `${fmtNum(bb.en)} → ${fmtNum(ba.en)}`],
        ["zh", `${fmtNum(bb.zh)} → ${fmtNum(ba.zh)}`],
        ["fails", `${bb.fails ?? "—"} → ${ba.fails ?? "—"}`],
      ]);
    }
    if ($("cycle-gate")) {
      $("cycle-gate").innerHTML = kvHtml([
        ["decision", s.decision],
        ["lint_ok", s.lint_ok],
        ["tdd", s.tdd && s.tdd.passed != null ? s.tdd.passed : "—"],
        ["base_sha", s.base_sha],
        ["commit_sha", s.commit_sha],
        ["layers", Array.isArray(s.layers) ? s.layers.join(", ") : s.layers],
      ]);
    }
    if ($("cycle-cluster")) $("cycle-cluster").textContent = s.cluster_summary || "—";
    if ($("cycle-log")) $("cycle-log").textContent = s.log_excerpt || "—";
  }

  function renderFailClusters(s) {
    const el = $("fail-clusters");
    if (!el) return;
    if (!s || !s.cluster_summary) {
      el.innerHTML = '<p class="muted-text small">选择 session 查看簇摘要</p>';
      return;
    }
    const pre = document.createElement("pre");
    pre.textContent = s.cluster_summary;
    el.innerHTML = "";
    el.appendChild(pre);
  }

  async function renderPrompt(s) {
    const empty = $("prompt-empty");
    const content = $("prompt-content");
    if (!s) {
      if (empty) empty.classList.remove("hidden");
      if (content) content.classList.add("hidden");
      return;
    }

    const path = s.prompt_path || "";
    let body = s.prompt || s.prompt_text || null;

    // Session may only store path; show path + available fields
    if (!body && s.cluster_summary) {
      body =
        `# Cycle ${s.iter} prompt\n\n` +
        `path: ${path || "(none)"}\n\n` +
        `## Cluster\n${s.cluster_summary}\n\n` +
        `## Layers\n${Array.isArray(s.layers) ? s.layers.join(", ") : s.layers || "—"}\n\n` +
        `## Skill\nskills/chem-tdd-skill/\n`;
    }

    if (!body && !path) {
      if (empty) empty.classList.remove("hidden");
      if (content) content.classList.add("hidden");
      return;
    }

    if (empty) empty.classList.add("hidden");
    if (content) content.classList.remove("hidden");
    if ($("prompt-path")) $("prompt-path").textContent = path || "(embedded)";
    if ($("prompt-body")) $("prompt-body").textContent = body || `(prompt file: ${path})`;
  }

  function renderBench(s) {
    const empty = $("bench-empty");
    const content = $("bench-content");
    if (!s || (!s.bench_before && !s.bench_after)) {
      if (empty) empty.classList.remove("hidden");
      if (content) content.classList.add("hidden");
      return;
    }
    if (empty) empty.classList.add("hidden");
    if (content) content.classList.remove("hidden");
    const side = (b) =>
      kvHtml([
        ["dual", fmtNum(b.dual)],
        ["en", fmtNum(b.en)],
        ["zh", fmtNum(b.zh)],
        ["fails", b.fails],
      ]);
    if ($("bench-before")) $("bench-before").innerHTML = side(s.bench_before || {});
    if ($("bench-after")) $("bench-after").innerHTML = side(s.bench_after || {});
  }

  function renderDualChart() {
    const svg = $("dual-chart");
    const empty = $("chart-empty");
    if (!svg) return;

    // Build series from dualHistory; if empty, try sessions with dual fields
    let points = state.dualHistory.slice();
    if (!points.length) {
      points = state.sessions
        .map((s) => ({ iter: s.iter, dual: dualFromSession(s) }))
        .filter((p) => p.dual != null)
        .sort((a, b) => a.iter - b.iter);
    }

    if (points.length < 1) {
      svg.innerHTML = "";
      if (empty) empty.classList.remove("hidden");
      return;
    }
    if (empty) empty.classList.add("hidden");

    const W = 600;
    const H = 160;
    const pad = 16;
    const xs = points.map((p) => p.iter);
    const ys = points.map((p) => p.dual);
    const minX = Math.min(...xs);
    const maxX = Math.max(...xs);
    const minY = Math.min(0, ...ys);
    const maxY = Math.max(1, ...ys);
    const dx = maxX - minX || 1;
    const dy = maxY - minY || 1;

    const coords = points.map((p) => {
      const x = pad + ((p.iter - minX) / dx) * (W - pad * 2);
      const y = H - pad - ((p.dual - minY) / dy) * (H - pad * 2);
      return [x, y];
    });

    const poly = coords.map(([x, y]) => `${x.toFixed(1)},${y.toFixed(1)}`).join(" ");
    const dots = coords
      .map(
        ([x, y], i) =>
          `<circle cx="${x.toFixed(1)}" cy="${y.toFixed(1)}" r="3" fill="#22C55E"><title>#${points[i].iter}: ${pct(points[i].dual)}</title></circle>`
      )
      .join("");

    svg.innerHTML =
      `<polyline fill="none" stroke="#22C55E" stroke-width="2" points="${poly}" />` + dots;
  }

  /* ---------- Controls ---------- */

  async function onStart() {
    setBusy(true);
    try {
      const r = await api(API.start, { method: "POST" });
      appendLog("loop.start", r);
      await fetchState();
      await fetchSessions();
    } catch (err) {
      appendLog("error", { message: err.message });
    } finally {
      setBusy(false);
    }
  }

  async function onPause() {
    setBusy(true);
    try {
      const r = await api(API.pause, { method: "POST" });
      appendLog("loop.pause", r);
      await fetchState();
    } catch (err) {
      appendLog("error", { message: err.message });
    } finally {
      setBusy(false);
    }
  }

  async function onStop() {
    const ok = window.confirm("确认停止 Agent 循环？将写入 STOP 并等待当前轮结束。");
    if (!ok) return;
    const force = window.confirm("是否强制停止（force join）？\n确定 = force，取消 = 普通 stop");
    setBusy(true);
    try {
      const r = await api(API.stop, {
        method: "POST",
        body: JSON.stringify({ force: !!force }),
      });
      appendLog("loop.stop", r);
      await fetchState();
    } catch (err) {
      appendLog("error", { message: err.message });
    } finally {
      setBusy(false);
    }
  }

  /* ---------- Namer ---------- */

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
    const input = $("smiles-input");
    const errEl = $("namer-error");
    const btn = $("btn-name");
    const smiles = (input && input.value || "").trim();
    if (errEl) {
      errEl.hidden = true;
      errEl.textContent = "";
    }
    if (!smiles) {
      if (errEl) {
        errEl.hidden = false;
        errEl.textContent = "请输入 SMILES";
      }
      return;
    }
    if (btn) btn.disabled = true;
    try {
      const result = await api(API.name, {
        method: "POST",
        body: JSON.stringify({ smiles }),
      });
      showNamerResult(result);
      state.namerHistory.unshift({
        smiles,
        en: result.en,
        zh: result.zh,
        success: result.success,
        time_ms: result.time_ms,
      });
      state.namerHistory = state.namerHistory.slice(0, 20);
      renderNamerHistory();
      appendLog("name", { smiles, en: result.en, zh: result.zh, success: result.success });
    } catch (err) {
      if (errEl) {
        errEl.hidden = false;
        errEl.textContent = err.message;
      }
      appendLog("error", { message: "name: " + err.message });
    } finally {
      if (btn) btn.disabled = false;
    }
  }

  function showNamerResult(result) {
    const box = $("namer-result");
    if (!box) return;
    box.classList.remove("hidden");
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

  /* ---------- Logs / SSE ---------- */

  function appendLog(type, payload) {
    const stream = $("log-stream");
    if (!stream) return;
    const ts = new Date().toISOString().slice(11, 19);
    const line = document.createElement("span");
    line.className = "log-line";
    const data =
      typeof payload === "string" ? payload : JSON.stringify(payload ?? {}, null, 0);
    line.innerHTML =
      `<span class="log-ts">${ts}</span> ` +
      `<span class="log-type">${escapeHtml(type)}</span> ` +
      escapeHtml(data);
    stream.appendChild(line);
    stream.appendChild(document.createTextNode("\n"));

    // cap lines
    while (stream.childNodes.length > 800) {
      stream.removeChild(stream.firstChild);
    }

    const auto = $("logs-autoscroll");
    if (!auto || auto.checked) {
      stream.scrollTop = stream.scrollHeight;
    }
  }

  function setSseBanner(show) {
    const el = $("sse-banner");
    if (!el) return;
    el.classList.toggle("hidden", !show);
  }

  function connectSSE() {
    let es;
    try {
      es = new EventSource(API.events);
    } catch (err) {
      setSseBanner(true);
      return;
    }

    es.onopen = () => {
      state.sseConnected = true;
      setSseBanner(false);
      appendLog("sse", { message: "connected" });
    };

    es.onerror = () => {
      state.sseConnected = false;
      setSseBanner(true);
    };

    // Generic message handler + named events
    function handleEvent(ev) {
      let data = {};
      try {
        data = JSON.parse(ev.data);
      } catch (_) {
        data = { raw: ev.data };
      }
      const type = data.type || ev.type || "message";
      const payload = data.payload != null ? data.payload : data;

      if (type === "log.line" || type === "log_line") {
        appendLog("log.line", payload);
      } else if (type === "cycle_end" || type === "cycle.updated" || type === "cycle_start") {
        appendLog(type, payload);
        fetchState();
        fetchSessions();
      } else if (type === "loop.state" || type === "state") {
        if (payload && typeof payload === "object") renderState({ ...state.loop, ...payload });
        else fetchState();
        appendLog(type, payload);
      } else {
        appendLog(type, payload);
      }
    }

    es.onmessage = handleEvent;
    // Listen for typed SSE events if server uses event: field
    ["log.line", "cycle_start", "cycle_end", "cycle.updated", "loop.state", "bench.done", "message"].forEach(
      (name) => {
        es.addEventListener(name, handleEvent);
      }
    );
  }

  /* ---------- Tabs ---------- */

  function switchTab(name) {
    document.querySelectorAll(".tab").forEach((t) => {
      const on = t.dataset.tab === name;
      t.classList.toggle("active", on);
      t.setAttribute("aria-selected", on ? "true" : "false");
    });
    document.querySelectorAll(".tab-panel").forEach((p) => {
      const on = p.id === "panel-" + name;
      p.classList.toggle("active", on);
      if (on) p.removeAttribute("hidden");
      else p.setAttribute("hidden", "");
    });
  }

  /* ---------- Init ---------- */

  function bind() {
    $("btn-start") && $("btn-start").addEventListener("click", onStart);
    $("btn-pause") && $("btn-pause").addEventListener("click", onPause);
    $("btn-stop") && $("btn-stop").addEventListener("click", onStop);
    $("namer-form") && $("namer-form").addEventListener("submit", onName);
    $("btn-clear-history") &&
      $("btn-clear-history").addEventListener("click", () => {
        state.namerHistory = [];
        renderNamerHistory();
      });
    $("btn-clear-logs") &&
      $("btn-clear-logs").addEventListener("click", () => {
        if ($("log-stream")) $("log-stream").textContent = "";
      });
    $("btn-copy-prompt") &&
      $("btn-copy-prompt").addEventListener("click", async () => {
        const text = $("prompt-body") ? $("prompt-body").textContent : "";
        try {
          await navigator.clipboard.writeText(text || "");
          appendLog("ui", { message: "prompt copied" });
        } catch (err) {
          appendLog("error", { message: "copy failed: " + err.message });
        }
      });

    document.querySelectorAll(".tab").forEach((tab) => {
      tab.addEventListener("click", () => switchTab(tab.dataset.tab));
    });
  }

  async function init() {
    bind();
    renderNamerHistory();
    await fetchState();
    await fetchSessions();
    connectSSE();
    // Fallback poll every 2s for state; sessions every 4s
    state.pollTimer = setInterval(() => {
      fetchState();
    }, 2000);
    setInterval(fetchSessions, 4000);
    appendLog("ui", { message: "ChemAgent Console ready" });
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", init);
  } else {
    init();
  }
})();
