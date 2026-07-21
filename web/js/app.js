/* ChemAgent Namer — client */
(function () {
  "use strict";

  const API = {
    name: "/api/v1/name",
  };

  const state = {
    namerHistory: [],
    liveNameEnabled: false,
    ketcherReady: false,
    lastNamedSmiles: "",
    nameReqSeq: 0,
    ketcherBridge: null,
    liveDebounceTimer: null,
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
        if (input) input.value = smiles;
      }
      await runName(smiles, { fromLive: true });
    }, 700);
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
        if (!state.liveNameEnabled && state.liveDebounceTimer) {
          clearTimeout(state.liveDebounceTimer);
          state.liveDebounceTimer = null;
        }
      });
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
