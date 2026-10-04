/* OpenSMILES Console — Ketcher bridge for Namer */
(function (global) {
  "use strict";

  const DEFAULT_SRC = "/vendor/ketcher/index.html";
  const READY_POLL_MS = 200;
  // Cold load pulls ~23MB main + ~11MB worker; allow slower disks/networks.
  const READY_TIMEOUT_MS = 60000;

  function shouldSkipLiveName(smiles, lastNamed) {
    const s = (smiles || "").trim();
    if (!s) return true;
    return s === (lastNamed || "");
  }

  function nextReqSeq(current) {
    const n = Number(current) || 0;
    return n + 1;
  }

  /** Append cache-bust query so iframe navigation cannot reuse a broken HTML shell. */
  function withCacheBust(url) {
    const base = String(url || DEFAULT_SRC);
    const sep = base.indexOf("?") >= 0 ? "&" : "?";
    return base + sep + "v=" + Date.now();
  }

  function createBridge(opts) {
    const iframe = opts.iframe;
    const src = opts.src || DEFAULT_SRC;
    const onReady = opts.onReady || function () {};
    const onError = opts.onError || function () {};
    const onChange = opts.onChange || function () {};

    let ready = false;
    let ketcher = null;
    let pollTimer = null;
    let changeUnsub = null;
    let destroyed = false;
    let lastPolledSmiles = null;
    let smilesPollTimer = null;
    let loadHandler = null;
    let errorHandler = null;
    let failed = false;

    function isReady() {
      return ready && !!ketcher;
    }

    function getKetcherFromFrame() {
      try {
        return iframe.contentWindow && iframe.contentWindow.ketcher;
      } catch (_) {
        return null;
      }
    }

    function clearPoll() {
      if (pollTimer) {
        clearInterval(pollTimer);
        pollTimer = null;
      }
    }

    function attachChange(k) {
      if (!k) return;
      const handler = function () {
        onChange();
      };
      try {
        if (k.changeEvent) {
          if (typeof k.changeEvent.add === "function") {
            k.changeEvent.add(handler);
            changeUnsub = function () {
              if (typeof k.changeEvent.remove === "function") k.changeEvent.remove(handler);
            };
            return;
          }
          if (typeof k.changeEvent.subscribe === "function") {
            const sub = k.changeEvent.subscribe(handler);
            changeUnsub = function () {
              if (sub && typeof sub.unsubscribe === "function") sub.unsubscribe();
            };
            return;
          }
        }
      } catch (_) {
        /* fall through to poll */
      }
      // Fallback: poll getSmiles while bridge is alive
      smilesPollTimer = setInterval(async function () {
        if (destroyed || !isReady()) return;
        try {
          const s = await k.getSmiles();
          const t = (s || "").trim();
          if (t !== lastPolledSmiles) {
            lastPolledSmiles = t;
            onChange();
          }
        } catch (_) {
          /* ignore */
        }
      }, 1000);
      changeUnsub = function () {
        if (smilesPollTimer) {
          clearInterval(smilesPollTimer);
          smilesPollTimer = null;
        }
      };
    }

    function markReady(k) {
      if (destroyed || ready) return;
      ketcher = k;
      ready = true;
      failed = false;
      clearPoll();
      attachChange(k);
      onReady(k);
    }

    function fail(err) {
      if (destroyed || ready) return;
      failed = true;
      clearPoll();
      // Allow a later init()/retry to run again.
      try {
        delete iframe.dataset.ketcherInit;
      } catch (_) {
        iframe.dataset.ketcherInit = "";
      }
      onError(err);
    }

    function startPolling() {
      clearPoll();
      const t0 = Date.now();
      pollTimer = setInterval(function () {
        if (destroyed) return;
        const k = getKetcherFromFrame();
        if (k) {
          markReady(k);
          return;
        }
        if (Date.now() - t0 > READY_TIMEOUT_MS) {
          fail(new Error("Ketcher 加载超时（可点「重试」或硬刷新；首次需下载约 35MB）"));
        }
      }, READY_POLL_MS);
    }

    function onFrameLoad() {
      if (destroyed) return;
      const win = iframe.contentWindow;
      if (!win) {
        fail(new Error("Ketcher iframe 不可用"));
        return;
      }
      const existing = getKetcherFromFrame();
      if (existing) {
        markReady(existing);
        return;
      }
      function onKetcherReady() {
        const k = getKetcherFromFrame();
        if (k) markReady(k);
      }
      try {
        win.addEventListener("ketcher-ready", onKetcherReady);
      } catch (_) {
        /* ignore */
      }
      startPolling();
    }

    function detachFrameHandlers() {
      if (loadHandler) {
        try {
          iframe.removeEventListener("load", loadHandler);
        } catch (_) {}
        loadHandler = null;
      }
      if (errorHandler) {
        try {
          iframe.removeEventListener("error", errorHandler);
        } catch (_) {}
        errorHandler = null;
      }
    }

    function init() {
      if (destroyed) return;
      if (ready) return;
      // Already loading
      if (iframe.dataset.ketcherInit === "1" && !failed) return;

      failed = false;
      ready = false;
      ketcher = null;
      clearPoll();
      if (changeUnsub) {
        try {
          changeUnsub();
        } catch (_) {}
        changeUnsub = null;
      }
      detachFrameHandlers();

      iframe.dataset.ketcherInit = "1";
      loadHandler = onFrameLoad;
      errorHandler = function () {
        fail(new Error("Ketcher 资源加载失败"));
      };
      iframe.addEventListener("load", loadHandler);
      iframe.addEventListener("error", errorHandler);

      // Cache-bust the iframe document itself (fetch no-cache does not apply to iframe.src).
      const busted = withCacheBust(src);
      // Probe that vendor shell exists; then navigate iframe with busted URL.
      fetch(src, { method: "GET", cache: "no-store" })
        .then(function (res) {
          if (!res.ok) throw new Error("HTTP " + res.status);
          return res.text();
        })
        .then(function (html) {
          if (html.indexOf("window.global") < 0) {
            throw new Error(
              "Ketcher shell 缺少 global shim；请运行 bash tools/build_ketcher_vendor.sh 或检查 server/web/vendor/ketcher/index.html"
            );
          }
          iframe.src = busted;
        })
        .catch(function (err) {
          fail(err);
        });
    }

    async function getSmiles() {
      if (!isReady()) return "";
      try {
        const s = await ketcher.getSmiles();
        const trimmed = (s || "").trim();
        // Workaround for Indigo bug: highly symmetric aromatics (e.g. benzene)
        // occasionally produce duplicated SMILES like "C1C=CC=CC=1.C1C=CC=CC=1".
        // When every fragment is identical, keep only the first one.
        if (trimmed.includes(".")) {
          const parts = trimmed.split(".");
          if (parts.length > 1 && parts.every(function (p) { return p === parts[0]; })) {
            return parts[0];
          }
        }
        return trimmed;
      } catch (_) {
        return "";
      }
    }

    async function setMolecule(smiles) {
      if (!isReady()) throw new Error("Ketcher 未就绪");
      await ketcher.setMolecule(smiles || "");
    }

    async function clear() {
      if (!isReady()) return;
      try {
        if (ketcher.editor && typeof ketcher.editor.clear === "function") {
          ketcher.editor.clear();
          return;
        }
      } catch (_) {
        /* fall through */
      }
      await ketcher.setMolecule("");
    }

    function destroy() {
      destroyed = true;
      clearPoll();
      if (smilesPollTimer) clearInterval(smilesPollTimer);
      if (changeUnsub) {
        try {
          changeUnsub();
        } catch (_) {}
      }
      detachFrameHandlers();
      ready = false;
      ketcher = null;
      try {
        delete iframe.dataset.ketcherInit;
      } catch (_) {}
    }

    return {
      init: init,
      isReady: isReady,
      getSmiles: getSmiles,
      setMolecule: setMolecule,
      clear: clear,
      destroy: destroy,
    };
  }

  global.ChemNamerKetcher = {
    createBridge: createBridge,
    shouldSkipLiveName: shouldSkipLiveName,
    nextReqSeq: nextReqSeq,
    DEFAULT_SRC: DEFAULT_SRC,
  };
})(typeof window !== "undefined" ? window : globalThis);
