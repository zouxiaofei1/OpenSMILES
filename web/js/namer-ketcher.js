/* ChemAgent Console — Ketcher bridge for Namer */
(function (global) {
  "use strict";

  const DEFAULT_SRC = "/vendor/ketcher/index.html";
  const READY_POLL_MS = 200;
  const READY_TIMEOUT_MS = 20000;

  function shouldSkipLiveName(liveEnabled, smiles, lastNamed) {
    if (!liveEnabled) return true;
    const s = (smiles || "").trim();
    if (!s) return true;
    return s === (lastNamed || "");
  }

  function nextReqSeq(current) {
    const n = Number(current) || 0;
    return n + 1;
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
      if (pollTimer) {
        clearInterval(pollTimer);
        pollTimer = null;
      }
      attachChange(k);
      onReady(k);
    }

    function startPolling() {
      const t0 = Date.now();
      pollTimer = setInterval(function () {
        if (destroyed) return;
        const k = getKetcherFromFrame();
        if (k) {
          markReady(k);
          return;
        }
        if (Date.now() - t0 > READY_TIMEOUT_MS) {
          clearInterval(pollTimer);
          pollTimer = null;
          onError(new Error("Ketcher 加载超时"));
        }
      }, READY_POLL_MS);
    }

    function onFrameLoad() {
      if (destroyed) return;
      const win = iframe.contentWindow;
      if (!win) {
        onError(new Error("Ketcher iframe 不可用"));
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

    function init() {
      if (destroyed) return;
      if (iframe.dataset.ketcherInit === "1") return;
      iframe.dataset.ketcherInit = "1";
      iframe.addEventListener("load", onFrameLoad);
      iframe.addEventListener("error", function () {
        onError(new Error("Ketcher 资源加载失败"));
      });
      fetch(src, { method: "GET", cache: "no-cache" })
        .then(function (res) {
          if (!res.ok) throw new Error("HTTP " + res.status);
          iframe.src = src;
        })
        .catch(function (err) {
          onError(err);
        });
    }

    async function getSmiles() {
      if (!isReady()) return "";
      try {
        const s = await ketcher.getSmiles();
        return (s || "").trim();
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
      if (pollTimer) clearInterval(pollTimer);
      if (smilesPollTimer) clearInterval(smilesPollTimer);
      if (changeUnsub) {
        try {
          changeUnsub();
        } catch (_) {}
      }
      ready = false;
      ketcher = null;
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
