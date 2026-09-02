/* Namer page: Ketcher bridge + live naming + history. */
import { $, state, api, API, escapeHtml } from "./core.js";

export const LIVE_NAME_DEBOUNCE_MS = 20;
// 载入引发的 change 全部安静多久后解除静默。必须覆盖 setMolecule 的布局
// worker 异步 change; 太小则一次迟到 change 会泄漏写回输入框。
const KETCHER_MUTE_TAIL_MS = 250;

/* Ketcher onChange 静默(引用计数版)。
   「输入框是权威来源」: Ketcher 会把同一分子重排成不同原子序的 SMILES(如
   kekulé→芳香式), 若据此用 getSmiles() 写回输入框并命名, 会覆盖用户输入、
   引入另一套 atom-ids/locants 索引。键入/点「从 SMILES 载入」都会触发
   setMolecule → onChange, 需要把该 change 静默掉。
   但一个载入的 change 是异步的、且可能跨越多次并发载入(连续点击), 固定
   超时(旧的 +100ms)会被并发载入击穿: 先到载入的超时把还没结束的后到载入
   的静默提前关掉。因此用 引用计数(进行中的载入) + "距最后一次被吞 change
   的安静期" 判定结束, 只有最后一次载入结束且再无新 change 才解除静默。 */
function beginKetcherMute() {
  state.ketcherMuted = true;
  state.ketcherLoadCount = (state.ketcherLoadCount || 0) + 1;
  clearKetcherMuteTail();
}

function endKetcherMute() {
  state.ketcherLoadCount = Math.max(0, (state.ketcherLoadCount || 0) - 1);
  if (state.ketcherLoadCount === 0) armKetcherMuteTail();
}

function armKetcherMuteTail() {
  clearKetcherMuteTail();
  state.ketcherMuteTimer = setTimeout(function () {
    state.ketcherMuteTimer = null;
    // 仅当已无进行中的载入才解除: 并发载入时由 count 把关, 先到载入的结束
    // 不会提前解除后到载入的静默。
    if ((state.ketcherLoadCount || 0) === 0) state.ketcherMuted = false;
  }, KETCHER_MUTE_TAIL_MS);
}

function clearKetcherMuteTail() {
  if (state.ketcherMuteTimer) {
    clearTimeout(state.ketcherMuteTimer);
    state.ketcherMuteTimer = null;
  }
}

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

export function ensureKetcher(forceRetry) {
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
  // 输入框是权威来源: Ketcher 会把同一分子重排成不同原子序的 SMILES
  // (如 kekulé→芳香式), 用它命名/渲染会让 atom-ids/locants 索引与
  // 用户输入的 SMILES(以及 /debug 对同一输入的结果)不一致。
  // 因此输入框有内容时优先采用, 仅在画布手绘(输入框为空)时才回退到 Ketcher。
  const typed = (($("smiles-input") && $("smiles-input").value) || "").trim();
  if (typed) return typed;
  let fromEditor = "";
  if (state.ketcherBridge && state.ketcherBridge.isReady()) {
    fromEditor = await state.ketcherBridge.getSmiles();
  }
  if (fromEditor) {
    const input = $("smiles-input");
    if (input) input.value = fromEditor;
    return fromEditor;
  }
  return "";
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
    if (result.success) {
      updateNamerSvg("namer-locants", API.locantsSvg, smiles, seq, state.namerOrient);
      updateNamerSvg("namer-atom-ids", API.atomIdsSvg, smiles, seq);
    }
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
  // 载入键入 SMILES 到 Ketcher 会触发 onChange: 此时结构并未真正变化,
  // 若据此用 getSmiles()(重排后的串)重命名会覆盖用户输入、引入另一套原子序。
  // 该 change 被静默吞掉, 但它是"还有载入在异步收尾"的信号 → 顺延静默期。
  if (state.ketcherMuted) {
    armKetcherMuteTail();
    return;
  }
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
      if (input) {
        state.suppressSmilesLive = true;
        input.value = smiles;
      }
    }
    await runName(smiles, { fromLive: true });
  }, LIVE_NAME_DEBOUNCE_MS);
}

function scheduleLiveSmilesName() {
  if (!state.liveNameEnabled) return;
  if (state.suppressSmilesLive) {
    state.suppressSmilesLive = false;
    return;
  }
  if (state.liveSmilesTimer) clearTimeout(state.liveSmilesTimer);
  state.liveSmilesTimer = setTimeout(async () => {
    state.liveSmilesTimer = null;
    if (!state.liveNameEnabled) return;
    const smiles = (($("smiles-input") && $("smiles-input").value) || "").trim();
    if (!smiles) return;
    const CK = window.ChemNamerKetcher;
    if (CK && CK.shouldSkipLiveName(true, smiles, state.lastNamedSmiles)) return;
    // Load SMILES into Ketcher. 静默其 onChange(见 scheduleLiveName), 否则
    // Ketcher 会把它重排成不同原子序的 SMILES 写回输入框并重命名, 导致
    // atom-ids/locants 用另一套索引(与用户输入及 /debug 不一致)。
    // begin/end 配对维持引用计数: 只 await setMolecule(载入引发的 change 都
    // 紧随其后), 命名本身不改动画布, 不必延长静默。
    beginKetcherMute();
    try {
      if (state.ketcherBridge && state.ketcherBridge.isReady()) {
        try {
          await state.ketcherBridge.setMolecule(smiles);
        } catch (_) {
          /* ignore ketcher load errors during live input */
        }
      }
    } finally {
      endKetcherMute();
    }
    // Run naming
    await runName(smiles, { fromLive: true });
  }, LIVE_NAME_DEBOUNCE_MS);
}

export function renderNamerHistory() {
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
  // Force re-trigger reveal animation
  box.classList.remove("hidden");
  box.style.animation = "none";
  box.offsetHeight; // force reflow
  box.style.animation = "";
  const badge = $("namer-success");
  if (badge) {
    badge.textContent = result.success ? "success" : "failed";
    badge.className = "badge " + (result.success ? "ok" : "fail");
  }
  if ($("namer-time")) $("namer-time").textContent = Math.round(result.time_ms || 0) + " ms";
  if ($("namer-en")) $("namer-en").textContent = result.en || "—";
  if ($("namer-zh")) $("namer-zh").textContent = result.zh || "—";
  if ($("namer-source")) $("namer-source").textContent = result.source || "—";
  // 命名失败时清空并隐藏结构图,避免旧结构残留
  if (!result.success) {
    hideNamerSvg("namer-locants");
    hideNamerSvg("namer-atom-ids");
  }
}

function hideNamerSvg(boxId) {
  const box = $(boxId);
  if (box) {
    box.hidden = true;
    box.innerHTML = "";
  }
}

async function updateNamerSvg(boxId, endpoint, smiles, seq, orient) {
  const box = $(boxId);
  if (!box) return;
  box.hidden = true;
  box.innerHTML = "";
  if (!smiles) return;
  const payload = { smiles };
  if (orient !== undefined) payload.orient = orient;
  try {
    const res = await api(endpoint, {
      method: "POST",
      body: JSON.stringify(payload),
    });
    if (seq !== state.nameReqSeq) return;
    if (res && res.ok && res.svg) {
      box.innerHTML = res.svg;
      box.hidden = false;
    }
  } catch (_) {
    /* 结构图失败静默,不干扰命名结果 */
  }
}

export function bindNamer() {
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
        // 引用计数 + 安静期判定解除静默, 连续点击(并发载入)不会互相击穿。
        // 见 KETCHER_MUTE_TAIL_MS 处注释。
        beginKetcherMute();
        await state.ketcherBridge.setMolecule(smiles);
      } catch (err) {
        setNamerError(err.message || "载入结构失败");
      } finally {
        endKetcherMute();
      }
    });
  $("btn-clear-ketcher") &&
    $("btn-clear-ketcher").addEventListener("click", async () => {
      if (!state.ketcherBridge || !state.ketcherBridge.isReady()) return;
      try {
        await state.ketcherBridge.clear();
        // Also clear the SMILES input
        const input = $("smiles-input");
        if (input) input.value = "";
        setNamerError("");
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
      if (!state.liveNameEnabled) {
        if (state.liveDebounceTimer) {
          clearTimeout(state.liveDebounceTimer);
          state.liveDebounceTimer = null;
        }
        if (state.liveSmilesTimer) {
          clearTimeout(state.liveSmilesTimer);
          state.liveSmilesTimer = null;
        }
      }
    });

  // Live SMILES input: when user edits the SMILES text box, auto-name + load into Ketcher
  $("smiles-input") &&
    $("smiles-input").addEventListener("input", () => {
      scheduleLiveSmilesName();
    });

  // L4 编号图「水平行」开关: 切换后按新取向重新生成结构图
  $("namer-orient-toggle") &&
    $("namer-orient-toggle").addEventListener("change", (ev) => {
      state.namerOrient = !!(ev.target && ev.target.checked);
      if (!state.lastNamedSmiles) return;
      const seq = ++state.nameReqSeq;
      updateNamerSvg("namer-locants", API.locantsSvg, state.lastNamedSmiles, seq, state.namerOrient);
    });

  renderNamerHistory();
}
