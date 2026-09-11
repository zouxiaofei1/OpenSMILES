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

/* 输入可能是 merged_benchmark id(chebi-1 / tiers-1), 判断其形态。
   真正的 SMILES 不可能整串匹配 ^[A-Za-z0-9]+-\d+$; 非 id 无需网络往返。 */
function isBenchmarkId(text) {
  return /^[A-Za-z0-9]+-\d+$/.test(text);
}

/* 把可能为 merged_benchmark id 的输入解析成该行的权威 smiles 串再往下走:
   命中 → 真实 SMILES; 否则原样返回(坏输入照旧交给 namer 报「无法解析」)。
   这样命名 / SVG / history 与直接输入那条 SMILES 完全一致, 不引入第二套索引。 */
async function resolveTextToSmiles(text) {
  if (!isBenchmarkId(text)) return text;
  try {
    const r = await api(API.resolveName, {
      method: "POST",
      body: JSON.stringify({ text }),
    });
    if (r && r.ok && r.kind === "id" && r.smiles) return r.smiles;
  } catch (_) {
    /* 网络失败按原样当 SMILES 走, 不阻断命名 */
  }
  return text;
}

async function resolveSmilesForName() {
  // 输入框是权威来源: Ketcher 会把同一分子重排成不同原子序的 SMILES
  // (如 kekulé→芳香式), 用它命名/渲染会让 atom-ids/locants 索引与
  // 用户输入的 SMILES(以及 /debug 对同一输入的结果)不一致。
  // 因此输入框有内容时优先采用, 仅在画布手绘(输入框为空)时才回退到 Ketcher。
  const typed = (($("smiles-input") && $("smiles-input").value) || "").trim();
  if (typed) return resolveTextToSmiles(typed);
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
      body: JSON.stringify({ smiles, engine: state.engine || "src" }),
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

/* 输入框 SMILES 变化 → 自动执行「从 SMILES 载入」(同步画板结构)。画板同步独立于
   「实时命名」开关: 关掉开关后输入变化仍会载入画板, 只是不再自动跑命名。 */
function scheduleLiveSmilesName() {
  if (state.suppressSmilesLive) {
    // 本次 input 事件来自画布写回输入框(scheduleLiveName): 输入框与画板本就一致,
    // 不必再 setMolecule 回灌(runName 已由写回路径执行)。
    state.suppressSmilesLive = false;
    return;
  }
  if (state.liveSmilesTimer) clearTimeout(state.liveSmilesTimer);
  state.liveSmilesTimer = setTimeout(async () => {
    state.liveSmilesTimer = null;
    const typed = (($("smiles-input") && $("smiles-input").value) || "").trim();
    if (!typed) return;
    // 用户可输 merged_benchmark id(chebi-1 / tiers-1): 先解析成真实 SMILES,
    // 之后载入画板 / 命名 / SVG 都用这一个权威串, 与直接输入该 SMILES 一致。
    const smiles = await resolveTextToSmiles(typed);
    if (!smiles) return;
    const CK = window.ChemNamerKetcher;
    // 1) 画板同步(不依赖「实时命名」开关)。载入会让 Ketcher 把结构重排成不同
    //    原子序的 SMILES 并触发 onChange, 需静默(见 scheduleLiveName 注释), 否则会
    //    写回输入框并引入另一套 atom-ids/locants 索引。begin/end 配对引用计数。
    if (state.ketcherBridge && state.ketcherBridge.isReady()) {
      let needLoad = true;
      try {
        // 画板已是该结构(含 Ketcher 重排后的等价串)则跳过, 避免无谓 setMolecule
        const cur = await state.ketcherBridge.getSmiles();
        needLoad = (cur || "") !== smiles;
      } catch (_) {
        needLoad = true;
      }
      if (needLoad) {
        beginKetcherMute();
        try {
          await state.ketcherBridge.setMolecule(smiles);
        } catch (_) {
          /* 逐键输入时中间态(如 "CC(")可能解析失败: 忽略, 等下一次输入 */
        } finally {
          endKetcherMute();
        }
      }
    }
    // 2) 命名仅当「实时命名」开启时执行
    if (!state.liveNameEnabled) return;
    if (CK && CK.shouldSkipLiveName(true, smiles, state.lastNamedSmiles)) return;
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
  // 新命名/新分子先清掉旧 PubChem 结果,避免残留给上一个结构
  hidePubchem();
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
  const engTag = $("namer-engine-tag");
  if (engTag) {
    if (result.engine) {
      engTag.textContent = result.engine;
      engTag.hidden = false;
    } else {
      engTag.hidden = true;
    }
  }
  if ($("namer-time")) $("namer-time").textContent = Math.round(result.time_ms || 0) + " ms";
  // 命中 merged_benchmark 且相似度达阈值时, 引擎名就地高亮与基准名的差异
  const gd = result.success ? result.gold_diff : null;
  setDiffHtml("namer-en", result.en, diffSide(gd, "en", "pred"));
  setDiffHtml("namer-zh", result.zh, diffSide(gd, "zh", "pred"));
  if ($("namer-source")) $("namer-source").textContent = result.source || "—";
  // 基准答案(merged_benchmark): 命中才显示,失败/未命中一律隐藏
  renderNamerGold(result.success ? result.gold : null, result);
  // 命名失败时清空并隐藏结构图,避免旧结构残留
  if (!result.success) {
    hideNamerSvg("namer-locants");
    hideNamerSvg("namer-atom-ids");
  }
}

/* 结果卡内显示 merged_benchmark 命中记录(gold): 无命中/传入 null 时隐藏并清空。
   result.gold_diff 里相似度达阈值的语言, gold 名称直接在原有行内逐字符高亮与引擎名
   的差异; 未达阈值(或无 gold_diff)照旧显示纯文本。 */
function renderNamerGold(gold, result) {
  const box = $("namer-gold");
  if (!box) return;
  if (!gold) {
    box.hidden = true;
    return;
  }
  const gd = (result && result.gold_diff) || null;
  setDiffHtml("namer-gold-en", gold.en, diffSide(gd, "en", "gold"));
  setDiffHtml("namer-gold-zh", gold.zh, diffSide(gd, "zh", "gold"));
  renderGoldSims(gd);
  const tier = $("namer-gold-tier");
  if (tier) {
    if (gold.tier !== undefined && gold.tier !== null) {
      tier.hidden = false;
      tier.textContent = "tier " + gold.tier;
    } else {
      tier.hidden = true;
    }
  }
  if ($("namer-gold-source")) {
    const src = gold.source || "";
    $("namer-gold-source").textContent = src ? "· " + src : "";
  }
  box.hidden = false;
}

/* 逐字符差异片段 → HTML: 一致片段原样, 差异片段用 <mark> 高亮。 */
function diffSegmentsHtml(segs) {
  return (segs || [])
    .map((s) => {
      const text = escapeHtml(s.t);
      return s.same ? text : '<mark class="gold-diff-mark">' + text + "</mark>";
    })
    .join("");
}

/* 取某语言某一侧的差异片段: 仅当后端判定相似度达阈值(gold_diff.show)时返回, 否则返回
   null, 让调用方回落到纯文本(差异过大时逐字标红没有信息量)。 */
function diffSide(goldDiff, lang, side) {
  const d = goldDiff && goldDiff[lang];
  return d && d.show ? d[side] : null;
}

/* 把名称写进结果行的 dd: 有差异片段则逐字符渲染(差异高亮), 否则纯文本。 */
function setDiffHtml(id, text, segs) {
  const el = $(id);
  if (!el) return;
  el.innerHTML = segs && segs.length ? diffSegmentsHtml(segs) : escapeHtml(text || "—");
}

/* 相似度条: 引擎名 ↔ 基准名的逐字符相似度(0-100), 有基准名的语言各一条。 */
function renderGoldSims(goldDiff) {
  const host = $("namer-gold-sims");
  if (!host) return;
  const parts = [];
  ["en", "zh"].forEach((lang) => {
    const d = goldDiff && goldDiff[lang];
    if (!d) return;
    const pct = Math.max(0, Math.min(100, Math.round((d.similarity || 0) * 100)));
    parts.push(
      '<span class="gold-sim">' +
        '<span class="gold-sim-lang">' + lang + "</span>" +
        '<span class="gold-sim-track">' +
        '<span class="gold-sim-fill" style="width:' + pct + '%"></span>' +
        "</span>" +
        '<span class="gold-sim-pct">' + pct + "%</span>" +
        "</span>"
    );
  });
  host.innerHTML = parts.join("");
  host.hidden = parts.length === 0;
}

/* 隐藏并清空 PubChem IUPAC 结果块(换分子 / 查询失败 / 命名失败时调用)。 */
function hidePubchem() {
  const box = $("namer-pubchem");
  if (box) box.hidden = true;
}

/* 渲染 PubChem 2.1.1 查询结果: 命中显示 iupac + CID 链接;未收录给提示。 */
function renderPubchem(res) {
  const box = $("namer-pubchem");
  if (!box) return;
  const note = $("namer-pubchem-note");
  const link = $("namer-pubchem-link");
  if (res && res.ok) {
    const nameEl = $("namer-pubchem-name");
    if (nameEl) nameEl.textContent = (res.iupac || "").trim() || "—";
    if (link) {
      if (res.cid) {
        link.hidden = false;
        link.textContent = "PubChem CID " + res.cid;
        link.href = res.url || "https://pubchem.ncbi.nlm.nih.gov/";
      } else {
        link.hidden = true;
      }
    }
    if (note) {
      let msg = "";
      if (!res.found) msg = "PubChem 未收录可用的 2.1.1 IUPAC 名称";
      else if (res.racemic) msg = "精确立体无记录，已按非立体（外消旋）形式匹配";
      note.hidden = !msg;
      note.textContent = msg;
    }
    box.hidden = false;
  } else {
    // 后端明确失败(ok:false / 无法解析 / 网络错误)→ 块隐藏,错误走 namer-error
    hidePubchem();
    if (note) note.hidden = true;
    if (res && res.error) setNamerError(res.error);
  }
}

async function onPubchem(ev) {
  ev.preventDefault();
  const smiles = await resolveSmilesForName();
  if (!smiles) {
    setNamerError("请输入或绘制 SMILES 再查 PubChem");
    return;
  }
  const panel = $("namer-result");
  // 结果卡尚未展示或命中的分子不是当前输入 → 先命名,让面板对应本分子(顺带刷新 gold)
  if (!panel || panel.classList.contains("hidden") || state.lastNamedSmiles !== smiles) {
    await runName(smiles, { fromLive: false });
  }
  const btn = $("btn-pubchem");
  const origLabel = btn ? btn.textContent : "";
  if (btn) {
    btn.disabled = true;
    btn.textContent = "查询中…";
  }
  try {
    const res = await api(API.pubchemIupac, {
      method: "POST",
      body: JSON.stringify({ smiles }),
    });
    renderPubchem(res);
  } catch (err) {
    hidePubchem();
    setNamerError(err.message || String(err));
  } finally {
    if (btn) {
      btn.disabled = false;
      btn.textContent = origLabel;
    }
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

/* 引擎开关(src/v2/v3/ml): 设 active 态 + 持久化。仅影响 /name 的 en/zh 输出;
   gold/PubChem/L4 编号/Atom 索引等结构图仍走原(src)管线, 与引擎无关。
   ml 为本地神经网络(仅英文), 每次请求 CPU 生成约几百 ms。 */
const ENGINE_VALUES = ["src", "v2", "v3", "ml"];

function applyEngineUI(engine) {
  state.engine = ENGINE_VALUES.indexOf(engine) >= 0 ? engine : "src";
  document.querySelectorAll(".namer-engine-btn").forEach(function (b) {
    b.classList.toggle("active", b.getAttribute("data-engine") === state.engine);
  });
  try {
    localStorage.setItem("namer.engine", state.engine);
  } catch (_) {
    /* localStorage 不可用时仅本次会话生效 */
  }
}

function bindEngineSwitch() {
  const btns = document.querySelectorAll(".namer-engine-btn");
  if (!btns.length) return;
  // 恢复上次选择; 首次访问缺省 src
  let stored = "src";
  try {
    const v = localStorage.getItem("namer.engine");
    if (ENGINE_VALUES.indexOf(v) >= 0) stored = v;
  } catch (_) {
    /* 忽略读取失败 */
  }
  applyEngineUI(stored);
  btns.forEach(function (b) {
    b.addEventListener("click", async function () {
      const next =
        ENGINE_VALUES.indexOf(b.getAttribute("data-engine")) >= 0
          ? b.getAttribute("data-engine")
          : state.engine;
      if (next === state.engine) return;
      applyEngineUI(next);
      // 当前已有分子(输入框/画板/上次命名)则立刻用新引擎重命名, 便于即时对照
      if (state.lastNamedSmiles) {
        const smiles = await resolveSmilesForName();
        if (smiles) runName(smiles, { fromLive: false });
      }
    });
  });
}

export function bindNamer() {
  bindEngineSwitch();
  $("namer-form") && $("namer-form").addEventListener("submit", onName);
  $("btn-pubchem") &&
    $("btn-pubchem").addEventListener("click", onPubchem);
  $("btn-clear-history") &&
    $("btn-clear-history").addEventListener("click", () => {
      state.namerHistory = [];
      renderNamerHistory();
    });
  $("btn-load-smiles") &&
    $("btn-load-smiles").addEventListener("click", async () => {
      const typed = (($("smiles-input") && $("smiles-input").value) || "").trim();
      if (!typed) {
        setNamerError("请输入 SMILES 再载入画板");
        return;
      }
      if (!state.ketcherBridge || !state.ketcherBridge.isReady()) {
        setNamerError("Ketcher 未就绪");
        return;
      }
      try {
        setNamerError("");
        // 输入可能为 merged_benchmark id → 先解析成该行真实 SMILES 再载入画板
        const smiles = await resolveTextToSmiles(typed);
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
