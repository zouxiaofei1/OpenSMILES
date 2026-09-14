/* Debug page: Pipeline 分层（/api/v1/debug）与 打印调试（/api/v1/debug-print）
   合并为单一页面，通过顶部开关切换模式，共享同一个 SMILES 输入框；批量打印模式
   忽略输入框，直接跑整个数据集。 */
import { $, escapeHtml, api, API } from "./core.js";

var debugDebounceTimer = null;
var DEBUG_DEBOUNCE_MS = 400;
var debugMode = "pipeline"; // "pipeline" | "print" | "batch"
var streamCtrl = null;      // 打印/批量流的 AbortController(两种模式互斥, 共用一个)
var streamSeq = 0;          // 每次 run 递增，用序号丢弃过期回调

function setDebugMode(mode) {
  debugMode = mode;
  stopStream(); // 离开当前模式即中止进行中的流
  var segBtns = document.querySelectorAll(".debug-mode-btn");
  for (var i = 0; i < segBtns.length; i++) {
    segBtns[i].classList.toggle("active", segBtns[i].getAttribute("data-debug-mode") === mode);
  }
  var pipeline = document.getElementById("debug-output");
  var printPanel = document.getElementById("print-debug-panel");
  var batchPanel = document.getElementById("batch-print-panel");
  if (pipeline) pipeline.hidden = mode !== "pipeline";
  if (printPanel) printPanel.hidden = mode !== "print";
  if (batchPanel) batchPanel.hidden = mode !== "batch";

  // 批量模式跑的是整个数据集, SMILES 输入无意义: 禁用并改为提示文案
  var input = document.getElementById("debug-smiles");
  if (input) {
    input.disabled = mode === "batch";
    input.placeholder = mode === "batch"
      ? "批量模式：不需要 SMILES，点 Run 跑整个 merged_benchmark"
      : "Enter SMILES, e.g. c1ccccc1C(=O)O or CC(=O)OC1CCCCC1";
  }

  // 批量模式不能一进来就自动开跑(整表要几十秒), 必须点 Run; 其余模式沿用实时执行
  if (mode !== "batch" && input && input.value.trim()) runDebug();
}

export function scheduleLiveDebug() {
  // 输入/更改停顿 DEBUG_DEBOUNCE_MS 后自动执行，无需点 Run；批量模式除外。
  // pipeline 走 /debug；print 走流式，新一轮 runPrintStream 会自动中止旧子进程。
  if (debugMode === "batch") return;
  if (debugDebounceTimer) clearTimeout(debugDebounceTimer);
  debugDebounceTimer = setTimeout(function () {
    debugDebounceTimer = null;
    runDebug();
  }, DEBUG_DEBOUNCE_MS);
}

export async function runDebug() {
  if (debugMode === "batch") {
    await runBatchStream();
    return;
  }

  var input = document.getElementById("debug-smiles");
  if (!input) return;
  var smiles = input.value.trim();
  if (!smiles) return;

  if (debugMode === "print") {
    await runPrintStream(smiles);
    return;
  }

  var btn = document.getElementById("debug-run");
  if (btn) btn.disabled = true;

  try {
    var res = await api(API.debug, { method: "POST", body: JSON.stringify({ smiles: smiles }) });
    renderDebug(res);
  } catch (err) {
    var out = document.getElementById("debug-output");
    if (out) {
      out.innerHTML = '<div class="debug-empty"><h2>Error</h2><p>' + escapeHtml(String(err)) + '</p></div>';
    }
  } finally {
    if (btn) btn.disabled = false;
  }
}

/* 打印调试 · 流式：读取 text/plain 响应，子进程 stdout 每分块一到就交给 onChunk
   增量更新显示，内容一变化即“实时”呈现，而不是等整次跑完。
   onChunk(text, err)：正常分块给 text；出错给 err（此时已到流尾）。 */
function stopStream() {
  if (streamCtrl) {
    try { streamCtrl.abort(); } catch (_) { /* ignore */ }
    streamCtrl = null;
  }
}

async function streamText(url, payload, onChunk) {
  stopStream(); // 新一轮 run 先取消仍在跑的旧流
  var mySeq = ++streamSeq;
  var ctrl = new AbortController();
  streamCtrl = ctrl;
  try {
    var res = await fetch(url, {
      method: "POST",
      headers: { "Content-Type": "application/json", Accept: "text/plain" },
      body: JSON.stringify(payload),
      signal: ctrl.signal,
    });
    if (!res.ok || !res.body) {
      var errText = "";
      try { errText = await res.text(); } catch (_) { /* ignore */ }
      onChunk("", "Error: " + (errText || (res.status + " " + res.statusText)));
      return;
    }

    var reader = res.body.getReader();
    var dec = new TextDecoder("utf-8");
    for (;;) {
      var chunk = await reader.read();
      if (chunk.done) break;
      if (mySeq !== streamSeq) { // 已被新一轮/切模式取消，丢弃
        try { reader.cancel(); } catch (_) { /* ignore */ }
        return;
      }
      onChunk(dec.decode(chunk.value, { stream: true }), "");
    }
    onChunk(dec.decode(), "");
  } catch (err) {
    if (mySeq !== streamSeq) return; // 被取消，正常路径已处理
    if (err && err.name === "AbortError") return;
    onChunk("", "Error: " + (err && err.message ? err.message : String(err)));
  } finally {
    if (mySeq === streamSeq) streamCtrl = null;
  }
}

async function runPrintStream(smiles) {
  var out = document.getElementById("print-debug-output");
  if (!out) return;
  out.classList.remove("has-error");
  out.textContent = "";

  var text = "";
  await streamText(API.debugPrintStream, { smiles: smiles }, function (chunk, err) {
    if (err) {
      out.classList.add("has-error");
      out.textContent = text + err;
      return;
    }
    text += chunk;
    out.textContent = text;
    // 原本贴底则跟随滚动；用户上翻查看时不被强制拉底
    if (out.scrollHeight - out.scrollTop - out.clientHeight < 24) {
      out.scrollTop = out.scrollHeight;
    }
  });
  if (!text && !out.classList.contains("has-error")) out.textContent = "(无输出)";
}

/* ── 批量打印模式：整表跑一遍，输出按分子分组 ──
   batchText 是累积的原始全文；无过滤时增量追加(便宜)，有过滤时按行筛选后整段重画
   (节流)，命中行前补上所属分子的分组头，才知道是哪个分子触发的该 print。 */
var batchText = "";
var batchPaintedLen = 0;     // 已画出的字符数（无过滤时用于增量追加）
var batchPaintedFilter = ""; // 上次画图用的过滤词，变了就整段重画
var batchLastPaint = 0;
var batchHits = 0;
var batchRunning = false;
var batchDone = false;
var BATCH_PAINT_MS = 150;

function batchFilterValue() {
  var el = document.getElementById("batch-filter");
  return el ? el.value.trim().toLowerCase() : "";
}

/* 分组头形如 "===== [12/3976] CCC(=O)N ====="，用于给命中行补出分子上下文。 */
function isBatchBanner(line) {
  return /^\s*=+ \[\d+\/\d+\]/.test(line);
}

/* 按行过滤批量输出：命中行连同它所属分子的分组头一起留下。 */
function filterBatchText(text, filter) {
  var lines = text.split("\n");
  var kept = [];
  var banner = "";
  var hits = 0;
  for (var i = 0; i < lines.length; i++) {
    var line = lines[i];
    if (isBatchBanner(line)) { banner = line; continue; }
    if (line.toLowerCase().indexOf(filter) !== -1) {
      if (banner) { kept.push(banner); banner = ""; }
      kept.push(line);
      hits++;
    }
  }
  return { text: kept.length ? kept.join("\n") : "（没有匹配「" + filter + "」的行）", hits: hits };
}

/* 取最近的进度（后端每 250 个分子打一行 "##### 进度 i/N"）。分组头是稀疏的，不能
   当进度用；只扫尾部避免全文正则。 */
function batchProgress(text) {
  var m = text.slice(-20000).match(/进度 (\d+)\/(\d+)/);
  return m ? m[1] + "/" + m[2] : "";
}

function updateBatchStat() {
  var el = document.getElementById("batch-filter-stat");
  if (!el) return;
  var parts = [];
  if (batchRunning) parts.push("运行中");
  else if (batchDone) parts.push("已完成");
  var progress = batchProgress(batchText);
  if (progress) parts.push(progress);
  if (batchFilterValue()) parts.push("命中 " + batchHits + " 行");
  el.textContent = parts.join(" · ");
}

/* 重画批量输出（force=true 立刻画，不受节流限制）。 */
function paintBatch(force) {
  var out = document.getElementById("batch-print-output");
  if (!out || !batchText) return;
  var filter = batchFilterValue();
  var now = Date.now();
  if (filter) {
    if (!force && filter === batchPaintedFilter && now - batchLastPaint < BATCH_PAINT_MS) return;
    var r = filterBatchText(batchText, filter);
    out.textContent = r.text;
    batchHits = r.hits;
  } else if (force || batchPaintedFilter || batchPaintedLen > batchText.length) {
    out.textContent = batchText;
  } else {
    out.textContent += batchText.slice(batchPaintedLen);
  }
  batchPaintedLen = batchText.length;
  batchPaintedFilter = filter;
  batchLastPaint = now;
  updateBatchStat();
}

async function runBatchStream() {
  var out = document.getElementById("batch-print-output");
  if (!out) return;
  out.classList.remove("has-error");
  batchText = "";
  batchPaintedLen = 0;
  batchPaintedFilter = "";
  batchHits = 0;
  batchRunning = true;
  batchDone = false;
  out.textContent = "";
  updateBatchStat();

  await streamText(API.debugPrintBatchStream, {}, function (chunk, err) {
    if (err) {
      out.classList.add("has-error");
      batchText += "\n" + err + "\n";
    } else {
      batchText += chunk;
    }
    if (!batchText) return;
    paintBatch(false);
    if (out.scrollHeight - out.scrollTop - out.clientHeight < 24) {
      out.scrollTop = out.scrollHeight;
    }
  });

  batchRunning = false;
  batchDone = true;
  if (!batchText && !out.classList.contains("has-error")) out.textContent = "(无输出)";
  paintBatch(true);
}

function renderDebug(data) {
  var out = document.getElementById("debug-output");
  if (!out) return;
  if (!data.success) {
    out.innerHTML = '<div class="debug-empty"><h2>Error</h2><p>' + escapeHtml(data.error || "Unknown error") + '</p></div>';
    return;
  }

  var layers = data.layers;
  var layerDefs = [
    { key: "L0_preprocess",   icon: "l0", name: "Layer 0 — Preprocess",     desc: "SMILES → Mol (RDKit sanitization, kekulization)" },
    { key: "L1_analyze",      icon: "l1", name: "Layer 1 — Analyze",         desc: "Mol → Info dict (FGs, ring systems, unsaturation)" },
    { key: "L2_parent",       icon: "l2", name: "Layer 2 — Parent Select",   desc: "Info → Parent candidates → Selected parent hydride" },
    { key: "L3_substituents", icon: "l3", name: "Layer 3 — Substituents",    desc: "Parent + Info → Substituent extraction + Coverage ledger" },
    { key: "L4_numbering",    icon: "l4", name: "Layer 4 — Numbering",       desc: "Parent + Substituents → Numbered dict (locants, chain orient)" },
    { key: "L5_assemble",     icon: "l5", name: "Layer 5 — Assembly",        desc: "Numbered dict → NameResult (en + zh + stereo)" },
  ];

  var html = '<div class="summary-row">';
  html += '<span class="chip info">Total: <b>' + data.total_time_ms + ' ms</b></span>';
  for (var i = 0; i < layerDefs.length; i++) {
    var ld = layerDefs[i];
    var l = layers[ld.key];
    if (l) {
      html += '<span class="chip">' + ld.key.replace("_"," ").replace("_"," ") + ': <b>' + (l.time_ms || "?") + ' ms</b></span>';
    }
  }
  html += '</div>';

  html += '<div class="layer-list">';
  for (var j = 0; j < layerDefs.length; j++) {
    var ld2 = layerDefs[j];
    var l2 = layers[ld2.key];
    if (!l2) continue;
    var hasError = !!l2.error;
    // 默认 open: 分层始终展开，用户无需每次手动点开（仍可点标题临时收起）
    html += '<div class="layer-card open' + (hasError ? ' error' : '') + '" data-layer="' + ld2.key + '">';
    html += '<div class="layer-head" onclick="ChemNamerDebug.toggleLayer(this)">';
    html += '<div class="layer-icon ' + ld2.icon + '">' + ld2.key[1] + '</div>';
    html += '<div class="layer-info"><div class="layer-name">' + ld2.name + '</div>';
    html += '<div class="layer-desc">' + ld2.desc + '</div></div>';
    html += '<div class="layer-meta">';
    if (hasError) {
      html += '<span class="badge fail">ERROR</span>';
    }
    html += '<span class="layer-time">' + (l2.time_ms || "?") + ' ms</span>';
    html += '<svg class="layer-chevron" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><polyline points="6 9 12 15 18 9"/></svg>';
    html += '</div></div>';
    html += '<div class="layer-content">' + renderDebugLayer(ld2.key, l2) + '</div>';
    html += '</div>';
  }
  html += '</div>';

  out.innerHTML = html;
}

function renderDebugLayer(key, l) {
  if (l.error) {
    return '<div style="color:var(--color-destructive-soft);font-family:var(--mono);">' + escapeHtml(l.error) + '</div>';
  }
  switch (key) {
    case "L0_preprocess":   return renderDebugL0(l);
    case "L1_analyze":      return renderDebugL1(l);
    case "L2_parent":       return renderDebugL2(l);
    case "L3_substituents": return renderDebugL3(l);
    case "L4_numbering":    return renderDebugL4(l);
    case "L5_assemble":     return renderDebugL5(l);
    default:                return '<pre class="json-block">' + JSON.stringify(l, null, 2) + '</pre>';
  }
}

function renderDebugL0(l) {
  var html = '<dl class="debug-kv">' +
    '<dt>SMILES</dt><dd>' + escapeHtml(l.smiles || "") + '</dd>' +
    '<dt>Atoms</dt><dd>' + l.num_atoms + '</dd>' +
    '<dt>Bonds</dt><dd>' + l.num_bonds + '</dd>' +
    '<dt>Composition</dt><dd>' + escapeHtml(l.composition || "") + '</dd>';
  if (l.atoms || l.bonds) {
    html += '<dt>raw graph</dt><dd><details class="mol-details"><summary>show atoms / bonds</summary><pre class="json-block">' +
      escapeHtml(JSON.stringify({ smiles: l.smiles, atoms: l.atoms, bonds: l.bonds }, null, 2)) + '</pre></details></dd>';
  }
  html += '</dl>';
  return html;
}

function renderDebugL1(l) {
  var html = '<dl class="debug-kv">';
  html += '<dt>n_carbons</dt><dd>' + l.n_carbons + '</dd>';
  html += '<dt>n_ring_systems</dt><dd>' + l.n_ring_systems + '</dd>';
  html += '<dt>n_rings</dt><dd>' + l.n_rings + '</dd>';
  html += '<dt>molecule</dt><dd>' + (l.mol ? (l.mol.num_atoms + ' atoms, ' + l.mol.num_bonds + ' bonds' +
    '<details class="mol-details"><summary>show raw Mol graph</summary><pre class="json-block">' +
    escapeHtml(JSON.stringify(l.mol, null, 2)) + '</pre></details>') : '—') + '</dd>';
  html += '</dl>';

  // Ring systems
  if (l.ring_systems && l.ring_systems.length) {
    html += '<div class="debug-sub"><div class="debug-sub-title">Ring Systems (' + l.ring_systems.length + ')</div>';
    for (var i = 0; i < l.ring_systems.length; i++) {
      var rs = l.ring_systems[i];
      var hets = (rs.hetero_atoms || []).map(function (h) { return 'idx=' + h.idx + '(Z=' + h.Z + ')'; }).join(', ');
      html += '<div class="candidate-card">';
      html += '<span class="badge">' + rs.topology + '</span> ';
      html += '<span class="chip">' + rs.n_rings + ' rings</span> ';
      html += '<span class="chip">' + rs.n_atoms + ' atoms</span> ';
      html += '<span class="chip ' + (rs.is_aromatic_mancude ? 'accent' : '') + '">mancude=' + rs.is_aromatic_mancude + '</span>';
      if (hets) html += '<span class="chip info">hetero: ' + escapeHtml(hets) + '</span>';
      html += '<div style="margin-top:4px"><span class="atom-list">atoms=[' + (rs.atom_ids || []).slice(0, 30).join(',') + (rs.atom_ids && rs.atom_ids.length > 30 ? '...' : '') + ']</span></div>';
      html += '</div>';
    }
    html += '</div>';
  }

  // Rings
  if (l.rings && l.rings.length) {
    html += '<div class="debug-sub"><div class="debug-sub-title">SSSR Rings (' + l.rings.length + ')</div>';
    html += '<table class="subst-table"><thead><tr><th>#</th><th>Size</th><th>Atoms</th></tr></thead><tbody>';
    for (var j = 0; j < l.rings.length; j++) {
      var r = l.rings[j];
      var ids = r.atom_ids || [];
      html += '<tr><td>' + j + '</td><td>' + ids.length + '-membered</td><td class="atom-list">[' + ids.join(',') + ']</td></tr>';
    }
    html += '</tbody></table></div>';
  }

  // Functional Groups: presence derived from every FG entry list (non-empty = present) + each list with its content
  html += '<div class="debug-sub"><div class="debug-sub-title">Functional Groups</div>';
  var fgMap = {
    carboxyls:'COOH', hydroxyls:'OH', aldehydes:'CHO', double_bonds:'C=C', triple_bonds:'C≡C',
    acyl_chlorides:'COCl', amides:'CON', amines:'NH2/NH', anhydrides:'(CO)2O', esters:'COOR',
    ketones:'C=O', nitriles:'C≡N', phosphates:'OPO3', thiols:'SH', acyls:'C(=O)-*',
    radicals:'*', demoted_carboxyls:'COOH↓', demoted_nitriles:'C≡N↓',
  };
  // FG entry lists (dynamic: any list key in the info dict that is not structural)
  var STRUCT_KEYS = { ring_systems: 1, rings: 1, carbon_ids: 1, mol: 1 };
  var fgListKeys = Object.keys(l).filter(function (k) {
    return Array.isArray(l[k]) && !STRUCT_KEYS[k];
  });
  html += '<div class="summary-row">';
  for (var i = 0; i < fgListKeys.length; i++) {
    var k = fgListKeys[i];
    if (l[k].length) {
      html += '<span class="chip accent" title="' + escapeHtml(k) + '">' + escapeHtml(fgMap[k] || k) + '</span>';
    }
  }
  html += '</div>';
  var off = fgListKeys.filter(function (k) { return !l[k].length; });
  if (off.length) {
    html += '<div class="fg-flags-off">' + off.map(function (k) {
      return '<span class="fg-flag-off" title="' + escapeHtml(k) + '">' + escapeHtml(fgMap[k] || k) + '</span>';
    }).join('') + '</div>';
  }

  var entryKeys = fgListKeys.filter(function (k) { return l[k].length; });
  if (entryKeys.length) {
    html += '<div class="fg-entries"><table class="subst-table"><thead><tr><th>FG</th><th>Count</th><th>Entries</th></tr></thead><tbody>';
    for (var j = 0; j < entryKeys.length; j++) {
      var ek = entryKeys[j];
      html += '<tr><td>' + escapeHtml(ek) + '</td><td>' + l[ek].length + '</td><td class="atom-list">' +
        escapeHtml(JSON.stringify(l[ek])) + '</td></tr>';
    }
    html += '</tbody></table></div>';
  }

  // Raw info dict fallback: everything else (carbon_ids, fg_inventory, ...)
  html += '<details class="mol-details"><summary>show raw L1 info dict (' + Object.keys(l).length + ' keys)</summary>' +
    '<pre class="json-block">' + escapeHtml(JSON.stringify(l, null, 2)) + '</pre></details>';

  return html;
}

function renderDebugL2(l) {
  var html = '<dl class="debug-kv">';
  html += '<dt>n_candidates</dt><dd>' + l.n_candidates + '</dd>';
  html += '<dt>selected kind</dt><dd style="color:var(--color-accent-soft);font-weight:600;">' + escapeHtml(l.selected && l.selected.kind || '—') + '</dd>';
  html += '</dl>';

  if (l.candidates && l.candidates.length) {
    html += '<div class="debug-sub"><div class="debug-sub-title">All Candidates</div>';
    for (var i = 0; i < l.candidates.length; i++) {
      var c = l.candidates[i];
      var sel = c.is_selected;
      html += '<div class="candidate-card' + (sel ? ' selected' : '') + '">';
      html += '<div class="candidate-kind">' + escapeHtml(c.kind || '?');
      if (sel) html += '<span class="sel-tag">SELECTED</span>';
      html += '</div>';
      html += '<pre class="json-block inline">' + JSON.stringify(c, null, 2) + '</pre>';
      html += '</div>';
    }
    html += '</div>';
  }

  return html;
}

function renderDebugL3(l) {
  var html = '<dl class="debug-kv">';
  html += '<dt>substituents</dt><dd>' + l.n_substituents + '</dd>';
  html += '<dt>coverage complete</dt><dd><span class="status-dot ' + (l.coverage && l.coverage.complete ? 'ok' : 'fail') + '"></span>' + (l.coverage && l.coverage.complete ? 'YES' : 'NO') + '</dd>';
  if (l.coverage && l.coverage.gap && l.coverage.gap.length) {
    html += '<dt>gap atoms</dt><dd class="atom-list">[' + l.coverage.gap.join(',') + ']</dd>';
  }
  if (l.coverage && l.coverage.overlap && l.coverage.overlap.length) {
    html += '<dt>overlap atoms</dt><dd class="atom-list">[' + l.coverage.overlap.join(',') + ']</dd>';
  }
  html += '<dt>owned atoms</dt><dd class="atom-list">[' + (l.owned_atoms||[]).slice(0,40).join(',') + ((l.owned_atoms||[]).length > 40 ? '...' : '') + ']</dd>';
  html += '</dl>';

  if (l.substituents && l.substituents.length) {
    html += '<div class="debug-sub"><div class="debug-sub-title">Substituent List</div>';
    html += '<table class="subst-table"><thead><tr><th>#</th><th>EN</th><th>ZH</th><th>Kind</th><th>Attach</th><th>Paren</th><th>Atoms</th></tr></thead><tbody>';
    for (var i = 0; i < l.substituents.length; i++) {
      var s = l.substituents[i];
      html += '<tr>' +
        '<td>' + i + '</td>' +
        '<td style="font-weight:500;">' + escapeHtml(s.en || '?') + '</td>' +
        '<td>' + escapeHtml(s.zh || '') + '</td>' +
        '<td>' + escapeHtml(s.kind || '') + '</td>' +
        '<td>' + (s.attach_idx !== undefined ? s.attach_idx : '—') + '</td>' +
        '<td>' + (s.paren ? 'yes' : '') + '</td>' +
        '<td class="atom-list">[' + (s.atoms||[]).join(',') + ']</td>' +
        '</tr>';
    }
    html += '</tbody></table></div>';
  }

  return html;
}

function renderDebugL4(l) {
  var html = '<dl class="debug-kv">';
  html += '<dt>parent kind</dt><dd style="color:var(--color-accent-soft);font-weight:600;">' + escapeHtml(l.parent_kind || '—') + '</dd>';
  html += '<dt>chain</dt><dd class="atom-list">[' + (l.chain||[]).join(',') + ']</dd>';
  html += '</dl>';

  if (l.locants && Object.keys(l.locants).length) {
    html += '<div class="debug-sub"><div class="debug-sub-title">Locants</div>';
    html += '<dl class="debug-kv">';
    var locKeys = Object.keys(l.locants);
    for (var i = 0; i < locKeys.length; i++) {
      var k = locKeys[i];
      html += '<dt>' + escapeHtml(k) + '</dt><dd>' + escapeHtml(String(l.locants[k])) + '</dd>';
    }
    html += '</dl></div>';
  }

  if (l.substituents && l.substituents.length) {
    html += '<div class="debug-sub"><div class="debug-sub-title">Numbered Substituents (' + l.substituents.length + ')</div>';
    html += '<table class="subst-table"><thead><tr><th>#</th><th>EN</th><th>Locant</th><th>Kind</th></tr></thead><tbody>';
    for (var j = 0; j < l.substituents.length; j++) {
      var s = l.substituents[j];
      html += '<tr>' +
        '<td>' + j + '</td>' +
        '<td style="font-weight:500;">' + escapeHtml(s.en || '?') + '</td>' +
        '<td>' + (s.locant !== undefined ? s.locant : '—') + '</td>' +
        '<td>' + escapeHtml(s.kind || '') + '</td>' +
        '</tr>';
    }
    html += '</tbody></table></div>';
  }

  return html;
}

function renderDebugL5(l) {
  var html = '';
  if (l.en || l.zh) {
    html += '<div style="display:flex;flex-direction:column;gap:var(--space-md);">';

    html += '<div style="background:var(--color-surface-2);border:1px solid rgba(34,197,94,0.3);border-radius:var(--radius);padding:var(--space-lg);">';
    html += '<div class="field-label">English Name</div>';
    html += '<div style="font-family:var(--mono);font-size:15px;font-weight:500;color:var(--color-foreground);word-break:break-all;">' + escapeHtml(l.en) + '</div>';
    html += '</div>';

    html += '<div style="background:var(--color-surface-2);border:1px solid rgba(56,189,248,0.3);border-radius:var(--radius);padding:var(--space-lg);">';
    html += '<div class="field-label">Chinese Name</div>';
    html += '<div style="font-family:var(--mono);font-size:15px;font-weight:500;color:var(--color-foreground);word-break:break-all;">' + escapeHtml(l.zh) + '</div>';
    html += '</div>';

    html += '<dl class="debug-kv">';
    html += '<dt>success</dt><dd>' + l.success + '</dd>';
    html += '<dt>source</dt><dd>' + escapeHtml(l.source || '') + '</dd>';
    html += '</dl>';

    html += '</div>';
  }
  if (l.meta && Object.keys(l.meta).length) {
    html += '<div class="debug-sub"><div class="debug-sub-title">Meta</div>';
    html += '<pre class="json-block inline">' + JSON.stringify(l.meta, null, 2) + '</pre></div>';
  }
  return html || '<span class="muted-text">(no output)</span>';
}

export function bindDebug() {
  // Debug page controls
  var debugRun = document.getElementById("debug-run");
  if (debugRun) debugRun.addEventListener("click", runDebug);
  var debugSmiles = document.getElementById("debug-smiles");
  if (debugSmiles) {
    debugSmiles.addEventListener("keydown", function (e) {
      if (e.key === "Enter") runDebug();
    });
    debugSmiles.addEventListener("input", scheduleLiveDebug);
  }
  // 模式开关：Pipeline 分层 / 打印输出 / 批量打印
  var modeBtns = document.querySelectorAll(".debug-mode-btn");
  for (var i = 0; i < modeBtns.length; i++) {
    modeBtns[i].addEventListener("click", function () {
      setDebugMode(this.getAttribute("data-debug-mode"));
    });
  }
  // 批量模式过滤框：只重画显示，原始输出仍留在 batchText 里，清空即恢复全文
  var batchFilter = document.getElementById("batch-filter");
  if (batchFilter) {
    batchFilter.addEventListener("input", function () { paintBatch(true); });
  }
}
