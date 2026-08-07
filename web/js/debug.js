/* Pipeline Debug page: run debug endpoint and render per-layer breakdown. */
import { $, escapeHtml, api, API } from "./core.js";

var debugDebounceTimer = null;
var DEBUG_DEBOUNCE_MS = 400;

export function scheduleLiveDebug() {
  if (debugDebounceTimer) clearTimeout(debugDebounceTimer);
  debugDebounceTimer = setTimeout(function () {
    debugDebounceTimer = null;
    runDebug();
  }, DEBUG_DEBOUNCE_MS);
}

export async function runDebug() {
  var input = document.getElementById("debug-smiles");
  if (!input) return;
  var smiles = input.value.trim();
  if (!smiles) return;

  var btn = document.getElementById("debug-run");
  if (btn) btn.disabled = true;

  try {
    var res = await api(API.debug, { method: "POST", body: JSON.stringify({ smiles: smiles }) });
    renderDebug(res);
  } catch (err) {
    var out = document.getElementById("debug-output");
    if (out) out.innerHTML = '<div class="debug-empty"><h2>Error</h2><p>' + escapeHtml(String(err)) + '</p></div>';
  } finally {
    if (btn) btn.disabled = false;
  }
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
    html += '<div class="layer-card' + (hasError ? ' error' : '') + '" data-layer="' + ld2.key + '">';
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
  return '<dl class="debug-kv">' +
    '<dt>SMILES</dt><dd>' + escapeHtml(l.smiles || "") + '</dd>' +
    '<dt>Atoms</dt><dd>' + l.num_atoms + '</dd>' +
    '<dt>Bonds</dt><dd>' + l.num_bonds + '</dd>' +
    '<dt>Composition</dt><dd>' + escapeHtml(l.composition || "") + '</dd>' +
    '</dl>';
}

function renderDebugL1(l) {
  var html = '<dl class="debug-kv">';
  html += '<dt>n_carbons</dt><dd>' + l.n_carbons + '</dd>';
  html += '<dt>n_ring_systems</dt><dd>' + l.n_ring_systems + '</dd>';
  html += '<dt>n_rings</dt><dd>' + l.n_rings + '</dd>';
  html += '<dt>molecule</dt><dd>' + (l.mol ? (l.mol.num_atoms + ' atoms, ' + l.mol.num_bonds + ' bonds') : '—') + '</dd>';
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

  // FG flags
  html += '<div class="debug-sub"><div class="debug-sub-title">Functional Groups</div>';
  var fgMap = [
    ['has_acid','COOH'],['has_alcohol','OH'],['has_alkene','C=C'],['has_alkyne','C≡C'],
    ['has_amide','CON'],['has_amine','NH2/NH'],['has_anhydride','(CO)2O'],['has_boronic','B(OH)2'],
    ['has_carbamate','OCON'],['has_carbonate','OCOO'],['has_ester','COOR'],['has_ether','C-O-C'],
    ['has_guanidine','N-C(=N)N'],['has_hydrazine','N-N'],['has_isocyanate','NCO'],
    ['has_isothiocyanate','NCS'],['has_ketone','C=O'],['has_nitrile','C≡N'],['has_nitro','NO2'],
    ['has_phosphate','OPO3'],['has_sulfide','C-S-C'],['has_sulfonamide','SO2N'],
    ['has_sulfonate','SO3R'],['has_sulfone','SO2'],['has_sulfonic_acid','SO3H'],
    ['has_sulfonyl_chloride','SO2Cl'],['has_sulfoxide','SO'],['has_thiol','SH'],['has_urea','NCON'],
  ];
  html += '<div class="summary-row">';
  for (var k = 0; k < fgMap.length; k++) {
    var f = fgMap[k];
    if (l[f[0]] === true) html += '<span class="chip accent">' + f[1] + '</span>';
  }
  html += '</div>';
  html += '<span class="small muted-text">hydroxyls=' + (l.hydroxyls||[]).length +
    ' amines=' + (l.amines||[]).length +
    ' carboxyls=' + (l.carboxyls||[]).length +
    ' esters=' + (l.esters||[]).length +
    ' ethers=' + (l.ethers||[]).length +
    ' ketones=' + (l.ketones||[]).length +
    ' amides=' + (l.amides||[]).length +
    ' dbl_bonds=' + (l.double_bonds||[]).length + '</span>';
  html += '</div>';

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
    debugSmiles.addEventListener("input", scheduleLiveDebug);
  }
}
