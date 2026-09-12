/* Code Analysis page: per-layer code-line stats (donut + expandable file list). */
import { $, state, escapeHtml, api, API } from "./core.js";

var CA_COLORS = ["#22c55e", "#38bdf8", "#f59e0b", "#a78bfa", "#f472b6", "#2dd4bf", "#f97316", "#94a3b8"];  // [6]=tools, [7]=核心包

// 深浅变体：amt>0 向白混合（变浅），amt<0 向黑混合（变深）
function caShade(hex, amt) {
  var c = parseInt(hex.slice(1), 16);
  var r = (c >> 16) & 255, g = (c >> 8) & 255, b = c & 255;
  if (amt >= 0) {
    r = Math.round(r + (255 - r) * amt);
    g = Math.round(g + (255 - g) * amt);
    b = Math.round(b + (255 - b) * amt);
  } else {
    amt = -amt;
    r = Math.round(r * (1 - amt));
    g = Math.round(g * (1 - amt));
    b = Math.round(b * (1 - amt));
  }
  return "#" + ((1 << 24) + (r << 16) + (g << 8) + b).toString(16).slice(1);
}

function ca$(id) {
  return document.getElementById(id);
}

function caSkeleton() {
  return (
    '<div class="ca-skeleton" aria-hidden="true">' +
    '<div class="skeleton" style="width:180px;height:180px;border-radius:50%;margin:0 auto"></div>' +
    "</div>" +
    '<div class="ca-list">' +
    Array.from({ length: 6 })
      .map(function () {
        return (
          '<div class="ca-row"><div class="skeleton" style="width:12px;height:12px;border-radius:50%"></div>' +
          '<div class="skeleton" style="width:80px;height:14px"></div>' +
          '<div class="skeleton" style="width:60px;height:14px"></div>' +
          '<div class="skeleton" style="flex:1;height:8px;border-radius:4px"></div></div>'
        );
      })
      .join("") +
    "</div>"
  );
}

export async function loadCodeAnalysis() {
  var donut = ca$("ca-donut");
  var list = ca$("ca-list");
  var summary = ca$("ca-summary");
  if (donut) donut.innerHTML = caSkeleton();
  if (list) list.innerHTML = "";
  if (summary) { summary.classList.remove("hidden"); summary.textContent = "统计中…"; }
  try {
    var data = await api(API.codeAnalysis);
    if (!data || !data.ok || !data.total) throw new Error("无数据");
    state.caData = data;
    if (summary) renderCaSummary(summary, data.total);
    if (donut) donut.innerHTML = renderCaDonut(data.layers, data.total.code);
    if (list) list.innerHTML = renderCaList(data.layers);
  } catch (err) {
    if (summary) {
      summary.classList.remove("hidden");
      summary.textContent = "加载失败: " + (err.message || String(err));
    }
    if (donut) donut.innerHTML = '<p class="muted-text">暂无统计数据</p>';
  }
}

function renderCaSummary(el, total) {
  el.innerHTML =
    '<div class="ca-stat"><span class="ca-stat-label">总文件</span><b class="mono">' + (total.file_count != null ? total.file_count : total.files) + "</b></div>" +
    '<div class="ca-stat"><span class="ca-stat-label">代码行</span><b class="mono">' + total.code.toLocaleString() + "</b></div>" +
    '<div class="ca-stat"><span class="ca-stat-label">注释行</span><b class="mono">' + total.comment.toLocaleString() + "</b></div>" +
    '<div class="ca-stat"><span class="ca-stat-label">总行数</span><b class="mono">' + total.lines.toLocaleString() + "</b></div>";
}

function renderCaDonut(layers, totalCode) {
  var size = 220, cx = 110, cy = 110, r = 78, w = 26;
  var circ = 2 * Math.PI * r;
  // 最大扇区从 12 点起，顺时针按占比降序
  var sorted = layers.slice().sort(function (a, b) { return b.code - a.code; });
  var acc = 0;
  var segs = sorted.map(function (l) {
    var frac = totalCode > 0 ? l.code / totalCode : 0;
    var len = frac * circ;
    var seg =
      '<circle cx="' + cx + '" cy="' + cy + '" r="' + r + '" fill="none" ' +
      'stroke="' + CA_COLORS[l.layer] + '" stroke-width="' + w + '" ' +
      'stroke-dasharray="' + len.toFixed(2) + ' ' + (circ - len).toFixed(2) + '" ' +
      'stroke-dashoffset="' + (-acc).toFixed(2) + '" ' +
      'transform="rotate(-90 ' + cx + ' ' + cy + ')" class="ca-seg">' +
      "<title>" + (l.label || ("Layer " + l.layer)) + ": " + l.code + " 行 (" + l.share + "%)</title></circle>";
    acc += len;
    return seg;
  });
  return (
    '<svg viewBox="0 0 ' + size + " " + size + '" class="ca-donut" role="img" aria-label="各 Layer 代码行占比">' +
    segs.join("") +
    '<text x="' + cx + '" y="' + (cy - 4) + '" class="ca-donut-total" text-anchor="middle">' + totalCode.toLocaleString() + "</text>" +
    '<text x="' + cx + '" y="' + (cy + 18) + '" class="ca-donut-cap" text-anchor="middle">code lines</text>' +
    "</svg>"
  );
}

function renderCaList(layers) {
  var maxCode = layers.reduce(function (m, l) { return Math.max(m, l.code); }, 0);
  return (
    '<div class="ca-list-head">' +
    '<span></span><span>Layer</span><span class="mono">文件</span><span class="mono">代码行</span><span>占比</span><span class="mono">%</span><span></span>' +
    "</div>" +
    layers
      .map(function (l) {
        var bw = maxCode > 0 ? Math.round(100 * l.code / maxCode) : 0;
        return (
          '<div class="ca-row ca-row-toggle" data-layer="' + l.layer + '" role="button" tabindex="0" aria-expanded="false">' +
          '<span class="ca-dot" style="background:' + CA_COLORS[l.layer] + '"></span>' +
          '<span class="ca-name">' + (l.label || ("Layer " + l.layer)) + "</span>" +
          '<span class="ca-files mono">' + l.file_count + "</span>" +
          '<span class="ca-lines mono">' + l.code.toLocaleString() + "</span>" +
          '<span class="ca-bar-wrap"><span class="ca-bar" style="width:' + bw + "%;background:" + CA_COLORS[l.layer] + '"></span></span>' +
          '<span class="ca-share mono">' + l.share + "%</span>" +
          '<span class="ca-chevron" aria-hidden="true"></span>' +
          "</div>" +
          '<div class="ca-file-sub hidden" data-file-layer="' + l.layer + '"></div>'
        );
      })
      .join("")
  );
}

function caFileRows(layer) {
  var l = (state.caData && state.caData.layers || []).find(function (x) { return x.layer === layer; });
  if (!l || !l.files || !l.files.length) return '<p class="muted-text small">无文件</p>';
  var color = CA_COLORS[layer];
  var shade = caShade(color, 0.6); // 浅色 = 总行数（含注释/空行）
  var maxTotal = l.files.reduce(function (m, f) { return Math.max(m, f.code + f.comment + f.blank); }, 0);
  return (
    '<div class="ca-file-row ca-file-head">' +
    '<span class="ca-file-name">文件名</span><span class="mono">代码行</span><span>分布</span><span class="mono">总行数</span>' +
    "</div>" +
    l.files
      .map(function (f) {
        var total = f.code + f.comment + f.blank;
        var codeW = maxTotal > 0 ? Math.round(100 * f.code / maxTotal) : 0;
        var otherW = maxTotal > 0 ? Math.round(100 * (total - f.code) / maxTotal) : 0;
        return (
          '<div class="ca-file-row">' +
          '<span class="ca-file-name mono">' + escapeHtml(f.name) + "</span>" +
          '<span class="ca-file-lines mono" style="color:' + color + '">' + f.code.toLocaleString() + "</span>" +
          '<span class="ca-file-bar-wrap">' +
          '<span class="ca-file-bar" style="width:' + codeW + "%;background:" + color + '"></span>' +
          '<span class="ca-file-bar ca-file-bar-shade" style="width:' + otherW + "%;background:" + shade + '"></span>' +
          "</span>" +
          '<span class="ca-file-total mono" style="color:' + shade + '">' + total.toLocaleString() + "</span>" +
          "</div>"
        );
      })
      .join("")
  );
}

function toggleCaRow(row) {
  var layer = parseInt(row.getAttribute("data-layer"), 10);
  var expanded = row.getAttribute("aria-expanded") === "true";
  row.setAttribute("aria-expanded", String(!expanded));
  row.classList.toggle("open", !expanded);
  var sub = document.querySelector('.ca-file-sub[data-file-layer="' + layer + '"]');
  if (sub) {
    if (!expanded && !sub.getAttribute("data-rendered")) {
      sub.innerHTML = caFileRows(layer);
      sub.setAttribute("data-rendered", "1");
    }
    sub.classList.toggle("hidden", expanded);
  }
}

function bindCaToggle() {
  var list = ca$("ca-list");
  if (!list) return;
  list.addEventListener("click", function (ev) {
    var row = ev.target.closest(".ca-row-toggle");
    if (row) toggleCaRow(row);
  });
  list.addEventListener("keydown", function (ev) {
    if (ev.key !== "Enter" && ev.key !== " ") return;
    var row = ev.target.closest(".ca-row-toggle");
    if (row) {
      ev.preventDefault();
      toggleCaRow(row);
    }
  });
}

function bindCaRefresh() {
  var btn = $("ca-refresh");
  if (!btn) return;
  btn.addEventListener("click", async function () {
    var st = $("ca-refresh-status");
    if (st) st.textContent = "刷新中…";
    await loadCodeAnalysis();
    if (st) st.textContent = "已刷新 " + new Date().toLocaleTimeString();
  });
}

export function bindCodeAnalysis() {
  bindCaToggle();
  bindCaRefresh();
}
