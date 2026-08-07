/* Wiki page: directory tree + markdown document viewer. */
import { $, state, escapeHtml, api, API } from "./core.js";

function renderWikiTree(entry, depth) {
  var html = "";
  if (entry.name !== "wiki") {
    html +=
      '<details class="wiki-dir" data-path="' + escapeHtml(entry.path) + '"' + (depth === 0 ? " open" : "") + ">" +
      '<summary style="padding-left:' + (depth * 14 + 4) + 'px">' + escapeHtml(entry.name) + "</summary>";
  }
  entry.files.forEach(function (f) {
    html += '<div class="wiki-file" data-path="' + escapeHtml(f.path) + '" style="padding-left:' + ((depth + 1) * 14 + 8) + 'px">' + escapeHtml(f.name) + "</div>";
  });
  entry.dirs.forEach(function (d) {
    html += renderWikiTree(d, depth + 1);
  });
  if (entry.name !== "wiki") html += "</details>";
  return html;
}

export async function loadWiki() {
  try {
    var data = await api(API.wiki);
    if (!data || !data.ok) throw new Error((data && data.error) || "加载失败");
    var treeEl = $("wiki-tree");
    if (treeEl) treeEl.innerHTML = renderWikiTree(data.tree, 0);
    loadWikiDoc("index.md");
  } catch (err) {
    var doc = $("wiki-doc");
    if (doc) doc.innerHTML = '<p class="muted-text">加载失败: ' + escapeHtml(err.message || err) + "</p>";
  }
}

async function loadWikiDoc(path) {
  // 解析相对路径：含 / 按 wiki 根相对；纯文件名按「当前文档同目录」优先、根兜底
  var cur = state.wikiCurrent || "";
  var curDir = cur.indexOf("/") !== -1 ? cur.substring(0, cur.lastIndexOf("/")) : "";
  var candidates = path.indexOf("/") !== -1
    ? [path]
    : (curDir ? [curDir + "/" + path, path] : [path]);

  var doc = $("wiki-doc");
  for (var i = 0; i < candidates.length; i++) {
    var data = null;
    try {
      data = await api(API.wikiDoc + "?path=" + encodeURIComponent(candidates[i]));
    } catch (e) {
      data = null;
    }
    if (data && data.ok) {
      state.wikiCurrent = data.path;
      var body = typeof marked !== "undefined" && marked.parse
        ? marked.parse(data.content)
        : escapeHtml(data.content);
      if (doc) doc.innerHTML = '<div class="wiki-title">' + escapeHtml(data.name) + "</div>" + body;
      document.querySelectorAll(".wiki-file.active").forEach(function (el) {
        el.classList.remove("active");
      });
      var sel = document.querySelector('.wiki-file[data-path="' + CSS.escape(data.path) + '"]');
      if (sel) sel.classList.add("active");
      return;
    }
  }
  if (doc) doc.innerHTML = '<p class="muted-text">无法打开文档: ' + escapeHtml(path) + "</p>";
}

export function bindWiki() {
  var tree = $("wiki-tree");
  if (tree) {
    tree.addEventListener("click", function (ev) {
      var f = ev.target.closest(".wiki-file");
      if (f) loadWikiDoc(f.getAttribute("data-path"));
    });
  }
  var doc = $("wiki-doc");
  if (doc) {
    doc.addEventListener("click", function (ev) {
      var a = ev.target.closest("a");
      if (!a) return;
      var href = (a.getAttribute("href") || "").trim();
      // 外部链接 / 锚点保留默认行为
      if (!href || /^[a-z]+:/i.test(href) || href.charAt(0) === "#") return;
      // 仅拦截指向 wiki 内 .md 的相对链接
      if (!/\.md(?:[#?]|$)/i.test(href)) return;
      ev.preventDefault();
      loadWikiDoc(href.split(/[#?]/)[0]);
    });
  }
}
