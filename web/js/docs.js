/* Docs page: docs/wiki 目录树 + IUPAC cn_translated 条目 共用查看器。
   侧边栏「Wiki / IUPAC」来源切换；加载与渲染分别搬自原 wiki.js / iupac.js。 */
import { $, state, escapeHtml, api, API } from "./core.js";

/* 来源元数据：各来源的列表/正文接口、默认文档、正文标题字段 */
const SOURCES = {
  wiki: {
    listURL: API.wiki,
    docURL: API.wikiDoc,
    defaultDoc: () => state.wikiCurrent || "index.md",
    setCurrent: (path) => { state.wikiCurrent = path; },
    titleField: "name", // wiki doc 用 data.name
    flat: false,
  },
  iupac: {
    listURL: API.iupac,
    docURL: API.iupacDoc,
    defaultDoc: () => state.iupacCurrent || "目录.md",
    setCurrent: (path) => { state.iupacCurrent = path; },
    titleField: "title", // iupac doc 用 data.title || data.name
    flat: true,
  },
};

/* Wiki：递归目录树 → <details>/<div class="wiki-file"> 字符串 */
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

/* IUPAC：扁平条目列表 → <div class="wiki-file"> 字符串 */
function renderIupacList(entries) {
  return entries.map(function (e) {
    return '<div class="wiki-file" data-path="' + escapeHtml(e.path) + '">'
      + escapeHtml(e.title || e.name) + "</div>";
  }).join("");
}

function setSourceUI(src) {
  var page = $("docs-page");
  if (page) page.setAttribute("data-src", src);
  var seg = page && page.querySelector(".wiki-source-seg");
  if (seg) {
    seg.querySelectorAll(".src-btn").forEach(function (b) {
      b.classList.toggle("active", b.getAttribute("data-src") === src);
    });
  }
}

/* 打开并渲染指定来源的一篇文档；来源列表不变，仅正文更新 */
async function loadDoc(src, path) {
  var meta = SOURCES[src];
  if (!meta) return;
  // wiki 解析相对路径：含 / 按 wiki 根相对；纯文件名按「当前文档同目录」优先、根兜底
  var candidates = [path];
  if (src === "wiki" && path.indexOf("/") === -1) {
    var cur = state.wikiCurrent || "";
    var curDir = cur.indexOf("/") !== -1 ? cur.substring(0, cur.lastIndexOf("/")) : "";
    candidates = curDir ? [curDir + "/" + path, path] : [path];
  }
  var doc = $("docs-doc");
  for (var i = 0; i < candidates.length; i++) {
    var data = null;
    try {
      data = await api(meta.docURL + "?path=" + encodeURIComponent(candidates[i]));
    } catch (e) {
      data = null;
    }
    if (data && data.ok) {
      meta.setCurrent(data.path);
      var body = typeof marked !== "undefined" && marked.parse
        ? marked.parse(data.content)
        : escapeHtml(data.content);
      var title = meta.titleField === "title"
        ? (data.title || data.name)
        : data.name;
      if (doc) doc.innerHTML = '<div class="wiki-title">' + escapeHtml(title) + "</div>" + body;
      document.querySelectorAll("#docs-tree .wiki-file.active").forEach(function (el) {
        el.classList.remove("active");
      });
      var sel = document.querySelector('#docs-tree .wiki-file[data-path="' + CSS.escape(data.path) + '"]');
      if (sel) sel.classList.add("active");
      return;
    }
  }
  if (doc) doc.innerHTML = '<p class="muted-text">无法打开文档: ' + escapeHtml(path) + "</p>";
}

/* 加载指定来源：拉列表/树 → 渲染到 #docs-tree → 打开默认（或记忆的）文档 */
export async function loadDocs(src) {
  var meta = SOURCES[src];
  if (!meta) return;
  state.docsSource = src;
  setSourceUI(src);
  var treeEl = $("docs-tree");
  try {
    var data = await api(meta.listURL);
    if (!data || !data.ok) throw new Error((data && data.error) || "加载失败");
    if (treeEl) {
      treeEl.innerHTML = meta.flat
        ? renderIupacList(data.entries)
        : renderWikiTree(data.tree, 0);
    }
    loadDoc(src, meta.defaultDoc());
  } catch (err) {
    var doc = $("docs-doc");
    if (doc) doc.innerHTML = '<p class="muted-text">加载失败: ' + escapeHtml(err.message || err) + "</p>";
  }
}

export function bindDocs() {
  // 来源切换
  var page = $("docs-page");
  if (page) {
    var seg = page.querySelector(".wiki-source-seg");
    if (seg) {
      seg.addEventListener("click", function (ev) {
        var b = ev.target.closest(".src-btn");
        if (!b) return;
        var src = b.getAttribute("data-src");
        if (src && src !== state.docsSource) loadDocs(src);
      });
    }
  }
  // 点击列表/树条目打开文档
  var tree = $("docs-tree");
  if (tree) {
    tree.addEventListener("click", function (ev) {
      var f = ev.target.closest(".wiki-file");
      if (f) loadDoc(state.docsSource, f.getAttribute("data-path"));
    });
  }
  // 正文内相对 .md 链接在当前来源内跳转；外部 / 锚点放行
  var doc = $("docs-doc");
  if (doc) {
    doc.addEventListener("click", function (ev) {
      var a = ev.target.closest("a");
      if (!a) return;
      var href = (a.getAttribute("href") || "").trim();
      if (!href || /^[a-z]+:/i.test(href) || href.charAt(0) === "#") return;
      if (!/\.md(?:[#?]|$)/i.test(href)) return;
      ev.preventDefault();
      loadDoc(state.docsSource, href.split(/[#?]/)[0]);
    });
  }
}
