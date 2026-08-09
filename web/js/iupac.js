/* IUPAC CN translation page: flat entry list + markdown document viewer.
   Layout mirrors the wiki page; entries are the files under
   docs/iupac/cn_translated (P-XX_中文翻译.md, 目录.md). */
import { $, state, escapeHtml, api, API } from "./core.js";

export async function loadIupac() {
  try {
    var data = await api(API.iupac);
    if (!data || !data.ok) throw new Error((data && data.error) || "加载失败");
    var listEl = $("iupac-list");
    if (listEl) {
      listEl.innerHTML = data.entries.map(function (e) {
        return '<div class="wiki-file" data-path="' + escapeHtml(e.path) + '">'
          + escapeHtml(e.title || e.name) + "</div>";
      }).join("");
    }
    if (state.iupacCurrent) {
      loadIupacDoc(state.iupacCurrent);
    } else {
      loadIupacDoc("目录.md");
    }
  } catch (err) {
    var doc = $("iupac-doc");
    if (doc) doc.innerHTML = '<p class="muted-text">加载失败: ' + escapeHtml(err.message || err) + "</p>";
  }
}

async function loadIupacDoc(path) {
  var doc = $("iupac-doc");
  var data = null;
  try {
    data = await api(API.iupacDoc + "?path=" + encodeURIComponent(path));
  } catch (e) {
    data = null;
  }
  if (data && data.ok) {
    state.iupacCurrent = data.path;
    var body = typeof marked !== "undefined" && marked.parse
      ? marked.parse(data.content)
      : escapeHtml(data.content);
    if (doc) doc.innerHTML = '<div class="wiki-title">' + escapeHtml(data.title || data.name) + "</div>" + body;
    document.querySelectorAll("#iupac-list .wiki-file.active").forEach(function (el) {
      el.classList.remove("active");
    });
    var sel = document.querySelector('#iupac-list .wiki-file[data-path="' + CSS.escape(data.path) + '"]');
    if (sel) sel.classList.add("active");
    return;
  }
  if (doc) doc.innerHTML = '<p class="muted-text">无法打开条目: ' + escapeHtml(path) + "</p>";
}

export function bindIupac() {
  var list = $("iupac-list");
  if (list) {
    list.addEventListener("click", function (ev) {
      var f = ev.target.closest(".wiki-file");
      if (f) loadIupacDoc(f.getAttribute("data-path"));
    });
  }
  var doc = $("iupac-doc");
  if (doc) {
    doc.addEventListener("click", function (ev) {
      var a = ev.target.closest("a");
      if (!a) return;
      var href = (a.getAttribute("href") || "").trim();
      // 外部链接 / 锚点保留默认行为
      if (!href || /^[a-z]+:/i.test(href) || href.charAt(0) === "#") return;
      // 仅拦截指向 cn_translated 内 .md 的相对链接
      if (!/\.md(?:[#?]|$)/i.test(href)) return;
      ev.preventDefault();
      loadIupacDoc(href.split(/[#?]/)[0]);
    });
  }
}
