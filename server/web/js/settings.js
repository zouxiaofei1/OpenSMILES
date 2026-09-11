/* Settings page: theme + model/API config. */
import { $, api, API } from "./core.js";

function settings$(id) {
  return document.getElementById("settings-" + id);
}

function applyTheme(theme) {
  var t = theme === "light" ? "light" : "dark";
  document.documentElement.setAttribute("data-theme", t);
  try { localStorage.setItem("chem-theme", t); } catch (_) {}
}

export function initTheme() {
  var t = null;
  try { t = localStorage.getItem("chem-theme"); } catch (_) {}
  if (t !== "light" && t !== "dark") t = "dark";
  applyTheme(t);
  var sel = settings$("theme");
  if (sel) sel.value = t;
}

function fillApiKeyHint(s) {
  var hint = settings$("api-key-hint");
  if (!hint) return;
  hint.textContent = s && s.api_key_configured
    ? "已配置（……" + (s.api_key_hint || "") + "），留空保持不变"
    : "未配置 API Key";
}

export async function loadSettings() {
  try {
    var data = await api(API.settings);
    if (!data || !data.ok) throw new Error((data && data.error) || "加载失败");
    var s = data.settings || {};
    if (settings$("model")) settings$("model").value = s.model || "";
    if (settings$("base-url")) settings$("base-url").value = s.base_url || "";
    if (settings$("concurrency")) settings$("concurrency").value = s.concurrency || 1;
    if (settings$("theme")) settings$("theme").value = s.theme || "dark";
    fillApiKeyHint(s);
    applyTheme(s.theme || "dark");
  } catch (err) {
    var st = settings$("status");
    if (st) st.textContent = "加载失败: " + (err.message || err);
  }
}

async function saveSettings() {
  var body = {
    model: settings$("model") ? settings$("model").value.trim() : "",
    base_url: settings$("base-url") ? settings$("base-url").value.trim() : "",
    concurrency: parseInt((settings$("concurrency") && settings$("concurrency").value) || "1", 10) || 1,
    theme: settings$("theme") ? settings$("theme").value : "dark",
  };
  var ak = settings$("api-key");
  if (ak && ak.value && ak.value.trim()) body.api_key = ak.value.trim();

  var st = settings$("status");
  if (st) { st.textContent = "保存中…"; st.className = "chip"; }
  try {
    var data = await api(API.settings, { method: "POST", body: JSON.stringify(body) });
    if (!data || !data.ok) throw new Error((data && data.error) || "保存失败");
    if (ak) ak.value = "";
    var s = data.settings || {};
    fillApiKeyHint(s);
    applyTheme(s.theme || "dark");
    if (st) { st.textContent = "已保存 " + new Date().toLocaleTimeString(); st.className = "chip accent"; }
  } catch (err) {
    if (st) { st.textContent = "保存失败: " + (err.message || err); st.className = "chip warn"; }
  }
}

export function bindSettings() {
  var saveBtn = document.getElementById("settings-save");
  if (saveBtn) saveBtn.addEventListener("click", saveSettings);
  var themeSel = settings$("theme");
  if (themeSel) {
    themeSel.addEventListener("change", function () {
      applyTheme(themeSel.value);
    });
  }
}
