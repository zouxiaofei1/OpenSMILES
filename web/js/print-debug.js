/* 打印调试 page: run scripts/debug.py via /api/v1/debug-print and show the
   captured stdout — i.e. every print() executed inside the naming path. */
import { $, escapeHtml, api, API } from "./core.js";

async function runPrintDebug() {
  var input = document.getElementById("print-debug-smiles");
  if (!input) return;
  var smiles = input.value.trim();
  if (!smiles) return;

  var out = document.getElementById("print-debug-output");
  if (!out) return;
  out.classList.remove("has-error");
  out.textContent = "正在执行 scripts/debug.py ...";

  try {
    var res = await api(API.debugPrint, { method: "POST", body: JSON.stringify({ smiles: smiles }) });
    if (res.error) {
      out.classList.add("has-error");
      out.textContent = res.error;
      return;
    }
    var txt = res.stdout || "";
    if (res.stderr) {
      txt += (txt ? "\n" : "") + "--- stderr (returncode=" + res.returncode + ") ---\n" + res.stderr;
      out.classList.add("has-error");
    }
    out.textContent = txt || "(无输出)";
  } catch (err) {
    out.classList.add("has-error");
    out.textContent = "Error: " + escapeHtml(String(err));
  }
}

export function bindPrintDebug() {
  var input = document.getElementById("print-debug-smiles");
  var btn = document.getElementById("print-debug-run");
  if (btn) btn.addEventListener("click", runPrintDebug);
  if (input) {
    input.addEventListener("keydown", function (e) {
      if (e.key === "Enter") runPrintDebug();
    });
  }
}
