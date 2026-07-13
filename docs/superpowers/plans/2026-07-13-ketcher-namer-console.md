# Ketcher Namer Console Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 在 ChemAgent Console 的 Namer tab 嵌入本地自托管 Ketcher 画板，支持画分子后调用 `POST /api/v1/name` 得到中英 IUPAC 名，并提供可选实时命名。

**Architecture:** 用一次性 Vite 壳工程（`tools/ketcher-shell/`）把 `ketcher-react` + `ketcher-standalone` 打成静态页，产出进 `web/vendor/ketcher/`（含 `index.html` + assets + `VERSION`）。Console 经同源 iframe 懒加载该页；父页 `KetcherBridge` 调用 `iframe.contentWindow.ketcher` 的 `getSmiles` / `setMolecule` / `changeEvent`。Namer 双栏布局；Name 优先取画板 SMILES；实时开关默认关、700ms debounce + `nameReqSeq` 防乱序。后端零改动。

**Tech Stack:** 现有静态 web（HTML/CSS/JS）、FastAPI 已挂载 `web/`、Ketcher 3.7.0（`ketcher-react` + `ketcher-standalone`）、Vite（仅 vendor 构建，不进入 Console 运行时）、Windows Git Bash。

## Global Constraints

- 形态：接入现有 Console Namer（非独立应用）
- Ketcher：本地自托管 → `web/vendor/ketcher/`；锁定版本写 `VERSION`
- 命名：手动 Name + 可选实时（默认关，debounce 700ms）
- 画板 → 输入框（Name/实时）；文本 → 画板仅「从 SMILES 载入」
- 复用 `POST /api/v1/name`；**默认不改** `server/*`
- Console 运行时**零 npm 构建**；vendor 由 `tools/ketcher-shell` 构建后提交
- vendor 缺失时文本 SMILES 命名仍可用
- 用户可见中文；代码标识符英文；`127.0.0.1:8765`
- 工作目录：`E:\dev\chem`
- venv：`source .venv/Scripts/activate`（本任务前端为主，一般不必）

---

## File Structure

```
tools/ketcher-shell/                 # NEW: 一次性构建壳（可 gitignore node_modules）
  package.json
  vite.config.js
  index.html
  src/main.jsx
  src/KetcherApp.jsx
  README.md                          # 如何 rebuild vendor
tools/build_ketcher_vendor.sh        # NEW: 构建并同步到 web/vendor/ketcher
web/vendor/ketcher/                  # NEW: 构建产物（提交）
  VERSION
  index.html
  assets/*
web/index.html                       # MOD: Namer 双栏 DOM
web/css/app.css                      # MOD: namer-layout / iframe
web/js/namer-ketcher.js              # NEW: KetcherBridge + 纯逻辑
web/js/app.js                        # MOD: 接线 Name / live / lazy tab
README.md                            # MOD: 一句说明
```

**不修改：** `server/*`、`src/namepredict/**`（除非静态托管实测失败再最小补丁）。

---

### Task 1: Ketcher vendor 壳 + 构建产物

**Files:**
- Create: `tools/ketcher-shell/package.json`
- Create: `tools/ketcher-shell/vite.config.js`
- Create: `tools/ketcher-shell/index.html`
- Create: `tools/ketcher-shell/src/main.jsx`
- Create: `tools/ketcher-shell/src/KetcherApp.jsx`
- Create: `tools/ketcher-shell/README.md`
- Create: `tools/build_ketcher_vendor.sh`
- Create: `web/vendor/ketcher/**`（构建输出）
- Create: `web/vendor/ketcher/VERSION`

**Interfaces:**
- iframe 页加载完成后：`window.ketcher` 为 Ketcher 实例（与 ketcher-core API 一致）
- `window.ketcher.getSmiles(): Promise<string>`
- `window.ketcher.setMolecule(structStr: string): Promise<void>`
- `window.ketcher.editor.clear(): void`（或等价清空；若无则 `setMolecule('')`）
- `window.ketcher.changeEvent`：可 `add(handler)` / 订阅结构变化（以 3.7 API 为准；若无则 bridge 轮询）
- `VERSION` 文件单行：`3.7.0`

- [ ] **Step 1: 创建壳 package.json**

```json
{
  "name": "ketcher-shell",
  "private": true,
  "version": "0.1.0",
  "type": "module",
  "scripts": {
    "build": "vite build"
  },
  "dependencies": {
    "ketcher-react": "3.7.0",
    "ketcher-standalone": "3.7.0",
    "react": "^18.3.1",
    "react-dom": "^18.3.1"
  },
  "devDependencies": {
    "@vitejs/plugin-react": "^4.3.4",
    "vite": "^5.4.11"
  }
}
```

- [ ] **Step 2: vite.config.js（base 指向 /vendor/ketcher/）**

```js
import { defineConfig } from "vite";
import react from "@vitejs/plugin-react";

export default defineConfig({
  plugins: [react()],
  base: "/vendor/ketcher/",
  build: {
    outDir: "dist",
    emptyOutDir: true,
    assetsDir: "assets",
  },
});
```

- [ ] **Step 3: index.html + React 入口**

`tools/ketcher-shell/index.html`:

```html
<!DOCTYPE html>
<html lang="en">
  <head>
    <meta charset="UTF-8" />
    <meta name="viewport" content="width=device-width, initial-scale=1.0" />
    <title>Ketcher</title>
    <style>
      html, body, #root { height: 100%; margin: 0; }
    </style>
  </head>
  <body>
    <div id="root"></div>
    <script type="module" src="/src/main.jsx"></script>
  </body>
</html>
```

`tools/ketcher-shell/src/main.jsx`:

```jsx
import React from "react";
import { createRoot } from "react-dom/client";
import "ketcher-react/dist/index.css";
import { KetcherApp } from "./KetcherApp.jsx";

createRoot(document.getElementById("root")).render(
  <React.StrictMode>
    <KetcherApp />
  </React.StrictMode>
);
```

`tools/ketcher-shell/src/KetcherApp.jsx`:

```jsx
import React, { useCallback } from "react";
import { Editor } from "ketcher-react";
import { StandaloneStructServiceProvider } from "ketcher-standalone";

const structServiceProvider = new StandaloneStructServiceProvider();

export function KetcherApp() {
  const onInit = useCallback((ketcher) => {
    window.ketcher = ketcher;
    window.dispatchEvent(new Event("ketcher-ready"));
  }, []);

  return (
    <div style={{ width: "100%", height: "100%" }}>
      <Editor
        staticResourcesUrl="."
        structServiceProvider={structServiceProvider}
        onInit={onInit}
      />
    </div>
  );
}
```

若 `Editor` 的 `onInit` prop 名称在 3.7 不同，查 `node_modules/ketcher-react` 类型/README，改成实际回调（常见：`onInit`）。

- [ ] **Step 4: 构建脚本**

`tools/build_ketcher_vendor.sh`:

```bash
#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SHELL_DIR="$ROOT/tools/ketcher-shell"
OUT="$ROOT/web/vendor/ketcher"
VER="3.7.0"

cd "$SHELL_DIR"
if [[ ! -d node_modules ]]; then
  npm install
fi
npm run build

rm -rf "$OUT"
mkdir -p "$OUT"
cp -R "$SHELL_DIR/dist/." "$OUT/"
echo "$VER" > "$OUT/VERSION"
echo "Wrote $OUT (ketcher $VER)"
```

`tools/ketcher-shell/README.md`：说明需 Node 18+、`bash tools/build_ketcher_vendor.sh` 重生成 vendor。

- [ ] **Step 5: 运行构建并检查产物**

```bash
bash tools/build_ketcher_vendor.sh
test -f web/vendor/ketcher/index.html
test -f web/vendor/ketcher/VERSION
cat web/vendor/ketcher/VERSION
# Expected: 3.7.0
ls web/vendor/ketcher/assets | head
```

若 `npm install` 失败（registry/网络），改用可用镜像后重试；勿提交 `tools/ketcher-shell/node_modules`。

- [ ] **Step 6: .gitignore（若需要）**

确认 `tools/ketcher-shell/node_modules/` 被 ignore（可加一行到根 `.gitignore`）：

```
tools/ketcher-shell/node_modules/
tools/ketcher-shell/dist/
```

**不要** ignore `web/vendor/ketcher/`。

- [ ] **Step 7: Commit**

```bash
git add tools/ketcher-shell tools/build_ketcher_vendor.sh web/vendor/ketcher .gitignore
git commit -m "$(cat <<'EOF'
feat(web): vendor Ketcher 3.7 standalone shell

Add Vite shell and build script; commit static assets under
web/vendor/ketcher for Console iframe embedding.

Co-Authored-By: Claude <noreply@anthropic.com>
EOF
)"
```

---

### Task 2: Namer 双栏 DOM

**Files:**
- Modify: `web/index.html`（`#panel-namer` 整段）

**Interfaces:**
- DOM ids（后续 JS 依赖）：
  - `namer-layout`
  - `ketcher-frame`（iframe，初始 **无 src** 或 `src="about:blank"`）
  - `ketcher-fallback`（vendor 失败提示，默认 hidden）
  - `btn-load-smiles`、`btn-clear-ketcher`
  - `live-name-toggle`、`live-name-hint`
  - 保留：`namer-form`、`smiles-input`、`btn-name`、`namer-error`、`namer-result`、`namer-en/zh/source/time`、`namer-history`

- [ ] **Step 1: 替换 panel-namer 结构**

将 `web/index.html` 中 `<!-- Namer -->` … `</section>` 换成：

```html
          <!-- Namer -->
          <section id="panel-namer" class="tab-panel" role="tabpanel" aria-labelledby="tab-namer" hidden>
            <div id="namer-layout" class="namer-layout">
              <div class="card namer-structure-card">
                <div class="card-head row-between">
                  <h3 class="card-title">Structure</h3>
                  <div class="namer-structure-actions">
                    <button type="button" id="btn-load-smiles" class="btn btn-secondary btn-sm" disabled>
                      从 SMILES 载入
                    </button>
                    <button type="button" id="btn-clear-ketcher" class="btn btn-ghost btn-sm" disabled>
                      清空画板
                    </button>
                  </div>
                </div>
                <div class="ketcher-wrap">
                  <iframe
                    id="ketcher-frame"
                    class="ketcher-frame"
                    title="Ketcher molecule editor"
                    src="about:blank"
                  ></iframe>
                  <div id="ketcher-fallback" class="ketcher-fallback" hidden>
                    <p class="empty-title">Ketcher 未安装</p>
                    <p class="empty-sub">请将构建产物放到 <code>web/vendor/ketcher/</code>，或运行 <code>bash tools/build_ketcher_vendor.sh</code>。文本 SMILES 命名仍可用。</p>
                  </div>
                </div>
              </div>
              <div class="namer-side">
                <div class="card namer-card">
                  <h3 class="card-title">SMILES Namer</h3>
                  <form id="namer-form" class="namer-form" autocomplete="off">
                    <label class="field-label" for="smiles-input">SMILES</label>
                    <div class="namer-row">
                      <input
                        id="smiles-input"
                        name="smiles"
                        class="input mono"
                        type="text"
                        placeholder="e.g. CCO"
                        required
                        spellcheck="false"
                        aria-describedby="namer-error live-name-hint"
                      />
                      <button type="submit" id="btn-name" class="btn btn-primary" aria-label="Name SMILES">
                        <svg class="icon" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true">
                          <path d="M10 2v7.527a2 2 0 0 1-.211.896L4.72 20.55a1 1 0 0 0 .9 1.45h12.76a1 1 0 0 0 .9-1.45l-5.069-10.127A2 2 0 0 1 14 9.527V2" />
                          <path d="M8.5 2h7" />
                          <path d="M7 16h10" />
                        </svg>
                        <span>Name</span>
                      </button>
                    </div>
                    <label class="live-toggle">
                      <input type="checkbox" id="live-name-toggle" />
                      <span>实时命名</span>
                    </label>
                    <p id="live-name-hint" class="muted-text small">结构变化后约 0.7s 自动命名（默认关）</p>
                    <p id="namer-error" class="field-error" role="alert" hidden></p>
                  </form>
                  <div id="namer-result" class="namer-result hidden" aria-live="polite">
                    <div class="result-head">
                      <span id="namer-success" class="badge">—</span>
                      <span id="namer-time" class="muted-text small mono"></span>
                    </div>
                    <dl class="kv">
                      <div><dt>en</dt><dd id="namer-en" class="mono">—</dd></div>
                      <div><dt>zh</dt><dd id="namer-zh" class="mono">—</dd></div>
                      <div><dt>source</dt><dd id="namer-source" class="mono">—</dd></div>
                    </dl>
                  </div>
                </div>
                <div class="card">
                  <div class="card-head row-between">
                    <h3 class="card-title">Recent (20)</h3>
                    <button type="button" id="btn-clear-history" class="btn btn-ghost btn-sm">清空</button>
                  </div>
                  <ul id="namer-history" class="history-list">
                    <li class="muted-text small">尚无命名记录</li>
                  </ul>
                </div>
              </div>
            </div>
          </section>
```

在 `</body>` 前、现有 `app.js` **之前**增加：

```html
  <script src="/js/namer-ketcher.js"></script>
```

（`app.js` 保持最后加载。）

- [ ] **Step 2: 静态检查**

```bash
rg -n "namer-layout|ketcher-frame|live-name-toggle|namer-ketcher" web/index.html
```

Expected: 均有匹配。

- [ ] **Step 3: Commit**

```bash
git add web/index.html
git commit -m "$(cat <<'EOF'
feat(web): Namer dual-pane DOM for Ketcher

Add structure card iframe, live-name toggle, and load/clear controls.

Co-Authored-By: Claude <noreply@anthropic.com>
EOF
)"
```

---

### Task 3: Namer 布局 CSS

**Files:**
- Modify: `web/css/app.css`（追加在文件末尾或 namer 区块附近）

- [ ] **Step 1: 追加样式**

```css
/* --- Namer + Ketcher --- */
.namer-layout {
  display: grid;
  grid-template-columns: minmax(0, 1.4fr) minmax(280px, 1fr);
  gap: var(--space-xl);
  align-items: start;
}

.namer-side {
  display: flex;
  flex-direction: column;
  gap: var(--space-xl);
  min-width: 0;
}

.namer-structure-card {
  min-width: 0;
}

.namer-structure-actions {
  display: flex;
  gap: var(--space-sm);
  flex-wrap: wrap;
}

.ketcher-wrap {
  position: relative;
  width: 100%;
  height: 480px;
  border: 1px solid var(--color-border);
  border-radius: var(--radius);
  overflow: hidden;
  background: #fff;
}

.ketcher-frame {
  width: 100%;
  height: 100%;
  border: 0;
  display: block;
  background: #fff;
}

.ketcher-fallback {
  position: absolute;
  inset: 0;
  display: flex;
  flex-direction: column;
  align-items: center;
  justify-content: center;
  padding: var(--space-2xl);
  text-align: center;
  background: var(--color-muted);
  color: var(--color-foreground);
}

.live-toggle {
  display: inline-flex;
  align-items: center;
  gap: var(--space-sm);
  font-size: 13px;
  color: var(--color-foreground);
  cursor: pointer;
  user-select: none;
}

.live-toggle input {
  width: 16px;
  height: 16px;
  accent-color: var(--color-accent);
}

@media (max-width: 959px) {
  .namer-layout {
    grid-template-columns: 1fr;
  }
  .ketcher-wrap {
    height: 420px;
  }
}
```

- [ ] **Step 2: Commit**

```bash
git add web/css/app.css
git commit -m "$(cat <<'EOF'
style(web): Namer dual-pane and Ketcher frame layout

Co-Authored-By: Claude <noreply@anthropic.com>
EOF
)"
```

---

### Task 4: `namer-ketcher.js` Bridge + 纯逻辑

**Files:**
- Create: `web/js/namer-ketcher.js`

**Interfaces:**
- 全局：`window.ChemNamerKetcher = { createBridge, shouldSkipLiveName, nextReqSeq }`
- `createBridge({ iframe, src, onReady, onError, onChange })` →
  `{ init(), isReady(), getSmiles(), setMolecule(s), clear(), destroy() }`
- `shouldSkipLiveName(liveEnabled, smiles, lastNamed) -> boolean`
- `nextReqSeq(current) -> number`

- [ ] **Step 1: 实现模块**

```js
/* ChemAgent Console — Ketcher bridge for Namer */
(function (global) {
  "use strict";

  const DEFAULT_SRC = "/vendor/ketcher/index.html";
  const READY_POLL_MS = 200;
  const READY_TIMEOUT_MS = 20000;

  function shouldSkipLiveName(liveEnabled, smiles, lastNamed) {
    if (!liveEnabled) return true;
    const s = (smiles || "").trim();
    if (!s) return true;
    return s === (lastNamed || "");
  }

  function nextReqSeq(current) {
    const n = Number(current) || 0;
    return n + 1;
  }

  function createBridge(opts) {
    const iframe = opts.iframe;
    const src = opts.src || DEFAULT_SRC;
    const onReady = opts.onReady || function () {};
    const onError = opts.onError || function () {};
    const onChange = opts.onChange || function () {};

    let ready = false;
    let ketcher = null;
    let pollTimer = null;
    let changeUnsub = null;
    let destroyed = false;

    function isReady() {
      return ready && !!ketcher;
    }

    function getKetcherFromFrame() {
      try {
        return iframe.contentWindow && iframe.contentWindow.ketcher;
      } catch (_) {
        return null;
      }
    }

    function attachChange(k) {
      if (!k || !k.changeEvent) return;
      const handler = function () {
        onChange();
      };
      try {
        if (typeof k.changeEvent.add === "function") {
          k.changeEvent.add(handler);
          changeUnsub = function () {
            if (typeof k.changeEvent.remove === "function") k.changeEvent.remove(handler);
          };
        } else if (typeof k.changeEvent.subscribe === "function") {
          const sub = k.changeEvent.subscribe(handler);
          changeUnsub = function () {
            if (sub && typeof sub.unsubscribe === "function") sub.unsubscribe();
          };
        }
      } catch (_) {
        /* optional */
      }
    }

    function markReady(k) {
      if (destroyed || ready) return;
      ketcher = k;
      ready = true;
      if (pollTimer) {
        clearInterval(pollTimer);
        pollTimer = null;
      }
      attachChange(k);
      onReady(k);
    }

    function startPolling() {
      const t0 = Date.now();
      pollTimer = setInterval(function () {
        if (destroyed) return;
        const k = getKetcherFromFrame();
        if (k) {
          markReady(k);
          return;
        }
        if (Date.now() - t0 > READY_TIMEOUT_MS) {
          clearInterval(pollTimer);
          pollTimer = null;
          onError(new Error("Ketcher 加载超时"));
        }
      }, READY_POLL_MS);
    }

    function onFrameLoad() {
      if (destroyed) return;
      const win = iframe.contentWindow;
      if (!win) {
        onError(new Error("Ketcher iframe 不可用"));
        return;
      }
      const existing = getKetcherFromFrame();
      if (existing) {
        markReady(existing);
        return;
      }
      function onKetcherReady() {
        const k = getKetcherFromFrame();
        if (k) markReady(k);
      }
      try {
        win.addEventListener("ketcher-ready", onKetcherReady);
      } catch (_) {
        /* ignore */
      }
      startPolling();
    }

    function init() {
      if (destroyed) return;
      if (iframe.dataset.ketcherInit === "1") return;
      iframe.dataset.ketcherInit = "1";
      iframe.addEventListener("load", onFrameLoad);
      iframe.addEventListener("error", function () {
        onError(new Error("Ketcher 资源加载失败"));
      });
      // probe vendor existence via fetch then set src
      fetch(src, { method: "GET", cache: "no-cache" })
        .then(function (res) {
          if (!res.ok) throw new Error("HTTP " + res.status);
          iframe.src = src;
        })
        .catch(function (err) {
          onError(err);
        });
    }

    async function getSmiles() {
      if (!isReady()) return "";
      try {
        const s = await ketcher.getSmiles();
        return (s || "").trim();
      } catch (_) {
        return "";
      }
    }

    async function setMolecule(smiles) {
      if (!isReady()) throw new Error("Ketcher 未就绪");
      await ketcher.setMolecule(smiles || "");
    }

    async function clear() {
      if (!isReady()) return;
      try {
        if (ketcher.editor && typeof ketcher.editor.clear === "function") {
          ketcher.editor.clear();
          return;
        }
      } catch (_) {
        /* fall through */
      }
      await ketcher.setMolecule("");
    }

    function destroy() {
      destroyed = true;
      if (pollTimer) clearInterval(pollTimer);
      if (changeUnsub) {
        try {
          changeUnsub();
        } catch (_) {}
      }
      ready = false;
      ketcher = null;
    }

    return {
      init: init,
      isReady: isReady,
      getSmiles: getSmiles,
      setMolecule: setMolecule,
      clear: clear,
      destroy: destroy,
    };
  }

  global.ChemNamerKetcher = {
    createBridge: createBridge,
    shouldSkipLiveName: shouldSkipLiveName,
    nextReqSeq: nextReqSeq,
    DEFAULT_SRC: DEFAULT_SRC,
  };
})(typeof window !== "undefined" ? window : globalThis);
```

- [ ] **Step 2: 浏览器控制台快速自检（可选，server 起来后）**

```js
// 在 Console 页
ChemNamerKetcher.shouldSkipLiveName(false, "CCO", "") === true
ChemNamerKetcher.shouldSkipLiveName(true, "CCO", "CCO") === true
ChemNamerKetcher.shouldSkipLiveName(true, "CCO", "CC") === false
ChemNamerKetcher.nextReqSeq(0) === 1
```

- [ ] **Step 3: Commit**

```bash
git add web/js/namer-ketcher.js
git commit -m "$(cat <<'EOF'
feat(web): KetcherBridge for Console Namer

Expose createBridge, live-name skip helper, and request seq helper.

Co-Authored-By: Claude <noreply@anthropic.com>
EOF
)"
```

---

### Task 5: 接线 `app.js`（lazy init、Name、载入/清空、实时）

**Files:**
- Modify: `web/js/app.js`

**Interfaces:**
- 扩展 `state`：`liveNameEnabled`, `ketcherReady`, `lastNamedSmiles`, `nameReqSeq`, `ketcherBridge`, `liveDebounceTimer`
- `ensureKetcher()`：首次进 Namer tab 调用
- `resolveSmilesForName()`：画板优先，否则输入框
- `runName(smiles, { fromLive })`：POST + seq + history（live 失败不炸 UI）
- 修改 `onName`、`switchTab`、`bind`

- [ ] **Step 1: 扩展 state**

在 `const state = { ... }` 内追加字段：

```js
    liveNameEnabled: false,
    ketcherReady: false,
    lastNamedSmiles: "",
    nameReqSeq: 0,
    ketcherBridge: null,
    liveDebounceTimer: null,
```

- [ ] **Step 2: 增加 namer/ketcher 辅助函数**（放在 `/* ---------- Namer ---------- */` 区块内，`renderNamerHistory` 之前）

```js
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

  function showKetcherFallback(show) {
    const fb = $("ketcher-fallback");
    const frame = $("ketcher-frame");
    if (fb) fb.hidden = !show;
    if (frame) frame.style.visibility = show ? "hidden" : "visible";
  }

  function ensureKetcher() {
    if (!window.ChemNamerKetcher) {
      showKetcherFallback(true);
      return;
    }
    if (state.ketcherBridge) {
      state.ketcherBridge.init();
      return;
    }
    const iframe = $("ketcher-frame");
    if (!iframe) return;
    state.ketcherBridge = window.ChemNamerKetcher.createBridge({
      iframe,
      src: window.ChemNamerKetcher.DEFAULT_SRC,
      onReady: () => {
        state.ketcherReady = true;
        setKetcherControlsEnabled(true);
        showKetcherFallback(false);
        appendLog("ketcher", { message: "ready" });
      },
      onError: (err) => {
        state.ketcherReady = false;
        setKetcherControlsEnabled(false);
        showKetcherFallback(true);
        appendLog("error", { message: "ketcher: " + (err && err.message ? err.message : String(err)) });
      },
      onChange: () => {
        scheduleLiveName();
      },
    });
    state.ketcherBridge.init();
  }

  async function resolveSmilesForName() {
    let fromEditor = "";
    if (state.ketcherBridge && state.ketcherBridge.isReady()) {
      fromEditor = await state.ketcherBridge.getSmiles();
    }
    if (fromEditor) {
      const input = $("smiles-input");
      if (input) input.value = fromEditor;
      return fromEditor;
    }
    return (($("smiles-input") && $("smiles-input").value) || "").trim();
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
      appendLog("name", { smiles, en: result.en, zh: result.zh, success: result.success, live: fromLive });
    } catch (err) {
      if (seq !== state.nameReqSeq) return;
      setNamerError(err.message || String(err));
      appendLog("error", { message: "name: " + (err.message || String(err)) });
    } finally {
      if (!fromLive && btn) btn.disabled = false;
    }
  }

  function scheduleLiveName() {
    if (!state.liveNameEnabled) return;
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
        if (input) input.value = smiles;
      }
      await runName(smiles, { fromLive: true });
    }, 700);
  }
```

- [ ] **Step 3: 替换 `onName`**

```js
  async function onName(ev) {
    ev.preventDefault();
    const smiles = await resolveSmilesForName();
    await runName(smiles, { fromLive: false });
  }
```

- [ ] **Step 4: 修改 `switchTab`**

在 `switchTab` 末尾增加：

```js
    if (name === "namer") ensureKetcher();
```

- [ ] **Step 5: 扩展 `bind`**

在现有 namer-form / clear-history 绑定旁增加：

```js
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
          await state.ketcherBridge.setMolecule(smiles);
        } catch (err) {
          setNamerError(err.message || "载入结构失败");
        }
      });

    $("btn-clear-ketcher") &&
      $("btn-clear-ketcher").addEventListener("click", async () => {
        if (!state.ketcherBridge || !state.ketcherBridge.isReady()) return;
        try {
          await state.ketcherBridge.clear();
        } catch (err) {
          setNamerError(err.message || "清空失败");
        }
      });

    $("live-name-toggle") &&
      $("live-name-toggle").addEventListener("change", (ev) => {
        state.liveNameEnabled = !!(ev.target && ev.target.checked);
        if (!state.liveNameEnabled && state.liveDebounceTimer) {
          clearTimeout(state.liveDebounceTimer);
          state.liveDebounceTimer = null;
        }
      });
```

- [ ] **Step 6: 手工冒烟**

```bash
source .venv/Scripts/activate
uvicorn server.app:app --host 127.0.0.1 --port 8765
```

浏览器打开 `http://127.0.0.1:8765/` → Namer：

1. iframe 出现 Ketcher（或 fallback）
2. 画 CCO → Name → 有 en/zh
3. 输入框贴 SMILES → 从 SMILES 载入 → 画板有结构
4. 实时 on/off 行为符合设计

- [ ] **Step 7: Commit**

```bash
git add web/js/app.js
git commit -m "$(cat <<'EOF'
feat(web): wire Ketcher into Namer name flow

Lazy-init iframe, prefer editor SMILES on Name, load/clear,
and optional 700ms live naming with request sequencing.

Co-Authored-By: Claude <noreply@anthropic.com>
EOF
)"
```

---

### Task 6: README + 验收收尾

**Files:**
- Modify: `README.md`（Console / API 小节加 2–4 行）

- [ ] **Step 1: README 补充**

在 Console URL 说明附近增加：

```markdown
Namer tab supports drawing structures via self-hosted Ketcher
(`web/vendor/ketcher/`, version in `VERSION`). Rebuild vendor:

```bash
bash tools/build_ketcher_vendor.sh
```

Requires Node 18+. Text SMILES naming works even if vendor is missing.
```

- [ ] **Step 2: 按 spec 验收清单跑一遍**

1. 画 CCO → Name  
2. 粘贴 SMILES → 载入  
3. 实时 on/off  
4. 临时改 iframe 路径或移走 vendor 验证 fallback + 文本 Name  

- [ ] **Step 3: Commit**

```bash
git add README.md
git commit -m "$(cat <<'EOF'
docs: note Ketcher vendor in Console Namer

Co-Authored-By: Claude <noreply@anthropic.com>
EOF
)"
```

---

## Spec coverage (self-review)

| Spec 要求 | Task |
|-----------|------|
| 接入 Console Namer 双栏 | 2, 3 |
| 本地自托管 vendor | 1 |
| 手动 Name，画板优先 SMILES | 5 |
| 实时默认关 + 700ms debounce | 5 |
| 从 SMILES 载入 / 清空画板 | 5 |
| vendor 缺失降级 | 4, 5 |
| 不改 name API | 全任务 |
| README | 6 |
| lazy init on Namer tab | 5 |
| seq 防乱序 / lastNamed 去重 | 4, 5 |

## Placeholder / consistency

- 版本统一 **3.7.0**（shell deps + VERSION + 脚本）
- API：`getSmiles` / `setMolecule` / `changeEvent` / `editor.clear`
- 若 3.7 `Editor` 无 `onInit`，Task 1 内查包文档改 prop，不改变父页 Bridge 契约
- 无 TBD 步骤

## 风险处理（实现时）

- npm 装 ketcher 失败：换 registry；或降到已知可装的 2.28.x 并同步改 VERSION
- `staticResourcesUrl` 路径：若图标 404，改为 `base` 相对路径或 `import.meta.env.BASE_URL`
- `changeEvent` 不可用：Bridge 内 fallback 为 1s 轮询 `getSmiles` 对比（仅 live 开启时），实现时若需要可在 Task 5 的 `onChange` 路径补上，仍保持 debounce 700ms
