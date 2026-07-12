# Console API Key + Full Bench Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 本机 ChemAgent Console 可配置 Anthropic API Key（落盘 secrets、无 Key 禁止 Start、真 pi 注入 env），并在 Bench 页一键跑全量 benchmark 且 SSE 分阶段日志。

**Architecture:** `agent_loop/secrets.py` 读写 `agent_loop/memory/secrets.json`；Settings/Bench 路由挂 FastAPI；`LoopController.start` 门禁 + `mock_pi=False` + `PiRunner` 子进程 env；`BenchController` 单飞线程跑 `python -m benchmarks.benchmark --json`（全量）并经 `EventBus` 发 `bench.stage`/`bench.done`；web 顶栏 Settings + Bench 工具条。

**Tech Stack:** Python 3.11+、FastAPI、pytest、现有 EventBus/SSE、静态 web（HTML/CSS/JS）。

## Global Constraints

- 仅 Anthropic：`ANTHROPIC_API_KEY`；secrets 路径 `agent_loop/memory/secrets.json` 且 gitignore
- GET/日志/SSE/session **永不**含完整 Key；hint 最多末 4 位
- 未配置 Key → Start HTTP 400 `missing_api_key` + 前端弹窗
- 有 Key → 真 pi（`mock_pi=False`），子进程 env 注入 Key
- 手动 Bench **始终全量** `data/merged_benchmark.json`（`limit=None`）
- 同时只允许 1 个手动 bench；与 loop **允许并行**
- 不改 benchmark 金标与计分逻辑；不改 namepredict 规则业务
- 单函数/方法体 ≤10 行；pathlib；Windows Git Bash；venv：`source .venv/Scripts/activate`
- 用户可见中文；代码标识符英文；Console 文档约定 `127.0.0.1`
- 工作目录：`E:\dev\chem`

---

## File Structure

```
agent_loop/secrets.py          # NEW: load/save/clear/status
agent_loop/pi_runner.py        # MOD: env injection for ANTHROPIC_API_KEY
server/routes_settings.py      # NEW
server/routes_bench.py         # NEW
server/bench_runner.py         # NEW: BenchController
server/deps.py                 # MOD: start gate, real pi, secrets path
server/app.py                  # MOD: include routers
web/index.html                 # MOD: settings modal + bench toolbar
web/css/app.css                # MOD: styles
web/js/app.js                  # MOD: settings + start gate + bench SSE
.gitignore                     # MOD: secrets.json
tests/unit/test_secrets.py
tests/unit/test_settings_api.py
tests/unit/test_bench_api.py
tests/unit/test_pi_runner.py   # MOD: env assertion
README.md                      # MOD: short notes
```

---

### Task 1: secrets 模块 + gitignore

**Files:**
- Create: `agent_loop/secrets.py`
- Create: `tests/unit/test_secrets.py`
- Modify: `.gitignore`

**Interfaces:**
- `DEFAULT_SECRETS_PATH = Path("agent_loop/memory/secrets.json")`
- `class SecretsStore:`
  - `__init__(self, path: Path | None = None)`
  - `is_configured(self) -> bool`
  - `get_key(self) -> str | None`  # stripped non-empty or None
  - `set_key(self, key: str) -> None`  # empty/whitespace → clear
  - `clear_key(self) -> None`
  - `public_status(self) -> dict`  # `{"configured": bool, "hint": str|None}` hint = last 4 chars if configured

- [ ] **Step 1: Write failing tests**

```python
# tests/unit/test_secrets.py
from pathlib import Path
from agent_loop.secrets import SecretsStore

def test_set_get_configured(tmp_path: Path):
    store = SecretsStore(tmp_path / "secrets.json")
    assert store.is_configured() is False
    assert store.public_status() == {"configured": False, "hint": None}
    store.set_key("sk-ant-test-1234")
    assert store.is_configured() is True
    assert store.get_key() == "sk-ant-test-1234"
    st = store.public_status()
    assert st["configured"] is True
    assert st["hint"] == "1234"
    assert "sk-ant" not in (st["hint"] or "")

def test_clear_and_empty_set(tmp_path: Path):
    store = SecretsStore(tmp_path / "secrets.json")
    store.set_key("abcd")
    store.clear_key()
    assert store.get_key() is None
    store.set_key("  ")
    assert store.is_configured() is False

def test_corrupt_file_unconfigured(tmp_path: Path):
    path = tmp_path / "secrets.json"
    path.write_text("{not json", encoding="utf-8")
    store = SecretsStore(path)
    assert store.is_configured() is False
```

- [ ] **Step 2: Run tests — expect FAIL**

```bash
source .venv/Scripts/activate
pytest -q tests/unit/test_secrets.py
```
Expected: import error or fail

- [ ] **Step 3: Implement `agent_loop/secrets.py`**

```python
"""Local Anthropic API key store (never log full key)."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

DEFAULT_SECRETS_PATH = Path("agent_loop/memory/secrets.json")
_KEY = "anthropic_api_key"


class SecretsStore:
    def __init__(self, path: Path | None = None) -> None:
        self.path = Path(path or DEFAULT_SECRETS_PATH)

    def is_configured(self) -> bool:
        return self.get_key() is not None

    def get_key(self) -> str | None:
        data = self._read()
        raw = data.get(_KEY)
        if not isinstance(raw, str):
            return None
        key = raw.strip()
        return key or None

    def set_key(self, key: str) -> None:
        text = (key or "").strip()
        if not text:
            self.clear_key()
            return
        self._write({_KEY: text})

    def clear_key(self) -> None:
        self._write({})

    def public_status(self) -> dict[str, Any]:
        key = self.get_key()
        if key is None:
            return {"configured": False, "hint": None}
        return {"configured": True, "hint": key[-4:] if len(key) >= 4 else key}

    def _read(self) -> dict[str, Any]:
        if not self.path.is_file():
            return {}
        try:
            data = json.loads(self.path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            return {}
        return data if isinstance(data, dict) else {}

    def _write(self, data: dict[str, Any]) -> None:
        self.path.parent.mkdir(parents=True, exist_ok=True)
        text = json.dumps(data, indent=2, ensure_ascii=False) + "\n"
        self.path.write_text(text, encoding="utf-8")
```

Keep each method body ≤10 lines (split helpers if needed).

- [ ] **Step 4: Update `.gitignore`**

Add line:
```
agent_loop/memory/secrets.json
```

- [ ] **Step 5: Run tests — expect PASS**

```bash
pytest -q tests/unit/test_secrets.py
```

- [ ] **Step 6: Commit**

```bash
git add agent_loop/secrets.py tests/unit/test_secrets.py .gitignore
git commit -m "feat: local Anthropic secrets store with public status"
```

---

### Task 2: Settings API + Start 门禁 + 真 pi env

**Files:**
- Create: `server/routes_settings.py`
- Create: `tests/unit/test_settings_api.py`
- Modify: `server/app.py`
- Modify: `server/deps.py` (`LoopController.start` / `_make_agent`)
- Modify: `server/routes_loop.py` (propagate 400)
- Modify: `agent_loop/pi_runner.py`
- Modify: `tests/unit/test_pi_runner.py` (env assert)
- Modify: `tests/unit/test_loop_controller.py` if start behavior changes

**Interfaces:**
- Settings: `GET/POST /api/v1/settings/api-key`
- `LoopController.start() -> dict` may include error; or raise / return status that route maps to 400
- Prefer: `start()` returns `{"status": "...", "error": "missing_api_key", "message": "请先配置 Anthropic API Key"}` and route uses `JSONResponse` 400 when error == missing_api_key
- `PiRunner._run_pi`: pass `env={**os.environ, "ANTHROPIC_API_KEY": key}` when key provided
- Add optional `PiRunner.__init__(self, api_key: str | None = None)` or `run(..., env: dict | None=None)` — prefer constructor `api_key: str | None = None` stored and merged in `_run_pi`
- `_make_agent`: if not secrets configured, never called from start; when called: `LoopConfig(mock_pi=False)`, `PiRunner(api_key=store.get_key())`

- [ ] **Step 1: Settings API tests**

```python
# tests/unit/test_settings_api.py
from pathlib import Path
from fastapi.testclient import TestClient
import server.deps as deps
from agent_loop.secrets import SecretsStore
from server.app import app

def test_api_key_roundtrip(tmp_path: Path, monkeypatch):
    store = SecretsStore(tmp_path / "secrets.json")
    monkeypatch.setattr(deps, "get_secrets", lambda: store)
    c = TestClient(app)
    r = c.get("/api/v1/settings/api-key")
    assert r.status_code == 200
    assert r.json()["configured"] is False
    r2 = c.post("/api/v1/settings/api-key", json={"api_key": "sk-ant-secret-9999"})
    assert r2.status_code == 200
    body = r2.json()
    assert body["configured"] is True
    assert body["hint"] == "9999"
    assert "secret" not in body.get("hint", "")
    # GET must not leak full key
    g = c.get("/api/v1/settings/api-key").json()
    assert "sk-ant-secret" not in str(g)

def test_start_without_key_400(tmp_path: Path, monkeypatch):
    store = SecretsStore(tmp_path / "secrets.json")
    monkeypatch.setattr(deps, "get_secrets", lambda: store)
    # ensure controller uses same store gate
    c = TestClient(app)
    r = c.post("/api/v1/loop/start")
    assert r.status_code == 400
    assert r.json()["error"] == "missing_api_key"
```

Wire `get_secrets()` in deps returning process-wide `SecretsStore()`.

- [ ] **Step 2: Run — expect FAIL**

```bash
pytest -q tests/unit/test_settings_api.py
```

- [ ] **Step 3: Implement routes_settings + app include**

```python
# server/routes_settings.py
from fastapi import APIRouter
from pydantic import BaseModel, Field
from server.deps import get_secrets

router = APIRouter(prefix="/api/v1/settings", tags=["settings"])

class ApiKeyBody(BaseModel):
    api_key: str = Field(default="")

@router.get("/api-key")
def get_api_key():
    return get_secrets().public_status()

@router.post("/api-key")
def post_api_key(body: ApiKeyBody):
    get_secrets().set_key(body.api_key)
    return get_secrets().public_status()
```

In `server/app.py` include `settings_router` **before** StaticFiles mount.

- [ ] **Step 4: Start gate in deps + routes_loop**

```python
# deps LoopController.start / _start_locked:
def _start_locked(self) -> dict[str, str]:
    if not get_secrets().is_configured():
        return {
            "status": self._status,
            "error": "missing_api_key",
            "message": "请先配置 Anthropic API Key",
        }
    ...
```

```python
# routes_loop.py
from fastapi.responses import JSONResponse

@router.post("/start")
def loop_start():
    result = get_controller().start()
    if result.get("error") == "missing_api_key":
        return JSONResponse(status_code=400, content=result)
    return result
```

`_make_agent`:
```python
def _make_agent(self) -> AgentLoop:
    key = get_secrets().get_key()
    cfg = LoopConfig(cwd=Path(".").resolve(), mock_pi=False)
    runner = PiRunner(api_key=key)
    return AgentLoop(cfg, pi_runner=runner, bus=self.bus, store=self.store)
```

Check `AgentLoop.__init__` accepts `pi_runner` — if yes use it; else set attribute after construct per existing API.

- [ ] **Step 5: PiRunner env injection**

```python
class PiRunner:
    def __init__(self, api_key: str | None = None) -> None:
        self.api_key = api_key

    def _run_pi(...):
        return subprocess.run(
            ["pi", "-p"],
            input=prompt,
            cwd=cwd,
            env=self._env(),
            **self._pi_run_kwargs(timeout_s),
        )

    def _env(self) -> dict[str, str]:
        env = dict(os.environ)
        if self.api_key:
            env["ANTHROPIC_API_KEY"] = self.api_key
        return env
```

Update unit test to assert `env` contains key when monkeypatching `subprocess.run`.

- [ ] **Step 6: Full unit pass**

```bash
pytest -q tests/unit/test_secrets.py tests/unit/test_settings_api.py tests/unit/test_pi_runner.py tests/unit/test_loop_controller.py tests/unit/test_api_name.py
pytest -q
```

- [ ] **Step 7: Commit**

```bash
git add agent_loop/pi_runner.py server tests .gitignore
git commit -m "feat: settings API, start requires API key, pi env inject"
```

---

### Task 3: BenchController + API + SSE stages

**Files:**
- Create: `server/bench_runner.py`
- Create: `server/routes_bench.py`
- Create: `tests/unit/test_bench_api.py`
- Modify: `server/deps.py` (get_bench_controller, bus)
- Modify: `server/app.py`

**Interfaces:**
- `class BenchController:`
  - `__init__(self, bus: EventBus, data_path: Path | None = None, cwd: Path | None = None)`
  - `start() -> dict`  # ok or error busy
  - `status_payload() -> dict`  # status, stages list, report
- Stages list items: `{"stage": str, "message": str, "progress": float|None}`
- status: `idle|running|done|error`
- POST `/api/v1/bench/run` → 200 `{status:running}` or **409** if running
- GET `/api/v1/bench/status` → status_payload
- Worker: subprocess `sys.executable -m benchmarks.benchmark --data <path> --json` (no --limit)
- Publish: `bus.publish("bench.stage", {...})`, `bus.publish("bench.done", report_summary)`, optional `bus.publish("log.line", {"line": message})`
- Report summary extract dual/en/zh/fails/n from JSON (read actual keys from `benchmarks/benchmark.py` finalize — use whatever keys exist, e.g. nested `overall` or top-level; implement `_summarize(report) -> dict` with defensive gets)

- [ ] **Step 1: Write tests with injectable runner**

```python
# tests/unit/test_bench_api.py
from fastapi.testclient import TestClient
from server.app import app
import server.deps as deps

def test_bench_run_and_status(monkeypatch):
    # Replace controller with fast fake that completes immediately
    ...
```

Minimum cases:
1. POST run when idle → 200, status becomes running or done (fake)
2. POST run when running → 409
3. GET status shape has status/stages/report
4. Unit test BenchController emits stages via bus subscriber list

For speed: `BenchController` accepts optional `run_fn: Callable[[], dict]` used instead of subprocess in tests.

- [ ] **Step 2: Implement bench_runner.py**

Skeleton:

```python
class BenchController:
    def start(self) -> dict[str, str]:
        with self._lock:
            if self._status == "running":
                return {"status": "running", "error": "busy"}
            self._status = "running"
            self._stages = []
            self._report = None
            t = threading.Thread(target=self._worker, daemon=True)
            t.start()
            return {"status": "running"}

    def _emit(self, stage: str, message: str, progress=None):
        item = {"stage": stage, "message": message, "progress": progress}
        self._stages.append(item)
        self.bus.publish("bench.stage", item)
        self.bus.publish("log.line", {"line": message})

    def _worker(self):
        try:
            self._emit("load", "[bench] 加载数据集 ...")
            self._emit("score", "[bench] 计分中…")
            report = self._run_fn() if self._run_fn else self._subprocess_bench()
            self._emit("summarize", "[bench] 汇总结果…")
            summary = self._summarize(report)
            self._report = summary
            self._emit("done", f"[bench] 完成 dual=...")
            self.bus.publish("bench.done", summary)
            self._status = "done"
        except Exception as exc:
            self._emit("error", f"[bench] 失败: {exc}")
            self._status = "error"
```

Route:

```python
@router.post("/run")
def bench_run():
    result = get_bench().start()
    if result.get("error") == "busy":
        return JSONResponse(status_code=409, content=result)
    return result
```

- [ ] **Step 3: pytest**

```bash
pytest -q tests/unit/test_bench_api.py
pytest -q
```

- [ ] **Step 4: Commit**

```bash
git add server tests
git commit -m "feat: full-benchmark runner API with staged SSE events"
```

---

### Task 4: Web UI — Settings + Start gate + Bench toolbar

**Files:**
- Modify: `web/index.html`
- Modify: `web/css/app.css`
- Modify: `web/js/app.js`

**Interfaces (frontend):**
- `GET/POST /api/v1/settings/api-key`
- `POST /api/v1/loop/start` handles 400
- `POST /api/v1/bench/run`, `GET /api/v1/bench/status`
- SSE events: `bench.stage`, `bench.done`, `log.line`

- [ ] **Step 1: HTML**

Topbar: before controls, button `#btn-settings` (gear SVG) `aria-label="设置 API Key"`.

Modal/drawer `#settings-modal`:
- status text `#settings-key-status`
- input `#settings-key-input` type=password
- buttons save / clear / close
- help text about secrets.json

Bench panel top:
```html
<div class="bench-toolbar">
  <button id="btn-bench-run" class="btn btn-primary">跑全量测试</button>
  <span id="bench-run-status">空闲</span>
</div>
<div class="card">
  <h3>手动全量 KPI</h3>
  <dl id="bench-manual-kpi" class="kv"></dl>
</div>
<pre id="bench-stage-log" class="log-view" aria-live="polite"></pre>
<button id="btn-bench-log-clear" type="button">清空日志</button>
<!-- existing session before/after below -->
```

- [ ] **Step 2: JS**

- `loadKeyStatus()` on boot and after save
- Save/clear POST
- `onStartClick`: await key status; if !configured → `alert("请先配置 Anthropic API Key")` + open modal; return
- else existing git confirm then POST start; if 400 missing_api_key same alert
- Bench: click → POST run; on 409 alert 测试正在进行
- SSE handler: if type bench.stage append log; if bench.done fill KPI
- Poll GET `/api/v1/bench/status` every 2s while running (backup if SSE miss)
- Disable `#btn-bench-run` while running

- [ ] **Step 3: CSS**

Modal overlay, toolbar spacing, log mono max-height scroll; Dark OLED vars.

- [ ] **Step 4: Manual smoke**

```bash
source .venv/Scripts/activate
uvicorn server.app:app --host 127.0.0.1 --port 8765
# curl settings + bench status
curl -s http://127.0.0.1:8765/api/v1/settings/api-key
curl -s -X POST http://127.0.0.1:8765/api/v1/loop/start  # expect 400 without key
```

- [ ] **Step 5: Commit**

```bash
git add web
git commit -m "feat: console settings UI and full bench run toolbar"
```

---

### Task 5: README + 回归

**Files:**
- Modify: `README.md`

- [ ] **Step 1: README 补充**

在 Console 节增加：

```markdown
### API Key（Anthropic）

在 Console 顶栏「设置」中粘贴 `ANTHROPIC_API_KEY`，保存在本机
`agent_loop/memory/secrets.json`（已 gitignore）。未配置时无法 Start。
Start 将使用真实 pi，并把 Key 注入子进程环境变量。

### 手动全量 Bench

Bench 页「跑全量测试」对 `data/merged_benchmark.json` 全量计分，
阶段日志经 SSE 推送。全量可能较慢。
```

- [ ] **Step 2: 全量回归**

```bash
source .venv/Scripts/activate
pytest -q
python tools/structure_lint.py --root src/namepredict
```

- [ ] **Step 3: Commit**

```bash
git add README.md
git commit -m "docs: console API key and full bench quick notes"
```

---

## Spec Coverage Checklist

| Spec 项 | Task |
|---------|------|
| secrets.json + gitignore | 1 |
| public_status 无明文 | 1, 2 |
| GET/POST settings API | 2 |
| Start 无 Key 400 + 前端弹窗 | 2, 4 |
| 真 pi + env 注入 | 2 |
| 全量 bench API + 409 | 3 |
| SSE bench.stage/done + log.line | 3, 4 |
| Bench UI 工具条 + KPI + 日志 | 4 |
| README | 5 |

## Placeholder Scan

无 TBD；CLI 使用已有 `benchmarks.benchmark --json`。

## Type Consistency

- `public_status`: `configured: bool`, `hint: str|None`
- start error: `error: "missing_api_key"`
- bench busy: `error: "busy"` → HTTP 409
- SSE payload: `stage`, `message`, `progress`
