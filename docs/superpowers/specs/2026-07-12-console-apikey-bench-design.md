# ChemAgent Console：Anthropic API Key 设置 + 手动全量 Bench

**日期**: 2026-07-12  
**状态**: 待实现  
**范围**: 本机 Console（`127.0.0.1`）— Key 持久化与 Start 门禁；Bench 页全量测试与分阶段日志  

---

## 1. 背景与目标

### 1.1 问题

- 真 agent 循环依赖 `pi`，需 `ANTHROPIC_API_KEY`；当前 Console Start 硬编码 `mock_pi=True`，且无 Key 配置 UI。
- Bench 页仅展示 session 摘要，无法从网页手动触发全量 benchmark，也无分阶段进度。

### 1.2 目标

1. **API Key（仅 Anthropic）**  
   - 网页保存 / 更新 / 清除  
   - 落盘本机文件（gitignore），接口**永不返回明文**  
   - **未配置 Key → 不能 Start**（前端弹窗 + 后端 400）  
   - 有 Key → 真 `PiRunner`，子进程注入 `ANTHROPIC_API_KEY`

2. **Bench 页手动全量跑**  
   - 按钮跑完整 `data/merged_benchmark.json`（无 limit）  
   - 分阶段日志经 SSE 推送  
   - 结果刷新「手动全量」KPI（dual / en / zh / fails / n）

### 1.3 非目标

- 多 provider、OAuth、云端多用户鉴权  
- 修改 benchmark 金标或计分逻辑  
- 多并发 bench 队列（同时只允许 1 个手动 bench）  
- 将 Key 写入 git、session JSON、progress、SSE payload  
- 强制改 uvicorn 绑定逻辑（文档仍约定 `127.0.0.1`）

### 1.4 成功标准

- 无 Key 点 Start → 弹窗，不启循环  
- 有 Key 点 Start → 真 pi + env 注入  
- Bench「跑全量测试」→ 阶段 log 可见 → 结束后有分数摘要  
- 重启 uvicorn 后 Key 仍在（文件持久化）  
- 相关 pytest 通过；`secrets.json` 被 gitignore  

---

## 2. 已确认决策

| # | 决策 | 选择 |
|---|------|------|
| 1 | Provider | 仅 Anthropic（`ANTHROPIC_API_KEY`） |
| 2 | Key 存储 | 本机文件 `agent_loop/memory/secrets.json`，gitignore |
| 3 | Start 与 Key | **必须已配置 Key 才能 Start**；否则弹窗报错 |
| 4 | Start 与 pi | 有 Key → 真 pi（`mock_pi=False`），不再 Console 死 mock |
| 5 | 手动 Bench 范围 | **始终全量** `data/merged_benchmark.json` |
| 6 | 阶段日志 | SSE `bench.stage` / `bench.done`（可选同步 `log.line`） |
| 7 | 与 loop 并发 | 手动 bench 与 loop **允许并行**（只读计分；占 CPU） |
| 8 | 架构方案 | 最小增量：现有 FastAPI + EventBus + 静态 web |

---

## 3. 架构与数据流

```
[Settings UI] ──POST──► /api/v1/settings/api-key ──► secrets.json
[Settings UI] ──GET───► configured + 末4位 hint（无明文）

[Start] ──► GET configured? ──否──► 弹窗 + 打开 Settings
                │是
                ▼
         POST /api/v1/loop/start
                │ 后端再校验 secrets
                ▼
         AgentLoop + PiRunner(env=ANTHROPIC_API_KEY)

[Bench 跑全量] ──POST──► /api/v1/bench/run
                ▼
         BenchController（单飞）
                │ load → score → summarize
                ▼
         EventBus: bench.stage / bench.done
                ▼
         SSE /api/v1/events → Bench 阶段日志 + KPI
```

---

## 4. Key 存储

| 项 | 约定 |
|----|------|
| 路径 | `agent_loop/memory/secrets.json` |
| 格式 | `{"anthropic_api_key": "sk-ant-..."}`；缺省或 `{}` = 未配置 |
| 清除 | 写空对象或删除键；POST `api_key: ""` 等价清除 |
| Git | `.gitignore` 增加该路径（或 `agent_loop/memory/secrets.json`） |
| 读取失败 / JSON 损坏 | 视为未配置；不抛到客户端细节 |

模块建议：`agent_loop/secrets.py`

- `is_configured() -> bool`  
- `get_key() -> str | None`（仅服务端内部）  
- `set_key(key: str) -> None`  
- `clear_key() -> None`  
- `public_status() -> {configured, hint?}`（hint = 末 4 位或省略）

---

## 5. API

### 5.1 Settings

| 方法 | 路径 | 请求 | 响应 |
|------|------|------|------|
| GET | `/api/v1/settings/api-key` | — | `{ "configured": bool, "hint": "xxxx" \| null }` |
| POST | `/api/v1/settings/api-key` | `{ "api_key": string }` | `{ "configured": bool, "hint": ... }`；空串 = 清除 |

**禁止**：任何响应、日志、SSE 包含完整 Key。

### 5.2 Loop Start 门禁

- `POST /api/v1/loop/start`：若未配置 Key → **HTTP 400**  
  `{ "error": "missing_api_key", "message": "请先配置 Anthropic API Key" }`  
- 已配置：`LoopConfig.mock_pi=False`；`PiRunner` 启动子进程时  
  `env = {**os.environ, "ANTHROPIC_API_KEY": key}`  
- 不把 Key 写入 prompt 文件或 session JSON  

### 5.3 手动 Bench

| 方法 | 路径 | 行为 |
|------|------|------|
| POST | `/api/v1/bench/run` | 启动全量 bench 后台任务；已在跑 → **409** |
| GET | `/api/v1/bench/status` | `{ "status": "idle\|running\|done\|error", "stages": [...], "report": null\|摘要 }` |

**report 摘要字段**（与现有 benchmark 报告对齐，至少）：`dual`（或 `acc_dual`）、en/zh 准确率、`fails` 数、`n`（样本数）。实现时以 `benchmarks.benchmark` 实际 JSON 键为准，前端做兼容映射。

### 5.4 SSE（`/api/v1/events`）

| event | payload |
|-------|---------|
| `bench.stage` | `{ "stage": "load\|score\|summarize\|done\|error", "message": str, "progress": number\|null }` |
| `bench.done` | 与 status.report 同形摘要 |
| `log.line`（可选） | 与 stage message 相同文案，便于 Logs 页 |

---

## 6. Bench 运行器

- 组件：`BenchController`（或 `server/bench_runner.py`），进程内单例，与 `LoopController` 类似。  
- **全量**：`limit=None`，数据路径默认 `data/merged_benchmark.json`（与 `LoopConfig.data_path` 一致时可复用）。  
- **冷导入**：优先 subprocess 调 `python -m benchmarks.benchmark --json`（无 limit），与 loop 子进程 bench 策略一致；阶段：  
  1. `load` — 开始 / 已加载（若 CLI 不拆阶段，则用「开始跑全量」「子进程进行中」「解析 JSON」近似）  
  2. `score` — 计分中  
  3. `summarize` — 汇总  
  4. `done` / `error`  
- 若现有 CLI 难以拆细 progress，允许：`load` → `score`（running）→ `done`/`error` 三级，但 message 必须中文可读。  
- 互斥：`running` 时第二次 POST → 409。  
- 与 loop：**允许并行**。

### 阶段文案（中文）

| stage | message 示例 |
|-------|----------------|
| load | `[bench] 加载数据集 data/merged_benchmark.json …` |
| load 完成 | `[bench] 已加载 N 条`（若可知 N） |
| score | `[bench] 计分中…` |
| summarize | `[bench] 汇总结果…` |
| done | `[bench] 完成 dual=…% en=… zh=… fails=…` |
| error | `[bench] 失败: <原因>`（原因不得含 Key） |

---

## 7. UI

### 7.1 Settings（顶栏）

- Start 控件左侧：齿轮按钮，打开 Settings 浮层/抽屉（Dark OLED）。  
- 字段：状态（已配置 ···xxxx / 未配置）、password 输入、保存、清除。  
- 说明：仅存本机 `agent_loop/memory/secrets.json`，不进 git，不写日志。  

### 7.2 Start

1. 点击 → 先 GET api-key；未配置 → 弹窗「请先配置 Anthropic API Key」并打开 Settings；**不**调 start。  
2. 已配置 → 保留现有 git 门禁确认文案。  
3. 后端 400 `missing_api_key` → 同样弹窗。  

### 7.3 Bench 页

- 顶部工具条（无 session 也可见）：  
  - **跑全量测试**（running 时 disabled + spinner）  
  - 状态：idle / 运行中 / 完成 / 失败  
  - 阶段日志 mono 滚动区（可清空）  
  - KPI 卡：「手动全量」dual / en / zh / fails / n  
- 下方保留 session Before/After 摘要。  

### 7.4 a11y

- 按钮 `aria-label`；可见 focus；`prefers-reduced-motion` 沿用；无 emoji 图标。  

---

## 8. 错误处理

| 场景 | 行为 |
|------|------|
| 未配置 Key Start | 前端弹窗；后端 400 |
| Key 无效 / pi 鉴权失败 | 现有 cycle gate revert；session/SSE 记失败，无 Key |
| secrets 损坏 | 视为未配置 |
| bench 已在跑 | 409 |
| 数据集缺失 | stage error，提示 merge_datasets |
| 清除 Key 时 loop 仍在跑 | 允许清除；**下一次** `PiRunner.run` 读盘生效 |

---

## 9. 测试

1. secrets：set → public configured true 且无明文；clear → false  
2. API：POST/GET settings；无 Key start → 400  
3. 有 Key（tmp_path secrets）start 不因 missing_api_key 失败（可用 mock 线程）  
4. bench：run → status/done 路径；二次 run → 409；可用 mock BenchController  
5. PiRunner 调用时 env 含 ANTHROPIC_API_KEY（subprocess mock 断言）  

手测：无 Key Start 弹窗；保存 Key 后 Start 走真 pi 路径；Bench 全量阶段 log + KPI。  

---

## 10. 文件清单

| 路径 | 作用 |
|------|------|
| `agent_loop/secrets.py` | 读/写/清除/状态 |
| `agent_loop/pi_runner.py` | 子进程 env 注入 |
| `server/deps.py` | Start 校验；真 pi；读 secrets |
| `server/routes_settings.py` | api-key 路由 |
| `server/routes_bench.py` | run + status |
| `server/bench_runner.py` | 全量 bench + 阶段事件 |
| `server/app.py` | 注册路由 |
| `web/index.html`, `web/css/app.css`, `web/js/app.js` | Settings + Bench UI |
| `.gitignore` | secrets.json |
| `tests/unit/test_secrets.py`, `test_settings_api.py`, `test_bench_api.py`（名可微调） | 单测 |
| `README.md` | Console 配 Key；Bench 手动全量 |

### 实现顺序

1. secrets 模块 + API + 测试  
2. Start 门禁 + PiRunner env + 去掉 Console 死 mock  
3. BenchController + SSE 阶段 + API  
4. 前端 Settings + Bench UI  
5. README / 手测  

---

## 11. 与既有约束的关系

- 不改金标、不计分逻辑  
- 生产命名代码约束不变  
- Console 仍本机；Key 仅服务真 pi  
- 命名/Namer 本身不需要 Key  

---

## 12. 验收清单

- [ ] `.gitignore` 含 secrets  
- [ ] GET api-key 无明文  
- [ ] 无 Key Start 前端弹窗 + 后端 400  
- [ ] 有 Key Start 使用真 PiRunner + env  
- [ ] Bench 全量 + 阶段 log + 结束 KPI  
- [ ] 二次 bench 409  
- [ ] pytest 绿  
