# Page Override: ChemAgent Console

> Overrides Master for the main operator console.  
> **Route:** `/`  
> **Purpose:** 控制 Agent 自循环、查看每轮上下文/session、在线测试 SMILES 命名、监控 benchmark 分数。

---

## Layout (Desktop-first, dense)

```
┌─ Top Bar (48px) ──────────────────────────────────────────────────────────┐
│ Logo · Status pill · dual% KPI · iter/K · [Start] [Pause] [Stop] · settings│
├─ Left Sidebar (280px) ──────────┬─ Main ──────────────────────────────────┤
│ Sessions                         │                                          │
│  · Live / History list           │  Tab: Overview | Cycle | Prompt | Bench  │
│  · Filter: running/done/failed   │       | Namer | Logs                     │
│  · Click → load cycle context    │                                          │
│                                  │  Content area (scroll)                   │
│ Fail clusters (compact)          │                                          │
│  · top failing features          │                                          │
└──────────────────────────────────┴──────────────────────────────────────────┘
```

- Desktop max content: full viewport width (ops dashboard, not 1200px marketing cap)
- Mobile (<768px): sidebar → drawer；Start/Pause/Stop 固定底栏
- Spacing: Master dense scale (`--space-*`)

---

## Primary Zones

### 1. Left: Session 列表
- 每项：`#iter` · 状态色点 · dual Δ · 主规则 ID · 相对时间
- 选中：左侧 accent bar（`#22C55E`）+ muted 背景
- >50 项虚拟列表；长文本 `truncate` + title tooltip
- Empty：`暂无循环 · 点击 Start 开始` + 次要说明

### 2. Main Tabs

| Tab | 内容 |
|-----|------|
| **Overview** | KPI：dual / en / zh / fails；dual 趋势折线；状态机当前步 |
| **Cycle** | 失败簇、主规则、Layer、TDD 红→绿、lint/pytest/bench、commit/revert |
| **Prompt** | 注入 pi 的完整 prompt（只读 + 复制）；skill 引用 |
| **Bench** | 分数表；失败样本（SMILES/gold/pred）；source·tier 分解 |
| **Namer** | SMILES → 本地 namepredict → en/zh/meta/耗时；最近 20 条 |
| **Logs** | 编排器 + pi 日志流（mono）；可暂停自动滚底 |

### 3. Controls
- **Start**：主 CTA（accent `#22C55E`）
- **Pause**：本轮结束后暂停（不杀进程）
- **Stop**：写 STOP + 可选强制 abort（confirm）
- 危险操作：destructive `#EF4444` + 确认对话框

---

## Status Semantics（色 + 文案 + 图标，不只靠颜色）

| 状态 | 色 | 文案 |
|------|----|------|
| idle | muted | 空闲 |
| running | accent + 脉冲点 | 运行中 |
| paused | amber | 已暂停 |
| gate_pass | green | 已提交 |
| gate_fail | red | 已回退 |
| error | red | 错误 |

---

## Charts
- dual 准确率：**折线**（时间 = 轮次）；单系列即可
- KPI 大数字始终可见（a11y / reduced-motion 时图表可静）
- 失败 feature：**横向条形** top-N（≤15）

---

## Real-time
- SSE 或 WebSocket：`loop.state` / `log.line` / `cycle.updated` / `bench.done`
- 断线顶栏警告 + 自动重连
- >300ms 操作：skeleton 或按钮 spinner；防连点 disable

---

## Namer 区
- 可见 label「SMILES」+ 输入 + 命名按钮（Enter 提交）
- 结果卡：en/zh mono；success badge；`time_ms`
- 错误在字段下方；提供重试
- 仅调用 API，不写缓存

---

## A11y / UX
- 对比度 ≥ 4.5:1；可见 focus ring
- 触控 ≥ 44px；icon 按钮 `aria-label`
- `prefers-reduced-motion`：关脉冲与图表入场动画
- 无 emoji 图标 → Lucide SVG
- `cursor-pointer` + 150–300ms hover

---

## Tech default
- FastAPI + `agent_loop` 同进程事件总线
- Frontend：`web/`（Tailwind dense dark；可先静态 SPA）
- API：`/api/v1/loop/*`, `/api/v1/sessions/*`, `/api/v1/name`, `/api/v1/events`

---

## Anti-patterns
- 默认浅色主题
- 侧栏塞完整 IUPAC 原文
- 浏览器内跑 RDKit
- 无确认的强制 Stop/reset
- 仅用颜色表示 gate 结果
