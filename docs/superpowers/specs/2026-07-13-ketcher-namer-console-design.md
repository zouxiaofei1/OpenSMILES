# ChemAgent Console：Namer 页嵌入 Ketcher 画板

**日期**: 2026-07-13  
**状态**: 待实现  
**范围**: Console Namer tab — 本地自托管 Ketcher + 现有 `POST /api/v1/name`

---

## 1. 背景与目标

### 1.1 问题

- Console Namer 仅支持手输 / 粘贴 SMILES，无法通过结构式编辑器绘制分子再命名。
- 已有 `POST /api/v1/name` 与暗色 Console UI，缺的是画板与 SMILES 的桥接。

### 1.2 目标

1. 在现有 **Namer** tab 嵌入 **本地自托管** Ketcher standalone 画板。
2. 从画板导出 SMILES，调用现有 namepredict API，展示 en / zh / source / time。
3. 支持 **手动 Name** + 可选 **实时命名**（默认关，debounce）。
4. 保留文本 SMILES 路径与 Recent(20) 历史；vendor 缺失时文本命名仍可用。

### 1.3 非目标

- 独立新应用 / 多页面产品化
- npm / Vite 等前端构建链路
- 修改 `POST /api/v1/name` 契约或命名算法
- SMILES ↔ 画板全双向自动同步（输入框改字不自动回写画板）
- MOL / InChI 作为主交互路径
- 用户偏好持久化、批量画板、云端 Ketcher

### 1.4 成功标准

- 画乙醇结构 → Name → 得到双语名（如 ethanol / 乙醇）
- 粘贴 SMILES →「从 SMILES 载入」→ 画板显示对应结构
- 实时开关 on：结构变化约 0.7s 后自动命名；相同 SMILES 不重复请求
- 实时开关 off：改结构不发请求
- 无 `web/vendor/ketcher` 时画板降级提示，文本 Name 仍可用
- 后端路由默认零改动

---

## 2. 已确认决策

| # | 决策 | 选择 |
|---|------|------|
| 1 | 前端形态 | 接入现有 Console Namer 页（非独立应用） |
| 2 | Ketcher 引入 | 本地自托管静态包 → `web/vendor/ketcher/` |
| 3 | 命名触发 | 手动 Name + 可选实时（默认关） |
| 4 | 架构方案 | 双栏布局；iframe 同源加载 standalone；零构建 |
| 5 | 画板 ↔ 文本 | 画板 → 输入框（Name / 实时）；文本 → 画板仅经「从 SMILES 载入」 |
| 6 | 后端 | 复用 `POST /api/v1/name`；默认不改 server |

---

## 3. 架构与数据流

```
用户在 Ketcher 画分子
        │
        ▼
┌───────────────────┐   getSmiles()    ┌──────────────────┐
│ ketcher-standalone│ ───────────────► │ SMILES 输入框     │
│ (本地 iframe)     │                  │ (可手改)          │
└───────────────────┘                  └────────┬─────────┘
        ▲                                       │
        │ setMolecule(smiles)                   │ 手动 Name 或 实时 debounce
        │ （「从 SMILES 载入」）                  ▼
        └────────────────────────────  POST /api/v1/name
                                       { "smiles": "..." }
                                                ▼
                                       展示 en / zh / source / time
                                       + 写入 Recent(20)
```

### 3.1 单元边界

| 单元 | 职责 | 依赖 |
|------|------|------|
| `web/vendor/ketcher/` | 自托管 standalone 静态资源 | 无 |
| Namer 页布局 | 双栏：编辑器 + 命名面板 | 现有 `app.css` 设计 token |
| `KetcherBridge`（`app.js` 内或 `namer-ketcher.js`） | iframe 初始化、get/set SMILES、结构变化订阅 | Ketcher iframe API |
| 现有 namer 表单逻辑 | 校验、POST name、结果与 history | `POST /api/v1/name` |
| FastAPI 静态托管 | 提供 `/vendor/ketcher/**` | 现有 `web/` mount |

### 3.2 后端

- **默认无变更**。现有 `routes_name.py`：`POST /api/v1/name` body `{ "smiles": str }`。
- 仅当静态 MIME / 挂载无法服务 vendor 时再补最小配置。

---

## 4. UI 布局与交互

沿用 ChemAgent Console 暗色 dense dashboard（背景 `#0F172A`，CTA `#22C55E`）。

### 4.1 布局

**宽屏（≥ 960px）**

```
┌─ Namer tab ─────────────────────────────────────────────────┐
│  ┌─ Structure card (≈ 58%) ─┐  ┌─ SMILES Namer (≈ 42%) ───┐ │
│  │  [从 SMILES 载入] [清空]  │  │  [SMILES 输入] [Name]    │ │
│  │  ┌────────────────────┐  │  │  ☐ 实时命名              │ │
│  │  │  Ketcher iframe    │  │  │  error / result          │ │
│  │  │  min-height ~420px │  │  │  en · zh · source · ms   │ │
│  │  └────────────────────┘  │  └──────────────────────────┘ │
│  └──────────────────────────┘  ┌─ Recent (20) ────────────┐ │
│                                │  历史 · 清空              │ │
│                                └──────────────────────────┘ │
└──────────────────────────────────────────────────────────────┘
```

**窄屏（< 960px）**：上下堆叠，画板在上，命名 + 历史在下。

### 4.2 交互

| 动作 | 行为 |
|------|------|
| 点 **Name** | 优先 `ketcher.getSmiles()`；空/失败则用输入框；trim 后 POST；结果 + history |
| 勾选 **实时命名** | 结构变化后 **700ms debounce**；与 `lastNamedSmiles` 相同则跳过；in-flight 用 seq 防乱序 |
| **从 SMILES 载入** | 输入框内容写入 Ketcher；非法结构显示 error |
| **清空画板** | 清空 Ketcher；**不**自动清结果与 history |
| 手改 SMILES 后 Name | 直接 API；不强制回写画板 |
| 空结构 | 提示「请绘制或输入 SMILES」 |
| API 失败 | `namer-error`；实时模式不打断画板 |

### 4.3 前端状态

| 字段 | 含义 | 默认 |
|------|------|------|
| `liveNameEnabled` | 实时命名开关 | `false` |
| `ketcherReady` | iframe / API 可用 | `false` |
| `lastNamedSmiles` | 去重用 | `""` |
| `namerHistory` | 最近 20 条（现有） | `[]` |
| `nameReqSeq` | 请求序号，防乱序 | `0` |

### 4.4 无障碍

- 实时开关：label + checkbox，`aria-describedby` 说明 debounce
- 结果区保持 `aria-live="polite"`
- 文案中文，与 Console 一致

---

## 5. Ketcher 资源与集成

### 5.1 静态资源

| 项 | 约定 |
|----|------|
| 目录 | `web/vendor/ketcher/` |
| 入口 | iframe `src="/vendor/ketcher/index.html"`（以实际包入口为准） |
| 版本 | 固定版本写入 `web/vendor/ketcher/VERSION` |
| 来源 | 官方 ketcher-standalone 发布物（实现时锁定具体版本号） |
| Git | **默认提交 vendor**（满足离线 / 本地自托管）；若体积不可接受，可改为 `tools/fetch_ketcher.sh` + gitignore，但需在 README 写清安装步骤 |

### 5.2 初始化

- **Lazy load**：用户**首次**切换到 Namer tab 时再设置 iframe `src`（或取消占位），避免拖慢 Overview。
- 同源托管，直接使用 `iframe.contentWindow.ketcher`（或该版本文档规定的全局对象）。
- 封装 `KetcherBridge`：`init` / `getSmiles` / `setMolecule` / `clear` / `onChange`，业务不散落 `contentWindow`。
- 具体方法名以所选 standalone 版本文档为准，在实现计划中对照 `VERSION` 锁定。

### 5.3 降级

| 情况 | 处理 |
|------|------|
| vendor 缺失 / iframe 失败 | 画板区提示「Ketcher 未安装」及路径说明；文本命名可用 |
| `getSmiles` 失败 | 回退输入框；实时模式跳过本次 |
| Ketcher 未就绪 | 「从 SMILES 载入」禁用；Name 仅用输入框 |

---

## 6. 代码改动一览

| 路径 | 变更 |
|------|------|
| `web/index.html` | Namer 双栏 DOM、iframe、实时开关、载入/清空 |
| `web/css/app.css` | `.namer-layout` grid、iframe 高度、窄屏堆叠 |
| `web/js/app.js` 和/或 `web/js/namer-ketcher.js` | Bridge、onName 扩展、实时 debounce、tab lazy init |
| `web/vendor/ketcher/**` | standalone + `VERSION` |
| `server/*` | 默认无变更 |
| `README.md` | 一句：Namer 支持 Ketcher；vendor 位置/版本 |

---

## 7. 错误处理

| 场景 | 行为 |
|------|------|
| 空画板 + 空输入 | 校验失败，提示绘制或输入 |
| 非法 SMILES / 命名失败 | badge failed + `namer-error` |
| 实时请求失败 | 写 error，不重试风暴 |
| 并发响应乱序 | 只应用最新 `nameReqSeq` 的结果 |

---

## 8. 验收清单

1. `uvicorn server.app:app --host 127.0.0.1 --port 8765` → 打开 Console → Namer  
2. 绘制 CCO → Name → en/zh 正确（或与当前 namer 行为一致）  
3. 粘贴 SMILES → 从 SMILES 载入 → 画板结构正确  
4. 实时 on：改结构 → ~0.7s 更新；重复相同结构不重复有效请求  
5. 实时 off：改结构不请求  
6. 移除或缺失 vendor 时文本 Name 仍可用  

---

## 9. 测试策略

- **手工验收**为主（画板 iframe 不适合纯 pytest 无头全覆盖）。
- 可选：对 `KetcherBridge` 的纯逻辑（debounce 去重、seq 乱序丢弃、空 SMILES 校验）抽可测函数，用现有前端习惯或极轻量单测；**不强制**引入前端测试框架。
- 回归：现有文本 Name + history 行为不变。

---

## 10. 实现顺序建议

1. 取得并 vendoring Ketcher standalone → `web/vendor/ketcher/` + `VERSION`  
2. Namer DOM 双栏 + CSS  
3. `KetcherBridge` + lazy init  
4. 接线 Name（画板优先 SMILES）  
5. 「从 SMILES 载入」/ 清空  
6. 实时开关 + debounce + seq  
7. 降级文案 + README  
8. 手工验收清单  

---

## 11. 风险与备注

- **体积**：standalone 可能数 MB～数十 MB；进仓需确认团队可接受。备选：fetch 脚本 + 本地安装说明。  
- **API 漂移**：不同 Ketcher 大版本 iframe API 可能不同；以 `VERSION` 锁定并用 Bridge 隔离。  
- **主题**：Ketcher 自带浅色 UI，嵌在暗色 Console 内可接受；本版不做深度换肤。  
