# General Ring System & Locant Engine Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 把 NamePredict 的杂环/稠环/碳环不饱和从「按 kind 堆 `_try_*` + orienter + locant 表」推进到 **RingSystemIR → ScaffoldSpec → NumberingPlan → locant API** 的通用管线；使同类骨架扩展数据化，堵住环己二烯等环多烯掉进 methane 的缺口，并提供 **饱和单环任意杂原子数/种类的替换命名生成器（EN `xa…azacycloalkane` / ZH `x杂环y烷`）**。每阶段 dual 可验收且不破层纯度。

**Architecture:** 环相关职责拆三层——**(A) 拓扑 IR**（环是什么）、**(B) 骨架匹配 ScaffoldSpec**（叫什么、编号策略是什么）、**(C) 位次引擎 NumberingPlan**（原子 → 位次标签 + P-14 择优）。kind 字符串降级为 scoring/L5 显示标签；真正驱动编号的是 Spec + Plan。L2 parent 与 L3 -yl 共用同一 Plan。

**饱和单杂环命名优先级（硬规则）：**

1. **Retained / Hantzsch–Widman 保留名**（oxolane、piperidine、morpholine…）— 命中则用现有 stem，**禁止**再生成 `1-oxacyclopentane` / `1-氧杂环戊烷` 覆盖金标。
2. **替换命名生成器**（本计划 Phase 1R）— 未命中保留表时，由拓扑 `(ring_size, hetero_sites[])` **纯函数**生成 stem，**不**为每个 (N,O,S)×尺寸再注册 kind。
3. 不在本计划内做完整系统稠合名（furo[2,3-b]…）；那是 Phase 4 可选远期。

**Tech Stack:** Python 3.12、RDKit（SSSR / 芳香性 / 键型）、pytest、现有 `SMILESNNamer` 端到端 dual、`benchmarks.benchmark_parallel`。

## Global Constraints

- 代码只能落在 `src/namepredict/layer0`–`layer5` 或入口点；其他目录只读（`prompt.txt`）。
- 单文件 ≤500 行；每个函数 ≤10 行；函数式、尽量单输入/单输出。
- 层边界：
  - L1：拓扑事实 + IR + fingerprint（不写名称）
  - L2：匹配 scaffold、门控 simple/FG、产出 parent（可带 `numbering`）
  - L3：取代基；环 -yl 消费 Plan
  - L4：`choose_numbering` / 兼容旧 chain 定向；**不**按 kind 无限加表
  - L5：只拼 stem + locant 字符串；不发明化学规则
- 禁止深度学习；禁止大缓存；禁止改 benchmark 数据。
- 每任务结束跑相关 unit tests；触及命名输出时跑 `python -m benchmarks.benchmark_parallel --data data/merged_benchmark.json --time`；dual 回退不得 >0.5%。
- 每任务一个改进点；review PASS 后再 commit；`workstate.md` 按既有日志格式记一条。
- IUPAC 依据：`docs/iupac/`（P-14 位次、P-22/P-25 retained、**P-22.2.3 skeletal replacement（x杂环y烷）**、P-29 取代基、P-31 不饱和）。
- **禁止** 新增 `_INDOLE_ORIENT_KINDS.append(...)` 式扩展；新骨架只加 Spec 数据。
- **禁止** 向 `sat_hetero._MONO` 无限 append 表外种类；表外必须走替换 stem 生成器。
- **禁止** L5 为每个新骨架写独立 locant 函数；只消费 `plan` 的 label。
- **禁止** 用替换名覆盖 retained/Hantzsch–Widman 金标（oxolane≠1-oxacyclopentane）。
- 兼容：parent 同时可带 `chain` + `numbering`；旧路径读 `chain`，新路径读 `numbering`；迁完再删特判。

---

## 现状基线（实施前必读）

| 区域 | 文件 | 现状 | 缺口 |
|---|---|---|---|
| 拓扑 | `layer1/ring_systems.py` | SSSR 稠合连通、杂原子、芳香 | 未驱动 L2；fingerprint 过粗 |
| Retained 表 | `layer2/retained_registry.py` | Phase-0 骨架 4 条 | 未取代 `_try_*`；无布局指纹 |
| 5+6 引擎 | `layer2/fused56.py` | Mono/Di13 Spec + chain 序 | 仅一类；未升格为通用 Spec |
| 注册 | `layer2/ring_producers.py` + `kind_registry.py` | ring try + KindMeta | kind 枚举驱动全链路 |
| 位次 | `layer4/numbering.py` | kind→orienter；`_INDOLE_LOCANTS` / `_NAPH_LOCANTS` | `_orient_indole` 直接 `return chain`；桥头硬编码 |
| 环多烯 | `parent_selector._is_polyene` / `ring_parent._is_cycloalkene_core` | 开链多烯 / 环单烯 | **环多烯全无** → 环己二烯 → methane |
| 饱和单杂环 | `layer2/sat_hetero.py` `_MONO` 硬表 ~11 条 | oxolane/piperidine/morpholine… | **表外尺寸/杂原子组合**（如 thiane、azepane、1,4-oxathiane）无 stem 生成器 |
| 环侧基 | `layer2/heteroaryl_sub.py` | 主要 pyridinyl | 稠杂环 -yl 未通用 |
| 拼名 | `layer5/benzene_names.py` / `assembler.py` | 部分 `*carboxylic` 特判 | stem 与 locant 未统一走 Plan |

**已知失败例（实施前可复现）：**

| SMILES | 期望 | 现状 |
|---|---|---|
| `C1=CCC=CC1` 等环己二烯 | cyclohexa-1,3-diene / 1,4- | parent=`alkane` → methane |
| `C1=CCCCC1` | cyclohexene | ✓ cycloalkene |
| `C=CC=C` | buta-1,3-diene | ✓ polyene |
| `C1CCSCC1` | thiane / 硫杂环己烷（或 1-thiacyclohexane） | 无 `_MONO` 项 → 失败/甲烷类 |
| `C1CCCCNC1` | azepane / 氮杂环庚烷 | 无表项 |
| `C1CSCCO1` | 1,4-oxathiane / 1,4-氧硫杂环己烷 | 无表项 |

**根因摘要：**

1. `cycloalkene` 要求 `_ring_double_count == 1`
2. `polyene` 要求 `not has_ring`
3. 环原子不进开链 → fallback 单碳 alkane
4. `sat_hetero._MONO` 只覆盖有限 (size, z_counts)；**无** 替换命名（P-22.2.3 skeletal replacement / 中文 x杂环y烷）生成路径

---

## 目标态文件结构

```
src/namepredict/
  layer1/
    ring_systems.py          # 已有 → 输出可被 IR 消费的 system dict
    ring_ir.py               # 新建：RingSystemIR / RingComponent / FusionEdge
    ring_fingerprint.py      # 新建：fingerprint 字符串（布局级）

  layer2/
    scaffold/
      __init__.py
      specs.py               # ScaffoldSpec / NumberingPolicy / LocantSite 注册
      match.py               # IR → ScaffoldHit（atom_map 角色）
      sub_rules.py           # SubRules（cap / allow / fg_block）读 Spec
      builders/
        fused56.py           # 从 layer2/fused56.py 迁入角色解析（薄）
        monohetero.py        # 单环杂芳 + 饱和 retained 薄封装
        carbocycle.py        # 环烷/环烯/环多烯
        sat_hetero_repl.py   # 新建：饱和单环替换命名 builder（表外兜底）
      sat_hetero_stem.py     # 新建：纯函数 stem 生成（EN/ZH x杂环y烷）
      hetero_a_names.py      # 新建：a-前缀表 oxa/thia/aza… + 中文 氧杂/硫杂/氮杂…
    fused56.py               # 过渡期保留；逐步变 re-export / 删
    sat_hetero.py            # 过渡期：_MONO 保留优先；未命中委托 repl builder
    retained_registry.py     # 过渡期：与 specs 对齐或删除
    ring_producers.py        # 逐步改为「对 match 结果建 parent」
    kind_registry.py         # 仍服务 scoring；stem 可从 Spec 同步

  layer4/
    locants/
      __init__.py
      plan.py                # NumberingPlan + locant() / locant_int()
      generate.py            # 按 NumberingPolicy.mode 出候选
      constraints.py         # P-14 约束栈（主 FG / 不饱和 / 取代基）
      engine.py              # choose_numbering()
    numbering.py             # thin wrapper：ring 走 engine，chain 走旧逻辑

  layer3/
    substituent_extractor.py # 环 -yl 逐步调 scaffold match + plan
    # heteroaryl_sub 变薄或删除重复 walk

  layer5/
    assembler.py / benzene_names.py  # 消费 plan labels；削减特判表
    # sat_hetero repl stem 由 parent.stem_* 或 Spec 带入，assembler 不解析拓扑

tests/unit/
  test_ring_ir.py
  test_ring_fingerprint.py
  test_scaffold_match.py
  test_numbering_plan.py
  test_locant_engine.py
  test_cyclopolyene.py           # 环己二烯 1,3 / 1,4
  test_scaffold_fused56_specs.py # indole/btz/box 数据化回归
  test_sat_hetero_stem.py        # 纯函数 stem：任意 Z 组合 + 尺寸
  test_sat_hetero_repl.py        # 端到端：表外 thiane/azepane/oxathiane
```

**层边界复核：**

- L1 `RingSystemIR`：只描述拓扑，不匹配名称
- L2 `ScaffoldHit`：匹配 + 门控 + parent dict；**替换 stem 在 L2 生成并写入 parent**（L5 只拼）
- L4 `NumberingPlan`：编号权威；`locant(atom)` 全项目唯一入口
- L5：`stem` + `locant` 字符串；1H- / join 规则保留；**禁止** L5 数环原子/猜杂原子

---

## 核心数据模型（接口契约）

### RingSystemIR（L1）

```python
# 概念字段（实现用 dataclass；函数 ≤10 行拆分）
RingComponent: sssr_idx, size, atom_ids, hetero[(idx,Z)...], aromatic
FusionEdge: a, b, shared, bond
RingSystemIR: atom_ids, components, fusions, topology, fingerprint
```

### ScaffoldSpec（L2 注册）

```python
ScaffoldSpec:
  id, naming_class, stem_en, stem_zh
  n_rings, ring ("hetero"|"carbo"), retained, fg_rank
  numbering: NumberingPolicy
  sub_rules: SubRules | None
  principal_slots: ...   # 可选 *carboxylic / *amine 变体
```

### NumberingPolicy

```python
mode:
  fixed_roles       # retained 稠环标准 path
  fixed_hetero      # 杂原子=1，两方向（单杂 / 保留单杂环）
  multi_hetero      # 多杂原子：a-前缀位次集最低 + 元素优先级（O>S>Se>N>…）
  carbocycle_free   # 环烷：旋转×反转
  poly_unsat        # 环多烯：ene 位次组最低
  naph_family       # 萘/喹啉族标准链集
  # fusion_systematic  # Phase 4 可选

standard_path: tuple[LocantSite, ...]   # fixed_roles
anchors: tuple[str, ...]
substitutable: frozenset[str]           # 角色名；桥头默认不可取代
```

### 替换命名 stem 契约（L2 纯函数，无 mol）

```python
# hetero_site: (locant_int, atomic_num)  — locant 已按 multi_hetero 择优后的编号
# ring_size: int  (3..N，首期 3–12，可配置上限)

def sat_hetero_stem_en(ring_size: int, sites: list[tuple[int, int]]) -> str:
    """e.g. (6, [(1,16)]) → '1-thiacyclohexane'
       (6, [(1,8),(4,16)]) → '1,4-oxathiane' 若走 HW 特例表；
       否则 '1-oxa-4-thiacyclohexane'
    """

def sat_hetero_stem_zh(ring_size: int, sites: list[tuple[int, int]]) -> str:
    """e.g. (6, [(1,16)]) → '1-硫杂环己烷' 或 '硫杂环己烷'（单杂省略 1- 按中文惯例）
       (7, [(1,7)]) → '氮杂环庚烷'
       (6, [(1,8),(4,16)]) → '1,4-氧硫杂环己烷'
    """
```

**元素 a-前缀表（首期必含，可扩）：**

| Z | EN a-prefix | ZH |
|---|---|---|
| 8 O | oxa | 氧杂 |
| 16 S | thia | 硫杂 |
| 7 N | aza | 氮杂 |
| 15 P | phospha | 磷杂 |
| 34 Se | selena | 硒杂 |
| 5 B | bora | 硼杂 |
| 14 Si | sila | 硅杂 |

**环尺寸词干（与开链 alkane 对齐，cyclo 前缀）：**

| size | EN cyclo…ane | ZH 环…烷 |
|---|---|---|
| 3 | cyclopropane | 环丙烷 |
| 4 | cyclobutane | 环丁烷 |
| 5 | cyclopentane | 环戊烷 |
| 6 | cyclohexane | 环己烷 |
| 7 | cycloheptane | 环庚烷 |
| 8 | cyclooctane | 环辛烷 |
| … | 至至少 12；更大用 环十三烷… 与 mult 前缀表 | |

**编号硬规则（P-22.2.3 / P-14，与 `multi_hetero` 对齐）：**

1. 杂原子整体位次集最低。
2. 同集时：按元素优先级（O > S > Se > Te > N > P > … 与 IUPAC a 次序）使较高优先级元素得较低位次。
3. 再同：取代基位次集 / 字母序。
4. 单杂原子：默认位次 1，名称中 **EN 可写 `1-thia…` 或省略**；**ZH 单杂常省略「1-」**（`硫杂环己烷`），多杂保留全部位次。
5. **Retained 优先**：`_MONO` / retained Spec 命中则 **不** 调用替换 stem（避免 oxolane 变 1-oxacyclopentane 破坏 dual）。

### ScaffoldSpec 对替换母体

```python
# 动态 Spec（非枚举每个分子）：
ScaffoldSpec(
  id="sat_hetero_repl",          # 或 id 含 fingerprint
  naming_class="replacement",
  stem_en=None, stem_zh=None,    # 运行时由 sat_hetero_stem_* 填 parent
  n_rings=1, ring="hetero", retained=False, fg_rank=0,
  numbering=NumberingPolicy(mode="multi_hetero" | "fixed_hetero", ...),
  sub_rules=与 sat_hetero 对齐的 simple 侧链门控,
)
# parent dict 额外字段：
#   kind: "sat_hetero_repl" | 或 retained kind
#   stem_en / stem_zh: 生成结果
#   hetero_sites: ((locant, Z), ...)
#   ring_size: int
#   numbering: NumberingPlan
```

### LocantSite

```python
role, label ("1"|"3a"|...), atom_selector, bridgehead, hetero
```

### NumberingPlan（L4 权威）

```python
scaffold_id
atom_order: list[int]
labels: list[str]                 # 与 atom_order 对齐
atom_to_label / label_to_atom
sub_atoms: frozenset[int]
constraints_applied: tuple[str, ...]
```

### API

```python
locant(plan, atom) -> str | None
locant_int(plan, atom) -> int | None   # 比较用；"3a" → 约定全序
choose_numbering(ir, spec, atom_map, parent_ctx) -> NumberingPlan
```

### 约束栈（软比较，字典序）

1. 主官能团位次集（含 virtual COOH/amine）
2. 不饱和位次集（单烯 / 多烯 tuple）
3. 取代基位次集
4. 取代基字母序（alpha key）
5. 稳定排序（标准方向 / 较小原子编号）

硬过滤：拓扑一致、锚点角色位次、1H 固定等。

---

## 阶段总览

| Phase | 名称 | dual 意图 | 行为变化 |
|---|---|---|---|
| **0** | IR + Plan 接口，零行为 | 持平 | 适配 indole/naph 为 Plan，结果不变 |
| **1** | 单环统一 + 环多烯 | **主涨**（环己二烯等） | 新 kind `cyclopolyene` 或扩展 cycloalkene |
| **1R** | **x杂环y烷替换生成器** | **中涨**（表外饱和单杂环） | 任意 n 杂 + 尺寸 → EN/ZH stem；retained 不回归 |
| **2** | Retained 稠环数据化 | 持平/微涨，可扩展 | 删 `_INDOLE_LOCANTS` 等特判 |
| **3** | 侧链对称（ring -yl） | 中涨（杂芳侧基 + 饱和杂环 -yl） | heteroaryl_sub 通用化；repl -yl |
| **4** | 系统稠合名（可选） | 长尾 | 仅 retained 不够时 |

**执行原则：** 一次只做一个 Task；0 做完再开 1。Phase 1 是用户可感知的第一个正确性修复（环己二烯）。**Phase 1R 紧接 1.3/1.4 之后**：依赖 `multi_hetero`/`fixed_hetero` 编号与 carbocycle 词干表，但不依赖稠环 Phase 2。Phase 4 默认不做，除非 dual 长尾明确卡在系统名。

**Phase 1R 范围边界（写死，防膨胀）：**

| 在范围内 | 不在 1R 内（另开任务） |
|---|---|
| 单环、全饱和、环上杂原子 Z∈{O,S,N,P,Se,B,Si}（首期至少 O/S/N） | 芳香杂环（pyridine 等仍 retained） |
| 环尺寸 3–12（可配置；>12 可二期） | 稠合/螺/桥环 |
| 任意杂原子个数 1..⌊size/1⌋（门控：每原子价合理、无相邻违规可配置） | 不饱和杂环（oxolene 等） |
| 未取代 + 与现 `sat_hetero` 同级 simple 侧链（C 上 n-alkyl/halo） | N-烷基（现 sat_hetero 亦 out of scope） |
| Retained `_MONO` **优先** | 用替换名覆盖 oxolane/morpholine 金标 |
| EN skeletal replacement + ZH `x杂环y烷` | Hantzsch–Widman 新造词（除已有 retained） |

---

## Phase 0 — 接口与 IR（行为不变）

### Task 0.1: RingSystemIR 类型与构建

**Files:**
- Create: `src/namepredict/layer1/ring_ir.py`
- Create: `tests/unit/test_ring_ir.py`
- Modify: `src/namepredict/layer1/ring_systems.py`（可选：导出 IR 的薄封装，或 ring_ir 调 build_ring_systems）

- [ ] **Step 1: 写失败测试** — 苯、萘、吲哚 SMILES → IR 的 n_rings / topology / hetero 断言
- [ ] **Step 2: 实现 dataclass + `build_ring_ir(mol) -> list[RingSystemIR]`**
- [ ] **Step 3: 单测通过**；不改 namer 路径
- [ ] **Step 4: commit** — `feat(layer1): RingSystemIR from SSSR fusion graph`

### Task 0.2: fingerprint（布局级）

**Files:**
- Create: `src/namepredict/layer1/ring_fingerprint.py`
- Create: `tests/unit/test_ring_fingerprint.py`

- [ ] **Step 1: 测试** — indole vs benzofuran fingerprint 不同；benzene 稳定
- [ ] **Step 2: 实现** — `sizes + fusion + hetero layout + aromatic`（函数拆短）
- [ ] **Step 3: 单测通过 + commit** — `feat(layer1): ring layout fingerprint`

### Task 0.3: NumberingPlan + locant API

**Files:**
- Create: `src/namepredict/layer4/locants/plan.py`
- Create: `tests/unit/test_numbering_plan.py`

- [ ] **Step 1: 测试** — 给定 atom_order + labels，`locant` / `locant_int`；3a 全序
- [ ] **Step 2: 实现纯数据 Plan（无引擎）**
- [ ] **Step 3: commit** — `feat(layer4): NumberingPlan and locant API`

### Task 0.4: 适配现有稠环 chain → Plan（零行为）

**Files:**
- Create: `src/namepredict/layer4/locants/adapt.py`（从 parent.chain + kind 填 Plan）
- Modify: `src/namepredict/layer4/numbering.py` — `_sub_locant` 优先读 plan，否则旧逻辑
- Test: 现有 indole/benzothiazole/naphthalene unit 全绿；dual 持平

- [ ] **Step 1: 适配 indole 族 9 原子 + `_INDOLE_LOCANTS` 语义**
- [ ] **Step 2: 适配 naph 10 原子 + `_NAPH_LOCANTS` 语义**
- [ ] **Step 3: 相关 unit + 抽样 dual；回退则停**
- [ ] **Step 4: commit** — `refactor(layer4): adapt retained chains to NumberingPlan`

---

## Phase 1 — 单环统一 + 环多烯（主修复）

### Task 1.1: ScaffoldSpec 最小集 + carbocycle builder

**Files:**
- Create: `src/namepredict/layer2/scaffold/specs.py`
- Create: `src/namepredict/layer2/scaffold/builders/carbocycle.py`
- Create: `tests/unit/test_scaffold_carbocycle.py`

- [ ] **Step 1: 注册** `cycloalkane` / `cycloalkene` / `cyclopolyene` Spec（mode 对应 free / free+ene / poly_unsat）
- [ ] **Step 2: builder 从 IR 或 info 识别单碳环 + 双键数**
- [ ] **Step 3: 单测 core 识别；不接 namer 也可先测 builder**
- [ ] **Step 4: commit** — `feat(layer2): carbocycle ScaffoldSpecs`

### Task 1.2: locant generate + constraints（环）

**Files:**
- Create: `src/namepredict/layer4/locants/generate.py`
- Create: `src/namepredict/layer4/locants/constraints.py`
- Create: `src/namepredict/layer4/locants/engine.py`
- Create: `tests/unit/test_locant_engine.py`

- [ ] **Step 1: `carbocycle_free` 候选 = 现有 `_ring_candidates` 语义**
- [ ] **Step 2: `poly_unsat` 主键 = ene locant 有序 tuple，其次取代基**
- [ ] **Step 3: `choose_numbering` 单测：甲基环己烷最低位次；环己-1,3-二烯 vs 1,4**
- [ ] **Step 4: commit** — `feat(layer4): ring locant engine (free + poly_unsat)`

### Task 1.3: 接入 parent 选择 — 修复环己二烯

**Files:**
- Modify: `src/namepredict/layer2/parent_selector.py` 或新建 thin try 经 scaffold
- Modify: `src/namepredict/layer2/kind_registry.py` — `cyclopolyene` KindMeta
- Modify: `src/namepredict/layer2/candidates.py` / ring 路径
- Modify: `src/namepredict/layer5/assembler.py` — cyclo + diene 词干（复用 polyene diene 后缀）
- Create: `tests/unit/test_cyclopolyene.py`

**验收 SMILES（至少）：**

| SMILES | en（示意） |
|---|---|
| `C1=CC=CCC1` 或 `C1C=CC=CC1` | cyclohexa-1,3-diene |
| `C1=CCC=CC1` | cyclohexa-1,4-diene |
| `C1=CCCCC1` | cyclohexene（不回归） |
| `C=CC=C` | buta-1,3-diene（不回归） |

- [ ] **Step 1: 红测** — 二烯不得为 methane；成功且含 diene
- [ ] **Step 2: L2 产出 parent kind + double_bonds + chain/plan**
- [ ] **Step 3: L4 定向 + L5 `cyclohexa-1,3-diene` / 中文环己-1,3-二烯**
- [ ] **Step 4: unit + dual；回退 >0.5% 则修或收窄**
- [ ] **Step 5: commit** — `feat(namepredict): cyclopolyene parent (cyclohexadiene)`
- [ ] **Step 6: workstate 记一条**

### Task 1.4: 单环杂芳 fixed_hetero 迁入引擎（可选同 Phase）

**Files:**
- Modify: pyridine / furan 等 orient → `fixed_hetero` generate
- Test: 现有 pyridine unit 不回归

- [ ] 仅当 1.3 dual 稳定后做；可拆下期
- [ ] commit — `refactor(layer4): monohetero orient via locant engine`

---

## Phase 1R — 饱和单环 x杂环y烷替换生成器

> **用户可感知交付 #2（仅次于环己二烯）：** 表外饱和单杂环不再静默失败，任意（合理）杂原子数/种类/环尺寸 → EN 替换名 + ZH x杂环y烷。  
> **依赖：** Phase 0 Plan API；Phase 1.1–1.2 的 carbocycle 尺寸词干与 locant engine 更佳。**不依赖** Phase 2 稠环。  
> **与 1.4 关系：** 1.4 做 retained 杂芳 orient；1R 做饱和替换 stem。可并行，但 1R 的 `multi_hetero` 编号应复用 engine，禁止第二套旋转逻辑。

### Task 1R.1: a-前缀表 + 环尺寸词干 + 纯函数 stem

**Files:**
- Create: `src/namepredict/layer2/scaffold/hetero_a_names.py`
- Create: `src/namepredict/layer2/scaffold/sat_hetero_stem.py`
- Create: `tests/unit/test_sat_hetero_stem.py`

**纯函数输入/输出（无 RDKit / 无 mol）：**

```python
# sites: 已编号后的 (locant, Z)，locant 从 1 起，升序
sat_hetero_stem_en(6, [(1, 16)])           → "1-thiacyclohexane"  # 或项目统一省略规则
sat_hetero_stem_zh(6, [(1, 16)])           → "硫杂环己烷"
sat_hetero_stem_en(7, [(1, 7)])            → "1-azacycloheptane"
sat_hetero_stem_zh(7, [(1, 7)])            → "氮杂环庚烷"
sat_hetero_stem_en(6, [(1, 8), (4, 16)])   → "1-oxa-4-thiacyclohexane"
sat_hetero_stem_zh(6, [(1, 8), (4, 16)])   → "1,4-氧硫杂环己烷"  # 或 1-氧杂-4-硫杂环己烷（定一种）
sat_hetero_stem_en(5, [(1, 8), (3, 8)])    → "1,3-dioxacyclopentane"  # 仅当未走 retained
sat_hetero_stem_en(8, [(1, 7), (5, 7)])    → "1,5-diazacyclooctane"
sat_hetero_stem_zh(8, [(1, 7), (5, 7)])    → "1,5-二氮杂环辛烷"
```

**拼装规则（EN，IUPAC skeletal replacement 简化）：**

1. 各位次 + a-prefix 按位次升序连接：`1-oxa-4-thia`（同种多杂用 di/tri + 位次集：`1,4-dioxa`）。
2. 接 `cyclo` + 碳环烷词干（**总环原子数**，含杂原子）：hex 对应 6 元 → `cyclohexane` 去 `cyclo` 再拼？  
   **标准写法：** `1-thiacyclohexane` = `1-thia` + `cyclohexane`（替换名嵌在 cycloalkane 前）。
3. 同种杂原子倍数：`1,4-dioxacyclohexane`；异种：`1-oxa-4-thiacyclohexane`（每位次各写 a-prefix）。
4. 函数 ≤10 行拆分：`_a_prefix(z)` / `_mult_hetero_chunk(sites)` / `_cyclo_ane(size)` / join。

**拼装规则（ZH）：**

1. 位次（多杂必写；单杂默认省略「1-」）+ 杂名序列 + `环` + 尺寸词 + `烷`。
2. 同种：`1,4-二氧杂环己烷`；异种：优先紧凑 `1,4-氧硫杂环己烷`，若实现成本高可用 `1-氧杂-4-硫杂环己烷`（**单测锁定一种**）。
3. 尺寸词与 `stems` / cycloalkane 中文一致（丙/丁/戊/己/庚/辛/壬/癸/十一/十二）。

- [ ] **Step 1: 红测** — 上表至少 6 组 EN/ZH 断言；未知 Z → 明确失败（返回 None / raise 约定）
- [ ] **Step 2: 实现 a-前缀表 + size 词干 + stem 纯函数**
- [ ] **Step 3: 单测通过；无 namer 依赖**
- [ ] **Step 4: commit** — `feat(layer2): sat hetero replacement stem generator (x杂环y烷)`

### Task 1R.2: multi_hetero 编号候选

**Files:**
- Modify: `src/namepredict/layer4/locants/generate.py` — `multi_hetero` 模式
- Modify: `src/namepredict/layer4/locants/constraints.py` — 杂原子位次集 + 元素优先级
- Create/Modify: `tests/unit/test_locant_engine.py`

**候选生成：**

- 单环：旋转 × 反转（与 `carbocycle_free` 同骨架）。
- 软键：杂原子 locant 有序 tuple 最低；并列时较高优先级元素位次更低；再取代基。
- 单杂：收敛为 hetero=1（`fixed_hetero` 可视为 `multi_hetero` 退化）。

**验收：**

| 拓扑 | 期望编号要点 |
|---|---|
| 六元 1 个 S | S=1 |
| 六元 O 与 S 相对 1,4 | O=1, S=4（O>S） |
| 六元 两 N 相对 1,4 | 1,4-diaza；取代基再打破平局 |

- [ ] **Step 1: 单测** multi_hetero 择优
- [ ] **Step 2: 实现 generate + constraints**
- [ ] **Step 3: commit** — `feat(layer4): multi_hetero numbering mode`

### Task 1R.3: builder 接入 — retained 优先，表外走替换

**Files:**
- Create: `src/namepredict/layer2/scaffold/builders/sat_hetero_repl.py`
- Modify: `src/namepredict/layer2/sat_hetero.py` — `_kind_of` 未命中时不直接 None，或 ring_producers 增加 repl try
- Modify: `src/namepredict/layer2/kind_registry.py` — `sat_hetero_repl` KindMeta（fg_rank=0, ring=hetero, retained=False）
- Modify: `src/namepredict/layer2/candidates.py` / `ring_producers.py`
- Modify: `src/namepredict/layer5/assembler.py` — 消费 `parent["stem_en"]`/`stem_zh`（若无则走 kind 表）
- Create: `tests/unit/test_sat_hetero_repl.py`

**门控（与现 sat_hetero 对齐后略放宽尺寸/杂数）：**

1. 恰 1 个 SSSR 环；全非芳；环键全单键。
2. 环原子仅 C + 允许杂原子集合。
3. 价态合理（三配位 N 可有 H；O/S 二配位）。
4. 外侧 simple：无主 FG；C 上侧链与现 `_outside_ok` / alkyl-halo 规则一致；**N-烷基仍排除**（与现注释一致，二期再开）。
5. **若 `_MONO` / retained Spec 命中 → 走旧 kind，禁止 repl stem。**

**验收 SMILES（至少）：**

| SMILES | EN（示意） | ZH（示意） | 备注 |
|---|---|---|---|
| `C1CCSCC1` | 1-thiacyclohexane / thiane* | 硫杂环己烷 | 表外单 S |
| `C1CCCCNC1` | 1-azacycloheptane | 氮杂环庚烷 | 7 元 N |
| `C1CSCCO1` | 1-oxa-4-thiacyclohexane | 1,4-氧硫杂环己烷 | 异种双杂 |
| `C1COCCOC1` | 1,4-dioxacyclohexane | 1,4-二氧杂环己烷 | 若 dioxane retained 则**不得**走这条 |
| `C1CCOC1` | oxolane | 氧杂环戊烷 | **retained 不回归** |
| `C1CCNCC1` | piperidine | 哌啶 | **retained 不回归** |
| `C1COCCN1` | morpholine | 吗啉 | **retained 不回归** |
| `CC1CCSCC1` | n-methyl-1-thiacyclohexane 位次择优 | 甲基硫杂环己烷 | simple 侧链 |

\* 若日后把 thiane 收进 retained 表，则改走 retained；生成器仍作表外兜底。

- [ ] **Step 1: 红测** — 表外分子 dual/unit 成功；retained 三例名不变
- [ ] **Step 2: builder 从 IR/info 抽 ring_size + hetero atom_ids → Plan → sites → stem**
- [ ] **Step 3: parent 带 `kind=sat_hetero_repl`, `stem_en`, `stem_zh`, `numbering`**
- [ ] **Step 4: L5 优先 parent stem；unit + dual；回退 >0.5% 则收窄门控**
- [ ] **Step 5: commit** — `feat(namepredict): replacement names for sat monoheterocycles (x杂环y烷)`
- [ ] **Step 6: workstate 记一条**

### Task 1R.4: 与 monohetero Spec 对齐 / 收缩 _MONO 为数据（可选）

**Files:**
- Modify: `scaffold/builders/monohetero.py` 或 specs — retained sat 条目数据化
- `_MONO` 可保留为 retained 快速路径，或改为 Spec 列表

- [ ] retained 与 repl **同一 core 检测**，仅 stem 来源分叉
- [ ] commit — `refactor(layer2): sat_hetero retained vs replacement single core`

### Task 1R.5: 替换母体的 -yl（可并入 Phase 3）

- 附着点最低碳位次（杂原子优先规则不变）
- stem + `-yl` / `基`：`1-thiacyclohexan-3-yl` / `硫杂环己烷-3-基`
- 默认可延后到 Phase 3.1 通用 ring_yl，避免 1R 膨胀

---

## Phase 2 — Retained 稠环数据化

### Task 2.1: ScaffoldSpec 表达 fused56 角色

**Files:**
- Create: `src/namepredict/layer2/scaffold/builders/fused56.py`
- Modify: 逐步让 `benzothiazole.py` / `benzoxazole.py` / `benzofuran.py` 变薄到 Spec 引用
- Test: 现有 fused56 / btz / box unit

- [ ] **Step 1: 把 Fused56MonoSpec / Di13Spec 映射为 ScaffoldSpec + atom_map 角色**
- [ ] **Step 2: standard_path 写入 labels 1…3a…7a**
- [ ] **Step 3: 行为不变 dual**
- [ ] **Step 4: commit** — `refactor(layer2): fused56 as ScaffoldSpec builder`

### Task 2.2: indole / benzimidazole / quinoline 进 Spec

**Files:**
- Spec 数据 + match 角色（NH=1 等）
- Modify: `layer4/numbering.py` 删除 `_INDOLE_ORIENT_KINDS` 长列表依赖（改 scaffold_id / mode）
- Delete 路径：`_indole_sub_locant` 硬表在适配期后删除

- [ ] **Step 1: 每迁一个骨架一组回归测**
- [ ] **Step 2: `_orient_indole` 改为 engine 或 no-op（Plan 已唯一）**
- [ ] **Step 3: dual 持平**
- [ ] **Step 4: commit** — `refactor: retained fused scaffolds via NumberingPlan`

### Task 2.3: naph_family 模式

**Files:**
- `generate.py` naph 候选 = 现有 `naph_chains` 语义
- 删除 `_NAPH_LOCANTS` 重复

- [ ] unit + dual；commit — `refactor(layer4): naph_family numbering mode`

### Task 2.4: retained_registry 与 specs 对齐或废弃

- [ ] 单一注册源：`scaffold/specs.py`
- [ ] `kind_registry` stem 可从 Spec 同步 bootstrap
- [ ] commit — `refactor(layer2): single scaffold registry`

---

## Phase 3 — 侧链对称（ring -yl）

### Task 3.1: 通用 ring_yl_name

**Files:**
- Create: `src/namepredict/layer2/scaffold/yl.py` 或 `layer3/ring_yl.py`
- Modify: `substituent_extractor` 调用
- Test: pyridinyl 回归 + 1–2 个稠环 -yl（若 Spec 已有）

- [ ] 附着点优先最低位次（P-29）
- [ ] leaves 仍走 leaves registry，挂在 plan.sub_atoms
- [ ] commit — `feat: generic ring-yl from NumberingPlan`

### Task 3.2: 收缩 heteroaryl_sub 手写 walk

- [ ] 删除与 Plan 重复的 pyridine walk
- [ ] dual；commit — `refactor: heteroaryl_sub via scaffold plan`

---

## Phase 4 — 系统稠合名（可选，默认不做）

仅当 dual 长尾明确需要：

- fusion_systematic mode：基础组分 + 融合边指示
- 三环以上 fingerprint
- 单独计划文档拆出，不在本文件扩写任务细节

---

## 设计原则（防再次腐化）

1. 新骨架 = 新 Spec 数据，不是新 orienter 函数。
2. 角色 `atom_map` 在匹配期完成；编号期不重新猜「谁是 N」。
3. bridgehead 是 **label 属性**，不是 kind 属性。
4. Parent 与 -yl 共用 Plan；差别只在 constraint 优先级。
5. 编号失败 / 无 plan → 候选丢弃或 `n_unhandled`，**禁止** 静默 methane。
6. 复杂分子优先扩大 Spec 覆盖 + 放宽 `sub_rules`，不在 assembler 加 if。
7. `fused56` 只做角色解析 builder，不做全项目第二套引擎。
8. **饱和单杂环：retained 表优先；表外只走替换 stem 生成器，禁止再往 `_MONO` 无限 append。**
9. **替换 stem 是 L2 纯函数（size + sites → EN/ZH）；L5 禁止数环/认杂原子。**
10. **禁止** 用 `1-oxacyclopentane` 覆盖 oxolane 等 dual 金标保留名。

---

## 与旧路径对照（迁移检查表）

| 旧 | 新 |
|---|---|
| `_is_cycloalkene_core` count==1 | Spec `cycloalkene` / `cyclopolyene` by n_db |
| `_is_polyene` and not has_ring | 开链仍 polyene；环走 cyclopolyene |
| `_chain_atoms` in fused56 | atom_map + standard_path → Plan |
| `_orient_indole` return chain | fixed_roles 唯一/镜像候选 |
| `_INDOLE_LOCANTS` / `_NAPH_LOCANTS` | `plan.labels` |
| `_sub_locant(kind)` | `locant(plan, atom)` |
| `_INDOLE_ORIENT_KINDS` | `spec.numbering.mode` / scaffold_id |
| `heteroaryl_sub` walk | ring_yl + Plan |
| kind 特判 L5 carboxylic | `stem` + `plan.locant(attach)` |
| `sat_hetero._MONO` 未命中 → None | retained miss → `sat_hetero_repl` + `sat_hetero_stem_*` |
| 每新杂环加 kind + stems 表 | a-前缀表 + size 词干 + 位次 Plan（生成器） |
| 手写 thiane/azepane try | 拓扑 sites → 自动 stem |

---

## 验收总标准

| 项 | 标准 |
|---|---|
| 环己二烯 | 1,3 / 1,4 正确命名；非 methane |
| 环己烯 / 开链二烯 | 不回归 |
| 现有 fused56 族 | unit + dual 不回退 >0.5% |
| 扩展新 5+6 | 理论上只加 Spec + 薄注册，不改 engine |
| **表外饱和单杂环** | thiane / azepane / 1,4-oxathiane 等 EN+ZH 正确；非 methane |
| **retained 饱和杂环** | oxolane / piperidine / morpholine / dioxane… **名与 dual 不回归** |
| **生成器任意性** | 同函数覆盖 1..k 杂原子、O/S/N（+可扩 Z）、环尺寸 3–12，**不**为每个组合加 kind |
| 层纯度 | L1 无名称；L5 无拓扑猜测；stem 在 L2 |
| 代码度量 | 单文件 ≤500；函数 ≤10 行 |

---

## 建议执行顺序（agent）

```
0.1 → 0.2 → 0.3 → 0.4 → dual 持平
  → 1.1 → 1.2 → 1.3（环己二烯）→ dual 涨分
  → 1.4（可选，杂芳 fixed_hetero）
  → 1R.1（stem 纯函数）→ 1R.2（multi_hetero）→ 1R.3（接入 namer，x杂环y烷）
  → 1R.4（可选对齐）
  → 2.x 按骨架逐个迁
  → 3.x -yl（含 sat_hetero_repl -yl）
  → 4 仅按需
```

**第一个用户可感知交付：Task 1.3（环己二烯）。**  
**第二个用户可感知交付：Task 1R.3（任意 x杂环y烷 / 替换饱和单杂环）。**  
Phase 0 不得跳过：没有 Plan API 就继续堆 kind 表。  
Phase 1R 不得把 retained 金标改成替换名。

---

## 参考（对话结论摘要）

- 杂环/稠环与 L3–L5 经 kind **紧耦合**；同类 5+6 已有 fused56 半数据驱动，真正复杂分子需要 IR+Spec+Plan。
- 环己二烯失败是 **功能缺口**（环单烯 ∩ 开链多烯 = 空），不是解析 bug。
- **x杂环y烷：** 现 `sat_hetero._MONO` 是有限 retained 表；用户要求 **任意数量/种类** → 必须是 **拓扑驱动的替换命名生成器**（EN skeletal replacement + ZH x杂环y烷），不是再 append `_MONO`。
- Retained（oxolane/哌啶/吗啉…）与替换生成器 **分叉 stem、共用 core 检测与 Plan**；金标优先保留名。
- 一句话蓝图：

```text
ring_systems(IR) → fingerprint match → ScaffoldSpec
  → roles/atom_map → generate NumberingPlan candidates
  → P-14 constraint stack → locant(atom)
  → L2 parent / L3 -yl / L4 / L5

饱和单杂环旁路：
  IR mono sat hetero → retained Spec? 
    yes → 旧 stem（oxolane…）
    no  → multi_hetero Plan → sites → sat_hetero_stem_en/zh（x杂环y烷）
```
