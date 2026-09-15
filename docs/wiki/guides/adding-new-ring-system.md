# 新增环系指南

本文档分步说明如何为 NamePredict 流水线新增一个环系（ring system）——新的单杂环、fused heterocycle（稠杂环）、稠合碳环，或新的单环烃稠合零件。

---

## 前置知识

阅读本指南前，建议先了解以下架构文档：

- [[architecture/overview]] — 六层流水线总览
- [[architecture/layer1-analyzer]] — Layer1 环系事实（`info["ring_systems"]`）的产出
- [[architecture/layer2-parent-selector]] — Layer2 母体选择器（骨架枚举 + P-44 筛选 + 稠环拆解）
- [[architecture/layer4-numbering]] — Layer4 编号引擎（P-14.4 三层分派 + 稠环编号）
- [[architecture/layer5-name-assembly]] — Layer5 中英双语名称组装
- [[concepts/functional-group-priority]] — FG 优先级规则
- [[reference/core-data-contracts]] — info / parent / numbered dict 字段表

---

## 流程概览

新增一个环系的改动按顺序为：

1. **Layer1 环系检测与划分**（通常无需改动——SSSR + 稠合边连通分量自动给出环系）
2. **Layer2 骨架身份与注册**（**主要工作量**）：在 `ring_scaffold._TEMPLATES` 注册 SMILES 模板；作稠合零件时补 `fused`/`fused_prefix`/`fused_stem`；保留编号时在同条目加 `standard`；新 `naming_class` 还要在 `ring_expression_policy._POLICIES` 放行环内/环外主 FG
3. **Layer4 编号**（**通常无需改动**——`numbering_engine.orient_numbering` 的 P-14.4 管线自动适用，`standard` 登记项自动走固定编号）
4. **Layer5 名称组装**（通常无需改动——词干由 `ScaffoldSpec` 经 `KindMeta` 注入；稠合名由 `fused_namer.fused_parent_names` 现场组装）
5. **FG-环组合支持**（环 + 主 FG 的正交化路径，通常只需第 2 步的 `_POLICIES` 放行）

> 核心原则：**scaffold × kind × 数量全正交**——加新环只加 `_TEMPLATES` 模板 + 可选 `standard`，不枚举 cycloalcohol/naphthalenol/indolol 组合 kind。

---

## 第一步: 环系检测与划分（Layer1）

> **源:** `src/namepredict/layer1/ring_systems.py`（124 行）

环系产出入口是 `build_ring_systems`（`ring_systems.py:117`），流程为：

1. `_sssr`（`:19`）取 RDKit `GetRingInfo().AtomRings()` 的全部最小环，经 `memo.by_mol` 缓存；`sssr_rings`（`:23`）是供各层统一调用的公开访问器（与 `_sssr` 同一次记忆）
2. `_ring_pairs`（`:27`）单遍扫全部环对：共享 **≥2** 原子记为**稠合边**（`fused`），恰好 **1** 个原子记为**螺环对**（`spiro`）
3. `_components`（`:54`）只用**稠合边**做并查集连通分量，每个分量即一个环系
4. `_system_entry`（`:104`）→ `_system_dict`（`:89`）组装事实 dict

由此得到的划分语义：**桥环**（两环共享 ≥3 原子，是 ≥2 的子集）天然并入同一分量；**螺环**（共享 1 原子）**不合并**，各自成分量——`build_ring_systems` 不消费 `spiro` 列表。

环系 dict 的键（`_system_dict`，`:93`）：`atom_ids` / `sssr_indices` / `fusion_edges` / `n_rings` / `n_atoms` / `hetero_atoms`（`:64`，每项 `{"idx", "Z"}`）/ `is_aromatic_mancude`（`None`）/ `topology`（`None`）。

**写盘点**：`analyzer._ring_meta`（`analyzer.py:121`）写入 `info["ring_systems"]`。

**消费点**：`parent_skeleton._ring_candidates`（`parent_skeleton.py:101`，只读 `atom_ids` 构造骨架候选，此阶段**不识别 scaffold 身份**）、`ring_scaffold._generic_carbocycle`（`:443`）、`numbering_engine._fused_numbering`（`:275`）、`fused_namer.fused_parent_names`（`:122`）。

另有跨层共用的两个访问器：`sssr_rings`（`:23`）与 `kekulized`（`:10`，返回 Kekulé 化并清芳香标志的副本，**失败返回 `None`**，调用方须退回原 mol）。下游调用者：`parent_skeleton.py:15`、`ring_scaffold.py:11`、`principal_expression.py:15`、`fused_system.py:15`、`numbering_engine.py:9`、`indicated_hydrogen.py:6`、`fused_namer.py:7`。

大多数新环系在 Layer1 无需任何改动。只有当环系划分本身要突破「稠合边连通分量」时才需要扩展本文件。

---

## 第二步: 骨架身份与注册（Layer2）—— 主要工作量

> **源:** `src/namepredict/layer2/ring_scaffold.py`（457 行）、`ring_expression_policy.py`（45 行）、`fused_system.py`（223 行）

### 2.1 身份识别入口：两步

`ParentSkeleton` 只承载 `topology` / `atom_ids` / `covered_principal_ids`（`parent_skeleton.py:25`），**不含 `scaffold_id`**——骨架身份**延迟到表达阶段**识别。识别入口是 `resolve_ring_scaffold`（`ring_scaffold.py:450`），只有两步：

```
resolve_ring_scaffold(info, skeleton)
  ├─ 1. match_retained(info, skeleton.atom_ids)   # SMILES 模板子图同构精确覆盖 → sid
  │      └─ get_spec(sid).identity                 # sid → ScaffoldIdentity
  └─ 2. _generic_carbocycle(info, skeleton)        # 兜底身份
```

- `match_retained`（`:387`）在 `_match_with_map`（`:393`）之上取 sid：① 按元素签名（`_elem_sig`，`:316` / `_TEMPLATE_ELEM`，`:321`）预筛，② 用 `_Q`（`:193`）的模板做元素标注子图同构并要求原子集**精确相等**；③ 精确匹配全失败后，按**完全氢化骨架**（`_Q_H`，`:235`）再比对一遍（P-25.3.4，加氢衍生物识别）。`mancude_only=True` 时过滤掉 `fused=False` 的饱和保留名
- `_generic_carbocycle`（`:437`）：全碳 → `ScaffoldIdentity("carbocycle", "carbocycle", 1, "carbo")`；非全碳且环数 ≥2 → `ScaffoldIdentity("fused_hetero", ...)`；**非全碳单环 → `None`（显式失败）**，即未注册的单杂环无法命名，必须注册模板

支持 API：`get_spec`（`:366`）、`all_specs`（`:371`）、`numbering_scaffold_facts`（`:376`）、`standard_chain`（`:428`）、`locant_prefix`（`:420`）。

### 2.2 注册保留母体：唯一事实来源 `_TEMPLATES`

**新增保留环系只改 `ring_scaffold.py` 的 `_TEMPLATES` 一张表**（`:60`，当前 **83 条**）。每条是 `sid → dict`：

```python
# ring_scaffold.py:71（`_TEMPLATES` 内 pyridine 条目）
"pyridine":    {"smiles": "n1ccccc1",   "stem_en": "pyridine",    "stem_zh": "吡啶",   "naming_class": "monohetero", "fused": True, "fused_prefix": ("pyrido", "吡啶并")},
```

字段清单：

| 字段 | 必填 | 作用 |
|---|---|---|
| `smiles` | 是 | 模板查询分子；`_Q`（`:193`）import 期由此构建，`_spec_from_template`（`:343`）由它算 `n_rings`/`ring` |
| `stem_en` / `stem_zh` | 是 | 双语词干，经 `ScaffoldSpec` → `KindMeta` 注入 L5 |
| `naming_class` | 是 | 命名类，`_POLICIES`（`ring_expression_policy.py:18`）按它放行 FG 表达 |
| `fused` | 否 | `True` 才可作稠合命名的**母体组分**（`component_stem`:148、`match_fusion_component`:189 的 `mancude_only` 过滤） |
| `fused_prefix` | 否 | 作**附加组分**时的保留稠合前缀 `(en, zh)`，如 `("benzo","苯并")` / `("pyrido","吡啶并")` |
| `fused_stem` | 否 | 稠合前取用的词干（去指示氢），如 `indole`→`indolo`；`component_stem`（`:148`）优先取它 |
| `locant_prefix` | 否 | 位次/指示氢前缀（`"1,3-"`、`"1H-"`、`"10H-"` 等） |
| `prefix_nh_conditional` | 否 | `True` 时仅当环内存在未取代芳香 NH 才注入 `locant_prefix`（`kind_registry._ring_keeps_nh_prefix`，`:42`） |
| `standard` | 否 | 固定编号 `(labels, order)` 二元组，见 2.3 |

当前 83 条的分布（`naming_class` → 条数）：`monohetero` 53、`naph_family` 9、`fused56` 8、`xanthene` 2、`mono_carbo`/`anthra`/`phenanthrene`/`pyrene`/`chrysene`/`carbazole`/`acridine`/`phenothiazine`/`benzodioxole`/`purine`/`steroid` 各 1。

可选字段统计：`fused=True` 67 条、`fused_prefix` 27 条、`fused_stem` 3 条（`indene`/`indole`/`purine`，`:66`/`:107`/`:126`）、`locant_prefix` 41 条、`prefix_nh_conditional` 14 条、`standard` 36 条。

词干本身的形态约定：

- `[1,2]oxazolo` / `[1,3]oxazolo` / `[1,3]thiazolo` / `[1,2,4]triazolo` / `[1,3,5]triazino` 等**方括号必须手写**——表征结构的位次要按 P-25.3.1.3 放进方括号，通用的「去尾 e 加 o」不含方括号（`_component_prefix`，`fused_namer.py:11`），会拼出 `1,2,4-triazolo`
- 词干把加氢前缀嵌在词中时（如 `4,5-dihydro-1,3-thiazole` 的 `1,3-`，`:106`），`locant_prefix` 与词干内前缀重复，靠 `pack_parent_stem`（`kind_registry.py:58`）的 `_strip_locant_prefix`（`:53`）避免位次重复
- **表序即优先级**：`_match_with_map` 按 `_Q` 的 dict 序（即 `_TEMPLATES` 的书写序）首次命中即返回，而元素签名相同的条目会撞车。实例：6 元含氧 mancude 母体 `pyran`（`:75`）**必须排在** `oxane`（`:89`）**之前**，否则环内 C=C 被静默丢弃（`pyran-2-one` → `oxan-2-one`）。pyran 与 oxane 的元素签名都是 `{C:5, O:1}`

`_spec_from_template`（`:343`）由条目派生 `ScaffoldSpec`（`retained=True`，`n_rings`/`ring` 从 smiles 算，`numbering.standard_path` 取 `standard[0]`），`_ALL_SPECS`（`:360`）与 `_BY_ID`（`:363`）是它的派生视图。`kind_registry._load_from_scaffold_specs`（`kind_registry.py:76`）据此自动注册 `KindMeta` 词干——**`ScaffoldSpec` 是词干的唯一权威**，新条目无需在 L5 另写词表。

位置异构体（quinoline/isoquinoline、二嗪、二唑等）在元素标注的子图同构下天然区分，无需额外消解。

### 2.3 固定编号：条目内的 `standard = (labels, order)`

若新环的保留编号不是「按模板原子序从 1 起」，在**同一条目内**加 `standard` 二元组：

- `order`：模板原子按 locant 顺序的下标排列
- `labels`：对应的 locant 标签元组
- **两者长度均等于模板原子数**

`_STANDARD_LABELS`（`:195`）/ `_STANDARD_ORDERS`（`:198`）是它的派生视图（当前 36 条登记）。`standard_chain`（`:428`）据此把模板编号映射到分子原子，L4 `_fixed_numbering`（`numbering_engine.py:214`，调用点 `:220`）消费。

```python
# ring_scaffold.py:118（quinoline 条目）
"standard": (NAPH_LABELS, (4, 5, 6, 7, 8, 9, 0, 1, 2, 3)),
```

标签常量定义在 `_TEMPLATES` 之前：`FUSED56_LABELS`（`:48`）、`PURINE_LABELS`（`:49`）、`CARBAZOLE_LABELS`（`:50`）、`ACRIDINE_LABELS`（`:51`）、`PHENOTHIAZINE_LABELS`（`:52`）、`NAPH_LABELS`（`:53`）、`ANTHRACENE_LABELS`（`:54`）、`PHENANTHRENE_LABELS`（`:55`）、`PYRENE_LABELS`（`:56`）、`XANTHENE_LABELS`（`:57`）、`STEROID_LABELS`（`:58`）。

标签形态须与环系拓扑对应：`purine` 桥头碳得纯数字 1–9、无 3a/7a 字母位（P-25.3.3）；`anthracene` 中环 9/10 为全数字、桥头 4a/10a/8a/9a，与萘的中环碳得字母位（`NAPH_LABELS`）不同；`FUSED56_LABELS` 的 `1,2,3,3a,4,5,6,7,7a` 用于 5+6 稠合（`indole`/`indazole`/`benzofuran` 等）。

**import 期自校验**：`_validate_standard_fields`（`:203`）逐一检查——`order` 必须是 `0..n-1` 的排列、`labels` 长度须等于模板原子数，不符即 `ValueError` 中止导入（调用点 `:217`）。`order`/`labels` 与 `smiles` 同条目，改 SMILES 后编号不会静默错位。

未登记 `standard` 的稠环走 P-25.3.3 通用外周编号（`fused_numbering.number_fused_system`），位号形态可能不合保留编号；单杂环走 P-14.4 候选枚举。`TRADITIONAL_NUMBERING_IDS`（`constants.py:147`，8 条：`anthracene`/`phenanthrene`/`acridine`/`carbazole`/`purine`/`xanthene`/`thioxanthene`/`cyclopenta[a]phenanthrene`）中的骨架被 `_fused_numbering` 在 `:264` 直接短路，只能靠同条目的 `standard` 定编号。

### 2.4 单环烃附加组分：`_FUSION_CARBOCYCLES`

`_FUSION_CARBOCYCLES`（`:323`）登记 P-25.3.2.2.1 的单环烃附加零件共 **6 条**（`cyclopropane`/`cyclobutane`/`cyclopentane`/`cyclohexane`/`cycloheptane`/`cyclooctane`，前缀为环名删尾 `ne` 加 `a`：`cyclopropa`…`cycloocta`）。

它们**不是保留母体**，只作稠环拆解时的附加零件，因此**不入 `_TEMPLATES`**——入表会让单环骨架解析成保留名，破坏 P-31 单环通用路径；也不作母体组分（P-25.3.2.1.1）。查询分子由 `_cyclo_component_query`（`:333`）从环大小的纯碳 SMARTS 派生，存于 `_Q_CYCLO`（`:339`）/ `_CYCLO_ELEM`（`:340`）。

配套 API：

| API | 作用 |
|---|---|
| `fusion_carbocycle_prefix(sid)`（`:164`） | 取 (en, zh) 前缀，非该类组分返回 `None` |
| `omits_fusion_numbers(sid)`（`:170`） | 一级单环烃附加组分省略数字位次（P-25.3.8.1；`benzene` 亦特判为省略） |
| `match_fusion_carbocycle(info, atom_ids)`（`:175`） | 骨架精确等于某附加组分时返回 sid |
| `match_fusion_component(info, atom_ids)`（`:189`） | 稠环拆解的组分匹配：保留 mancude 母体优先（`match_retained(..., mancude_only=True)`），其次单环烃附加组分 |
| `retained_fusion_prefix(sid)`（`:156`） | 两者的统一入口（模板 `fused_prefix` → 单环烃前缀） |
| `component_stem(sid)`（`:148`） | 取稠合组分词干（`fused_stem` → 词干名） |

### 2.5 稠环拆解：`fused_system.py`

`decompose_fused_system`（`:217`）是公共入口：由环系的 `sssr_indices` + `fusion_edges` 递归拆成 `FusedNode` 树（`:24`，字段 `scaffold_id`/`atom_ids`/`ring_indices`/`fusion_shared`/`attached`/`fused_stem`/`fused_prefix`/`fused_omit_numbers`）。

- `_candidates_for`（`:53`）：从每个种子环出发增长式枚举保留母体候选（原子集去重、`_has_template_superset`〔`:37`〕超集剪枝）；种子必须本身已是母体候选或单环烃附加组分（`:60`）
- `_select_base`（`:81`）：按 P-25.3.2.4 准则 **(a)–(j)** 逐条收窄选母体组分。**(a)** `_key_a`（`:94`）用 `constants.P25_SENIOR`（`constants.py:41`，N 最优先）；(b) 环数更多；(c) 环大小降序最大；(d) 杂原子总数更多；(e) 杂原子种类更多；(f) `_key_f`（`:106`）用 `P145_SENIOR`；(g)–(j) 经 `_numbered_locants`（`:149`）调 **L4** 的 `preferred_orientations` + `number_fused_system` 得到位次表后，按位次元组收窄。并列时兜底取环集升序最小（`:144`）
- `_decompose`（`:192`）递归：选母体后，剩余环按 `_ring_components`（`:170`）的融合图连通分量递归为附加组分

**新环系作附加组分时，只需要在 `_TEMPLATES` 条目上给 `fused_prefix`**；拆解本身是通用准则，无需改 `fused_system.py`。

### 2.6 FG-环组合的 kind：无需手动注册组合 kind

**无需**为 FG-环组合（pyridinol、cycloalkanol、naphthalenol）手动注册组合 kind——scaffold × FG 已正交化：

- **环 + 主 FG**：`express_ring_principal`（`principal_expression.py:232`）把 kind **收敛为 FG 类别**（alcohol/acid/amine/...，经 `_ring_kind`〔`:169`〕+ `_chain_kind`〔`:175`〕），词干由 scaffold 承载
- **苯 + FG 保留名**：`chain_engine._BENZENE_RETAINED`（`:444`，**10 个 kind**：alcohol/amine/radical/acid/sulfonic/ester/acyl/aldehyde/nitrile/amide）提供 phenol/aniline/benzoic acid/benzaldehyde/benzonitrile/benzamide/benzenesulfonic acid/benzoate/benzoyl/phenyl；其余 scaffold 走通用词干命名（naphthalen-1-ol / pyridine-3-carboxylic acid / imidazol-2-amine）
- **环内主 FG 的准入**：`ring_expression_policy._POLICIES`（`:18`，当前 **16 条**）按 `scaffold.naming_class` 白名单放行**环内**（`in_skeleton`）或**环外**（`exocyclic`）的主 FG 表达。`RingExpressionPolicy`（`:10`）只有三个字段：`naming_classes` / `group_class` / `relations`。`supports_ring_expression`（`:41`）的判定写入 parent 的 `typed_ring_expression_supported`（`principal_expression.py:203`），由 `parent_select._unsupported_typed_ring`（`parent_select.py:33`）拦截不支持者

**新 scaffold 若引入新 `naming_class` 且要支持环内酮/醇/胺（或环外酸/腈/酯/酰胺/酰基头），必须把该 `naming_class` 加进对应策略行**，否则候选被整体拦截成空输出。多元素集合行（如 `:27` 的 `{"monohetero","fused56","purine","carbazole","acridine","phenothiazine","benzodioxole","anthra","phenanthrene","pyrene"}` 环内酮行）就是这类集合的实例。

### 2.7 L2 ↔ L4 双向耦合（改表前必读）

两个模块互为依赖，构成**双向耦合**，两侧符号都是事实上的公共接口：

| 方向 | 符号 | 位置 |
|---|---|---|
| L2 → L4 | `narrow` | `fused_system.py:16`（模块级 import） |
| L2 → L4 | `fused_atoms` / `locant_key` | `fused_system.py:116` / `:117` |
| L2 → L4 | `preferred_orientations` / `number_fused_system` | `fused_system.py:158` / `:159` |
| L4 → L2 | `_Q` / `_Q_H` / `_hydrogenated` | `numbering_engine.py:202`（`_template_matches`） |
| L4 → L2 | `_Q` / `standard_chain` | `numbering_engine.py:220`（`_fixed_numbering`） |
| L4 → L2 | `_STANDARD_LABELS` | `numbering_engine.py:233`、`:240`、`:369`、`:370` |
| L4 → L2 | `get_spec` | `numbering_engine.py:267`（`_fused_numbering`） |
| L4 → L2 | `_STANDARD_ORDERS` | `numbering_engine.py:382`（`fused_component_numbering`） |
| L4 → L2 | `extra_indicated_atoms` | `numbering.py:112`（`_extra_indicated`） |

**约束**：① 通用枚举（L4）与稠环编号/稠环拆解（L4/L2）**共用同一套收窄语义** `narrow`/`narrow_by_senior`——改收窄语义会同时影响链/环/稠环三条路径；② 改 `_TEMPLATES` 的 `standard` 会同时改 L2 拆解选母体（准则 g–j）与 L4 固定编号，这正是单一事实来源的收益，但也意味着 L4 侧的 `_Q`/`_STANDARD_LABELS`/`standard_chain` 不可重命名而不改 L4。

---

## 第三步: 编号（Layer4）—— 通常无需改动

> **源:** `src/namepredict/layer4/numbering_engine.py`（401 行）、`fused_numbering.py`（203 行）、`fused_orientation.py`（364 行）

### 3.1 三层分派

`orient_numbering`（`:325`）按顺序尝试，返回 P-14.4 定向后的原子顺序（不适用返回 `None`）：

1. **`_fixed_numbering`**（`:214`）——保留模板的 `standard` 固定编号：`_template_matches`（`:200`）先直接子图匹配、失败退到氢化骨架（`_Q_H`），每次匹配经 `standard_chain` 转成候选链；有多个匹配时按 `_locant_key_of`（`:245`）用 `_STANDARD_LABELS` 的**标签**比较位次（而非链位置）——蒽的中环 10 位在模板链上排在 5 位之前，按链位置比较会把 10 位判成更低。`_fixed_key`（`:251`）按 (c) 后缀（principal/自由价，`:234-236`）→ (f) 前缀（取代基）→ 引用序收窄
2. **`_fused_numbering`**（`:260`）——P-25.3.3 稠环编号：`TRADITIONAL_NUMBERING_IDS`（`constants.py:147`）直接短路返回 `None`（`:264`）；适用条件（`:269-274`）为链全芳香、或 `spec.n_rings >= 2` 的已注册稠环、或未注册稠环（环数 ≥2）；由 `build_ring_systems` 取的环系须与 chain 原子集一致（`:278`），`sssr` 少于 2 环返回 `None`（`:283`）。随后组装收窄层序（`:293-309`）：**radical 位** → **principal 特征基团位**（P-14.4(c)）→ **`INDICATED_H` 哨兵层**（P-25.3.3.1.2(f)，指示氢位次插在 (c) 之后）→ **取代基位次集** → **hydro 加氢位**（P-31.2.2），末尾按取代基字母序（`:310`）破局。成功时写 `parent["numbering_scaffold"] = {"scaffold_id": "fused", "labels": tuple(labels)}`（`:317`）——**只有两个键**
3. **通用候选枚举**——杂环走 `_narrow_hetero_ring`（`:182`），碳环走 `_ring_cands`（`:33`，每原子 1 号位 × 双向 = 2n 个），链走正反 2 个（`:341-342`）。之后逐条收窄：principal 附着原子（`:345`）→ 多重键端点边位（`_unsat_bonds`:159 / `_bond_locants`:52）→ 取代基位次集（`:353`）→ 字母序（`_stem_loc_pairs`:25）→ CIP（`_chain_rs_codes`:124 / `_rs_locant_key`:137）

`_narrow_hetero_ring`（`:182`）的杂环元素序窄化（P-22.2.2.1.3）：(a)(b) 经 `narrow_by_senior`（`:86`）——全杂原子集最低位次，再按 `constants.P145_SENIOR`（`constants.py:42`，F>Cl>Br>I>O>S>…>N>…）逐元素收窄；(c) 环内带 H 或 3 价取代的 N 得低位（唑 NH=1），`float_hetero=True` 时跳过 (c)。

> **注意 `P25_SENIOR`（`constants.py:41`，N 最优先）与 `P145_SENIOR`（`:42`，F/O/S 先于 N）同源不同序**：`P25_SENIOR` 只用于 P-25.3.2.4 选稠环母体组分（`fused_system._select_base` 的 `_key_a`），编号低位次用 `P145_SENIOR`。

### 3.2 公共收窄原语

`narrow(cands, key_fn, *, reverse=False, skip_none=False)`（`:75`）与 `narrow_by_senior(cands, key_fn, heteros, by_z, *, skip_none=False)`（`:86`）是**公共 API**：L4 的链/环/稠环编号与 L2 的稠环拆解共用同一套收窄语义。`skip_none=True` 表示特征全缺时规则不适用、候选原样返回。

### 3.3 稠环编号原语（`fused_numbering.py` / `fused_orientation.py`）

- `fused_atoms(rings)`（`:14`）——出现在 ≥2 个环的稠合原子
- `_top_rings`（`:23`）/ `_top_atoms`（`:34`）——P-25.3.3.1.1 起点：y 聚行后取 x 最大环，环内起点**只在紧邻稠合原子的非稠合原子中取**（`:41`）
- `_boundary_walk`（`:69`）——沿外部边（`_exterior_edges`:59）单闭环行走外边界，鞋带面积强制顺时针
- `_assign_labels`（`:94`）——非稠合原子与稠合杂原子得下一个数字；稠合碳得「紧邻前数字 + a/b/c 递增」
- `number_fused_system(mol, rings, coords, sub_layers=None, alpha_subs=None)`（`:147`）——按 (a)(b) 杂原子 → (c) 稠合碳 → (d) 稠合杂原子 → `sub_layers` 逐层 → 字母序 → CIP 收窄；`sub_layers` 里的 `INDICATED_H`（`:11`）哨兵走 `_as_indicated`（`:174`，饱和带氢位为**奇数**个时才适用，`:176`）
- `Orientation`（`fused_orientation.py:20`）**只有三个字段**：`row`（水平行环序）/ `coords`（`(原子,x,y)` 元组）/ `quad`（四象限面积分数）；`coord_dict()`（`:27`）转回 dict。`preferred_orientations`（`:328`）枚举所有平局取向（含镜像），排序键 `(len(row), q1, -q3, above)`（`:350`）

### 3.4 指示氢与加氢通道

三个通道均以 `(mol, scaffold_id, match)` 为入参，判定写在 `_TEMPLATES` 模板上，**L4 不写特判**：

| API | 作用 |
|---|---|
| `mancude_ring_atoms(scaffold_id, match)`（`ring_scaffold.py:259`） | 模板 Kekulé 双键（`_kekule_double_atoms`:244，经 `kekulized`）所在的**整个不饱和环**映射到分子：该集合内的 C=C 由母体氢化物名隐含（P-31.1.2）。L2 `principal_expression._implied_ring_atoms`（`:308`）据此排除环内隐含不饱和键 |
| `extra_indicated_atoms(mol, scaffold_id, match)`（`:268`） | 模板同位无 H 而分子有 H 的芳香杂环原子（P-58.2.1）：如 1H-喹啉-4-酮的 N1。**仅对稠合母体（模板 ≥2 环）生效**（`:271`）；L4 `numbering._extra_indicated`（`numbering.py:110`）取它写入指示氢 |
| `hydrogenated_atoms(mol, scaffold_id, match)`（`:284`） | 模板 Kekulé 位在分子中已饱和（P-31.2.2）的加氢位；L2 `principal_expression._scaffold_fields`（`:184`）写入 `hydro_atoms`，L4 `numbering.number`（`numbering.py:82`）经 `hydro_prefix`（`:62`）生成 hydro 前缀。环杂原子失双键新增的 H 归指示氢，不计入 hydro 计数（`:306`） |

位次换算统一走 `locant_calc.atom_locant`（`locant_calc.py:33`，公开别名），排序键走 `locant_calc.locant_key`——**不要在下游自己算位次**。

若新环系的保留名不隐含某些环内双键或需要额外的指示氢，只需核对模板的 Kekulé 形式与上述三通道的判定，无需改 L4。

---

## 第四步: 名称组装（Layer5）

> **源:** `src/namepredict/layer5/assembler.py`、`fused_namer.py`（131 行）、`chain_engine.py`（529 行）

### 4.1 保留名环系：词干注入

`_TEMPLATES` 的 `stem_en`/`stem_zh` 经 `ScaffoldSpec` → `KindMeta`（`kind_registry._load_from_scaffold_specs`，`kind_registry.py:76`）→ `pack_parent_stem`（`:58`）写入 parent 的 `stem_en`/`stem_zh`。L5 侧在 `assembler._names_for`（`:301`）消费：

```python
# assembler.py:311-326
if kind == "alkane" and parent.get("fused_tree") and sid != "benzene":
    return _parent_stem_names(numbered)      # 未注册稠环：L5 已组好稠合 base 名
if sid == "benzene" and kind == "alkane":
    return ("benzene", "苯")
if parent.get("stem_en") and parent.get("stem_zh"):   # 稠环/杂环 scaffold 词干注入
    entry = replace(entry, stem=(parent["stem_en"], parent["stem_zh"]), coda="",
                    omit_rule=lambda n, loc, omit: bool(omit), aromatic=(sid == "benzene"))
elif sid == "carbocycle":
    entry = replace(entry, cyclic=True, ene_loc_omit=True, omit_rule=...)
sc_variant = _BENZENE_RETAINED[kind] if sid == "benzene" and kind in _BENZENE_RETAINED \
    else (entry.variant or {}).get(sid)
```

英文词尾 `e` 的省略由 `chain_engine._elide_parent_e`（`:224`）按后缀首字母判定（P-60.2(a)，调用点 `:238`、`:390`）——diol/dione/carbaldehyde 等辅音开头后缀保留 `e`。

### 4.2 稠合名：`fused_namer.fused_parent_names`

未注册稠环（或已注册但走稠合组装的环系）由 `assembler._ensure_fused_stem`（`:250`）经 `fused_namer.fused_parent_names`（`:112`）现场组装 base 名，再补指示氢前缀（`_indicated_h_prefix`，`:244`）。

`fused_parent_names` 的路径：`node.attached` 为空 → 返回 `None`（单节点保留名由 `_parent_stem_names` 处理）；否则由 `_collect_attached`（`:99`）递归收集附加组分前缀。单级组装在 `_fused_one`（`:66`）：

- 母体与附加组分**共用同一稠合共享原子集编号**——`_component_numbering`（`:21`）把 `shared` 作取代基传入 L4 `fused_component_numbering`（`numbering_engine.py:376`），母体外周再剔除出现在 ≥3 环的 perifused 中心（`_inner_atoms`:35）
- 融合描述符 = 附加组分共享原子位次 + 母体外周侧字母（`_fusion_letter`:44 / `_fusion_numbers`:56），形如 `[2,3-b]`
- `fused_omit_numbers=True` 的一级单环烃附加组分只出 `[b]`（`:78`）

**整名替换**：组装出的名字命中 `constants.RETAINED_FUSION_ALIASES`（`constants.py:186`，7 条：`benzo[c]furan`→`2-benzofuran`、`benzo[c]pyrrole`→`isoindole`、`benzo[b]benzofuran`→`dibenzofuran`、`benzo[b]quinoxaline`→`phenazine`、`benzo[a]indene`→`fluorene`、`benzo[d]1,2-oxazole`→`1,2-benzoxazole`、`benzo[b]anthracene`→`tetracene`）时整名替换（`fused_namer.py:130`）。

### 4.3 特殊命名需求

- 环外主基的后缀改写：`chain_engine._exo_ring_spec`（`:298`）消费 `constants.EXO_RING_SUF`（`constants.py:161`，6 类：acid/aldehyde/ester/amide/nitrile/acyl）；`_ring_prefix_located`（`:291`）判定主基位次是否必带
- 前缀构建规则：`layer5/assembler_prefixes.py`
- 苯环专属保留名：`chain_engine._BENZENE_RETAINED`（`:444`）——新增苯系保留名加在这里（key 为 kind）
- kind → 词尾分派：`chain_engine._KIND_TABLE`（`:457`，15 entry）

---

## 第五步: FG-环组合支持

| 情形 | 承担路径 |
|---|---|
| 苯单取代保留名（10 kind） | `chain_engine._BENZENE_RETAINED`（`:444`），经 `assembler._names_for:323` 写入 `variant` |
| 稠环/杂环 + 主 FG | `express_ring_principal`（`principal_expression.py:232`）收敛为 FG 类别 kind，词干由 parent 的 `stem_en`/`stem_zh` 注入 |
| 环外酸/醛/酯/酰胺/腈/酰基头 | `chain_engine._exo_ring_spec`（`:298`）+ `constants.EXO_RING_SUF`（`:161`）→ `cyclohexanecarboxylic acid` / `cyclohexanecarbaldehyde` / `cyclohexanecarboxamide` / `cyclohexanecarbonitrile` / `cyclohexanecarbonyl` / `cyclohexanecarboxylate`；多取代拼倍数词且位次必带（P-65.2.2 / P-66.6.1.1.3），苯单取代回落保留名 |
| 环内 keto/alcohol/amine 与环外 acid/nitrile 的准入 | `ring_expression_policy.supports_ring_expression`（`:41`）按 `naming_class` 白名单判定，见 2.6 |

**这些路径都与 scaffold 身份正交**——新增 scaffold 只需保证 `naming_class` 落在 `_POLICIES` 里，无需为每个 FG 组合另加代码。

---

## 工作示例: 新增吡啶（pyridine）

### Step 1 — 环系检测

吡啶的六元芳环（5C + 1N）由 `_sssr`（`ring_systems.py:19`）自动检出，单环自成一个分量（`_components`:54）。Layer1 无改动。

### Step 2 — Layer2 骨架识别

pyridine 走**数据驱动**识别，全部在 `ring_scaffold.py`：

1. **模板注册**：`_TEMPLATES` 的 pyridine 条目（`:71`）——`smiles="n1ccccc1"`、`stem_en="pyridine"`、`stem_zh="吡啶"`、`naming_class="monohetero"`、`fused=True`、`fused_prefix=("pyrido","吡啶并")`。无 `standard`（单杂环走 P-14.4 候选枚举）。`_spec_from_template`（`:343`）派生 `ScaffoldSpec`（`retained=True`，`n_rings=1`，`ring="hetero"`）
2. **子图同构匹配**：`resolve_ring_scaffold`（`:450`）→ `match_retained`（`:387`）的元素标注子图同构精确覆盖 → `get_spec("pyridine").identity`，`scaffold_id="pyridine"`
3. **词干注册**：`_load_from_scaffold_specs`（`kind_registry.py:76`）自动把 `("pyridine","吡啶")` 注册为 `KindMeta`

### Step 3 — 编号

无需 orienter，也不需要 `standard`。`orient_numbering`（`:325`）在 `_fixed_numbering`/`_fused_numbering` 均返回 `None` 后，走通用环枚举 `_ring_cands`（`:33`）；因环内含非碳原子，分派到 `_narrow_hetero_ring`（`:182`）：`narrow_by_senior`（`:86`）先取全杂原子集最低位次、再按 `P145_SENIOR` 逐元素收窄，使唯一的 N 定在 1 号位。

### Step 4 — 命名

`pack_parent_stem`（`kind_registry.py:58`）把 `stem_en`/`stem_zh` 写入 parent，`assembler._names_for`（`:301`）在 `:315` 注入词干 `("pyridine","吡啶")` 到 `_KIND_TABLE["alkane"]`（`chain_engine.py:457`）的 entry。FG 变体（pyridinecarboxylic acid 等）由 `express_ring_principal`（`:232`）收敛为 FG 类别后走同一链引擎。

### Step 5 — 测试

- `c1ccncc1` — pyridine / 吡啶
- `c1ccncc1C(=O)O` — pyridine-2-carboxylic acid（吡啶-2-甲酸）
- `c1ccncc1N` — pyridin-2-amine（吡啶-2-胺）
- `c1ccncc1O` — pyridin-2-ol / 吡啶-2-醇（验证 `_POLICIES` 环内醇/酮放行）

---

## 文件改动清单

| Layer | 文件 | 改动内容 | 是否必需 |
|-------|------|---------|---------|
| L1 | `layer1/ring_systems.py` | 仅当环系划分要突破「稠合边连通分量」时扩展（入口 `build_ring_systems`:117） | 否 |
| L2 | `layer2/ring_scaffold.py` | **主改动**：在 `_TEMPLATES`（`:60`）新增 `{smiles, stem_en, stem_zh, naming_class}`；作附加组分补 `fused_prefix`，作母体组分补 `fused`；词干含加氢前缀时补 `locant_prefix`（+ `prefix_nh_conditional`）；非默认编号时在同条目加 `standard = (labels, order)`（标签常量定义在 `:48-58`） | **是** |
| L2 | `layer2/ring_expression_policy.py` | 新 `naming_class` 要支持环内酮/醇/胺或环外酸/腈时，把该 `naming_class` 加进 `_POLICIES`（`:18`）对应策略行 | 视 FG 而定 |
| L2 | `layer2/fused_system.py` | 通常无需改动（拆解走 `_select_base`:81 的通用准则 (a)–(j)） | 否 |
| L2 | `layer2/kind_registry.py` | 通常无需改动（`_load_from_scaffold_specs`:76 自动注册词干） | 否 |
| L4 | `layer4/numbering_engine.py` | 通常无需改动（三层分派 + `standard` 固定编号自动适用） | 否 |
| L4 | `layer4/fused_numbering.py` | 通常无需改动（无 `standard` 的稠环走 P-25.3.3 外周编号） | 否 |
| L4 | `layer4/fused_orientation.py` | 通常无需改动（优选取向由环几何通用决定） | 否 |
| L5 | `layer5/assembler.py` | 通常无需改动（词干在 `_names_for`:315 自动注入；稠合名走 `_ensure_fused_stem`:250） | 否 |
| L5 | `layer5/fused_namer.py` | 通常无需改动 | 否 |
| L5 | `layer5/chain_engine.py` | 仅当新增苯环专属保留名（`_BENZENE_RETAINED`:444）或新 kind（`_KIND_TABLE`:457）时 | 否 |
| — | `constants.py` | 仅当新增稠合组装名 → 保留名映射（`RETAINED_FUSION_ALIASES`:186）或环外主基后缀（`EXO_RING_SUF`:161）时 | 否 |
| — | `tests/` | 添加 SMILES 测试用例 | **是** |

---

## 常见陷阱

1. **`_TEMPLATES` 表序即优先级，元素签名相同的条目会撞车。** `_match_with_map`（`:393`）按 dict 序返回首个命中。`pyran`（`:75`）必须排在 `oxane`（`:89`）之前，否则 6 元含氧 mancude 环的环内 C=C 被静默丢弃。
2. **`standard` 写错在 import 期就炸。** `_validate_standard_fields`（`:203`，调用点 `:217`）要求 `order` 是 `0..n-1` 的排列、`labels` 长度等于模板原子数，不符即 `ValueError`。
3. **`standard` 的 `labels` 形态要与环系拓扑对应**。「稠环中环碳得字母位」不是通例：`anthracene` 中环 9/10 是全数字，`NAPH_LABELS` 的中环碳才是字母位；`purine` 桥头得纯数字 1–9。照抄邻表的标签元组会出错名。
4. **单环烃附加组分不要塞进 `_TEMPLATES`。** `_FUSION_CARBOCYCLES`（`:323`）成员只作稠合零件；入 `_TEMPLATES` 会让单环骨架解析成保留名，破坏 P-31 单环通用路径。
5. **`fused=False` 的饱和保留名不能作稠合组分。** `match_fusion_component`（`:189`）与 `_match_with_map`（`:400`、`:410`）都按 `mancude_only` 过滤——`pyrrolidine`/`piperidine`/`oxolane` 等只作母体，不作稠合零件。反之要当附加组分就必须 `fused=True` + `fused_prefix`。
6. **`_Q_H` 不含 `mono_carbo`。** `:237` 显式排除 `naming_class == "mono_carbo"`（苯），故苯不参与氢化骨架匹配——加氢苯衍生物不会因氢化回退误判为 benzene。
7. **新 `naming_class` 不加进 `_POLICIES` → 空输出。** 环内酮/醇/胺和环外酸/腈的 typed 表达由 `supports_ring_expression`（`:41`）白名单放行，不支持者被 `parent_select._unsupported_typed_ring`（`parent_select.py:33`）整体拦截，候选数为 0 而不是报错。
8. **`P25_SENIOR` 与 `P145_SENIOR` 同源不同序，勿混用。** `constants.py:41` 的 `P25_SENIOR`（N 最优先）只用于选稠环母体组分；`constants.py:42` 的 `P145_SENIOR`（F/O/S 先于 N）用于「哪个杂原子得低位次」。
9. **`kekulized` 失败返回 `None`，调用方必须降级。** `_kekule_double_atoms`（`:251`）退回空双键集、`indicated_hydrogen.saturated_ring_atoms`（`indicated_hydrogen.py:10`，`:17` 调 `kekulized`）退回原 mol、`principal_expression._kekule_ring_dbs`（`:291`，`:295` 调 `kekulized`）放弃补环内 C=C——新环系的模板若无法 Kekulé 化，加氢/指示氢判定会静默退化。
10. **`labels` 与链长不符时 L4 静默回落链序号。** `_component_labels`（`:363`）在 `numbering_scaffold` 与 `_STANDARD_LABELS` 都取不到匹配长度的标签时返回 `1..n` 纯数字——位号形态错了但不会失败，容易漏检。
11. **`TRADITIONAL_NUMBERING_IDS`（`constants.py:147`）里的骨架只能靠 `standard` 定编号。** `_fused_numbering`（`:264`）对它们直接返回 `None`，未登记 `standard` 就会掉到通用外周编号。
12. **螺环不合并。** `build_ring_systems`（`:117`）只用稠合边做连通分量，`_ring_pairs`（`:27`）产出的 `spiro` 列表不被消费——螺环两侧各自是一个独立环系。
13. **别在下游自己算位次。** 位次换算统一走 `locant_calc.atom_locant`（`locant_calc.py:33`），排序走 `locant_key`。
14. **改 L2/L4 公共面要两侧同步。** `_Q`/`_Q_H`/`_hydrogenated`/`_STANDARD_LABELS`/`_STANDARD_ORDERS`/`standard_chain`/`get_spec`/`extra_indicated_atoms` 被 L4 直接 import；`narrow`/`narrow_by_senior`/`fused_atoms`/`locant_key`/`preferred_orientations`/`number_fused_system` 被 L2 直接 import（见 2.7）。

---

## 相关文档

- [[architecture/layer1-analyzer]] — Layer1 官能团与环系事实
- [[architecture/layer2-parent-selector]] — Layer2 母体选择器完整架构
- [[architecture/layer4-numbering]] — Layer4 编号引擎（P-14.4 三层分派 + 稠环编号）
- [[architecture/layer5-name-assembly]] — Layer5 名称组装
- [[concepts/functional-group-priority]] — FG 优先级与 IUPAC P-44 规则
- [[reference/core-data-contracts]] — parent dict / numbered dict 字段表
- [[guides/adding-new-functional-group]] — 新增官能团指南（thiol 工作示例）
