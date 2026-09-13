# 新增环系指南

本文档分步说明如何为 NamePredict 流水线新增一个环系（ring system）——无论是新的 fused heterocycle（稠杂环）、bridged system（桥环）、spiro system（螺环），还是新的单杂环。

---

## 前置知识

阅读本指南前，建议先了解以下架构文档：

- [[architecture/overview]] — 六层流水线总览
- [[architecture/layer2-parent-selector]] — Layer2 母体选择器的核心架构（候选生成 + 评分排序）
- [[architecture/layer4-numbering]] — Layer4 编号引擎（P-14.4 候选管线，无需注册 orienter）
- [[architecture/layer5-name-assembly]] — Layer5 中英双语名称组装
- [[concepts/functional-group-priority]] — FG 优先级规则

---

## 流程概览

新增一个环系涉及以下改动，按顺序为：

1. **Layer1 环检测**（通常无需改动，SSSR 已自动识别）
2. **Layer2 骨架识别**（主要工作量——在 `ring_scaffold.py` 的 `_TEMPLATES` 注册模板；新 `naming_class` 还要在 `ring_expression_policy._POLICIES` 放行环内/环外主 FG）
3. **Layer4 编号**（**通常无需改动**——`numbering_engine` 的 P-14.4 候选管线自动适用）
4. **Layer5 名称组装**（添加 `_KIND_TABLE` 词干 / variant）
5. **FG-环组合支持**（如环上有 COOH/CHO/CN/OH/NH2 等 FG 变体）

> 核心原则：**scaffold × kind × 数量全正交**——加新环只加 `_TEMPLATES` 模板 + 词干表，不枚举 cycloalcohol/naphthalenol/indolol 组合 kind。

---

## 第一步: 环检测前提（Layer1）

大多数环系已被 Layer1 的 `ring_systems.py` 通过 **SSSR（Smallest Set of Smallest Rings）+ Union-Find 融合** 自动检测，无需额外工作。

> **源:** `src/namepredict/layer1/ring_systems.py`

`sssr_rings`（`ring_systems.py:14`）利用 RDKit 的 `GetRingInfo().AtomRings()` 获取所有最小环，然后 `build_ring_systems`（`:289`）通过 fusion graph（共享 >=2 原子的边）将稠合环合并为 ring system。桥环（共享 >=3 原子）和螺环（共享 1 原子）也被识别。结果由 `analyzer._ring_meta`（`analyzer.py:380`）写入 `info["ring_systems"]`；L2 消费其中的 `atom_ids`/`sssr_indices`/`fusion_edges`（`parent_skeleton._ring_candidates`，`:111`；`fused_system.decompose_fused_system`，`:225`）。

如果你的环系需要**超越 SSSR 的特殊识别**，则需扩展 ring_systems.py 或创建新的检测模块。

---

## 第二步: 母体选择（Layer2）—— 主要工作量

### 2.1 骨架识别入口

母体选择统一走 P-44 规则管线（`rule_driven_parent_candidates`，`principal_parent.py:45`），环骨架身份由 **`ring_scaffold.py` 的 `resolve_ring_scaffold`**（`ring_scaffold.py:477`）在表达阶段解析：

```
resolve_ring_scaffold(info, skeleton)
  ├─ 1. get_identity(skeleton.scaffold_id)  # 模板 id 直接命中
  ├─ 2. match_retained(info, atom_ids)      # SMILES 模板子图同构精确覆盖
  └─ 3. _generic_carbocycle                 # 全碳非保留环 → ScaffoldIdentity("carbocycle",...)
```

支持 API：`get_identity`（`:383`）、`match_retained`（`:410`）、`_generic_carbocycle`（`:464`）。

### 2.2 注册保留母体（唯一事实来源 `_TEMPLATES`）

环骨架身份由 `resolve_ring_scaffold`（`ring_scaffold.py:477`）识别：① `get_identity(skeleton.scaffold_id)` 显式命中 → ② **SMILES 模板子图同构**（`match_retained`，`:410`）→ ③ 兜底 `_generic_carbocycle`（`:464`，全碳非保留环 → carbocycle，非全碳多环 → `fused_hetero`，非全碳单环 → None 显式失败）。

**新增保留环系（如 quinoline）只改 `ring_scaffold.py` 的 `_TEMPLATES` 一张表**（`ring_scaffold.py:71`，当前 83 条）——加一条 `{smiles, stem_en, stem_zh, naming_class}`：

```python
# ring_scaffold.py:129（`_TEMPLATES` 内 quinoline 条目）
"quinoline":    {"smiles": "c1ccc2ncccc2c1", "stem_en": "quinoline", "stem_zh": "喹啉", "naming_class": "naph_family", "fused": True, "standard": (NAPH_LABELS, (4, 5, 6, 7, 8, 9, 0, 1, 2, 3))},
```

可选字段：

- `locant_prefix`（如 `"1,3-"`/`"1H-"`，位置异构体与指示氢前缀）与 `prefix_nh_conditional`（NH 存在时才加 `1H-`）；词干本身已把 locant 前缀嵌在词中（如 `4,5-dihydro-1,3-thiazole` 的 `1,3-`）时，`kind_registry.pack_parent_stem`（`kind_registry.py:65`）的 `_embeds_locant_prefix`（`:60`）让前缀留在组分名原位，不会前移或剥离
- `fused`（True 才可作稠合命名零件）、`fused_stem`（稠合前缀取去指示氢的组分名，如 `indole`→`indolo`）、`fused_prefix`（附加组分的保留稠合前缀，如 `benzo`/`吡啶并`）。表征结构的位次（杂原子位置）在稠合名中按 P-25.3.1.3 置于方括号内，且 P-25.3.2.1.2 强制 isoxazole/oxazole/thiazole 在稠合名中改用 Hantzsch-Widman 名，故 `[1,2]oxazolo`/`[1,3]oxazolo`/`[1,3]thiazolo`/`[1,2,4]triazolo`/`[1,3,5]triazino` 一律手写方括号——通用「去尾 e 加 o」不含方括号，会拼出 `1,2,4-triazolo` 而非 `[1,2,4]triazolo`
- `standard`（固定编号，见下）

当前 83 条按命名类分布：`mono_carbo` 1 条（`benzene`）、`naph_family` 9 条（`naphthalene`/`quinoline`/`isoquinoline`/`quinazoline`/`quinoxaline`/`cinnoline`/`chromene`/`isochromene`/`pteridine`）、`fused56` 8 条（`indole`/`indazole`/`benzimidazole`/`benzofuran`/`benzothiophene`/`benzothiazole`/`benzoxazole`/`indene`）、`monohetero` 53 条（mancude 杂芳环 `furan`/`thiophene`/`pyrrole`/`pyridine`/`pyran`/`imidazole`/`pyrazole`/`oxazole`/`thiazole`/`isoxazole`/`triazole`/`tetrazole`/`triazine`/`dioxine`/`oxadiazole124`/`oxadiazole134`/`oxadiazole125`/`thiadiazole134`/`thiadiazole124`/`triazine124`/`tetrazine1245`/`thiazole12`/`thiazine13`；饱和/部分饱和杂环 `pyrrolidine`/`piperidine`/`morpholine`/`piperazine`/`oxolane`/`oxane`/`oxirane`/`aziridine`/`oxetane`/`azetidine`/`thiolane`/`thiane`/`dioxolane`/`dioxane`/`trioxane`/`oxazolidine`/`imidazolidine`/`thiazolidine`/`dihydrofuran`/`dihydropyran`/`dihydropyrrole`/`dihydroimidazole`/`dihydrothiazole`/`oxepane`/`azepane`/`oxazepane`/`thiazepane`），以及各单条命名类（`anthracene`/`phenanthrene`/`pyrene`/`chrysene`/`carbazole`/`acridine`/`phenothiazine`/`benzodioxole`/`purine`/`xanthene`/`thioxanthene`/`cyclopenta[a]phenanthrene`）。

其中 67 条带 `fused=True`、27 条带 `fused_prefix`、3 条带 `fused_stem`（`indene`/`indole`/`purine`）、41 条带 `locant_prefix`、36 条带 `standard`。

`pyran` 是「表序即优先级」的实例：6 元含氧 mancude 母体须排在 `oxane` 之前，否则环内 C=C 被静默丢弃（`pyran-2-one` → `oxan-2-one`）。`benzofuran`/`benzothiophene` 另带 `locant_prefix="1-"`。

`_spec_from_template`（`ring_scaffold.py:359`）自动派生 `ScaffoldSpec`（n_rings/ring 从 smiles 算，retained=True），`all_specs()`（`:393`）/`get_spec()`（`:388`）/`get_identity()`（`:383`）均由此派生；`kind_registry._load_from_scaffold_specs`（`kind_registry.py:83`）据此自动注册 KindMeta 词干（bootstrap 唯一一步，ScaffoldSpec 是词干权威）。位置异构体（quinoline/isoquinoline、二嗪、二唑等）在元素标注的子图同构下天然区分，无需额外消解。

#### 固定编号：条目内的 `standard = (labels, order)`

若新环的保留编号不是"按模板原子序从 1 起"，在**同一条目内**加 `standard = (labels, order)` 二元组：`order` 是模板原子按 locant 顺序的下标排列，`labels` 为对应的 locant 标签元组，**两者长度均等于模板原子数**。`_STANDARD_LABELS`（`ring_scaffold.py:206`）/`_STANDARD_ORDERS`（`:209`）是它的派生视图（当前 36 条登记），`standard_chain`（`:451`）据此把模板编号映射到分子原子（L4 `numbering_engine` 的 `_fixed_numbering`，`numbering_engine.py:204` 调用）。`_fixed_numbering` 比较候选位次时按键内标签而非链位置（`_locant_key_of`，`numbering_engine.py:236`）：蒽的中环 10 位在模板链上排在 5 位之前，按链位置比较会把 10 位判成更低。

import 期 `_validate_standard_fields()`（`ring_scaffold.py:214`）逐一校验：`order` 必须是 `0..n-1` 的排列、`labels` 长度须等于模板原子数，不符即 `ValueError` 中止导入。`order`/`labels` 与 `smiles` 同条目，改 SMILES 后编号不会静默错位。标签序列常量定义在 `_TEMPLATES` 之前（`FUSED56_LABELS`:59、`PURINE_LABELS`:60、`CARBAZOLE_LABELS`:61、`ACRIDINE_LABELS`:62、`PHENOTHIAZINE_LABELS`:63、`NAPH_LABELS`:64、`ANTHRACENE_LABELS`:65、`PHENANTHRENE_LABELS`:66、`PYRENE_LABELS`:67、`XANTHENE_LABELS`:68、`STEROID_LABELS`:69）。标签形态与环系拓扑对应：`purine` 桥头碳得纯数字 1–9、无 3a/7a 字母位（P-25.3.3）；`anthracene` 中环 9/10 为全数字、桥头 4a/10a/8a/9a，与萘的中环碳得字母位（`NAPH_LABELS`）不同。未登记 `standard` 的稠环走 P-25.3.3 通用外周编号（对称碳环 naphthalene、`chrysene` 等），位号形态可能不合保留编号；单杂环走 P-14.4 候选枚举（pyrrole/pyridine 等）。

#### 单环烃附加组分：`_FUSION_CARBOCYCLES`

`_FUSION_CARBOCYCLES`（`ring_scaffold.py:339`）登记 P-25.3.2.2.1 的单环烃附加零件共 6 条（`cyclopropane`/`cyclobutane`/`cyclopentane`/`cyclohexane`/`cycloheptane`/`cyclooctane`，饱和环名删尾 'ne' 得 `cyclopropa`/`cyclohexa` 前缀）。它们**不是保留母体**，只作稠合拆解的附加零件，故不入 `_TEMPLATES`（入表会让单环骨架解析成保留名，破坏 P-31 单环通用路径）；也不是母体组分（P-25.3.2.1.1）。配套 API：

| API | 作用 |
|---|---|
| `fusion_carbocycle_prefix(sid)`（`:175`） | 取 (en, zh) 前缀，非该类组分返回 None |
| `omits_fusion_numbers(sid)`（`:181`） | 一级单环烃附加组分省略数字位次（P-25.3.8.1） |
| `match_fusion_carbocycle(info, atom_ids)`（`:186`） | 骨架精确等于某附加组分时返回 sid |
| `match_fusion_component(info, atom_ids)`（`:200`） | 稠环拆解的组分匹配：保留 mancude 母体优先，其次单环烃附加组分 |

`retained_fusion_prefix`（`:167`）是两者的统一入口（模板 `fused_prefix` → 单环烃前缀），`component_stem`（`:159`）取稠合组分词干（`fused_stem` → 词干名）。

稠环拆解（`layer2/fused_system.py`，`_select_base`（`:85`）走 P-25.3.2.4 准则 (a)–(j)，入口 `decompose_fused_system`（`:225`））经 `match_fusion_component` 取种子环与母体候选。

> 无 `_TOPOLOGY` 五元组表与手写的 `_ALL_SPECS`；`_ALL_SPECS`（`:376`）由 `_TEMPLATES` 派生，查询 API 为 `get_spec`（`:388`）/`get_identity`（`:383`）/`all_specs`（`:393`）。

### 2.3 FG-环组合的 kind（无需手动注册组合 kind）

**无需**为 FG-环组合（如 pyridinol、cycloalkanol、naphthalenol）手动注册组合 kind——scaffold × FG 已正交化。环上带 FG 时：

- **环 + 主 FG**：`express_ring_principal`（`principal_expression.py:256`）把 kind **收敛为 FG 类别**（alcohol/acid/amine/...，经 `_ring_kind`，`:176`），词干由 scaffold 承载（`pack_parent_stem`〔`kind_registry.py:65`〕按 scaffold_id 注入）；苯/饱和环/稠环/杂环一律平等
- **苯 + FG 保留名**：chain_engine `_KIND_TABLE` 各 entry 的 `variant` 提供 phenol/benzoic acid/aniline/benzaldehyde/benzonitrile/benzamide 等（苯专属）；其余 scaffold 走通用词干命名（naphthalen-1-ol / pyridine-3-carboxylic acid / imidazol-2-amine）
- 命名时 `chain_engine._KIND_TABLE` 的 entry 按 `scaffold_id` 注入词干（`assembler._ring_stem`〔`:27`〕从 parent 的 stem_en/zh 取完整词干；结尾 e 的省略由 chain_engine `_elide_parent_e`〔`:271`〕按后缀首字母判定，P-60.2(a)——diol/dione/carbaldehyde 等辅音开头后缀保留 e）
- **环内主 FG 的额外准入**：`layer2/ring_expression_policy.py` 的 `_POLICIES`（`:18`）按 `scaffold.naming_class` 白名单放行环内（`in_skeleton`）或环外（`exocyclic`）的主 FG 表达，`supports_ring_expression`（`:41`）判定结果写入 `typed_ring_expression_supported`，由 `principal_parent._unsupported_typed_ring`（`:28`）拦截不支持者。**新 scaffold 若引入新 `naming_class` 且要支持环内酮/醇/胺（或环外酸/腈），必须把该 `naming_class` 加进对应策略行**，否则候选被整体拦截成空输出

---

## 第三步: 编号（Layer4）—— 通常无需改动

**`numbering_engine.orient_numbering`**（`:316`，P-14.4 候选管线）按 kind 无关的统一规则收窄：

1. 枚举候选编号：链正反（2 个）/ 环每原子 1 号位 × 双向（2n 个）
2. 稠环先走固定编号或外周编号：`_fixed_numbering`（`:204`，`standard` 登记的模板编号经 `standard_chain` 映射）→ `_fused_numbering`（`:251`，P-25.3.3 优选取向 + 外周字母位）
3. 杂环走 `_narrow_hetero_ring`（`:168`）元素序窄化定起点，其余按 principal FG 附着原子 → 多重键 → 取代基位次集逐条收窄，末尾按字母序与 CIP（`_rs_locant_key`，`:123`）破局

因此**新增环系无需写 orienter**。唯一需要保证的是 L2 正确注入 `scaffold_id`。固定编号事实（稠环 `3a`/`4a` 标签）由 `numbering_scaffold_facts`（`ring_scaffold.py:398`）基于 `_TEMPLATES` 的 `standard` 生成，`kind_registry._attach_numbering_scaffold`（`kind_registry.py:33`）把它与 `numbering_scaffold_required` 一并写入 parent，L4 经 `standard_chain`（`ring_scaffold.py:451`）映射 `standard` 登记的固定编号。

杂环的位次收窄在 `_narrow_hetero_ring`（`numbering_engine.py:168`）：(a) 全杂原子集最低位次 → (b) 按 `constants.P145_SENIOR`（`constants.py:42`，F>Cl>Br>I>O>S>…>N>…）逐元素收窄 → (c) 同元素 N 中带 H/3 价取代者得低位（唑 NH=1）。**注意与 `P25_SENIOR`（`:41`，N 最优先）同源不同序**：`P25_SENIOR` 只用于 P-25.3.2.4 选稠环母体组分，编号低位次用 `P145_SENIOR`。

稠环的收窄层序由 `_fused_numbering`（`numbering_engine.py:251`）组装：后缀 (c) principal → `INDICATED_H` 哨兵层（指示氢位次最小化，`fused_numbering.py:10`）→ 取代基前缀；层级混为一集会把前缀位次判到指示氢之前。稠合点集合另作 `sub_layers` 传入 `fused_component_numbering`（`numbering_engine.py:370`），供多环无固定编号组分的镜像对（苯并咪唑 N1/N3 互换）逐层收窄；`fused_numbering._top_atoms`（`fused_numbering.py:33`）的外周行走起点只在该环紧邻稠合原子的非稠合原子里取。

#### 指示氢通道（保留母体）

保留母体模板隐含的环内不饱和度与 H 分布由以下通道补足，均以 (mol, scaffold_id, match) 为入参：

| API | 作用 |
|---|---|
| `mancude_ring_atoms(scaffold_id, match)`（`ring_scaffold.py:274`） | 模板 mancude（Kekulé 双键）**整个不饱和环**位映射到分子：该集合内的 C=C 由母体氢化物名隐含（P-31.1.2），不得再写成 -ene/-yne；L2 `principal_expression._implied_ring_atoms`（`:336`）据此排除环内隐含不饱和键 |
| `extra_indicated_atoms(mol, scaffold_id, match)`（`:283`） | 保留母体名未隐含、而分子中该芳香杂环位带 H 的原子（P-58.2.1）：模板同位无 H 而分子有 H，如 1H-喹啉-4-酮的 N1。仅对稠合母体（≥2 环）生效；L4 `numbering._extra_indicated`（`layer4/numbering.py:115`）取它写入 `indicated_h_forced` |
| `hydrogenated_atoms(mol, scaffold_id, match)`（`:300`） | 模板 Kekulé 位在分子中已饱和（P-31.2.2）的加氢位，L2 `principal_expression._scaffold_fields`（`:229`）写入 `hydro_atoms`；L4 `numbering.number`（`layer4/numbering.py:84`）据此经 `hydro_prefix`（`:63`）生成 hydro 前缀；环杂原子失双键新增的 H 归指示氢，不计入 hydro 计数 |

若新环系的保留名不隐含某些环内双键或需要额外的指示氢，按上述通道核对 `mancude_ring_atoms`/`extra_indicated_atoms` 的模板判定即可，无需在 L4 写特判。

> **源:** `src/namepredict/layer4/numbering_engine.py`, `src/namepredict/layer2/ring_scaffold.py:398`

> 注：无 `NumberingPlan` 机制与 `locants/` 子包——编号完全走候选枚举 + `chain.index + 1`，不依赖固定编号 plan。

---

## 第四步: 名称组装（Layer5）

### 4.1 保留名环系

对于保留名环系（苯、吡啶、萘、吲哚等），词干在 `ring_scaffold.py` 的 `_TEMPLATES`（`stem_en`/`stem_zh`）定义。`chain_engine._KIND_TABLE` 的 entry 在 `_names_for` 中按 `scaffold_id` 运行时替换：

```python
# assembler.py `_names_for`（:372；词干/variant 注入在 :389-:404）
if kind == "alkane" and parent.get("fused_tree") and sid != "benzene":
    return _parent_stem_names(numbered)   # 未注册稠环：L5 fused_namer 已组好 base 名
if sid == "benzene" and kind == "alkane":
    return ("benzene", "苯")
ring_stem = _ring_stem(numbered)          # 从 parent 的 stem_en/stem_zh 派生
if ring_stem:                             # 稠环/杂环 scaffold 词干
    entry = replace(entry, stem=ring_stem, coda="",
                    omit_rule=lambda n, loc, omit: bool(omit), aromatic=(sid == "benzene"))
elif sid == "carbocycle":
    entry = replace(entry, cyclic=True, ene_loc_omit=True, omit_rule=...)
sc_variant = (entry.variant or {}).get(sid)   # 苯环保留名等 scaffold 专属特例
```

未注册稠环（`scaffold_id` 不在 `_TEMPLATES`）由 `assembler._ensure_fused_stem`（`:352`）经 `layer5/fused_namer.fused_parent_names`（`fused_namer.py:126`）现场组装稠合 base 名，再补指示氢前缀。

### 4.2 特殊命名需求

如果你的环系需要特殊的命名处理，则可能需要：

- 在 `assembler.py` 的 `_names_for`（`:372`）添加 worker（环外主基已有通用 worker `_exocyclic_ring_names`〔`:78`〕+ 后缀表 `EXO_RING_SUF`〔`constants.py:179`〕，新环系无需再加）
- 在 `assembler_prefixes.py` 中注册前缀构建规则
- 在 `chain_engine.py` 的 `_KIND_TABLE` entry 添加 `variant`（苯环保留名特例）

---

## 第五步: FG-环组合支持

FG-环组合的命名由以下路径承担：

- **苯系保留名**（benzoic/phenol/aniline/benzaldehyde/benzonitrile/benzamide/benzoate）：`chain_engine._KIND_TABLE` 各 entry 的 `variant["benzene"]`（要求 `scaffold_id=="benzene"` 且 `multiplicity==1`）
- **稠环/杂环 FG 收敛**：`express_ring_principal`（`principal_expression.py:256`）收敛为 FG 类别，词干由 `_ring_stem`（`assembler.py:27`）注入
- **环外酸**（`facts.relation == "exocyclic"`）：`assembler._exocyclic_ring_names`（`:78`）→ `cyclohexanecarboxylic acid`；多羧酸（multiplicity≥2）拼 …-di/tricarboxylic acid 且位次必带（P-65.2.2），苯单酸回落 benzoic acid 保留名
- **环外醛**（`facts.relation == "exocyclic"`）：`assembler._exocyclic_ring_names`（`assembler.py:78`）→ `cyclohexanecarbaldehyde`；多醛 → -dicarbaldehyde（P-66.6.1.1.3）；苯单醛仍走 benzaldehyde 保留名
- **环外酰胺**（`facts.relation == "exocyclic"`）：`assembler._exocyclic_ring_names` → `cyclohexanecarboxamide`
- **环外腈**（`facts.relation == "exocyclic"`）：`assembler._exocyclic_ring_names` → `cyclohexanecarbonitrile`
- **环外酰基头**（`facts.relation == "exocyclic"`）：`assembler._exocyclic_ring_names` → `cyclohexanecarbonyl`
- **环外酯**（`facts.relation == "exocyclic"`）：`assembler._exocyclic_ring_names` → `cyclohexanecarboxylate`
- **环内主 FG 的准入**：环内 keto/alcohol/amine 与环外 acid/nitrile 是否可表达由 `ring_expression_policy.supports_ring_expression`（`:41`）按 `naming_class` 白名单判定（见 2.3）

---

## 文件改动清单

| Layer | 文件 | 改动内容 |
|-------|------|---------|
| L1 | `layer1/ring_systems.py` | 特殊环检测（通常无需改动；SSSR + 融合由 `build_ring_systems`:289 承担） |
| L2 | `layer2/ring_scaffold.py` | **主改动**：在 `_TEMPLATES`（`:71`）新增 `{smiles, stem_en, stem_zh, naming_class}`（+ 可选 `locant_prefix`/`prefix_nh_conditional`，作稠合命名零件时再补 `fused`/`fused_stem`/`fused_prefix`）；非默认编号时在同条目加 `standard = (labels, order)` |
| L2 | `layer2/ring_expression_policy.py` | 新 `naming_class` 要支持环内酮/醇/胺或环外酸/腈时，把该 `naming_class` 加进 `_POLICIES`（`:18`）对应策略行 |
| L2 | `layer2/kind_registry.py` | 通常无需改动（`_load_from_scaffold_specs`:83 自动注册） |
| L2 | `layer2/fused_system.py` | 通常无需改动（稠环拆解走 `_select_base`:85 的通用准则） |
| L4 | `layer4/numbering_engine.py` | 通常无需改动（P-14.4 候选管线 + `standard` 固定编号自动适用） |
| L4 | `layer4/fused_numbering.py` | 通常无需改动（无 `standard` 的稠环走 P-25.3.3 外周编号） |
| L5 | `layer5/assembler.py` | `_ring_stem`（`:27`）词干注入 / `_names_for`（`:372`）worker |
| L5 | `layer5/fused_namer.py` | 通常无需改动（未注册稠环的稠合名组装，`fused_parent_names`:126） |
| L5 | `layer5/chain_engine.py` | `_KIND_TABLE`（`:421`）的 `variant`（苯环保留名特例） |
| L5 | `layer5/assembler_prefixes.py` | 前缀构建规则（如需） |
| — | `tests/` | 添加 SMILES 测试用例 |

---

## 工作示例: 新增吡啶（pyridine）支持

以下是 NamePredict 中吡啶支持的实现路径（实际代码），理解它即可掌握新增保留名单杂环的完整流程。

### Step 1 — 环检测

吡啶的六元芳环（5C + 1N）被 SSSR 自动检测。未对 Layer1 做额外改动。

### Step 2 — Layer2 骨架识别

pyridine 走**数据驱动**识别，全部在 `ring_scaffold.py`：

1. **模板注册**：`_TEMPLATES` 的 `pyridine` 条目（`ring_scaffold.py:82`，`stem_en="pyridine"`, `stem_zh="吡啶"`, `naming_class="monohetero"`, `fused=True`, `fused_prefix=("pyrido","吡啶并")`），`_spec_from_template`（`:359`）派生 ScaffoldSpec（retained=True）
2. **子图同构匹配**：`match_retained`（`:410`）用 SMILES 模板做元素标注的子图同构，`resolve_ring_scaffold`（`:477`）标注 `scaffold_id="pyridine"`

### Step 3 — 编号

无需 orienter。`numbering_engine.orient_numbering`（`:316`）的 P-14.4 候选管线对环枚举 2n 个候选；`_narrow_hetero_ring`（`numbering_engine.py:168`）对杂环做元素序窄化（P-22.2.2.1.3：全杂原子最低位次→逐元素序→唑 NH=1，先行于 principal），使单杂环的唯一杂原子/吡啶的 N 定于 1 号位。

### Step 4 — 命名

`assembler._names_for`（`:372`）查 `_KIND_TABLE`（`alkane` entry，`chain_engine.py:438`），`_ring_stem`（`assembler.py:27`）从 parent 的 `stem_en`/`stem_zh` 注入词干 `("pyridine","吡啶")`（`assembler.py:393`）。FG 变体（pyridinecarboxylic 等）由 L2 `express_ring_principal`（`principal_expression.py:256`）收敛为 FG 类别后同样走 chain_engine。

### Step 5 — 测试

SMILES 测试用例示例：
- `c1ccncc1` — pyridine / 吡啶
- `c1ccncc1C(=O)O` — pyridine-2-carboxylic acid（吡啶甲酸 / 吡啶-2-甲酸）
- `c1ccncc1N` — pyridin-2-amine（吡啶-2-胺）

---

## 常见问题

### Q: 有 ScaffoldSpec 但没有骨架识别路径？

`ScaffoldSpec` 仅提供 stem 和编号策略的元数据；分子仍须在 `resolve_ring_scaffold`（`ring_scaffold.py:477`）的识别路径（`get_identity`:383 直查 / `match_retained`:410 子图同构 / `_generic_carbocycle`:464 兜底）中被命中，其骨架才会使用该 scaffold。

### Q: 为什么有些环系没有独立的 .py 文件？

部分简单环系（如 furan, thiophene, pyrrole）不设独立模块——它们经 `ring_scaffold.py` 的 `_TEMPLATES` 批量注册，共享相同的编号策略，仅杂原子与 stem 不同。

### Q: 多环稠合系统如何处理？

多环稠合系统（如 naphthalene, indole）经 `ring_systems.py` 的 ring_system 分组识别（`build_ring_systems`，`:289`），骨架身份由 `resolve_ring_scaffold`（`ring_scaffold.py:477`）解析（`match_retained` 子图同构）。编号方面，L2 注入 `numbering_scaffold` 事实（`{scaffold_id, labels, relative_stereo}`，由 `numbering_scaffold_facts`:398 生成、`kind_registry._attach_numbering_scaffold`:33 写入），L4 校验其存在性；`labels` 与链长不符时 L4 退回链序号编号。

### Q: ScaffoldSpec 和 KindMeta 的关系？

`ScaffoldSpec` 是 stem 和编号的单一权威来源（Single Authority）。`kind_registry.py` 在 bootstrap 的唯一一步（`_load_from_scaffold_specs`，`:83`）读取所有 `ScaffoldSpec` 并转换为 `KindMeta` 注册（`_bootstrap`，`:94`）。ScaffoldSpec 未给出 stem 的骨架跳过注册。

> **源:** `src/namepredict/layer2/kind_registry.py:83`

---

## 相关文档

- [[architecture/layer2-parent-selector]] — Layer2 母体选择器完整架构
- [[architecture/layer4-numbering]] — Layer4 编号引擎（P-14.4 候选管线）
- [[architecture/layer5-name-assembly]] — Layer5 名称组装
- [[concepts/functional-group-priority]] — FG 优先级与 IUPAC P-44 规则
