# 指南：新增环系 (Adding a Ring System)

适用范围：`src/namepredict` 命名引擎；以 `layer2/ring_scaffold.py` 的 `_TEMPLATES` 为唯一事实来源，符号名与源码一致。

本指南以**含硫饱和六元杂环（thiane 类）**为主线，走完「识别 → 注册 → 固定编号 → 命名 → 与 FG 组合」五步。

## 落地路径总览

| 步骤 | 承担模块 | 新增 thiane 时的动作 |
|---|---|---|
| 1 环检测与划分 | `layer1/ring_systems.py` | 无 |
| 2 骨架注册 | `layer2/ring_scaffold.py` | `_TEMPLATES` 加一条 |
| 3 固定编号 | `layer4/`：`numbering_engine`、`fused_numbering`、`fused_orientation`、`ring_geometry` | 无（走 P-14.4 候选枚举） |
| 4 指示氢 | `layer4/numbering_engine._nh_sites` | 无 |
| 5 名称组装 | `layer5/chain_engine.py`、`fused_namer.py` | 无 |
| 6 与 FG 组合 | `layer2/`：`principal_expression`、`parent_skeleton`、`ring_expression_policy` | 无（`monohetero` 已在策略内） |

核心原则：**scaffold × kind × 数量正交**——新环只加模板条目，不新增 kind，也不枚举「环 + 某 FG」的组合 kind。

识别有四条路径，只有第一条需注册模板：**保留模板**（`match_retained` → `_TEMPLATES`）、**P-25 稠合**（`decompose_fused_system`）、**P-23 桥环**（`try_bridged_scaffold`）、**生成式词干**。后三条免注册，见 2.5。

## 第 1 步 L1：环检测与划分

`layer1/ring_systems.py` 提供跨层公共原语：`sssr_rings(mol)` 返回 SSSR 最小环原子序列（经 `tools/memo.by_mol` 缓存）；`kekulized(mol)` 返回 Kekulé 化并清芳香标志的副本，**失败返回 `None`**，调用方须退回原 mol；`build_ring_systems(mol)` 是环系划分入口——`_ring_pairs` 扫环对（共享 ≥2 原子记稠合边 `fusion_edges`，恰好 1 个记螺环对），`_components` 只用**稠合边**做并查集分量。

环系 dict 字段：`atom_ids` / `sssr_indices` / `fusion_edges` / `n_rings` / `n_atoms` / `hetero_atoms` / `is_aromatic_mancude` / `topology`；`sssr_indices` 与 `fusion_edges` 是 L2 拆解与 L4 编号的输入。新增单环或普通稠环在 L1 无需改动。

## 第 2 步 L2：骨架注册（主要工作量）

### 2.1 唯一注册表：`_TEMPLATES`

`_TEMPLATES` 是保留母体 SMILES 模板的唯一注册表；`get_spec(scaffold_id)` 返回派生的 `ScaffoldSpec`，`all_specs()` 返回全部。thiane 条目：

```python
"thiane": {"smiles": "C1CCSCC1", "stem_en": "thiane", "stem_zh": "四氢噻喃",
           "naming_class": "monohetero", "fused": True,
           "fused_prefix": ("thiopyrano", "噻喃并")},
```

必填字段：`smiles`（模板查询分子；import 期由此构建子图同构查询与元素签名）、`stem_en` / `stem_zh`（双语词干，经 `ScaffoldSpec` 注入 L5）、`naming_class`（决定 FG 准入策略与稠环拆解行为）。

可选字段：`fused`（才可作稠合命名的**母体组分**）、`fused_prefix`（作**附加组分**时的保留稠合前缀，由 `retained_fusion_prefix` 消费）、`fused_stem`（稠合前取用的词干；缺省时 `component_stem` 走 `_bracketed_stem`，见 2.2）、`locant_prefix`（词干内位次前缀）、`prefix_nh_conditional`（仅环内有未取代 NH 时才注入）、`standard`（固定编号，见 2.3）。

按保留名登记的条目有 `picene`、`indolizine`、`pyrrolizine`、`pyrazolidine`、`thiadiazolidine124`；`azepane` 另带 `fused_stem`（`("azepine", "氮杂卓")`），作稠合母体时取 mancude 词干。

### 2.2 `fused_prefix` 的语义（新增饱和杂环必读）

稠合组分的名不能用饱和环名去尾加「并」：须取对应的 **mancude 名**（P-25.3.1.2.1 / P-25.3.1.2.3）。`fused_prefix` 缺失时 `fused_namer._component_prefix` 退回「去尾 `e` 加 `o`」，把 `thiane` 拼成 `thiano`、中文拼成「四氢噻喃并」。现有取值：`oxane → pyrano`、`thiolane → thieno`、`thiane → thiopyrano`；可能作稠合组分的饱和杂环须同时给 `fused=True` 与 `fused_prefix`。

组分名（`component_stem`）还会经 `_bracketed_stem` 定形（P-25.3.5）：`locant_prefix` 为纯数字位次（`1-` / `1,3-` / `1,2,4-`）且 `stem_en` 不以数字开头时写作方括号形式（`[1]benzofuran` / `[1,3]benzothiazole`）；`1H-` 型前缀与词干已含位次者（`1,2-oxazole`）不加。`RETAINED_FUSION_ALIASES` 的键须按此形态书写。

### 2.3 固定编号登记：`standard`

单杂环默认走 P-14.4 候选枚举，无需 `standard`；只有保留编号不能由枚举得到时，才在**同一条目内**加 `(labels, order)`（`order` 是模板原子按 locant 顺序的下标排列，`labels` 是对应标签元组，长度都等于模板原子数）。import 期由 `_validate_standard_fields` 校验，不符即 `ValueError`。

`locant_prefix` 是词干里的位次前缀（如 `dioxolane` 的 `"1,3-"`），`standard` 是模板原子到 locant 的映射，两者是两件事；词干本身含加氢或位次前缀（`4,5-dihydro-...`）的环须同步 `locant_prefix`。

词干已含字面位次的部分不饱和杂环（`dihydrofuran` / `dihydropyran` / `dihydropyrrole` / `dihydroimidazole` / `dihydrothiazole`），字面位次即固定编号，条目内登记 `standard`（P-14.4(a)/(b)），否则枚举给出别的取向。

### 2.4 碳环稠合前缀：`_FUSION_CARBOCYCLES`

单环烃附加组分表（P-25.3.2.2.1），字段 `smiles` / `prefix_en` / `prefix_zh`，成员 `cyclopropane` 至 `cyclooctane`，由 `fusion_carbocycle_prefix` 与 `match_fusion_carbocycle` 消费。六元碳环的稠合前缀是保留前缀 `benzo` / `苯并`（非 `cyclohexa`）；这些成员只作稠环拆解的附加零件，不入 `_TEMPLATES`。

### 2.5 免注册的路径（P-23）

**桥环**：`bridged_system.decompose_bridged_system(info, system)` 按扩展 von Baeyer 规则自动拆解——桥头取环系内度 ≥3 的原子，去桥头分量逐块成桥，桥头直键记为零原子桥（P-23.1.2）；候选经 P-23.2.1～P-23.2.6.2.5 收窄后产出 `BridgedNode`，`scaffold_identity()` 给出 `ScaffoldIdentity("bridged", "bridged", n_rings, ring)`。门槛（`try_bridged_scaffold`，环系 ≥2 个 SSSR 环时才调用）：「稠合命名法不适用」，或「`fused_hetero` 且稠环树盖不住全部环 / `_p25_names_ok` 拼不出稠合词干」。身份写在 `bridged_node` / `bridged_nodes`。

**生成式词干**：`assembler._ensure_generated_stem` 对词干为空、>10 元、`scaffold_id == "carbocycle"`、`kind == "alkane"` 的环，取「`prefix_from_chain` 的 'a' 前缀 + 环烷」（`1-azacyclotridecane`）；≤10 元环仍走 Hantzsch-Widman 保留名。无模板命中的单环由 `_generic_carbocycle` 归为 `carbocycle` 身份。

## 第 3 步 L4：固定编号

### 3.1 单环：P-14.4 候选枚举

`numbering_engine.orient_numbering` 逐级尝试：`_fixed_numbering` 与 `_fused_numbering` 都不适用时，杂环进 `_narrow_hetero_ring`、碳环进 `_ring_cands`；桥环（带 `bridged_nodes`）在入口即交 `bridged_numbering`。收窄原语 `narrow(cands, key_fn, *, reverse=False, skip_none=False)`（按 key 取最优值；`skip_none=True` 表示特征全缺时规则不适用）与 `narrow_by_senior(...)`（先取全杂原子集最低位次，再按 `constants.P145_SENIOR` 逐元素）是跨层公共 API，L2 的 `_select_base`、`narrow_candidates` 与 L4 桥环裁决都复用。

> `constants.P25_SENIOR` 与 `P145_SENIOR` 同源不同序：前者选稠环母体组分（`_select_base` 准则 (a)），后者定杂原子低位次，勿混用。

### 3.2 稠环：定向 → 几何 → 分层编号

- 定向：`fused_orientation.preferred_orientations(mol, rings, fusion_edges)` 枚举全部平局取向（含镜像）；`_layout` 摆放行外环时按「已有已摆放邻居」反复递推，直到无可推进，返回 `Orientation`（`row` / `coords` / `quad`）。
- 几何：`ring_geometry.py` 给模板与合法性判定——`regular_polygon`、`RING_TEMPLATES`（覆盖 3–98 元环）、`ring_shape_template`（5/7 元变形模板）、`rigid_fit`、`overlap_area`，阈值 `DEFORM_MAX` / `OVERLAP_FRAC`。
- 编号：`fused_numbering.number_fused_system(...)` 对每个取向分层收窄：**(a)(b) 杂原子 → (c) 稠合碳 → (d) 稠合杂原子 → `sub_layers` → 字母序 → CIP**；`fused_atoms` 给出出现在 ≥2 环的稠合原子。

`_fused_numbering` 组装的 `sub_layers` 顺序即分层优先级，传统编号骨架（`constants.TRADITIONAL_NUMBERING_IDS`）交由 `standard` 固定编号。

桥环编号不由模板注册决定：parent 带 `bridged_node` 时 `_fused_numbering` 返回 `None`，编号改由 `layer4/bridged_numbering.bridged_numbering` 在 L2 下传的并列候选间裁决——按 P-23.3.2 收窄杂原子位次后依次施加 P-14.4(c)(f)(g)，候选渲染特征全同才取首个，否则不可判定返回 `None`；选中者写回 `parent["bridged_node"]`。

### 3.3 指示氢

`numbering_engine._nh_sites(heteros, mol)` 给出 P-22.2.2.1.4 的指示氢候选位（环内可带 H 的 N），由 `_narrow_hetero_ring` 收窄；`float_hetero=True` 时跳过该层，稠环侧由 `fused_numbering.INDICATED_H` 哨兵承担。稠环层级顺序：**principal 特征基团 → 指示氢 → 加氢位次（`hydro_atoms`）→ 取代基位次集 → 字母序平局。**

## 第 4 步 L5：名称组装

| 环节 | 承担符号 |
|---|---|
| 词干注入 | `assembler._ensure_parent_stem`：桥环走 `_ensure_bridged_stem` → `bridged_namer.bridged_parent_names`（环数词头 + 描述符 + 烃词干 + `prefix_from_chain` 的 'a' 前缀，另给 `stem_bare_en/zh`）；否则 `_ensure_fused_stem` → `_ensure_generated_stem` |
| 词干落地 | `stem_en` / `stem_zh` 经 `ScaffoldSpec` 与 `kind_registry._load_from_scaffold_specs` 注册为 `KindMeta`，由 `pack_parent_stem` 写入 parent，`_names_for` 注入 `chain_engine._KIND_TABLE` 词条 |
| 不饱和渲染 | `_names_for` 对带 `stem_bare_en` 且 kind 为 `bridged` / `alkane` 的母体改走链引擎，en/yne 位次由此保留；`_chain_unsat_fields` 对 `carbocycle` / `bridged` 身份补回 Kekulé 环内 C=C |
| 苯单取代保留名 | `chain_engine._BENZENE_RETAINED`（key 为 kind） |
| 环外主基后缀 | `constants.EXO_RING_SUF` 经 `chain_engine._exo_ring_spec` 改写；位次是否必带由 `_ring_prefix_located` 判定 |
| 稠合成名 | `fused_namer.fused_parent_names`（`_ensure_fused_stem` 调用）→ `_collect_attached` / `_fused_one` / `_component_prefix`；整名命中 `RETAINED_FUSION_ALIASES` 时替换为保留名 |
| 前缀构建 | `layer5/assembler_prefixes.py` |

新增 thiane 只需 `stem_en` / `stem_zh` 两个字段，L5 无需改动。

## 第 5 步 与 FG 组合

`principal_expression.express_ring_principal` 决定环 + 主 FG 的 kind：命中 FG 类别时经 `_ring_kind` / `_chain_kind` 收敛为 FG 类别 kind，词干仍由 scaffold 承载；否则回落 `_resolved_ring_kind` 的结构 kind（`bridged` 即由此给出）。scaffold、稠环树与桥环候选由 `_ring_scaffold_and_nodes` 一次算好并写入 parent。

环骨架是否被接受有两道闸：

1. **骨架筛选谓词**：`parent_skeleton.select_principal_skeletons` 依次施加 `keep_max_principal_coverage` → `keep_p44_1_2`（拓扑不唯一时先经 `keep_senior_atom`）→ `keep_p44_2` 或 `keep_p44_3` → `keep_p44_4_unsaturation`。
2. **环表达策略**：`ring_expression_policy.supports_ring_expression` 按 `naming_class` + FG 类别 + 关系（`in_skeleton` / `exocyclic`）白名单判定，结果写入 `typed_ring_expression_supported`；`parent_select._unsupported_typed_ring` 依 `scaffold_identity` 与该字段判环酮是否被排除。

新 `naming_class` 若要支持环内酮/醇/胺或环外酸/腈，须加进 `_POLICIES` 的对应策略行。

## 改动清单

| 改动点 | 文件 | 是否派生自动 |
|---|---|---|
| 环系划分 / 词干注册 | `layer1/ring_systems.py`、`kind_registry._load_from_scaffold_specs` | 是 |
| 模板条目 `smiles` / `stem_en` / `stem_zh` / `naming_class` | `layer2/ring_scaffold.py` 的 `_TEMPLATES` | **否，唯一的手写改动** |
| `ScaffoldSpec` / `ScaffoldIdentity` | `ring_scaffold._spec_from_template` | 是 |
| `standard` / `fused_prefix` / `fused_stem` / `_FUSION_CARBOCYCLES` | 同条目内 | 否（仅枚举得不到的编号 / 可作稠合组分的环） |
| 稠环 / 桥环拆解 | `layer2/fused_system.py`、`layer2/bridged_system.py` | 是（`_select_base` 准则 (a)–(j) 与 P-23 规则通用） |
| L4 单环 / 稠环编号 | `numbering_engine.py`、`fused_numbering.py`、`fused_orientation.py`、`ring_geometry.py` | 是（P-14.4 枚举 / 分层收窄） |
| L4 桥环编号裁决 | `layer4/bridged_numbering.py` | 是 |
| L4 指示氢 | `numbering_engine._nh_sites` | 是 |
| L5 词干 / 后缀 / 稠合名 | `chain_engine.py`、`fused_namer.py` | 是 |
| L5 桥环 / 生成式词干 | `assembler.py`、`bridged_namer.py`、`skeleton_replacement.py` | 是 |
| 苯环保留名 / 环外主基后缀 / FG 准入 | `_BENZENE_RETAINED`、`constants.EXO_RING_SUF`、`ring_expression_policy._POLICIES` | 否（仅苯环 / 环外表达 / 新 `naming_class`） |
| 测试用例 | `tests/` | 否 |

## 常见陷阱

1. **诱导覆盖未校验**：模板匹配要求原子集诱导子图键数等于模板键数（`_is_induced_match`），否则真子图模板（萘嵌进三环骨架）会静默丢环；氢化骨架比对同样施加。
2. **孤立饱和环配上芳香名**：`_isolated_saturated_ring` 判「不与他环稠合、Kekulé 视图环内全单键」的单环，此类环只匹配全单键模板（P-31.2）。
3. **漏 `fused_prefix`**：饱和环作稠合组分时走「去尾 e 加 o」兜底；六元碳环的稠合前缀是 `benzo` / `苯并`，不是 `cyclohexa`。
4. **模板与真实环的取代模式不匹配**：`match_retained` / `_match_with_map` 要求模板原子集与骨架原子集**精确相等**，失配即退化为通用命名。
5. **`standard` 与 `locant_prefix` 不同步**：两者须对齐，否则位次整体错位；`labels` 与链长不符时 `_component_labels` 静默回落成 `1..n`。
6. **指示氢未收窄**：`_nh_sites` 的结果不进 `_narrow_hetero_ring` 的收窄层时，`1H-` / `10H-` 前缀随之出错。
7. **`naming_class` 未进 `_POLICIES`**：环内酮/醇/胺与环外酸/腈的 typed 表达被 `_unsupported_typed_ring` 拦截。
8. **`P25_SENIOR` 与 `P145_SENIOR` 混用**：前者选稠环母体组分，后者定杂原子低位次。单环烃附加组分不要写进 `_TEMPLATES`，入表会让单环骨架解析成保留名。

## 验收用例（thiane 类）

| SMILES | EN | ZH |
|---|---|---|
| `C1CCSCC1` | thiane | 四氢噻喃 |
| `C1CCSCC1O` | thian-3-ol | 四氢噻喃-3-醇 |
| `C1CCSCC1C(=O)O` | thiane-3-carboxylic acid | 四氢噻喃-3-羧酸 |
| `c1cc2c(cn1)SCCC2` | 3,4-dihydro-2H-thiopyrano[2,3-c]pyridine | 3,4-二氢-2H-噻喃并[2,3-c]吡啶 |
| `c1cc2c(cn1)SCC2` | 2,3-dihydrothieno[2,3-c]pyridine | 2,3-二氢噻吩并[2,3-c]吡啶 |

后两行验证 `thiane → thiopyrano` 与 `thiolane → thieno` 的 `fused_prefix` 生效。

## 相关文档

- [[architecture/layer1-analyzer]]
- [[architecture/layer2-parent-selector]]
- [[architecture/layer4-numbering]]
- [[architecture/layer5-name-assembly]]
- [[reference/core-data-contracts]]
- [[guides/adding-new-functional-group]]
