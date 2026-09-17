# 指南：新增环系 (Adding a Ring System)

适用范围：`src/namepredict` 命名引擎；以 `layer2/ring_scaffold.py` 的 `_TEMPLATES` 为唯一事实来源，符号名与当前源码一致。

本指南以**含硫饱和六元杂环（thiane 类）**为主线，走完「识别 → 注册 → 固定编号 → 命名 → 与 FG 组合」五步。

## 落地路径总览

| 步骤 | 承担模块 | 新增 thiane 时的动作 |
|---|---|---|
| 1 环检测与划分 | `layer1/ring_systems.py` | 无 |
| 2 骨架注册 | `layer2/ring_scaffold.py` | `_TEMPLATES` 加一条 |
| 3 固定编号 | `layer4/numbering_engine.py`、`fused_numbering.py`、`fused_orientation.py`、`ring_geometry.py` | 无（走 P-14.4 候选枚举） |
| 4 指示氢 | `layer4/numbering_engine.py` 的 `_nh_sites` | 无 |
| 5 名称组装 | `layer5/chain_engine.py`、`fused_namer.py` | 无 |
| 6 与 FG 组合 | `layer2/principal_expression.py`、`parent_skeleton.py`、`ring_expression_policy.py` | 无（`monohetero` 已在策略内） |

核心原则：**scaffold × kind × 数量正交**。新环只加模板条目（外加可选的 `standard` 与 `fused_prefix`），不新增 kind，也不枚举「环 + 某 FG」的组合 kind。

## 第 1 步 L1：环检测与划分

`layer1/ring_systems.py` 提供跨层公共原语：

| 符号 | 作用 |
|---|---|
| `sssr_rings(mol)` | 统一环访问器，返回 SSSR 最小环原子序列（经 `tools/memo.by_mol` 缓存） |
| `kekulized(mol)` | 返回 Kekulé 化并清芳香标志的副本；**失败返回 `None`**，调用方须退回原 mol |
| `build_ring_systems(mol)` | 环系划分入口：`_ring_pairs` 单遍扫环对（共享 ≥2 原子记为稠合边 `fusion_edges`，恰好 1 个记为螺环对），`_components` 只用**稠合边**做并查集连通分量 |

每个环系 dict 的字段：`atom_ids` / `sssr_indices` / `fusion_edges` / `n_rings` / `n_atoms` / `hetero_atoms` / `is_aromatic_mancude` / `topology`。

其中 `sssr_indices` 与 `fusion_edges` 是 L2 稠环拆解与 L4 稠环编号的输入。新增单环或普通稠环在 L1 无需改动。

## 第 2 步 L2：骨架注册（主要工作量）

### 2.1 唯一注册表：`_TEMPLATES`

`layer2/ring_scaffold.py` 的 `_TEMPLATES` 是保留母体 SMILES 模板的唯一注册表；查询入口 `get_spec(scaffold_id)` 返回派生出的 `ScaffoldSpec`，`all_specs()` 返回全部。

thiane 条目：

```python
"thiane": {"smiles": "C1CCSCC1", "stem_en": "thiane", "stem_zh": "四氢噻喃",
           "naming_class": "monohetero", "fused": True,
           "fused_prefix": ("thiopyrano", "噻喃并")},
```

字段清单：

| 字段 | 必填 | 作用 |
|---|---|---|
| `smiles` | 是 | 模板查询分子；import 期由此构建子图同构查询与元素签名 |
| `stem_en` / `stem_zh` | 是 | 双语词干，经 `ScaffoldSpec` 注入 L5 |
| `naming_class` | 是 | 命名类，决定 FG 准入策略与稠环拆解行为 |
| `fused` | 否 | 才可作稠合命名的**母体组分**；`match_fusion_component` 的 `mancude_only` 按它过滤 |
| `fused_prefix` | 否 | 作**附加组分**时的保留稠合前缀 `(en, zh)`，由 `retained_fusion_prefix` 消费 |
| `fused_stem` | 否 | 稠合前取用的词干（去指示氢），`component_stem` 优先取它 |
| `locant_prefix` | 否 | 词干内位次前缀（`"1H-"`、`"1,3-"` 等） |
| `prefix_nh_conditional` | 否 | `True` 时仅当环内存在未取代 NH 才注入 `locant_prefix` |
| `standard` | 否 | 固定编号 `(labels, order)`，见 2.3 |

### 2.2 `fused_prefix` 的语义（新增饱和杂环必读）

稠合组分的名不能用饱和环名去尾加「并」：须取对应的 **mancude 名**（P-25.3.1.2.1 / P-25.3.1.2.3）。`fused_prefix` 缺失时 `fused_namer._component_prefix` 退回「去尾 `e` 加 `o`」的通用规则，把 `thiane` 拼成 `thiano`、中文拼成「四氢噻喃并」。

现有取值：

| scaffold_id | `fused_prefix` |
|---|---|
| `oxane` | `("pyrano", "吡喃并")` |
| `thiolane` | `("thieno", "噻吩并")` |
| `thiane` | `("thiopyrano", "噻喃并")` |

判定：新增饱和杂环若可能作稠合组分，必须同时给 `fused=True` 与 `fused_prefix`。

### 2.3 固定编号登记：`standard`

单杂环默认走 P-14.4 候选枚举，无需 `standard`。只有当保留编号不能由枚举得到时，才在**同一条目内**加 `(labels, order)`：`order` 是模板原子按 locant 顺序的下标排列，`labels` 是对应标签元组，两者长度都等于模板原子数。import 期由 `_validate_standard_fields` 校验，不符即 `ValueError` 中止导入。

`locant_prefix` 与 `standard` 是两件不同的事：前者是词干里的位次前缀（如 `dioxolane` 的 `"1,3-"`），后者是模板原子到 locant 的映射。新增环的词干本身若含加氢或位次前缀（`4,5-dihydro-...`），须同步 `locant_prefix`。

### 2.4 碳环稠合前缀：`_FUSION_CARBOCYCLES`

单环烃附加组分表（P-25.3.2.2.1），字段 `smiles` / `prefix_en` / `prefix_zh`，由 `fusion_carbocycle_prefix` 与 `match_fusion_carbocycle` 消费，成员为 `cyclopropane` 至 `cyclooctane`。

注意六元碳环的稠合前缀是保留前缀 `benzo` / `苯并`，而不是 `cyclohexa` / `环己并`。

这些成员只作稠环拆解的附加零件，不入 `_TEMPLATES`。

## 第 3 步 L4：固定编号

### 3.1 单环：P-14.4 候选枚举

`numbering_engine.orient_numbering` 逐级尝试，`_fixed_numbering` 与 `_fused_numbering` 都不适用时走通用枚举：杂环进 `_narrow_hetero_ring`，碳环进 `_ring_cands`。收窄原语是跨层公共 API：

| 原语 | 语义 |
|---|---|
| `narrow(cands, key_fn, *, reverse=False, skip_none=False)` | 按 key 取最优值收窄；`skip_none=True` 表示特征全缺时规则不适用、候选原样返回 |
| `narrow_by_senior(cands, key_fn, heteros, by_z, *, skip_none=False)` | 先取全杂原子集最低位次，再按 `constants.P145_SENIOR` 逐元素收窄 |

L2 的稠环拆解 `fused_system._select_base` 复用 `narrow`，因此这两条原语是跨层公共面。

> `constants.P25_SENIOR` 与 `constants.P145_SENIOR` 同源不同序：前者用于「选哪个组分当母体」（`_select_base` 准则 (a)），后者用于「哪个杂原子得低位次」，勿混用。

### 3.2 稠环：定向 → 几何 → 分层编号

- 定向：`fused_orientation.preferred_orientations(mol, rings, fusion_edges)` 枚举全部平局取向（含镜像），返回 `Orientation`（字段 `row` / `coords` / `quad`）。
- 几何：`ring_geometry.py` 提供模板与合法性判定——`regular_polygon` / `RING_TEMPLATES`、`ring_shape_template`（5 元与 7 元环的变形模板）、`rigid_fit`、`clip_polygon`、`overlap_area`，阈值常量 `DEFORM_MAX` / `OVERLAP_FRAC`。
- 编号：`fused_numbering.number_fused_system(mol, rings, coords, sub_layers=None, alpha_subs=None)` 对每个取向的候选分层收窄：**(a)(b) 杂原子 → (c) 稠合碳 → (d) 稠合杂原子 → `sub_layers` → 字母序 → CIP**；辅助原语 `fused_atoms` 给出出现在 ≥2 环的稠合原子。

`_fused_numbering` 组装的 `sub_layers` 顺序即分层优先级，传统编号骨架在 `constants.TRADITIONAL_NUMBERING_IDS` 内直接交由 `standard` 固定编号。

### 3.3 指示氢

`numbering_engine._nh_sites(heteros, mol)` 给出 P-22.2.2.1.4 的指示氢候选位（环内可带 H 的 N），由 `_narrow_hetero_ring` 收窄；`float_hetero=True` 时跳过该层。稠环侧由哨兵 `fused_numbering.INDICATED_H` 承担。

稠环完整层级顺序：

**principal 特征基团 → 指示氢 → 加氢位次（`hydro_atoms`）→ 取代基位次集 → 字母序平局。**

## 第 4 步 L5：名称组装

| 环节 | 承担符号 |
|---|---|
| 词干 | `_TEMPLATES` 的 `stem_en` / `stem_zh` 经 `ScaffoldSpec` 与 `kind_registry._load_from_scaffold_specs` 注册为 `KindMeta`，再由 `pack_parent_stem` 写入 parent，`assembler._names_for` 把 `parent["stem_en"]` / `parent["stem_zh"]` 注入 `chain_engine._KIND_TABLE` 的词条 |
| 苯单取代保留名 | `chain_engine._BENZENE_RETAINED`（key 为 kind） |
| 环外主基后缀 | `constants.EXO_RING_SUF` 经 `chain_engine._exo_ring_spec` 改写；主基位次是否必带由 `_ring_prefix_located` 判定 |
| 稠合成名 | `fused_namer.fused_parent_names`（由 `assembler._ensure_fused_stem` 调用）→ `_collect_attached` / `_fused_one` / `_component_prefix`；整名命中 `constants.RETAINED_FUSION_ALIASES` 时替换为保留名 |
| 前缀构建 | `layer5/assembler_prefixes.py` |

新增 thiane 只需 `stem_en` / `stem_zh` 两个字段，L5 无需改动。

## 第 5 步 与 FG 组合

`principal_expression.express_ring_principal` 决定环 + 主 FG 的 kind：命中 FG 类别时经 `_ring_kind` / `_chain_kind` 收敛为 FG 类别 kind，词干仍由 scaffold 承载；否则回落 `_resolved_ring_kind` 的结构 kind。

环骨架是否被接受有两道闸：

1. **骨架筛选谓词**（`parent_skeleton.py`）：入口 `select_principal_skeletons` 依次施加 `keep_max_principal_coverage` → `keep_p44_1_2`（拓扑不唯一时先经 `keep_senior_atom`）→ `keep_p44_2` 或 `keep_p44_3` → `keep_p44_4_unsaturation`。
2. **环表达策略**：`ring_expression_policy.supports_ring_expression` 按 `naming_class` + FG 类别 + 关系（`in_skeleton` / `exocyclic`）白名单判定，结果写入 parent 的 `typed_ring_expression_supported`；不支持者被 `parent_select._unsupported_typed_ring` 整体拦截，候选数为 0 而非报错。

新 scaffold 若引入新 `naming_class` 且要支持环内酮/醇/胺或环外酸/腈/酰基头，必须把该 `naming_class` 加进 `ring_expression_policy._POLICIES` 的对应策略行。

## 改动清单

| 改动点 | 文件 | 是否派生自动 |
|---|---|---|
| 环系划分 | `layer1/ring_systems.py` | 是（`build_ring_systems` 通用） |
| 模板条目 `smiles` / `stem_en` / `stem_zh` / `naming_class` | `layer2/ring_scaffold.py` 的 `_TEMPLATES` | **否，唯一的手写改动** |
| `ScaffoldSpec` / `ScaffoldIdentity` | `ring_scaffold._spec_from_template` | 是 |
| 固定编号映射 `standard` | 同条目内 | 否（仅保留编号非枚举可得的环） |
| 稠合前缀 `fused_prefix` | 同条目内 | 否（仅可能作稠合组分的饱和杂环） |
| 单环烃附加组分 | `ring_scaffold._FUSION_CARBOCYCLES` | 否 |
| 词干注册 | `layer2/kind_registry.py` 的 `_load_from_scaffold_specs` | 是 |
| 稠环拆解 | `layer2/fused_system.py` | 是（`_select_base` 准则 (a)–(j) 通用） |
| L4 单环编号 | `layer4/numbering_engine.py` | 是（P-14.4 枚举） |
| L4 稠环编号 | `fused_numbering.py`、`fused_orientation.py`、`ring_geometry.py` | 是 |
| L4 指示氢 | `numbering_engine._nh_sites` | 是 |
| L5 词干 / 后缀 / 稠合名 | `chain_engine.py`、`fused_namer.py` | 是 |
| 苯环专属保留名 | `chain_engine._BENZENE_RETAINED` | 否（仅苯环保留名） |
| 环外主基后缀 | `constants.EXO_RING_SUF` | 否（仅新 FG 类别） |
| FG 准入 | `layer2/ring_expression_policy.py` 的 `_POLICIES` | 否（仅新 `naming_class`） |
| 测试用例 | `tests/` | 否 |

## 常见陷阱

1. **漏 `fused_prefix`**：饱和环作稠合组分时走「去尾 e 加 o」兜底，把「四氢噻吩并」「四氢噻喃并」拼进稠合名。`thiane` 作附加组分必须给 `("thiopyrano", "噻喃并")`。
2. **六元碳环误用 `cyclohexa`**：`_FUSION_CARBOCYCLES` 里六元环的前缀是保留前缀 `benzo` / `苯并`。
3. **模板 SMILES 与真实环的取代模式不匹配**：`match_retained` / `_match_with_map` 要求模板原子集与骨架原子集**精确相等**，失配则 `get_spec` 取不到，环退化为通用命名（未注册的单杂环直接失败）。
4. **`standard` 与 `locant_prefix` 不同步**：位次前缀挂在词干里而 `standard` 另算映射时，两者须对齐，否则位次整体错位。此外 `labels` 与链长不符时，`numbering_engine._component_labels` 静默回落成 `1..n` 纯数字，产出错名而不报错。
5. **指示氢未收窄**：`_nh_sites` 的结果不进 `_narrow_hetero_ring` 的收窄层时，杂环编号会在同分位次里任选，`1H-` / `10H-` 前缀随之出错。
6. **`naming_class` 未进 `_POLICIES`**：环内酮/醇/胺与环外酸/腈的 typed 表达被 `parent_select._unsupported_typed_ring` 整体拦截，产出空候选而不是报错。
7. **`P25_SENIOR` 与 `P145_SENIOR` 混用**：前者选稠环母体组分，后者定杂原子低位次。
8. **单环烃附加组分不要写进 `_TEMPLATES`**：它们只作稠合零件；入表会让单环骨架解析成保留名，破坏单环通用路径。
9. **`kekulized` 失败返回 `None`**：调用方须降级退回原 mol，否则加氢位与指示氢判定会静默退化。

## 验收用例（thiane 类）

| SMILES | EN | ZH |
|---|---|---|
| `C1CCSCC1` | thiane | 四氢噻喃 |
| `C1CCSCC1O` | thian-3-ol | 四氢噻喃-3-醇 |
| `C1CCSCC1N` | thian-3-amine | 四氢噻喃-3-胺 |
| `C1CCSCC1C(=O)O` | thiane-3-carboxylic acid | 四氢噻喃-3-羧酸 |
| `C1CCSCC1C=O` | thiane-3-carbaldehyde | 四氢噻喃-3-甲醛 |
| `C1CCSCC1C#N` | thiane-3-carbonitrile | 四氢噻喃-3-甲腈 |
| `c1cc2c(cn1)SCCC2` | 3,4-dihydro-2H-thiopyrano[2,3-c]pyridine | 3,4-二氢-2H-噻喃并[2,3-c]吡啶 |
| `c1cc2c(cn1)SCC2` | 2,3-dihydrothieno[2,3-c]pyridine | 2,3-二氢噻吩并[2,3-c]吡啶 |

最后两行分别验证 `thiane → thiopyrano` 与 `thiolane → thieno` 的 `fused_prefix` 生效。

## 相关文档

- [[architecture/layer1-analyzer]]
- [[architecture/layer2-parent-selector]]
- [[architecture/layer4-numbering]]
- [[architecture/layer5-name-assembly]]
- [[reference/core-data-contracts]]
- [[guides/adding-new-functional-group]]
