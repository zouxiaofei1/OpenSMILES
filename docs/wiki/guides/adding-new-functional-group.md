# 添加新官能团 (Adding a New Functional Group)

> **扩展指南** | 关联: [[architecture/overview]], [[architecture/layer1-analyzer]], [[architecture/layer2-parent-selector]], [[architecture/layer3-substituents]], [[architecture/layer4-numbering]], [[architecture/layer5-name-assembly]], [[concepts/functional-group-priority]]

---

## 概述

本指南详细说明如何向 NamePredict 命名管线添加一个新的 functional group（官能团，FG）类型。整个流程跨越 L1（检测与登记）、L2（母体选择）、L5（名称组装）三层；**L4 通常无需改动**（`numbering_engine.orient_numbering` 的 P-14.4 候选管线自动适用，位次记录 `locant_calc._FG_LOCANTS` 由 `FgSpec.locant_kind` 派生）。

我们以 **thiol（硫醇，R-SH）** 作为工作示例——它是 `FG_SPECS` 中的一条 `FgSpec`（`p41=17`、`path=(2,)`、`chain=True`、`multi=True`、`locant_kind="sh"`，`fg_registry.py:57`），代表"检测于 analyzer.py 内联 + `chain=True` 的 SUFFIX 母体 + `_KIND_TABLE` 命名"的标准路径。

> 新增 FG 的核心是三步：**L1 检测 → L1 `FG_SPECS` 登记 `FgSpec` → L5 `_KIND_TABLE` 添加 `_Chain` spec**。下游表（`PRINCIPAL_REGISTRY`、`_CHAIN_FG`/`_MULTI_FG`、`locant_calc._FG_LOCANTS`、`assembler_prefixes._KEEP_LOCANT_KINDS`、`stereo._RS_KINDS`）全部由 `FG_SPECS` 派生，无需逐个登记；也无 `_no_fgs` 互斥检查、L4 orienter 注册、组合 kind。

在开始之前，建议先阅读 [[concepts/functional-group-priority]] 了解 FG 优先级体系，以及 [[architecture/overview]] 了解 6 层管线架构。

---

## 步骤 1: Layer 1 — 检测 (Detection)

Layer 1 的职责是从 RDKit `Mol` 对象中检测官能团，输出结构化的 atom index 条目，并汇总为带类型的 `FunctionalGroupInventory`。

### 检测位置

- **简单 FG**（如 thiol, alcohol, ketone）：直接在 `analyzer.py` 中内联实现，通常只需要 3 个函数：候选原子判断、条目构建、收集函数（参照 `_thiol_entries`，`analyzer.py:222`）。
- **较复杂 FG**（如 acyl halide）：在 `src/namepredict/layer1/{fg_name}.py` 建独立模块，导出 `{fg}_entries(mol) -> list[dict]` 函数（`layer1/acyl_halide.py:acyl_halide_entries:50`），由 `analyzer._detect_parts` 内的包装函数（`_acyl_chloride_entries`，`analyzer.py:283`）接入。`phosphate_entries`（`analyzer.py:98`）属同类检测器，直接内联于 `analyzer.py`。
- **共享羰基原语**：涉及 C=O 检测（酸/酯/酰胺/醛/酮）时，把通用判断函数放进 `layer1/_carbonyl_common.py` 供复用（`_has_double_bonded_o`/`_has_acid_o_neighbor`/`_amide_n_info`/`_ester_alkoxy_of` 等），`analyzer.py` 与 `layer1/acyl_halide.py` 均由该模块导入。
- **环内羰基**：环内单碳羰基与环内零碳羰基由 `_is_ketone_carbon`（`analyzer.py:129`）统一按环酮 -one 母体处理（内酰胺/环酮/内酯/硫代内酯、N-酰基环胺；环脲/环碳酸酯）。「环内杂原子」集来自 `constants.RING_HETERO`（`constants.py:36`，N/O/S）——新 FG 若让羰基经新的环内杂原子闭合成环，需同步该集合，否则该羰基不会被判定为环酮母体。

### 注册到 info dict

L1 的 FG 出口是单一 `info["fg_inventory"]`（`FunctionalGroupInventory`），由 `analyzer._collect_fgs`（`analyzer.py:477`）产出；info dict 另带 `double_bonds`/`triple_bonds`（`:439`，结构事实，独立于 FG 通道）、`carbon_ids`/`n_carbons`（`:484`）与环事实 `rings`/`n_rings`/`has_ring`/`ring_systems`/`n_ring_systems`（`_ring_meta`，`:380`）。**无 `has_*` 布尔键、无 `_fg_parts`/`_fg_bools`**：FG 存在性一律由清单内容判定（`inventory.occurrences(...)`，`functional_group_inventory.py:45`）。

以 thiol 为例，新增 FG 需要在 L1 落地：

1. 在 `analyzer.py` 定义候选原子判断 + 条目构建 + 收集函数（参照 `_thiol_entries`，`:222`）：

```python
def _is_thiol_s(atom) -> bool:      # S 原子 + 恰好一个碳邻居（R-SH 模式）
    ...

def _thiol_entry(atom) -> dict:     # 中心原子与周边原子索引
    return {"center_idx": atom.GetIdx(), "surr_idx": [_carbon_neighbor(atom).GetIdx()]}

def _thiol_entries(mol: Mol) -> list[dict]:
    return [_thiol_entry(a) for a in mol.GetAtoms() if _is_thiol_s(a)]
```

2. 在 `_detect_parts`（`analyzer.py:462`）返回的 dict 中加一条以 `FgSpec.list_key` 为键的条目列表（如 `"thiols": _thiol_entries(mol)`）——键名须与 `FG_SPECS` 的 `list_key` 逐字一致，`functional_group_inventory._LIST_CLASSES`（`:53`）按该键反查 FG 类别。
3. 在 `fg_registry.FG_SPECS`（`fg_registry.py:32`）加一条 `FgSpec`（字段与示例见步骤 2.1）；其 `p41` 同时供 `analyzer._arbitrate_parts`（`:447`）做 P-41 组合 FG 抑制。
4. 在 `functional_group_inventory.py` 的 `FunctionalGroupClass`（`:10`）加枚举成员，值须等于 `FgSpec.fg`（`_LIST_CLASSES`/`_ANCHOR_KEYS` 已由 `FG_SPECS` 派生，`:53`/`:55`，无需手写）。
5. 若条目原子集**不是**「`center_idx` ∪ `surr_idx`」的通用形态（如 phosphate 要把 P 的全部 O 邻居纳入），在 `functional_group_inventory.FG_ATOM_FNS`（`:88`）登记类别专属特征原子函数 `fn(mol, payload) -> set[int]`（现例 `_phosphate_atoms`，`:73`）。
6. 若新 FG 需要压制组合羰基 FG（酸/酯/酰卤/酰胺/醛/酸酐）或腈，把 `fg → parts 键` 加进 `constants.FG_PARTS_KEY`（`constants.py:114`）——该表是 `_arbitrate_parts` 判定「哪些 FG 在场」的唯一来源，未登记的 FG 不参与压制。

---

## 步骤 2: Layer 2 — 母体选择 (Parent Selection)

### 2.1 注册 FgSpec

主官能团优先级的单一权威是 `src/namepredict/layer1/fg_registry.py` 的 `FG_SPECS`（`fg_registry.py:32`，14 条 `FgSpec`，顺序遵循 IUPAC P-41）：`p41` 是 P-41 类号、`path` 是 P-43 优先级路径，二者合成 `PrincipalPriority`（`principal.py:17`）；`PRINCIPAL_REGISTRY`（`principal.py:48`）经 `_spec_from_fg`（`principal.py:39`）由它派生，`select_principal_group`（`principal.py:70`）据此取优先级最小的 SUFFIX 类。主官能团等级**不存于 `kind_registry`**——kind_registry 只承载 scaffold 词干（`KindMeta`，由 `ring_scaffold.ScaffoldSpec` bootstrap 注入）。

1. 在 `FG_SPECS` 加一条 `FgSpec`。字段清单见 `fg_registry.py:10`（`FgSpec` 定义）；链式 SUFFIX FG 通常只需 `fg`/`list_key`/`p41`/`anchors`/`chain=True`，按需再加 `path`、`multi`（数量后缀）、`rs`（R/S）、`keep_locant`（取代基位次保留）、`locant_kind`/`locant_source`（L4 位次记录）、`parent_anchor_fields`（单/复锚点字段）、`oh_parent`/`nh2_parent`/`oxo_parent`（醇/胺/酮母体抑制）。thiol 的声明是 `FgSpec("thiol", "thiols", p41=17, path=(2,), anchors=("surr_idx",), chain=True, multi=True, rs=True, locant_kind="sh")`（`fg_registry.py:57`）。
2. `list_key` 必须与 `analyzer._detect_parts` 返回 dict 的键一致；`anchors` 列出的 payload 键即该 FG 的母体锚点（`functional_group_inventory._ANCHOR_KEYS`，`:55`）。`expr` 默认 `"suffix"`；`principal_spec()`（`principal.py:58`）只放行 SUFFIX 表达，`"prefix_only"`/`"legacy_compat"` 取值不作主官能团（当前 `FG_SPECS` 全部为 SUFFIX）。

   磷酸即按此登记：`FgSpec("phosphate", "phosphates", p41=9, path=(1,), anchors=("p_idx",), chain=True)`（`fg_registry.py:41`）——`path=(1,)` 让同类（酯）内 `path=()` 的羧酸酯优先，磷酸随之降级为 `phosphonooxy` 前缀（P-67.1.5.1，前缀名注册于 `tools/anchored_table.py:81`）。

派生视图（新增 FG 由 `FG_SPECS` 自动进入，无需逐个登记）：`chain_fgs()`/`multi_fgs()`/`srs_fgs()`/`keep_locant_fgs()`（`fg_registry.py:63`/`:68`/`:73`/`:78`），消费方为 `principal_expression._CHAIN_FG`（`:45`）/`_MULTI_FG`（`:46`）、`locant_calc._FG_LOCANTS`（`:140`）、`assembler_prefixes._KEEP_LOCANT_KINDS`（`:42`）、`stereo._RS_KINDS`（`stereo.py:107`）。

### 2.2 母体接线（按表达权限分档）

**SUFFIX 档**：`principal_expression.py` 的 `_CHAIN_FG`（`:45`）限定可链式表达的 FG——由 `fg_registry.chain_fgs()` 派生，即 `FG_SPECS` 中 `chain=True` 的全部条目（radical/acyl/acid/phosphate/ester/acyl_halide/amide/nitrile/aldehyde/ketone/alcohol/thiol/amine）。`_chain_kind`（`:71`）按类别与多重度返回 kind：NONE（纯烃）在 count==0 时返回 `"alkane"`；ACYL / RADICAL 前支特判返回 `"acyl"`/`"radical"`（P-65.1.7.2 / P-29）；其余 `_CHAIN_FG` 成员若在 `_MULTI_FG` 内，任意 count≥1 恒返回 `group_class.value`，不在 `_MULTI_FG` 内则仅 count==1 放行。

```python
# FgSpec 设 chain=True 即进 _CHAIN_FG；再设 multi=True 即进 _MULTI_FG；_chain_kind 返回 group_class.value
```

> 设 `multi=True` 的 FG（acid/ester/amide/ketone/alcohol/thiol/amine）任意数量恒用基团名，数量由 `principal_expression_facts.multiplicity`（`principal_expression.py:119`）承载；未设的（nitrile/aldehyde/acyl_halide/phosphate）仅单基（count≠1 → None）。ALDEHYDE 开链亦仅单基，但环骨架外环 -CHO 例外（`_ring_kind`，`:176`，放行多醛，母体 kind 仍 `aldehyde`，供 L5 `_exocyclic_ring_names` 拼 -carbaldehyde/-dicarbaldehyde，P-66.6.1.1.3）。二酮由 L5 chain_engine `mult_ok` 生成式命名，L2 不产 `dione` kind。

**非 SUFFIX 表达**（`expr="prefix_only"`/`"legacy_compat"`）：`principal_spec()`（`principal.py:58`）只放行 SUFFIX 表达，这些取值不参与主官能团选择——纯醚类分子会落入纯烃兜底。`FgSpec.compat` 经 `_spec_from_fg` 投影为 `PrincipalFeatureSpec.compatibility_rank`（`principal.py:35`），当前无消费方。

### 2.3 互斥检查（无需改动）

无 `_no_fgs` 互斥谓词。互斥由 `select_principal_group`（`principal.py:70`）的结构性单选择实现：只选单个最高优先级主官能团，低优先级 FG 一律成为取代基——新增 FG 无需配置互斥 keys。

组合 FG（酸/酯/酰卤/酰胺/醛/酸酐 + 腈）另有一层 P-41 仲裁：`analyzer._arbitrate_parts`（`:447`）按 `_SUPPRESSIBLE`（`:443`）判定——存在更高优先级 FG 时组合 FG 整组退出主基团。`carboxyls`/`nitriles` 是叶型降级（`_LEAF_DEMOTED`，`:444`）：条目留在清单中并置 `demoted=True`（`functional_group_inventory.FunctionalGroupOccurrence.demoted`，`:37`），其中心碳由 `parent_skeleton._demoted_leaf_carbons`（`:53`）排除出主链；其余组合 FG 整组清空、羰基碳降级为 oxo 前缀。参与这层仲裁的 FG 必须登记在 `constants.FG_PARTS_KEY`（`constants.py:114`）。

### 2.4 添加 ownership 逻辑

ownership 由 `src/namepredict/layer2/parent_ownership.py:_kind_fg_atoms`（`:12`）通用推导：以 `principal_occurrences` 的锚点中落在母体链上的原子作种子（锚点全在骨架外时——如苯甲酸的羧基——改取与骨架相邻的锚点），再把种子的、属于该 FG 特征原子集（`occurrence.characteristic_atoms`）的邻居并入。因此**含杂原子的 FG 通常无需改动**：thiol 的 S 是锚点碳（`surr_idx`）的邻居且属特征原子集，自动进入 `owned_atoms`。

只有特征原子不在锚点邻域时才需要额外机制。phosphate 的 P 位于 O 臂之外，故在 `functional_group_inventory.FG_ATOM_FNS`（`:88`）登记 `_phosphate_atoms`（`:73`），把 P 中心与其全部 O 邻居并入特征原子，`_kind_fg_atoms` 便能沿锚点把它收进母体；O–R 臂的碳留在边界外由 L3 提取。

`compute_owned_atoms`（`:32`）取「链 ∪ `_kind_fg_atoms`」得到最终所有权集合，`finalize_parent_ownership`（`:37`）把它固化为不可变 `owned_atoms` frozenset。

若新 FG 需要向 L5 传递整名所需的计数或盐元数据，还要在 `principal_expression.py` 加 producer 分支（如 `_chain_phosphate_fields`，`:346`，含盐门控）并在 L5 `_names_for` 加 worker 分支（见步骤 5.1 方式 B）。

---

## 步骤 3: Layer 3 — 取代基提取 (Substituent Extraction)

如果新 FG 可以作为取代基前缀出现（大多数 FG 都可以），需要：

1. **保留取代名注册**：在 `src/namepredict/tools/anchored_table.py` 的 `_REGISTRY`（`:24`）中添加 `RetainedSubstituent` 条目（字段 `en`/`zh`/`anchored`/`paren`，`:16`），锚定 canonical-SMILES 写在 `anchored` 字段（如 `"sulfanyl": RetainedSubstituent("sulfanyl", "巯基", anchored=("*S",))`，`:62`）。取代基查表的入口是 `anchored_lookup`（`:157`）；`_WHOLE_ONLY_KEYS`（`:146`）中的单原子杂原子键（`*O`/`*[O]`/`*N`）只供整分子自由基命名，不参与取代基查表，避免游离 O/N 被误作侧链。import 期 `_build_anchor_index`（`:110`）会核对 anchored 键的 canonical 形式并校验跨条目唯一性，非 canonical 或重复键立即 `ValueError`。
2. **前缀命名**：由 `anchored_table` 查表或 `SubstituentNamer` 兜底完成，通常无需额外代码
3. **括号规则**：复合前缀（词干带取代）按 PIN P-16.5.1.1 加围栏，简单前缀免括（`_radical_yl_from_sub`，`as_substituent.py:66`，判据为 `meta.parent_substituent_count > 0`）；若取代名以 `oxy`/`sulfanyl` 结尾且前端 R 为简单取代基，`_obridge_front_simple`（`as_substituent.py:47`）会免去围栏（P-63.2.1/.2.2，gold/ChEBI 平铺式）——前端是否简单由命名后端 retained→recursive 判定，无法判定时保守加括号

---

## 步骤 4: Layer 4 — 编号 (Numbering) —— 通常无需改动

**`numbering_engine.orient_numbering`**（`:316`，P-14.4 候选管线）按 kind 无关的统一规则收窄：

1. 枚举候选编号（链正反 2 个 / 环每原子 1 号位 × 双向 2n 个）
2. 杂环先按 `_narrow_hetero_ring`（`:168`）元素序定起点；稠环走 `_fixed_numbering`（`:204`）/`_fused_numbering`（`:251`）
3. 按 P-14.4 逐条收窄：principal FG 附着原子最低位次集 → 多重键 → 取代基位次集 → 字母序 / CIP（`_rs_locant_key`，`:123`）破局

**新增 FG 无需注册 orienter**——收窄用的 principal 附着原子取自 `parent["principal_expression_facts"]`（`_principal_atoms`，`:139`），只要 L2 正确写入 facts 即自动生效。

FG 位次记录由 `locant_calc.py` 的 `_fg_locants`（`:142`）产出，表 `_FG_LOCANTS`（`:140`）由 `FG_SPECS` 中设了 `locant_kind` 的条目派生；取哪些原子由 `_locants_for`（`:132`）按 `FgSpec.locant_source` 决定：

- `"attachment"`（默认）：取 `principal_expression_facts` 中该 FG 类别的骨架内附着原子
- `"attachment_exocyclic"`：仅当主基以环外方式表达（`facts.relation == "exocyclic"`）时取附着原子（ester/amide/nitrile/aldehyde 用）
- `"anchor_field"`：取 `parent_anchor_fields[0]` 指定的固定 locant 1 字段（radical 用）

位次省略由 `omit_locants.omit_fg_locant`（`omit_locants.py:10`）承担，`locant_calc._omit_for`（`:108`）按 `_FG_GROUP`（`:105`，locant_kind → FG 类别）分派——新 FG 若需要环单 FG / C1–C2 的位次省略，把 `locant_kind → 类别名` 加进 `_FG_GROUP`。若新 FG 的位次须由「固定 locant 1」的扁平字段承载（如 `radical_c_idx`/`acyl_c_idx`），再把其类别加进 `principal_expression._SEMANTIC_ANCHOR_FGS`（`:55`）。

> **源:** `src/namepredict/layer4/numbering_engine.py`

---

## 步骤 5: Layer 5 — 名称组装 (Name Assembly)

### 5.1 命名实现（当前机制）

Layer5 的母体命名以 `chain_engine.py` 的 **`_KIND_TABLE`**（`:421`）链引擎为主（数据驱动）：

**方式 A：链状/环状 FG → 在 `_KIND_TABLE` 添加 `_Chain` spec。** `_Chain`（`chain_engine.py:216`）是声明式配置：`kind`/`en_suf`/`zh_suf`/`coda`、`fg`（对应 `fg_locants` 记录的 kind）、`need`、`omit_rule`、`mult_ok`、`ene_seg`/`yne_seg`、`variant` 等。thiol 的 spec 在 `chain_engine.py:481`，声明 `en_suf="thiol"`、`zh_suf="硫醇"`、`coda="ane"`、`fg="sh"`、`mult_ok=True`、`ene_seg=("ene","烯")`。`_chain_names`（`:332`）统一渲染词干、位次与省略规则——新增此类 FG 通常只需一条 spec，无需修改派发逻辑。

**方式 B：特殊拼接 → 在 `_names_for` 添加 worker 分支**（`assembler.py:372`，分支顺序即优先级）。现有 worker：`kind == "phosphate"`（`:374` → `layer5/phosphate.py:phosphate_names`，`:111`，整名自组装，不经 `_KIND_TABLE`）、`kind in EXO_RING_SUF` 的环外主基（`:378` → `_exocyclic_ring_names`，`assembler.py:78`，后缀表 `EXO_RING_SUF`〔`constants.py:179`〕驱动——加新环外主基只需在表里加一行）、`kind == "radical"` 的杂原子锚点自由基（`:382` → `_mononuclear_radical_names`，`assembler.py:292`）。

P 锚点母体自由基是「L2 词干 + L5 拼装」的接线范例：L2 `principal_expression._mononuclear_radical`（`:398`）按锚点元素查 `constants.MONONUCLEAR_BY_ELEMENT`（`constants.py:136`）取 free 名，氧化态/自由价键级词干由 `SULFUR_STEM_BY_OXO`（`:137`）/`PHOSPHORUS_STEM_BY_OXO`（`:138`，P-67.1.4.1.1.2/.6）/`NITROGEN_STEM_BY_FREE_DOUBLE`（`:139`）按 =O 数或自由价键级改写（P 上多于 1 个 =O 时明确失败），L5 由 `_mononuclear_radical_names` 调用 `_phosphoryl_sub_names`（`assembler.py:148`）拼取代基。

### 5.2 kind 收敛（无需改动）

kind 收敛在 L2 `_chain_kind`（`principal_expression.py:71`），L5 直接按 FG 类别 kind 查 `_KIND_TABLE`（无 `typed_kinds` 模块）。链式 FG 只需在 FgSpec 设 `chain=True`。若新 FG 在某个 scaffold 上有保留名（苯环的 phenol/aniline/benzoic acid 式），在 `_KIND_TABLE` 对应 entry 的 `variant` 中登记 `{scaffold_id: {multiplicity: {覆盖字段}}}`（苯环条目如 `chain_engine.py:426`）：`_names_for` 在 `assembler.py:402` 按 `scaffold_id` 取该条目，`_chain_names` 在 `:342`/`:344` 按 multiplicity 取覆盖字段。环/稠环 scaffold 的通用词干则由 `assembler._ring_stem`（`:27`）从 parent 的 `stem_en`/`stem_zh` 注入为 `stem=(...)` 覆盖（`assembler.py:393`）。

### 5.3 接线到组装流水线

`_names_for` 返回 `(en, zh)` 后，由 `assemble`（`assembler.py:562`）统一完成后续变换——`join_hydro_prefix`（`:571`，hydro + 指示氢前缀）、`join_ring_cation_suffix`（`:572`，环内 N⁺/O⁺ 的 `-ium`）、`join_kind_name`（`:574`，前缀/酯拼接）、`join_anion_names`（`:578`，阴离子）、`join_rs_prefix`（`:579`，R/S）、`join_metal_salt_names`（`:580`，金属盐）——无需为单个 FG 手动接线。若新 FG 自行组装盐形态（如 phosphate），还需在 `namer._apply_salt_suffix`（`namer.py:192`）加跳过条件（现按 `result.meta["parent_kind"] == "phosphate"` 跳过，`:196`），避免通用盐后缀重复追加。

---

## 测试

每个新 FG 必须通过以下测试验证：

1. **简单 FG 作为母体**：仅含该 FG 的简单分子。例如 `CCS`（ethanethiol / 乙硫醇）
2. **复杂 FG 作为母体**：含不同链长的分子
3. **FG 作为取代基**：该 FG 与更高优先级 FG 共存（如含 -COOH 和 -SH 的分子）
4. **与其他取代基共存**：在基本 FG 上叠加 alkyl/halogen 等取代基，验证编号和位次正确
5. **中英双语**：所有测试用例必须验证 en 和 zh 两种输出
6. **Coverage Ledger**：无 gap（无遗漏）、无 overlap（无冲突）——新增 FG 不应破坏已有 FG 的命名结果

---

## 需修改文件汇总表

| Layer | 文件 | 添加内容 |
|:---:|---|------|
| L1 | `analyzer.py:_detect_parts`（`:462`） | 条目列表键（键名 = `FgSpec.list_key`）+ 判断/条目/收集函数 |
| L1 | `layer1/{fg}.py` | 独立检测器 `{fg}_entries(mol)`（仅较复杂 FG；现例 `acyl_halide.py`） |
| L1 | `layer1/_carbonyl_common.py` | 共享羰基原语（如涉及 C=O） |
| L1 | `layer1/fg_registry.py:FG_SPECS`（`:32`） | `FgSpec` 声明（`fg`/`list_key`/`p41`/`path`/`expr`/`anchors`/`parent_anchor_fields`/`chain`/`multi`/`rs`/`keep_locant`/`locant_kind`/`locant_source`/`oh_parent`/`nh2_parent`/`oxo_parent`）；派生 `PRINCIPAL_REGISTRY`、`_CHAIN_FG`、`_MULTI_FG`、`_FG_LOCANTS`、`_RS_KINDS`、`_KEEP_LOCANT_KINDS` |
| L1 | `layer1/functional_group_inventory.py` | `FunctionalGroupClass` 枚举成员（`:10`）；特征原子非通用形态时加 `FG_ATOM_FNS`（`:88`）。`_LIST_CLASSES`/`_ANCHOR_KEYS` 自 `FG_SPECS` 派生 |
| L1 | `constants.py:FG_PARTS_KEY`（`:114`） | 参与 P-41 组合 FG 压制时登记 `fg → parts 键` |
| L2 | `principal_expression.py` | 一般无需改动（`_CHAIN_FG`/`_MULTI_FG` 自 `FG_SPECS` 派生）。仅当需新 kind、新固定 locant 1 字段或新 producer 字段时改 `_chain_kind`（`:71`）/`_ring_kind`（`:176`）/`_SEMANTIC_ANCHOR_FGS`（`:55`）/worker 分支（如 `_chain_phosphate_fields`，`:346`） |
| L2 | `parent_ownership.py:_kind_fg_atoms`（`:12`） | 一般无需改动（锚点邻域自动收敛）；特征原子不在锚点邻域时改走 `FG_ATOM_FNS` |
| L2 | `kind_registry.py` | 无需改动——只注册 scaffold 词干，FG 主基团等级不存于此 |
| L3/tools | `tools/anchored_table.py:_REGISTRY`（`:24`） | 保留取代名 `RetainedSubstituent`（含 `anchored` 锚定键与 `paren`） |
| L3 | `layer3/as_substituent.py` | 一般无需改动（括号规则由 `parent_substituent_count` / `_obridge_front_simple` 自动判定） |
| L4 | `layer4/locant_calc.py` | 一般无需改动（`_FG_LOCANTS` 自 `locant_kind` 派生）；位次省略需在 `_FG_GROUP`（`:105`）登记 `locant_kind → 类别` |
| L4 | `layer4/omit_locants.py` | 省略位次规则（如需；`numbering_engine` 无需注册 orienter） |
| L5 | `chain_engine.py:_KIND_TABLE`（`:421`） | `_Chain` spec（链式/环式 FG 的通用词干、位次与数量后缀） |
| L5 | `assembler.py:_names_for`（`:372`） | worker 分支（仅特殊拼接） |
| L5 | `layer5/{fg}.py` | 整名 worker（如需，如 `phosphate.py`） |
| L5 | `constants.py:EXO_RING_SUF`（`:179`） | 新环外主基的后缀行（如涉及） |
| Namer | `namer.py:_apply_salt_suffix`（`:192`） | 自行组装盐形态时的跳过条件（如需） |

---

## 常见陷阱

1. **`list_key` 与 `_detect_parts` 键不一致 / 漏加枚举成员**：`FgSpec.list_key` 与 `_detect_parts`（`analyzer.py:462`）的键必须逐字一致，且 `FunctionalGroupClass` 要有值为 `FgSpec.fg` 的成员，否则 `_LIST_CLASSES`（`functional_group_inventory.py:53`）构建即 `KeyError`；只在 `_detect_parts` 加键而漏登 `FG_SPECS`，清单里就没有该类别，L2 母体选择会漏掉它。
2. **`p41`/`path` 设置错误**：`p41` 是 P-41 类号（`0` = 非主官能团，`PRINCIPAL_REGISTRY` 只收 `p41 != 0` 的条目），`path` 是同类内的 P-43 路径；二者共同决定主官能团选择。只有 `expr="suffix"` 的条目才作主官能团。
3. **ownership 遗漏**：如果 FG 的杂原子不在锚点邻域（如 phosphate 的 P），必须把它纳入特征原子集（`FG_ATOM_FNS`）——遗漏会导致这些原子被 Layer 3 误识别为"未被 parent 拥有"的取代基碎片。
4. **数量派生 kind 不要造**：数量由 `facts.multiplicity` 承载（FgSpec 设 `multi=True` 后由 L5 `_generated_mult_fields` 生成数量后缀）——不要注册 `dithiol`/`trithiol` 之类的新 kind。
5. **中文命名约定不一致**：中文 stems 应遵循 IUPAC 中文命名规范（CCS 规则）。
6. **整名 worker 与盐后缀重复**：自行组装盐形态的 FG（phosphate）必须同时在 `namer._apply_salt_suffix`（`namer.py:192`）跳过通用盐后缀，否则金属盐会叠加两次。
7. **新 FG 漏进 `FG_PARTS_KEY`**：`_arbitrate_parts`（`analyzer.py:447`）只按 `constants.FG_PARTS_KEY`（`:114`）判定「哪些 FG 在场」，未登记的 FG 不参与组合 FG 压制。

---

## 相关页面

- [[architecture/layer1-analyzer]] -- Layer 1 官能团检测器架构
- [[architecture/layer2-parent-selector]] -- Layer 2 母体选择器架构
- [[architecture/layer3-substituents]] -- Layer 3 取代基提取
- [[architecture/layer4-numbering]] -- Layer 4 编号引擎
- [[architecture/layer5-name-assembly]] -- Layer 5 名称组装
- [[architecture/overview]] -- 6 层管线架构总览
- [[concepts/functional-group-priority]] -- FG 优先级体系（`p41`/`path`）详解
