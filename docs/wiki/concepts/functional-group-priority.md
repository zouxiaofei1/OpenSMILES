# 官能团优先层级 (Functional Group Priority Hierarchy)

> **核心概念** | 关联: [[architecture/overview]], [[architecture/layer1-analyzer]], [[architecture/layer2-parent-selector]], [[architecture/layer5-name-assembly]], [[concepts/atom-ownership]]

---

## 概述

IUPAC 有机命名中，一个分子可能同时含有多种官能团（Functional Group, FG），例如羟基酸（含 -COOH 和 -OH）、氨基酮（含 -NH&#8322; 和 >C=O）、氰基酯（含 -CN 和 -COOR）等。按照 IUPAC P-41 规则，这些官能团之间存在严格的**优先顺序**：优先级最高的官能团成为 **principal characteristic group（主特征基团）**，以母体后缀（suffix）表达；较低优先级的官能团退化为取代基前缀（prefix）。

NamePredict 将这一优先级体系的**单一事实来源放在 `fg_registry.FG_SPECS`**（`FG_SPECS` 于 `fg_registry.py:32`）的 `p41` / `path` 两个字段上，由 `principal.py` 的 `PRINCIPAL_REGISTRY`（`principal.py:48`）派生为 `PrincipalPriority(p41_class, p43_path)`（`principal.py:16`，`order=True` 可比较），`select_principal_group`（`principal.py:70`）按其**升序取最小**选出主基团，并在三层流水线中接力使用：

1. **Layer 1（`analyzer.py`）**: 检测 FG 时采用排他性优先级（carboxyl 排斥 ester/amide/anhydride），并以 `_arbitrate_parts`（`analyzer.py:447`）按 p41 做 P-41 主基团仲裁
2. **Layer 2（`parent_skeleton.py` + `parent_selector.py`）**: 先由 `select_principal_skeletons`（`parent_skeleton.py:253`）做骨架级 P-44 筛选，再由 `_rank_candidates`（`parent_selector.py:32`）按 P-44 键排序、`_reorder_p45_2`（`parent_selector.py:44`）按 P-45.2.1 重排并取并列组
3. **Layer 5（`chain_engine.py`）**: 直接查 `_KIND_TABLE`（`chain_engine.py:421`）分派后缀（kind 收敛在 L2 `_chain_kind`，无 `typed_kinds` 模块），高优先级 FG 获得后缀，低优先级 FG 转为前缀

> 磷酸/磷酸酯（`phosphate`）由 `layer1/analyzer.py` 的 `phosphate_entries`（`analyzer.py:98`）检测，并在 `FG_SPECS` 登记为 `FgSpec("phosphate", "phosphates", p41=9, path=(1,), anchors=("p_idx",), chain=True)`（`fg_registry.py:41`），已进入优先级体系；其余 12 个扩展 FG（sulfoxide/sulfone/sulfonate/sulfonamide/sulfonic_acid/sulfonyl_chloride/boronic/carbamate/carbonate/urea/guanidine/hydrazine）在 layer1 中没有检测器，不进入优先级体系。

---

## 官能团优先级表（`FgSpec.p41` / `FgSpec.path`）

以下表格展示 NamePredict 中实现的官能团优先级，遵循 IUPAC P-41 顺序。**`p41` 数值越小优先级越高**（`select_principal_group` 取 `min`），`path` 为同 `p41` 类内的 P-43 决胜键。表格依据 `src/namepredict/layer1/fg_registry.py` 的 `FG_SPECS`（`fg_registry.py:32`）逐条生成，经 `principal.py` 的 `PRINCIPAL_REGISTRY`（`principal.py:48`）投影为 `PrincipalFeatureSpec`。

| `p41` | `path` | 官能团种类 (kind) | 英文后缀示例 | 中文后缀示例 | 说明 |
|:---:|:---:|---|---|---|---|
| **1** | `()` | `radical` | -yl | -基 | 自由价（P-29/P-31）；自由基压制所有组合羰基 FG |
| **1** | `()` | `acyl` | -oyl | -酰基 | 酸衍生的酰基残基头（P-65.1.7.2） |
| **7** | `(1,)` | `acid` | -oic acid | -酸 | 羧酸 |
| **8** | `()` | `anhydride` | -oic anhydride | -酸酐 | 登记保留；L1 不产 occurrence，不参与候选生成 |
| **9** | `()` | `ester` | -oate | -酸酯 | 与 phosphate 同 `p41=9`，`path=()` 使其优先（P-43） |
| **9** | `(1,)` | `phosphate` | -phosphate / phosphoric acid | -磷酸 / -磷酸酯 | 整名由 L5 `phosphate_names`（`layer5/phosphate.py:111`）组装（P 中心无碳词干） |
| **10** | `()` | `acyl_halide` | -oyl fluoride/chloride/bromide/iodide | -酰氟/酰氯/酰溴/酰碘 | 后缀随实际卤素（`parent.hal_z`）变（P-65.5） |
| **11** | `()` | `amide` | -amide | -酰胺 | 酰胺 |
| **14** | `()` | `nitrile` | -nitrile | -腈 | C&#8801;N |
| **15** | `()` | `aldehyde` | -al / carbaldehyde | -醛 / 甲醛 | 醛基 |
| **16** | `()` | `ketone` | -one | -酮 | 酮，含环内单碳羰基（内酰胺/内酯/硫代内酯、N-酰基环胺，`_is_ketone_carbon` 按「单碳/零碳邻居 + 环内杂原子 N/O/S」判定）（`dione` 由 L5 chain_engine `mult_ok` 生成式命名产生） |
| **17** | `(1,)` | `alcohol` | -ol | -醇 | 羟基；与 thiol 同 `p41=17`，`path=(1,)` 使其优先（P-43） |
| **17** | `(2,)` | `thiol` | -thiol | -硫醇 | 巯基 |
| **19** | `()` | `amine` | -amine | -胺 | 氨基 |

> **源:** `src/namepredict/layer1/fg_registry.py` 的 `FG_SPECS`（`p41`/`path`）→ `src/namepredict/layer2/principal.py` 的 `PRINCIPAL_REGISTRY`（P-41 链状 FG）

> 注 1：`FgSpec` 另有 `compat` 字段（`fg_registry.py:17`），14 条登记项均未赋值因而恒为 `0`，`PrincipalFeatureSpec.compatibility_rank`（`principal.py:35`）随之恒为 `0`，不承载任何排序。
>
> 注 2：`p41=1` 的 `radical` 与 `acyl` 两行键完全相同，`select_principal_group` 的 `min` 在并列时取集合迭代顺序首位——实际互斥由 L1 `_arbitrate_parts`（`analyzer.py:447`）与 L2 所有权判定各自收敛。
>
> 注 3：数量派生 kind（`diacid`/`polycarboxylic`/`diol`/`triol`/`diamine`/`triamine`/`tetraamine`）与组合 kind（`cycloalcohol`/`cycloketone`/`cycloamine`/`cycloalkane_polycarboxylic` 等）不存在——链式 acid/alcohol/amine/ketone 对任意主基团数 kind 恒为基团名，数量由 `principal_expression_facts.multiplicity` 承载。苯/杂环保留名（benzoic/phenol/aniline 等）由 L5 chain_engine `variant` 提供，不占独立 kind 行。

---

## KindMeta 数据结构

`kind_registry.py:7-16` 定义了 `KindMeta` 数据类，作为 scaffold 母体种类的统一元数据容器：

```python
@dataclass(frozen=True)
class KindMeta:
    kind: str              # 母体种类标识符（如 "benzene", "naphthalene"）
    en: str | None         # 英文 stem 名称
    zh: str | None         # 中文 stem 名称
    ring: str              # 环系类型: "none" | "hetero" | "carbo"
    n_rings: int           # 环的数量
    retained: bool         # 是否为 IUPAC 保留名（retained name）
```

名称为 `None` 的 scaffold 不登记（`kind_registry.py:88`）。所有 scaffold 母体种类通过模块级 `_bootstrap()` 函数统一注册（`kind_registry.py:94`），注册**只有一步**：`_load_from_scaffold_specs()`（`kind_registry.py:83`，从 `ring_scaffold.all_specs()` 读取）。**`ring_scaffold.py` 的 `_TEMPLATES`** 是 **stem 的最终权威来源**。链式 FG kind 不预先注册——`KindMeta` 不承载主官能团等级，等级由 `FgSpec.p41` → `PrincipalPriority` 承载。

---

## 三层 FG 生命周期

### Layer 1: 检测（Detect）—— 排他性优先级

`src/namepredict/layer1/analyzer.py` 在检测 FG 时，必须处理**化学上的包含关系**：一个 carbon 如果已经是 carboxyl carbon，就不能同时被识别为 ester 或 amide。为实现这一点，Layer 1 采用"排他性检测"模式：

- **Carboxyl 优先于 ester/amide/anhydride**: `_is_carboxyl_carbon()`（`analyzer.py:107`）匹配 C(=O)OH 或 C(=O)O⁻ 模式（`_has_double_bonded_o` + `_has_acid_o_neighbor`），`_is_ester_carbon`（`analyzer.py:176`）与 `_is_amide_carbon`（`analyzer.py:117`）显式判定酸性氧邻居为 False 后才匹配。
- **酸酐桥氧从酯氧中排除**: `_is_ester_alkoxy_o()`（`analyzer.py:155`）先检查 `_is_anhydride_bridge_o()`（`_carbonyl_common.py:56`）。
- **环内酯并入环母体**: `_is_ester_carbon` 经 `_is_lactone_carbon`（`analyzer.py:167`）排除酯氧在环内者，交 `_is_ketone_carbon` 作环酮（2H-chromen-2-one / 2-benzofuran-1-one / 1,3-dioxolan-2-one）。
- **环内 N 不作酰胺**: `_amide_n_info`（`_carbonyl_common.py:104`）排除环内 N（P-66.1.1），故 N-酰基环胺（1-(pyrrolidin-1-yl)ethanone）的羰基改由 `_is_ketone_carbon`（`analyzer.py:129`）按"单碳邻居 + 环内杂原子"判为酮母体——环内杂原子由 `_has_ring_hetero_neighbor`（`analyzer.py:125`）取 `constants.RING_HETERO`（N/O/S），故环内 O/S 与羰基同环的内酯/硫代内酯与环内 N 的 N-酰基环胺同归环酮（`oxolan-2-one`/`thiolan-2-one`/`1-pyrrolidin-1-ylethanone`）；零碳邻居的环内羰基（环脲/环碳酸酯型）同样作环酮。
- **醛/酮边界**: `_is_aldehyde_carbon`（`analyzer.py:190`）要求羰基带 H（无 H 的 N-酰基/环酮/内酰胺不作醛），环内羰基遂由 `_is_ketone_carbon` 作环酮、以 `-one` 后缀表达（P-66.6.1）——两函数互斥，同一羰基不会被酮/醛双计。
- **锚定酰基头优先于自由基**: `_is_acyl_head`（`analyzer.py:399`）识别带 `*` 锚点、带 =O、恰 1 个单键碳邻居且无其它重邻居的羰基碳（P-65.1.7.2），其条目进 `acyls`；`_radical_entries`（`analyzer.py:427`）以该碳集为排除集，避免醛→酮误降级。

carbamate/carbonate/urea/guanidine/sulfoxide/sulfone 等 12 个扩展 FG 在 layer1 中没有检测器，因此无跨模块排他链（carbamate 排除 ester、urea 排除 amide 等）——排他检测只发生在 analyzer.py 内部的核心羰基族与酰卤模块。共享的羰基检测原语集中在 `_carbonyl_common.py`，酰卤检测独立在 `layer1/acyl_halide.py`（`acyl_halide_entries`，`acyl_halide.py:50`）。磷酸走 `analyzer.py` 的 `phosphate_entries`（`analyzer.py:98`）/ `_one_phosphate`（`analyzer.py:44`）：识别 P(=O)(O)₃ 中心（恰好 1 个 =O、3 个单键 O，整分子重原子须全部落在中心与臂内），产出 `{p_idx, n_oh, n_om, n_arms}`，不与羰基族互斥。

所有检测到的 FG 经 `_detect_parts`（`analyzer.py:462`）汇总、`_arbitrate_parts`（`analyzer.py:447`）按 P-41 仲裁后，由 `build_inventory`（`functional_group_inventory.py:115`）收敛为**一个** `FunctionalGroupInventory`（14 个 `FunctionalGroupClass`），随 `double_bonds`/`triple_bonds` 一起构成 info dict 传递给 Layer 2。P-41 仲裁把被更高优先级 FG 压制的组合羰基 FG 整组清空（其羰基碳降级入 `ketones` 作 oxo 前缀候选），叶型降级（`_LEAF_DEMOTED = ("carboxyls", "nitriles")`，`analyzer.py:444`）则在 occurrence 上打 `demoted` 标记、其碳排除出主链（P-61.1.3 carboxy/cyano）。

> **源:** `src/namepredict/layer1/analyzer.py`, `src/namepredict/layer1/_carbonyl_common.py` | 详情见 [[architecture/layer1-analyzer]]

### Layer 2: 评分与选择（Score & Select）

L2 选择分**骨架级**与**候选级**两段，评分键分散在 `parent_skeleton.py` 与 `parent_selector.py`：

**骨架级（P-44，`parent_skeleton.py`）**: 主官能团先由 `select_principal_group`（`principal.py:70`）选出（入口 `select_principal_parent_skeletons`，`principal_parent.py:21`），再由 `select_principal_skeletons(info, occurrences)`（`parent_skeleton.py:253`）在候选骨架上依次施加——`enumerate_principal_skeletons`（`parent_skeleton.py:266`）枚举开链/环骨架 → `keep_max_principal_coverage`（`parent_skeleton.py:198`）取覆盖主基团最多者 → `keep_p44_1_2`（`parent_skeleton.py:147`，混合拓扑时环优先 + `_SENIOR_ATOMS` 最优先元素）→ 开链走 `keep_p44_3`（`parent_skeleton.py:169`，键 `p44_3_key`：杂原子数、原子数、元素计数）/ 环走 `keep_p44_2`（`parent_skeleton.py:191`，键 `p44_2_key`：含杂原子、N 计数、senior、环数）→ `keep_p44_4_unsaturation`（`parent_skeleton.py:240`，键：多重键数、双键数）。主官能团类别由 `_chain_kind` / `_ring_kind`（`principal_expression.py:71` / `:176`）投影为母体 kind。

**候选级（`parent_selector.py`）**: `_rank_candidates`（`parent_selector.py:32`）按 P-44 键 `_p44_1_1`（`parent_selector.py:26`）降序排列。该键由 `principal_key`（`parent_selector.py:21`）从 `P44Facts`（`parent_selector.py:7-11`）取出：

```python
@dataclass(frozen=True, order=True)
class P44Facts:
    principal_group_class: int    # 主官能团类等级
    principal_group_count: int    # 主官能团实例数
```

`principal_key` 以 `P44Facts(None, parent["principal_group_count"])` 构造，类等级维当前恒为 `None`，故排序的实际判别量是 `principal_group_count`。随后 `_reorder_p45_2`（`parent_selector.py:44`）按 **P-45.2.1** 前缀取代基团数（`_p45_2_prefix_count`，`parent_selector.py:37`，`iter_claims` 枚举所得）降序稳定重排，`tied=True` 时只返回并列最大组。`select_parent`（`parent_selector.py:69`）即返回该并列组（`list[dict]`）。

P-44 评分并列时不直接取首位：`namer` 对并列组每个候选各跑一次 L3–L5（上限 `_MAX_TIED_CANDIDATES=4`，`namer.py:107`），再由 `_best_hit`（`namer.py:154`）裁决——P-44.1.1 后缀位次集合（`suffix_locant_set`，`candidate_keys.py:7`）已分胜负时保持候选顺序，仍并列时才取 P-45.2.2 前缀位次集合（`prefix_locant_set`，`candidate_keys.py:22`）最小者。

> **源:** `src/namepredict/layer2/parent_skeleton.py:253`, `src/namepredict/layer2/parent_selector.py:26-74`, `src/namepredict/layer2/principal.py:70` | 详情见 [[architecture/layer2-parent-selector]]

### Layer 5: 后缀分派（Suffix Dispatch）

Layer5 由 `_names_for`（`assembler.py:372`）查 **`chain_engine._KIND_TABLE`**（`chain_engine.py:421`，14 个 `_Chain` spec：`alcohol`/`ketone`/`alkane`/`acid`/`sulfonic`/`ester`/`acyl`/`thiol`/`amine`/`aldehyde`/`nitrile`/`amide`/`acyl_halide`/`radical`）渲染词干、不饱和段、位次与环前缀（kind 收敛在 L2 `_chain_kind`，无 `typed_kinds` 模块；其中 `sulfonic` 一行只在本表登记，L2 无产生它的 FG 类别）。`acyl_halide` 一行的实际 spec 由 `_ACYL_HALIDE_BY_HAL`（`chain_engine.py:419`，按 F/Cl/Br/I 键控，`_ac_hal_chain` 于 `chain_engine.py:402` 生成）按 `parent.hal_z` 覆盖。`mult_ok` 生成式由 `_generated_mult_fields`（`chain_engine.py:294`）按 `principal_expression_facts.multiplicity`（`_parent_multiplicity`，`chain_engine.py:52`）派生数量后缀（alcohol→diol/triol/tetraol，amine→diamine/triamine/tetraamine，acid→dioic acid），`variant`（`_Chain` 字段，`chain_engine.py:248`；在 `_chain_names`，`chain_engine.py:332` 应用）仅作 scaffold 特例覆盖（苯 → phenol/benzoic acid/aniline 等）。环外（exocyclic）FG 走 worker——`_exocyclic_ring_names`（`assembler.py:78`）按后缀表 `EXO_RING_SUF`（`constants.py:179`）拼 …-carboxylic acid/…-carboxylate/…-carboxamide/…-carbonitrile/…-carbaldehyde/…-carbonyl（多羧酸 → -dicarboxylic acid，环外二醛 → -dicarbaldehyde，P-66.6.1.1.3；酰基头 -carbonyl，furan-2-carbonyl，P-65.1.7.2；苯单取代 → 保留名回落 chain_engine variant）——否则落 `_parent_stem_names`（`assembler.py:410`）回退。单环环烷/环烯的环外系统名共用 `_ring_carbocycle_stem`（`assembler.py:47`，环烯/另带前缀取代时后缀 locant 显式）。

`phosphate` 是例外分支：`_names_for` 首条即 `kind == "phosphate"`（`assembler.py:374`），不经 `_KIND_TABLE`，直接调 `layer5/phosphate.py` 的 `phosphate_names`（`phosphate.py:111`）组装整名——P 中心无碳词干，用 L2 `_chain_phosphate_fields`（`principal_expression.py:346`）注入的 `n_oh`/`n_om`/`n_arms`/`salt_meta` 与 `o_side` 臂（`constants.ESTER_O_SIDE_KINDS` 含 `"phosphate"`，`constants.py:156`；判定于 `claim_extract.py:81`）拼装；`namer._apply_salt_suffix`（`namer.py:192`）对 `parent_kind == "phosphate"` 跳过通用金属盐后缀（盐形态已在整名内处理）。杂原子锚点自由基另走 `_mononuclear_radical_names`（`assembler.py:292`），按其氧化态/自由价键级选定的单核氢化物词干取名，不经 `_KIND_TABLE`。

调度是分层的：L1 收敛 occurrence → L2 收敛 kind → L5 查链引擎 → 命中后直接返回。这确保了被选为母体的 principal FG 获得后缀，而劣后 FG 在 Layer 3 中被转为取代基前缀（如 hydroxy-、oxo-、amino-）。

> **源:** `src/namepredict/layer5/chain_engine.py:421`, `src/namepredict/layer5/assembler.py:372` | 详情见 [[architecture/layer5-name-assembly]]

---

## 互斥排除（P-41 互斥约束）

IUPAC P-41 规定某些 FG 之间不能作为母体共存。NamePredict 通过 **`principal.py` 的 `select_principal_group`**（`principal.py:70`）结构性实现互斥：在未降级（`not entry.demoted`）的 occurrence 里筛出 `principal_spec`（`principal.py:58`，只认 `PrincipalExpression.SUFFIX`）非空的类，按 `PrincipalPriority` 取 `min` 选单个最高优先级主官能团，低优先级 FG 一律成为取代基，不需要逐候选互斥检查。L1 侧另有一道 `_arbitrate_parts`（`analyzer.py:447`）P-41 仲裁：`_SUPPRESSIBLE`（`analyzer.py:443`，组合羰基 FG + 腈）在存在更高优先级 FG 时整组退出主基团。

> 无 `fg_helpers.py` 与 `_no_fgs(info, keys)` 互斥谓词——互斥由 `select_principal_group` 结构性单选择实现。

---

## 注册体系：kind_registry 作为词干注册中心

`src/namepredict/layer2/kind_registry.py` 是**词干与环元数据的注册查询中心**：保留 scaffold 的 stem/ring/n_rings/retained 由 `ring_scaffold.py` 的 `_TEMPLATES` 派生（bootstrap 唯一一步 `_load_from_scaffold_specs`，`kind_registry.py:83`），**不承载主官能团等级**（等级单一权威在 `fg_registry.FG_SPECS.p41/path` → `principal.PRINCIPAL_REGISTRY`）。公共 API：

| API | 功能 |
|---|---|
| `kind_registry.get(kind)`（`kind_registry.py:20`） | 查 `KindMeta`（未注册返回 None） |
| `kind_registry.parent_names(kind)`（`kind_registry.py:25`） | 返回 `(en_stem, zh_stem)` 元组 |
| `kind_registry.pack_parent_stem(parent, mol)`（`kind_registry.py:65`） | 为候选补齐 `stem_en`/`stem_zh`（含五元杂环 locant 前缀终态化）与 `numbering_scaffold`/`numbering_scaffold_required`；由 `parent_selector._finalize_ranked`（`parent_selector.py:57`）调用 |

新增 FG 种类时，若走 principal typed 管线（acid/alcohol/ketone/amine 等），在 `fg_registry.FG_SPECS` 登记 `p41`/`path`/`chain` 等字段，并在 `principal_expression.py` 的 `_CHAIN_FG`（`principal_expression.py:45`，由 `fg_registry.chain_fgs()` 派生）覆盖的集合内即可。评分、候选收集与组装逻辑不变，仍保持**开放-封闭原则**。

> **源:** `src/namepredict/layer2/kind_registry.py`

---

## 总结

官能团优先级是 NamePredict 命名正确性的基石。从 Layer 1 的排他性检测、Layer 2 的主基团等级与骨架评分、Layer 5 的后缀分派，整个系统严格遵循 IUPAC P-41 和 P-44 的层级化思维：

1. **检测层保证"一碳一 FG"**——排他性模式使得同一碳原子不会被重复识别为多个 FG
2. **评分层保证"优先级强制"**——`select_principal_group` 按 `PrincipalPriority` 取最小选出唯一的 principal FG，骨架级 P-44 规则在其上收敛候选
3. **组装层保证"劣后 FG 降级"**——低优先级 FG 从候选后缀降为取代基前缀
4. **互斥层保证"互斥正确"**——`select_principal_group` 的单选择结构性排除无法共存的 FG 组合
