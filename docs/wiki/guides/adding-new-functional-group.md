# 添加新官能团 (Adding a New Functional Group)

> **扩展指南** | 关联: [[architecture/overview]], [[architecture/layer1-analyzer]], [[architecture/layer2-parent-selector]], [[architecture/layer3-substituents]], [[architecture/layer4-numbering]], [[architecture/layer5-name-assembly]], [[concepts/functional-group-priority]], [[concepts/atom-ownership]], [[reference/core-data-contracts]]

---

## 概述

本指南说明如何向命名管线添加一个新的 functional group（官能团，FG）。工作示例是 **thiol（硫醇，R–SH）**：它是链式 SUFFIX 主官能团，走「L1 局部 SMARTS 检测 → L2 母体接线 → L5 链引擎命名」的标准路径，**L4 无需任何改动**。

新增一个 FG 要动的**表项共 5 处**，其中 ①②③④ 在 L1、⑤ 在 L5：

| # | 位置 | 作用 |
|:---:|---|---|
| ① | `layer1/fg_local_smarts.FG_SMARTS`（`fg_local_smarts.py:26`） | 加一条 `(FG 键, SMARTS)`，模式首原子即中心原子 |
| ② | `layer1/functional_group_inventory.FunctionalGroupClass`（`functional_group_inventory.py:10`） | 加枚举成员，值 = ① 的 FG 键 |
| ③ | `layer1/fg_registry.FG_SPECS`（`fg_registry.py:20`） | 登记 `FgSpec`（优先级 / 锚点 / 位次来源） |
| ④ | `layer1/analyzer._LOCAL_ENTRY_FGS`（`analyzer.py:157`） | 把 FG 键加进元组，条目才由通用路径产出 |
| ⑤ | `layer5/chain_engine._KIND_TABLE`（`chain_engine.py:457`） | 加 `_Chain` spec（该 FG 作母体时的词干与后缀） |

**①②③④ 缺一不可**，漏一处即静默失效或 import 期报错：

- 缺 ① → 无匹配中心原子；缺 ② → import 期 `ValueError`（`_ANCHOR_KEYS` 构建）；缺 ③ → 清单里没有该类别，L2 视同分子中无此 FG；缺 ④ → SMARTS 有命中但不产出条目。
- 缺 ⑤ → 该 FG 作母体时 L5 查不到 spec：环骨架回落环词干（FG 后缀静默丢失），开链母体则以 `unsupported` 失败（`assembler.py:510`）。详见「常见陷阱」5。

**`analyzer.py` 只在两类非局部判据下才需要改**：条目需要与另一 FG 的中心原子联动（酰基头 heads 联动），或需要对整分子做纯度/连通性校验（磷酸臂回接）——见 §1.4。

**L2 / L3 / L4 的表全部由 `FG_SPECS` 与 `FunctionalGroupClass` 现算**，通常无需逐条登记：

- `layer2/principal.PRINCIPAL_REGISTRY`（`principal.py:44`）——主官能团等级
- `layer2/principal_expression._chain_kind`（`principal_expression.py:70`）——kind 收敛
- `layer4/locant_calc._FG_LOCANTS`（`locant_calc.py:146`）——FG 位次记录
- `layer5/stereo._RS_KINDS`（`stereo.py:141`）——R/S 适用类别

开始之前建议先读 [[concepts/functional-group-priority]]（`p41`/`path` 体系）与 [[architecture/overview]]（6 层管线）。

---

## 步骤 1: Layer 1 — 检测 (Detection)

Layer 1 的职责是枚举各类官能团条目，汇总为带类型的 `FunctionalGroupInventory`。**检测本身是表驱动的**：SMARTS 表 + 局部命中枚举，非局部判据才回落到 `analyzer` 后置。

### 1.1 加 SMARTS 表项

`FG_SMARTS`（`fg_local_smarts.py:26`）是 `tuple[tuple[str, str], ...]`，共 17 条，每条 = `(FG 键, SMARTS)`。

```python
# 硫醇：S 上恰一个 H，与一个碳相连（R–S–H）；首原子 S 即中心原子
("thiol", "[#16;X2;H1]~[#6]"),          # fg_local_smarts.py:46
```

规则：

- **模式首原子 = 该 FG 的中心原子**；`match_local_fg`（`:68`）把它当作条目的 `center_idx`。
- **同名多条取并集**：`_compile`（`:54`）按 FG 键归并编译；`match_local_fg` 用 `by_center.setdefault(m[0], m)` 按中心原子去重，同一中心原子只保留首个命中，故同名多条不会重复计数。现例：`ketone` 三条（`:40`/`:41`/`:42`，按碳邻居数与环内杂原子分支）、`amine` 三条（`:48`/`:49`/`:50`，按取代度分支）。
- 输出 `{FG 键: [匹配元组, ...]}`，元组按中心原子索引升序；**除首元素外，匹配元组的其余元素不被通用条目路径消费**（`_fg_entry` 只读 `t[0]`，周边原子由 RDKit 邻居关系现算）。
- SMARTS 写错时**编译期立即失败**：`_compile`（`:60`）对无法解析的模式 `raise ValueError`。
- 涉及羰基时复用模块级私有原语，不要另起一份：`_ACID_O`（`:11`）、`_NOT_ACID`（`:13`）、`_NOT_ACYCLIC_ESTER`（`:15`）、`_NOT_HALO`（`:17`）、`_ONE_C`（`:19`）、`_RING_HET`（`:21`，环内杂原子，与 `constants.RING_HETERO`〔`constants.py:36`〕同义）、`_O_PHOS`（`:23`）。
- 非局部判据不写进 SMARTS：酰基头的 heads 联动、磷酸的臂回接与整分子纯度都留在 `analyzer` 后置（§1.4）。

### 1.2 登记 `FunctionalGroupClass`

`FunctionalGroupClass`（`functional_group_inventory.py:10`）是 `str, Enum`，现有 13 个实类 + `NONE = "alkane"`（`:25`，纯烃占位，不参与 FG 检测）。

加一个成员，**值必须与 `FgSpec.fg` 及 `FG_SMARTS` 的 FG 键逐字一致**：

```python
THIOL = "thiol"          # functional_group_inventory.py:23
```

必须加的原因：`_ANCHOR_KEYS`（`:54`）在 import 期、`_one`（`:106`）在构建期都做 `FunctionalGroupClass(key)`，缺成员直接 `ValueError`。

### 1.3 登记 `FgSpec`

`FgSpec`（`fg_registry.py:10`）是 `frozen dataclass`，7 个字段：

| 字段 | 默认 | 含义 |
|---|---|---|
| `fg` | 必填 | `FunctionalGroupClass` 值；同时是 L1 条目键与 L4 `fg_locants` 记录的 kind |
| `p41` | `0` | P-41 主官能团等级，**0 = 非主官能团**；越小优先级越高 |
| `path` | `()` | P-43 同类内优先级路径，与 `p41` 合成 `PrincipalPriority` |
| `expr` | `"suffix"` | 表达类型；`PrincipalExpression` 枚举只有 `suffix` 一个成员（`principal.py:23`） |
| `anchors` | `()` | occurrence payload 中承载**母体连接原子**的键；空 = 不收集锚点 |
| `parent_anchor_fields` | `None` | `(单数, 复数)` 锚点字段名，仅供固定 locant 1 的扁平字段使用 |
| `locant_source` | `"attachment"` | L4 位次原子来源：`attachment` / `attachment_exocyclic` / `anchor_field` |

thiol 的声明（`fg_registry.py:32`），与上一行的 alcohol（`:31`）对照：

```python
FgSpec("alcohol", p41=17, path=(1,), anchors=("surr_idx",)),
FgSpec("thiol",   p41=17, path=(2,), anchors=("surr_idx",)),
```

- `p41=17` 同属 P-41 第 17 类，`path` 决定醇优先于硫醇（`PrincipalPriority` 可排序，`min()` 取胜）。
- `anchors=("surr_idx",)` 因为中心原子是杂原子 S，与母体相连的是它的碳邻居。**锚点键选取规则**：中心原子即连接原子时用 `("center_idx",)`（acid/ester/amide/nitrile/aldehyde/ketone/acyl_halide），连接原子在中心周边时用 `("surr_idx",)`（alcohol/thiol/amine）。胺的多臂体现在 `surr_idx` 一次给出多个碳锚点（`functional_group_inventory.py:54` 注释，P-62.2）。
- `locant_source` 保持默认 `"attachment"`（取 `principal_expression_facts` 中该 FG 类别的骨架内附着原子）。另两个取值的适用者：`attachment_exocyclic`（仅环外表达时取，ester/amide/nitrile/aldehyde）、`anchor_field`（取 `parent_anchor_fields[0]`，radical/acyl）。
- `p41=0` 的 FG 不进 `PRINCIPAL_REGISTRY`（`principal.py:44` 过滤 `sp.p41`），也不进 `_arbitrate_parts` 的 p41 表（`analyzer.py:143`）。

### 1.4 让条目进入 `_detect_parts`（**最易漏的一步**）

`_detect_parts`（`analyzer.py:167`）的返回 dict 是清单构建的唯一输入。**FG_SMARTS 命中本身不产出条目**：通用路径由 `_local_entries`（`:161`）按 `_LOCAL_ENTRY_FGS`（`:157`）这个硬编码元组逐个 FG 键取值。

```python
_LOCAL_ENTRY_FGS = ("acid", "alcohol", "ester", "amide", "ketone", "amine", "thiol",
                    "nitrile", "acyl_halide")          # analyzer.py:157
```

**新 FG 必须把 FG 键加进该元组**，否则 `_local_entries` 不产出该键 → `build_inventory`（`functional_group_inventory.py:114`）取到空序列 → 清单里该类别零 occurrence，L2 视同分子中没有这个官能团，命名静默回落到更低优先级的母体。元组顺序不影响清单产出顺序（清单顺序随 `FG_SPECS`），仅为便于对照。

加进去后条目由 `_fg_entry`（`analyzer.py:95`）统一产出：

```python
{"center_idx": int, "surr_idx": [int, ...]}
```

`surr_idx` 由 `_surr_idx`（`:89`）算：中心原子的全部重原子邻居；中心是**非环碳**时排除环内邻居（避免环上取代基被当成 FG 周边）。

**需要非局部校验的 FG 不走这条路径**，而是在 `_detect_parts` 里显式加键做后置过滤。现有四条：

| FG 键 | 处理 | 原因 |
|---|---|---|
| `radical`（`:174`） | `_radical_entry`（`:131`），并排除 acyl 头 | 酰基残基 R–C(=O)–\* 的羰基碳不算自由基中心 |
| `aldehyde`（`:176`） | `_fg_entry`，并排除 acyl 头 | 同上，酰基头优先作 acyl |
| `acyl`（`:170`） | `_fg_entry`，**须先于 aldehyde/radical 计算 heads** | 提供 heads 集合 |
| `phosphate`（`:177`） | `phosphate_entries`（`:83`）→ `_phosphate_entry`（`:58`） | 需臂单点回接 + 整分子纯度校验（重原子 = core ∪ 臂，排除焦磷酸/臂间成环） |

`phosphate` 是**唯一不用 `center_idx`/`surr_idx` 形态的类别**：条目为 `{"p_idx", "n_oh", "n_om"}`。对应的 `FgSpec.anchors` 写 `("p_idx",)`（`fg_registry.py:24`），特征原子另由 `FG_ATOM_FNS` 定义（§1.5）。

契约：`_detect_parts` 的键集必须落在 `FG_SPECS` 的 `fg` 集合内（`tests/unit/test_architecture_contracts.py:154` 的 `test_unsaturation_stays_out_of_fg_channel`）。双键/三键不走 FG 通道，由 `_collect_fgs`（`analyzer.py:181`）写入 `double_bonds`/`triple_bonds`。

### 1.5 特征原子（默认无需改动）

`FunctionalGroupOccurrence.characteristic_atoms`（`functional_group_inventory.py:33`）默认取「`center_idx` ∪ `surr_idx`」（`center_surr_atoms`，`:78`）。thiol 无需额外配置。

**只有特征原子不在锚点邻域时**才在 `FG_ATOM_FNS`（`:87`）登记类别专属函数 `fn(mol, payload) -> set[int]`。现仅一条：`"phosphate": _phosphate_atoms`，把 P 中心与其全部 O 邻居（含 O–R 桥氧）并入特征原子。

### 1.6 P-41 组合 FG 仲裁（仅组合羰基 / 腈相关才需要）

`_arbitrate_parts`（`analyzer.py:141`）的 parts dict **键是 FG 类别值**（`acid`/`ester`/`acyl_halide`/…），三张元素表也都是 FG 键：

- `_SUPPRESSIBLE`（`:136`）：可被更高优先级 FG 整体压制的组合 FG（acid/ester/acyl_halide/amide/nitrile/aldehyde）
- `_PRESENCE_SKIP`（`:137`）：不参与存在性判定的 FG（phosphate——`p41` 与酯同为 9，纳入会改写压制结果）
- `_LEAF_DEMOTED`（`:138`）：降级为「前缀叶」的 FG（acid → carboxy、nitrile → cyano，P-61.1.3）；叶条目留在清单中并置 `demoted=True`（`FunctionalGroupOccurrence.demoted`，`functional_group_inventory.py:36`），其中心碳由 `parent_skeleton._demoted_leaf_carbons`（`parent_skeleton.py:43`）排除出主链；其余组合 FG 整组清空、羰基碳留在链内作 oxo/formyl 前缀

present 集由 `p41 = {sp.fg: sp.p41 for sp in FG_SPECS if sp.p41}` 现算（`:143`），所以新 FG 只要 `p41 != 0` 就自动参与压制判定，无需登记。只有**新 FG 本身是「应被更高优先级 FG 压制的组合 FG」**时，才把它的 FG 键加进 `_SUPPRESSIBLE`（并按需加进 `_LEAF_DEMOTED`）。

---

## 步骤 2: Layer 2 — 母体接线

### 2.1 优先级与 kind（无需改动）

`PRINCIPAL_REGISTRY`（`principal.py:44`）由 `FG_SPECS` 派生：

```python
PRINCIPAL_REGISTRY = {FG(sp.fg): _spec_from_fg(sp) for sp in FG_SPECS if sp.p41}
```

`select_principal_group`（`:60`）取 `PrincipalPriority`（`:17`，`(p41_class, p43_path)`，`order=True`）最小的类，低优先级 FG 一律成为取代基。**无互斥谓词**——互斥由「只选单个最高优先级类」这一结构性单选择实现，新增 FG 无需配置互斥 keys。

kind 由 `principal_expression._chain_kind`（`principal_expression.py:70`）返回 `FunctionalGroupClass.value`：

- `NONE`（无 FG）在 count == 0 时返回 `"alkane"`，否则 `None`
- `ACYL` 前支返回 `"acyl"`；`RADICAL` 返回 `"radical"`
- 其余类别：count ≥ 1 即返回 `group_class.value`

环骨架走 `_ring_kind`（`:169`）：scaffold 存在且主 FG 是注册类别时复用 `_chain_kind`；否则回落 `_resolved_ring_kind`（`:158`）——保留 scaffold 取其 `id`，碳环与未注册芳香稠环收敛 `"alkane"`，`benzene` 也收敛 `"alkane"`（苯环身份由 `scaffold_id` 承载）。

**因此新增 FG 不需要在 L2 写任何 kind 分派**：`FunctionalGroupClass` 成员一加，`_chain_kind` 就能为它产出 kind，L5 按该 kind 查 `_KIND_TABLE`。

### 2.2 数量与环外表达（由表驱动，无需改动代码）

- **多重度**由 occurrence 个数承载（`PrincipalExpressionFacts.multiplicity`，`principal_expression.py:36`；L5 侧 `_parent_multiplicity`，`chain_engine.py:46` 读它）。多基后缀由 `_generated_mult_fields`（`chain_engine.py:247`）自动生成——**不要为多重度另造 kind**。
- **环外主基**：若新 FG 能以环上附着方式作主基（`facts.relation == "exocyclic"`），在 `constants.EXO_RING_SUF`（`constants.py:161`）加一行 `group_class → (singular, plural)`；`chain_engine._exo_ring_spec`（`:298`）据此改写 spec 的词尾/词干/位次规则，六个环外类别（acid/aldehyde/ester/amide/nitrile/acyl）共用该管线。`plural=None` 表示该类无多取代系统名。
- **保留 scaffold 单取代保留名**：苯环走 `_BENZENE_RETAINED`（`chain_engine.py:444`，消费点 `assembler._names_for:323`）；其他 scaffold 走 `_Chain.variant` 的 `{scaffold_id: {multiplicity: 覆盖字段}}`（`chain_engine.py:203`）。

### 2.3 ownership（无需改动）

ownership 由 `parent_select._kind_fg_atoms`（`parent_select.py:63`）通用推导：

1. 种子 = 主基团 occurrence 的 `parent_anchors` ∩ 骨架；
2. 锚点全在骨架外时（苯甲酸的羧基）改取「与骨架相邻的锚点」；
3. 把种子的、属于该 FG `characteristic_atoms` 的邻居并入。

`finalize_parent_ownership`（`:83`）取「链 ∪ 上述原子」固化为不可变 `owned_atoms` frozenset，L3 据此划分母体与边界外的取代基。

**含杂原子的 FG 通常无需改动**：thiol 的 S 是锚点碳（`surr_idx`）的邻居且属特征原子集，自动进入 `owned_atoms`。只有特征原子不在锚点邻域（phosphate 的 P）才需 §1.5 的 `FG_ATOM_FNS`。详见 [[concepts/atom-ownership]]。

### 2.4 需要额外字段时的挂钩点

| 需求 | 挂钩点 |
|---|---|
| 位次由「固定 locant 1」的扁平字段承载 | `_SEMANTIC_ANCHOR_FGS`（`principal_expression.py:54`）放行类别 + `_semantic_anchor_fields`（`:57`）写字段 + `FgSpec.parent_anchor_fields` |
| L5 整名需要新的计数 / 盐元数据字段 | 在 `express_chain_principal`（`:408`）的 `kind ==` 分支加 producer，如 `_chain_phosphate_fields`（`:317`，含盐门控，门控不过返回 `None` 使候选作废） |
| 主官能团降级 / 表达标志 | `_expression_flags`（`:88`）、`_charge_state`（`:95`） |

---

## 步骤 3: Layer 3 — 取代基提取（该 FG 作前缀时）

当新 FG 不是主官能团时，它会以取代基前缀出现。需要做的事只有一件：**登记保留取代名**。

`tools/anchored_table.py` 的 `_REGISTRY`（`:24`）是保留名的注册表，条目为 `RetainedSubstituent`（`:16`，字段 `en`/`zh`/`anchored`/`paren`）。`anchored` 是 canonical SMILES 且以 `*` 哑原子标注连接点：

```python
"sulfanyl": RetainedSubstituent("sulfanyl", "巯基", anchored=("*S",)),   # anchored_table.py:61
```

- import 期 `_build_anchor_index`（`:109`）会核对 anchored 键的 canonical 形式并校验跨条目唯一性，**非 canonical 或重复键立即 `ValueError`**。
- 查表入口 `anchored_lookup`（`:156`）；`_WHOLE_ONLY_KEYS`（`:145`，`*O`/`*[O]`/`*N`）中的单原子杂原子键只供整分子自由基命名，不参与取代基查表，避免游离 O/N 被误作侧链。
- 前缀的括号规则由 `layer3/as_substituent` 与 `_obridge_front_simple` 判定，通常无需为单个 FG 配置。
- 前缀的词序、倍增与位次由 `assembler_prefixes._prefix_for`（`assembler_prefixes.py:296`）/`_build_prefix`（`:279`）统一处理。
- claim 槽位到取代基 kind 的映射在 `constants.CLAIM_KIND`（`constants.py:141`，3 键：`amine_n`/`ring_c`/`chain_c`），消费点 `layer3/substituent_extractor.py:14`；缺省槽位回落 `"side"`。新 FG 一般不涉及该表。

验证示例：`CC(O)CS` → `1-sulfanylpropan-2-ol` / 1-巯基丙-2-醇（thiol 让位于 OH，走 `sulfanyl` 前缀）。

---

## 步骤 4: Layer 4 — 编号 (Numbering) —— 通常无需改动

`numbering_engine.orient_numbering`（`numbering_engine.py:325`）的 P-14.4 候选管线按 kind 无关的统一规则收窄：枚举候选编号 → 杂环/稠环定起点 → 逐条按 principal 附着原子最低位次集 → 多重键 → 取代基位次集 → 字母序 / CIP 破局。

**新增 FG 无需注册 orienter**：收窄用的 principal 附着原子取自 `parent["principal_expression_facts"]`，只要 L2 正确写入 facts 就自动生效。

FG 位次记录由 `locant_calc._fg_locants`（`locant_calc.py:148`）产出，记录为 `{kind, locants, omit}`，`kind` 即 FG 类别值；表 `_FG_LOCANTS`（`:146`）由 `FG_SPECS` 现算：

```python
_FG_LOCANTS = tuple((sp.fg, sp) for sp in FG_SPECS if sp.fg is not None)   # locant_calc.py:146
```

取哪些原子由 `_locants_for`（`:138`）按 `FgSpec.locant_source` 决定（§1.3）。

**唯一可能需要动的地方是位次省略**：`_omit_for`（`:115`）按 `_FG_GROUP`（`:112`，`{记录 kind → principal_expression_facts 类别}`）分派到 `omit_locants.omit_fg_locant`（`omit_locants.py:5`，判据为「环状单环且无取代 → 省；否则 `pos == 1 and n_carbons <= 2`」）。

```python
_FG_GROUP = {"oh": "alcohol", "amine": "amine", "ketone": "ketone", "sh": "thiol"}   # locant_calc.py:112
```

**要为新 FG 开环单 FG / C1–C2 省略，就加一条：键 = `FgSpec.fg`（记录 kind），值 = `principal_expression_facts` 的类别名。** 未落入 `_FG_GROUP` 的 kind 恒返回 `omit=False`（L5 侧 `omit_rule` 仍兜底端位省略）。

---

## 步骤 5: Layer 5 — 名称组装 (Name Assembly)

### 5.1 在 `_KIND_TABLE` 加 `_Chain` spec（方式 A）

`_KIND_TABLE`（`chain_engine.py:457`）现有 15 entry：alcohol / ketone / alkane / acid / sulfonic / ester / phosphate / acyl / thiol / amine / aldehyde / nitrile / amide / acyl_halide / radical。

`_Chain`（`:174`）是声明式配置规格：`kind`/`en_suf`/`zh_suf`/`coda`，关键字段另有

- `fg`——对应 `fg_locants` 记录的 kind，**必须等于 `FgSpec.fg`**；
- `need`——FG 位次数要求（`need=1` 单 FG；多 FG 由数量机制生成）；
- `no_loc`——`"plain"`（无位次输出普通名）/ `"none"`（返回 None，如 ketone）；
- `omit_rule`——`(n, loc, omit) -> bool`，True 表示省略位次；
- `unsat_polyol`——多 FG 词干支持烯/炔插入；`ene_seg`/`yne_seg`/`ene_base`/`yne_suf`——不饱和段形态；
- `variant`——`{scaffold_id: {multiplicity: 覆盖字段}}`；`plain_maps`/`plain_fn`——俗名/保留名；`wrap`——整体包裹（E/Z）；`plain_hook`——词尾由分子计数直接决定（磷酸）。

thiol 的 spec（`:496`）：

```python
"thiol": _Chain(kind="thiol", en_suf="thiol", zh_suf="硫醇", coda="ane",
                fg="thiol", need=1, omit_rule=_omit_term_locant,
                ez_ene=ez_for_parent, ez_ene_multi=ez_for_parent,
                zh_full=True, unsat_polyol=True,
                ene_seg=("ene", "烯"), yne_seg=("yne", "炔")),
```

- `coda="ane"` + `en_suf="thiol"` → 词干 `ethan` + `e` 省略（`_elide_parent_e`，`:224`）+ `thiol` = `ethanethiol`；
- `need=1` + `omit_rule=_omit_term_locant`（`:8`，`omit or loc is None or (loc == 1 and n <= 2)`）→ `CCS` 出 `ethanethiol`，不出 `ethane-1-thiol`；
- `unsat_polyol=True` + `ene_seg=("ene", "烯")` → 不饱和链走段式引擎（`_chain_enyne`，`:121`）。

`_chain_names`（`:328`）统一渲染词干、位次与省略规则；**新增此类 FG 通常只需一条 spec，不动派发逻辑**。环式 FG 复用同一 entry：`_names_for`（`assembler.py:315`）在 parent 带 `stem_en`/`stem_zh` 时注入 `stem=(...)` 并把 `coda` 置空，环词干由 `kind_registry.pack_parent_stem`（`kind_registry.py:58`）在 L2 `_finalize_ranked`（`parent_select.py:118`）阶段写入。

C1/C2 英文保留名放 `constants.CHAIN_RETAINED`（`constants.py:178`）；scaffold 专属保留名放 `_Chain.variant`（苯环单取代另见 `_BENZENE_RETAINED`，`chain_engine.py:444`）。

### 5.2 特殊拼接才改 `_names_for`（方式 B）

`_names_for`（`assembler.py:301`）当前只有 3 个分支：

| 行 | 分支 |
|---|---|
| `:304` | `kind == "radical"` 且 parent 带 `radical_anchor_element` → `_mononuclear_radical_names`（`:190`，杂原子锚点母体自由基，P 的词干按 =O 数改写） |
| `:307` | `kind == "acyl_halide"` 时按 `parent.hal_z` 换成对应卤素的 spec（`_ACYL_HALIDE_BY_HAL`，`chain_engine.py:418`） |
| `:306`–`:330` | 其余一律查 `_KIND_TABLE`（`:306`），命中即交 `_chain_names` 并返回；未命中回落 `_parent_stem_names`（`:332`） |

**整名组装不在 `_names_for` 分派，而在 `join_kind_name`（`assembler.py:427`）**：`kind in ESTER_O_SIDE_KINDS`（`constants.py:143`，ester/phosphate）时走 O-侧臂拼接（`join_ester_name` / `join_phosphate_name`，`:279`）。磷酸的母体词尾仍由 `_KIND_TABLE["phosphate"]` 的 `plain_hook=_phosphate_tail`（`chain_engine.py:428`）产出，O-侧烷氧臂在拼装期并入。

### 5.3 接线到组装流水线（无需改动）

`assemble`（`assembler.py:500`）统一完成后续变换，无需为单个 FG 手动接线：

| 行 | 变换 |
|---|---|
| `:511` | `join_hydro_prefix`——hydro 与指示氢前缀 |
| `:512` | `join_ring_cation_suffix`——环内 N⁺/O⁺ 的 `-ium` |
| `:514` | `join_kind_name`——前缀 / 酯 / 磷酸拼接 |
| `:518` | `join_anion_names`——阴离子 |
| `:519` | `join_ez_prefix`——母体外挂双键 E/Z |
| `:520` | `join_rs_prefix`——R/S（`stereo._RS_KINDS`，`stereo.py:141` = 全部 `FunctionalGroupClass` 值 + `alkane`，新 FG 自动纳入） |

若新 FG 自行组装盐形态，需在 `namer._apply_salt_suffix`（`namer.py:172`）加跳过条件（现按 `result.meta["parent_kind"] == "phosphate"` 跳过，`:176`）。

---

## 测试

新增 FG 至少覆盖四类用例，双语（en/zh）都要断言：

1. **简单 FG 作母体**：`CCS` → `ethanethiol` / 乙硫醇
2. **不同链长 / 不饱和**：`CCCS`、`C=CCS` 等，验证词干、位次与段式不饱和
3. **FG 作取代基**：与更高优先级 FG 共存，如 `CC(O)CS` → `1-sulfanylpropan-2-ol`
4. **与其他取代基共存**：叠加 alkyl / halogen，验证编号与位次
5. **契约测试**：`tests/unit/test_architecture_contracts.py:154` 会校验 `_detect_parts` 键集 ⊆ `FG_SPECS` 的 `fg` 集；`tests/unit/test_coverage_ledger.py` 校验无 gap / 无 overlap——新增 FG 不应改写已有 FG 的命名结果

按主题落 `tests/unit/`（现例 `test_ether_sulfide.py` 覆盖 thiol/ether 负例）。pytest 用并行运行，benchmark 全量耗时较长，不要反复跑。

---

## 需修改文件汇总表

| Layer | 文件 | 函数 / 表（锚点） | 添加内容 |
|:---:|---|---|---|
| L1 | `layer1/fg_local_smarts.py` | `FG_SMARTS`（`:26`） | `(FG 键, SMARTS)` 表项；首原子 = 中心原子，同名多条取并集 |
| L1 | `layer1/functional_group_inventory.py` | `FunctionalGroupClass`（`:10`） | 枚举成员，值 = FG 键（缺则 import / 构建期 `ValueError`） |
| L1 | `layer1/fg_registry.py` | `FG_SPECS`（`:20`） | `FgSpec`（7 字段：`fg`/`p41`/`path`/`expr`/`anchors`/`parent_anchor_fields`/`locant_source`） |
| L1 | `layer1/analyzer.py` | `_LOCAL_ENTRY_FGS`（`:157`） | FG 键（走通用 `center_idx`/`surr_idx` 条目路径时**必须**加） |
| L1 | `layer1/analyzer.py` | `_detect_parts`（`:167`） | **仅非局部判据**：显式条目键 + 后置过滤函数 |
| L1 | `layer1/functional_group_inventory.py` | `FG_ATOM_FNS`（`:87`） | 仅当特征原子非「`center_idx` ∪ `surr_idx`」（现例 `_phosphate_atoms`，`:72`） |
| L1 | `layer1/analyzer.py` | `_SUPPRESSIBLE`（`:136`）/ `_LEAF_DEMOTED`（`:138`） | 仅当新 FG 是应被更高优先级 FG 压制的组合 FG |
| L2 | — | `PRINCIPAL_REGISTRY`（`principal.py:44`）、`_chain_kind`（`principal_expression.py:70`） | **无需改动**（由 `FG_SPECS` + `FunctionalGroupClass` 现算） |
| L2 | `layer2/parent_select.py` | `_kind_fg_atoms`（`:63`） | **无需改动**（锚点邻域自动收敛；特征原子不在邻域时改走 `FG_ATOM_FNS`） |
| L2 | `layer2/principal_expression.py` | `_chain_phosphate_fields`（`:317`）等 kind 分支 | 仅当 L5 整名需要新计数 / 盐元数据字段 |
| L2 | `layer2/principal_expression.py` | `_SEMANTIC_ANCHOR_FGS`（`:54`） | 仅当位次由固定 locant 1 的扁平字段承载 |
| L2 | `layer2/kind_registry.py` | — | **无需改动**（只注册 scaffold 词干，FG 主基团等级不存于此） |
| L3 | `tools/anchored_table.py` | `_REGISTRY`（`:24`） | `RetainedSubstituent`（保留取代名 + `anchored` 锚定键 + `paren`） |
| L3 | `layer3/as_substituent.py` | — | **无需改动**（括号规则自动判定） |
| L4 | `layer4/locant_calc.py` | `_FG_LOCANTS`（`:146`）/ `_locants_for`（`:138`） | **无需改动**（由 `FG_SPECS` 现算） |
| L4 | `layer4/locant_calc.py` | `_FG_GROUP`（`:112`） | 位次省略需登记「记录 kind → facts 类别」（键用 FG 类别值） |
| L5 | `layer5/chain_engine.py` | `_KIND_TABLE`（`:457`） | `_Chain` spec（`fg` 必须 = `FgSpec.fg`） |
| L5 | `layer5/chain_engine.py` | `_BENZENE_RETAINED`（`:444`）/ `_Chain.variant` | 保留 scaffold 单取代保留名（如需） |
| L5 | `constants.py` | `EXO_RING_SUF`（`:161`） | 环外主基系统名后缀行（如需） |
| L5 | `constants.py` | `CHAIN_RETAINED`（`:178`） | C1/C2 英文保留名（如需） |
| L5 | `layer5/assembler.py` | `_names_for`（`:301`） | worker 分支（仅特殊拼接；整名拼接在 `join_kind_name`，`:427`） |
| Namer | `namer.py` | `_apply_salt_suffix`（`:172`） | 自行组装盐形态时的跳过条件（如需） |

---

## 常见陷阱

1. **漏登记 `_LOCAL_ENTRY_FGS`**：`FG_SMARTS` + `FG_SPECS` + 枚举都加了，但没把 FG 键加进 `_LOCAL_ENTRY_FGS`（`analyzer.py:157`）→ `_detect_parts` 不产出该键 → 清单里零 occurrence → L2 视同无此 FG，命名静默回落到更低优先级母体（常表现为纯烃兜底名），且无任何报错。
2. **四处命名不一致**：`FG_SMARTS` 键、`FgSpec.fg`、`FunctionalGroupClass` 值、`_Chain.fg` 必须同字。L4 按 `sp.fg` 记 kind（`locant_calc.py:146`），L5 按 `_Chain.fg` 查记录（`chain_engine._fg_locant`，`:56`）——不一致则位次丢失，或 `no_loc="none"` 时直接返回 None 使候选作废。
3. **`anchors` 取错键**：`anchors` 决定 `parent_anchors`，进而决定 ownership 种子（`parent_select.py:63`）与 L4 附着位次（`_locants_for`）。中心原子即连接原子用 `center_idx`，连接原子在中心周边（醇/硫醇/胺）用 `surr_idx`。取错会让特征原子进不了 `owned_atoms`，被 L3 当成边界外的取代基碎片。
4. **`p41` / `path` 设错**：`p41` 越小优先级越高，`0` = 不作主官能团（`PRINCIPAL_REGISTRY` 与 `_arbitrate_parts` 的 p41 表都不收 0）。同类内用 `path` 排序（alcohol `(1,)` 压 thiol `(2,)`）。`expr` 只有 `"suffix"` 一个合法值（`PrincipalExpression`，`principal.py:23`），写别的值 import 期 `ValueError`。
5. **只注册 L5 不算接线**：`_KIND_TABLE["sulfonic"]`（`chain_engine.py:479`）与 `_BENZENE_RETAINED["sulfonic"]`（`:449`）都是惰性条目——没有 `FunctionalGroupClass` 成员、也没有 L2 生产者能产出 kind `sulfonic`，该 entry 永不被查到。**只加 L5 entry 不会带来任何命名能力**；顺序必须是 L1 检测 → L1 登记 → L2 自动收敛 → L5 spec。
6. **不要为数量造 kind**：多重度由 `facts.multiplicity` 承载（§2.2），`_generated_mult_fields`（`chain_engine.py:247`）自动生成数量后缀。注册 `dithiol`/`trithiol` 之类的新 kind 不会有任何调用方。
7. **位次省略的键名**：`_FG_GROUP`（`locant_calc.py:112`）的键是**记录 kind**（= `FgSpec.fg` = FG 类别值），值才是 `principal_expression_facts` 的类别名。用 `FgSpec` 里不存在的字段名作键，省略规则永远不会命中。
8. **特征原子非通用形态时漏改 `FG_ATOM_FNS`**：中心原子不在锚点邻域（如磷酸的 P）会被 L3 当成边界外碎片，母体名缺该原子、碎片被误命名。
9. **`_SUPPRESSIBLE` 键必须真实存在**：`_arbitrate_parts`（`analyzer.py:148`）做 `out[fg]`，登记了 `_detect_parts` 不产出的键会直接 `KeyError`。
10. **整名 worker 与盐后缀重复**：自行组装盐形态的 FG 必须同时在 `namer._apply_salt_suffix`（`namer.py:172`）跳过通用盐后缀，否则金属盐叠加两次。
11. **中文命名约定**：中文词干须遵循中国化学会命名原则（CCS）；双语输出一律同时验证。

---

## 相关页面

- [[architecture/layer1-analyzer]] -- Layer 1 官能团检测（`FG_SMARTS` 表 + 清单构建）
- [[architecture/layer2-parent-selector]] -- Layer 2 母体选择与主基团接线
- [[architecture/layer3-substituents]] -- Layer 3 取代基提取与锚定表
- [[architecture/layer4-numbering]] -- Layer 4 编号引擎与位次记录
- [[architecture/layer5-name-assembly]] -- Layer 5 链引擎与名称组装
- [[architecture/overview]] -- 6 层管线架构总览
- [[concepts/functional-group-priority]] -- FG 优先级体系（`p41` / `path`）
- [[concepts/atom-ownership]] -- 母体所有权（`owned_atoms`）
- [[reference/core-data-contracts]] -- info dict / parent dict 字段契约
- [[guides/adding-new-ring-system]] -- 添加新环系（骨架扩展指南）
