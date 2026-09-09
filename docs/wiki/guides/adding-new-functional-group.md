# 添加新官能团 (Adding a New Functional Group)

> **扩展指南** | 关联: [[architecture/overview]], [[architecture/layer1-analyzer]], [[architecture/layer2-parent-selector]], [[architecture/layer3-substituents]], [[architecture/layer4-numbering]], [[architecture/layer5-name-assembly]], [[concepts/functional-group-priority]]

---

## 概述

本指南详细说明如何向 NamePredict 命名管线添加一个新的 functional group（官能团，FG）类型。整个流程跨越 L1（检测）、L2（母体选择）、L5（名称组装）三层；**L4 通常无需改动**（`numbering_engine` 的 P-14.4 候选管线自动适用）。

我们以 **thiol（硫醇，R-SH）** 作为工作示例——它是已实现的中优先级（`fg_rank=4`）SUFFIX 类链状 FG，代表"检测于 analyzer.py 内联 + SUFFIX 档母体 + `_KIND_TABLE` 命名"的标准路径。

> 新增 FG 的核心是三步：**L1 检测 → L2 `PRINCIPAL_REGISTRY`/`_CHAIN_FG` → L5 `_KIND_TABLE`**。无 `_no_fgs` 互斥检查、L4 orienter 注册、组合 kind。

在开始之前，建议先阅读 [[concepts/functional-group-priority]] 了解 FG 优先级体系，以及 [[architecture/overview]] 了解 6 层管线架构。

---

## 步骤 1: Layer 1 — 检测 (Detection)

Layer 1 的职责是从 RDKit `Mol` 对象中检测官能团，输出结构化的 atom index 条目和布尔标志。

### 检测位置

- **简单 FG**（如 thiol, ether, sulfide）：直接在 `analyzer.py` 中内联实现，通常只需要 3 个函数：候选原子判断、条目构建、收集函数。
- **较复杂 FG**（如 isocyanate, acyl halide, phosphate）：可创建独立模块 `src/namepredict/layer1/{fg_name}.py`，导出 `{fg}_entries(mol) -> list[dict]` 函数。当前有 `isocyanate.py`、`acyl_halide.py`、`phosphate.py` 三个独立模块。
- **共享羰基原语**：如果新 FG 涉及 C=O 检测（酸/酯/酰胺/醛/酮），把通用判断函数放进 `_carbonyl_common.py` 供复用。

### 注册到 info dict

以 thiol 为例（`analyzer.py` 内联）：

1. 定义候选原子判断函数 + 条目构建函数 + 收集函数：

```python
def _is_thiol_s(atom) -> bool:
    # S 原子 + 恰好一个 H（R-SH 模式）
    ...

def _thiol_entry(atom) -> dict:
    return {"s_idx": atom.GetIdx(), "c_idx": ...}

def _thiol_entries(mol: Mol) -> list[dict]:
    return [_thiol_entry(a) for a in mol.GetAtoms() if _is_thiol_s(a)]
```

2. 在 `_fg_parts`（`analyzer.py:536`）中添加条目列表键（如 `"thiols": _thiol_entries(mol)`）；若该键不在 `_fg_lists` 的 5 个核心键内，还要同步加进 `_fg_more_lists`（`analyzer.py:440`）——`phosphates` 即如此加入
3. 在 `_fg_bools` 的 `_FG_BOOL_MORE_KEYS`（`analyzer.py:20`）中添加布尔标志（如 `"has_phosphate"`）
4. 如果新 FG 需要类型化的 `FunctionalGroupClass`，在 `functional_group_inventory.py` 中添加：枚举成员 + `_LIST_CLASSES` 键映射 + `_ANCHOR_KEYS` 锚点键（`_LIST_CLASSES`/`_ANCHOR_KEYS` 均自 `FG_SPECS` 派生，`functional_group_inventory.py:68`/`:71`；thiol 已存在：`THIOL` + `"thiols"` → `s_idx`）

---

## 步骤 2: Layer 2 — 母体选择 (Parent Selection)

### 2.1 注册 fg_rank

`fg_rank` 的单一权威是 `src/namepredict/layer1/fg_registry.py` 的 `FG_SPECS`（`FgSpec.compat`，遵循 IUPAC P-41 顺序；`FG_SPECS` 于 `fg_registry.py:36`），`PRINCIPAL_REGISTRY`（`principal.py:51`）由它派生。主官能团等级**不存于 `kind_registry`**——kind_registry 仅注册 scaffold 词干（`KindMeta`）；旧式 parent 的等级兜底由 `parent_candidate._kind_rank`（`parent_candidate.py:20`）按 kind → `legacy_rank`（`principal.py:67`）实时投影。新增 chain FG 只需：

1. 在 `fg_registry.FG_SPECS` 加一条 `FgSpec` 声明 `p41`/`compat`/`expr`/`anchors`/`parent_anchor_fields`/`chain`/`multi`/`rs`/`locant_kind` 等字段。SUFFIX 档 `expr="suffix"`（默认，如 `FgSpec("thiol", "thiols", p41=17, path=(2,), compat=4, ...)`，`fg_registry.py:76`）；非 SUFFIX 类（如 sulfide/ether）用 `expr="legacy_compat"/"prefix_only"`。需在 L4 `fg_locants` 记主基团位次的设 `locant_kind`（aldehyde/acid/amide 等即如此），支持 R/S 的设 `rs=True`。`PRINCIPAL_REGISTRY` 经 `_spec_from_fg`（`principal.py:39`）自动派生，无需手写 `_suffix`。

   磷酸即按此登记：`FgSpec("phosphate", "phosphates", p41=9, path=(1,), compat=10, anchors=("p_idx",), parent_anchor_fields=("p_idx","p_idxs"), chain=True)`（`fg_registry.py:51`）——`path=(1,)` 让同类（酯）内 `path=()` 的羧酸酯优先，磷酸随之降级为 `phosphonooxy` 前缀（P-67.1.5.1，前缀名注册于 `tools/anchored_table.py:77`）。

### 2.2 母体接线（按表达权限分档）

**SUFFIX 档**：`principal_expression.py` 的 `_CHAIN_FG`（`:46`）限定可链式表达的 FG——由 `fg_registry.chain_fgs()` 派生，含 acid/ester/acyl_halide/amide/nitrile/aldehyde/ketone/alcohol/thiol/amine/**acyl**/**phosphate**（锚定酰基残基/环外酰基头由 `_chain_kind` 特判返回 "acyl"，P-65.1.7.2；radical 亦在 `_chain_kind` 前支特判，不在 `chain_fgs()`），`_chain_kind`（`:62`）按类别与多重度返回 kind：

```python
# _CHAIN_FG 含该 FG 类即自动支持；_chain_kind 返回 group_class.value
```

> ACID/ALCOHOL/AMINE/KETONE 任意 count≥1 恒返回基团名（`_MULTI_FG`，数量由 `principal_expression_facts.multiplicity` 承载）；ESTER/AMIDE/NITRILE/ACYL_HALIDE 仅单基（count≠1 → None）；ALDEHYDE 开链亦仅单基，但环骨架外环 -CHO 例外（`_ring_kind` 放行多醛，母体 kind 仍 `aldehyde`，供 L5 `_exocyclic_aldehyde_names` 拼 -carbaldehyde/-dicarbaldehyde，P-66.6.1.1.3）。`dione` 不再由 L2 产生，二酮由 L5 chain_engine `mult_ok` 生成式命名。新增 SUFFIX FG 若属同类，把其类别加进 `_CHAIN_FG` 即可。

**LEGACY_COMPAT 档**（sulfide/isocyanate/isothiocyanate 等非 SUFFIX 类）：仅取 `compatibility_rank`（按 FG 枚举经 `legacy_rank` 查询，`principal.py:64`），但**不参与主官能团选择**（`principal_spec()` 只放行 SUFFIX 表达）——纯醚/纯硫醚类分子会落入纯烃兜底。

### 2.3 互斥检查（无需改动）

无 `_no_fgs` 互斥谓词（无 `fg_helpers.py`）。互斥由 `select_principal_group` 的结构性单选择实现：P-44 只选单个最高优先级主官能团，低优先级 FG 一律成为取代基。新增 FG 无需配置互斥 keys。

### 2.4 添加 ownership 逻辑

如果新 FG 包含不在 chain 中的 heteroatom（如 thiol 的 S），需要在 `src/namepredict/layer2/parent_ownership.py` 中添加 `_kind_fg_atoms` 的分支。以 thiol 为例（`_thiol_fg_atoms` `:171-183`）：

```python
def _thiol_fg_atoms(mol: Mol, parent: dict) -> set[int]:
    """Thiol: attachment carbon + SH sulfur."""
    c_idx = parent.get("sh_c_idx")
    if c_idx is None:
        return set()
    out = {int(c_idx)}
    for n in mol.GetAtomWithIdx(int(c_idx)).GetNeighbors():
        if n.GetAtomicNum() == 16:
            out.add(n.GetIdx())
    return out
```

然后在 `_kind_fg_atoms`（`:254-277`）的 parts tuple 中添加调用。

磷酸的做法（`_phosphate_fg_atoms`，`parent_ownership.py:242`）：`kind == "phosphate"` 时把 P 中心与其全部 O 邻居纳入 `owned_atoms`，O–R 臂的碳留在边界外由 L3 提取。若新 FG 需要向 L5 传递整名所需的计数或盐元数据，还要在 `principal_expression.py` 加 producer 分支（如 `_chain_phosphate_fields`，`:279`，含盐门控）并在 L5 `_names_for` 加 worker 分支（见步骤 5.1 方式 B）。

---

## 步骤 3: Layer 3 — 取代基提取 (Substituent Extraction)

如果新 FG 可以作为取代基前缀出现（大多数 FG 都可以），需要：

1. **保留取代名注册**：在 `src/namepredict/tools/anchored_table.py` 的 `_REGISTRY` 中添加 `RetainedSubstituent` 条目，锚定 canonical-SMILES 写在条目 `anchored` 字段（如 `"sulfanyl": RetainedSubstituent(..., anchored=("*S",))`）。只有**非保留**取代基才在 `ANCHOR_TABLE` 内联 `(en, zh, paren, kind)`。构建期校验会检查 anchored 键的 canonical 形式与唯一性，`*C=C-C` 类死条目会立即报错。
2. **前缀命名**：由 `anchored_table` 查表或 `SubstituentNamer` 兜底完成，通常无需额外代码
3. **括号规则**：复合前缀（词干带取代）按 PIN P-16.5.1.1 加围栏，简单前缀免括（`_radical_yl_from_sub`，`as_substituent.py:95`）；若取代名以 `oxy`/`sulfanyl` 结尾且前端 R 为简单取代基，`_obridge_front_simple`（`as_substituent.py:48`）会免去围栏（P-63.2.1/.2.2，gold/ChEBI 平铺式）——前端是否简单由命名后端 retained→recursive 判定，无法判定时保守加括号

---

## 步骤 4: Layer 4 — 编号 (Numbering) —— 通常无需改动

**`numbering_engine.orient_numbering`**（P-14.4 候选管线）不再按 kind 枚举 orienter：

1. 枚举候选编号（链正反/环 2n）
2. 按 P-14.4 逐条收窄：principal FG 最低位次集 → 多重键 → 取代基位次集

**新增 FG 无需注册 orienter**——只要 parent dict 携带正确的 FG 定位字段，引擎自动把主官能团放到最低位次。FG 位次由 `locant_calc.py` 的 `_fg_locants` 产出（`_FG_LOCANTS` 由 `fg_registry.FgSpec.locant_kind` 派生）：需记主基团位次的 FG 在 FgSpec 设 `locant_kind` 并在 `locant_calc._LOCANT_FNS` 补位次函数（如 `_aldehyde_fg_locants`/`_acid_fg_locants`）；若新 FG 的位次需要省略规则，在 `omit_locants.py` 中添加。

> **源:** `src/namepredict/layer4/numbering_engine.py`

---

## 步骤 5: Layer 5 — 名称组装 (Name Assembly)

### 5.1 命名实现（当前机制）

Layer5 的母体命名以 `chain_engine.py` 的 **`_KIND_TABLE`** 链引擎为主（数据驱动）：

**方式 A：链状/环状 FG → 在 `_KIND_TABLE` 添加 `_Chain` spec。** 例如 `thiol`（`chain_engine.py`）声明 `en_suf="thiol"`、`zh_suf="硫醇"`、`fg="sh"` 等。`_chain_names` 统一渲染词干、位次与省略规则——新增此类 FG 通常只需一行 spec，无需修改派发逻辑。

**方式 B：特殊拼接 → 在 `_names_for` 添加 worker 分支**（先于 chain_engine 查表尝试）。如 `_exocyclic_acid_names`（`assembler.py:78`）、`_exocyclic_aldehyde_names`（`assembler.py:208`，外环 -CHO）、`_mononuclear_radical_names`（`assembler.py:315`，杂原子锚点自由基）、`kind == "phosphate"` 分支（`assembler.py:370` → `layer5/phosphate.py:phosphate_names`，整名自组装，不经 `_KIND_TABLE`）等。

### 5.2 kind 收敛（无需改动）

kind 收敛在 L2 `_chain_kind`，L5 直接按 FG 类别 kind 查 `_KIND_TABLE`（无 `typed_kinds` 模块）。链式 FG 只需加进 `_CHAIN_FG`；若新 FG 是环/苯组合保留名，在 `_KIND_TABLE` 对应 entry 的 `variant` 中登记 `{"benzene": {1: {plain_fn=...}}}` 特例。

### 5.3 接线到组装流水线

`_names_for` 返回 `(en, zh)` 后，由 `assemble`（`assembler.py:493`）统一完成后续变换——`join_kind_name` 拼接前缀、`maybe_anion_names` 阴离子、`apply_rs_prefix` 立体前缀、`maybe_metal_salt_names` 盐后缀——无需为单个 FG 手动接线。若新 FG 自行组装盐形态（如 phosphate），还需在 `namer._apply_salt_suffix`（`namer.py:242`）加跳过条件，避免通用盐后缀重复追加。

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
| L1 | `analyzer.py:_fg_parts`/`_fg_bools` | 检测逻辑 + entries 函数 + 列表键 + 布尔标志 |
| L1 | `analyzer.py` 或新建 `{fg}.py` | 独立检测器（较复杂 FG） |
| L1 | `_carbonyl_common.py` | 共享羰基原语（如涉及 C=O） |
| L1 | `functional_group_inventory.py` | `FunctionalGroupClass` + `_LIST_CLASSES` + `_ANCHOR_KEYS`（如需类型化类别） |
| L2 | `fg_registry.FG_SPECS` → `principal.py:PRINCIPAL_REGISTRY`（派生） | `FgSpec`（`p41`/`compat`/`expr`/`chain`/`multi`/`rs`/`locant_kind`；非 SUFFIX 用 `legacy_compat`） |
| L2 | `kind_registry.py` | 一般无需改动——只注册 scaffold 词干，FG 主基团等级不存于此 |
| L2 | `principal_expression.py:_CHAIN_FG` | 链式 FG 类别集合（SUFFIX 类） |
| L2 | `parent_ownership.py:_kind_fg_atoms` | FG heteroatom ownership 函数 |
| L2 | `principal_expression.py` | producer 字段（整名计数/盐门控等，如需；如 `_chain_phosphate_fields`:279） |
| L3/tools | `tools/anchored_table.py` | 前缀名注册（如需） |
| L4 | `omit_locants.py` | 省略位次规则（如需；numbering_engine 无需注册） |
| L5 | `chain_engine.py:_KIND_TABLE` | `_Chain` spec |
| L5 | `assembler.py:_names_for` | worker 分支（特殊拼接） |
| L5 | `layer5/{fg}.py` | 整名 worker（如需，如 `phosphate.py`） |
| Namer | `namer.py:_apply_salt_suffix` | 自行组装盐形态时的跳过条件（如需） |

---

## 常见陷阱

1. **忘记添加 `has_*` 布尔键**：如果 `_fg_bools` 中没有映射，L2 无法感知该 FG，母体选择会漏掉它。
2. **fg_rank 设置错误**：fg_rank 决定评分优先级与 suffix 资格。SUFFIX 类才能成为主官能团。
3. **ownership 遗漏**：如果 FG 包含 heteroatom（如 S, O, N），必须确保这些原子被标记为 `owned_atoms`。遗漏会导致这些原子被 Layer 3 误识别为"未被 parent 拥有"的取代基碎片。
4. **数量派生 kind 不要造**：数量由 `multiplicity` 承载，`_KIND_TABLE` 的 `variant` 字段切换后缀——不要注册 `dithiol`/`trithiol` 之类的新 kind。
5. **中文命名约定不一致**：中文 stems 应遵循 IUPAC 中文命名规范（CCS 规则）。
6. **整名 worker 与盐后缀重复**：自行组装盐形态的 FG（phosphate）必须同时在 `namer._apply_salt_suffix` 跳过通用盐后缀，否则金属盐会叠加两次。
7. **L1 列表键漏进 `_fg_more_lists`**：新列表键若不在 `_fg_lists` 的 5 个核心键内，`_fg_lists` 不会把它带进 info dict（`phosphates` 即靠 `_fg_more_lists` 收录）。

---

## 相关页面

- [[architecture/layer1-analyzer]] -- Layer 1 官能团检测器架构
- [[architecture/layer2-parent-selector]] -- Layer 2 母体选择器架构
- [[architecture/layer3-substituents]] -- Layer 3 取代基提取
- [[architecture/layer4-numbering]] -- Layer 4 编号引擎
- [[architecture/layer5-name-assembly]] -- Layer 5 名称组装
- [[architecture/overview]] -- 6 层管线架构总览
- [[concepts/functional-group-priority]] -- FG 优先级体系与 fg_rank 详解
