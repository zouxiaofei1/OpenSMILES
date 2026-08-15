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
- **较复杂 FG**（如 isocyanate, acyl halide）：可创建独立模块 `src/namepredict/layer1/{fg_name}.py`，导出 `{fg}_entries(mol) -> list[dict]` 函数。当前仅剩 `isocyanate.py` 与 `acyl_halide.py` 两个独立模块。
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

2. 在 `_fg_parts`（`analyzer.py:362-373`）中添加条目列表键（如 `"thiols": _thiol_entries(mol)`）
3. 在 `_fg_bools`（`analyzer.py:352-359`）中添加布尔标志（如 `"has_thiol"`）
4. 如果新 FG 需要类型化的 `FunctionalGroupClass`，在 `functional_group_inventory.py` 中添加：枚举成员 + `_LIST_CLASSES` 键映射 + `_ANCHOR_KEYS` 锚点键（thiol 已存在：`THIOL` + `"thiols"` → `s_idx`）

---

## 步骤 2: Layer 2 — 母体选择 (Parent Selection)

### 2.1 注册 fg_rank

`fg_rank` 的单一权威是 `src/namepredict/layer2/principal.py` 的 `PRINCIPAL_REGISTRY`（FG → `PrincipalFeatureSpec.compatibility_rank`，遵循 IUPAC P-41 顺序）。`kind_registry` 通过 `_principal_rank`（`kind_registry.py:18`）把 kind 字符串转 FG 枚举后取 `legacy_rank` 实时投影，无预注册表（`_KIND_CLASS`/`_load_chain_fg` 不存在）。新增 chain FG 只需：

1. 在 `PRINCIPAL_REGISTRY` 中注册该 FG 的 `PrincipalFeatureSpec`。SUFFIX 档用 `_suffix(p41_class, rank, *path)`（如 `FG.THIOL: _suffix(17, 4, 2)`）；非 SUFFIX 类（如 sulfide/ether）用 `PrincipalFeatureSpec(..., PrincipalExpression.LEGACY_COMPAT/PREFIX_ONLY, rank)`

### 2.2 母体接线（按表达权限分档）

**SUFFIX 档**（acid/ester/amide/nitrile/aldehyde/ketone/alcohol/thiol/amine）：`principal_expression.py` 的 `_CHAIN_FG`（`:41`）限定 8 类可链式表达的 FG，`_chain_kind`（`:65`）按类别与多重度返回 kind：

```python
# _CHAIN_FG 含该 FG 类即自动支持；_chain_kind 返回 group_class.value
```

> ACID/ALCOHOL/AMINE/KETONE 任意 count≥1 恒返回基团名（`_MULTI_FG`，数量由 `principal_expression_facts.multiplicity` 承载）；ESTER/AMIDE/NITRILE/ALDEHYDE 仅单基（count≠1 → None）。`dione` 不再由 L2 产生，二酮由 L5 chain_engine `mult_ok` 生成式命名。新增 SUFFIX FG 若属同类，把其类别加进 `_CHAIN_FG` 即可。

**LEGACY_COMPAT 档**（sulfide 等非 SUFFIX 类）：`compatibility_rank` 经 `_principal_rank` 投影为 `fg_rank`，但**不参与主官能团选择**（`principal_spec()` 只放行 SUFFIX）——纯醚/纯亚砜类分子会落入纯烃兜底。

### 2.3 互斥检查（无需改动）

无 `_no_fgs` 互斥谓词（无 `fg_helpers.py`）。互斥由 `select_principal_group` 的结构性单选择实现：P-44 只选单个最高优先级主官能团，低优先级 FG 一律成为取代基。新增 FG 无需配置互斥 keys。

### 2.4 添加 ownership 逻辑

如果新 FG 包含不在 chain 中的 heteroatom（如 thiol 的 S），需要在 `src/namepredict/layer2/parent_ownership.py` 中添加 `_kind_fg_atoms` 的分支。以 thiol 为例（`_thiol_fg_atoms` `:145-154`）：

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

然后在 `_kind_fg_atoms`（`:203-224`）的 parts tuple 中添加调用。

---

## 步骤 3: Layer 3 — 取代基提取 (Substituent Extraction)

如果新 FG 可以作为取代基前缀出现（大多数 FG 都可以），需要：

1. **保留取代名注册**：在 `src/namepredict/tools/anchored_table.py` 的 `_REGISTRY` 中添加 `RetainedSubstituent` 条目，锚定 canonical-SMILES 写在条目 `anchored` 字段（如 `"sulfanyl": RetainedSubstituent(..., anchored=("*S",))`）。只有**非保留**取代基才在 `ANCHOR_TABLE` 内联 `(en, zh, paren, kind)`。构建期校验会检查 anchored 键的 canonical 形式与唯一性，`*C=C-C` 类死条目会立即报错。
2. **前缀命名**：由 `anchored_table` 查表或 `SubstituentNamer` 兜底完成，通常无需额外代码

---

## 步骤 4: Layer 4 — 编号 (Numbering) —— 通常无需改动

**`numbering_engine.orient_numbering`**（P-14.4 候选管线）不再按 kind 枚举 orienter：

1. 枚举候选编号（链正反/环 2n）
2. 按 P-14.4 逐条收窄：principal FG 最低位次集 → 多重键 → 取代基位次集

**新增 FG 无需注册 orienter**——只要 parent dict 携带正确的 FG 定位字段，引擎自动把主官能团放到最低位次。FG 位次由 `locant_calc.py` 的 `_FG_LOCANTS` 表产出；若新 FG 的位次需要省略规则，在 `omit_locants.py` 中添加。

> **源:** `src/namepredict/layer4/numbering_engine.py`

---

## 步骤 5: Layer 5 — 名称组装 (Name Assembly)

### 5.1 命名实现（当前机制）

Layer5 的母体命名以 `chain_engine.py` 的 **`_KIND_TABLE`** 链引擎为主（数据驱动）：

**方式 A：链状/环状 FG → 在 `_KIND_TABLE` 添加 `_Chain` spec。** 例如 `thiol`（`chain_engine.py`）声明 `en_suf="thiol"`、`zh_suf="硫醇"`、`fg="sh"` 等。`_chain_names` 统一渲染词干、位次与省略规则——新增此类 FG 通常只需一行 spec，无需修改派发逻辑。

**方式 B：特殊拼接 → 在 `_names_for` 添加 worker 分支。** 如 `_exocyclic_acid_names`（`assembler.py:55-76`）、`_sulfide_names` 等。

### 5.2 kind 收敛（无需改动）

kind 收敛在 L2 `_chain_kind`，L5 直接按 FG 类别 kind 查 `_KIND_TABLE`（无 `typed_kinds` 模块）。链式 FG 只需加进 `_CHAIN_FG`；若新 FG 是环/苯组合保留名，在 `_KIND_TABLE` 对应 entry 的 `variant` 中登记 `{"benzene": {1: {plain_fn=...}}}` 特例。

### 5.3 接线到组装流水线

`_names_for` 返回 `(en, zh)` 后，由 `assemble`（`assembler.py:196`）统一完成后续变换——`join_kind_name` 拼接前缀、`maybe_anion_names` 阴离子、`apply_rs_prefix` 立体前缀、`maybe_metal_salt_names` 盐后缀——无需为单个 FG 手动接线。

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
| L2 | `principal.py:PRINCIPAL_REGISTRY` | `PrincipalFeatureSpec`（`_suffix` 或 `LEGACY_COMPAT`） |
| L2 | `kind_registry.py` | fg_rank 投影（`_principal_rank`） |
| L2 | `principal_expression.py:_CHAIN_FG` | 链式 FG 类别集合（SUFFIX 类） |
| L2 | `parent_ownership.py:_kind_fg_atoms` | FG heteroatom ownership 函数 |
| L3/tools | `tools/anchored_table.py` | 前缀名注册（如需） |
| L4 | `omit_locants.py` | 省略位次规则（如需；numbering_engine 无需注册） |
| L5 | `chain_engine.py:_KIND_TABLE` | `_Chain` spec |
| L5 | `assembler.py:_names_for` | worker 分支（特殊拼接） |

---

## 常见陷阱

1. **忘记添加 `has_*` 布尔键**：如果 `_fg_bools` 中没有映射，L2 无法感知该 FG，母体选择会漏掉它。
2. **fg_rank 设置错误**：fg_rank 决定评分优先级与 suffix 资格。SUFFIX 类才能成为主官能团。
3. **ownership 遗漏**：如果 FG 包含 heteroatom（如 S, O, N），必须确保这些原子被标记为 `owned_atoms`。遗漏会导致这些原子被 Layer 3 误识别为"未被 parent 拥有"的取代基碎片。
4. **数量派生 kind 不要造**：数量由 `multiplicity` 承载，`_KIND_TABLE` 的 `variant` 字段切换后缀——不要注册 `dithiol`/`trithiol` 之类的新 kind。
5. **中文命名约定不一致**：中文 stems 应遵循 IUPAC 中文命名规范（CCS 规则）。

---

## 相关页面

- [[architecture/layer1-analyzer]] -- Layer 1 官能团检测器架构
- [[architecture/layer2-parent-selector]] -- Layer 2 母体选择器架构
- [[architecture/layer3-substituents]] -- Layer 3 取代基提取
- [[architecture/layer4-numbering]] -- Layer 4 编号引擎
- [[architecture/layer5-name-assembly]] -- Layer 5 名称组装
- [[architecture/overview]] -- 6 层管线架构总览
- [[concepts/functional-group-priority]] -- FG 优先级体系与 fg_rank 详解
