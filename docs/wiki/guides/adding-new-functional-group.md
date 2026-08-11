# 添加新官能团 (Adding a New Functional Group)

> **扩展指南** | 关联: [[architecture/overview]], [[architecture/layer1-analyzer]], [[architecture/layer2-parent-selector]], [[architecture/layer3-substituents]], [[architecture/layer4-numbering]], [[architecture/layer5-name-assembly]], [[concepts/functional-group-priority]]

---

## 概述

本指南详细说明如何向 NamePredict 命名管线添加一个新的 functional group（官能团，FG）类型。整个流程跨越全部 5 个 layer，从前端检测到末端名称组装。我们以 **sulfoxide（亚砜，R-S(=O)-R'）** 作为完整的工作示例——它是一个中等优先级（`fg_rank=2`）、非末端的链状 FG，在 IUPAC P-63.3 中有明确定义，代码实现简洁清晰，非常适合作为新增 FG 的参考模板。

在开始之前，建议先阅读 [[concepts/functional-group-priority]] 了解 FG 优先级体系和互斥检查规则，以及 [[architecture/overview]] 了解 5 层管线架构。

---

## 步骤 1: Layer 1 — 检测 (Detection)

Layer 1 的职责是从 RDKit `Mol` 对象中检测官能团，输出结构化的 atom index 条目和布尔标志，供下游使用。

### 决策: 内联检测 vs. 独立模块

根据 FG 的复杂度选择合适的检测方式：

- **简单 FG**（如 sulfide, ether, thiol）：直接在 `analyzer.py` 中内联实现，通常只需要 3 个函数：候选原子判断、条目构建、收集函数。
- **复杂 FG**（如 carbamate, sulfoxide, boronic）：创建独立模块 `src/namepredict/layer1/{fg_name}.py`，导出 `{fg}_entries(mol: Mol) -> list[dict]` 函数，可额外导出 `is_{fg}_atom(atom) -> bool` 供其他检测器做排他判断。

sulfoxide 属于中等复杂度——有自己独特的原子环境约束（S, degree=3, 一个 =O, 两个 C），因此采用独立模块。

### 选项 A: analyzer.py 内联检测

适用于检测逻辑简单（<20 行）且不需要被其他检测器排他引用的 FG。模式如下：

1. 定义候选原子判断函数：

```python
def _is_{fg}_atom(atom) -> bool:
    # 检查原子序数、度数、键型、邻居类型
    ...
```

2. 定义条目构建函数：

```python
def _fg_entry(atom) -> dict:
    # 打包 atom.GetIdx() 及相关邻居索引
    return {"key_idx": atom.GetIdx(), ...}
```

3. 定义收集函数并注册到 `_collect_fgs`：

```python
def _{fg}_entries(mol: Mol) -> list[dict]:
    return [_fg_entry(a) for a in mol.GetAtoms() if _is_{fg}_atom(a)]
```

4. 在 `_FG_BOOL_MORE_KEYS`（`analyzer.py:5-20`）中添加映射：

```python
("has_{fg}", "{fg}s"),
```

5. 根据复杂度，将收集函数注册到 `_fg_parts_b_core`、`_p_fg_a`、`_p_fg_b1` 或 `_p_fg_b2` 中。

### 选项 B: 独立模块检测（sulfoxide 的模式）

创建 `src/namepredict/layer1/sulfoxide.py`（`sulfoxide.py:1-39`）：

```python
"""L1 detection of open-chain dialkyl sulfoxides (P-63.3)."""
from __future__ import annotations
from rdkit.Chem import BondType, Mol

def _is_sulfoxide_sulfur(atom) -> bool:
    """S with exactly one =O and two C neighbors (not sulfone/sulfide)."""
    if atom.GetAtomicNum() != 16 or atom.GetTotalNumHs() != 0:
        return False
    if atom.GetTotalDegree() != 3:
        return False
    if len(_dbl_o_nbs(atom)) != 1:
        return False
    return len(_c_nbs(atom)) == 2

def sulfoxide_entries(mol: Mol) -> list[dict]:
    return [_sulfoxide_entry(a) for a in mol.GetAtoms()
            if _is_sulfoxide_sulfur(a)]
```

检测逻辑的关键约束：

- `GetAtomicNum() == 16`：必须是硫原子
- `GetTotalDegree() == 3`：排除 sulfide（degree=2）和 sulfone（degree=4）
- 恰好一个双键 O：区分 sulfoxide 和 sulfone（两个双键 O）
- 恰好两个 C 邻居：确保是 R-S(=O)-R' 模式

### 注册到 analyzer.py

在 `analyzer.py` 中：

1. **导入并注册收集函数**（`analyzer.py:438-441`）：sulfoxide 被放入 `_p_fg_a` 中，与其他中等复杂度 FG（phosphate, carbamate, carbonate）同组：

```python
def _p_fg_a(mol: Mol) -> dict:
    from namepredict.layer1.sulfoxide import sulfoxide_entries
    return {"phosphates": phosphate_entries(mol), ...,
            "sulfoxides": sulfoxide_entries(mol)}
```

2. **添加 `has_*` 布尔键映射**（`analyzer.py:13`）：

```python
_FG_BOOL_MORE_KEYS = (
    ...
    ("has_sulfoxide", "sulfoxides"),
    ...
)
```

3. **确保键名在 `_fg_more_lists` 函数中出现**（`analyzer.py:405-414`），该函数定义了哪些键从 parts dict 传递到 info dict。

4. **如果新 FG 需要被更高优先级 FG 排他排除**（如 carboxyl 排除 ester），在高优先级检测器中添加对 `is_{fg}_atom()` 的调用。sulfoxide 在 fg_rank=2，低于大多数常见 FG，因此不需要被其他检测器排他排除——较低的 fg_rank 足以通过 Layer 2 的互斥检查（`_no_fgs`）将其降级为取代基。

---

## 步骤 2: Layer 2 — 母体选择 (Parent Selection)

Layer 2 是管线中最复杂的部分，负责根据检测到的 FG 集合选择母体结构（parent）。添加新 FG 需要触及多个文件。

### 2.1 注册 kind 和 fg_rank

`fg_rank` 的单一权威是 `src/namepredict/layer2/principal.py` 的 `PRINCIPAL_REGISTRY`
（FG → `PrincipalFeatureSpec.compatibility_rank`，遵循 IUPAC P-41 顺序）；`kind_registry.py`
的 `_KIND_CLASS` 表把 kind 字符串映射到 FG 枚举，`_load_chain_fg()` 遍历其 keys 注册。
新增 chain FG 分两步：

1. 在 `PRINCIPAL_REGISTRY` 中补充该 FG 的 `PrincipalFeatureSpec`，expression 用
   `LEGACY_COMPAT`（非 SUFFIX 类不参与主官能团选择，`principal_spec()` 只放行 SUFFIX）：

```python
FG.SULFOXIDE: PrincipalFeatureSpec(PrincipalPriority(41, (3,)), PrincipalExpression.LEGACY_COMPAT, 2),
```

2. 在 `kind_registry.py` 的 `_KIND_CLASS` 中补 kind→FG 映射：

```python
"sulfoxide": FG.SULFOXIDE,
```

`sulfoxide` 的 `fg_rank=2`（即 `compatibility_rank`），与 ether 和 sulfide 同级——在 IUPAC P-41 中，
它们都是最低优先级的含杂原子 FG。由于 `LEGACY_COMPAT` 档不参与主官能团选择（见 2.2），
NamePredict 中 sulfoxide 不会成为 principal group，仅作为 fg_rank 供评分与互斥使用。

### 2.2 母体接线（按表达权限分档）

新 FG 的母体接入方式取决于它属于 `PRINCIPAL_REGISTRY` 的哪一档（`principal.py:39-71`）。系统
不存在独立的 producer/try 函数模块——母体选择统一走 `rule_driven_parent_candidates` 的 typed
表达。

**SUFFIX 档（acid/ester/amide/nitrile/aldehyde/ketone/alcohol/thiol/amine）**：
在 `PRINCIPAL_REGISTRY` 用 `_suffix(p41_class, rank, *path)` 注册（如 `FG.THIOL: _suffix(17, 4, 2)`），
并在 `principal_expression.py` 的表达表中声明 kind 与字段：
- 开链：`_CHAIN_KINDS`（group_class → kind by multiplicity，如 acid→acid/diacid/polycarboxylic）
- 苯环 + 单 FG：`_RETAINED_RING_KINDS`（如 `ESTER → "benzoate"`）
- 环骨架：`_RING_FIELDS`（位次字段名）+ `_resolved_ring_kind`（scaffold 组合）

**LEGACY_COMPAT 档（sulfoxide/sulfone/carbamate/磷酸等非 SUFFIX 类）**：
在 `PRINCIPAL_REGISTRY` 注册 `PrincipalFeatureSpec`（`LEGACY_COMPAT` + `compatibility_rank`，
如 `FG.SULFOXIDE: PrincipalFeatureSpec(PrincipalPriority(41, (3,)), PrincipalExpression.LEGACY_COMPAT, 2)`），
并在 `kind_registry._KIND_CLASS`（`kind_registry.py:10-33`）补 kind→FG 映射。该档
`compatibility_rank` 经 `_principal_rank` 投影为 `fg_rank`，但**不参与主官能团选择**
（`principal_spec()` 只放行 SUFFIX）——纯醚/纯亚砜分子会落入纯烃兜底，属已知能力缺口。

### 2.3 互斥 keys（如果需要新组合）

如果新 FG 作为母体时不能与其他 FG 共存，通过 `_no_fgs(info, keys)` 传参即可（`fg_helpers.py:26`）：

```python
keys = ("has_acid", "has_ester", "has_alcohol", "has_amine", ...)
ok = _no_fgs(info, keys)
```

FG 特有的互斥 keys 可直接内联在该 FG 所在模块，不需要被其他模块共享。`_no_fgs` 是唯一的共享
互斥原语。

### 2.4 添加 ownership 逻辑

如果新 FG 包含不在 chain 中的 heteroatom（如 sulfoxide 的 S 和 O），需要在 `src/namepredict/layer2/parent_ownership.py` 中添加 `_kind_fg_atoms` 的分支，确保这些原子被标记为 `owned_atoms`。以 sulfoxide 为例，需要添加类似 `_sulfide_fg_atoms` 的函数：

```python
def _sulfoxide_fg_atoms(mol: Mol, parent: dict) -> set[int]:
    """Sulfoxide S + =O + both carbon arms."""
    if parent.get("kind") != "sulfoxide" or parent.get("s_idx") is None:
        return set()
    s_idx = int(parent["s_idx"])
    out = {s_idx}
    for n in mol.GetAtomWithIdx(s_idx).GetNeighbors():
        if n.GetAtomicNum() == 8:
            out.add(n.GetIdx())      # =O atom
        elif n.GetAtomicNum() == 6:
            out |= _ether_arm_atoms(mol, s_idx, n.GetIdx())  # carbon arms
    return out
```

然后将其添加到 `_kind_fg_atoms` 的 parts tuple 中（与 `_sulfide_fg_atoms` 模式一致）。

---

## 步骤 3: Layer 3 — 取代基提取 (Substituent Extraction)

Layer 3 负责将未被 parent 拥有的剩余原子识别为 substituent（取代基）。如果新 FG 可以作为取代基前缀出现（大多数 FG 都可以），需要在相关文件中注册。

### 3.1 保留取代名注册

如果新 FG 作为取代基时有 IUPAC 保留名（retained name），在 `src/namepredict/tools/anchored_table.py` 的 `_REGISTRY` 中添加 `RetainedSubstituent` 条目（该注册表同时承载 anchored 查表的保留名解析）。

sulfoxide 作为取代基的前缀名是 "-sulfinyl"（亚磺酰基），但目前 NamePredict 中 sulfoxide 的 fg_rank=2 意味着它几乎总是被更高优先级的 FG 压制为取代基前缀。如果将来需要显式注册 sulfinyl 前缀，在此处添加。

### 3.2 取代基命名逻辑

如果新 FG 作为取代基时有特殊的命名模式（如 -sulfinyl vs. 作为母体时的 -sulfoxide），在 `src/namepredict/layer3/substituent_extractor.py` 或 `substituent_namer.py` 中添加对应的提取和命名逻辑。

---

## 步骤 4: Layer 4 — 编号 (Numbering)

Layer 4 负责为 parent chain 分配位次编号。

### 4.1 添加 orienter 函数

在 `src/namepredict/layer4/numbering.py` 中，根据 FG 类型选择合适的 orienter：

- **末端 FG**（acid, aldehyde, nitrile, ester, amide）：使用 `_orient_to_terminal(chain, c_idx)`，确保 FG 碳在链的 1 号位。
- **非末端单 FG**（ketone, alcohol, amine, thiol）：使用 `_orient_by_single_fg(chain, parent, substituents, key)`，将 FG 碳放到尽可能低的位置。
- **双 FG**（diol, dione, diamine）：使用 `_orient_pair` 或自定义逻辑。

sulfoxide 不是 terminal FG——它的 S 原子位于链内，两侧各有一条碳臂。由于 sulfoxide parent 的 chain 已经是较长臂（由 `_make_parent` 选择），且 S 原子在 chain 的端点之一（通过 `chain[0]` 可达），默认的 `_orient_alkane` 即可处理——不需要特殊的 orienter。

但如果新 FG 需要特殊的 orienter（如将某个特征原子放在 1 号位），需要在 `_kind_orienters()` 字典（`numbering.py:327`）中注册。sulfoxide 不需要此项，因此跳过。

### 4.2 省略位次规则

如果新 FG 的位次在某些情况下可以省略（如 ethanamine 不需要标注 1-amine），在 `src/namepredict/layer4/omit_locants.py` 中添加规则。

sulfoxide 的对称命名（如 dimethyl sulfoxide）不需要位次编号，这一逻辑由 Layer 5 的对称名称字典处理，而非 omit_locants。

---

## 步骤 5: Layer 5 — 名称组装 (Name Assembly)

Layer 5 将编号后的 parent 和 substituent 信息组装为最终的中英双语名称。

### 5.1 命名实现（当前机制）

Layer5 的母体命名以 `assembler.py` 的 `_KIND_TABLE` 链引擎为主（数据驱动），特殊拼接模式在 `_names_for` 中添加 worker 分支：

**方式 A：链状/环状 FG → 在 `_KIND_TABLE` 添加 `_Chain` spec。** 例如 `acid`（`assembler.py:363-369`）声明 `en_suf="oic acid"`、`zh_suf="酸"`、不饱和段 `ene_base=("enoic acid","烯酸")`/`yne_suf=("ynoic acid","炔酸")`、E/Z 包装 `ez_ene=_ez_prefix`。醇/酮/胺/醛/腈/烯/炔/环烷等 29 种链/环 kind 均由 `_chain_names` 统一渲染词干、不饱和段、位次与环前缀——新增此类 FG 通常只需一行 spec，无需修改派发逻辑。

**方式 B：特殊拼接 → 在 `_names_for` 添加 worker 分支。** 例如硫醚（`assembler.py:482-483`）：

```python
if kind == "sulfide":
    return _sulfide_names(numbered)
```

`_sulfide_names`（`:49`）按对称/不对称臂长分派 `_sym_sulfide_names`/`_asym_sulfide_names`，返回 `(en, zh)` 或 `None`（None 触发 fallback）。

**复用 stems**：烷基词干从 `stems.py` 的 `ALKANE_EN/ZH`、`ester_alkyl_*` 等表取，避免重复定义。对称/不对称模式（如 `dimethyl sulfide`/`二甲基硫醚` vs `ethyl methyl sulfide`/`乙基甲基硫醚`）复用同一批 alkyl 词干。

### 5.2 接线到组装流水线

`_names_for` 返回 `(en, zh)` 后，由 `assemble`（`assembler.py:510`）统一完成后续变换——`join_kind_name` 拼接前缀、`maybe_anion_names` 阴离子、`apply_rs_prefix` 立体前缀、`maybe_metal_salt_names` 盐后缀——无需为单个 FG 手动接线。

如果新 FG 需要专属拼接（如酯/苯甲酸酯在 alkyl 与 acyl 之间插入前缀），在 `benzene_names.py` 增加 `join_{fg}_name`，并在 `join_kind_name`（`benzene_names.py:189-199`）添加 kind 分支。

---

## 测试

每个新 FG 必须通过以下测试验证：

### 基础测试用例

1. **简单 FG 作为母体**：仅含该 FG 的简单分子。例如 `CS(=O)C`（dimethyl sulfoxide）验证 sulfoxide 的基本命名。
2. **复杂 FG 作为母体**：含不同臂长的分子。例如 `CS(=O)CC`（ethyl methyl sulfoxide）验证不对称命名。
3. **FG 作为取代基**：该 FG 与更高优先级 FG 共存。例如含有 -COOH 和 -S(=O)- 的分子，验证 sulfoxide 降级为 sulfinyl 前缀。
4. **与其他取代基共存**：在基本 FG 上叠加 alkyl/halogen 等取代基，验证编号和位次正确。
5. **中英双语**：所有测试用例必须验证 en 和 zh 两种输出。

### Coverage Ledger 验证

- **无 gap**（无遗漏）：该 FG 的所有合理分子变体都应被正确命名。
- **无 overlap**（无冲突）：新 FG 不应破坏已有 FG 的命名结果（通过现有测试回归验证）。
- **dual coverage**：验证新 FG 的 kind 在 principal 表达表（`_CHAIN_KINDS` / `_RETAINED_RING_KINDS`）或 `kind_registry` 中的注册不会导致候选选择冲突。

### 集成测试

- 验证新 FG 作为 **parent**（母体）时的命名链完整。
- 验证新 FG 作为另一 FG 的 **substituent**（取代基）时的命名链完整。

---

## 需修改文件汇总表

| Layer | 文件 | 添加内容 |
|:---:|---|------|
| L1 | `analyzer.py` 或新建 `{fg}.py` | 检测逻辑 + entries 函数 |
| L1 | `analyzer.py:_FG_BOOL_MORE_KEYS` | `has_*` 布尔键映射 |
| L1 | `analyzer.py:_p_fg_a/b1/b2` 或 `_fg_parts_b_core` | 注册 entries 收集函数 |
| L2 | `principal.py:PRINCIPAL_REGISTRY` | `PrincipalFeatureSpec`（`LEGACY_COMPAT` + `compatibility_rank`，或 `_suffix`） |
| L2 | `kind_registry.py:_KIND_CLASS` | kind→FG 映射 |
| L2 | `principal_expression.py`（SUFFIX 类） | `_CHAIN_KINDS` / `_RETAINED_RING_KINDS` / `_RING_FIELDS` 表达表 |
| L2 | `fg_helpers.py` | `_no_fgs` 互斥检查（如需要） |
| L2 | `parent_ownership.py:_kind_fg_atoms` | FG heteroatom ownership 函数 |
| L3 | `substituent_extractor.py` / `substituent_namer.py` | 取代基提取/命名规则（如果需要） |
| L3/tools | `tools/anchored_table.py` | 前缀名注册（如果需要） |
| L4 | `numbering.py:_kind_orienters()` | orienter 函数（如果需要） |
| L4 | `omit_locants.py` | 省略位次规则（如果需要） |
| L5 | 新建 `layer5/{fg}_names.py` 或加进 `assembler.py` | 命名模块 / worker（en + zh stems） |
| L5 | `assembler.py:_KIND_TABLE` 或 `_names_for` | `_Chain` spec 或 worker 分支 |

---

## sulfoxide 的当前实现位置（参考）

sulfoxide 功能在当前源码中的位置：

| 文件 | 功能 |
|---|---|
| `src/namepredict/layer1/sulfoxide.py` | 检测 S(=O) 环境 |
| `src/namepredict/layer1/analyzer.py` | 注册 has_sulfoxide 布尔键 + entries 收集 |
| `src/namepredict/layer2/principal.py` | `PRINCIPAL_REGISTRY` 的 `FG.SULFOXIDE`（LEGACY_COMPAT, compatibility_rank=2） |
| `src/namepredict/layer2/kind_registry.py` | `_KIND_CLASS` 的 kind→FG 映射 + bootstrap 注册 |
| `src/namepredict/layer5/assembler.py` | `_KIND_TABLE` 链引擎 + `_sulfide_names` 等 worker（sulfoxide 专属命名模块 `sulfur_names.py` 已删除，L5 目前不输出 sulfoxide 名称） |

---

## 常见陷阱

1. **忘记添加 `has_*` 布尔键映射**：这是最常见的遗漏。如果 `_FG_BOOL_MORE_KEYS` 中没有映射，Layer 2 的 `_no_fgs` 互斥检查无法正确排除该 FG，导致互斥规则失效。

2. **fg_rank 设置错误**：fg_rank 不仅影响评分优先级，还影响互斥规则语义。如果将 sulfoxide 的 fg_rank 设为 6（与 ketone 同级），它就会错误地在含酮分子中成为 principal FG。

3. **ownership 遗漏**：如果 FG 包含 heteroatom（如 S, O, N），必须确保这些原子被标记为 `owned_atoms`。遗漏会导致这些原子被 Layer 3 误识别为"未被 parent 拥有"的取代基碎片。

4. **注册顺序/表达表错误**：新 FG 加入 principal 表达表（`_CHAIN_KINDS` / `_RETAINED_RING_KINDS`）时需按 fg_rank 语义排布，避免高优先级 FG 被低优先级 FG 抢先；修改后运行 dual coverage 测试。

5. **arm_pair 中的 `_arm_ok` 检查遗漏**：对于双臂 FG（sulfide, ether, sulfoxide），必须调用 `_arm_ok(mol, arm, hetero_idx)` 确保每条臂都不含环或其他 heteroatom。遗漏此检查会导致含环结构的分子被错误地选择为 open-chain parent。

6. **中文命名约定不一致**：中文 stems 应遵循 IUPAC 中文命名规范（CCS 规则）。例如 sulfoxide 的对称命名用"二甲基亚砜"而非"二甲基亚磺酰"，因为中文化学命名传统上对称 dialkyl sulfoxide 使用"亚砜"后缀。

---

## 相关页面

- [[architecture/layer1-analyzer]] -- Layer 1 官能团检测器架构
- [[architecture/layer2-parent-selector]] -- Layer 2 母体选择器架构
- [[architecture/layer3-substituents]] -- Layer 3 取代基提取
- [[architecture/layer4-numbering]] -- Layer 4 编号引擎
- [[architecture/layer5-name-assembly]] -- Layer 5 名称组装
- [[architecture/overview]] -- 6 层管线架构总览
- [[concepts/functional-group-priority]] -- FG 优先级体系与 fg_rank 详解
