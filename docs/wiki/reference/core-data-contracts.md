# 核心数据契约 (Core Data Contracts)

> 本文档定义 NamePredict 6 层流水线中流转的 **全部关键数据结构** 的字段规格。每个数据结构均标注源文件与行号，便于代码导航。

---

## 数据流全景

```mermaid
graph TD
    SMILES["SMILES 输入"] --> L0["Layer 0<br/>预处理 / 盐拆分"]
    L0 -->|"Mol (organic)"| L1["Layer 1<br/>analyze(mol)"]
    L1 -->|"info dict"| L2["Layer 2<br/>select_parent(info)"]
    L2 -->|"parent dict"| L3["Layer 3<br/>extract_substituents(info, parent)"]
    L3 -->|"list[subst dict]"| L4["Layer 4<br/>number(parent, subs)"]
    L4 -->|"numbered dict"| L5["Layer 5<br/>assemble(numbered)"]
    L5 -->|"NameResult"| OUT["输出: {en, zh}"]

    L2 -->|"iter_claims(mol, owned_atoms)"| CB["ClaimedBlock"]
    CB -->|"try_name(mol, claim)"| SN["SubstituentName"]
    SN -->|"build_coverage_ledger()"| CL["CoverageLedger"]
    CL -->|"complete?"| NAMER["namer.py<br/>_ok_result()"]

    style SMILES fill:#e1f5fe
    style OUT fill:#c8e6c9
    style CL fill:#fff9c4
    style NAMER fill:#fff9c4
```

**流程说明**：

1. SMILES 经 Layer 0 预处理后得到有机部分的 `Mol` 对象
2. Layer 1 分析 Mol，产出 **info dict** —— 分子中所有官能团和环系的"一次性全景快照"
3. Layer 2 基于 info 选择母体，产出 **parent dict**；同时通过 `iter_claims()` 将 parent 未覆盖的原子切分为 **ClaimedBlock** 列表
4. Layer 3 为每个 ClaimedBlock 命名，产出 **SubstituentName**；同时从 info + parent 提取内置取代基，产出 **subst dict** 列表；通过 **CoverageLedger** 验证全分子覆盖
5. Layer 4 对 parent + subs 进行编号，产出 **numbered dict** —— 带完整位次信息的增强 parent + 带 locant 的 subs
6. Layer 5 将 numbered dict 组装为最终中英双语 **NameResult**

---

## 1. info dict（L1 输出）

**全称**: Functional Group Information Dictionary（官能团信息字典）

**产出**: `analyze(mol)` in `src/namepredict/layer1/analyzer.py:484-485`

**构建**: `_info(mol, carbons, fgs)` at `analyzer.py:481-483`，由三层合并：
- `base`: `mol`, `carbon_ids`, `n_carbons`
- `fgs`: 所有 FG 条目列表 + 布尔标志（`_collect_fgs` at `analyzer.py:478-480`）
- `_ring_meta(mol)`: 环系元信息（`analyzer.py:393-400`）

**输入方**: Layer 2 (`select_parent(info)`), Layer 3 (`extract_substituents(info, parent)`), Layer 4 (间接通过 parent), Layer 5 (间接)

### 基础字段

| 字段 | 类型 | 说明 |
|---|---|---|
| `mol` | `rdkit.Chem.Mol` | 原始 RDKit 分子对象引用（有机部分，已去盐） |
| `carbon_ids` | `list[int]` | 分子中所有碳原子的 atom index |
| `n_carbons` | `int` | 碳原子总数（= `len(carbon_ids)`） |

### 官能团条目列表（FG entry lists）

每个 FG 列表为 `list[dict]`，每项是一个 dict，其字段因 FG 类型而异。以下列出全部列表键名：

| 键名 | 条目 dict 典型字段 | 来源 |
|---|---|---|
| `hydroxyls` | `o_idx`, `c_idx` | `analyzer.py:_hydroxyl_entries` |
| `carboxyls` | `c_idx`, `anion` | `analyzer.py:_carboxyl_entries` |
| `esters` | 酯键原子索引 | `analyzer.py:_ester_entries` |
| `amides` | 酰胺键原子索引 | `analyzer.py:_amide_entries` |
| `ketones` | `c_idx` | `analyzer.py:_ketone_entries` |
| `aldehydes` | `c_idx` | `analyzer.py` |
| `amines` | `n_idx` | `analyzer.py` |
| `nitriles` | `c_idx`, `n_idx` | `analyzer.py` |
| `double_bonds` | 双键原子对 | `analyzer.py` |
| `triple_bonds` | 三键原子对 | `analyzer.py` |
| `acyl_chlorides` | — | `analyzer.py` |
| `anhydrides` | — | `analyzer.py` |
| `thiols` | `s_idx` | `analyzer.py` |
| `ethers` | `o_idx` | `analyzer.py` |
| `sulfides` | `s_idx` | `analyzer.py` |
| `nitros` | — | `analyzer.py` |
| `phosphates` | — | `layer1/phosphate.py` |
| `phosphonics` | — | `layer1/phosphate.py` |
| `carbamates` | — | `layer1/carbamate.py` |
| `carbonates` | — | `layer1/carbonate.py` |
| `sulfoxides` | — | `layer1/sulfoxide.py` |
| `isocyanates` | — | `layer1/isocyanate.py` |
| `isothiocyanates` | — | `layer1/isocyanate.py` |
| `ureas` | — | `layer1/urea.py` |
| `hydrazines` | — | `layer1/hydrazine.py` |
| `guanidines` | — | `layer1/guanidine.py` |
| `sulfonamides` | — | `layer1/sulfonamide.py` |
| `sulfonates` | — | `layer1/sulfonate.py` |
| `sulfonyl_chlorides` | — | `layer1/sulfonyl_chloride.py` |
| `sulfonic_acids` | — | `layer1/sulfonic_acid.py` |
| `sulfones` | — | `layer1/sulfone.py` |
| `boronics` | — | `layer1/boronic.py` |

### 布尔标志（Boolean flags）

对应 `_fg_bools(lists)` at `analyzer.py:420-427`，每个 `has_*` 标志 = `bool(对应的 FG 列表)`：

| 标志 | 对应列表 | 标志 | 对应列表 |
|---|---|---|---|
| `has_alcohol` | `hydroxyls` | `has_acid` | `carboxyls` |
| `has_ester` | `esters` | `has_amide` | `amides` |
| `has_ketone` | `ketones` | `has_aldehyde` | `aldehydes` |
| `has_amine` | `amines` | `has_nitrile` | `nitriles` |
| `has_alkene` | `double_bonds` | `has_alkyne` | `triple_bonds` |
| `has_acyl_chloride` | `acyl_chlorides` | `has_anhydride` | `anhydrides` |
| `has_thiol` | `thiols` | `has_ether` | `ethers` |
| `has_sulfide` | `sulfides` | `has_nitro` | `nitros` |
| `has_phosphate` | `phosphates` | `has_phosphonic` | `phosphonics` |
| `has_carbamate` | `carbamates` | `has_carbonate` | `carbonates` |
| `has_sulfoxide` | `sulfoxides` | `has_isocyanate` | `isocyanates` |
| `has_isothiocyanate` | `isothiocyanates` | `has_urea` | `ureas` |
| `has_hydrazine` | `hydrazines` | `has_guanidine` | `guanidines` |
| `has_sulfonamide` | `sulfonamides` | `has_sulfonate` | `sulfonates` |
| `has_sulfonyl_chloride` | `sulfonyl_chlorides` | `has_sulfonic_acid` | `sulfonic_acids` |
| `has_sulfone` | `sulfones` | `has_boronic` | `boronics` |

> 共计 **31 个布尔标志**（5 个 core + 26 个 from `_FG_BOOL_MORE_KEYS` at `analyzer.py:5-20`）。

### 环系元信息

来自 `_ring_meta(mol)` at `analyzer.py:393-400`：

| 字段 | 类型 | 说明 |
|---|---|---|
| `rings` | `list[dict]` | 每个环的原子索引信息，如 `{"atom_ids": [0,1,2,3,4,5]}` |
| `n_rings` | `int` | 环的数量 |
| `has_ring` | `bool` | 是否含环（= `n_rings > 0`） |
| `ring_systems` | `list[dict]` | 环系（fused ring systems），由 `build_ring_systems(mol)` 产出 |
| `n_ring_systems` | `int` | 环系的数量 |

---

## 2. parent dict（L2 输出）

**全称**: Parent Structure Dictionary（母体结构字典）

**产出**: `select_parent(info)` at `src/namepredict/layer2/parent_selector.py:539-541`；候选列表由 `iter_parent_candidates(info)` at `parent_selector.py:532-536` 生成。

**构建**: 每个候选由 `_collect_candidates(info)` 生成，经 `_finalize_ranked(info, cands)` 排序并注入 `owned_atoms`（frozenset）。

### 通用字段

| 字段 | 类型 | 说明 |
|---|---|---|
| `kind` | `str` | 母体类型标识符，如 `"alcohol"`, `"acid"`, `"pyridine"`, `"benzene"`, `"alkane"` 等 |
| `chain` | `list[int]` | 有序的母体骨架原子索引（主链或环系统） |
| `n_carbons` | `int` | 母体中碳原子数量 |
| `owned_atoms` | `frozenset[int]` | 母体拥有的全部原子索引（由 `finalize_parent_ownership` 设置，不可变） |
| `stem_en` | `str` | 英文 stem 名称（如 `"ethanol"`, `"benzoic acid"`） |
| `stem_zh` | `str` | 中文 stem 名称（如 `"乙醇"`, `"苯甲酸"`） |
| `mol` | `rdkit.Chem.Mol` | RDKit Mol 引用（由 `pack_parent_stem` 附加） |

### FG 专属字段（按 kind 选择性存在）

以下字段根据母体类型（`kind`）选择性出现，并非所有 parent dict 都包含全部字段：

| 字段 | 类型 | 出现条件（kind） | 说明 |
|---|---|---|---|
| `cooh_c_idx` | `int` / `list[int]` | `acid`, `polycarboxylic` 等 | 羧基碳的原子索引 |
| `aldehyde_c_idx` | `int` | `aldehyde` 系 | 醛基碳的原子索引 |
| `ketone_c_idx` | `int` | `ketone` 系 | 酮羰基碳的原子索引 |
| `oh_c_idx` | `int` | `alcohol` 系 | 羟基所连碳的原子索引 |
| `amine_c_idx` | `int` | `amine` 系 | 胺基所连碳的原子索引 |
| `nitrile_c_idx` | `int` | `nitrile` 系 | 腈基碳的原子索引 |
| `ester_c_idx` | `int` | `ester` 系 | 酯基碳的原子索引 |
| `amide_c_idx` | `int` | `amide` 系 | 酰胺基碳的原子索引 |
| `sh_c_idx` | `int` | `thiol` 系 | 硫醇所连碳的原子索引 |
| `double_bond` | `tuple[int,int]` | `alkene` 系 | 双键原子对 |
| `triple_bond` | `tuple[int,int]` | `alkyne` 系 | 三键原子对 |
| `numbering_scaffold` | `dict` / `None` | fused ring 系 | 稠环编号骨架事实（来自 L2 scaffold specs） |
| `numbering_scaffold_required` | `bool` | fused ring 系 | 是否必须提供 numbering scaffold |
| `scaffold_id` | `str` | fused ring 系 | 骨架标识符（如 `"naphthalene"`） |

---

## 3. ClaimedBlock（L2 -> L3 桥接类型）

**全称**: Ownership-Only Claimable Side Block

**定义**: `src/namepredict/layer2/claimable_block.py:20-25`

**用途**: 描述 parent 未覆盖的一个"外侧"原子组件（side block），仅记录拓扑归属信息，不含名称。Layer 3 基于 ClaimedBlock 生成 `SubstituentName`。

```python
@dataclass(frozen=True)
class ClaimedBlock:
    slot: SideSlot          # 附着位置类型
    attach_parent: int      # 母体侧附着点 atom index
    root: int               # 取代基侧根原子 atom index
    atoms: frozenset[int]   # 该侧块的所有原子索引（不可变集合）
```

### SideSlot 枚举

定义于 `claimable_block.py:12-17`：

| 值 | 含义 |
|---|---|
| `CHAIN_C` | 附着于母体链上的碳原子 |
| `RING_C` | 附着于母体环上的碳原子 |
| `AMIDE_N` | 附着于酰胺氮（母体拥有羰基，侧块从氮延伸） |
| `ETHER_O` | 附着于醚氧 |
| `OTHER` | 其他类型（如杂原子附着） |

### 产出函数

`iter_claims(mol, owned_atoms)` at `claimable_block.py:152-159`：遍历所有外侧重原子组件，为每个组件确定 canonical edge `(attach_parent, root)` 和 slot，返回排序后的 `list[ClaimedBlock]`。

---

## 4. SubstituentName（L3 内部类型）

**定义**: `src/namepredict/layer3/substituent_namer.py:14-19`

**用途**: ClaimedBlock 经过命名后的结果，包含中英文名称和书写规则。

```python
class SubstituentName:
    claim: ClaimedBlock        # 来源声明块
    en: str                    # 英文取代基名称（如 "methyl", "chloro"）
    zh: str                    # 中文取代基名称（如 "甲基", "氯"）
    requires_parentheses: bool # 命名时是否需要括号（如复合取代基）
    backend: str               # 命名后端: "retained" / "rooted_tree" / "recursive"
```

**产出**: `SubstituentNamer.name(mol, claim)` at `substituent_namer.py:328-333`，依次尝试 retained -> rooted_tree -> recursive 三个后端，返回第一个成功的结果。

---

## 5. Substituent dict（L3 输出 -> L4 输入）

**全称**: Substituent Dictionary（取代基字典）

**产出**: `extract_substituents(info, parent)` at `src/namepredict/layer3/substituent_extractor.py:434-444`

**用途**: 每个 dict 描述一个待编号的取代基。Layer 4 通过 `_with_locants()` 向每个 substit dict 注入 `locant` 字段。

### 内置取代基字段（L3 产出）

| 字段 | 类型 | 说明 |
|---|---|---|
| `kind` | `str` | 取代基类型（如 `"alkyl"`, `"alkoxy"`, `"halo"`, `"hydroxy"` 等） |
| `en` | `str` | 英文名称（如 `"methyl"`, `"chloro"`, `"hydroxy"`） |
| `zh` | `str` | 中文名称（如 `"甲基"`, `"氯"`, `"羟基"`） |
| `attach_idx` | `int` | 附着于母体链上的 atom index（用于编号定位） |
| `atoms` | `list[int]` | 取代基占用的所有原子索引 |
| `paren` | `bool` | 命名时是否需要括号包裹 |
| `n_carbons` | `int` | 取代基中的碳原子数（影响排序优先级） |
| `backend` | `str` | 命名后端标识 |

### L4 注入字段

| 字段 | 类型 | 注入方 | 说明 |
|---|---|---|---|
| `locant` | `int` / `str` | `_with_locants()` at `numbering.py:472-473` | 该取代基在母体链上的位次编号 |

---

## 6. CoverageLedger（L3 -> Namer 桥接）

**全称**: Heavy-Atom Coverage Ledger（重原子覆盖台账）

**定义**: `src/namepredict/layer3/coverage.py:12-21`

**用途**: 验证 parent + named claims 是否完整覆盖分子的所有重原子（排除氢），无遗漏（gap）、无重叠（overlap）。

```python
@dataclass(frozen=True)
class CoverageLedger:
    owned_atoms: frozenset[int]                # parent 拥有的原子
    named_claims: tuple[SubstituentName, ...]  # 已命名的取代基列表
    gap: frozenset[int]                         # 未被任何方覆盖的重原子
    overlap: frozenset[int]                     # 被多方同时覆盖的重原子

    @property
    def complete(self) -> bool:
        return not self.gap and not self.overlap
```

**构建**: `build_coverage_ledger(mol, owned_atoms=owned, names=names)` at `coverage.py:51-66`

**使用**: `namer.py:_ledger_complete()` at `namer.py:63-65`。只有当 `complete == True` 时，命名结果才被视为有效（`coverage_complete: True` 写入 NameResult.meta）。

---

## 7. numbered dict（L4 输出 -> L5 输入）

**全称**: Numbered Structure Dictionary（编号结构字典）

**产出**: `number(parent, substituents)` at `src/namepredict/layer4/numbering.py:519-526`

**构建**: `_pack(oriented, subs_with_locants)` at `numbering.py:510-518`

### 结构

```python
{
    "parent":       {...},    # 增强的 parent dict（带 oriented chain + numbering plan）
    "substituents": [...],    # list[subst dict]，每项已注入 locant
    # 以下为展平到顶层的 locant 字段（来自 _fg_locants + _pack）:
    "oh_locant":              int | None,
    "oh_locants":             list[int] | None,
    "omit_oh_locant":         bool,
    "amine_locant":           int | None,
    "amine_locants":          list[int] | None,
    "omit_amine_locant":      bool,
    "sh_locant":              int | None,
    "omit_sh_locant":         bool,
    "ketone_locant":          int | None,
    "ketone_locants":         list[int] | None,
    "omit_ketone_locant":     bool,
    "cooh_locants":           list[int] | None,
    "ene_locant":             int | None,
    "ene_locants":            list[int] | None,
    "omit_ene_locant":        bool,
    "yne_locant":             int | None,
    "omit_yne_locant":        bool,
    # 来自 polycarboxylic / cyclo_relative_stereo:
    "stem_en":                str,   # (polycarboxylic 特定)
    "stem_zh":                str,   # (polycarboxylic 特定)
    "relative_stereo_prefix": str,
    # 由 namer._ok_result 注入:
    "name_mode":              str,   # "general" | "screened" | ...
}
```

### locant 字段详解

`_fg_locants()` at `numbering.py:500-509` 通过聚合三个子函数计算所有位次信息：

| 子函数 | 产出字段 |
|---|---|
| `_oh_am_locants()` at `numbering.py:483-491` | `oh_locant`, `oh_locants`, `omit_oh_locant`, `amine_locant`, `amine_locants`, `omit_amine_locant` |
| `_sh_locants()` at `numbering.py:492-494` | `sh_locant`, `omit_sh_locant` |
| `_unsat_locants()` at `numbering.py:474-482` | `ene_locant`, `ene_locants`, `omit_ene_locant`, `yne_locant`, `omit_yne_locant` |
| (inline in `_fg_locants`) | `ketone_locant`, `ketone_locants`, `omit_ketone_locant`, `cooh_locants` |

---

## 8. NameResult（L5 输出 / 最终 API 返回）

**全称**: Named Result Dataclass

**定义**: `src/namepredict/types.py:6-13`

**用途**: NamePredict 的最终对外返回值，包含中英文 IUPAC 名称及元信息。

```python
@dataclass
class NameResult:
    en: str                # 英文 IUPAC 名称
    zh: str                # 中文 IUPAC 名称
    success: bool          # 命名是否成功
    source: str = "iupac"  # 命名来源标识
    time_ms: float = 0.0   # 累计处理时间（毫秒）
    meta: dict = {}        # 元信息字典
```

### meta 字段

由 `_ok_result()` at `src/namepredict/namer.py:68-73` 填充：

| 字段 | 类型 | 来源 | 说明 |
|---|---|---|---|
| `parent_chain` | `list[int]` | `_chain_meta()` at `namer.py:31-33` | 母体链的原子索引列表 |
| `parent_kind` | `str` | `_chain_meta()` | 母体类型标识符 |
| `depth` | `int` | `_ok_result()` 参数 | 递归深度（一般化合物为 0） |
| `coverage_complete` | `bool` | 硬编码为 `True` | 仅当 CoverageLedger.complete 时调用 |
| `salt` | `str` | (如有盐) | 盐部分的名称 |
| `reason` | `str` | (如失败) | 失败原因（如 `"parse"`, `"no_candidate"`） |

### 实例化路径

- **成功路径**: `namer.py:68-70` → `_ok_result()` → `assemble(numbered, time_ms=...)` → `NameResult(en=..., zh=..., success=True, ...)`
- **失败路径**: `namer.py:20-21` → `_fail()` → `NameResult(en="", zh="", success=False, reason="...")`

---

## 数据流转汇总表

| 阶段 | 输入 | 产出 | 产出类型 | 关键函数 |
|---|---|---|---|---|
| L0 -> L1 | SMILES | `Mol` (organic) | `rdkit.Chem.Mol` | `preprocess_smiles()` |
| L1 | `Mol` | **info dict** | `dict` | `analyze(mol)` |
| L2 | info dict | **parent dict** | `dict` | `select_parent(info)` |
| L2 -> L3 | parent (owned_atoms) | **ClaimedBlock** | `dataclass` | `iter_claims(mol, owned_atoms)` |
| L3 | ClaimedBlock | **SubstituentName** | `dataclass` | `SubstituentNamer.name()` |
| L3 | info + parent | **subst dict** | `list[dict]` | `extract_substituents(info, parent)` |
| L3 -> Namer | owned + names | **CoverageLedger** | `dataclass` | `build_coverage_ledger()` |
| L4 | parent + subst dicts | **numbered dict** | `dict` | `number(parent, substituents)` |
| L5 | numbered dict | **NameResult** | `dataclass` | `assemble(numbered)` |

---

## 相关页面

- [[architecture/layer0-preprocessor]] — 预处理与盐拆分
- [[architecture/layer1-analyzer]] — FG 分析器详解
- [[architecture/layer2-parent-selector]] — 母体选择器详解（核心）
- [[architecture/layer3-substituents]] — 取代基提取与命名
- [[architecture/layer4-numbering]] — 编号与定位符分配
- [[architecture/layer5-name-assembly]] — 中英双语名称组装
- [[concepts/parent-selection]] — 母体选择规则与优先级
- [[concepts/functional-groups]] — 官能团分类与优先级表
- [[concepts/iupac-rules]] — IUPAC 蓝皮书规则映射
- [[concepts/bilingual-naming]] — 中英双语命名约定
- [[index]] — Wiki 首页
