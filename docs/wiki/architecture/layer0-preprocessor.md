# Layer0: 预处理 (Preprocessor)

> **管线位置:** 第 0 层 / 6 层 | **源文件:** 5 个 `.py` (335 行) | **最后更新:** 2026-09-10

---

## 概述

Layer0 是 NamePredict 6 层命名管线的入口层，负责将外部输入的 SMILES 字符串转换为内部可操作的分子表示（RDKit `Mol` 对象），并执行盐（salt）检测与解离。`preprocess` 除解析与消毒外，还会做**立体初步指派**、**酰胺烯醇互变异构归一化**（把非芳香中性的 `C(OH)=N` 位点归一为酮式 `C(=O)-NH`，见 §1.1）与**酸性质子重定位**（把质子化强酸 OH 与去质子化弱酸阴离子共存时的负电荷收敛到最强酸位，见 §1.2）。它是整个管线中唯一与"原始文本"打交道的层，也是错误处理的第一道防线。

**管线中的位置：**

```
SMILES 字符串
    ↓
┌──────────────┐
│  Layer0 预处理 │  ← 当前层
└──────────────┘
    ↓  RDKit Mol
┌──────────────┐
│  Layer1 分析  │
└──────────────┘
    ↓  ... → Layer5 组装 → NameResult
```

**职责边界：**

| 职责 | 说明 |
|------|------|
| SMILES 解析与归一化 | 将文本输入转换为 `rdkit.Chem.Mol` 对象：空白校验 → `MolFromSmiles(sanitize=False)` → `SanitizeMol` → `AssignStereochemistry`（立体初步指派）→ `normalize_amide_tautomer`（酰胺烯醇 → 酮式）→ `normalize_acid_charge`（酸性质子重定位） |
| 输入验证 | 检测空字符串、纯空白、无效 SMILES 语法，及消毒/归一化异常 |
| 盐解离 | 识别碱性金属盐（Li/Na/K）和盐酸盐（HCl），分离有机片段 |
| 盐元数据生成 | 生成中英双语盐名称元数据，供 Layer5 拼接 |

**非职责（由下游层承担）：**

- 不执行任何官能团识别（属 Layer1）
- 不进行母体选择或命名决策（属 Layer2）
- 不对有机片段本身做任何化学分析

**输入/输出类型：**

```
preprocess(smiles: str) -> Mol | None
dissociate_salt(mol: Mol) -> tuple[Mol, dict]
```

`preprocess` 返回 `None` 表示输入无法解析，管线在此终止并返回失败 `NameResult`。`dissociate_salt` 的 `dict` 为空 `{}` 表示分子不是可识别的简单盐，管线继续以原始分子运行；非空时包含 `salt` 元数据，合并到最终 `NameResult.meta` 中。

---

## 核心逻辑

### 1. SMILES 解析 (`preprocess`)

`preprocess` 函数是 NamePredict 管线入口中最关键的函数之一，负责把原始 SMILES 规整为可供下游各层直接消费的分子形态。其内部为五步：

1. **空白校验**：检查输入是否为 `None`、空字符串或纯空白字符串。任何无法形成有效文本的输入均立即返回 `None`，避免将空值传递给 RDKit 造成异常。
2. **解析与消毒**：调用 `Chem.MolFromSmiles(str(smiles).strip(), sanitize=False)` 先做语法解析，成功后再以 `Chem.SanitizeMol(mol)` 完成消毒。RDKit 的 `MolFromSmiles` 在遇到无效 SMILES 语法时返回 `None`（而非抛出异常），这与 NamePredict 的错误处理策略一致——所有解析失败都统一为"静默失败"，由上层 `_pipeline` 函数捕获并生成 `success=False` 的 `NameResult`，其中 `meta={"reason": "parse"}` 标记失败原因。
3. **立体初步指派**：调用 `Chem.AssignStereochemistry(mol, force=True, cleanIt=False, flagPossibleStereoCenters=True)`，在进入 L1 前即对手性中心 / 双键预先指派 R/S、E/Z，供 L4 位次与 L5 立体描述符使用。
4. **酰胺烯醇互变异构归一化**：调用 `normalize_amide_tautomer(mol)`（`tautomer.py`），把分子内非芳香中性的 `C(OH)=N` 烯醇位点归一化为酮式 `C(=O)-NH`（见下 §1.1），使下游面对统一的酰胺形态。
5. **酸性质子重定位**：调用 `normalize_acid_charge(mol)`（`charge.py`），把同片段内质子化强酸（默认羧酸 `C(=O)OH`）与去质子化弱酸位（酚氧/烯醇氧/酰胺氧、去质子化氮）共存时的质子搬到弱酸位，使等价负电荷收敛到最强酸（见下 §1.2）——让既有 carboxylate→-oate 的 L1–L5 能力作用于输入侧电荷错位的结构（如 chebi-433）。

以上任一步抛出异常（如消毒失败）都会返回 `None`。

**错误分类**：Layer0 能检测两类输入错误：
- **空白输入**：`""`、`None`、`"   "` -- 在进入 RDKit 之前即被拦截
- **语法无效**：如 `"CCOO"` 中碳的五价——由 RDKit `MolFromSmiles` 返回 `None` 拦截

> **源:** `src/namepredict/layer0/preprocessor.py:11-25`

#### 1.1 酰胺烯醇互变异构归一化 (`tautomer.py`)

模块 `src/namepredict/layer0/tautomer.py`（65 行）专门做互变异构规范化：

- 位点判定：`_is_amide_enol_o`（`tautomer.py:13`）要求与碳**单键**相连、中性、仅带隐氢的羟基 O；`_is_amide_enol_n`（`tautomer.py:23`）要求与碳**双键**相连、非芳香、中性、无显式 H 的亚胺 N。遍历排除芳香 C 后由 `_amide_enol_sites(mol)`（`tautomer.py:33`）汇总全部 `(c_idx, n_idx, o_idx)` 位点。
- 改写方式：`normalize_amide_tautomer(mol)`（`tautomer.py:48`）对每个位点把 C=N 双键降为单键、C–O 单键升为 C=O 双键，**只改键级**并靠 RDKit 隐氢重算完成质子迁移——不增删重原子、不改原子序；改写后再次 `SanitizeMol` + `AssignStereochemistry`。
- 保守跳过：带电 N、O⁻ 阴离子、显式 `[H]`、硫类似物（C=S 等）位点一律不处理，不做质子化 / 阴离子改写。
- 返回约定：无位点，或改写后消毒失败时，原样返回输入的 `mol`。

#### 1.2 酸性质子重定位 (`charge.py`)

模块 `src/namepredict/layer0/charge.py`（129 行）专门处理输入侧**电荷错位**结构——同一片段内去质子化弱酸位（酚氧/烯醇氧/酰胺 O⁻、去质子化 N⁻）与质子化强酸位（默认羧酸 `C(=O)OH`）共存时，把质子从强酸搬到弱酸位：等价负电荷收敛到最强酸，使既有 carboxylate→-oate 的 L1–L5 命名能力作用于此类输入（如 chebi-433）。**只改 FormalCharge 与 H 记账，不增删重原子**，走 RWMol 原子级编辑，改后 `SanitizeMol` + `AssignStereochemistry`；含 `*` dummy 或改动失败时保守跳过。

- 强弱判定：`_acid_kind_of_oh`（`charge.py:36`）识别中性含 H 的 O 是否质子化强酸 OH（carboxyl：O–H 连 `C(=O)`；phospho：连 `P(=O)`；sulfo：连 `S(=O)n`；醇/酚/酯 O 无成酸中心天然返回 None）；`_is_weak_anion`（`charge.py:62`）识别去质子化弱酸位（-1 电荷、O/N、非强酸共轭碱、邻接无 +1 内平衡写法）。
- 方向性保证：供体只取 `_DONOR_KIND` 中酸类（`charge.py:14`，默认只开 `"carboxyl"`，磷酸/磺酸 donor 经实测不贡献修复只扩大 blast radius）、受体只取弱酸阴离子 → 单调收敛，每次搬走一对后不再进入候选，终止于无配对。搬运序以酸 kind 优先级（`_KIND_PRIO` `charge.py:15`）加 canonical rank 最小保证确定性。
- 搬运方式：`_relocate_proton`（`charge.py:73`）把强酸 OH 质子搬到弱酸受体——受体变中性 +1 显式 H，供体变 -1 减 1 H；**原子级记账而非重写 SMILES**，保立体中心不随邻居重排翻转；`SanitizeMol` 失败返回 `None` 保守跳过。
- 返回约定：无改动返回原 `mol` 对象；同片段成对质子全部搬运（上界 `mol.GetNumAtoms()` 次循环，每次消耗一对 donor/acceptor）。

### 2. 盐检测与解离 (`dissociate_salt`)

`salt.py` 实现了对简单盐结构的检测和有机片段的分离，是 Layer0 中逻辑最复杂的模块。其设计遵循 IUPAC 功能类命名（functional class nomenclature）中对盐的处理规则：将盐视为"阳离子（cation）修饰有机母体名称"的模式，例如 "sodium benzoate"（苯甲酸钠）。

#### 2.1 碎片分割

函数入口 `dissociate_salt(mol)`（`salt.py:95`）先做**计数早退**：调用不带 `asMols` 的 `Chem.GetMolFrags(mol)`（`salt.py:104`），它只返回各碎片的原子索引元组，不重建分子；片段数 < 2 时直接返回 `(mol, {})`。这一步是性能关键——`asMols=True, sanitizeFrags=True` 会为每个片段重建 `Mol` 并重新 sanitize（实测 0.003ms vs 0.103ms），而 `_name_mol` 每次递归都会调用本函数，绝大多数分子只有 1 个片段。

片段数 ≥ 2 时才真正分割：`Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)` 返回可独立操作的碎片分子，`sanitizeFrags=True` 确保每个碎片都是有效的化学子结构；分割后若片段少于 2 个同样返回原分子与空元数据。

#### 2.2 碱性金属盐检测 (`_alkali_en`)

支持的碱性金属阳离子限于三种：**锂 (Li, Z=3)**、**钠 (Na, Z=11)**、**钾 (K, Z=19)**。检测条件严格：
- 碎片必须为单原子（`GetNumAtoms() == 1`）
- 该原子的形式电荷（formal charge）必须为 `+1`

这排除了中性金属原子（电荷 0）和多原子阳离子（如铵根 NH4+）。注意：钙 Ca²⁺、镁 Mg²⁺ 等多价阳离子不在当前支持范围内，这是有意为之的设计简化——多价盐需要更复杂的化学计量处理（如 "calcium dibenzoate" vs "sodium benzoate"）。

> **源:** `src/namepredict/layer0/salt.py:15-23`（`_alkali_en`；金属名表 `_ALKALI_EN`:11 / `_METAL_ZH`:12）

#### 2.3 盐酸盐检测 (`_is_hcl_frag`)

盐酸盐（hydrochloride）的检测覆盖两种 RDKit 可能产生的表示形式：
- **游离氯离子** `[Cl-]`：原子序数 17，形式电荷 -1，无氢原子
- **中性 HCl 分子**：原子序数 17，形式电荷 0，带 1 个氢原子

两种表示均被识别为同一类反离子。这一设计处理了 SMILES 输入的不同写法（如 `CCO.Cl` 产生 HCl 碎片、`CCO.[Cl-]` 产生 Cl⁻ 碎片），确保两种写法都能正确生成 "hydrochloride / 盐酸盐" 的元数据。

> **源:** `src/namepredict/layer0/salt.py:33-43`

#### 2.4 水分子的静默忽略 (`_is_water`)

在盐的 SMILES 表示中，水合水（water of hydration）以独立碎片形式出现（如 `O` 表示 H₂O）。`_is_water` 检测单原子氧（Z=8）、电荷 0、带 2 个氢原子的碎片，并在碎片分类阶段将其静默丢弃——不计入金属列表、不计入有机碎片列表、不干扰 HCl 计数。这意味着 NamePredict 自动忽略水合水，不对其命名（如不输出 "monohydrate"）。

#### 2.5 碎片分类与仲裁 (`_partition` + `_from_frags`)

碎片分类由 `_bucket_frag` 完成，将每个碎片分入三类之一：金属阳离子列表、有机碎片列表、或 HCl 计数。`_partition` 累加所有碎片后，`_from_frags` 执行仲裁逻辑：

1. **单有机碎片约束**：仅当 `len(organics) == 1` 时才进行盐识别。多个有机碎片的情况（如混合物、共晶、多组分盐）返回空元数据。
2. **优先级规则**：碱性金属盐优先于盐酸盐。如果同时检测到金属阳离子和 HCl，当前实现拒绝识别（返回 `None`），这是对复杂混合盐的保守处理。
3. **同种金属约束**：多个金属阳离子碎片必须为同一种金属（`len(set(metals)) == 1`），否则拒绝识别。这处理了混合碱盐（如 LiNa 混合盐）的罕见情况。
4. **单 HCl 约束**：仅支持恰好 1 个 HCl 反离子。多个 HCl 的情况（如二盐酸盐 `dihydrochloride`）不在当前支持范围内。

> **源:** `src/namepredict/layer0/salt.py:58-91`（`_partition` / `_from_frags`）

#### 2.6 盐元数据结构

成功识别盐后，返回的元数据字典包含中英双语字段，供 Layer5 在名称组装时使用：

| 盐类型 | 元数据键 | 示例值 |
|--------|---------|--------|
| 碱性金属盐 | `metal` | `"sodium"` |
| 碱性金属盐 | `metal_zh` | `"钠"` |
| 碱性金属盐 | `n_metal` | `1` |
| 盐酸盐 | `acid_salt` | `"hydrochloride"` |
| 盐酸盐 | `acid_salt_zh` | `"盐酸盐"` |

这些元数据在 `_name_mol` 中通过 `result.meta["salt"] = salt` 注入到最终的 `NameResult`，供 Layer5 组装器在生成中英双语名称时追加盐后缀。

### 3. 管线集成

Layer0 在 `namer.py` 的 `_pipeline` 函数中被调用（`namer.py:289`），是管线的第一个处理步骤。调用流程：

```
SMILESNNamer.name(smiles)              # namer.py:343
  → memo.begin_run()                    # namer.py:349 清空本次命名的中间结果记忆
  → _pipeline(smiles, t0)               # namer.py:289
      → preprocess(smiles)              # namer.py:291 → layer0/preprocessor.py:11
      → if None: _fail("parse")         # 解析失败终止
      → _name_mol(mol, ...)             # namer.py:259
          → dissociate_salt(mol)        # namer.py:272 → layer0/salt.py:95
          → analyze(organic)            # 进入 Layer1
          → _run_candidates(...)        # Layer2-5
          → if salt: result.meta["salt"] = salt  # 注入盐元数据
```

当 `preprocess` 返回 `None` 时，`_pipeline` 通过 `_fail` 生成 `NameResult(en="", zh="", success=False, meta={"reason": "parse"})` 并直接返回，不再进入后续层。这是 NamePredict 的快速失败（fail-fast）策略——在管线最前端拦截无效输入，避免下游层对空对象进行无效计算。

盐解离发生在 `_name_mol` 中（`namer.py:272`），在 Layer1 分析之前。这意味着 Layer1-5 始终处理的是解离后的纯有机片段，保证了各层代码无需关心盐的存在，实现了关注点分离。

`SMILESNNamer.name` 在进入 `_pipeline` 前调用 `memo.begin_run()`（`namer.py:349`）清空本次命名的**中间结果记忆**（`src/namepredict/tools/memo.py`）：它按分子对象记忆同一次命名内恒定、会被各层反复计算的量（环感知、CIP 标签、完全氢化骨架、锚定子分子等），纯消除重复计算、不改变任何返回值，跨分子不共享。同一机制在 L1 有调用点（见 [[architecture/layer1-analyzer]] §7.1）。另有 `src/namepredict/tools/rdkit_fast.py` 在包导入时由 `namepredict/__init__.py:5` 安装补丁，把 `Chem.Mol.GetAtoms/GetBonds` 换成索引循环、去掉 RDKit 生成器的每项包装开销。

> **源:** `src/namepredict/namer.py:289` / `:349`

---

## 文件清单

| 文件 | 行数 | 说明 |
|------|------|------|
| `src/namepredict/layer0/__init__.py` | 7 | 包入口，导出 `preprocess` 和 `dissociate_salt` 两个公共接口 |
| `src/namepredict/layer0/preprocessor.py` | 25 | SMILES 预处理：`preprocess(smiles) -> Mol | None`，空白校验 + 解析消毒 + 立体初步指派 + 酰胺烯醇归一化 + 酸性质子重定位 |
| `src/namepredict/layer0/tautomer.py` | 65 | 酰胺烯醇互变异构归一化：非芳香中性 `C(OH)=N` → `C(=O)-NH`，`normalize_amide_tautomer` 供 preprocess 复用 |
| `src/namepredict/layer0/charge.py` | 129 | 酸性质子重定位：同片段质子化强酸 + 去质子化弱酸位共存时收敛负电荷到最强酸，`normalize_acid_charge` 供 preprocess 复用 |
| `src/namepredict/layer0/salt.py` | 109 | 盐解离引擎：检测 Li/Na/K 金属盐和 HCl 盐酸盐，返回有机片段与双语元数据；`dissociate_salt`:95 先以 `GetMolFrags(mol)` 计数早退（>:104）再按需重建碎片 |

---

## 数据流图

### 主流程

```mermaid
flowchart TD
    A["SMILES 字符串输入"] --> B{"preprocess()"}
    B -->|"空/无效 SMILES"| C["返回 None"]
    C --> D["_pipeline: _fail('parse')"]
    D --> E["NameResult(success=false)"]

    B -->|"有效 SMILES"| P2["MolFromSmiles(sanitize=False)<br/>+ SanitizeMol<br/>+ AssignStereochemistry(立体初步指派)"]
    P2 --> P3["normalize_amide_tautomer<br/>(C(OH)=N → C(=O)-NH)"]
    P3 --> P3b["normalize_acid_charge<br/>(负电荷收敛到最强酸)"]
    P3b --> F["RDKit Mol 对象(酮式酰胺、电荷收敛)"]
    F --> G{"dissociate_salt()"}
    G --> H{"GetMolFrags(mol)<br/>碎片数 ≥ 2?<br/>(计数早退, 不重建碎片)"}

    H -->|"否 (单一分子)"| I["返回 (原始 mol, {})"]
    I --> J["进入 Layer1 分析"]

    H -->|"是"| K{"碎片分类\n_bucket_frag × N"}

    K --> L["碱性金属?\nLi+/Na+/K+"]
    K --> M["HCl/Cl-?"]
    K --> N["H₂O?"]
    K --> O["有机片段"]

    L -->|"是"| P["加入 metals 列表"]
    M -->|"是"| Q["递增 n_hcl 计数"]
    N -->|"是"| R["静默丢弃"]
    O --> R2["加入 organics 列表"]

    P --> S{"_from_frags 仲裁"}
    Q --> S
    R2 --> S

    S -->|"单有机 + 同种金属"| T["生成金属盐元数据\n{metal, metal_zh}"]
    S -->|"单有机 + 1 HCl"| U["生成盐酸盐元数据\n{acid_salt, acid_salt_zh}"]
    S -->|"不符合规则"| V["返回 (原始 mol, {})"]

    T --> W["返回 (有机 mol, salt_meta)"]
    U --> W
    W --> J
    V --> J

    style A fill:#e1f5fe
    style E fill:#ffcdd2
    style J fill:#c8e6c9
    style C fill:#ffcdd2
```

### 管线集成上下文

```mermaid
flowchart LR
    subgraph Layer0
        P["preprocess()"] --> S["dissociate_salt()"]
    end

    subgraph 管线
        SMILES["SMILES"] --> P
        P -->|Mol| S
        S -->|"organic Mol + salt_meta"| L1["Layer1\nanalyze()"]
        L1 --> L2["Layer2\nparent select"]
        L2 --> L3["Layer3\nsubstituents"]
        L3 --> L4["Layer4\nnumbering"]
        L4 --> L5["Layer5\nassemble()"]
    end

    subgraph 输出
        L5 --> NR["NameResult"]
        S -.->|"salt_meta 注入"| NR
    end

    style Layer0 fill:#fff3e0,stroke:#ff9800
    style NR fill:#c8e6c9
```

---

## 对外接口

### `preprocess(smiles: str) -> Mol | None`

SMILES 字符串到 RDKit 分子对象的转换函数；返回的 `Mol` 已完成消毒、立体初步指派、酰胺烯醇互变异构归一化与酸性质子重定位。返回 `None` 表示输入无效或归一化失败。

| 参数 | 类型 | 说明 |
|------|------|------|
| `smiles` | `str` | 输入的 SMILES 字符串 |

| 返回值 | 说明 |
|--------|------|
| `rdkit.Chem.Mol` | 消毒后、立体已指派、完成酰胺烯醇归一化与酸性质子收敛的分子对象 |
| `None` | 空输入、无效 SMILES，或解析/消毒/归一化抛异常 |

**调用者:** 仅 `namer.py:_pipeline` 直接调用。

### `dissociate_salt(mol: Mol) -> tuple[Mol, dict]`

检测并解离简单盐结构，返回有机片段与盐元数据。

| 参数 | 类型 | 说明 |
|------|------|------|
| `mol` | `rdkit.Chem.Mol` | 输入分子（可能含盐碎片） |

| 返回值 | 说明 |
|--------|------|
| `tuple[Mol, dict]` | `(organic_mol, salt_meta)`；未经识别的盐返回 `(原始mol, {})` |

**`salt_meta` 字典字段（当非空时）：**

| 键 | 类型 | 出现条件 | 说明 |
|----|------|---------|------|
| `metal` | `str` | 碱性金属盐 | 英文金属名 (`"sodium"`) |
| `metal_zh` | `str` | 碱性金属盐 | 中文金属名 (`"钠"`) |
| `n_metal` | `int` | 碱性金属盐 | 金属阳离子数量 |
| `acid_salt` | `str` | 盐酸盐 | 英文酸盐名 (`"hydrochloride"`) |
| `acid_salt_zh` | `str` | 盐酸盐 | 中文酸盐名 (`"盐酸盐"`) |

**调用者:** `namer.py:_name_mol`。

---

## 相关页面

- [[architecture/overview]] — 6 层架构总览与层间数据流
- [[architecture/layer1-analyzer]] — 下一层：官能团分析器（接收 Layer0 输出的 Mol）
- [[architecture/layer5-name-assembly]] — 盐元数据的最终消费方：中英双语名称拼接
- [[reference/core-data-contracts]] — `NameResult` 及 `meta` 字段结构
- [[index]] — Wiki 首页
