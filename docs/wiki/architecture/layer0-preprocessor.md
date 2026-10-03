# Layer0: 预处理 (Preprocessor)

> **源文件:** `layer0/` 5 个 `.py` (549 行) | **对外接口:** `preprocess` / `dissociate_salt` / `pair_unique_ions`

## 概述

本层把外部 SMILES 字符串转成内部可用的 RDKit `Mol`，并完成同位素采集、互变归一、酸性质子重定位与盐解离。层内共享常量取自 `constants.py`。

```
SMILES 字符串
   ↓ Layer0 预处理（当前层）
RDKit Mol
   ↓ Layer1 分析 → … → Layer5 组装
NameResult
```

| 接口 | 签名 | 说明 |
|---|---|---|
| `preprocessor.preprocess` | `(smiles: str) -> Mol \| None` | 清洗解析 + 同位素采集 + 归一化；空输入或异常返回 `None` |
| `salt.dissociate_salt` | `(mol: Mol) -> tuple[Mol, dict]` | 分离有机片段与盐元数据；`dict` 为空表示非简单盐 |
| `salt.pair_unique_ions` | `(frags) -> tuple[Mol, Mol, int, int] \| None` | 唯一可配对的有机阴阳离子及份数，供多片段 P-77 二元盐名 |

非职责：不识别官能团（L1）、不选母体、不做位次与名称组装（L2–L5）。

本层是全管线唯一接触原始文本之处，也是错误处理的第一道防线：空串、纯空白、非法 SMILES 语法、消毒或归一化抛出的异常，一律收敛为 `None`，由上层转成解析失败结果。

## 核心逻辑

### 同位素采集（`preprocessor.py`）

`_capture_isotopes` 在清零同位素前把信息记到重原子属性上，供 L5 渲染同位素名：

- 非氢核素：属性 `ISO_NUCLIDE_PROP`（`isoNuclide`）写在核素原子自身，值形如 `13C`/`18F`/`123I`。
- 氘/氚：按质量数经 `ISO_H_PROPS` 映射为计数属性 `iso2H`/`iso3H`，累加到其唯一重原子邻居（计数为 0 时不写属性）。

`_strip_isotopes` 随后对每个原子 `SetIsotope(0)` 清零——命名管线用重原子属性渲染同位素名，同位素标记留在原子上会让锚定键匹配失败。

### 互变异构归一化（`tautomer.py`）

入口 `normalize_amide_tautomer(mol)` 是双向归一化器：把烯醇/烯硫醇式拉回酮式、环外亚胺式拉回胺式，并把酰胺肟拉回羟基亚胺式。依次取三类位点——`_enol_candidates`（烯醇）→ `_amidine_candidates`（亚胺）→ `_amidoxime_candidates`（酰胺肟），后两类排除与已选烯醇位点同中心碳者；三类皆空直接返回原 mol。有位点在 `RWMol` 上就地改键后 `SanitizeMol`，失败整体回退，成功再 `AssignStereochemistry`。可调判据为文件顶部的 `ENOL_X`(O/S) 与 `N_VALENCE`(3)。

- 烯醇侧 `_enol_candidates`：中心碳任取；可迁移杂原子由 `_is_amide_enol_x` 判——`ENOL_X`、中性、与碳单键、度为 1、带 ≥1 H（隐氢与显式 H 记账都算）。受体氮由 `_is_amide_enol_n` 判——中性、无显式 H、度+氢数未达 `N_VALENCE`、与碳成双键或芳香键。
- 亚胺侧 `_amidine_candidates`：中心碳在环内且带环外亚胺氮（`_is_imine_n`：中性、非碳邻居一律排除、可再容纳 1 个 H、不在环内），受体为环内氨基氮（`_is_amidine_n`：在环内、带 H、单键或芳香键连碳）。
- 酰胺肟侧 `_amidoxime_candidates`（P-66.4.4）：中心碳恰连两个氮、无其它杂原子，其一为环外亚胺氮、另一为 `_is_hydroxyamino_n` 判定的羟基氨基氮（中性、非环、带 1 H、单键连碳、另连锁 `ENOL_X` 且带 ≥1 H）。
- `_assign` 按候选受体氮数由少到多贪心配对，一个氮只服务一个位点。
- `_normalize_amide` / `_normalize_amidine` 就地搬迁 H 与键级；`_normalize_amide` 跳过羟基氢数不为 1 的位点，`_normalize_amidine` 同时承担亚胺与酰胺肟两类位点。显式 H 记账必须保留——芳环 N 的隐氢受芳香性约束不会自动补，只能 `SetNumExplicitHs` 显式加，否则 Kekulé 奇偶不匹配。

### 预处理链（`preprocess`）

空白校验 → `MolFromSmiles(sanitize=False)` → `_strip_isotopes`（内含 `_capture_isotopes`）→ `SanitizeMol` → `AssignStereochemistry` → 互变归一 → 酸性质子重定位 → 二次互变归一。第二次是因为重定位把酰胺 O⁻ 变成中性 `C(OH)=N`，此时才出现新的烯醇位。

### 电荷归一（`charge.py`）

`normalize_acid_charge` 把质子从强酸根供体搬到受体（`_relocate_proton`，只做电荷/H 记账）；消毒失败、SMILES 不变或无同片段受体即停，含通配原子（`*`）时返回原 mol。

- 供体：`_acid_kind_of_oh` + `ACID_CENTERS` 判——中性 OH 所连非芳香中心带 ≥1 双键氧，酸类须在 `DONOR_KIND`（`carboxyl`/`phospho`/`sulfo`）内，表键序即酸强度序。
- 受体分两类：`_is_weak_anion` 判定的去质子化弱酸位（元素限 `ACCEPTOR_Z`，现仅 `O`）为首选；仅当整分子净电荷为负且无弱酸位时，才把 `_is_base_n` 判定的中性脂肪胺氮当受体（保净电荷内的两性离子式，P-73.1.2）。
- `_is_weak_anion` 的排除项：邻接 +1 电荷者（硝基/N-氧化物等内平衡写法）、已是强酸共轭碱的 O（羧酸/磷酸/磺酸根）、以及环系邻域带强酸的芳环氧负离子。后者经 `_ring_system_of`（并含共享原子的稠环系）→ `_ring_bears_strong_acid`（环原子及其 2 键邻域内有中性 `-COOH`/`-SO3H`/`-PO(OH)2`）判定：有此酸时保留酚 O⁻ 作 `-olate` 母体（酸降为 `carboxy`/`sulfo` 前缀）；但 `_ring_has_nh`（环系内有带 H 的环氮）为真时不保留，按内酰胺/酮式互变优先搬质子（金标 4-oxo-1H 式）。
- 确定性与循环：供体按 `DONOR_KIND` 序优先、再取 canonical rank 最小排序，受体取同片段内 rank 最小者；循环逐次消耗一对供体/受体，上界为原子数。

`N⁻` 不在 `ACCEPTOR_Z` 内，按 P-72.2.2.2 保留为 `azanide` 母体，不作质子受体。

### 盐解离（`salt.py`）

`dissociate_salt` 先按片段数过滤（单片段直接视为无盐），再由 `_from_frags` 识别反离子。`_frag_role` 对单片段给出三种角色及其中英盐名，非盐离子与有机片段返回 `None`：

| 角色 | 判据 | 名称来源 |
|---|---|---|
| `metal` | 规范 SMILES 命中 `POLY_CATION`；或单原子、电荷 ≥ +1、元素在 `METAL_ION_EN` | `POLY_CATION`；`METAL_ION_EN`/`METAL_ION_ZH`（P-71.2） |
| `halide` | 规范 SMILES 命中 `POLY_ANION`；或单原子、元素在 `HALO_Z`、电荷 -1 且无 H | `POLY_ANION`；`HALIDE_EN`/`HALIDE_ZH` |
| `hx` | 单原子、元素在 `HALO_Z`、电荷 0 且带 1 H | `HALIDE_HX_EN`/`HALIDE_HX_ZH`（P-71.3） |

`POLY_CATION` = `{[NH4+] → azanium/铵}`；`POLY_ANION` = `{[OH-] → hydroxide/氢氧化物, O=[N+]([O-])[O-] → nitrate/硝酸盐, [O-][Cl+3]([O-])([O-])[O-] → perchlorate/高氯酸盐, [N-]=[N+]=[N-] → azide/叠氮化物, F[P-](F)(F)(F)(F)F → hexafluorophosphate/六氟磷酸盐}`，键为片段规范 SMILES。多原子阳离子并入 `metal` 通道、多原子阴离子并入 `halide` 通道。

三种反离子各有额外的有机侧约束，不满足时返回 `(mol, {})`：

| 反离子 | 有机侧约束 | 元数据键 |
|---|---|---|
| `metal` | 金属种类唯一、不与阴离子/HX 共存、有机片段为阴离子（净电荷 < 0） | `metal`/`metal_zh`/`n_metal`/`n_org` |
| `halide` | 种类唯一、不与 HX 共存；有机物可为阳离子、中性碱或阴离子（P-71.2）；单一 `Cl⁻` + 净电荷 0 且无内盐的有机碱改按盐酸盐命名（P-71.3） | `halide`/`halide_zh`/`n_org`（Cl⁻ 特例则 `acid_salt`/`acid_salt_zh`/`n_org`） |
| `hx` | 种类唯一、有机片段中性 | `acid_salt`/`acid_salt_zh`/`n_org` |

有机片段须至少一个且彼此一致（按非立体 SMILES 去重），份数记入 `n_org`。份数 >1 时经 `_mult_word` 加 `MULT_EN`/`MULT_ZH` 数量前缀（两表现覆盖 1–9999）。

`pair_unique_ions` 面向多片段组分名：净电荷 > 0 的片段归阳离子、< 0 的归阴离子，要求两侧物种各唯一且整体电荷守恒，返回 `(阳离子片段, 阴离子片段, n_阳, n_阴)`；否则 `None`。`namer._name_components` 据此拼 P-77.1.1 二元盐名，`namer._name_component` 也用 `_frag_role` 直接为单个反离子片段取名。

### 简单分子查表（`constants.SIMPLE_MOLECULES`）

单元素、单质与极简单分子（重原子 ≤ `SIMPLE_MAX_HEAVY`=3）无母体链/环，走不到 L2–L5，故在命名入口由 `namer._simple_species` 查表取 IUPAC 保留名。键为 L0 电荷归一后的规范 SMILES，值为 `(en, zh)`，覆盖单元素分子（`[H][H]`/`O=O`/`O=[O+][O-]`）、母体氢化物（`O`/`N`/`P`/`S`/`Si` 等第 15–17 族）、简单氧化物/硫化物/氮化物（`NO`/`C#O`/`O=C=O`/`O=S=O` 等）、简单卤素化合物与单质金属（由 `ELEMENT_METAL_NAMES` 生成 `[X]` 键，P-21 保留名）。

## 与其他层的契约

- 输出：消毒并完成互变/电荷归一的 `Mol`；有机片段原子索引不变，同位素信息以原子属性保留——非氢核素 `ISO_NUCLIDE_PROP`（字符串），氘/氚 `ISO_H_PROPS` 派生的计数属性——供 L5 渲染同位素名。
- 盐元数据含中英双语盐名，键随反离子类型而定：金属/多原子阳离子盐为 `metal`/`metal_zh`/`n_metal`/`n_org`，卤素/多原子阴离子盐为 `halide`/`halide_zh`/`n_org`，氢卤酸（含 `Cl⁻`+中性碱特例）为 `acid_salt`/`acid_salt_zh`/`n_org`。
- 消费方：`namer` 经 `preprocess` 取分子、经 `dissociate_salt` 分离有机片段，并把盐元数据注入 L1 的 `info["salt"]`，供 L2 盐门控（写 `salt_meta`）与 L5 拼盐名后缀；`namer._apply_salt_suffix` 经 `stems.join_metal_salt_names` 拼后缀，含氧酸中心母体自带盐组装故跳过。多片段且无简单盐时，`namer._name_components` 经 `pair_unique_ions` 走 P-77 二元盐名，其余组分按字母序以 `; ` 连接。
- 两个接口相互独立：`preprocess` 只吃 SMILES 做解析与归一化，`dissociate_salt` 只吃已消毒的 `Mol` 做片段拆分；调用方先预处理再解离。盐解离不改变有机片段的原子索引，故 L1 之后的层无需感知盐的存在。
- 本层不保留任何输入文本，下游一律基于 `Mol` 与原子索引工作。盐数据契约详见 [[reference/core-data-contracts]]。

相关页面：[[architecture/layer1-analyzer]]、[[architecture/layer5-name-assembly]]、[[reference/core-data-contracts]]。
