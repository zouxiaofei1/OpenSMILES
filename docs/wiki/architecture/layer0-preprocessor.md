# Layer0: 预处理 (Preprocessor)

> **源文件:** `layer0/` 5 个 `.py` (398 行) | **对外接口:** `preprocess` / `dissociate_salt`

## 概述

本层把外部 SMILES 字符串转成内部可用的 RDKit `Mol`，并完成互变归一、酸性质子重定位与盐解离。层内共享常量取自 `constants.py`。

```
SMILES 字符串
   ↓ Layer0 预处理（当前层）
RDKit Mol
   ↓ Layer1 分析 → … → Layer5 组装
NameResult
```

| 接口 | 签名 | 说明 |
|---|---|---|
| `preprocessor.preprocess` | `(smiles: str) -> Mol \| None` | 清洗解析 + 归一化；空输入或异常返回 `None` |
| `salt.dissociate_salt` | `(mol: Mol) -> tuple[Mol, dict]` | 分离有机片段与盐元数据；`dict` 为空表示非简单盐 |

非职责：不识别官能团（L1）、不选母体、不做位次与名称组装（L2–L5）。

本层是全管线唯一接触原始文本之处，也是错误处理的第一道防线：空串、纯空白、非法 SMILES 语法、消毒或归一化抛出的异常，一律收敛为 `None`，由上层转成解析失败结果。

## 核心逻辑

### 互变异构归一化（`tautomer.py`）

入口 `normalize_amide_tautomer(mol)` 是双向归一化器：既把烯醇/烯硫醇式拉回酮式，也把环外亚胺式拉回胺式。先 `_assign(_enol_candidates(mol))` 取烯醇位点，再 `_assign(_amidine_candidates(mol))` 取亚胺位点，并排除与烯醇位点同中心碳者；无位点直接返回原 mol。有位点在 `RWMol` 上就地改键后 `SanitizeMol`，失败整体回退，成功再 `AssignStereochemistry`。可调判据为文件顶部的 `ENOL_X`(O/S) 与 `N_VALENCE`(3)。

- 烯醇侧 `_enol_candidates`：中心碳任取；可迁移杂原子由 `_is_amide_enol_x` 判——`ENOL_X`、中性、与碳单键、度为 1、带 ≥1 H。受体氮由 `_is_amide_enol_n` 判——中性、无显式 H、度+氢数未达 `N_VALENCE`、与碳成双键或芳香键。
- 亚胺侧 `_amidine_candidates`：中心碳在环内且带环外亚胺氮（`_is_imine_n`：中性、非碳邻居一律排除、可再容纳 1 个 H、不在环内），受体为环内氨基氮（`_is_amidine_n`：在环内、带 H、单键或芳香键连碳）。
- `_assign` 按候选受体氮数由少到多贪心配对，一个氮只服务一个位点。
- `_normalize_amide` / `_normalize_amidine` 就地搬迁 H 与键级；`_normalize_amide` 跳过羟基氢数不为 1 的位点。显式 H 记账必须保留——芳环 N 的隐氢受芳香性约束不会自动补，只能 `SetNumExplicitHs` 显式加，否则 Kekulé 奇偶不匹配。

### 预处理链（`preprocess`）

空白校验 → `MolFromSmiles(sanitize=False)` → `_strip_isotopes` → `SanitizeMol` → `AssignStereochemistry` → 互变归一 → 酸性质子重定位 → 二次互变归一。第二次是因为重定位把酰胺 O⁻ 变成中性 `C(OH)=N`，此时才出现新的烯醇位。`_strip_isotopes` 清除同位素标记：命名管线不产出同位素名，保留会让锚定键匹配失败。

### 电荷归一（`charge.py`）

`normalize_acid_charge` 把质子从强酸根供体搬到受体（`_relocate_proton`，只做电荷/H 记账），消毒失败或 SMILES 不变即停。供体由 `_acid_kind_of_oh` + `ACID_CENTERS` 判：中性 OH 所连非芳香中心带 ≥1 双键氧，酸类须在 `DONOR_KIND` 内，表键序即酸强度序。受体分两类：`_is_weak_anion` 判定的去质子化弱酸位（元素限 `ACCEPTOR_Z`，邻接 +1 电荷者与强酸共轭碱除外）为首选；仅当整分子净电荷为负且无弱酸位时，才把 `_is_base_n` 判定的中性脂肪胺氮当受体（保净电荷内的两性离子式，P-73.1.2）。供体按 `DONOR_KIND` 序优先、再取 canonical rank 最小排序，受体取同片段内 rank 最小者，以保证结果确定。循环逐次消耗一对供体/受体，上界为原子数；含通配原子（`*`）或无同片段受体时返回原 mol。

### 盐解离（`salt.py`）

`dissociate_salt` 先按片段数过滤（单片段直接视为无盐），再由 `_from_frags` 识别单原子反离子。`_frag_role` 给出三种角色：电荷 ≥ +1 且元素在 `METAL_ION_EN` 者为金属阳离子（P-71.2）；元素在 `HALO_Z` 内、电荷 -1 且无 H 者为卤素阴离子 X⁻；电荷 0 带 1 H 者为中性卤化氢 HX。其余片段归入有机侧。

三种反离子互斥，且各有额外的有机侧约束：金属须种类唯一、不与卤素反离子共存、有机片段为阴离子；卤素阴离子须种类唯一、不与 HX 共存、有机片段为阳离子；HX 须种类唯一、有机片段中性（P-71.3）。有机片段须至少一个且彼此一致（按非立体 SMILES 去重），份数记入 `n_org`。中英盐名由 `METAL_ION_EN`/`METAL_ION_ZH`、`HALIDE_EN`/`HALIDE_ZH`、`HALIDE_HX_EN`/`HALIDE_HX_ZH` 查得，份数 >1 时经 `_mult_word` 加 `MULT_EN`/`MULT_ZH` 数量前缀。不满足时返回 `(mol, {})`。

## 与其他层的契约

- 输出：消毒并完成互变/电荷归一、清除同位素与立体指派的 `Mol`；盐元数据含中英双语盐名，金属盐为 `metal`/`metal_zh`/`n_metal`/`n_org`，卤化物盐为 `halide`/`halide_zh`/`n_org`，氢卤酸盐为 `acid_salt`/`acid_salt_zh`/`n_org`。
- 消费方：`namer` 经 `preprocess` 取分子、经 `dissociate_salt` 分离有机片段，并把盐元数据注入 L1 的 `info["salt"]`，供 L2 盐门控（写 `salt_meta`）与 L5 拼盐名后缀；`namer._apply_salt_suffix` 经 `stems.join_metal_salt_names` 拼后缀，含氧酸中心母体自带盐组装故跳过。
- 两个接口相互独立：`preprocess` 只吃 SMILES 做解析与归一化，`dissociate_salt` 只吃已消毒的 `Mol` 做片段拆分；调用方先预处理再解离。盐解离不改变有机片段的原子索引，故 L1 之后的层无需感知盐的存在。
- 本层不保留任何输入文本，下游一律基于 `Mol` 与原子索引工作。盐数据契约详见 [[reference/core-data-contracts]]。
