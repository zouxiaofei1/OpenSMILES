# Layer0: 预处理 (Preprocessor)

> **源文件:** `layer0/` 5 个 `.py` (355 行) | **对外接口:** `preprocess` / `dissociate_salt`

## 概述

本层把外部 SMILES 字符串转成内部可用的 RDKit `Mol`，并完成盐解离。层内共享常量取自 `constants.py`。

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

`tautomer.py` 的可调判据集中在文件顶部三个常量：`ENOL_X`（可迁移杂原子 O/S）、`N_VALENCE`（中性氮价态上限）、`MIN_RING_HETERO`（芳香中心碳归一所需的最少环内杂原子数）。

本层是全管线唯一接触原始文本之处，也是错误处理的第一道防线：空串、纯空白、非法 SMILES 语法、消毒或归一化抛出的异常，一律收敛为 `None`，由上层转成解析失败结果。

## 核心逻辑

### 互变异构归一化（`tautomer.py`）

唯一入口 `normalize_amide_tautomer(mol)` 是双向归一化器：既把烯醇/烯硫醇式拉回酮式，也把环外亚胺式拉回胺式。先 `_assign(_enol_candidates(mol))` 取烯醇位点，再 `_assign(_amidine_candidates(mol))` 取亚胺位点，并用集合差排除与烯醇位点同中心碳者；无位点直接返回原 mol；有位点在 `RWMol` 上就地改键后 `SanitizeMol`，失败整体回退，成功再 `AssignStereochemistry`。

- 烯醇侧 `_enol_candidates`：中心碳非芳香，或虽是芳香但落在含杂原子环内（`_in_poly_hetero_ring`，阈值 `MIN_RING_HETERO`）。可迁移杂原子 X 由 `_is_amide_enol_x` 判：`ENOL_X`(O/S)、中性、与碳单键、带 ≥1 H。受体氮由 `_is_amide_enol_n` 判：中性、无显式 H、度+氢数未达 `N_VALENCE`(3)、与碳成双键或芳香键。
- 亚胺侧 `_amidine_candidates`：中心碳在环内且带环外亚胺氮（`_is_imine_n`：中性、非碳邻居一律排除、可再容纳 1 个 H、不在环内），受体为环内氨基氮（`_is_amidine_n`：在环内、带 H、单键或芳香键连碳）。
- `_assign`：按候选受体氮数由少到多贪心配对，一个氮只服务一个位点。
- `_normalize_amide` / `_normalize_amidine`：就地搬迁 H 与键级。显式 H 记账必须保留——芳环 N 的隐氢受芳香性约束不会自动补，只能 `SetNumExplicitHs` 显式加，否则 Kekulé 奇偶不匹配。
- 无位点、或 `SanitizeMol` 抛异常时返回原 `mol`；成功时对结果重跑 `AssignStereochemistry`。

### 预处理链（`preprocess`）

空白校验 → `MolFromSmiles(sanitize=False)` → `SanitizeMol` → `AssignStereochemistry` → 互变归一 → 酸性质子重定位 → 二次互变归一。第二次是因为重定位把酰胺 O⁻ 变成中性 `C(OH)=N`，此时才出现新的烯醇位。

### 电荷归一（`charge.py`）

`normalize_acid_charge` 依 `ACID_CENTERS` 查表判成酸中心，把质子从强酸根供体搬到同片段内的弱酸受体（`_relocate_proton`，元素限 `ACCEPTOR_Z`），消毒失败或 SMILES 不变即停。字典键序即酸强度序。供体按 `DONOR_KIND` 序优先、再取 canonical rank 最小排序，受体取同片段内 rank 最小者，以保证结果确定。循环逐次消耗一对供体/受体，上界为原子数；含通配原子（`*`）或跨片段配对时直接返回原 mol。

### 盐解离（`salt.py`）

`dissociate_salt` 先按片段数过滤（单片段直接视为无盐），再由 `_from_frags` 识别单原子 +1 碱金属盐（`ALKALI_EN`，金属种类须唯一）与 HCl 盐（`_is_hcl_frag`：中性带 1 H 的 Cl 或 `[Cl-]`）。其余片段归入有机侧，只有恰好一个有机片段时才判为简单盐，否则返回 `(mol, {})`。碱金属盐时不接受 HCl 共抗衡离子；中英盐名由 `METAL_ZH` 等表查得。

## 与其他层的契约

- 输出：消毒并完成互变/电荷归一与立体指派的 `Mol`；盐元数据含中英双语盐名，碱金属盐为 `metal`/`metal_zh`/`n_metal`，HCl 盐为 `acid_salt`/`acid_salt_zh`。
- 消费方：`namer` 经 `preprocess` 取分子、经 `dissociate_salt` 分离有机片段，并把盐元数据注入 L1 的 `info["salt"]`，供 L5 拼接盐名后缀。
- 两个接口相互独立：`preprocess` 只吃 SMILES 做解析与归一化，`dissociate_salt` 只吃已消毒的 `Mol` 做片段拆分；调用方先预处理再解离。盐解离不改变有机片段的原子索引，故 L1 之后的层无需感知盐的存在。
- 本层不保留任何输入文本，下游一律基于 `Mol` 与原子索引工作。
