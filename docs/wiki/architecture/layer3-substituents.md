# Layer3: 取代基提取 (Substituent Extraction)

> **源文件:** 7 `.py`（525 行） | **对外接口:** `extract_substituents(info, parent, *, cache)`、`iter_claims` / `ClaimedBlock` / `SideSlot`、`SubstituentNamer` / `SubstituentName` / `SubstituentBackend`、`name_as_substituent`、`build_anchor_submol`、`build_coverage_ledger` / `CoverageLedger`

## 概述

为母体所有权边界之外的重原子组分逐一命名，产出取代基列表供 L4 编号与 L5 组装。入口 `extract_substituents` 只做汇总与参数推导：claim 的拓扑枚举在 `claimable_block`，命名在 `substituent_namer`。

**输入:**
- `info["mol"]`（RDKit `Mol`）、`info["root_ctx"]`（`(root_mol, to_root)`，供递归取代基回根分子校正 R/S）
- `parent`: 只读 `owned_atoms`（为 `None` 即返回 `[]`）、`kind`、`o_idx`、`thio_side`
- `cache: CommonNameCache | None`

**输出:** `list[dict]`，每项含 `kind`（`halo`/`alkyl`/`n_block`/`side`）、`n_carbons`、`attach_idx`（母体侧连接原子）、`atoms`（升序索引）、`en`/`zh`、`paren`（是否需围栏）；O 侧臂另带 `o_side=True`。

**职责边界:** 不选母体、不编号、不拼名，只负责边界外的组分拓扑与命名。命名失败的 claim 在 `_append_named` 处静默跳过。

## claim 枚举与槽位

`iter_claims(mol, owned_atoms)` 是唯一提取路径：`_unique_components` 用 `side_roots` / `cut_block` 切出并去重外部重原子组分，`_try_claim` 逐个建 claim，最后按 `(attach_parent, root, slot.value)` 排序返回 `ClaimedBlock`。

- `_try_claim` 取最小 `(attach_parent, root)` 边为连接点；`_has_dbl_o_edge` 排除「双键 O 连到非碳 owned 原子」的组分（主官能团成分不 claim），环内 S 的 `=O` 例外——它无主 FG 承接，须按 oxo 前缀 claim。
- `claim_block` 要求 `attach_parent` 落在 owned 内，且组分对外连接点唯一（`_attach_parents_of` 计数为 1），否则放弃。
- `SideSlot` 四值 `CHAIN_C` / `RING_C` / `AMINE_N` / `OTHER`，由 `derive_slot` 仅从连接原子角色导出：`_is_amine_n` 认非芳香、非环员、电中性的 N 为 `AMINE_N`（环 N 用环上位次定位，N⁺ 作母体阳离子不写 N- 前缀），碳按是否成环分 `RING_C`/`CHAIN_C`。
- `CLAIM_KIND` 把槽位映射为取代基 kind：`amine_n → n_block`，`ring_c`/`chain_c → alkyl`，未登记槽位取 `side`。
- `sub_from_named` 先用 `NAME_KIND`（保留名暗含 kind，如卤素 → `halo`），再回退 `CLAIM_KIND`；命中 `N_PREFIX_KINDS` 但附着原子为环员时改判 `ring_c`，改走环上位次而非 N- 前缀。

## 双后端命名

`SubstituentNamer` 按序尝试后端取首个命中；`SubstituentBackend` Protocol 只要求 `name` 属性与 `try_name(mol, claim)`。

- `RetainedBackend`（`name="retained"`）: 调 `tools.anchored_table.anchored_lookup(mol, claim.atoms, claim.root)` 查 canonical-SMILES 锚定表，命中保留叶子名。
- `RecursiveBackend`（`name="recursive"`）: 调 `name_as_substituent` 走有界递归 cut → free-name → yl_form，携带 `cache` 与 `root_ctx`。

`_named` 把后端命中的 `(en, zh, paren)` 封成 `SubstituentName`，双语名任一为空即视为未命中。

**递归实现:** `submol_build.build_anchor_submol` 复制组分原子与内部键，在连接原子处挂 dummy 锚点（`*`，键型取与母体间的真实键型，`=CH2` 这类双键叶用双键），照搬子图内双键 E/Z 后消毒；连接点不在 `atoms` 内或消毒失败则返回 `None`。`_radical_yl_from_sub` 命中后用 `_fix_rs_with_real` 在原始根分子上重算 R/S（连接点被 dummy 顶替后 CIP 会翻转），位次改用整体编号标签，单原子链省略位次号。片段自由基名写入 `cache`，但立体随宿主根重算，不跨根共享。

**查表:** `tools.anchored_table._REGISTRY` 是叶子保留表，`RetainedSubstituent(en, zh, anchored=(), paren=False)` 存双语名、锚定键与是否需括号；`anchored_key` 由 `memo` 缓存，`_ANCHOR_INDEX` 建键索引，`_WHOLE_ONLY_KEYS` 中的键由 L3 跳过。registry 未命中且给出 `attach_old` 时，`anchored_lookup` 回退 `_alkoxycarbonyl`：`_ester_o_side` 判定块恰为「连接点羰基碳 + 一个双键 O + 一个单键 O + 前端」时，把 `-C(=O)-O-R` 收成 `<R>oxycarbonyl` / `<R>氧羰基`（前端名须以 `oxy` 结尾，中文侧经 `zh_bridge_root` 拼接），并返回 `paren=False`。

## 围栏规则

`as_substituent._radical_yl_from_sub` 在 free-name 命中后判定 `requires_parentheses`（`name_as_substituent` 在连接点非重原子时直接返回 `None`，避免 dummy 锚点自复制致无限递归）。

- **复合前缀必括**（P-16.5.1.1）: `hit.meta["parent_substituent_count"] > 0` 即需括号；免括白名单为裸 `"phenyl"`、`SIMPLE_ALKOXY_NO_PAREN`、`AMIDO_RETAINED_EN`（P-66.1.1.4.3）、`SIMPLE_BRIDGE_YL_NO_PAREN`。
- **已自含围栏**: `hit.meta["bridge_self_enclosed"]` 为真表示前端名自带围栏（如 `(4-甲氧基苯基)磺酰基`），则不对整名加括号。
- **简单前端 + O/S 桥**（P-63.2.1/.2.2）: 名以 `oxy`/`sulfanyl` 结尾且需括号时，由 `_obridge_front_simple` 把前端装成 `ClaimedBlock`（`slot=OTHER`）交给同一条 `SubstituentNamer` 链路（retained → recursive）判简单性，为真则免括。
- **双原子桥**: 名以 `DIATOMIC_BRIDGE_YL`（`diazenyl` / `disulfanyl`）结尾时跳过上述 O/S 桥免括，围栏由 L5 依前端定形。

## O-侧臂标记

`_append_named(mol, claim, namer, out, *, o_side=False, side_z=8)`: 命名成功后，仅当 `o_side` 为真**且** `claim.attach_parent` 原子的元素序数等于 `side_z` 时，才给结果 dict 加 `o_side=True`。标记按连接原子元素判定，`side_z` 由调用方给出。

`extract_substituents` 为所有 claim 一次算好两个参数：

- `o_side = parent.get("kind") in ESTER_O_SIDE_KINDS or parent.get("o_idx") is not None`；`ESTER_O_SIDE_KINDS = {ester, phosphate, phosphonate, sulfate, sulfonate}`。苯甲酸酯类（苯 base + ester FG）无 kind 表项，靠 `o_idx` 识别。
- `side_z = 16 if parent.get("thio_side") else 8`: 硫代酯的酯侧臂元素为 S（P-65.6.3.3.7.1），故判据由 O 换成硫。

L5 以 `o_side` 挑出侧臂拼酯/含氧酸整名，见 [[architecture/layer5-name-assembly]]。

## 覆盖台账

`build_coverage_ledger(mol, *, owned_atoms, names)` 返回 `CoverageLedger(owned_atoms, named_claims, gap, overlap)`，只统计重原子：

- `gap = 全部重原子 − owned_atoms`
- `overlap` = 归属次数 > 1 的原子（owned 记 1 次，每个命名 claim 再各加 1 次）
- `complete` 属性 = 无 gap 且无 overlap

本层不据此门控主管线：消费方是调试接口 `server/backend/routes_debug` 与契约测试。

## 与其他层的契约

- **L2:** 消费 `owned_atoms`（为 `None` 即不提取）、`kind`、`o_idx`、`thio_side`；边界由 L2 定，见 [[architecture/layer2-parent-selector]] 与 [[concepts/atom-ownership]]。
- **L4/L5:** 每项以 `attach_idx` 声明挂点、`atoms` 声明覆盖的重原子，供 [[architecture/layer4-numbering]] 定位次；`kind` 选 L5 前缀通道（`n_block` 走 N- 前缀，`alkyl` 走位次通道），`paren` 定围栏，`o_side` 交 L5 拼侧臂。双语名分列 `en`/`zh`，见 [[concepts/bilingual-naming]] 与 [[reference/core-data-contracts]]。
- **工具层:** 依赖 `tools.block_cut.side_roots` / `cut_block` 与 `tools.anchored_table`；`tools.re.SUB_LOCANT_RE` 为 L5 与本层共用的位次段判据。
