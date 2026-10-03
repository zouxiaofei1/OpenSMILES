# 原子归属跟踪 (Atom Ownership)

> **概念层级:** 跨层核心概念 | **链路:** L1 出事实 → L2 合边界 → L3 消耗 | **核心数据:** `owned_atoms: frozenset[int]`、`FunctionalGroupOccurrence`、`ClaimedBlock`、`CoverageLedger`

原子归属划定"母体 / 取代基"的分界线：分子里每个重原子（非氢原子）或者由母体声明所有权，或者落在边界外、作为一个连通块被 L3 claim。归属结果 `owned_atoms` 由 L2 冻结为 `frozenset[int]`，L3 只读它做几何切分，不回写。

## 三段式

| 层 | 职责 | 入口 |
|----|------|------|
| L1 | 为每条 occurrence 出特征原子集与母体锚点 | `build_inventory` → `_one` |
| L2 | 把锚点对齐骨架，合成不可变 `owned_atoms` | `_kind_fg_atoms` → `finalize_parent_ownership` |
| L3 | 在边界外枚举并命名取代基 | `extract_substituents` → `iter_claims` |

`owned_atoms` 同时是 L2 候选排序的基准：`_p45_2_prefix_count` 以边界外的 claim 个数作 P-45.2.1 键，`_reorder_p45_2` 据此降序、只保留并列最大组；`select_parent` 再对磷酸候选施加 `_reorder_oxo_ester_side` 定向（P-67.1.3，取糖 / 多元醇侧），返回该组。归属边界越紧（边界外组分越多），计数键越高。

## L1：事实来源

`FunctionalGroupInventory` 由 `FG_SPECS` 投影出两张表：

- `_FG_KEYS` —— FG 类别键，顺序即注册顺序，决定清单条目的产出顺序；
- `_ANCHOR_KEYS` —— 只收录声明了 `anchors` 的类别，取值即 occurrence payload 里承载锚点的键名。胺的锚点是多臂锚点（`surr_idx` 上的全部碳邻居），经 `_indices` 展平为索引集合。

特征原子由 `_characteristic_atoms` 给出，分两条路：

- **通用规则** `center_surr_atoms`：`{center_idx} ∪ surr_idx`；
- **例外表** `FG_ATOM_FNS`：登记六类。`oxoacid` 与 `sulfonamide` 共用 `_oxoacid_atoms` = 锚点碳 + 中心 P/S（payload 的 `oxo_z`）+ 中心的非碳邻居，中心的**碳臂不进特征原子**，留给链或取代基侧；`cation` / `azanide` / `heterane` 共用 `_cation_atoms`，只占中心自身，周边原子全退为取代基 / 前缀；`nitrile` 用 `_nitrile_atoms`，只取腈碳与三键氮，R 侧连接原子不入。

以磺酰胺 `CS(=O)(=O)N` 为例：`center_idx` 是锚点碳、`oxo_z` 是 S，特征原子为 `{C, S, =O, =O, N}`。N 与锚点碳相隔两跳，靠邻接扩展取不到——这正是含氧酸族走整取的原因。`center_surr_atoms` 同时是通用规则的唯一实现：未在 `FG_ATOM_FNS` 登记例外的类别一律退回它。

各类别的锚点键与特征原子口径：

| FG 类别 | anchors key | 特征原子 |
|---------|-------------|----------|
| radical / acyl | `center_idx` | `center_surr_atoms` |
| acid | `center_idx` | `{C, =O, OH}` |
| ester | `center_idx` | `{C, =O, O_single}`，不含烷氧臂碳 |
| acyl_halide / amide / aldehyde / ketone / thione | `center_idx` | `center_surr_atoms`，即羰基（或 C=S / 酰卤）中心与其周边原子 |
| nitrile | `center_idx` | `_nitrile_atoms`：腈碳与三键氮，R 侧连接原子不入 |
| alcohol / thiol / amine | `surr_idx` | `center_surr_atoms`；锚点是连接碳，O / S / N 由一跳扩展进入 |
| oxoacid / sulfonamide | `center_idx` | `_oxoacid_atoms` |
| cation / azanide / heterane | `center_idx` | `_cation_atoms`：仅中心原子 |

锚点键决定"哪条原子是母体入口"：羧酸族取羰基碳，醇 / 硫醇 / 胺取连接碳而非中心杂原子。

事实的聚合发生在清单构建期：`build_inventory` 遍历 `_FG_KEYS`，为每条检测条目调用 `_one`，一次算出 `parent_anchors` 与 `characteristic_atoms` 并固化进 `FunctionalGroupOccurrence`（字段为 `id` / `group_class` / `characteristic_atoms` / `parent_anchors` / `payload` / `demoted`）。`demoted` 标记 P-41 仲裁降级为前缀叶的条目，`inventory.occurrences()` 按类查询时排除它们，被降级的羧基 / 氰基因此不参与主基团竞争，其碳也退为边界外组分。

扩展一个新 FG 类别只需三处落点，不必触碰 Layer2 的归属逻辑：`FG_SPECS` 加一条带 `anchors` 的注册项、检测表加对应模式、特征原子不符合通用规则时在 `FG_ATOM_FNS` 加例外。

## L2：几何边界合成

`_kind_fg_atoms(parent, mol)` 读 `parent["principal_expression_facts"]`（缺失即返回空集）与 `parent["principal_occurrences"]`（全量主官能团 occurrence，不限本骨架覆盖），按三态处理：

1. **取 seed**：`seeds = anchors ∩ chain`；锚点全在骨架外时，改取与骨架原子相邻的锚点。该兜底承接的是外环官能团——苯甲酸 `c1ccccc1C(=O)O` 的羧基碳不在环骨架上，锚点经邻接命中后成为 seed。
2. **`facts.group_class ∈ OXO_FG_CLASSES = {oxoacid, sulfonamide}`**：只取本骨架覆盖的 occurrence（`parent["covered_principal_ids"]`）的 `characteristic_atoms` 并入 seed，**不做邻接扩展**。这两类的中心 P/S 不是碳骨架成员，没有可外借的臂，跨两跳的氧 / 氮必须整取。
3. **其余类别**：沿 seed 的邻居做**一跳**扩展，只并入既落在全体 occurrence 特征原子集内、又不是锚点的邻居，不回代。

`finalize_parent_ownership(parent, mol)` 折叠出 `owned_atoms`：

- 若 `parent["kind"] ∈ OXO_CENTER_KINDS`（`phosphate` / `phosphonate` / `sulfate` / `boronic`）：只取主官能团特征原子，**不含 `chain`**——链仅供编号，臂一律退为取代基。
- 否则取 `frozenset(_chain_atoms(parent) | _kind_fg_atoms(parent, mol))`，即骨架原子与主官能团特征原子的并集。
- 结果带缓存语义：`owned_atoms` 已是 `frozenset` 时直接返回原 dict，不重算。

两处调用点：`select_parent` 经 `_finalize_ranked`（先 `pack_parent_stem`）固化每个候选；`namer._prepare_candidate` 在进入 L3 前再兜一次。

### 骨架原子

`_chain_atoms(parent)` 取 `parent["chain"]` 转成集合，该字段由表达阶段写为骨架原子 id 列表，**环骨架的环原子同样在列**——苯酚的 `chain` 覆盖整个环，苯甲酸乙酯的 `chain` 覆盖环与羰基碳。骨架原子是所有权的基础层：任何候选至少拥有自己的骨架，开链母体与环系母体在这一点上无差别。候选未能选出骨架时 `chain` 为空，此时提前失败，不进入取代基提取。

实测对照（`owned_atoms` 的取值随 `kind` 与骨架而变）：

| SMILES | kind | chain | owned_atoms |
|--------|------|-------|-------------|
| `CCO` | alcohol | `[0, 1]` | `{0, 1, 2}` |
| `CC(=O)C` | ketone | `[0, 1, 3]` | `{0, 1, 2, 3}` |
| `CC(=O)C(C)=O` | ketone | 单个乙酰基骨架 | `{0, 1, 2, 3, 4, 5}`，两条 occurrence 的特征原子全并 |
| `c1ccccc1C(=O)O` | acid | `[0..5]` | `{0..8}`（环 + 羧基 C、=O、OH） |
| `CCOC(=O)c1ccccc1` | ester | 环 `[5..10]` | 环 + 羰基 C、=O、酯 O；乙氧臂在边界外 |
| `OP(=O)(O)O` | phosphate | `[1]` | `{0, 1, 2, 3, 4}`（P 与 4 个 O） |
| `COP(=O)(O)OC` | phosphate | `[2]` | `{1, 2, 3, 4, 5}`，两个甲氧臂在边界外 |
| `CP(=O)(O)O` | phosphonate | `[0, 1]` | `{1, 2, 3, 4}`——链不入所有权，甲基碳退为取代基 |
| `CS(=O)(=O)N` | sulfonamide | `[0]` | `{0, 1, 2, 3, 4}` |
| `CCN(CC)CC` | amine | `[0, 1]` | `{0, 1, 2}`，另两条乙基臂在边界外 |

### 环骨架与环内官能团

环骨架的环原子同样进 `chain`，因此环上的杂原子随 `chain` 一并进入所有权。环内单碳羰基的锚点即环内羰基碳、本身在 `chain` 内直接成为 seed，其环内杂原子邻居由 `chain` 侧纳入：硫代内酯 `O=C1CCCS1` → `chain=[1..5]`、`owned_atoms={0..5}`，环内 S 属于母体。

环内零碳羰基（内酯）走同一路径：`C1COC(=O)C1` → `chain=[0,1,2,3,5]`、`owned_atoms={0..5}`，环内酯氧一并归母体。环外部分（如 N-酰基环胺的吡咯烷环、内酰胺 N 上的烷基）留在边界之外，由 Layer3 作为取代基提取——`O=C1CCCN1C` 的 N-甲基即落 `other` 槽位。

### facts 字段与全量 occurrence

`PrincipalExpressionFacts` 记录 `group_class` / `multiplicity` / `relation` / `occurrence_ids` 与三个原子集：`characteristic_atoms`（全体 occurrence 特征原子的并集）、`anchor_atoms`（原锚点）、`attachment_atoms`（骨架内附着原子，骨架外的锚点取其骨架内邻居）。`_kind_fg_atoms` 只消费 occurrence 上的两个原子集与 `parent["covered_principal_ids"]`，不重算邻接关系。

母体 dict 由 `_parent_dict` 组装，同时写入两个 occurrence 视图：`covered_principal_ids` 是本骨架覆盖的条目 id 元组，`principal_occurrences` 则是**全量**主官能团 occurrence。后者取全量而非本骨架子集，是因为未被本骨架覆盖的同级基团仍属母体——漏掉其异原子会被 L3 切成假前缀。`OXO_FG_CLASSES` 分支则反向只认 `covered_principal_ids`，故 `CP(=O)(O)O` 的甲基碳虽在 `chain` 内也不进所有权。

### 边界的反向使用：P-45.2.1 计数

`owned_atoms` 不只是 L3 的输入，也是 L2 候选排序的基准。`_p45_2_prefix_count` 就地取用 L3 的 `iter_claims`，以边界外的 claim 个数作 P-45.2.1 键；`_reorder_p45_2` 按 `(计数降序, 缩合含氧酸链内桥氧数降序, 原序)` 稳定重排，第二键由 `_condensed_rank` 给出（`oxo_kind` 为 `phosphate` 或 `sulfate` 时取链中中心原子的 `analyzer._oxo_bridge_arms`，其余恒 0）。`tied=True` 时只保留计数最高的并列组。`select_parent` 随后对整组施加 `_reorder_oxo_ester_side`：当候选全为 `oxo_kind == phosphate` 时，按 `_oxo_ester_side_score`（臂上芳香含氮环数升序、臂内氧数降序）定向，把糖 / 多元醇侧选作母体（P-67.1.3）。归属边界越紧，P-45.2.1 计数键越高，越可能在并列中胜出。

## L3：消耗与槽位

`extract_substituents(info, parent)` 在所有权边界外取 claim：

- `iter_claims(mol, owned)` 遍历边界外的重原子连通组分，结果按 `(attach_parent, root, slot.value)` 排序，输出顺序确定；
- `_owned_edges` 列出组分与母体之间的全部 `(owned 原子, 组分原子)` 边，`_try_claim` 与 `claim_block` 要求这些边的 owned 端只有唯一一个（多回接点即桥连，不成 claim），并取最小边作 `(attach_parent, root)`；
- `_has_dbl_o_edge` 滤掉"经双键连到 owned 内非碳重原子"的组分（环内 S / P 的 `=O` 无主 FG 承接，作例外放行），主官能团成分不切成假羟基侧链；
- 命名返回 `None` 的 claim 被静默跳过，其原子只能体现为 gap。

`claimable_block` 提供两种结构：`ClaimedBlock`（`slot` / `attach_parent` / `root` / `atoms`，冻结 dataclass）与 `SideSlot`（4 值：`chain_c` / `ring_c` / `amine_n` / `other`）。槽位到取代基 kind 由 `constants.CLAIM_KIND` 映射（`amine_n` → `n_block`，`ring_c` / `chain_c` → `alkyl`），未登记槽位退回 `"side"`。

O- 侧臂不走槽位通路：`extract_substituents` 先算 `o_side`（`kind ∈ ESTER_O_SIDE_KINDS` 或存在 `o_idx` 字段），`_append_named` 仅当 claim 的 `attach_parent` 原子序数等于 `side_z` 时打标；`side_z = 16 if parent["thio_side"] else 8`，硫代酯的酯侧臂因此按 S 识别。实测 `CCOC(=O)CC` 的乙氧臂、`COP(=O)(O)OC` 的两个甲氧臂、`CC(=O)SC` 的甲硫臂都落 `other` 槽位并带 `o_side`，kind 取 `side`。

### 槽位判据

`derive_slot` 只看母体侧连接原子的角色：`_is_amine_n` 判定非芳香、非环员的 N 才是 `amine_n`；碳原子按是否成环分 `ring_c` / `chain_c`；其余杂原子落 `other`。环员 N 与芳香 N 因此走 `other`——内酰胺 N 上的甲基、环胺的环外臂都按位次或普通前缀处理，不用 N- 前缀。

kind 落定后还有一道改写：`sub_from_named(named, mol, parent)` 在 kind 属 `N_PREFIX_KINDS`（`n_block` / `n_alkyl`）且附着原子是环员、或附着原子落在母体骨架 `chain` 内时，改用 `ring_c` 映射的 `alkyl`，即改以数字位次定位。实测乙酰苯胺 `CC(=O)Nc1ccccc1` 的苯基臂落 `amine_n` → `n_block`；`CCN(CC)CC` 两条链外乙基臂同样落 `amine_n`。

### claim 枚举

`iter_claims` 的组分来源是 `side_roots`（母体原子的外部重原子邻居）逐一起 `cut_block`，割出的连通块按原子集去重（`_unique_components`），再对每个组分试建 claim。`cut_block` 拒绝穿越传入的 owned 集，因此 claim 原子按构造位于边界外、彼此不重叠。`extract_substituents` 实际传入 `cut_owned = owned_atoms ∪ _amidine_n_owned ∪ _hydrazide_n_owned`：脒 / 胍母体的亚胺氮与酰肼远端 N 临时并入切分边界，其上的臂从该 N 外侧键起切，避免被吸进母体或丢失。

`_side_arm_claim` 再把切在桥杂原子（O / S）上、非链内的侧臂并回该桥原子、改挂到链上原子，使硫代酯等的支臂拿到可定位次。命中的 claim 经 `SubstituentNamer` 取到名字后，由 `sub_from_named` 封装为取代基 dict：`kind`、`n_carbons`（块内碳数）、`attach_idx`（即 `attach_parent`）、`atoms`（排序后的原子索引）、`en` / `zh` / `paren`，O- 侧臂再补 `o_side` 标记，酰肼两端 N 的取代基由 `_mark_hydrazide_primes` 补 `n_prime`（0 = N′）。归属信息到此转化为下游可消费的字段，`atoms` 是后续覆盖台账的唯一输入面。

## 覆盖完整性

`build_coverage_ledger(mol, owned_atoms=..., names=...)` 给出两个诊断集：

- `gap = heavy − owned_atoms`：口径只扣所有权、不扣 claims，未被覆盖的重原子即缺口；
- `overlap = {i | counts[i] > 1}`：`counts` 以 `owned_atoms` 起算，再叠加每个 `SubstituentName.claim.atoms`。

管线内 `namer._prepare_candidate` 以 `names=[]` 调用，`gap` 因此等价于"`owned_atoms` 是否覆盖全部重原子"，`overlap` 恒为空集。该布尔值作为候选三元组的第三项返回，`_try_phase` 解包而不消费，成功候选在 `meta` 标注 `fallback = "no_coverage_gate"`。真正的前置失败条件只有两条：`owned_atoms` 为空，或既无 `chain` 又无环。

带真实 `names` 的台账在诊断路径构建：`server/backend/routes_debug.py` 用取代基原子集重建 `SubstituentName` 后调 `build_coverage_ledger`，把 `coverage` 与排序后的 `owned_atoms` 一并输出。测试侧以 `names=[]` 校验并列母体组内至少有一个候选覆盖完整。

gap 随母体类别而变：`CC(=O)C(C)=O` 的酮母体覆盖全部 6 个重原子，`complete` 为真；无 FG 归属规则的分子则留缺口——`CCOCC`（醚）选出 `kind=alkane`、`owned_atoms={0, 1}`，醚氧与另一端碳落 gap；`CSCCO` 的醇母体只拥有 `{2, 3, 4}`，S 与甲基硫臂落 gap。缺口由 L3 的 claim 在整名层面补齐，台账只作诊断。

### 所有权的不变性

`owned_atoms` 的类型是 `frozenset[int]`，`finalize_parent_ownership` 对其幂等：字段已是 frozenset 时返回同一 dict 对象，不重新求值。候选组内各个候选拥有各自的 `owned_atoms`（记录的是自身母体范围），组内同分不代表边界相同。提取取代基不修改母体边界，该契约由类型约束与测试共同固定。

### 归属边界与命名实体

归属只解决"哪些原子属于谁"，不解决"叫什么"。母体侧的名字由 L4 编号与 L5 词干组装承担，边界外的 claim 由取代基命名通道产出词条；两侧在整名阶段合流。因此 `owned_atoms` 的正确性判据是几何的：重原子被恰好一次声明，而命名是否成立由后续层独立决定。

相关页面：[[architecture/layer1-analyzer]]、[[architecture/layer2-parent-selector]]、[[architecture/layer3-substituents]]、[[architecture/layer4-numbering]]、[[reference/core-data-contracts]]。
