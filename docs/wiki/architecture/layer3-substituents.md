# Layer3: 取代基提取器 (Substituent Extractor)

> **最后更新:** 2026-09-04 | **源文件:** 9 `.py` (1020 行) | **公开 API:** `extract_substituents(info, parent, *, name_mode, cache) -> list[dict]`

## 概述

Layer3 是 NamePredict 6 层流水线中的第三层，负责从母体结构中识别、切除并命名所有非母体原子作为取代基（substituents）。其核心职责：给定 layer1 的分析结果 `info` 和 layer2 选出的母体 `parent`，提取所有附加在母体链/环上的侧链和官能团，为每个取代基生成中英双语名称，并通过覆盖台账（Coverage Ledger）验证所有重原子均已覆盖、无重叠冲突。

**输入:**
- `info: dict` — layer1 的分析结果，包含 `mol` (RDKit Mol 对象)、羟基列表 `hydroxyls`、酮基列表 `ketones`、硝基列表 `nitros`、胺基列表 `amines`、醚列表 `ethers` 等
- `parent: dict` — layer2 选出的母体，包含 `kind`（母体类型）、`chain`（主链原子索引列表）、`owned_atoms`（母体声明拥有的所有原子 frozenset）、`ring`（环原子）、以及各种 FG 定位符
- `name_mode: str` — 命名模式（`"general"` 保留名 / `"pin"` PIN 首选名）

**输出:**
- `list[dict]` — 取代基列表，每个元素包含 `kind`（取代基类别）、`en`/`zh`（中英名称）、`attach_idx`（连接点原子索引）、`atoms`（取代基包含的原子列表）、`paren`（是否需要括号）、`n_carbons`（碳原子数）、`backend`（命名后端标识）等字段

**在流水线中的位置:**

```
layer1 analyze → layer2 select_parent → layer3 extract_substituents → layer4 number → layer5 assemble
```

> **源:** `src/namepredict/namer.py` — `extract_substituents` 由 `_assemble_candidate` 调用，随后 `_ledger_complete` 用 `build_coverage_ledger` 验证覆盖完整性。

## 核心逻辑

### 提取主流程 (extract_substituents)

`extract_substituents` 是 layer3 的入口函数。主流程为 **三段流水线**：核心 FG 提取 + 锚定查表烷基 + 覆盖补全（`extract_claimed_sides`）。

```python
def extract_substituents(info, parent, *, name_mode="general", cache=None) -> list:
    mol, chain = info["mol"], parent.get("chain") or []
    base = (
       _extract_core_subs(info, parent)                       # 卤素 + OH + NH2 + 氧代
       + _extract_alkyls_no_aryl(mol, chain, name_mode=name_mode)  # anchored_table 查表烷基
    )
    owned = parent.get("owned_atoms")
    if owned:
        base = [_with_full_atoms(mol, owned, s) for s in base]   # side_atoms 全连通补全
    return base + extract_claimed_sides(info, parent, base, name_mode=name_mode, cache=cache)
```

> **源:** `src/namepredict/layer3/substituent_extractor.py:219`

**第一阶段: 核心官能团提取 (`_extract_core_subs`)** — 按固定顺序收集卤素（F/Cl/Br/I，直接连接在链碳上）、羟基（仅当母体不是醇/酚类）、氨基（仅当母体不是胺/苯胺类，委托 `amino_side.py`）、氧代（仅当母体不是酮/醌类）。每个提取器都检查母体类型 (`_PARENT_OH_KINDS`, `_PARENT_NH2_KINDS`, `_PARENT_OXO_KINDS`) 以及 principal-expression 的附着位（`_principal_attachments`），当母体本身以该官能团为主要官能团时，链上相应基团视为母体骨架而非取代基。

卤素提取有特殊过滤 (`_filter_fg_halos`)：对酰氯/酰溴母体排除官能团自身的卤素；对功能性醚（如 HFIP）跳过醚臂上的氟原子。

> **源:** `src/namepredict/layer3/substituent_extractor.py:198`（`_extract_core_subs`）

**第二阶段: 烷基侧链提取 (`_extract_alkyls_no_aryl`)** — 遍历母体主链上的每个碳原子，通过 `_side_starts`（`carbon_neighbors`，来自 `tools/chain`）识别链碳的界外邻居。对每个侧链起点调用 `_one_alkyl` → `_one_anchored_alkyl`：用 `tools.block_cut.side_atoms` 计算全非母体连通分量，再经 `tools.anchored_table.anchored_entry` 查表。**仅 `kind == "alkyl"` 的条目被认领**（纯碳侧链）；cyano/nitroso 等杂叶留给 claim 补全器。未命中则静默回退（由 `extract_claimed_sides` 的 SubstituentNamer 兜底）。

> **源:** `src/namepredict/layer3/substituent_extractor.py:91, 189`（`_one_anchored_alkyl` / `_extract_alkyls_no_aryl`）

**第三阶段: 声明侧链补全 (`extract_claimed_sides`)** — 覆盖补全机制（`claim_extract.py`）。遍历 `iter_claims`（`claimable_block.py`）生成的所有 `ClaimedBlock`（ownership 边界处的未覆盖原子块），对每个未被前面阶段覆盖的 claim 调用 `SubstituentNamer.name()` 命名。**仅原子已被覆盖的 claim 跳过**（`_should_skip`；酰胺 N 不再特判跳过）。这是 layer3 的"安全网"——任何逃过前两段的侧链原子（芳基、烷氧基、N-侧链、复杂分支烷基、杂环等）最终在这里被捕获命名。

> **源:** `src/namepredict/layer3/claim_extract.py:94`

### 保留取代基查表 (anchored_table, tools/)

核心机制是**锚定 canonical-SMILES 查表**。`submol_build.build_anchor_submol` 在取代基连接位打上 dummy 原子 `*`，使 canonical SMILES 同时编码形状与附着点（`*C(C)C` 异丙基 vs `*CCC` 正丙基；`*c1ccc(Cl)cc1` 4-氯苯基）。`tools/anchored_table.py`（201 行）维护**单一取代基注册表 `_REGISTRY`**（45 条 `RetainedSubstituent`）——基础取代基与 IUPAC 保留取代基统一入表，各带 `anchored` 锚定键；原内联 `ANCHOR_TABLE` 已并入 registry，**长链烷基（C5+）/环烷基/苯基/卤代烷基不再查表**，改走完整递归/radical 命名管线：

| 类别 | 条目示例 |
|------|---------|
| 基础取代基（halo/leaf） | `*F` fluoro、`*Cl` chloro、`*[N+](=O)[O-]` nitro、`*N=C=O` isocyanato、`*N=C=S` isothiocyanato |
| 线性 n-烷基 C1-C4 | `*C` methyl … `*CCCC` butyl |
| 支链/烯/炔保留基 | `*C(C)(C)C` tert-butyl、`*C(C)C` isopropyl、`*CC=C` allyl、`*CC#C` propargyl、`*CC(C)(C)C` neopentyl |
| 含杂原子保留基（leaf） | `*OC` methoxy、`*O` hydroxy、`*SC` methylsulfanyl、`*S(=O)(=O)O` sulfo、`*C=O` formyl、`*C#N` cyano、`*N` amino、`*N=N` diazenyl 等 |
| 芳烷基 | `*Cc1ccccc1` benzyl |

构建期 `_build_anchor_index`（`anchored_table.py:105`）把各条目 `anchored` 键 canonical 化后反查生成 `_ANCHOR_INDEX`；`anchored_entry`（`:170`）/`anchored_lookup`（`:182`）经 `_table_hit`（`:162`）命中即解析为 `(en, zh, paren, kind)`。registry 键经 `resolve_name`（`:128`）解析——其 general/pin 分派已注释停用，恒返回系统名（`systematic_en`/`systematic_zh`，如 isopropyl → propan-2-yl）。`kind` 分 `alkyl` / `aryl` / `halo` / `leaf` 四类——**提取器只认领 `alkyl`**，其余类别留给命名器/补全器。构建期校验 anchored 键的 canonical 形式与唯一性。registry 现有 **45 条**（新增 `formyl` 甲酰、anchored=`("*C=O",)`；删除 `methylsulfinyl` 甲亚磺酰基；`methylsulfonyl` 的 `systematic_en` 由 `methanesulfonyl` 修正为 `methylsulfonyl`）。

> **源:** `src/namepredict/tools/anchored_table.py`

### 取代基命名引擎 (SubstituentNamer)

当直接提取器无法命名或 claim 补全器认领的取代基出现时，`SubstituentNamer` 承担命名职责。它采用有序后端链（`_default_backends`，`substituent_namer.py:92-94`），首个命中即返回：

1. **RetainedBackend** — 保留名/锚定查表。`_try_anchored_lookup` 用 `tools.anchored_table.anchored_lookup` 查 `_REGISTRY`（`_ANCHOR_INDEX` 反查），命中即返回。

2. **RecursiveBackend**（`substituent_namer.py:65`）— 有界递归切割命名。`as_substituent.name_as_substituent` 将 claim 原子从母分子切出为 submol，作为独立分子跑完整 L1-L5 管道，再转 P-29 -yl 形式。递归深度上限 `max_depth=4`。`RecursiveBackend` 与 `SubstituentNamer`（`substituent_namer.py:97`）现接收 `root_ctx`（本层块原子 → 原始根分子索引映射，`claim_extract.py:101` 自 `info.get("root_ctx")` 注入 namer），贯穿传给 `name_as_substituent`，供 radical 取代基回根分子校正 R/S。

> 命名由 retained 与 recursive 两个后端承担（无 RootedTreeBackend）。

> **源:** `src/namepredict/layer3/substituent_namer.py:52-109`

### 通用 cut→free-name→yl 管道 (as_substituent / submol_build)

`as_substituent.py`（114 行）承载 RecursiveBackend 的底层 cut→free-name→yl 管道，核心是**按连接点类型分派的锚定 radical 优先路径**，并经 `root_ctx` 回根分子校正 R/S：

- `submol_build.py` 提供 `build_cut_submol`（诱导子分子 + attach 处 H 封端）与 `build_anchor_submol`（attach 打 dummy `*`，供 anchored SMILES 使用），定义 `CutSubmol` 数据类（`atom_map`/`inv_map`/`attach_new`/`attach_old`）。`_copy_bonds` 只重建键型会丢双键奇偶，故两条构建路径都调 `_carry_alkene_stereo`（`submol_build.py:43`）把诱导子图内 C=C 的 E/Z 标签照搬进子分子——子分子只是命名替身，标签取原分子即真实立体，不随配基被切/被 `*` 顶替而重判；两端引用邻接被保留则映射到新索引，被切配基由连在 sp2 端的 dummy 顶替（`_add_anchor`（`:117`）现返回 dummy 索引供其使用），无法唯一解析（如 H 封端把 sp2 端变成非手性 CH2）则跳过留无立体。
- `_radical_yl_from_sub`（`as_substituent.py:51`）处理**碳连接点**：`build_anchor_submol` 打 `*` 后经完整管线自由命名——L1 检测自由基（p41=1）、L2 选 radical 主基团、L4 锚定位次、L5 输出 `{locants}yl`；canonical SMILES 作缓存键与 `_name_mol` 共享 `CommonNameCache`。**缓存语义**：只缓存片段自身的自由基名（`copy.copy(hit)`，不经 `_canonical_result` 规范化、保留锚定链）——立体随宿主根分子变化，不能跨根共享，故 fresh 与 cache 命中都按当前根分子再过一次 R/S 校正。仅当 `parent_kind == "radical"` 时采用（理论必达，防御性判断）。
- `_fix_rs_with_real`（`as_substituent.py:21`）做 **radical 取代基 R/S 回根分子重算**：糖苷/醚类（如 O-C 糖-糖）异头位被 `*` dummy 顶替后，锚定子结构的 CIP 会随配基翻转，须回到完整根分子重算才正确（与 ChEBI 一致）。命中 meta 的 `parent_kind == "radical"` 时，把其 `parent_chain` 索引经 `block_root_order`（根上下文映射）映到根分子索引，用 `layer5.stereo` 的 `_cip_on_chain` 分别对锚定子结构与根分子取 R/S，不一致则以 `_with_rs` 重写 en/zh。
- **括号判定改由 meta 接口**（不再反编译名字，P-16.5.1.1）：词干是否复合由自由基命名的母体取代基数给出——`composite = meta.parent_substituent_count > 0`（namer 在 meta 注入母体取代基数），`need_paren = composite and en not in ("phenyl", *_SIMPLE_ALKOXY_NO_PAREN)`——复合前缀必括、未取代简单基免括。
- `_yl_from_sub`（`:90`）按连接点原子序数分派：碳 → `_radical_yl_from_sub`；非碳路径暂返回 None（留 H 封端 free-name + `free_to_yl` 的扩展位）。
- `name_as_substituent`（`:105`）为入口，`atoms` 强制 frozenset 后委托 `_yl_from_sub`；新增 `root_ctx` 参数贯穿到子路径（与 `_name_mol` 的 `root_ctx` 同为 `(root_mol, to_root)`）。
- yl 转换由 `tools/free_to_yl.free_to_yl` 完成（`as_substituent.py` 直接 `from namepredict.tools.free_to_yl import free_to_yl`，无 `layer3/yl_form.py`），处理官能团后缀到前缀的特殊转换：醇→烷氧基 (P-63.2.2)、硫醇→烷硫基 (P-63.2.1)、伯胺→烷氨基 (P-62.2)。

> 注：芳基/苯环侧链不再走独立 `_arene_yl_from_sub`——苯基等经 `_radical_yl_from_sub` 锚定自由基管线统一命名（苯 variant 保留名）。

> **源:** `src/namepredict/layer3/as_substituent.py`, `src/namepredict/layer3/submol_build.py`

### 氨基取代基 (amino_side)

`amino_side.py` 处理非胺母体上的氨基取代基，主入口 `_extract_aminos`（由 `substituent_extractor._extract_core_subs` 调用）：母体本身以胺为主官能团时跳过（`_principal_attachments` 检查 principal expression 的胺附着位）；否则遍历 `info["amines"]`，链上氨基直接 `_make_amino`，链外仲氨基经 `_sec_amino_off_chain` 组装为带括号的 `({name})amino` 形式。**酰胺母体的 N-烷基由 layer1 检测 + layer5 组装处理，不走 layer3**；胺母体的 N 端取代基（P-62.2 N- 前缀）由 `claimable_block` 的 `AMINE_N` slot 识别 + L5 `_N_PREFIX_KINDS` 组装（8584795）。`claimable_block._is_amine_n`（`claimable_block.py:68`）现排除**环员 N**（`atom.IsInRing()`）：饱和环内 N（哌啶/氮杂环等）是环杂原子、按环位次定位，不走 N- 前缀，避免 ring-N 取代基误用 N- 前缀。

> **源:** `src/namepredict/layer3/amino_side.py`

### 覆盖台账 (Coverage Ledger)

覆盖台账是 layer3 与上层质量控制的桥梁。`build_coverage_ledger` 接收分子、母体拥有的原子集合和已命名的取代基列表，输出 `CoverageLedger`：`owned_atoms`、`named_claims`、`gap`（遗漏）、`overlap`（冲突），`complete` = `not gap and not overlap`。在 `namer.py` 中，覆盖完整的候选优先被采用（Pass 1），仅当全部候选不完整时才回退到部分覆盖结果（Pass 2，标记 `fallback: no_coverage_gate`）。

> **源:** `src/namepredict/layer3/coverage.py`

### 取代基注册表 (anchored_table._REGISTRY)

`tools/anchored_table.py` 维护集中式取代基注册表 (`_REGISTRY`)，基础取代基与 IUPAC 2013 蓝皮书 P-29/P-57/P-61-P-68 条目统一入表。每条 `RetainedSubstituent`（L26）含保留英文名/中文名、系统名、IUPAC 推荐级别（PIN/GENERAL/NOT_RECOMMENDED）、`anchored` 锚定键、`paren`、`kind`。`resolve_name(key, name_mode)`（`:128`）是 `anchored_key` 的解析后端：general/pin 分派已注释停用，恒返回系统名（`systematic_en`/`systematic_zh`）。

> **源:** `src/namepredict/tools/anchored_table.py`

## 数据流图

```mermaid
flowchart TD
    subgraph L2["Layer2 输出"]
        PARENT["parent dict<br/>chain / owned_atoms / kind"]
        INFO["info dict<br/>mol / FG lists / ethers / amines"]
    end

    subgraph EXTRACT["substituent_extractor.py"]
        direction TB
        CORE["_extract_core_subs<br/>卤素 + OH + NH2 + 氧代"]
        ALKYL["_extract_alkyls_no_aryl<br/>anchored_table 查表烷基"]
        FULL["_with_full_atoms<br/>side_atoms 全连通补全"]
        CLAIM["extract_claimed_sides<br/>claimable_block → SubstituentNamer"]
    end

    subgraph NAMER["substituent_namer.py (二后端)"]
        direction TB
        RETAINED["RetainedBackend<br/>anchored_table 查表"]
        RECURSE["RecursiveBackend<br/>cut→free-name→yl"]
    end

    subgraph LEDGER["coverage.py"]
        BUILD["build_coverage_ledger<br/>owned_atoms + named_claims → gap / overlap"]
    end

    PARENT --> CORE
    INFO --> CORE
    PARENT --> ALKYL
    ALKYL --> FULL
    CORE --> FULL

    FULL --> CLAIM
    CLAIM -- "未覆盖 claim" --> RETAINED
    RETAINED -- "未命中" --> RECURSE

    RETAINED --> SUBST["list[dict]"]
    RECURSE --> SUBST
    CLAIM --> SUBST

    SUBST --> BUILD
```

```mermaid
flowchart LR
    subgraph CUT["cut → free-name → yl 管道 (as_substituent.py)"]
        BC["build_cut_submol<br/>诱导子分子 + H 封端"]
        YL["_yl_from_sub<br/>连接点类型分派"]
        RAD["_radical_yl_from_sub<br/>锚定 * → radical 主基团管线"]
        FY["tools.free_to_yl<br/>-yl 转换"]
    end

    subgraph ANCHOR["锚定查表 (tools/anchored_table.py)"]
        BA["build_anchor_submol<br/>attach 打 dummy *"]
        AT["_REGISTRY<br/>45 条 + _ANCHOR_INDEX"]
    end

    BC --> YL
    YL --> RAD
    RAD --> FY
    BA --> AT
    AT -->|"anchored_key"| REG["tools/anchored_table<br/>resolve_name"]
```

## 文件清单

| 文件 | 行数 | 描述 |
|------|------|------|
| `__init__.py` | 6 | 公开 API 导出：`extract_substituents` |
| `substituent_extractor.py` | 231 | **主提取器**。三段流水线：`_extract_core_subs` + `_extract_alkyls_no_aryl`（anchored 查表）+ `extract_claimed_sides`。含 `alkyl_alpha_key`（P-14.5 字母序排序键，剥前导位次/立体组 `(2S,3R)-`/`(E)-`/整体 `[` 后按实质词干排，L3/L5 共享）。 |
| `substituent_namer.py` | 109 | **命名引擎**。有序后端链：Retained(anchored) → Recursive（接收 `root_ctx`）。 |
| `as_substituent.py` | 114 | **cut→free-name→yl 管道**。`_yl_from_sub` 碳连接点 → `_radical_yl_from_sub`（锚定 * radical 管线），`_fix_rs_with_real` 回根分子校正 R/S，括号经 meta 接口判定；共享 CommonNameCache。 |
| `submol_build.py` | 152 | **子分子构建**。`build_cut_submol` / `build_anchor_submol` / `_carry_alkene_stereo`（E/Z 迁移）/ `CutSubmol`。 |
| `claim_extract.py` | 101 | **声明侧链补全**。`extract_claimed_sides` 遍历 `iter_claims`，对未覆盖 claim 调 `SubstituentNamer`（`root_ctx` 自 `info.get("root_ctx")` 注入）。 |
| `claimable_block.py` | 194 | `ClaimedBlock` / `SideSlot`(CHAIN_C/RING_C/AMIDE_N/**AMINE_N**/ETHER_O/OTHER) / `iter_claims`。`_is_amine_n` 排除环员 N。 |
| `amino_side.py` | 42 | 氨基取代基：伯氨基 + 仲氨基（`_extract_aminos`）。 |
| `coverage.py` | 71 | **覆盖台账**。`build_coverage_ledger` 计算 gap/overlap。 |

> 备注：layer3 只有以上 9 个模块。`side_facts.py`/`aryl_sub.py`/`yl_form.py` 均不存在——`carbon_neighbors` 位于 `tools/chain.py`，yl 转换在 `tools/free_to_yl.py`。

### tools/ 层无关工具（layer3 消费）

| 文件 | 行数 | 描述 |
|------|------|------|
| `tools/anchored_table.py` | 201 | **锚定 canonical-SMILES 查表 + 取代基注册表**。`_REGISTRY`（45 条统一 registry）+ `_ANCHOR_INDEX` 锚定反查 + `anchored_entry`/`anchored_lookup`/`pick_root`/`resolve_name`/`anchored_whole_mol`。核心机制。 |
| `tools/block_cut.py` | 82 | 母体边界块切割：`side_atoms`（全连通分量）/`cut_block`/`side_roots`。 |
| `tools/chain.py` | 43 | 碳链行走原语 `_carbon_neighbors`/`_longest_from`/`carbon_neighbors`，L2/L3 共享。 |
| `tools/free_to_yl.py` | 194 | **-yl 转换**（layer-agnostic）。 |

> 备注：`tools/alkoxy_side.py` 不存在（酯 O 侧拓扑由 L2 principal_expression + L5 `join_ester_name` 承担）。

## 对外接口

### 主入口

```python
def extract_substituents(
    info: dict,
    parent: dict,
    *,
    name_mode: str = "general",
    cache: CommonNameCache | None = None,
) -> list[dict]:
```

- **info**: layer1 分析结果，必须包含 `"mol"` 键（RDKit Mol 对象）以及各类 FG 列表
- **parent**: layer2 母体选择结果，必须包含 `"chain"` 和 `"owned_atoms"` 键
- **name_mode**: `"general"` 使用保留/通用名；`"pin"` 使用 IUPAC 首选名
- **cache**: 递归子结构命名的共享缓存（`CommonNameCache`），加速跨 cut 命名
- **返回**: 取代基 dict 列表，共同键包括 `kind`, `en`, `zh`, `attach_idx`, `atoms`, `paren`

### 内部命名类

```python
class SubstituentNamer:
    def __init__(self, backends=None, *, name_mode="general", cache=None)
    def name(self, mol, claim: ClaimedBlock, *, depth=0) -> SubstituentName | None
```

`SubstituentName` 数据类包含 `claim`, `en`, `zh`, `requires_parentheses`, `backend` 字段。

### 覆盖台账

```python
def build_coverage_ledger(mol, *, owned_atoms, names) -> CoverageLedger
```

`CoverageLedger` 属性: `owned_atoms`, `named_claims`, `gap`, `overlap`, `complete`.

## 相关页面

- [[architecture/layer2-parent-selector]] — Layer2 母体选择器，提供 `parent` dict（含 `owned_atoms` 和 `chain`）；L2 与 L3 互不调用，共享仅经 `tools/`
- [[architecture/layer1-analyzer]] — Layer1 官能团分析器，提供 `info` dict
- [[architecture/layer4-numbering]] — Layer4 编号引擎，消费 layer3 输出的取代基列表
- [[architecture/layer5-name-assembly]] — Layer5 名称组装，最终拼接母体名和取代基前缀（-yl 转换在 `tools/free_to_yl`）
- [[architecture/overview]] — 系统架构概览，6 层流水线总览
- [[concepts/bilingual-naming]] — 中英双语命名约定
- [[concepts/atom-ownership]] — ClaimedBlock / owned_atoms / CoverageLedger 原子归属模型
