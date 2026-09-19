# NamePredict Wiki

> **源码规模:** `src/namepredict` 共 61 个 `.py` / 9,278 行（含 `tools/`、`constants.py`、`namer.py`，排除 `__pycache__`）
> **文档页面:** architecture 7 篇 + concepts 3 篇 + guides 2 篇 + reference 1 篇

## 项目概述

**NamePredict** 是规则驱动的 SMILES → IUPAC 双语命名引擎：输入 SMILES 字符串，输出中英双语系统名，封装为 `types.NameResult`（`en` / `zh` / `success` / `source` / `time_ms` / `meta`）。对外入口为 `namer.SMILESNNamer.name`。

判据取自 IUPAC 建议（P-14 编号、P-41 主基团、P-44 母体选择、P-45 排序、P-6x 各类官能团等）。设计原则是**表驱动优先、代码驱动兜底**：能落表的判据（SMARTS、词干、优先级、词尾、保留名）一律进表，代码只负责分派与拼接。

核心依赖 RDKit（`Mol` 为层间传递的分子表示）；Python 3.11+。

## 架构一览

```
SMILES 输入
    │
    ▼  L0  layer0/   预处理            → mol + salt 元数据
    │      解析消毒 · 互变异构归一 · 酸性质子重定位 · 盐解离
    │
    ▼  L1  layer1/   官能团分析         → info{mol, fg_inventory, salt, …}
    │      FG_SMARTS 表驱动局部检测 + analyzer 非局部后处理
    │
    ▼  L2  layer2/   母体选择           → parent{chain, kind, owned_atoms, …}
    │      P-44 骨架筛选 · P-41 主基团 · P-45.2 排序 · kind 正交化
    │
    ▼  L3  layer3/   取代基提取         → subst[{kind, en, zh, attach_idx, paren, …}]
    │      claim 枚举 · 锚定查表/递归命名 · 覆盖台账
    │
    ▼  L4  layer4/   编号与位次         → numbered{parent, substituents, numbering_scaffold}
    │      P-14.4 三层分派编号引擎 · 指示氢与加氢前缀
    │
    ▼  L5  layer5/   名称组装           → NameResult{en, zh, success, meta}
    │      链引擎取名 · 前缀/桥/立体/盐拼接
    │
    ▼  {en: "ethanol", zh: "乙醇"}
```

## 架构文档

| 页面 | 内容 |
|------|------|
| [[architecture/overview]] | 架构总览：项目定位、六层流水线与逐层职责、`namer.py` 协调器、跨层设计模式表、61 文件目录树与页内导航 |
| [[architecture/layer0-preprocessor]] | Layer0 预处理：`preprocess` 六步链、`tautomer` 双向互变异构归一化（烯醇/烯硫醇式 → 酮式、环外亚胺式 → 胺式）、`charge` 按 `ACID_CENTERS` 查表重定位酸性质子（字典键序即酸强度序）、`salt` 盐解离与双语盐元数据 |
| [[architecture/layer1-analyzer]] | Layer1 官能团分析器：`fg_local_smarts.FG_SMARTS` 表驱动局部检测、`analyzer` 非局部后处理（含氧酸 `oxo_kind` 归一为重点，P/S 中心含 S–O–S 多硫酸链）、`fg_registry.FG_SPECS` 跨层元数据与 P-41 酸式升类判据 `oxoacid_is_acid`、`functional_group_inventory` 单一出口 |
| [[architecture/layer2-parent-selector]] | Layer2 母体选择器：`parent_select` 层唯一门面（P-44 编排 / 候选收集 / 原子归属 / P-45.2 排序）、`parent_skeleton` 筛选谓词、`principal` 主基团选择与 `_effective_priority`、`principal_expression` kind 正交化与含氧酸字段、`ring_scaffold` 模板表与诱导覆盖校验、`fused_system` 稠环树与 `bridged_system` P-23 桥环拆解 |
| [[architecture/layer3-substituents]] | Layer3 取代基提取器：claim 枚举为唯一提取路径（`iter_claims` / `ClaimedBlock` / `SideSlot`）、锚定查表 + 递归命名双后端、`_alkoxycarbonyl` 回退、围栏白名单与双原子桥例外、`side_z` 驱动的 O-侧臂标记、覆盖台账 |
| [[architecture/layer4-numbering]] | Layer4 编号与位次：P-14.4 三层分派（桥环裁决 → fixed → fused → 通用候选枚举）、`narrow`/`narrow_by_senior` 收窄原语、`bridged_numbering` 桥环编号裁决、稠环编号与定向、`indicated_hydrogen` 护栏与指示氢收窄、位次省略与 `alpha_order_key` 字母数字序 |
| [[architecture/layer5-name-assembly]] | Layer5 名称组装：`_ensure_parent_stem` 词干注入分派（桥环 `bridged_namer` / 生成式 / 稠环 / C1 保留名）、`skeleton_replacement` 骨架置换 'a' 前缀、`chain_engine._KIND_TABLE` 链引擎、`_OXO_TAIL` 含氧酸词尾表、`join_kind_name` 三条拼接路径、酯与硫代酯语序、环外主基命名、前缀组装与 `oxo_arm_fence` 围栏、指示氢前缀、立体化学 |

## 核心概念

| 页面 | 内容 |
|------|------|
| [[concepts/functional-group-priority]] | 官能团优先级体系：`FgSpec.p41`/`p43_path` 为单一权威（`PrincipalPriority` 取最小）、类 7 的含氧酸漂移判据、L1 声明 → L2 消费 → L5 出词尾的三层接力、压制与降级集合、P-45.2 排序与缩合含氧酸（P/S）次序 |
| [[concepts/atom-ownership]] | 原子归属跟踪：L1 出 FG 锚点与特征原子事实、L2 `finalize_parent_ownership` 合成 `owned_atoms` 几何边界、L3 在边界内消耗；含氧酸的「特征原子全归主基团」与中心母体只留 FG 原子两条规则、覆盖完整性诊断 |
| [[concepts/bilingual-naming]] | 中英双语命名：词表分层（`constants.py` → `stems.py` → `_KIND_TABLE`）、kind → 词尾分派、`_OXO_TAIL` 含氧酸词尾表、桥环 von Baeyer 名与骨架置换前缀的双语构形、C1 保留名族与 N 位次口径、酯/硫代酯/盐/中心母体的语序反转、中文侧特有处理、围栏与 `alpha_order_key` 排序 |

## 参考

| 页面 | 内容 |
|------|------|
| [[reference/core-data-contracts]] | 核心数据合约：info dict、FG inventory 与含氧酸 payload、parent dict（含 `bridged_node`/`bridged_nodes`/`stem_bare_*`/`locant_kind`/`subs_consumed`）、`BridgedNode`/`BridgeSegment`、numbered dict、`ClaimedBlock`/`SideSlot`、subst dict、`NameResult` 的字段表与端到端数据流图 |

## 扩展指南

| 页面 | 内容 |
|------|------|
| [[guides/adding-new-functional-group]] | 新增官能团指南：以碳锚定磺酸族为例走「检测 → 归一 → 注册 → 母体化 → 命名」五个落点（`FG_SMARTS` → `analyzer._oxo_kind`（含 S–O–S 多硫酸链）→ `FG_SPECS` → `_WHOLE_FG_ATOMS`/`OXO_CENTER_KINDS` → `_KIND_TABLE`/`_OXO_TAIL`/C1 保留名），附派生自动/手动对照表与常见陷阱 |
| [[guides/adding-new-ring-system]] | 新增环系指南：以 thiane 为例覆盖 L1 环检测 → `ring_scaffold._TEMPLATES` 注册（含 `fused_prefix` 的 mancude 语义、`_FUSION_CARBOCYCLES` 的 benzo 特例与诱导覆盖/饱和单环两重匹配校验）→ L4 固定编号与指示氢 → L5 命名 → 与 FG 组合；并给出免注册的两条路径（P-23 桥环与 >10 元生成式词干），附改动清单与陷阱 |

## 外部资源

- `docs/iupac/` — IUPAC 蓝皮书原文（中英两版）
  - `docs/iupac/cn_translated/` — Blue Book 中文译本（规则为英文原文）
  - `docs/iupac/cn/` — 中国化学会《有机化合物命名原则》（中文规则）
- `tools/BlueBookV2.pdf` — 带图原文
