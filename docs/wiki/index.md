# NamePredict Wiki

> 最后更新: 2026-09-13 | 源文件: 64 `.py` / 8,880 行（`src/namepredict` 全量，含 `tools/`、`constants.py`、`namer.py`，排除 `__pycache__`）| Wiki 页面: 15

## 项目概述

**NamePredict** 是一个基于规则的 SMILES -> IUPAC 双语命名引擎。输入 SMILES 字符串，输出中英双语 IUPAC 名称。

版本: `0.1.0` | 核心依赖: RDKit | Python 3.10+

## 架构一览

```
SMILES 输入
    |
[layer0] 预处理 -> 分子解析、立体指派、互变异构/电荷归一、盐解离
    |
[layer1] 官能团分析 -> FG_SPECS 统一检测 -> fg_inventory + 环系事实 -> info dict
    |
[layer2] 母体选择 -> P-44 规则管线选主官能团/骨架、kind 正交化 <- 核心
    |
[layer3] 取代基提取 -> claim 枚举 + 锚定查表 + 递归命名、覆盖台账
    |
[layer4] 编号与定位符 -> P-14.4 三层分派编号引擎、位次分配
    |
[layer5] 名称组装 -> 中英双语拼接 <- 输出层
    |
{en: "ethanol", zh: "乙醇"}
```

## 架构文档

| 页面 | 内容 |
|------|------|
| [[architecture/overview]] | 架构总览：6 层流水线、数据流图、跨层设计模式、namer.py 协调器与 P-44.1.1/P-45.2.2 并列候选裁决 |
| [[architecture/layer0-preprocessor]] | Layer0 预处理：SMILES 解析 + 立体初步指派 + 酰胺烯醇互变异构归一化（电荷重定位后二次归一）+ 酸性质子收敛（charge.py 按 `ACID_CENTERS` 查表识别成酸中心）、盐检测与解离、盐元数据注入 |
| [[architecture/layer1-analyzer]] | Layer1 官能团分析器：`FG_SPECS`（14 条 `FgSpec`）驱动的统一检测接口（条目形态 `{center_idx, surr_idx}`）、`fg_inventory` 单一出口、P-41 仲裁 `_arbitrate_parts`（组合羰基压制 + carboxy/cyano 降级叶）、环系拓扑与桥/螺/稠合识别 |
| [[architecture/layer2-parent-selector]] | Layer2 母体选择器：P-44 规则驱动主链管线 + P-45.2.1 前缀取代基打平（`select_parent` 返回并列组 `list[dict]`）、kind 正交化（`ring_scaffold._TEMPLATES` 83 条模板 / 36 条登记固定编号，唯一来源）、稠环拆解（fused_system）、杂原子锚点词干、owned_atoms 归属 |
| [[architecture/layer3-substituents]] | Layer3 取代基提取器：claim 枚举为唯一提取路径（`iter_claims`/`ClaimedBlock`/`SideSlot`）、二后端命名（anchored 查表 `_REGISTRY` 70 条 + 递归）、覆盖台账 |
| [[architecture/layer4-numbering]] | Layer4 编号：P-14.4 三层分派（fixed → fused → 普用候选枚举）+ 稠环编号（fused_orientation/fused_numbering/ring_geometry）、指示氢与加氢前缀（`indicated_h_locants` + `hydro_prefix`）、FG 位次表 `_FG_LOCANTS`（由 `FG_SPECS` 投影）、locant 省略（FG/烯/炔）、并列候选裁决键 |
| [[architecture/layer5-name-assembly]] | Layer5 名称组装：组装流水线、`chain_engine._KIND_TABLE` 链引擎（14 entry 含 `sulfonic`、逐卤素 `acyl_halide`）、`_chain_enyne` 统一烯/炔段、环外主基系统名与稠合名组装、`free_to_yl`（单核母体氢化物 → -yl）、N- 前缀、立体化学（E/Z + R/S 含稠环字母位次） |

## 核心概念

| 页面 | 内容 |
|------|------|
| [[concepts/functional-group-priority]] | 官能团优先级体系：`FgSpec.p41`/`path` 为单一权威（`PrincipalPriority` 取最小）、L1 `_arbitrate_parts` → L2 骨架级/候选级排序 → L5 `_KIND_TABLE` 三层接力 |
| [[concepts/atom-ownership]] | 原子归属跟踪：L1 出 FG 锚点/特征原子事实、L2 `_kind_fg_atoms` 做几何边界合成（`owned_atoms`）、L3 消耗边界；`ClaimedBlock`/`SideSlot`、gap/overlap 覆盖完整性 |
| [[concepts/bilingual-naming]] | 中英双语命名约定：词表集中在 `constants.py`、词干派生在 `stems.py`、酯类/盐类语序反转、`zh_stem`/`zh_num` 中文数字、递归命名的双语缓存与宿主根传播 |

## 参考

| 页面 | 内容 |
|------|------|
| [[reference/core-data-contracts]] | 核心数据合约：info dict（11 键 + root_ctx/salt）、parent dict、numbered dict、NameResult、ClaimedBlock、SubstituentName、CoverageLedger、subst dict -- 数据结构完整字段表与数据流图 |

## 扩展指南

| 页面 | 内容 |
|------|------|
| [[guides/adding-new-functional-group]] | 新增官能团指南：以硫醇 (thiol) 为工作示例，覆盖 L1 检测与 `FgSpec` 登记 → L2 母体接线 → L5 `_KIND_TABLE` 注册，含派生表清单、文件改动清单与常见陷阱 |
| [[guides/adding-new-ring-system]] | 新增环系指南：以吡啶 (pyridine) 为工作示例，覆盖 L1 环检测 → `ring_scaffold._TEMPLATES` 注册 → 固定编号 → L5 命名 → FG-环组合，含文件改动清单 |

## 外部资源

- `docs/ARCHITECTURE.md` -- 原始架构设计文档
- `docs/RULES.md` -- IUPAC 规则实现状态跟踪
- `docs/iupac/` -- IUPAC 蓝皮书原文（中英两版）
