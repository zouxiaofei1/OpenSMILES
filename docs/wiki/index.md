# NamePredict Wiki

> 最后更新: 2026-08-11 | 源文件: 111 `.py` / 12,374 行 | Wiki 页面: 15

## 项目概述

**NamePredict** 是一个基于规则的 SMILES -> IUPAC 双语命名引擎。输入 SMILES 字符串，输出中英双语 IUPAC 名称。

版本: `0.1.0` | 核心依赖: RDKit | Python 3.10+

## 架构一览

```
SMILES 输入
    |
[layer0] 预处理 -> 分子解析、盐检测
    |
[layer1] 官能团分析 -> 识别 FG、环指纹  -> info dict
    |
[layer2] 母体选择 -> P-44 规则管线选主官能团/骨架、纯烃兜底、保留名 <- 核心
    |
[layer3] 取代基提取 -> anchored 查表 + claim 补全、覆盖台账
    |
[layer4] 编号与定位符 -> 位次分配、方向决策
    |
[layer5] 名称组装 -> 中英双语拼接 <- 输出层
    |
{en: "ethanol", zh: "乙醇"}
```

## 架构文档

| 页面 | 内容 |
|------|------|
| [[architecture/overview]] | 架构总览：6 层流水线、数据流图、跨层设计模式、namer.py 协调器 |
| [[architecture/layer0-preprocessor]] | Layer0 预处理：SMILES 解析、盐检测与解离、盐元数据注入 |
| [[architecture/layer1-analyzer]] | Layer1 官能团分析器：34 种 FG 检测、环系拓扑、排他性优先级 |
| [[architecture/layer2-parent-selector]] | Layer2 母体选择器：P-44 规则驱动主链管线、纯烃表达、苯环保留名、kind_registry、评分与门控 |
| [[architecture/layer3-substituents]] | Layer3 取代基提取器：anchored 查表主流程 + claim 补全、SubstituentNamer 有序后端（retained/rooted-tree/recursive）、侧链拓扑事实层 (side_facts)、覆盖台账 |
| [[architecture/layer4-numbering]] | Layer4 编号：链/环方向决策、约束驱动编号引擎、omit_locants 规则 |
| [[architecture/layer5-name-assembly]] | Layer5 名称组装：5 步流水线、双语词干表、FG 专属命名模块、立体化学 |

## 核心概念

| 页面 | 内容 |
|------|------|
| [[concepts/functional-group-priority]] | 官能团优先级体系：P-41 降序排列 (fg_rank 0-13)、kind_registry 单一权威注册中心、互斥检查规则、L1/L2/L5 三层 FG 生命周期 |
| [[concepts/atom-ownership]] | 原子归属跟踪：owned_atoms 计算 (链骨架+FG 杂原子)、ClaimedBlock 桥接类型与 SideSlot 枚举、gap/overlap 覆盖完整性、14 种 FG 的异原子归属规则 |
| [[concepts/bilingual-naming]] | 中英双语命名约定：en/zh 元组惯例、词干表 (stems.py) 作为单一权威、酯类/盐类的语序反转、zh_stem 转换、zh_num 中文数字生成 |

## 参考

| 页面 | 内容 |
|------|------|
| [[reference/core-data-contracts]] | 核心数据合约：info dict、parent dict、numbered dict、NameResult、ClaimedBlock、SubstituentName、CoverageLedger、subst dict -- 8 种数据结构的完整字段表与数据流图 |

## 扩展指南

| 页面 | 内容 |
|------|------|
| [[guides/adding-new-functional-group]] | 新增官能团指南：以亚砜 (sulfoxide) 为工作示例，跨越 L1-L5 的完整注册模式 (检测->母体选择->取代基->编号->组装)，含文件改动清单和常见陷阱 |
| [[guides/adding-new-ring-system]] | 新增环系指南：以吡啶 (pyridine) 为工作示例，跨越 L1-L5 的完整流程 (环检测->ScaffoldSpec 注册->编号骨架->命名->FG-环组合)，含文件改动清单 |

## 外部资源

- `docs/ARCHITECTURE.md` -- 原始架构设计文档
- `docs/RULES.md` -- IUPAC 规则实现状态跟踪
- `docs/iupac/` -- IUPAC 蓝皮书原文（中英两版）
