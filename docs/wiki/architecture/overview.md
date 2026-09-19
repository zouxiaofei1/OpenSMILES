# 架构总览 (Architecture Overview)

> **源码规模:** src/namepredict 共 61 个 `.py` / 9,278 行

## 项目定位

NamePredict 是规则驱动的 SMILES → IUPAC 双语命名引擎：输入 SMILES 字符串，输出 `{en, zh}` 中英双语系统名，封装为 `types.NameResult`（`en` / `zh` / `success` / `source` / `time_ms` / `meta`）。判据取自 IUPAC 建议（P-14 编号、P-41 主基团、P-44 母体选择、P-45 排序、P-6x 各类官能团等），原则是表驱动优先、代码驱动兜底：能落表的判据（SMARTS、词干、优先级、词尾、保留名）一律进表，代码只做分派与拼接。对外入口为 `namer.SMILESNNamer.name`。

## 六层流水线

```
SMILES
  │
  ▼  L0  layer0/   预处理             5 文件 /   355 行   → mol + salt
  │
  ▼  L1  layer1/   官能团分析          6 文件 /   733 行   → info{mol, fg_inventory, …}
  │
  ▼  L2  layer2/   母体选择           11 文件 / 2,418 行   → parent{chain, kind, owned_atoms, …}
  │
  ▼  L3  layer3/   取代基提取          7 文件 /   520 行   → subst[{kind, en, zh, attach_idx, paren}]
  │
  ▼  L4  layer4/   编号与位次         11 文件 / 1,719 行   → numbered{parent, locants, …}
  │
  ▼  L5  layer5/   名称组装            9 文件 / 2,463 行
  │
  ▼  NameResult{en, zh, success, meta}
```

- **L0 预处理 `layer0/`**：`preprocess` 清洗并解析 SMILES；`tautomer` 把次要互变体归一为优势式；`charge` 重定位酸性质子；`salt` 拆出有机母体与盐元数据。
- **L1 官能团分析 `layer1/`**：`fg_local_smarts.FG_SMARTS` 表驱动局部检测，`analyzer` 补非局部后处理（含氧酸 `oxo_kind` 归一、P-41 降级仲裁），产出 `FunctionalGroupInventory`。
- **L2 母体选择 `layer2/`**：P-44 逐级筛骨架 + P-41 选主基团 + P-45.2 按前缀数排序 + `kind` 正交化，并定下 scaffold 身份（保留模板 / 稠环树 / 桥环节点）。
- **L3 取代基提取 `layer3/`**：`claimable_block.iter_claims` 在归属边界内枚举 claim，`substituent_namer` 先查锚定表再递归命名，`coverage` 出覆盖台账。
- **L4 编号与位次 `layer4/`**：`numbering_engine.orient_numbering` 按 P-14.4 三层分派定序（桥环另走 P-23 裁决分支），`numbering` 计算位次、指示氢与加氢前缀。
- **L5 名称组装 `layer5/`**：注入母体词干（稠环 / 桥环 / 生成式 / C1 保留名）后由链引擎取名，再拼前缀、O-侧臂、桥、立体与盐后缀，输出 `NameResult`。

## 协调器 namer.py

`SMILESNNamer` 持有 `tools.common_names.CommonNameCache`（默认 20000 条），`name` 先查缓存，命中即返回，未命中走 `_pipeline` 并只回写成功结果；`_canonical_result` 把 `meta.parent_chain` 转成 RDKit 规范原子序，使缓存命中的元数据与首算一致。`_name_mol` 编排 L0→L5：盐解离 → `analyze` → `select_parent` → `extract_substituents` → `number` → `assemble`。

盐后缀在组装之后由 `_apply_salt_suffix` 经 `stems.join_metal_salt_names` 补齐；含氧酸中心母体（`OXO_CENTER_KINDS`）自带盐组装，故跳过。P-44.1.1 并列组的候选逐个跑一遍 L3–L5，裁决键为 `(p44_1_1_key, p45_2_2_key)`：P-44.1.1 未决时取 P-45.2.2 前缀位次最小者。失败经 `_fail` 返回，`meta.reason` 记录原因（`parse` / `no_assemblable_candidate`），覆盖门控回退标 `meta.fallback`。诊断用的链路元数据一并落在 `meta`：`parent_chain` / `parent_kind` / `parent_labels` / `bridge_self_enclosed` / `parent_substituent_count` / `salt`。

## 跨层设计模式

| 模式 | 承载点 | 要点 |
| --- | --- | --- |
| 表驱动单一事实来源 | `layer1/fg_registry.FG_SPECS`；`layer5/chain_engine._KIND_TABLE`；`constants.py` | FG 类别、kind 词尾、跨层词表各只有一处登记，下游表由之派生 |
| kind 正交化 | L2 `principal` / `kind_registry` → L4 / L5 | 主基团类别压成 `kind` 字符串，L4/L5 只认 `kind`、不认 FG 枚举 |
| 含氧酸合一管线 | L1 `analyzer` → `fg_registry.oxoacid_is_acid` → L5 `chain_engine._OXO_TAIL` | P 与 S 中心共用一条链，含 S–O–S 多硫酸链 |
| 多环拆解择一 | L2 `fused_system` / `bridged_system` → L4 / L5 | 稠合与桥环互斥：桥环仅在稠合名不适用或拼不出词干时接管 `scaffold_id` |
| 原子归属 | L2 `parent_select.finalize_parent_ownership` → L3 `claimable_block` / `coverage` | `owned_atoms` 划定 L3 枚举边界，台账以 `gap = 重原子 − owned_atoms` 兜底 |
| P-41 优先级接力 | `FgSpec.p41` → L2 `principal._effective_priority` / `parent_skeleton` → L5 | 优先级声明一次，骨架筛选与出词尾共用 |
| 前缀排序键统一 | `tools/re.alpha_order_key` | P-14.5 字母数字序唯一实现，编号与组装两侧共用 |
| 围栏由产生方定形 | L3 → L5 `assembler_prefixes` | L3 决定取代基是否自带围栏，L5 只做渲染与必要的整体围栏 |

- **`FG_SPECS` 与 `_KIND_TABLE`**：`FgSpec` 一条记录同时描述 FG 枚举值、L1 检测列表键、锚点字段、P-41 等级、P-43 路径与表达类型；L2 的 `PRINCIPAL_REGISTRY`、`functional_group_inventory._FG_KEYS` 与 `_ANCHOR_KEYS` 全部由它派生，扩展官能团只需加一条表项。命名端同理：`_KIND_TABLE` 每个 `kind` 一条 `_Chain` 规格（词尾、保留名 variant、位次省略规则）。
- **含氧酸合一管线**：`FG_SMARTS` 按「元素 + 双键氧数 + 臂型」逐条登记 P/S 中心，`analyzer` 把臂角色归一为 `oxo_kind`（`phosphate` / `phosphonate` / `sulfate` / `sulfonic` / 磺酸酯 / `sulfonamide` / `sulfonyl_chloride`）；`oxoacid_is_acid` 判是否按 P-41 类 7 酸式参与主基团竞争，否则落类 9 酯或类 11 磺酰胺，`_OXO_TAIL` 按 `(oxo_kind, n_oh)` 出词尾。
- **前缀排序键统一**：`tools/re.alpha_order_key` 循环剥除 N- 前缀、斜体前缀、位次、立体描述符与外层括号，先比字母序列、再比首字母前位次；`numbering_engine` 的 P-14.4(f) 平局与 L5 前缀排序共用同一键。
- **围栏由产生方定形**：`SubstituentName.requires_parentheses` 随取代基一路传到 L5，L5 仅在渲染前缀时套用（`assembler_prefixes._enclose`、`oxo_arm_fence`），并对整名补必要围栏，不在组装侧二次拆分。
- **多环拆解择一**：一个环系要么走 P-25 稠合、要么走 P-23 桥环，不并存。`principal_expression._ring_scaffold_and_nodes` 一次算好两者，桥环命中即用 `BridgedNode.scaffold_identity()` 覆盖 scaffold 身份；桥环本身只交出**并列最优候选**（`bridged_nodes`），编号交给 L4 按 P-23.3.2 与 P-14.4 裁决后写回 `bridged_node`，L5 只读 `descriptor` 与上标位次对。

## 文件组织

```
src/namepredict/                （各层 __init__.py 只重导出该层入口）
├── __init__.py                 包装入口，安装 rdkit 迭代器加速
├── namer.py                    顶层协调器：缓存 + L0–L5 编排
├── constants.py                跨层常量：原子序数、倍数前缀、词表
├── types.py                    NameResult 数据类
├── layer0/
│   ├── preprocessor.py         preprocess：清洗解析 + 互变归一链
│   ├── tautomer.py             酮式/内酰胺/硫酮/胺式归一
│   ├── charge.py               酸性质子重定位与电荷归一
│   └── salt.py                 盐解离与 salt 元数据
├── layer1/
│   ├── fg_local_smarts.py      FG_SMARTS 局部检测表（中心原子）
│   ├── analyzer.py             非局部后处理、oxo_kind 归一、降级仲裁
│   ├── fg_registry.py          FG_SPECS：FG 跨层元数据唯一注册处
│   ├── functional_group_inventory.py  类型化清单与 occurrence
│   └── ring_systems.py         环系拓扑（L1/L4/L5 共用）
├── layer2/
│   ├── parent_select.py        select_parent 入口 + owned_atoms 固化
│   ├── parent_skeleton.py      P-44 拓扑优先骨架枚举与逐级窄化
│   ├── principal.py            P-41 / P-43 主基团选择
│   ├── principal_expression.py 为选中骨架标注主基团表达式 facts
│   ├── ring_expression_policy.py  环骨架类型化表达的能力策略
│   ├── chain_walk.py           脂肪族碳链行走
│   ├── fused_system.py         P-25.3.2.4 稠环组分树
│   ├── bridged_system.py       P-23 扩展 von Baeyer 桥环拆解
│   ├── ring_scaffold.py        骨架规格 + 保留 SMILES 模板表
│   └── kind_registry.py        kind 词干注册与词干打包
├── layer3/
│   ├── substituent_extractor.py  extract_substituents 入口
│   ├── claimable_block.py      claim 枚举、slot 推导与块切割
│   ├── substituent_namer.py    命名后端：锚定表 → 递归
│   ├── coverage.py             重原子覆盖台账
│   ├── submol_build.py         诱导切割子分子（H 封端）
│   └── as_substituent.py       子分子 → free-name → -yl
├── layer4/
│   ├── numbering.py            number 入口 + 指示氢/加氢前缀
│   ├── numbering_engine.py     P-14.4 定向编号引擎
│   ├── locant_calc.py          附着原子 → 链上位次
│   ├── candidate_keys.py       P-44.1.1 / P-45.2.2 候选比较键
│   ├── fused_numbering.py      稠环编号：外周序、稠合位次、CIP 破局
│   ├── bridged_numbering.py    P-23 桥环编号裁决
│   ├── fused_orientation.py    稠环优选取向
│   ├── ring_geometry.py        环几何坐标与变形阈值
│   ├── indicated_hydrogen.py   指示氢位次
│   └── omit_locants.py         FG 位次省略规则
├── layer5/
│   ├── assembler.py            assemble 入口：取名 → 前缀 → 立体 → 盐
│   ├── chain_engine.py         链式词干引擎与 _KIND_TABLE
│   ├── assembler_prefixes.py    前缀分组与双语渲染、围栏
│   ├── stems.py                碳数词干表与金属盐后缀
│   ├── stereo.py               E/Z 与 CIP R/S 前缀
│   ├── fused_namer.py          稠合名组装
│   ├── bridged_namer.py        P-23 桥环母体名（环数词头 + 描述符）
│   └── skeleton_replacement.py 骨架置换 'a' 前缀
└── tools/
    ├── anchored_table.py       取代基锚定 canonical-SMILES 查表
    ├── re.py                   名称文本辅助与 alpha_order_key
    ├── common_names.py         SMILES → NameResult 内存缓存
    ├── memo.py                 单次命名内的中间结果记忆
    ├── block_cut.py            母体边界块切割
    ├── rdkit_fast.py           RDKit 迭代器加速补丁
    └── chain.py                开链最长路径搜索
```

## 页内导航

- [[architecture/layer0-preprocessor]] — 预处理、互变归一、质子重定位与盐解离
- [[architecture/layer1-analyzer]] — FG 检测表与 `FunctionalGroupInventory`
- [[architecture/layer2-parent-selector]] — P-44 骨架筛选、P-41 主基团与 P-45.2 排序
- [[architecture/layer3-substituents]] — claim 枚举、锚定查表与递归命名
- [[architecture/layer4-numbering]] — P-14.4 编号引擎与位次表达
- [[architecture/layer5-name-assembly]] — 链引擎、前缀组装与立体化学
- [[concepts/functional-group-priority]] — 官能团优先级体系
- [[concepts/atom-ownership]] — `owned_atoms` 原子归属机制
- [[concepts/bilingual-naming]] — 中英双语命名策略
- [[reference/core-data-contracts]] — 核心数据合约
- [[guides/adding-new-functional-group]] / [[guides/adding-new-ring-system]] — 扩展指南
