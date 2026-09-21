# 架构总览 (Architecture Overview)

> **源码规模:** src/namepredict 共 62 个 `.py` / 10,581 行（含各层 `__init__.py`）

## 项目定位

NamePredict 是规则驱动的 SMILES → IUPAC 双语命名引擎：输入 SMILES 字符串，输出 `{en, zh}` 双语系统名，封装为 `types.NameResult`（`en` / `zh` / `success` / `source` / `time_ms` / `meta`）。判据取自 IUPAC 建议（P-14 编号、P-23 桥环、P-24 螺环、P-25 稠环、P-41 主基团、P-44 母体选择、P-45 排序、P-6x 各类官能团、P-71 盐、P-73 阳离子等），原则是表驱动优先、代码驱动兜底：能落表的判据（SMARTS、词干、优先级、词尾、保留名、锚定键）一律进表，代码只做分派与拼接。对外入口为 `namer.SMILESNNamer.name`。

## 六层流水线

```
SMILES
  │
  ▼  L0  layer0/   预处理与解盐      5 文件 /   398 行  → mol + salt
  │
  ▼  L1  layer1/   官能团分析        6 文件 /   793 行  → info{mol, fg_inventory, rings, …}
  │
  ▼  L2  layer2/   母体选择         12 文件 / 2,957 行  → 候选组[{chain, kind, owned_atoms, scaffold_id, …}]
  │
  ▼  L3  layer3/   取代基提取        7 文件 /   525 行  → subst[{kind, en, zh, attach_idx, requires_parentheses}]
  │
  ▼  L4  layer4/   编号与位次       10 文件 / 1,942 行  → numbered{parent, locants, ene/yne, hydro, …}
  │
  ▼  L5  layer5/   名称组装         10 文件 / 2,893 行
  │
  ▼  NameResult{en, zh, success, meta}

顶层 namer.py / constants.py / types.py 4 文件 / 551 行；tools/ 8 文件 / 522 行，为各层共用的拓扑与文本工具。
```

- **L0 预处理 `layer0/`**：`preprocess` 清洗字符串并解析（`_strip_isotopes` 去掉同位素标记，否则锚定键失配）；互变链为 `normalize_amide_tautomer` → `normalize_acid_charge` → 再跑一次 `normalize_amide_tautomer`（电荷重定位后才会出现可归一的酰胺烯醇位）。`dissociate_salt` 按 `METAL_ION_EN` / `HALIDE_ZH` / `HALIDE_HX_EN` 三类反离子拆出有机母体与盐元数据。
- **L1 官能团分析 `layer1/`**：`fg_local_smarts.FG_SMARTS` 逐条登记局部检测，模式首原子即中心原子、同名多条取并集，新增官能团只加表项。`analyze` 补非局部判据：`_oxo_kind` 归一含氧酸臂型，`oxoacid_entries` / `boronic_entries` 处理 P/S/B 中心，`_arbitrate_parts` 做 P-41 降级仲裁，`_drop_claimed_cations` 回收已被其它检测器命中的阳离子。出口为 `info`（`carbon_ids`、`fg_inventory`、`double_bonds` / `triple_bonds`、`rings` / `ring_systems`）。
- **L2 母体选择 `layer2/`**：`select_parent` = `rule_driven_parent_candidates` 收集 → `_reorder_p45_2` 按 P-45.2.1 前缀数与 P-67.2.1 缩合含氧酸桥氧数重排 → 返回并列最优候选组。骨架侧 `parent_skeleton` 按 P-44.1.2 / 44.2 / 44.3 / 44.4 逐级窄化；`principal_expression` 为选中骨架标注主基团 facts，并把环系一次性分派到稠合（P-25）、桥环（P-23）、螺环（P-24）三族之一；`finalize_parent_ownership` 固化 `owned_atoms`；`kind_registry.pack_parent_stem` 打包母体词干。
- **L3 取代基提取 `layer3/`**：`extract_substituents` 以 `claimable_block.iter_claims` 为唯一提取路径，`substituent_namer` 先查 `anchored_table` 的锚定 canonical-SMILES 表、未命中回退递归命名（`as_substituent.name_as_substituent`），未命名成功的 claim 静默跳过。`coverage.build_coverage_ledger` 产出重原子覆盖台账，供 Debug 视图与契约测试消费。
- **L4 编号与位次 `layer4/`**：`number` 取 `orient_numbering` 的定向结果后算位次、指示氢（`indicated_hydrogen`）与指示氢/加氢前缀（`_extra_indicated` / `_extra_hydrogenated`）。`orient_numbering` 按族分派：组分式螺环 → 螺环 → 桥环 → 固定编号骨架 → 稠环 → 通用 P-14.4 枚举；`locant_calc` 出附着原子位次、FG 位次省略规则与并列候选比较键。
- **L5 名称组装 `layer5/`**：`assemble` 先经 `_ensure_parent_stem` 注入母体词干（组分式螺环 → 螺环 → 桥环 → 未注册稠环 → 大环生成式），再由 `_names_for` 走链引擎取名，依次接 `join_hydro_prefix`、`join_ring_cation_suffix`、`join_kind_name`（前缀渲染）、`join_anion_names`、`join_ez_prefix`、`join_rs_prefix`，输出 `NameResult`。

## 协调器 namer.py

`SMILESNNamer` 持有 `tools.common_names.CommonNameCache`（默认 20000 条）：`name` 先查缓存，命中即返回；未命中则 `memo.begin_run` 开一次运行级记忆，走 `_pipeline` 并只回写成功结果。`_canonical_result` 把 `meta.parent_chain` 重写为 RDKit `CanonicalRankAtoms` 序，使缓存命中的元数据与首算一致。

`_name_mol` 编排主体：`dissociate_salt` → `analyze` → `_run_candidates` → `_apply_salt_suffix`。`info["root_ctx"]` 携带整分子根与原子映射，供锚定碎片回指宿主；`info["salt"]` 是磷酸母体生产者与盐门控的来源。带跨分子缓存时顶层整分子另开一份 2000 条运行内缓存（`*` 片段名只在同根内共享），锚定碎片无盐则沿用传入缓存。

`_run_candidates` 分两段：`select_parent` 的候选逐个 `finalize_parent_ownership` 后 `extract_substituents`，再由 `_assemble_candidate` 跑 `number` + `assemble`，`meta` 写入裁决键 `p44_1_1_key` / `p45_2_2_key`（取自 `locant_calc.suffix_locant_set` / `prefix_locant_set`）。`_best_hit` 只在此二键上裁决：后缀位次键不一致说明 P-44.1.1 已决，保持候选原序取首选；一致时取 P-45.2.2 前缀位次最小者。全部失败则 `_fail` 返回 `no_assemblable_candidate`。

`_subs_for_numbering` 只放行 O 侧臂、链上附着与 `N_PREFIX_KINDS` 三类取代基；`_remap_attach` 把不落在母体链上的附着原子改挂到 `ring_attach_idx` 或单锚点官能团锚点。

盐后缀在组装之后由 `_apply_salt_suffix` 经 `stems.join_metal_salt_names` 补齐；`OXO_CENTER_KINDS` 中心母体自带盐组装，故跳过。成功结果另带诊断元数据：`parent_chain` / `parent_kind` / `parent_labels`（取自 `numbering_scaffold.labels` 的整体编号标签，如稠环桥头 `3a`、`6a`）/ `bridge_self_enclosed` / `parent_substituent_count` / `salt`。失败经 `_fail` 返回，`meta.reason` 记 `parse` 或 `no_assemblable_candidate`。

## 跨层设计模式

| 模式 | 承载点 | 要点 |
| --- | --- | --- |
| 表驱动单一事实来源 | `layer1/fg_registry.FG_SPECS`；`layer5/chain_engine._KIND_TABLE`；`constants.py` | FG 类别、kind 词尾、跨层词表各只有一处登记，下游表由之派生 |
| kind 正交化 | L2 `principal` / `kind_registry` → L4 / L5 | 主基团类别压成 `kind` 字符串，L4/L5 只认 `kind`、不认 FG 枚举 |
| 含氧酸合一管线 | L1 `analyzer` → `fg_registry.oxoacid_is_acid` → L5 `chain_engine._OXO_TAIL` | P/S/B 中心共用一条链，含 S–O–S 多硫酸链 |
| 多环拆解择一 | L2 `principal_expression` / `fused_system` / `bridged_system` / `spiro_system` → L4 / L5 | 稠合、桥环、螺环三族互斥，一次识别即定 `scaffold_id` |
| 候选枚举与裁决分离 | L2 `*_system` → L4 `*_numbering` → L5 `*_namer` | L2 只交并列最优候选，L4 裁决并写回 `*_node`，L5 只读 |
| 原子归属 | L2 `parent_select.finalize_parent_ownership` → L3 `claimable_block` | `owned_atoms` 划定 L3 的 claim 枚举边界 |
| P-41 优先级接力 | `FgSpec.p41` → L2 `principal._effective_priority` / `parent_skeleton` → L5 | 优先级声明一次，骨架筛选与出词尾共用 |
| 排序键统一 | `tools/re.alpha_order_key`、`layer4/locant_calc.locant_key` | P-14.5 字母数字序与位次序各只有一份实现，编号与组装两侧共用 |
| 上层不反向依赖 | L5 `spiro_namer` / `bridged_namer` | node 按鸭子类型读 `descriptor` / `components` / `links`，L5 不 import L2 |
| 围栏由产生方定形 | L3 → L5 `assembler_prefixes` | L3 决定取代基是否自带围栏，L5 只做渲染与必要的整体围栏 |

- **`FG_SPECS` 与 `_KIND_TABLE`**：`FgSpec` 一条记录同时描述 FG 枚举值、L1 检测列表键、锚点字段、P-41 等级、P-43 路径与表达类型；L2 的 `PRINCIPAL_REGISTRY`、`functional_group_inventory._FG_KEYS` 与 `_ANCHOR_KEYS` 全部由它派生，扩展官能团只需加一条表项。命名端同理：`_KIND_TABLE` 每个 `kind` 一条 `_Chain` 规格（词尾、保留名 variant、位次省略规则）。
- **含氧酸合一管线**：`FG_SMARTS` 按「元素 + 双键氧数 + 臂型」逐条登记 P/S/B 中心，`analyzer` 把臂角色归一为 `oxo_kind`（`phosphate` / `phosphonate` / `sulfate` / `sulfonic` / `sulfonamide` / `boronic` 等）；`oxoacid_is_acid` 判是否按 P-41 类 7 酸式参与主基团竞争，否则落类 9 酯或类 11 磺酰胺，`_OXO_TAIL` 按 `(oxo_kind, n_oh)` 出词尾。
- **多环拆解择一**：一个环系要么走 P-25 稠合、P-23 桥环、要么走 P-24 螺环，不并存。`principal_expression._ring_scaffold_and_nodes` 一次算好三族候选：桥环命中即用 `BridgedNode.scaffold_identity()` 覆盖 scaffold 身份；螺环按 `SPIRO_SCAFFOLDS` 分 `mono_spiro` 与 `fused_bridged_spiro`，全单环组分走 P-24.2 螺描述符，含多环组分走 P-24.5~24.7 组分式命名。
- **候选枚举与裁决分离**：桥环与螺环各自只交出并列最优候选（`bridged_nodes` / `spiro_nodes` / `fbs_nodes`），编号交 L4 按 P-23.3.2 / P-24.2 / P-24.5 裁决后写回 `bridged_node` / `spiro_node` / `fbs_node`；命中即不再下落，L5 只读描述符、组分与上标位次对。
- **排序键统一**：`alpha_order_key` 循环剥除 N- 前缀、斜体前缀、位次、立体描述符与外层括号，先比字母序列、再比首字母前位次；`locant_key` 把数字位次按数值排、字母尾与撇号作次级键。`numbering_engine` 的 P-14.4(f) 平局与 L5 前缀排序共用前者，组分式螺环的多撇号排序用后者。
- **围栏由产生方定形**：`SubstituentName.requires_parentheses` 随取代基一路传到 L5，L5 仅在渲染前缀时套用（`assembler_prefixes._enclose`、`oxo_arm_fence`），并对整名补必要围栏，不在组装侧二次拆分。

## 文件组织

```
src/namepredict/
├── __init__.py  包装入口：安装 rdkit 迭代器加速
├── namer.py  顶层协调器：缓存 + L0–L5 编排
├── constants.py  跨层常量：原子序数、盐/阳离子词表、倍数前缀、保留名别名
├── types.py  NameResult 数据类
├── layer0/
│   ├── __init__.py  L0 入口重导出
│   ├── preprocessor.py  preprocess：清洗、解析与互变归一链
│   ├── tautomer.py  酮式/内酰胺/硫酮/胺式/脒式归一
│   ├── charge.py  酸性质子重定位与电荷归一
│   └── salt.py  dissociate_salt：反离子识别与盐元数据
├── layer1/
│   ├── __init__.py  L1 入口重导出
│   ├── fg_local_smarts.py  FG_SMARTS 局部检测表（中心原子）
│   ├── analyzer.py  analyze：非局部后处理、oxo_kind 归一、降级仲裁
│   ├── fg_registry.py  FG_SPECS：FG 跨层元数据唯一注册处
│   ├── functional_group_inventory.py  类型化清单与 occurrence
│   └── ring_systems.py  build_ring_systems：环系拓扑（L2/L4/L5 共用）
├── layer2/
│   ├── __init__.py  L2 入口重导出
│   ├── parent_select.py  select_parent 入口 + owned_atoms 固化
│   ├── parent_skeleton.py  P-44 拓扑优先骨架枚举与逐级窄化
│   ├── principal.py  P-41 / P-43 主基团选择
│   ├── principal_expression.py  主基团表达式 facts 与环族分派
│   ├── ring_expression_policy.py  环骨架类型化表达的能力策略
│   ├── chain_walk.py  脂肪族碳链行走
│   ├── fused_system.py  P-25.3.2.4 稠环组分树
│   ├── bridged_system.py  P-23 扩展 von Baeyer 桥环拆解
│   ├── spiro_system.py  P-24 螺环拆解：组分划分与螺描述符候选
│   ├── ring_scaffold.py  骨架规格 + 保留 SMILES 模板表 + 环解析
│   └── kind_registry.py  kind 词干注册与词干打包
├── layer3/
│   ├── __init__.py  L3 入口重导出
│   ├── substituent_extractor.py  extract_substituents 入口
│   ├── claimable_block.py  claim 枚举、slot 推导与块切割
│   ├── substituent_namer.py  命名后端：锚定表 → 递归
│   ├── coverage.py  重原子覆盖台账（Debug 与测试消费）
│   ├── submol_build.py  诱导切割子分子（H 封端）
│   └── as_substituent.py  子分子 → free-name → -yl
├── layer4/
│   ├── __init__.py  L4 入口重导出
│   ├── numbering.py  number 入口 + 指示氢/加氢前缀
│   ├── numbering_engine.py  orient_numbering：P-14.4 枚举与逐族分派
│   ├── locant_calc.py  附着原子 → 位次、位次省略规则与候选比较键
│   ├── fused_numbering.py  稠环编号：外周序、稠合位次、CIP 破局
│   ├── fused_orientation.py  稠环优选取向
│   ├── bridged_numbering.py  P-23 桥环编号裁决
│   ├── spiro_numbering.py  P-24 螺环编号：描述符候选与逐组分编号
│   ├── ring_geometry.py  环几何坐标与变形/重叠阈值
│   └── indicated_hydrogen.py  指示氢位次
├── layer5/
│   ├── __init__.py  L5 入口重导出
│   ├── assembler.py  assemble 入口：取名 → 前缀 → 阴离子/立体后缀
│   ├── chain_engine.py  链式词干引擎与 _KIND_TABLE
│   ├── assembler_prefixes.py  前缀分组与双语渲染、围栏
│   ├── stems.py  碳数词干表、金属盐与阴离子后缀
│   ├── stereo.py  E/Z（P-91.2）与 CIP R/S（P-92）前缀
│   ├── fused_namer.py  稠合名组装
│   ├── bridged_namer.py  P-23 桥环母体名（环数词头 + 描述符）
│   ├── spiro_namer.py  P-24 螺环母体名（螺描述符 / 组分式）
│   └── skeleton_replacement.py  骨架置换 'a' 前缀
└── tools/
    ├── __init__.py  跨层拓扑工具包标记
    ├── anchored_table.py  取代基锚定 canonical-SMILES 查表
    ├── re.py  名称文本辅助、alpha_order_key 与 SUB_LOCANT_RE
    ├── common_names.py  SMILES → NameResult 内存缓存
    ├── memo.py  单次命名内的中间结果记忆
    ├── block_cut.py  母体边界块切割
    ├── rdkit_fast.py  RDKit 迭代器加速补丁
    └── chain.py  开链最长路径搜索
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
