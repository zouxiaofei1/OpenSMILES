# 最差批次错误模式报告（SMILES → IUPAC 命名器）

## 口径与范围

- 数据：`data/merged_benchmark.json` 全量 3989 条（gold 为人工核校的双语名）与现役 `src/` 命名器逐条比对。
- 入选条件：`sim_max < 0.4`（英/中两侧取较高者，difflib SequenceMatcher 比值）**或**任一侧输出为空。
- 入选规模：**375 条**（占全量 9.4%），其中完全无输出 196 条；来源 chebi 328 条、smiles_tiers 47 条。
- 结构规模：重原子数 4–383（中位 30）；含 ≥2 个环的 331 条；含 ≥2 个环系统的 212 条。
- 排序：全表按 RDKit 复杂度（重原子数 → 环数 → 分子量）**从简单到复杂**排序；报告中每处 id 列表均沿用该顺序。
- 分析方法：375 条 round-robin 均分为 8 片，每片由一个独立 agent 逐条比对 gold/pred 并归纳错误模式，要求同一模式内机制一致；责任模块由 agent 在 `src/` 中 Grep/Read 或最小复现验证后给出，非推测。
- 归并结果：16 个主要模式 + 7 个零散尾类，375 条全部归类，无遗漏。

## 模式总览

| # | 模式 | 条数 | 其中无输出 | 占比 |
|---|---|---:|---:|---:|
| 1 | 多环母体词干组装失败（稠合/桥/螺未细分） | 78 | 78 | 20.8% |
| 2 | coverage 门控放宽：残缺局部片段名被当作成功输出 | 65 | 0 | 17.3% |
| 3 | 稠合环系拆解退化/稠合词干失败 | 38 | 35 | 10.1% |
| 4 | 杂原子中心（P/B）被当作碳命名 | 29 | 0 | 7.7% |
| 5 | 桥环/螺环/笼状母体缺 von Baeyer 词干 | 29 | 27 | 7.7% |
| 6 | 未注册单环/杂环母体缺失（HW 模板与保留名注册表缺项） | 20 | 19 | 5.3% |
| 7 | 离子/电荷/反离子未体现（季铵、阴离子、叶立德） | 20 | 1 | 5.3% |
| 8 | 链式同种主基团多重性不可表达（二酯/二酰胺/酸酐/偕二醇） | 18 | 18 | 4.8% |
| 9 | 多环骨架母体名生成错误（伪稠合/拓扑不符/位次缺失） | 17 | 0 | 4.5% |
| 10 | 取代基递归命名残缺（块内原子静默丢弃） | 9 | 0 | 2.4% |
| 11 | 母体选择与主基优先级错误 | 9 | 3 | 2.4% |
| 12 | 寡糖/糖苷母体单糖选取与位次、立体描述错误 | 8 | 0 | 2.1% |
| 13 | 母体候选为空（骨架表达失败） | 8 | 8 | 2.1% |
| 14 | 酰胺/内酰胺互变异构未归一 | 8 | 0 | 2.1% |
| 15 | 磺酸/硫酸（氢）酯主基团未登记，降级为前缀 | 6 | 0 | 1.6% |
| 16 | 大环母体/环几何模板缺失（含编号引擎 KeyError） | 5 | 5 | 1.3% |
| 17 | 环上取代基位次编号方向错误 | 2 | 0 | 0.5% |
| 18 | G gold 存疑 | 1 | 0 | 0.3% |
| 19 | 取代基前缀字母序与保留前缀未命中 | 1 | 0 | 0.3% |
| 20 | 多组分体系只命名单一组分 | 1 | 0 | 0.3% |
| 21 | 寡核苷酸链无命名支持 | 1 | 1 | 0.3% |
| 22 | 稠环氢化程度/指示氢表达错误 | 1 | 0 | 0.3% |
| 23 | 金属有机/配位化合物无母体 | 1 | 1 | 0.3% |

## 1. 多环母体词干组装失败（稠合/桥/螺未细分）

**条数** 78（无输出 78）｜**出现分片** 01, 02, 05, 06, 07

- **判定特征**：预测两侧皆空（内部 reason=no_assemblable_candidate）。分子含 ≥2 环的非苯环多环骨架（稠合杂环、桥环、螺环、大环），gold 用 bicyclo/tricyclo/spiro[…] 或稠合母体名。诊断可见母体带 fused_tree 但 stem_en/stem_zh 为空。本类由 5 个分片各自独立归出，未再细分稠合与桥环两条子路径。
- **根因**：L5 的母体词干只有两条来源：L2 注入的保留 scaffold 词干，或 `_ensure_fused_stem → fused_parent_names` 对 fused_tree 现场组装稠合 base 名。后者只覆盖「相邻环共享一条边」的邻位稠合拓扑——既没有 von Baeyer 桥环/多环（P-23）与螺环（P-24）的组装分支，也未注册青霉烷/头孢烯/安甾等桥环保留名；同时它要求根组分已有词干且 attached 非空，而 L2 `_decompose` 会把匹配不到保留模板的环静默丢弃，使 attached 变空。任一环节缺失即返回 None → assemble 走 `_unsupported`，全部候选作废、整分子无输出，且管线不降级、不报错。
- **IUPAC 依据**：P-23（von Baeyer 双环及多环）、P-24（螺环）、P-25.3.2（稠合名组装与稠合描述符）、P-25.3.3（稠环编号）、P-101.1（保留母体名）
- **责任模块**：src/namepredict/layer5/fused_namer.py:fused_parent_names；src/namepredict/layer5/assembler.py:_ensure_fused_stem；src/namepredict/layer2/fused_system.py:_decompose

**代表样例**（序号为复杂度排名，越小的分子越简单）

- 25 chebi-2415 | 1,3,3-trimethylbicyclo[2.2.1]heptane | (空)
- 105 chebi-2477 | (1R,4S,5R,8S,9R,12S,13R)-1,5,9-trimethyl-11,14,15,16-tetraoxatetracyclo[10.3.1.04,13.08,13]hexadecan-10-one | (空)
- 185 chebi-84 | (6R,7R)-3-[(5-methyl-1,3,4-thiadiazol-2-yl)sulfanylmethyl]-8-oxo-7-[[2-(tetrazol-1-yl)acetyl]amino]-5-thia-1-azabicyclo[4.2.0]oct-2-ene-2-carboxylic acid | (空)

**全部 78 个 id（按复杂度升序）**

```text
  chebi-2415 chebi-885 tiers-66561 tiers-44245 chebi-3151 tiers-17996 tiers-85289 chebi-1796
  chebi-82 chebi-1147 chebi-2221 chebi-1256 tiers-68910 chebi-474 chebi-1053 chebi-2517
  chebi-1703 chebi-1081 chebi-1588 chebi-2811 tiers-129945 chebi-2477 chebi-3297 chebi-1872
  chebi-403 chebi-2143 chebi-342 chebi-1356 chebi-2311 chebi-1484 chebi-1908 chebi-2406
  chebi-3266 chebi-3078 chebi-1305 chebi-434 chebi-606 chebi-2022 chebi-2999 chebi-2041
  chebi-1389 chebi-1248 chebi-1086 chebi-129 chebi-84 chebi-1240 chebi-526 chebi-2714
  chebi-64 tiers-131559 tiers-148964 chebi-1017 chebi-568 chebi-1295 chebi-665 chebi-1985
  chebi-1905 chebi-3262 chebi-2651 chebi-1186 chebi-3295 chebi-924 chebi-1593 chebi-788
  chebi-2789 chebi-2386 chebi-2952 chebi-473 chebi-156 chebi-1896 chebi-1466 chebi-2839
  chebi-2923 chebi-948 chebi-2954 chebi-1665 chebi-2013 chebi-3004
```

## 2. coverage 门控放宽：残缺局部片段名被当作成功输出

**条数** 65（无输出 0）｜**出现分片** 01, 02, 03, 04, 05, 06, 07, 08

- **判定特征**：pred 非空但只描述分子的一小块（常为 gold 名的 1/5~1/8），形态是 acetate / methanol / ethane / methyl / formamide 这类小片段官能团名；`meta.fallback == "no_coverage_gate"`、`coverage_complete=False`（少数被误判为 True）。顶层命名与取代基递归层同样生效。
- **根因**：L3 coverage ledger 要求母体 owned_atoms 与全部取代基原子覆盖整分子。当分子里存在既非母体、又无法提取为可命名取代基的原子（未支持的大环/稠环、糖链、季铵中心、磷酸二酯链、非盐多组分）时，所有候选的 coverage 均不完整，`_complete_hit` 全落空，`_try_phase` 便无条件退到 `_partial_hit`——该分支完全绕过 coverage 门控直接组装，未覆盖原子被静默丢弃，结果仍返回 `success=True`。取代基递归命名（`_radical_yl_from_sub → _name_mol`）内部同样触发该回退。输出前没有任何环节校验名称与分子原子集的对应关系。
- **IUPAC 依据**：P-1.1 与 P-14（名称须完整描述整个分子）、P-44.1（母体选择须覆盖全部骨架原子）、P-44.2（环系与链的母体优先级）、P-45.2（前缀位次完整性）
- **责任模块**：src/namepredict/namer.py:_partial_hit；src/namepredict/namer.py:_run_candidates；src/namepredict/layer3/coverage.py:build_coverage_ledger

**代表样例**（序号为复杂度排名，越小的分子越简单）

- 81 tiers-53012 | 5-piperidin-4-yl-3-pyrazin-2-yl-1,2,4-oxadiazole | pyrazine
- 153 chebi-424 | (5Z)-7-[(3R,4S)-3-[(1E,3S)-3,7-dihydroxyoct-1-enyl]-2-oxabicyclo[3.1.1]heptan-4-yl]hept-5-enoate | (5Z)-hept-5-enoate
- 345 chebi-2492 | N,N-diethylethanamine 5-[2-[2,3-di(hexadecanoyloxy)propoxy-hydroxyphosphoryl]oxyethylsulfamoyl]-2-(3-oxa-23-aza-9-azoniaheptacyclo[17.7.1...]octacosa-…-16-yl)benzenesulfonate | 3-[[2-[[4-(…pyrido[3,2-i]…

**全部 65 个 id（按复杂度升序）**

```text
  chebi-2334 chebi-3285 tiers-80645 tiers-53012 chebi-1338 chebi-2616 chebi-2136 chebi-1979
  chebi-2427 chebi-2327 chebi-1196 chebi-2242 chebi-424 chebi-2938 chebi-3085 chebi-2674
  chebi-2351 chebi-2128 chebi-782 tiers-140247 chebi-968 chebi-2791 chebi-3226 chebi-17
  chebi-2184 tiers-112188 chebi-406 chebi-1330 chebi-2709 chebi-2426 chebi-258 chebi-1197
  chebi-488 chebi-2453 chebi-1187 chebi-520 chebi-332 chebi-1661 chebi-2333 chebi-2176
  chebi-2963 chebi-3130 chebi-1461 chebi-1082 chebi-249 chebi-1599 chebi-2702 chebi-2665
  chebi-1857 chebi-3147 chebi-3275 chebi-1066 chebi-1227 chebi-599 chebi-1034 chebi-68
  chebi-3029 chebi-3083 chebi-2586 chebi-2321 chebi-1077 chebi-1918 chebi-2492 chebi-3221
  chebi-2149
```

## 3. 稠合环系拆解退化/稠合词干失败

**条数** 38（无输出 35）｜**出现分片** 02, 03, 04, 08

- **判定特征**：两种表现：(a) 输出空名——环系含 ≥2 个 SSSR 环且 L2 已产出 fused_tree，但拆解树退化为无附加组分的单节点，或根组分无保留词干；(b) 有输出但稠合描述符残缺——出现 `pyrido[-a]pyrimidin`、`[,-a]`、`[-c]` 这类空位次，或母体覆盖的环数少于分子实际环数（4 环体系只写出 3 环，9 环体系只覆盖 3 环）。
- **根因**：L2 `_decompose` 按「相邻环共享 2 个原子」建稠合边，且只收容能精确匹配保留模板的组分（`match_fusion_component`）；未注册的稠环组分（氮杂䓬、二氧杂环、薁型 5-7 稠环、phenalenone 角环、四环素型角稠环等）匹配失败即被 `if node is not None` 静默丢弃，拆解树退化为单节点；而 `compute_owned_atoms` 仍把这些环记为已归属，coverage 门控无法察觉。L5 `_fused_one` 只校验稠合字母可取就拼接 `[位次-字母]`，不校验位次集合（共享边两原子未同时落在子组分编号链时 `_fusion_numbers` 返回空元组），于是产出空位次残名。
- **IUPAC 依据**：P-25.3.2（稠合命名与稠合描述符位次-字母）、P-25.3.2.4（母体组分选择准则）、P-25.3.3（稠环编号）、P-25.2.1（保留稠环组分名）、P-23.2.5（无法稠合命名时改用 von Baeyer）
- **责任模块**：src/namepredict/layer2/fused_system.py:_decompose；src/namepredict/layer5/fused_namer.py:_fused_one

**代表样例**（序号为复杂度排名，越小的分子越简单）

- 178 chebi-1692 | [(4S,5R,6S,8S,10R)-10-[(R)-(2,4-dioxo-1H-pyrimidin-6-yl)-hydroxymethyl]-5-methyl-2,11-diaza-12-azoniatricyclo[6.3.1.04,12]dodec-1(12)-en-6-yl] sulfate | 2-hydroxy-6-[hydroxy-[(2aS,3R,4S,5aS,7R)-3-methyl…
- 27 chebi-3137 | 6,6-dimethyl-2-methylidenebicyclo[3.1.1]heptane | (空)
- 35 chebi-1503 | (1R,3R,4S,5S)-4-methyl-1-propan-2-ylbicyclo[3.1.0]hexan-3-ol | (空)

**全部 38 个 id（按复杂度升序）**

```text
  chebi-3137 chebi-1503 tiers-65526 chebi-1919 tiers-59592 tiers-96844 tiers-7355 tiers-18169
  chebi-298 chebi-846 tiers-124039 chebi-1459 tiers-30123 chebi-1864 chebi-1560 chebi-2696
  chebi-263 chebi-967 chebi-2409 chebi-2257 chebi-289 chebi-604 chebi-2589 chebi-1567
  chebi-2283 chebi-1417 chebi-162 chebi-1692 tiers-96196 chebi-929 chebi-1501 chebi-1636
  chebi-1512 chebi-1135 chebi-3119 chebi-1488 chebi-1238 chebi-2062
```

## 4. 杂原子中心（P/B）被当作碳命名

**条数** 29（无输出 0）｜**出现分片** 02, 03, 04, 05, 06, 07, 08

- **判定特征**：pred 把含 P 的片段写成碳醚/醇——出现 `hydroxymethoxy`、`hydroxy(ethoxy)methoxy`、`[hydroxy(phosphonooxy)methoxy]` 这类串，或把 B(OH)₂ 写成 `dihydroxymethyl`；gold 对应位置是 phosphonooxy / oxidophosphoryl / …phosphate / boronic acid。典型分子是含 P–O–P 焦磷酸桥、磷酸二酯桥或多磷酸链的核苷酸、辅酶、磷脂。
- **根因**：L1 `phosphate.py:_one_phosphate` 末尾的「整分子纯度」校验要求 core ∪ 臂 = 全部重原子（注释明写用于排除 P–O–P），凡含 P–O–P 或多磷酸链的分子，每个 P 都判不过而被拒收为磷酸官能团。L3 取代基侧，`tools/anchored_table.py` 只为 P(=O)(OH)₂ 型单酯登记了 phosphono/phosphonooxy（`*OP(=O)(O)O` / `*OP(=O)([O-])O`），磷酸二酯与焦磷酸没有任何条目。未被识别的 P 于是落入通用碳骨架路径，P 上的 =O/OH 变成碳上的羟基、桥氧变成甲氧基，名称少一个 P、多一个碳。硼酸同理：元素表里有 B，但 FG_SPECS 无硼酸条目，B 被当作四价碳纳入烃链。
- **IUPAC 依据**：P-67.1.3（单核非碳酸含氧酸的酯，如 methyl dihydrogen phosphate）、P-67.1.4（多核含氧酸与焦磷酸）、P-67.1.5.1（磷酸降级为 phosphonooxy 等前缀）、P-68.1（硼酸）、P-68.3（第 15 族元素母体氢化物 phosphane/phosphoryl）
- **责任模块**：src/namepredict/layer1/phosphate.py:_one_phosphate；src/namepredict/tools/anchored_table.py:_REGISTRY；src/namepredict/layer3/as_substituent.py:name_as_substituent

**代表样例**（序号为复杂度排名，越小的分子越简单）

- 290 chebi-3120 | S-[2-[3-[[(2R)-4-[[[(2R,3S,4R,5R)-5-(6-aminopurin-9-yl)-4-hydroxy-3-phosphonooxyoxolan-2-yl]methoxy-hydroxyphosphoryl]oxy-hydroxyphosphoryl]oxy-2-hydroxy-3,3-dimethylbutanoyl]amino]propanoylamino]ethyl]…
- 243 chebi-990 | [(2R,3R,4R,5S,6R)-3-acetamido-4,5-dihydroxy-6-(hydroxymethyl)oxan-2-yl] [[(2R,3S,4R,5R)-5-(2,4-dioxopyrimidin-1-yl)-3,4-dihydroxyoxolan-2-yl]methoxy-hydroxyphosphoryl] hydrogen phosphate | N-[(2R,3R,4R,5…
- 339 chebi-3287 | ...(2Z,6Z,...,38E)-3,7,11,15,19,23,27,31,35,39,43-undecamethyltetratetraconta-2,6,10,14,18,22,26,30,34,38,42-undecaenoxy]phosphoryl] phosphate | ...(2Z,...)-3,7,...-undecaenyloxy]methoxy]methoxy]-5-(1-h…

**全部 29 个 id（按复杂度升序）**

```text
  tiers-43418 tiers-8543 chebi-571 chebi-369 chebi-2204 chebi-1169 chebi-990 chebi-1729
  chebi-3120 chebi-2534 chebi-1394 chebi-2951 chebi-902 chebi-2054 chebi-2304 chebi-1280
  chebi-3098 chebi-1676 chebi-3287 chebi-2624 chebi-1830 chebi-183 chebi-766 chebi-559
  chebi-2411 chebi-3244 chebi-1645 chebi-352 chebi-1544
```

## 5. 桥环/螺环/笼状母体缺 von Baeyer 词干

**条数** 29（无输出 27）｜**出现分片** 04, 07, 08

- **判定特征**：pred 为空（或只剩能命名的小侧链）。环系含桥头原子（两环共享 ≥3 原子，或某原子属 ≥3 个 SSSR 环）或螺原子（两环共享 1 原子），gold 直接用 bicyclo/tricyclo/tetracyclo/pentacyclo[…] 或 spiro[…] 作母体。最小复现：bicyclo[2.2.1]heptane、bicyclo[2.2.2]octane、camphor、penam 型 O=C1CC2SCCCN12 一律无输出。
- **根因**：当前管线只有稠合名（P-25）组装器，完全没有 von Baeyer 桥环/螺环母体生成器。L1 `ring_systems` 虽算出 bridgeheads/bridge_lengths，L2 `_decompose` 仍按「共享 2 原子」建稠合边；桥环（共享 1 或 ≥3 原子）无稠合边 → 剩余环各自成为独立分量、又匹配不到保留模板 → 拆解树退化为单节点。L5 `fused_parent_names` 遇 `not node.attached` 或根无词干立即返回 None → unsupported → 空输出。原实现 `layer2/scaffold/bridged_parent.py`（c040df0 引入）已在 b823cc7 被删除，当前 src 中已无任何 bicyclo/tricyclo 生成代码。
- **IUPAC 依据**：P-23.2（von Baeyer 桥环母体氢化物与编号：P-23.2.2 双环、P-23.2.5 三环、P-23.2.6 多环）、P-23.3（杂环 von Baeyer，骨架置换 a 前缀）、P-24.2（螺环母体氢化物）
- **责任模块**：src/namepredict/layer5/fused_namer.py:fused_parent_names；src/namepredict/layer2/fused_system.py:_decompose

**代表样例**（序号为复杂度排名，越小的分子越简单）

- 28 chebi-903 | 2,6,6-trimethylbicyclo[3.1.1]hept-2-ene | (空)
- 44 tiers-6342 | N-methyl-1,4-dioxaspiro[4.5]decan-8-amine | (空)
- 180 tiers-95650 | N-[3-(2-phenylmorpholin-4-yl)propyl]adamantane-1-carboxamide | (空)

**全部 29 个 id（按复杂度升序）**

```text
  chebi-1786 chebi-971 chebi-826 chebi-903 chebi-2147 chebi-2153 tiers-6342 chebi-1520
  chebi-594 chebi-1438 chebi-503 chebi-2303 chebi-226 chebi-2537 chebi-479 chebi-3288
  tiers-95650 chebi-676 chebi-3053 chebi-3208 chebi-125 chebi-1115 chebi-3064 chebi-14
  chebi-2925 chebi-2687 chebi-210 chebi-1854 chebi-1763
```

## 6. 未注册单环/杂环母体缺失（HW 模板与保留名注册表缺项）

**条数** 20（无输出 19）｜**出现分片** 01, 02, 03, 07, 08

- **判定特征**：pred 两侧皆空（或只剩环外小片段如甲醇/苯酚）。主骨架是 3–10 元杂单环（异噻唑酮、1,2,4-噁二唑、1,3,4-噁二唑、1,3-噻嗪、1,4-氧氮杂环庚烷、1,2,4-三嗪、1,2,4,5-四嗪等）或 7 元氧杂环内酯；诊断 `select_parent_tied` 返回空列表，即连候选母体都没有。对照：morpholine、pyrrolidin-2-one 可正常命名。
- **根因**：母体环候选只来自 `layer2/ring_scaffold.py:_TEMPLATES` 中登记了模板与词干的保留名杂环（`_BY_ID` 约 69 条），没有 Hantzsch–Widman（P-22.2.2，3–10 元）生成器；未命中时 `resolve_ring_scaffold` 返回 None，`_generic_ring_kind` 又只接受全碳环或 ≥2 个芳香环的稠环，其余一律 None。于是环骨架被整体判为不可表达，`express_ring_principal` 返回 None、候选集为空 → no_assemblable_candidate。若分子另有小官能团，锚点只能拉出那条小链，环整体在名称中丢失。
- **IUPAC 依据**：P-22.2.2（3–10 元杂单环 Hantzsch–Widman 名）、P-22.2.1 与表 2.2（杂单环母体保留名）、P-22.2.4（≥11 元 mancude 杂单环）、P-25.2（稠环保留名）、P-44.1（母体选择）
- **责任模块**：src/namepredict/layer2/principal_expression.py:_generic_ring_kind；src/namepredict/layer2/ring_scaffold.py:_TEMPLATES；src/namepredict/layer2/principal_expression.py:express_ring_principal

**代表样例**（序号为复杂度排名，越小的分子越简单）

- 9 chebi-2752 | 2-methyl-1,2-thiazol-3-one | (空)
- 73 chebi-868 | 3-phenyl-5-thiophen-2-yl-1,2,4-oxadiazole | (空)
- 98 chebi-164 | 3,6-bis(2-chlorophenyl)-1,2,4,5-tetrazine | (空)

**全部 20 个 id（按复杂度升序）**

```text
  chebi-2752 tiers-5188 tiers-58291 chebi-1403 chebi-3066 chebi-531 chebi-1615 chebi-432
  chebi-868 chebi-93 chebi-164 tiers-26894 chebi-3257 chebi-2251 chebi-1978 chebi-1020
  chebi-510 chebi-1146 chebi-435 chebi-2488
```

## 7. 离子/电荷/反离子未体现（季铵、阴离子、叶立德）

**条数** 20（无输出 1）｜**出现分片** 01, 02, 03, 04, 05, 07

- **判定特征**：SMILES 含形式电荷（[NH3+]/[N+]/[O-]/[S-]/[B-]）或独立反离子（.[Na+]/.[Cl-]），而 pred 用中性措辞：amino/氨基 代替 azanium/铵，hydroxy 代替 oxidophenyl/-olate，hydroxysulfonyl 代替 sulfonate/硫酸酯；或把带电中心整块丢弃（季铵 N+ 及其 4 个取代基整体消失）。探针：同一 SMILES 去掉 .[Na+] 后可正常命名。
- **根因**：电荷在 L1 之后就基本失去表达通道。(1) `fg_registry.FG_SPECS` 里 quaternary_ammonium 既无 anchors 也无前缀通道，`analyzer._quaternary_ammonium_entries` 只登记不消费；(2) `analyzer._amine_degree` 不校验 FormalCharge，[NH3+] 按 H 数被当作一级胺写成 amino；(3) `analyzer._hydroxyl_entries` 只收中性羟基氧，[O-] 走普通氧路径写成 hydroxy，而 anion 标志只在 carboxyls/phosphates 的 payload 上传递，`stems.maybe_anion_names` 只处理羧酸根；(4) 硫酸/磺酸/硼酸类目在官能团清单里不存在，只能退回锚定前缀表当普通取代基；(5) `layer0/charge.py:normalize_acid_charge` 的「负电荷收敛到最强酸位」策略在羧酸与弱酸阴离子共存时把羧酸质子搬到弱酸位点，得到与输入相反的电荷分布。
- **IUPAC 依据**：P-71（阴离子命名与 -olate/-ide 后缀）、P-72.2.2（-oate/-ate/-sulfonate 等阴离子后缀）、P-73.1/P-73.2（-ium/-ylium 阳离子后缀，azanium/ammonium）、P-73.4（azonia 等置换式阳离子前缀）、P-74.1（离子化合物与电荷表达）、P-41 表 4.1 第 4–6 类（阴离子/两性离子/阳离子优先于酸）
- **责任模块**：src/namepredict/layer1/fg_registry.py:FG_SPECS；src/namepredict/layer1/analyzer.py:_amine_degree；src/namepredict/layer0/charge.py:normalize_acid_charge；src/namepredict/layer5/stems.py:maybe_anion_names

**代表样例**（序号为复杂度排名，越小的分子越简单）

- 1 chebi-1152 | hydrogen carbonotrithioate | methanethiol
- 41 chebi-402 | (2-aminophenyl) sulfate | 2-(hydroxysulfonyloxy)aniline
- 297 chebi-1358 | dipotassium (2Z)-2-[(2E,4E,6E)-7-[5-carboxy-3,3-dimethyl-1-(4-sulfonatobutyl)indol-1-ium-2-yl]hepta-2,4,6-trienylidene]-3-ethyl-1,1-dimethylbenzo[e]indole-6,8-disulfonate | dipotassium 2-[7-[3-ethyl-6-(…

**全部 20 个 id（按复杂度升序）**

```text
  chebi-1152 chebi-2055 chebi-2628 chebi-1282 chebi-895 chebi-3086 chebi-402 chebi-384
  chebi-3280 chebi-1684 chebi-2707 chebi-648 chebi-3 chebi-378 tiers-67934 chebi-2503
  chebi-1916 chebi-2867 chebi-1358 chebi-982
```

## 8. 链式同种主基团多重性不可表达（二酯/二酰胺/酸酐/偕二醇）

**条数** 18（无输出 18）｜**出现分片** 01, 02, 03, 04, 07

- **判定特征**：pred 两侧皆空（reason=no_assemblable_candidate）。分子为开链，同一主官能团出现 2 次（二酯、二酰胺、二元酰卤），或属酸酐类（R-C(=O)-O-C(=O)-R），或同一碳上带两个羟基（偕二醇）。最小复现：oxalyl dibromide、acetyl acetate、NC(=O)CC(N)=O、CC(O)(O)C 全部无输出，而 NC(=O)CC(=O)O（酸+酰胺）与 OCCO 正常。
- **根因**：L2 `principal_expression._chain_kind` 用 fg_registry 的 multi 标志派生出白名单 `_MULTI_FG`（仅 acid/ketone/alcohol/amine/thiol 允许 count ≥2）；酯、酰胺、酰卤一旦出现 2 次即返回 None，酸酐则根本没注册进 chain FG 表。骨架枚举本身找到了含两个基团的链，但 `express_chain_principal` 返回 None，候选集为空 → 直接判 no_assemblable_candidate，未回退到 P-31.1.4 的倍增后缀（-dioate/-diamide）或酸酐类名表达。偕二醇是同一缺陷的另一形态：两个 -ol 落在同一锚原子，L5 多重后缀（-diol）对重复位次返回失败并落到 unsupported。
- **IUPAC 依据**：P-31.1.4（倍增后缀 -dioate/-diamide）、P-65.1.2（二酸保留名 ethanedioic acid/oxalic acid）、P-65.7.1（对称酸酐）、P-66.1.1.1.1.1（二元酰胺 pentanediamide）、P-66.1.1（后缀多重化 -diol）、P-65.6.3.3.2（单一酸组分的多元酯 butanedioate）
- **责任模块**：src/namepredict/layer2/principal_expression.py:_chain_kind；src/namepredict/layer2/principal_expression.py:express_chain_principal；src/namepredict/layer5/assembler.py:assemble

**代表样例**（序号为复杂度排名，越小的分子越简单）

- 17 chebi-2703 | dimethyl (2E)-but-2-enedioate | (空)
- 161 chebi-9 | dodecanoyl dodecanoate | (空)
- 129 tiers-128583 | N-(1,4-dioxaspiro[4.5]decan-3-ylmethyl)-N'-(3-fluorophenyl)oxamide | (空)

**全部 18 个 id（按复杂度升序）**

```text
  tiers-9441 tiers-93 tiers-180 chebi-2703 tiers-186 tiers-61915 tiers-1630 chebi-862
  tiers-128583 chebi-9 tiers-91385 tiers-107887 tiers-129440 tiers-131305 chebi-1363 chebi-47
  chebi-2398 chebi-2382
```

## 9. 多环骨架母体名生成错误（伪稠合/拓扑不符/位次缺失）

**条数** 17（无输出 0）｜**出现分片** 01, 03, 04, 05, 06, 07, 08

- **判定特征**：pred 有输出，但母体骨架名与分子拓扑不符：用稠合名 + dodecahydro/tetradecahydro 表达桥环、螺环或迫稠合体系（gold 用 bicyclo/tricyclo/spiro），或稠合描述符出现空位次 `[-a]`、非法位次 0，或隐含环原子数与分子不符（3 个六元环只有 13 个环原子却拼出隐含 14 原子的稠合名），或环系整体丢失退化成开链烷烃名（undecane）。
- **根因**：L2 `_select_base`/`_candidates_for` 只做「增长式子结构模板匹配」，不校验匹配到的母体是否与骨架同构、也不要求共享边为外周邻边；桥环的每个环都能各自命中六元模板（个别案例 4 个环全被认成 pyridine），拼成看似合理的稠合名后用 hydro 前缀补氢掩盖。螺连接（共享 1 原子）也被当作稠合处理，螺原子同时被算作稠合桥头。L1 `_compute_bridged_info`/`_topology` 对共享原子（peri）稠合不判为 bridged，仍标 topology=fused，于是按邻稠合拆组分。保留名覆盖不足（picene、xanthene 类缺失）时只能由稠合组分自组装，选出的基名与权威名不同，位次随之整体移位。
- **IUPAC 依据**：P-25.3.1（稠合命名只适用于邻位稠合体系）、P-23.2（不可邻稠合的体系用 von Baeyer）、P-24.2（螺环须用螺原子描述符）、P-25.3.2.1.3 与 P-25.3.2.4（保留名作优选稠合母体组分及组分选择准则）、P-25.1.1 与表 2.7（稠环保留名）
- **责任模块**：src/namepredict/layer2/fused_system.py:_select_base；src/namepredict/layer5/fused_namer.py:_fused_one；src/namepredict/layer1/ring_systems.py:_compute_bridged_info

**代表样例**（序号为复杂度排名，越小的分子越简单）

- 113 chebi-2594 | (1R,15R,17S,18S)-17-ethyl-3-aza-13-azoniapentacyclo[13.3.1.02,10.04,9.013,18]nonadeca-2(10),4,6,8-tetraene | 0-ethyl-3,4-dihydro-1H-pyrido[-b]indole
- 305 chebi-966 | (1R,2S,3S,5R,6S,7S,8R,10S,11S,14E,16Z,18R,19S,22R,24R,25S,26R,28S,30S,31R,33S)-…-4,12,27,29,32-pentaoxaoctacyclo[26.3.1.119,22.01,8.02,26.03,5.07,11.025,30]tritriaconta-14,16-dien-18-yl benzoate | …-hexa…
- 83 chebi-1236 | (1S,9S,10S)-7,15-diazatetracyclo[7.7.1.02,7.010,15]heptadec-2-ene | (3S,5S,14S)-3,4,5,8,10,11,12,13,14,15,16,17-dodecahydro-2H-pyrido[1,2-a]pyridine

**全部 17 个 id（按复杂度升序）**

```text
  tiers-29902 chebi-395 chebi-194 chebi-1236 chebi-1727 chebi-2576 chebi-2594 chebi-1698
  chebi-2061 chebi-822 chebi-197 tiers-121896 chebi-758 chebi-529 chebi-390 chebi-966
  chebi-3149
```

## 10. 取代基递归命名残缺（块内原子静默丢弃）

**条数** 9（无输出 0）｜**出现分片** 01, 03, 05, 06, 08

- **判定特征**：母体骨架大体正确，但某个前缀把 gold 中 20 原子以上的大块（多环、糖链、含 N+ 的胺块）描述成远小于它的名字（23 原子含四环只写 methyl，61 原子只写糖链，80 原子只写 8 原子片段）；或 `coverage_complete=True` 却少了一整个环。诊断：某取代基 claim 的原子数远大于其名字所覆盖的原子数。
- **根因**：L3 取代基命名走 `as_substituent._radical_yl_from_sub`——把 claim 原子集复制成带 `*` 锚点的子分子再递归调用 `_name_mol`。递归层内部同样遇到 coverage 不完整，于是在子分子层触发 `_partial_hit`，把「连接点那一两个原子」当结果返回（实测 23 原子块在 depth=1 返回 methyl）；返回路径只检查名字非空、不校验名字覆盖了 claim 的多少原子。另一种形态是 coverage 台账只比对 claim 原子集与母体原子集，不比对名字能否表达这些原子，于是磺酰胺块只渲染了磺酰基却被记为已命名；还有 `namer._subs_for_numbering` 只保留挂点在母体链上的取代基，挂点在另一取代基内部的支链被静默丢弃，而 coverage 校验集仍是 L3 全量取代基。
- **IUPAC 依据**：P-14.1（取代基须完整表达所辖原子）、P-22.1 与 P-29.2（取代基与自由基命名）、P-29.3（复合取代基前缀，如三甲基铵基）、P-44.3 与 P-45.2（母体与取代基的原子划分及前缀位次）
- **责任模块**：src/namepredict/layer3/as_substituent.py:_radical_yl_from_sub；src/namepredict/layer3/coverage.py:build_coverage_ledger；src/namepredict/namer.py:_subs_for_numbering

**代表样例**（序号为复杂度排名，越小的分子越简单）

- 209 chebi-3087 | 2-[[(1S,2S,5S,6S,7R,10R,11S)-5-hydroxy-5,7,11-trimethyl-12-oxo-13-oxatetracyclo[9.3.3.01,10.02,7]heptadecan-6-yl]methyl]-6-methoxycyclohexa-2,5-diene-1,4-dione | 2-methoxy-6-methylcyclohexa-2,5-diene-1,…
- 337 chebi-1548 | 2-[2-[[2-[(12S,19S,26Z,29S)-19-(2-carboxyethyl)-12-[(3,4-dimethyl-1H-indole-2-carbonyl)sulfanylmethyl]-…-nonazahexacyclo[30.2.1.18,11.115,18.122,25.02,7]octatriaconta-…-5-yl]-1,3-thiazole-4-carbonyl]ami…
- 321 chebi-2580 | (2S,3R,4S,5S,6R)-2-[(2S,3R,4S,5R,6R)-2-[(2R,3R,4R,5R,6R)-4,5-dihydroxy-2-(hydroxymethyl)-6-[(1R,2S,4S,5'S,…)-5',7,9,13-tetramethylspiro[5-oxapentacyclo[10.8.0.02,9.04,8.013,18]icosane-6,2'-piperidine]-1…

**全部 9 个 id（按复杂度升序）**

```text
  chebi-1907 chebi-1695 chebi-1283 chebi-1787 chebi-3087 chebi-151 chebi-2710 chebi-2580
  chebi-1548
```

## 11. 母体选择与主基优先级错误

**条数** 9（无输出 3）｜**出现分片** 01, 02, 05, 06, 07

- **判定特征**：构成式相同、仅母体名不同。(a) 分子同时含磷酸（酐）与硫酯/酯/酰胺时，pred 选磷酸作后缀、糖或核苷酸大块降为取代基，gold 反之；(b) O-C(=O)-O（碳酸酯）、N-C(=O)-N（脲）被写成 C1 保留名（formamide/formic acid）加杂原子前缀，gold 用 carbonate/urea 保留名；(c) 含 ≥3 个羧基时 pred 以二元酸为母体、其余写成 carboxy 前缀，gold 用能容纳最多后缀的链；(d) gold 用天然产物保留母体名（picene、corrin），pred 用稠合组分拼装；(e) 母体退化为单碳官能团碳（二硫代羧基碳、胍碳）或带电氮鎓小片段后无词干可拼 → 空输出。
- **根因**：多个独立缺口叠加。(1) `fg_registry` 给 phosphate 与 ester 同为 p41=9，磷酸靠 path 顺序胜出，而 P-41 表 4.1 要求羧酸衍生物（含硫代酸酯）优先于无机酸衍生物；(2) FunctionalGroupClass 里没有 carbonate/urea 类目，该羰基碳只按普通羧酸/酰胺处理，被选成 1 碳母体链，L5 `chain_engine._RETAINED` 对 n=1 直接给出 formic acid/formamide；(3) 链骨架比较先按最长链取，没有先应用 P-44.1.1「主特征基团作后缀数目最多者优先」；(4) `_RETAINED_FUSION_ALIASES` 只收了 benzo[c]furan/dibenzofuran 等 6 条，缺 picene 等天然产物保留名（src 中检索不到 picene/picen 登记）；(5) 官能团碳本身被当作母体氢化物，1 碳的 amine/thiol 在 L5 词干表中无表项。
- **IUPAC 依据**：P-41.1 与表 4.1（特性基团优选顺序；第 7a 类磺酸 > 第 11 类酰胺）、P-44.1 与 P-44.3（母体骨架须为链或环系，官能团碳本身不能作母体氢化物）、P-44.1.1（主特征基团作后缀数目最多者优先）、P-66.1.6.1.1.1（尿素 urea 保留名）、P-65.1.1（碳酸/碳酸酯保留名）、P-25.1.1 与 P-101（稠环保留名与天然产物母体氢化物）
- **责任模块**：src/namepredict/layer2/principal.py:select_principal_group；src/namepredict/layer2/parent_skeleton.py:select_principal_skeletons；src/namepredict/layer5/chain_engine.py:_RETAINED；src/namepredict/layer5/fused_namer.py:_RETAINED_FUSION_ALIASES

**代表样例**（序号为复杂度排名，越小的分子越简单）

- 289 chebi-203 | S-[2-[3-[[(2R)-4-[[[(2R,3S,4R,5R)-5-(6-aminopurin-9-yl)-4-hydroxy-3-phosphonooxyoxolan-2-yl]methoxy-hydroxyphosphoryl]oxy-hydroxyphosphoryl]oxy-2-hydroxy-3,3-dimethylbutanoyl]amino]propanoylamino]ethyl] …
- 10 chebi-2007 | [hydroxy(oxido)phosphoryl] carbonate | phosphonatooxyformate
- 170 tiers-118482 | 1-[1-(2,3-dihydro-1,4-benzodioxin-6-yl)-5-oxopyrrolidin-3-yl]-3-(3-methylpyridin-2-yl)urea | 1-(3-methylpyridin-2-ylamino)-N-(5-oxopyrrolidin-3-yl)formamide

**全部 9 个 id（按复杂度升序）**

```text
  chebi-709 chebi-2007 chebi-983 chebi-3281 chebi-1914 tiers-118482 chebi-2826 chebi-203
  chebi-2611
```

## 12. 寡糖/糖苷母体单糖选取与位次、立体描述错误

**条数** 8（无输出 0）｜**出现分片** 02, 03, 05, 06, 07, 08

- **判定特征**：pred 与 gold 都是覆盖完整的糖链长名（常 >1000 字符）、糖残基种类与个数一致，但选中的母体糖环不是同一个（pred 用带苷元的环、gold 用还原端），导致该环所有取代位次与 R/S 描述符整体移位；或同一糖环上两个支链的位次互换（5 位与 6 位对调）。对照：单糖、二糖、三糖直线链可正确命名。
- **根因**：每个糖环在 L2 被当作独立 RING_SYSTEM 候选，`parent_selector.select_parent_tied` 只按通用骨架评分与取代基数（P-44.1.1 后缀位次集合、P-45.2.2 前缀位次集合）裁决，没有实现 P-102.7.2.2 的低聚糖专用母体规则（有游离半缩醛基的低聚糖按 glycosyl[glycosyl]n glycose 命名、glycose 即还原端作母体），于是按枚举顺序或评分选到了苷元所在的那一环。母体一环之差让整条链的糖基残基方向与位次、异头立体描述全部错位。另一形态是糖环的起始位次与编号方向交给 `orient_numbering` 的「旋转/翻转候选 + locant_key 最小化」裁决，对多手性多臂糖环权衡失当。
- **IUPAC 依据**：P-102.7.2.2（有游离半缩醛基低聚糖的母体规则）、P-102.6.1.1.1（糖基残基 -yl/-yloxy）、P-102.3.4.2（异头中心 α/β 构型）、P-44.1.1 与 P-44.2（母体与环系优先级）、P-34（环系编号方向与起始位次）
- **责任模块**：src/namepredict/layer2/parent_selector.py:select_parent_tied；src/namepredict/layer4/numbering_engine.py:orient_numbering；src/namepredict/layer2/parent_skeleton.py:select_principal_skeletons

**代表样例**（序号为复杂度排名，越小的分子越简单）

- 346 chebi-2988 | N-[(3R,4R,5S,6R)-5-[(2S,3R,4R,5S,6R)-3-acetamido-5-[(2S,3S,4S,5R,6R)-4-…-6-[[…]oxymethyl]-3,5-dihydroxyoxan-2-yl]oxy-4-hydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-2,4-dihydroxy-6-(hydroxymethyl)oxan-3-yl]ace…
- 322 chebi-1840 | N-[(2S,3R,4R,5R,6R)-2-[(2R,3S,4R,5R,6R)-6-[(2R,3R,4R,5R,6S)-5-acetamido-2-(hydroxymethyl)-4-…-6-[(2R,3S,4R,5S,6S)-4,5,6-trihydroxy-2-(hydroxymethyl)oxan-3-yl]oxyoxan-3-yl]oxy-4,5-dihydroxy-2-(hydroxymet…
- 371 chebi-388 | ...(2R,3S,4S,5S,6S)-5-[(...)-氧基]-6-[[(...)-氧基]-...]... | 同一糖环上 5-/6- 两个支链互换（-6-[[...]-5-[(...)]）

**全部 8 个 id（按复杂度升序）**

```text
  chebi-2536 chebi-2747 chebi-1840 chebi-2988 chebi-649 chebi-2813 chebi-650 chebi-388
```

## 13. 母体候选为空（骨架表达失败）

**条数** 8（无输出 8）｜**出现分片** 05, 06

- **判定特征**：pred 两侧皆空且 reason=no_assemblable_candidate；L2 已枚举出骨架（`select_principal_skeletons` 返回非空），但每个骨架在表达层都返回 None，`_collect_candidates` 候选数为 0。分子特征是含未注册命名类的杂环（1,2,3-三唑、1,2,4-噻二唑、1,2,3-氧硫嗪二氧化物）或链上两个同种主官能团。
- **根因**：缺口在 L2「骨架 → 母体」的表达层而非枚举层：`principal_expression` 的表达表按（命名类 × 主官能团）组织，未注册环系不能被表达成可用 class（环内酮还会被 `_unsupported_typed_ring` 整体拦截），链上两个同种主官能团也无法收敛成单一母体骨架，两路都返回空列表 → L4/L5 无输入。属能力表缺项，不是覆盖或编号问题。
- **IUPAC 依据**：P-44.1 与 P-44.3（母体结构选择）、P-22.2.2（Hantzsch–Widman 杂环母体）、P-23 与 P-25（多环母体）
- **责任模块**：src/namepredict/layer2/principal_parent.py:rule_driven_parent_candidates；src/namepredict/layer2/principal_parent.py:_express_selected

**代表样例**（序号为复杂度排名，越小的分子越简单）

- 45 chebi-2761 | gold: 7-hydroxy-7-methyl-4-prop-1-en-2-yloxepan-2-one | pred: (空)
- 109 tiers-96746 | gold: 4-chloro-6-methoxy-7-[2-(1,2,3-triazol-1-yl)ethoxy]quinazoline | pred: (空)
- 149 chebi-1904 | gold: bis(2-ethylhexyl) hexanedioate | pred: (空)

**全部 8 个 id（按复杂度升序）**

```text
  chebi-318 chebi-2761 tiers-42304 tiers-96746 chebi-1904 tiers-98623 tiers-111140 chebi-1738
```

## 14. 酰胺/内酰胺互变异构未归一

**条数** 8（无输出 0）｜**出现分片** 01, 02, 05, 06

- **判定特征**：输入把酰胺写成亚氨酸/亚氨酸根形态（C([O-])=N、N=C(O)、N=C([O-])），或芳香环上的内酰胺位点（嘌呤/嘧啶的 2-氨基-6-酮 ↔ 2-亚氨基-6-羟基）；pred 保留 -OH 形式，出现 `1-hydroxyethylideneamino`、`6-hydroxy-2-imino-…-purin-9-yl`、`1-hydroxy-1-(methylimino)ethane` 这类亚氨酸名并平白多出一个 OH，gold 统一按酰胺/内酰胺（-amido、-one 后缀）。最小复现：[O-]C(C)=NC → 1-hydroxy-1-(methylimino)ethane（应为 N-甲基乙酰胺）。
- **根因**：L0 `tautomer.normalize_amide_tautomer` 只把「非芳香、中性、O 为隐氢羟基且 N 无显式 H」的 C(OH)=N 位点改成 C(=O)-NH（`_is_amide_enol_n` 明确要求非芳香环），以阴离子 [O-] 或显式 H 书写的亚氨酸位点被守卫条件跳过，芳香环上的互变位点也不处理；`charge.normalize_acid_charge` 同样未把酰胺氧阴离子还原为 C=O。残留的亚胺-醇式让 L1 把 N=C 判成亚胺、把 O⁻ 当羟基补回，净多一个 O-H；芳香内酰胺位点还会让主特征基团被识别为醇而非内酰胺，母体从环系退化为多元醇链（P-43 后缀优先级中酰胺应高于醇）。
- **IUPAC 依据**：P-66.1.1（酰胺与亚氨酸的命名、伯酰胺后缀）、P-66.1.1.4.3（-amido 保留式）、P-15.5.3.4（互变异构式的规范取舍）、P-43（后缀优先级）、P-25.1.1（嘌呤保留母体名）
- **责任模块**：src/namepredict/layer0/tautomer.py:normalize_amide_tautomer；src/namepredict/layer0/charge.py:normalize_acid_charge；src/namepredict/layer2/principal_expression.py:_NITROGEN_STEM_BY_FREE_DOUBLE

**代表样例**（序号为复杂度排名，越小的分子越简单）

- 313 chebi-886 | [(2S,3S,4S)-2-[[(2S,3R)-2,3-dihydroxytetracosanoyl]amino]-3,4-dihydroxyoctadecyl] [(2R,3S,5R,6R)-2,3,4,5,6-pentahydroxycyclohexyl] phosphate | (2R,3S,5R,6R)-2,3,4,5,6-pentahydroxycyclohexyl (2S,3S,4S)-3,…
- 353 chebi-1018 | [(2R,3R,4R,5S,6R)-3-acetamido-4-[(2S,…)-…oxan-2-yl]oxy-5-hydroxy-6-(hydroxymethyl)oxan-2-yl] [oxido-[(2Z,6Z,…)-…undecaenoxy]phosphoryl] phosphate | (2R,3S,4S,5S,6R)-2-[…-5-(1-hydroxyethylideneamino)-2-(…
- 369 chebi-1772 | (2S,4S,5R,6R)-5-acetamido-2-[…]oxy-…-4-hydroxyoxane-2-carboxylate | (2S,4S,5R,6R)-2-[…]-4-hydroxy-5-(1-hydroxyethylideneamino)oxane-2-carboxylate

**全部 8 个 id（按复杂度升序）**

```text
  chebi-3169 chebi-3261 chebi-886 chebi-1018 chebi-1144 chebi-1772 chebi-2877 chebi-1429
```

## 15. 磺酸/硫酸（氢）酯主基团未登记，降级为前缀

**条数** 6（无输出 0）｜**出现分片** 02, 06, 07

- **判定特征**：SMILES 含 -SO₂-OH（可带负电）或 -O-SO₂-OH；pred 把它写成 sulfo / sulfooxy / hydroxysulfonyloxy 前缀并把母体选成胺、吡啶或苯胺，gold 用 -sulfonic acid / hydrogen sulfate 后缀。同族症状：偶氮 -N=N- 被写成 `iminoamino` 而非 diazenyl。
- **根因**：`fg_registry.FG_SPECS` 未登记 sulfonic acid / hydrogen sulfate 主官能团（只在 `tools/anchored_table.py` 里以 sulfo 前缀存在），磺酸永远不能作后缀或母体，P-41 表 4.1 中第 7a 类（磺酸）高于第 11 类（酰胺）的优先序被倒置；L5 的 `_KIND_TABLE` 与 chain_engine 也没有对应后缀词干，阴离子端还叠加电荷丢失。`-N=N-` 桥没有 diazenyl 条目，取代基提取把它拆成两个氮前缀，经 assembler 的 imine→imino 映射拼出 iminoamino。
- **IUPAC 依据**：P-42.1 与 P-42.3（酸的后缀优先顺序，-sulfonic acid 属第 7a 类）、P-67.1.3.1 与 P-67.1.3.2（非碳酸含氧酸的盐与酯，如 methyl hydrogen sulfate）、P-68.3.1.3.2（diazenyl 保留前缀）、P-41 表 4.1（第 7a 类 > 第 11 类）
- **责任模块**：src/namepredict/layer1/fg_registry.py:FG_SPECS；src/namepredict/tools/anchored_table.py；src/namepredict/layer3/substituent_namer.py

**代表样例**（序号为复杂度排名，越小的分子越简单）

- 42 chebi-695 | 2-morpholin-4-ylethanesulfonic acid | 4-(2-sulfoethyl)morpholine
- 58 chebi-1733 | decyl sulfate | 1-(hydroxysulfonyloxy)decane
- 74 chebi-187 | 2-butyloctyl hydrogen sulfate | 5-(sulfooxymethyl)undecane

**全部 6 个 id（按复杂度升序）**

```text
  chebi-1495 chebi-695 chebi-1733 chebi-187 chebi-1982 chebi-1797
```

## 16. 大环母体/环几何模板缺失（含编号引擎 KeyError）

**条数** 5（无输出 5）｜**出现分片** 01, 02, 03, 07

- **判定特征**：pred 全空。分子含 ≥12 元（例 16/18/30/36/39 元）大环或大环内酯、环肽骨架；内部不是「无候选」而是 `number()` 抛 KeyError，栈落在 `fused_orientation._layout` 的 `RING_TEMPLATES[n]`。
- **根因**：L4 编号前要算环变形度，`_ring_deform`/`_layout` 直接索引 `RING_TEMPLATES[n]`，而该表只登记 3–19 元环；环尺寸 ≥20 时 KeyError，被 `_assemble_candidate` 的 `except (ValueError, KeyError, TypeError)` 静默吞掉，所有候选被丢弃而不冒泡、不降级、不报错。同组另一些分子（12–18 元环）是母体侧缺口：`ring_scaffold._TEMPLATES` 只覆盖 3–6 元单环与常见稠合母体，缺 >10 元环的骨架置换式（a 前缀 + cyclo 词干）生成器。
- **IUPAC 依据**：P-22.2.4（≥11 元 mancude 杂单环）、P-22.2.3（骨架置换 a 命名杂单环）、P-23.2（大环 von Baeyer 命名）、P-25.3.3（环系编号）
- **责任模块**：src/namepredict/layer4/fused_orientation.py:_layout（入口 src/namepredict/layer4/numbering.py:number）；src/namepredict/layer2/ring_scaffold.py:_TEMPLATES

**代表样例**（序号为复杂度排名，越小的分子越简单）

- 361 chebi-3272 | (1S,6R,8R,...)-8,17,26,35,44,53-hexakis(6-aminopurin-9-yl)-3,12,21,30,39,48-hexaoxido-…-hexaphosphaheptacyclo[49.3.0.06,10.015,19.024,28.033,37.042,46]tetrapentacontane-9,18,27,36,45,54-hexol | (空)
- 218 chebi-1243 | (1S,11R,13S,14S,24R,26S)-13,26-dimethyl-2,15-dioxa-12,25-diazatricyclo[22.2.2.211,14]triacontane-3,16-dione | (空)
- 315 chebi-31 | (1S,15S,16R,17R,18S,19E,21E,25E,27E,29E,31E)-33-[(2S,3S,4S,5S,6R)-4-amino-3,5-dihydroxy-6-methyloxan-2-yl]oxy-1,3,4,7,9,11,17,37-octahydroxy-15,16,18-trimethyl-13-oxo-14,39-dioxabicyclo[33.3.1]nonatriacon…

**全部 5 个 id（按复杂度升序）**

```text
  chebi-1243 chebi-690 chebi-31 chebi-439 chebi-3272
```

## 尾类（零散缺陷）

条数少但机制独立的缺陷，逐条列出全部 id（按复杂度升序）。

### T17. 环上取代基位次编号方向错误（2 条）

- **判定特征**：pred 与 gold 是同骨架、同取代基集合，但环上取代基 locant 整体互换（糖环 2↔6、2↔3 对调），位次集合相同或高度重叠，立体描述符随之整体错位。
- **根因**：L4 `orient_numbering` 在多方向编号候选的收窄中，对「位次集合并列」的等价方向缺少与 IUPAC 一致的后续裁决准则（环上取代基按字母序/CIP 或异头碳优先取最低位次），选了与 gold 相反的方向；位次整体镜像后名称内部仍自洽，不触发任何失败门控，只在比对时表现为名称整体不同。
- **IUPAC 依据**：P-34.2 与 P-34.4（最低位次与并列情形的裁决）、P-102.3（糖类的环编号与异头位）
- **责任模块**：src/namepredict/layer4/numbering_engine.py:orient_numbering
- **全部 id**：chebi-225、chebi-2396

### T18. G gold 存疑（1 条）

- **判定特征**：pred 与 gold 骨架一致但母体名体系不同：gold 用 von Baeyer 名（pentacyclo[11.8.0.0³,¹¹.0⁴,⁹.0¹⁵,²⁰]henicosa），pred 用稠合名 benzo[b]naphtho[3,2-f]benzothiophene；分子 21 个环原子、5 环，pred 的稠合算式 9+4+8=21 与分子相符。
- **根因**：不是命名器能力缺陷：ChEBI 对全稠合多环芳烃采用 von Baeyer 形式，而 P-25.3 对邻稠合体系优先稠合命名，两者为等价合法名，difflib 字符串相似度无法识别等价名。另有若干条目的 gold_zh 存在多环符号直译损坏（gold_en 正确），未单列。
- **IUPAC 依据**：P-25.3 与 P-23.2（邻稠合体系优先稠合命名，von Baeyer 用于不可邻稠合者）
- **责任模块**：未定位（gold 来源与评分口径问题，非 src 缺陷）
- **全部 id**：chebi-1639

### T19. 取代基前缀字母序与保留前缀未命中（1 条）

- **判定特征**：pred 与 gold 描述同一结构、母体与位次相同，差异只在取代基前缀写法与引用顺序：gold 用保留前缀 anilino/苯胺基，pred 退回系统式 phenylamino/苯氨基；且复杂取代基与简单取代基的字母序颠倒。
- **根因**：`tools/anchored_table.py` 的 anilino 保留叶未命中，取代基命名退回 free-name 系统式；`assembler_prefixes._sorted_stems` 的字母序比较键未剔除位次与 N-斜体等不参与排序的成分，也未按取代基名首字母比较，前缀引用顺序与 P-14.5 不符。
- **IUPAC 依据**：P-14.5 与 P-14.5.1（前缀按字母（数字）序引用）、P-62.2.1.1.2（anilino 保留前缀优先于 phenylamino）
- **责任模块**：src/namepredict/layer5/assembler_prefixes.py:_sorted_stems；src/namepredict/tools/anchored_table.py
- **全部 id**：tiers-150173

### T20. 多组分体系只命名单一组分（1 条）

- **判定特征**：SMILES 含多个以 `.` 分隔的共价组分（非金属盐、非盐酸盐），pred 只有其中一个组分名（通常是优先级最高官能团所在的小组分），gold 按组分顺序列出全部。
- **根因**：多片段输入不走 salt 拆分路径：主官能团选择在整分子（含所有片段）上只挑出优先级最高者，母体骨架只从该组分生成，其余组分既不能作母体也不进入取代基体系，最终被当作 coverage 缺口丢弃。
- **IUPAC 依据**：P-24（多组分/多片段体系各组分须分别命名并按序排列）
- **责任模块**：src/namepredict/layer2/principal_parent.py:rule_driven_parent_candidates；src/namepredict/layer0/salt.py:dissociate_salt
- **全部 id**：chebi-2034

### T21. 寡核苷酸链无命名支持（1 条）

- **判定特征**：多核苷酸（SMILES 含多个 OP(=O)(O)OC 桥与多种核苷酸碱基，重原子数百），pred 为空；对照判据：单核苷酸（dAMP）、二核苷酸（d(ApA)）可正常命名，说明不是碱基、糖或磷酸本身的问题。
- **根因**：磷酸二酯桥既没有核酸链专用命名（P-102 糖基化规则在 >2 单元时不适用），也超出 `tools/anchored_table.py` 单磷酸前缀的覆盖；母体只能选一个核苷单元，其余单元作为取代基层层嵌套时在某处失败（coverage 不完整），而 `namer._run_candidates` 只有单一候选阶段、无降级路径，直接返回 no_assemblable_candidate。
- **IUPAC 依据**：P-102.6（糖苷与糖基残基）、P-67.1.3（磷酸的酯）、P-14.4（位次必需引用）
- **责任模块**：src/namepredict/namer.py:_run_candidates；src/namepredict/tools/anchored_table.py
- **全部 id**：chebi-2969

### T22. 稠环氢化程度/指示氢表达错误（1 条）

- **判定特征**：pred 与 gold 的母体骨架及甲基/羟基/糖取代位次一致，但 hydro 前缀的加氢位次与数量不同（dodecahydro-1H vs decahydro-1H,2H,13H），即两个名称描述的氢化程度与不饱和位点不同，相差 2 个 H。
- **根因**：L4 推导加氢原子集 `hydro_atoms` 的来源（保留模板比对，缺失时由 `_fallback_hydro_atoms` 从分子自身饱和环位反推）与 gold 采用的不饱和位点判定不一致；`hydro_prefix` 在加氢计数不合法或位次不在链内时把整套 hydro 位次退回指示氢（P-58.2.1），未真正加氢的环碳被写成 2H/13H，名称的氢化程度与输入结构不符。
- **IUPAC 依据**：P-31.2.2（加氢程度前缀 hydro）、P-58.2.1（指示氢）、P-14.3.4.5（完全氢化省略位次）
- **责任模块**：src/namepredict/layer4/numbering.py:number
- **全部 id**：chebi-2171

### T23. 金属有机/配位化合物无母体（1 条）

- **判定特征**：pred 为空；SMILES 含过渡金属原子（Ti/Fe/Mo…）或 π-配位碳环阴离子（如 [cH-]1[cH-][cH-][cH-][cH-]1），gold 为夹心/配位型名称。
- **根因**：金属原子既不是骨架碳也不在任何官能团清单里，配位的环戊二烯基阴离子也不构成可命名母体；主官能团为空（NONE）时环骨架走 `_generic_ring_kind`，芳香但未注册的单环在 n_rings<2 时返回 None → 候选集为空。缺少 P-69 的配位/夹心型（titanocene 类）母体生成路径。
- **IUPAC 依据**：P-69.2.4（对碳的多中心键有机金属基团）、P-69.2.6（不饱和分子配体）
- **责任模块**：src/namepredict/layer2/principal_expression.py:_generic_ring_kind
- **全部 id**：chebi-1246

## 附录：最差批次全量索引（按复杂度升序）

| # | id | 重原子 | 环 | 环系统 | sim_max | 输出 | 归属模式 |
|---:|---|---:|---:|---:|---:|---|---|
| 1 | chebi-1152 | 4 | 0 | 0 | 0.32 | 有 | 离子/电荷/反离子未体现（季铵、阴离子、叶立德） |
| 2 | chebi-2055 | 5 | 0 | 0 | 0.26 | 有 | 离子/电荷/反离子未体现（季铵、阴离子、叶立德） |
| 3 | tiers-9441 | 6 | 0 | 0 | 0.00 | 空 | 链式同种主基团多重性不可表达（二酯/二酰胺/酸酐/偕二醇） |
| 4 | tiers-93 | 7 | 0 | 0 | 0.00 | 空 | 链式同种主基团多重性不可表达（二酯/二酰胺/酸酐/偕二醇） |
| 5 | chebi-2628 | 7 | 0 | 0 | 0.33 | 有 | 离子/电荷/反离子未体现（季铵、阴离子、叶立德） |
| 6 | chebi-709 | 7 | 0 | 0 | 0.28 | 有 | 母体选择与主基优先级错误 |
| 7 | chebi-1282 | 7 | 0 | 0 | 0.26 | 有 | 离子/电荷/反离子未体现（季铵、阴离子、叶立德） |
| 8 | chebi-2334 | 7 | 0 | 0 | 0.33 | 有 | coverage 门控放宽：残缺局部片段名被当作成功输出 |
| 9 | chebi-2752 | 7 | 1 | 1 | 0.00 | 空 | 未注册单环/杂环母体缺失（HW 模板与保留名注册表缺项） |
| 10 | chebi-2007 | 8 | 0 | 0 | 0.39 | 有 | 母体选择与主基优先级错误 |
| 11 | tiers-5188 | 8 | 1 | 1 | 0.00 | 空 | 未注册单环/杂环母体缺失（HW 模板与保留名注册表缺项） |
| 12 | tiers-180 | 9 | 0 | 0 | 0.00 | 空 | 链式同种主基团多重性不可表达（二酯/二酰胺/酸酐/偕二醇） |
| 13 | chebi-895 | 9 | 0 | 0 | 0.36 | 有 | 离子/电荷/反离子未体现（季铵、阴离子、叶立德） |
| 14 | chebi-3285 | 9 | 1 | 1 | 0.24 | 有 | coverage 门控放宽：残缺局部片段名被当作成功输出 |
| 15 | chebi-1786 | 9 | 2 | 1 | 0.00 | 空 | 桥环/螺环/笼状母体缺 von Baeyer 词干 |
| 16 | chebi-971 | 9 | 3 | 1 | 0.00 | 空 | 桥环/螺环/笼状母体缺 von Baeyer 词干 |
| 17 | chebi-2703 | 10 | 0 | 0 | 0.00 | 空 | 链式同种主基团多重性不可表达（二酯/二酰胺/酸酐/偕二醇） |
| 18 | tiers-186 | 10 | 0 | 0 | 0.00 | 空 | 链式同种主基团多重性不可表达（二酯/二酰胺/酸酐/偕二醇） |
| 19 | tiers-61915 | 10 | 0 | 0 | 0.00 | 空 | 链式同种主基团多重性不可表达（二酯/二酰胺/酸酐/偕二醇） |
| 20 | tiers-43418 | 10 | 0 | 0 | 0.30 | 有 | 杂原子中心（P/B）被当作碳命名 |
| 21 | chebi-1907 | 10 | 1 | 1 | 0.39 | 有 | 取代基递归命名残缺（块内原子静默丢弃） |
| 22 | chebi-318 | 10 | 1 | 1 | 0.00 | 空 | 母体候选为空（骨架表达失败） |
| 23 | tiers-58291 | 10 | 1 | 1 | 0.23 | 有 | 未注册单环/杂环母体缺失（HW 模板与保留名注册表缺项） |
| 24 | chebi-826 | 10 | 2 | 1 | 0.00 | 空 | 桥环/螺环/笼状母体缺 von Baeyer 词干 |
| 25 | chebi-2415 | 10 | 2 | 1 | 0.00 | 空 | 多环母体词干组装失败（稠合/桥/螺未细分） |
| 26 | chebi-885 | 10 | 2 | 1 | 0.00 | 空 | 多环母体词干组装失败（稠合/桥/螺未细分） |
| 27 | chebi-3137 | 10 | 3 | 1 | 0.00 | 空 | 稠合环系拆解退化/稠合词干失败 |
| 28 | chebi-903 | 10 | 3 | 1 | 0.00 | 空 | 桥环/螺环/笼状母体缺 von Baeyer 词干 |
| 29 | tiers-66561 | 10 | 3 | 1 | 0.00 | 空 | 多环母体词干组装失败（稠合/桥/螺未细分） |
| 30 | chebi-1495 | 11 | 1 | 1 | 0.35 | 有 | 磺酸/硫酸（氢）酯主基团未登记，降级为前缀 |
| 31 | chebi-2147 | 11 | 2 | 1 | 0.00 | 空 | 桥环/螺环/笼状母体缺 von Baeyer 词干 |
| 32 | chebi-2153 | 11 | 2 | 1 | 0.00 | 空 | 桥环/螺环/笼状母体缺 von Baeyer 词干 |
| 33 | tiers-44245 | 11 | 2 | 1 | 0.00 | 空 | 多环母体词干组装失败（稠合/桥/螺未细分） |
| 34 | chebi-3151 | 11 | 2 | 1 | 0.00 | 空 | 多环母体词干组装失败（稠合/桥/螺未细分） |
| 35 | chebi-1503 | 11 | 2 | 1 | 0.00 | 空 | 稠合环系拆解退化/稠合词干失败 |
| 36 | chebi-1246 | 11 | 2 | 2 | 0.00 | 空 | 金属有机/配位化合物无母体 |
| 37 | tiers-29902 | 11 | 3 | 1 | 0.26 | 有 | 多环骨架母体名生成错误（伪稠合/拓扑不符/位次缺失） |
| 38 | chebi-983 | 12 | 0 | 0 | 0.00 | 空 | 母体选择与主基优先级错误 |
| 39 | chebi-3086 | 12 | 1 | 1 | 0.22 | 有 | 离子/电荷/反离子未体现（季铵、阴离子、叶立德） |
| 40 | chebi-1403 | 12 | 1 | 1 | 0.00 | 空 | 未注册单环/杂环母体缺失（HW 模板与保留名注册表缺项） |
| 41 | chebi-402 | 12 | 1 | 1 | 0.38 | 有 | 离子/电荷/反离子未体现（季铵、阴离子、叶立德） |
| 42 | chebi-695 | 12 | 1 | 1 | 0.39 | 有 | 磺酸/硫酸（氢）酯主基团未登记，降级为前缀 |
| 43 | tiers-65526 | 12 | 2 | 1 | 0.00 | 空 | 稠合环系拆解退化/稠合词干失败 |
| 44 | tiers-6342 | 12 | 2 | 1 | 0.00 | 空 | 桥环/螺环/笼状母体缺 von Baeyer 词干 |
| 45 | chebi-2761 | 13 | 1 | 1 | 0.00 | 空 | 母体候选为空（骨架表达失败） |
| 46 | tiers-17996 | 13 | 2 | 1 | 0.00 | 空 | 多环母体词干组装失败（稠合/桥/螺未细分） |
| 47 | tiers-85289 | 13 | 2 | 1 | 0.00 | 空 | 多环母体词干组装失败（稠合/桥/螺未细分） |
| 48 | chebi-1919 | 13 | 2 | 1 | 0.00 | 空 | 稠合环系拆解退化/稠合词干失败 |
| 49 | chebi-1796 | 13 | 2 | 1 | 0.00 | 空 | 多环母体词干组装失败（稠合/桥/螺未细分） |
| 50 | chebi-82 | 13 | 2 | 1 | 0.00 | 空 | 多环母体词干组装失败（稠合/桥/螺未细分） |
| 51 | tiers-59592 | 13 | 2 | 1 | 0.00 | 空 | 稠合环系拆解退化/稠合词干失败 |
| 52 | chebi-395 | 13 | 4 | 1 | 0.21 | 有 | 多环骨架母体名生成错误（伪稠合/拓扑不符/位次缺失） |
| 53 | chebi-1695 | 14 | 0 | 0 | 0.27 | 有 | 取代基递归命名残缺（块内原子静默丢弃） |
| 54 | chebi-3281 | 14 | 0 | 0 | 0.39 | 有 | 母体选择与主基优先级错误 |
| 55 | chebi-3066 | 14 | 1 | 1 | 0.00 | 空 | 未注册单环/杂环母体缺失（HW 模板与保留名注册表缺项） |
| 56 | chebi-531 | 14 | 1 | 1 | 0.00 | 空 | 未注册单环/杂环母体缺失（HW 模板与保留名注册表缺项） |
| 57 | chebi-1147 | 14 | 2 | 1 | 0.00 | 空 | 多环母体词干组装失败（稠合/桥/螺未细分） |
| 58 | chebi-1733 | 15 | 0 | 0 | 0.29 | 有 | 磺酸/硫酸（氢）酯主基团未登记，降级为前缀 |
| 59 | tiers-8543 | 15 | 2 | 2 | 0.32 | 有 | 杂原子中心（P/B）被当作碳命名 |
| 60 | chebi-384 | 15 | 2 | 1 | 0.37 | 有 | 离子/电荷/反离子未体现（季铵、阴离子、叶立德） |
| 61 | chebi-2221 | 15 | 2 | 1 | 0.00 | 空 | 多环母体词干组装失败（稠合/桥/螺未细分） |
| 62 | chebi-1914 | 15 | 2 | 1 | 0.00 | 空 | 母体选择与主基优先级错误 |
| 63 | chebi-1615 | 15 | 2 | 2 | 0.00 | 空 | 未注册单环/杂环母体缺失（HW 模板与保留名注册表缺项） |
| 64 | tiers-96844 | 15 | 3 | 1 | 0.00 | 空 | 稠合环系拆解退化/稠合词干失败 |
| 65 | chebi-1256 | 15 | 3 | 1 | 0.00 | 空 | 多环母体词干组装失败（稠合/桥/螺未细分） |
| 66 | tiers-1630 | 16 | 0 | 0 | 0.00 | 空 | 链式同种主基团多重性不可表达（二酯/二酰胺/酸酐/偕二醇） |
| 67 | chebi-432 | 16 | 1 | 1 | 0.00 | 空 | 未注册单环/杂环母体缺失（HW 模板与保留名注册表缺项） |
| 68 | tiers-7355 | 16 | 2 | 1 | 0.00 | 空 | 稠合环系拆解退化/稠合词干失败 |
| 69 | tiers-80645 | 16 | 2 | 1 | 0.18 | 有 | coverage 门控放宽：残缺局部片段名被当作成功输出 |
| 70 | tiers-42304 | 16 | 2 | 2 | 0.00 | 空 | 母体候选为空（骨架表达失败） |
| 71 | chebi-194 | 16 | 3 | 1 | 0.34 | 有 | 多环骨架母体名生成错误（伪稠合/拓扑不符/位次缺失） |
| 72 | tiers-18169 | 16 | 3 | 1 | 0.00 | 空 | 稠合环系拆解退化/稠合词干失败 |
| 73 | chebi-868 | 16 | 3 | 3 | 0.00 | 空 | 未注册单环/杂环母体缺失（HW 模板与保留名注册表缺项） |
| 74 | chebi-187 | 17 | 0 | 0 | 0.27 | 有 | 磺酸/硫酸（氢）酯主基团未登记，降级为前缀 |
| 75 | chebi-298 | 17 | 2 | 1 | 0.00 | 空 | 稠合环系拆解退化/稠合词干失败 |
| 76 | chebi-846 | 17 | 2 | 1 | 0.00 | 空 | 稠合环系拆解退化/稠合词干失败 |
| 77 | tiers-68910 | 17 | 2 | 1 | 0.00 | 空 | 多环母体词干组装失败（稠合/桥/螺未细分） |
| 78 | chebi-1283 | 17 | 2 | 2 | 0.38 | 有 | 取代基递归命名残缺（块内原子静默丢弃） |
| 79 | chebi-474 | 17 | 3 | 2 | 0.00 | 空 | 多环母体词干组装失败（稠合/桥/螺未细分） |
| 80 | tiers-124039 | 17 | 3 | 1 | 0.00 | 空 | 稠合环系拆解退化/稠合词干失败 |
| 81 | tiers-53012 | 17 | 3 | 3 | 0.29 | 有 | coverage 门控放宽：残缺局部片段名被当作成功输出 |
| 82 | chebi-1053 | 17 | 3 | 1 | 0.00 | 空 | 多环母体词干组装失败（稠合/桥/螺未细分） |
| 83 | chebi-1236 | 17 | 4 | 1 | 0.36 | 有 | 多环骨架母体名生成错误（伪稠合/拓扑不符/位次缺失） |
| 84 | chebi-1520 | 17 | 4 | 1 | 0.00 | 空 | 桥环/螺环/笼状母体缺 von Baeyer 词干 |
| 85 | chebi-2517 | 18 | 3 | 1 | 0.00 | 空 | 多环母体词干组装失败（稠合/桥/螺未细分） |
| 86 | chebi-1703 | 18 | 3 | 1 | 0.00 | 空 | 多环母体词干组装失败（稠合/桥/螺未细分） |
| 87 | chebi-93 | 18 | 3 | 3 | 0.00 | 空 | 未注册单环/杂环母体缺失（HW 模板与保留名注册表缺项） |
| 88 | chebi-1459 | 19 | 2 | 1 | 0.00 | 空 | 稠合环系拆解退化/稠合词干失败 |
| 89 | chebi-1081 | 19 | 2 | 1 | 0.00 | 空 | 多环母体词干组装失败（稠合/桥/螺未细分） |
| 90 | chebi-1338 | 19 | 3 | 2 | 0.35 | 有 | coverage 门控放宽：残缺局部片段名被当作成功输出 |
| 91 | tiers-30123 | 19 | 4 | 2 | 0.00 | 空 | 稠合环系拆解退化/稠合词干失败 |
| 92 | chebi-862 | 20 | 2 | 1 | 0.00 | 空 | 链式同种主基团多重性不可表达（二酯/二酰胺/酸酐/偕二醇） |
| 93 | chebi-2034 | 20 | 2 | 2 | 0.35 | 有 | 多组分体系只命名单一组分 |
| 94 | chebi-3169 | 20 | 2 | 1 | 0.40 | 有 | 酰胺/内酰胺互变异构未归一 |
| 95 | chebi-2616 | 20 | 3 | 1 | 0.39 | 有 | coverage 门控放宽：残缺局部片段名被当作成功输出 |
| 96 | chebi-1864 | 20 | 3 | 2 | 0.00 | 空 | 稠合环系拆解退化/稠合词干失败 |
| 97 | chebi-1588 | 20 | 3 | 2 | 0.00 | 空 | 多环母体词干组装失败（稠合/桥/螺未细分） |
| 98 | chebi-164 | 20 | 3 | 3 | 0.00 | 空 | 未注册单环/杂环母体缺失（HW 模板与保留名注册表缺项） |
| 99 | chebi-2136 | 20 | 3 | 2 | 0.25 | 有 | coverage 门控放宽：残缺局部片段名被当作成功输出 |
| 100 | chebi-1979 | 20 | 3 | 2 | 0.31 | 有 | coverage 门控放宽：残缺局部片段名被当作成功输出 |
| 101 | chebi-1727 | 20 | 4 | 1 | 0.30 | 有 | 多环骨架母体名生成错误（伪稠合/拓扑不符/位次缺失） |
| 102 | chebi-2811 | 20 | 4 | 1 | 0.00 | 空 | 多环母体词干组装失败（稠合/桥/螺未细分） |
| 103 | tiers-129945 | 20 | 4 | 2 | 0.00 | 空 | 多环母体词干组装失败（稠合/桥/螺未细分） |
| 104 | chebi-594 | 20 | 5 | 1 | 0.00 | 空 | 桥环/螺环/笼状母体缺 von Baeyer 词干 |
| 105 | chebi-2477 | 20 | 5 | 1 | 0.00 | 空 | 多环母体词干组装失败（稠合/桥/螺未细分） |
| 106 | tiers-26894 | 21 | 2 | 2 | 0.00 | 空 | 未注册单环/杂环母体缺失（HW 模板与保留名注册表缺项） |
| 107 | chebi-3280 | 21 | 3 | 1 | 0.21 | 有 | 离子/电荷/反离子未体现（季铵、阴离子、叶立德） |
| 108 | chebi-1560 | 21 | 3 | 1 | 0.00 | 空 | 稠合环系拆解退化/稠合词干失败 |
| 109 | tiers-96746 | 21 | 3 | 2 | 0.00 | 空 | 母体候选为空（骨架表达失败） |
| 110 | chebi-2576 | 21 | 5 | 1 | 0.30 | 有 | 多环骨架母体名生成错误（伪稠合/拓扑不符/位次缺失） |
| 111 | chebi-1639 | 21 | 5 | 1 | 0.15 | 有 | G gold 存疑 |
| 112 | chebi-1438 | 21 | 5 | 1 | 0.00 | 空 | 桥环/螺环/笼状母体缺 von Baeyer 词干 |
| 113 | chebi-2594 | 21 | 6 | 1 | 0.19 | 有 | 多环骨架母体名生成错误（伪稠合/拓扑不符/位次缺失） |
| 114 | chebi-3297 | 22 | 3 | 1 | 0.00 | 空 | 多环母体词干组装失败（稠合/桥/螺未细分） |
| 115 | chebi-2696 | 22 | 3 | 2 | 0.00 | 空 | 稠合环系拆解退化/稠合词干失败 |
| 116 | chebi-503 | 22 | 3 | 1 | 0.00 | 空 | 桥环/螺环/笼状母体缺 von Baeyer 词干 |
| 117 | chebi-1872 | 22 | 5 | 2 | 0.00 | 空 | 多环母体词干组装失败（稠合/桥/螺未细分） |
| 118 | chebi-1982 | 23 | 2 | 2 | 0.34 | 有 | 磺酸/硫酸（氢）酯主基团未登记，降级为前缀 |
| 119 | chebi-1684 | 23 | 2 | 2 | 0.00 | 空 | 离子/电荷/反离子未体现（季铵、阴离子、叶立德） |
| 120 | chebi-263 | 23 | 3 | 2 | 0.00 | 空 | 稠合环系拆解退化/稠合词干失败 |
| 121 | chebi-403 | 23 | 3 | 2 | 0.00 | 空 | 多环母体词干组装失败（稠合/桥/螺未细分） |
| 122 | chebi-2143 | 23 | 4 | 2 | 0.00 | 空 | 多环母体词干组装失败（稠合/桥/螺未细分） |
| 123 | chebi-967 | 23 | 4 | 1 | 0.00 | 空 | 稠合环系拆解退化/稠合词干失败 |
| 124 | chebi-2303 | 23 | 4 | 1 | 0.00 | 空 | 桥环/螺环/笼状母体缺 von Baeyer 词干 |
| 125 | chebi-2427 | 23 | 4 | 2 | 0.25 | 有 | coverage 门控放宽：残缺局部片段名被当作成功输出 |
| 126 | chebi-342 | 23 | 4 | 2 | 0.00 | 空 | 多环母体词干组装失败（稠合/桥/螺未细分） |
| 127 | chebi-1356 | 23 | 4 | 2 | 0.00 | 空 | 多环母体词干组装失败（稠合/桥/螺未细分） |
| 128 | chebi-226 | 23 | 5 | 1 | 0.00 | 空 | 桥环/螺环/笼状母体缺 von Baeyer 词干 |
| 129 | tiers-128583 | 24 | 3 | 2 | 0.00 | 空 | 链式同种主基团多重性不可表达（二酯/二酰胺/酸酐/偕二醇） |
| 130 | chebi-2311 | 24 | 3 | 2 | 0.00 | 空 | 多环母体词干组装失败（稠合/桥/螺未细分） |
| 131 | chebi-2409 | 24 | 3 | 2 | 0.00 | 空 | 稠合环系拆解退化/稠合词干失败 |
| 132 | chebi-2257 | 24 | 4 | 1 | 0.38 | 有 | 稠合环系拆解退化/稠合词干失败 |
| 133 | chebi-1484 | 24 | 4 | 1 | 0.00 | 空 | 多环母体词干组装失败（稠合/桥/螺未细分） |
| 134 | chebi-1908 | 24 | 4 | 1 | 0.00 | 空 | 多环母体词干组装失败（稠合/桥/螺未细分） |
| 135 | chebi-3257 | 24 | 4 | 3 | 0.00 | 空 | 未注册单环/杂环母体缺失（HW 模板与保留名注册表缺项） |
| 136 | chebi-2327 | 24 | 4 | 2 | 0.40 | 有 | coverage 门控放宽：残缺局部片段名被当作成功输出 |
| 137 | chebi-2406 | 24 | 5 | 1 | 0.00 | 空 | 多环母体词干组装失败（稠合/桥/螺未细分） |
| 138 | chebi-1196 | 25 | 2 | 1 | 0.24 | 有 | coverage 门控放宽：残缺局部片段名被当作成功输出 |
| 139 | chebi-289 | 25 | 3 | 2 | 0.00 | 空 | 稠合环系拆解退化/稠合词干失败 |
| 140 | chebi-604 | 25 | 3 | 2 | 0.00 | 空 | 稠合环系拆解退化/稠合词干失败 |
| 141 | chebi-3266 | 25 | 4 | 1 | 0.00 | 空 | 多环母体词干组装失败（稠合/桥/螺未细分） |
| 142 | chebi-3078 | 25 | 4 | 1 | 0.00 | 空 | 多环母体词干组装失败（稠合/桥/螺未细分） |
| 143 | chebi-1305 | 25 | 4 | 1 | 0.00 | 空 | 多环母体词干组装失败（稠合/桥/螺未细分） |
| 144 | chebi-2589 | 25 | 4 | 2 | 0.00 | 空 | 稠合环系拆解退化/稠合词干失败 |
| 145 | chebi-434 | 25 | 4 | 1 | 0.00 | 空 | 多环母体词干组装失败（稠合/桥/螺未细分） |
| 146 | chebi-606 | 25 | 4 | 2 | 0.00 | 空 | 多环母体词干组装失败（稠合/桥/螺未细分） |
| 147 | chebi-1567 | 25 | 4 | 2 | 0.00 | 空 | 稠合环系拆解退化/稠合词干失败 |
| 148 | chebi-2537 | 25 | 5 | 1 | 0.00 | 空 | 桥环/螺环/笼状母体缺 von Baeyer 词干 |
| 149 | chebi-1904 | 26 | 0 | 0 | 0.00 | 空 | 母体候选为空（骨架表达失败） |
| 150 | chebi-2242 | 26 | 2 | 1 | 0.18 | 有 | coverage 门控放宽：残缺局部片段名被当作成功输出 |
| 151 | chebi-2022 | 26 | 3 | 1 | 0.00 | 空 | 多环母体词干组装失败（稠合/桥/螺未细分） |
| 152 | chebi-2283 | 26 | 3 | 1 | 0.00 | 空 | 稠合环系拆解退化/稠合词干失败 |
| 153 | chebi-424 | 26 | 3 | 1 | 0.32 | 有 | coverage 门控放宽：残缺局部片段名被当作成功输出 |
| 154 | chebi-2999 | 26 | 4 | 1 | 0.00 | 空 | 多环母体词干组装失败（稠合/桥/螺未细分） |
| 155 | chebi-1417 | 26 | 4 | 1 | 0.00 | 空 | 稠合环系拆解退化/稠合词干失败 |
| 156 | chebi-479 | 26 | 4 | 1 | 0.00 | 空 | 桥环/螺环/笼状母体缺 von Baeyer 词干 |
| 157 | chebi-1698 | 26 | 5 | 1 | 0.36 | 有 | 多环骨架母体名生成错误（伪稠合/拓扑不符/位次缺失） |
| 158 | chebi-2061 | 26 | 5 | 1 | 0.38 | 有 | 多环骨架母体名生成错误（伪稠合/拓扑不符/位次缺失） |
| 159 | chebi-2041 | 26 | 6 | 1 | 0.00 | 空 | 多环母体词干组装失败（稠合/桥/螺未细分） |
| 160 | chebi-3288 | 26 | 6 | 1 | 0.00 | 空 | 桥环/螺环/笼状母体缺 von Baeyer 词干 |
| 161 | chebi-9 | 27 | 0 | 0 | 0.00 | 空 | 链式同种主基团多重性不可表达（二酯/二酰胺/酸酐/偕二醇） |
| 162 | chebi-1389 | 27 | 2 | 1 | 0.00 | 空 | 多环母体词干组装失败（稠合/桥/螺未细分） |
| 163 | chebi-2938 | 27 | 2 | 1 | 0.21 | 有 | coverage 门控放宽：残缺局部片段名被当作成功输出 |
| 164 | chebi-162 | 27 | 3 | 2 | 0.00 | 空 | 稠合环系拆解退化/稠合词干失败 |
| 165 | chebi-1248 | 27 | 3 | 2 | 0.00 | 空 | 多环母体词干组装失败（稠合/桥/螺未细分） |
| 166 | chebi-3085 | 27 | 3 | 2 | 0.24 | 有 | coverage 门控放宽：残缺局部片段名被当作成功输出 |
| 167 | chebi-1086 | 27 | 3 | 2 | 0.00 | 空 | 多环母体词干组装失败（稠合/桥/螺未细分） |
| 168 | chebi-571 | 27 | 3 | 2 | 0.40 | 有 | 杂原子中心（P/B）被当作碳命名 |
| 169 | tiers-91385 | 27 | 4 | 2 | 0.00 | 空 | 链式同种主基团多重性不可表达（二酯/二酰胺/酸酐/偕二醇） |
| 170 | tiers-118482 | 27 | 4 | 3 | 0.31 | 有 | 母体选择与主基优先级错误 |
| 171 | tiers-107887 | 27 | 4 | 2 | 0.00 | 空 | 链式同种主基团多重性不可表达（二酯/二酰胺/酸酐/偕二醇） |
| 172 | chebi-2674 | 27 | 5 | 1 | 0.22 | 有 | coverage 门控放宽：残缺局部片段名被当作成功输出 |
| 173 | chebi-822 | 27 | 7 | 1 | 0.38 | 有 | 多环骨架母体名生成错误（伪稠合/拓扑不符/位次缺失） |
| 174 | chebi-2351 | 28 | 1 | 1 | 0.14 | 有 | coverage 门控放宽：残缺局部片段名被当作成功输出 |
| 175 | tiers-129440 | 28 | 2 | 2 | 0.00 | 空 | 链式同种主基团多重性不可表达（二酯/二酰胺/酸酐/偕二醇） |
| 176 | chebi-369 | 28 | 2 | 2 | 0.38 | 有 | 杂原子中心（P/B）被当作碳命名 |
| 177 | chebi-129 | 28 | 4 | 1 | 0.00 | 空 | 多环母体词干组装失败（稠合/桥/螺未细分） |
| 178 | chebi-1692 | 28 | 4 | 2 | 0.32 | 有 | 稠合环系拆解退化/稠合词干失败 |
| 179 | chebi-2707 | 28 | 5 | 1 | 0.39 | 有 | 离子/电荷/反离子未体现（季铵、阴离子、叶立德） |
| 180 | tiers-95650 | 28 | 6 | 3 | 0.00 | 空 | 桥环/螺环/笼状母体缺 von Baeyer 词干 |
| 181 | chebi-197 | 28 | 7 | 1 | 0.40 | 有 | 多环骨架母体名生成错误（伪稠合/拓扑不符/位次缺失） |
| 182 | tiers-98623 | 29 | 3 | 3 | 0.00 | 空 | 母体候选为空（骨架表达失败） |
| 183 | chebi-2128 | 29 | 3 | 2 | 0.13 | 有 | coverage 门控放宽：残缺局部片段名被当作成功输出 |
| 184 | chebi-676 | 29 | 4 | 3 | 0.00 | 空 | 桥环/螺环/笼状母体缺 von Baeyer 词干 |
| 185 | chebi-84 | 29 | 4 | 3 | 0.00 | 空 | 多环母体词干组装失败（稠合/桥/螺未细分） |
| 186 | chebi-1240 | 29 | 5 | 1 | 0.00 | 空 | 多环母体词干组装失败（稠合/桥/螺未细分） |
| 187 | chebi-648 | 30 | 0 | 0 | 0.31 | 有 | 离子/电荷/反离子未体现（季铵、阴离子、叶立德） |
| 188 | chebi-782 | 30 | 2 | 1 | 0.20 | 有 | coverage 门控放宽：残缺局部片段名被当作成功输出 |
| 189 | tiers-111140 | 30 | 2 | 2 | 0.00 | 空 | 母体候选为空（骨架表达失败） |
| 190 | chebi-2204 | 30 | 2 | 2 | 0.39 | 有 | 杂原子中心（P/B）被当作碳命名 |
| 191 | chebi-526 | 30 | 3 | 2 | 0.00 | 空 | 多环母体词干组装失败（稠合/桥/螺未细分） |
| 192 | chebi-3053 | 30 | 3 | 1 | 0.00 | 空 | 桥环/螺环/笼状母体缺 von Baeyer 词干 |
| 193 | chebi-2714 | 30 | 3 | 2 | 0.00 | 空 | 多环母体词干组装失败（稠合/桥/螺未细分） |
| 194 | chebi-64 | 30 | 4 | 2 | 0.00 | 空 | 多环母体词干组装失败（稠合/桥/螺未细分） |
| 195 | tiers-140247 | 30 | 5 | 2 | 0.40 | 有 | coverage 门控放宽：残缺局部片段名被当作成功输出 |
| 196 | chebi-3 | 31 | 3 | 2 | 0.38 | 有 | 离子/电荷/反离子未体现（季铵、阴离子、叶立德） |
| 197 | tiers-131559 | 31 | 4 | 3 | 0.00 | 空 | 多环母体词干组装失败（稠合/桥/螺未细分） |
| 198 | tiers-148964 | 31 | 5 | 3 | 0.00 | 空 | 多环母体词干组装失败（稠合/桥/螺未细分） |
| 199 | chebi-1797 | 32 | 3 | 2 | 0.38 | 有 | 磺酸/硫酸（氢）酯主基团未登记，降级为前缀 |
| 200 | chebi-3208 | 32 | 4 | 2 | 0.00 | 空 | 桥环/螺环/笼状母体缺 von Baeyer 词干 |
| 201 | chebi-968 | 32 | 4 | 2 | 0.31 | 有 | coverage 门控放宽：残缺局部片段名被当作成功输出 |
| 202 | chebi-2251 | 32 | 5 | 5 | 0.00 | 空 | 未注册单环/杂环母体缺失（HW 模板与保留名注册表缺项） |
| 203 | chebi-1787 | 32 | 5 | 2 | 0.17 | 有 | 取代基递归命名残缺（块内原子静默丢弃） |
| 204 | chebi-125 | 32 | 5 | 1 | 0.00 | 空 | 桥环/螺环/笼状母体缺 von Baeyer 词干 |
| 205 | chebi-3261 | 33 | 3 | 3 | 0.38 | 有 | 酰胺/内酰胺互变异构未归一 |
| 206 | tiers-150173 | 33 | 4 | 3 | 0.37 | 有 | 取代基前缀字母序与保留前缀未命中 |
| 207 | chebi-1017 | 33 | 4 | 1 | 0.00 | 空 | 多环母体词干组装失败（稠合/桥/螺未细分） |
| 208 | tiers-96196 | 33 | 4 | 2 | 0.00 | 空 | 稠合环系拆解退化/稠合词干失败 |
| 209 | chebi-3087 | 33 | 5 | 2 | 0.39 | 有 | 取代基递归命名残缺（块内原子静默丢弃） |
| 210 | chebi-568 | 33 | 5 | 1 | 0.00 | 空 | 多环母体词干组装失败（稠合/桥/螺未细分） |
| 211 | chebi-929 | 33 | 5 | 1 | 0.00 | 空 | 稠合环系拆解退化/稠合词干失败 |
| 212 | chebi-1115 | 33 | 6 | 1 | 0.00 | 空 | 桥环/螺环/笼状母体缺 von Baeyer 词干 |
| 213 | chebi-378 | 34 | 1 | 1 | 0.37 | 有 | 离子/电荷/反离子未体现（季铵、阴离子、叶立德） |
| 214 | chebi-2791 | 34 | 2 | 2 | 0.15 | 有 | coverage 门控放宽：残缺局部片段名被当作成功输出 |
| 215 | tiers-131305 | 34 | 3 | 3 | 0.00 | 空 | 链式同种主基团多重性不可表达（二酯/二酰胺/酸酐/偕二醇） |
| 216 | tiers-121896 | 34 | 4 | 2 | 0.39 | 有 | 多环骨架母体名生成错误（伪稠合/拓扑不符/位次缺失） |
| 217 | chebi-1295 | 34 | 4 | 3 | 0.00 | 空 | 多环母体词干组装失败（稠合/桥/螺未细分） |
| 218 | chebi-1243 | 34 | 6 | 1 | 0.00 | 空 | 大环母体/环几何模板缺失（含编号引擎 KeyError） |
| 219 | tiers-67934 | 35 | 0 | 0 | 0.34 | 有 | 离子/电荷/反离子未体现（季铵、阴离子、叶立德） |
| 220 | chebi-3226 | 35 | 2 | 2 | 0.36 | 有 | coverage 门控放宽：残缺局部片段名被当作成功输出 |
| 221 | chebi-17 | 35 | 3 | 1 | 0.29 | 有 | coverage 门控放宽：残缺局部片段名被当作成功输出 |
| 222 | chebi-2184 | 35 | 4 | 2 | 0.19 | 有 | coverage 门控放宽：残缺局部片段名被当作成功输出 |
| 223 | chebi-665 | 35 | 4 | 1 | 0.00 | 空 | 多环母体词干组装失败（稠合/桥/螺未细分） |
| 224 | chebi-1501 | 35 | 4 | 3 | 0.00 | 空 | 稠合环系拆解退化/稠合词干失败 |
| 225 | chebi-1985 | 35 | 4 | 3 | 0.00 | 空 | 多环母体词干组装失败（稠合/桥/螺未细分） |
| 226 | chebi-1905 | 35 | 5 | 2 | 0.00 | 空 | 多环母体词干组装失败（稠合/桥/螺未细分） |
| 227 | chebi-1636 | 35 | 5 | 2 | 0.00 | 空 | 稠合环系拆解退化/稠合词干失败 |
| 228 | tiers-112188 | 35 | 5 | 3 | 0.40 | 有 | coverage 门控放宽：残缺局部片段名被当作成功输出 |
| 229 | chebi-3262 | 35 | 5 | 3 | 0.00 | 空 | 多环母体词干组装失败（稠合/桥/螺未细分） |
| 230 | chebi-2651 | 35 | 6 | 2 | 0.00 | 空 | 多环母体词干组装失败（稠合/桥/螺未细分） |
| 231 | chebi-1186 | 35 | 6 | 2 | 0.00 | 空 | 多环母体词干组装失败（稠合/桥/螺未细分） |
| 232 | chebi-406 | 36 | 2 | 1 | 0.38 | 有 | coverage 门控放宽：残缺局部片段名被当作成功输出 |
| 233 | chebi-3295 | 36 | 4 | 3 | 0.00 | 空 | 多环母体词干组装失败（稠合/桥/螺未细分） |
| 234 | chebi-924 | 36 | 4 | 2 | 0.00 | 空 | 多环母体词干组装失败（稠合/桥/螺未细分） |
| 235 | chebi-1512 | 36 | 4 | 3 | 0.00 | 空 | 稠合环系拆解退化/稠合词干失败 |
| 236 | chebi-3064 | 36 | 5 | 2 | 0.00 | 空 | 桥环/螺环/笼状母体缺 von Baeyer 词干 |
| 237 | chebi-1593 | 36 | 5 | 1 | 0.00 | 空 | 多环母体词干组装失败（稠合/桥/螺未细分） |
| 238 | chebi-788 | 37 | 3 | 2 | 0.00 | 空 | 多环母体词干组装失败（稠合/桥/螺未细分） |
| 239 | chebi-2503 | 38 | 2 | 2 | 0.32 | 有 | 离子/电荷/反离子未体现（季铵、阴离子、叶立德） |
| 240 | chebi-1169 | 38 | 3 | 3 | 0.36 | 有 | 杂原子中心（P/B）被当作碳命名 |
| 241 | chebi-2789 | 39 | 2 | 1 | 0.00 | 空 | 多环母体词干组装失败（稠合/桥/螺未细分） |
| 242 | chebi-1330 | 39 | 2 | 2 | 0.21 | 有 | coverage 门控放宽：残缺局部片段名被当作成功输出 |
| 243 | chebi-990 | 39 | 3 | 3 | 0.38 | 有 | 杂原子中心（P/B）被当作碳命名 |
| 244 | chebi-14 | 39 | 4 | 2 | 0.00 | 空 | 桥环/螺环/笼状母体缺 von Baeyer 词干 |
| 245 | chebi-2386 | 39 | 6 | 2 | 0.00 | 空 | 多环母体词干组装失败（稠合/桥/螺未细分） |
| 246 | chebi-2709 | 39 | 7 | 4 | 0.35 | 有 | coverage 门控放宽：残缺局部片段名被当作成功输出 |
| 247 | chebi-2925 | 40 | 2 | 1 | 0.29 | 有 | 桥环/螺环/笼状母体缺 von Baeyer 词干 |
| 248 | chebi-2687 | 41 | 7 | 4 | 0.31 | 有 | 桥环/螺环/笼状母体缺 von Baeyer 词干 |
| 249 | chebi-2426 | 42 | 5 | 2 | 0.17 | 有 | coverage 门控放宽：残缺局部片段名被当作成功输出 |
| 250 | chebi-2952 | 42 | 5 | 4 | 0.00 | 空 | 多环母体词干组装失败（稠合/桥/螺未细分） |
| 251 | chebi-1135 | 42 | 5 | 4 | 0.00 | 空 | 稠合环系拆解退化/稠合词干失败 |
| 252 | chebi-258 | 42 | 7 | 2 | 0.08 | 有 | coverage 门控放宽：残缺局部片段名被当作成功输出 |
| 253 | chebi-473 | 42 | 8 | 4 | 0.00 | 空 | 多环母体词干组装失败（稠合/桥/螺未细分） |
| 254 | chebi-156 | 43 | 6 | 3 | 0.00 | 空 | 多环母体词干组装失败（稠合/桥/螺未细分） |
| 255 | chebi-1363 | 44 | 0 | 0 | 0.00 | 空 | 链式同种主基团多重性不可表达（二酯/二酰胺/酸酐/偕二醇） |
| 256 | chebi-1978 | 44 | 2 | 2 | 0.00 | 空 | 未注册单环/杂环母体缺失（HW 模板与保留名注册表缺项） |
| 257 | chebi-1896 | 44 | 3 | 2 | 0.00 | 空 | 多环母体词干组装失败（稠合/桥/螺未细分） |
| 258 | chebi-1197 | 44 | 5 | 1 | 0.22 | 有 | coverage 门控放宽：残缺局部片段名被当作成功输出 |
| 259 | chebi-488 | 44 | 6 | 4 | 0.28 | 有 | coverage 门控放宽：残缺局部片段名被当作成功输出 |
| 260 | chebi-2453 | 44 | 6 | 2 | 0.28 | 有 | coverage 门控放宽：残缺局部片段名被当作成功输出 |
| 261 | chebi-1187 | 44 | 7 | 3 | 0.29 | 有 | coverage 门控放宽：残缺局部片段名被当作成功输出 |
| 262 | chebi-520 | 45 | 3 | 2 | 0.18 | 有 | coverage 门控放宽：残缺局部片段名被当作成功输出 |
| 263 | chebi-690 | 45 | 4 | 3 | 0.00 | 空 | 大环母体/环几何模板缺失（含编号引擎 KeyError） |
| 264 | chebi-758 | 45 | 6 | 2 | 0.33 | 有 | 多环骨架母体名生成错误（伪稠合/拓扑不符/位次缺失） |
| 265 | chebi-1466 | 45 | 8 | 1 | 0.00 | 空 | 多环母体词干组装失败（稠合/桥/螺未细分） |
| 266 | chebi-1020 | 46 | 2 | 2 | 0.00 | 空 | 未注册单环/杂环母体缺失（HW 模板与保留名注册表缺项） |
| 267 | chebi-510 | 46 | 4 | 3 | 0.00 | 空 | 未注册单环/杂环母体缺失（HW 模板与保留名注册表缺项） |
| 268 | chebi-3119 | 46 | 5 | 2 | 0.00 | 空 | 稠合环系拆解退化/稠合词干失败 |
| 269 | chebi-332 | 46 | 5 | 2 | 0.10 | 有 | coverage 门控放宽：残缺局部片段名被当作成功输出 |
| 270 | chebi-529 | 46 | 6 | 2 | 0.29 | 有 | 多环骨架母体名生成错误（伪稠合/拓扑不符/位次缺失） |
| 271 | chebi-1661 | 46 | 7 | 3 | 0.34 | 有 | coverage 门控放宽：残缺局部片段名被当作成功输出 |
| 272 | chebi-1488 | 47 | 5 | 2 | 0.00 | 空 | 稠合环系拆解退化/稠合词干失败 |
| 273 | chebi-1916 | 47 | 6 | 4 | 0.39 | 有 | 离子/电荷/反离子未体现（季铵、阴离子、叶立德） |
| 274 | chebi-1146 | 48 | 4 | 3 | 0.00 | 空 | 未注册单环/杂环母体缺失（HW 模板与保留名注册表缺项） |
| 275 | chebi-435 | 48 | 4 | 3 | 0.00 | 空 | 未注册单环/杂环母体缺失（HW 模板与保留名注册表缺项） |
| 276 | chebi-2333 | 48 | 5 | 1 | 0.22 | 有 | coverage 门控放宽：残缺局部片段名被当作成功输出 |
| 277 | chebi-390 | 48 | 7 | 2 | 0.26 | 有 | 多环骨架母体名生成错误（伪稠合/拓扑不符/位次缺失） |
| 278 | chebi-2176 | 50 | 4 | 3 | 0.30 | 有 | coverage 门控放宽：残缺局部片段名被当作成功输出 |
| 279 | chebi-2963 | 50 | 5 | 1 | 0.03 | 有 | coverage 门控放宽：残缺局部片段名被当作成功输出 |
| 280 | chebi-2488 | 51 | 3 | 3 | 0.00 | 空 | 未注册单环/杂环母体缺失（HW 模板与保留名注册表缺项） |
| 281 | chebi-3130 | 51 | 4 | 1 | 0.15 | 有 | coverage 门控放宽：残缺局部片段名被当作成功输出 |
| 282 | chebi-1461 | 51 | 6 | 5 | 0.31 | 有 | coverage 门控放宽：残缺局部片段名被当作成功输出 |
| 283 | chebi-151 | 52 | 2 | 2 | 0.19 | 有 | 取代基递归命名残缺（块内原子静默丢弃） |
| 284 | chebi-1082 | 52 | 4 | 3 | 0.35 | 有 | coverage 门控放宽：残缺局部片段名被当作成功输出 |
| 285 | chebi-2826 | 52 | 6 | 4 | 0.00 | 空 | 母体选择与主基优先级错误 |
| 286 | chebi-249 | 53 | 4 | 1 | 0.14 | 有 | coverage 门控放宽：残缺局部片段名被当作成功输出 |
| 287 | chebi-1729 | 54 | 3 | 2 | 0.32 | 有 | 杂原子中心（P/B）被当作碳命名 |
| 288 | chebi-2710 | 54 | 6 | 5 | 0.22 | 有 | 取代基递归命名残缺（块内原子静默丢弃） |
| 289 | chebi-203 | 55 | 3 | 2 | 0.35 | 有 | 母体选择与主基优先级错误 |
| 290 | chebi-3120 | 55 | 3 | 2 | 0.35 | 有 | 杂原子中心（P/B）被当作碳命名 |
| 291 | chebi-2867 | 55 | 6 | 5 | 0.30 | 有 | 离子/电荷/反离子未体现（季铵、阴离子、叶立德） |
| 292 | chebi-1238 | 55 | 9 | 2 | 0.21 | 有 | 稠合环系拆解退化/稠合词干失败 |
| 293 | chebi-2839 | 56 | 3 | 1 | 0.00 | 空 | 多环母体词干组装失败（稠合/桥/螺未细分） |
| 294 | chebi-2534 | 56 | 4 | 3 | 0.35 | 有 | 杂原子中心（P/B）被当作碳命名 |
| 295 | chebi-2536 | 56 | 5 | 5 | 0.34 | 有 | 寡糖/糖苷母体单糖选取与位次、立体描述错误 |
| 296 | chebi-2747 | 56 | 5 | 5 | 0.20 | 有 | 寡糖/糖苷母体单糖选取与位次、立体描述错误 |
| 297 | chebi-1358 | 56 | 5 | 2 | 0.37 | 有 | 离子/电荷/反离子未体现（季铵、阴离子、叶立德） |
| 298 | chebi-2923 | 57 | 3 | 1 | 0.00 | 空 | 多环母体词干组装失败（稠合/桥/螺未细分） |
| 299 | chebi-1394 | 57 | 3 | 2 | 0.40 | 有 | 杂原子中心（P/B）被当作碳命名 |
| 300 | chebi-225 | 57 | 4 | 4 | 0.20 | 有 | 环上取代基位次编号方向错误 |
| 301 | chebi-1599 | 57 | 6 | 1 | 0.03 | 有 | coverage 门控放宽：残缺局部片段名被当作成功输出 |
| 302 | chebi-2702 | 58 | 4 | 4 | 0.31 | 有 | coverage 门控放宽：残缺局部片段名被当作成功输出 |
| 303 | chebi-2611 | 58 | 7 | 3 | 0.38 | 有 | 母体选择与主基优先级错误 |
| 304 | chebi-2665 | 59 | 6 | 2 | 0.03 | 有 | coverage 门控放宽：残缺局部片段名被当作成功输出 |
| 305 | chebi-966 | 59 | 12 | 3 | 0.09 | 有 | 多环骨架母体名生成错误（伪稠合/拓扑不符/位次缺失） |
| 306 | chebi-1857 | 60 | 5 | 1 | 0.20 | 有 | coverage 门控放宽：残缺局部片段名被当作成功输出 |
| 307 | chebi-3147 | 60 | 9 | 2 | 0.19 | 有 | coverage 门控放宽：残缺局部片段名被当作成功输出 |
| 308 | chebi-2951 | 61 | 3 | 2 | 0.36 | 有 | 杂原子中心（P/B）被当作碳命名 |
| 309 | chebi-902 | 61 | 6 | 6 | 0.31 | 有 | 杂原子中心（P/B）被当作碳命名 |
| 310 | chebi-3275 | 62 | 5 | 1 | 0.10 | 有 | coverage 门控放宽：残缺局部片段名被当作成功输出 |
| 311 | chebi-1066 | 63 | 5 | 1 | 0.06 | 有 | coverage 门控放宽：残缺局部片段名被当作成功输出 |
| 312 | chebi-210 | 63 | 10 | 5 | 0.00 | 空 | 桥环/螺环/笼状母体缺 von Baeyer 词干 |
| 313 | chebi-886 | 64 | 1 | 1 | 0.40 | 有 | 酰胺/内酰胺互变异构未归一 |
| 314 | chebi-1227 | 64 | 4 | 4 | 0.01 | 有 | coverage 门控放宽：残缺局部片段名被当作成功输出 |
| 315 | chebi-31 | 65 | 3 | 2 | 0.00 | 空 | 大环母体/环几何模板缺失（含编号引擎 KeyError） |
| 316 | chebi-2171 | 66 | 7 | 4 | 0.34 | 有 | 稠环氢化程度/指示氢表达错误 |
| 317 | chebi-948 | 67 | 8 | 5 | 0.00 | 空 | 多环母体词干组装失败（稠合/桥/螺未细分） |
| 318 | chebi-599 | 68 | 8 | 8 | 0.30 | 有 | coverage 门控放宽：残缺局部片段名被当作成功输出 |
| 319 | chebi-439 | 72 | 5 | 5 | 0.00 | 空 | 大环母体/环几何模板缺失（含编号引擎 KeyError） |
| 320 | chebi-1034 | 72 | 10 | 5 | 0.28 | 有 | coverage 门控放宽：残缺局部片段名被当作成功输出 |
| 321 | chebi-2580 | 72 | 10 | 5 | 0.28 | 有 | 取代基递归命名残缺（块内原子静默丢弃） |
| 322 | chebi-1840 | 73 | 6 | 6 | 0.35 | 有 | 寡糖/糖苷母体单糖选取与位次、立体描述错误 |
| 323 | chebi-3149 | 73 | 9 | 5 | 0.36 | 有 | 多环骨架母体名生成错误（伪稠合/拓扑不符/位次缺失） |
| 324 | chebi-47 | 74 | 6 | 4 | 0.00 | 空 | 链式同种主基团多重性不可表达（二酯/二酰胺/酸酐/偕二醇） |
| 325 | chebi-2054 | 74 | 7 | 3 | 0.18 | 有 | 杂原子中心（P/B）被当作碳命名 |
| 326 | chebi-2304 | 74 | 7 | 3 | 0.17 | 有 | 杂原子中心（P/B）被当作碳命名 |
| 327 | chebi-68 | 75 | 4 | 2 | 0.19 | 有 | coverage 门控放宽：残缺局部片段名被当作成功输出 |
| 328 | chebi-3029 | 76 | 4 | 2 | 0.04 | 有 | coverage 门控放宽：残缺局部片段名被当作成功输出 |
| 329 | chebi-3083 | 76 | 6 | 4 | 0.03 | 有 | coverage 门控放宽：残缺局部片段名被当作成功输出 |
| 330 | chebi-2398 | 76 | 6 | 4 | 0.00 | 空 | 链式同种主基团多重性不可表达（二酯/二酰胺/酸酐/偕二醇） |
| 331 | chebi-1280 | 76 | 7 | 3 | 0.15 | 有 | 杂原子中心（P/B）被当作碳命名 |
| 332 | chebi-3098 | 77 | 7 | 3 | 0.16 | 有 | 杂原子中心（P/B）被当作碳命名 |
| 333 | chebi-1676 | 77 | 7 | 3 | 0.17 | 有 | 杂原子中心（P/B）被当作碳命名 |
| 334 | chebi-2586 | 77 | 8 | 6 | 0.07 | 有 | coverage 门控放宽：残缺局部片段名被当作成功输出 |
| 335 | chebi-2321 | 79 | 6 | 4 | 0.14 | 有 | coverage 门控放宽：残缺局部片段名被当作成功输出 |
| 336 | chebi-1077 | 85 | 2 | 2 | 0.14 | 有 | coverage 门控放宽：残缺局部片段名被当作成功输出 |
| 337 | chebi-1548 | 85 | 9 | 3 | 0.21 | 有 | 取代基递归命名残缺（块内原子静默丢弃） |
| 338 | chebi-2382 | 88 | 3 | 3 | 0.00 | 空 | 链式同种主基团多重性不可表达（二酯/二酰胺/酸酐/偕二醇） |
| 339 | chebi-3287 | 89 | 2 | 2 | 0.12 | 有 | 杂原子中心（P/B）被当作碳命名 |
| 340 | chebi-2624 | 90 | 3 | 3 | 0.35 | 有 | 杂原子中心（P/B）被当作碳命名 |
| 341 | chebi-2954 | 90 | 13 | 1 | 0.00 | 空 | 多环母体词干组装失败（稠合/桥/螺未细分） |
| 342 | chebi-1665 | 90 | 13 | 1 | 0.00 | 空 | 多环母体词干组装失败（稠合/桥/螺未细分） |
| 343 | chebi-982 | 91 | 5 | 5 | 0.19 | 有 | 离子/电荷/反离子未体现（季铵、阴离子、叶立德） |
| 344 | chebi-1918 | 92 | 6 | 5 | 0.03 | 有 | coverage 门控放宽：残缺局部片段名被当作成功输出 |
| 345 | chebi-2492 | 95 | 8 | 2 | 0.20 | 有 | coverage 门控放宽：残缺局部片段名被当作成功输出 |
| 346 | chebi-2988 | 101 | 8 | 8 | 0.26 | 有 | 寡糖/糖苷母体单糖选取与位次、立体描述错误 |
| 347 | chebi-1830 | 101 | 10 | 10 | 0.17 | 有 | 杂原子中心（P/B）被当作碳命名 |
| 348 | chebi-183 | 108 | 3 | 3 | 0.36 | 有 | 杂原子中心（P/B）被当作碳命名 |
| 349 | chebi-3221 | 108 | 4 | 2 | 0.09 | 有 | coverage 门控放宽：残缺局部片段名被当作成功输出 |
| 350 | chebi-2013 | 108 | 9 | 8 | 0.00 | 空 | 多环母体词干组装失败（稠合/桥/螺未细分） |
| 351 | chebi-1854 | 109 | 8 | 6 | 0.00 | 空 | 桥环/螺环/笼状母体缺 von Baeyer 词干 |
| 352 | chebi-766 | 110 | 2 | 2 | 0.36 | 有 | 杂原子中心（P/B）被当作碳命名 |
| 353 | chebi-1018 | 111 | 4 | 4 | 0.32 | 有 | 酰胺/内酰胺互变异构未归一 |
| 354 | chebi-2149 | 111 | 13 | 4 | 0.33 | 有 | coverage 门控放宽：残缺局部片段名被当作成功输出 |
| 355 | chebi-2062 | 113 | 10 | 8 | 0.00 | 空 | 稠合环系拆解退化/稠合词干失败 |
| 356 | chebi-2396 | 121 | 6 | 6 | 0.10 | 有 | 环上取代基位次编号方向错误 |
| 357 | chebi-3004 | 125 | 6 | 5 | 0.00 | 空 | 多环母体词干组装失败（稠合/桥/螺未细分） |
| 358 | chebi-649 | 130 | 12 | 12 | 0.37 | 有 | 寡糖/糖苷母体单糖选取与位次、立体描述错误 |
| 359 | chebi-559 | 132 | 4 | 4 | 0.19 | 有 | 杂原子中心（P/B）被当作碳命名 |
| 360 | chebi-1763 | 132 | 16 | 4 | 0.00 | 空 | 桥环/螺环/笼状母体缺 von Baeyer 词干 |
| 361 | chebi-3272 | 132 | 19 | 7 | 0.00 | 空 | 大环母体/环几何模板缺失（含编号引擎 KeyError） |
| 362 | chebi-1144 | 137 | 3 | 3 | 0.25 | 有 | 酰胺/内酰胺互变异构未归一 |
| 363 | chebi-2411 | 137 | 3 | 3 | 0.25 | 有 | 杂原子中心（P/B）被当作碳命名 |
| 364 | chebi-3244 | 141 | 18 | 14 | 0.35 | 有 | 杂原子中心（P/B）被当作碳命名 |
| 365 | chebi-2813 | 144 | 13 | 13 | 0.30 | 有 | 寡糖/糖苷母体单糖选取与位次、立体描述错误 |
| 366 | chebi-1738 | 150 | 8 | 6 | 0.00 | 空 | 母体候选为空（骨架表达失败） |
| 367 | chebi-650 | 154 | 12 | 12 | 0.39 | 有 | 寡糖/糖苷母体单糖选取与位次、立体描述错误 |
| 368 | chebi-1645 | 161 | 11 | 11 | 0.10 | 有 | 杂原子中心（P/B）被当作碳命名 |
| 369 | chebi-1772 | 167 | 8 | 8 | 0.25 | 有 | 酰胺/内酰胺互变异构未归一 |
| 370 | chebi-2877 | 171 | 6 | 6 | 0.17 | 有 | 酰胺/内酰胺互变异构未归一 |
| 371 | chebi-388 | 236 | 17 | 17 | 0.35 | 有 | 寡糖/糖苷母体单糖选取与位次、立体描述错误 |
| 372 | chebi-352 | 261 | 12 | 12 | 0.20 | 有 | 杂原子中心（P/B）被当作碳命名 |
| 373 | chebi-1429 | 264 | 5 | 4 | 0.18 | 有 | 酰胺/内酰胺互变异构未归一 |
| 374 | chebi-1544 | 278 | 12 | 12 | 0.36 | 有 | 杂原子中心（P/B）被当作碳命名 |
| 375 | chebi-2969 | 383 | 45 | 38 | 0.00 | 空 | 寡核苷酸链无命名支持 |
