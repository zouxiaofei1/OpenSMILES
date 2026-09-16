# SMILES → IUPAC 双语命名错误分析报告

- **数据**：`benchmarks/merged_benchmark.json`，3976 条（chebi 3188 + smiles_tiers 788；后者含中文金标 748 条）
- **被测**：`src/namepredict/` 当前 HEAD（工作区 clean）
- **判分**：`benchmarks/benchmark.py`，`normalize_en` 后字符串全等 —— 括号、连字符、位次写法都是硬判分项
- **方法**：全量跑一遍导出全部错例 → 按复杂度（重原子数 → 环数 → 碳数）升序切 6 块 → 6 个并行 subagent 逐块用 RDKit 核实结构 + 比对 IUPAC 规则（`docs/iupac/`）归纳根因 → 本报告按根因跨块合并
- **脚本与原始报告**：`tmp/`（git 忽略）

---

## 1. 评测基线

| 指标 | 值 |
|---|---|
| EN 准确率 | **61.95%**（2463 / 3976） |
| ZH 准确率 | **77.41%**（579 / 748） |
| dual 准确率 | 61.62%（2450 / 3976） |
| 错例总数 | **1526**（占 38.4%） |

错例构成：仅 EN 错 1357 条、仅 ZH 错 13 条、EN/ZH 双错 156 条。chebi 源全部不评中文（无中文金标），所以中文侧只能从 smiles_tiers 的 13 条单错里观察。

按复杂度分块后的错例分布（复杂度越高越难）：

| 分块 | HAC 范围 | 条数 |
|---|---:|---:|
| chunk_01 | 3 – 17 | 255 |
| chunk_02 | 17 – 24 | 255 |
| chunk_03 | 24 – 30 | 255 |
| chunk_04 | 30 – 36 | 255 |
| chunk_05 | 36 – 56 | 255 |
| chunk_06 | 56 – 383 | 251 |

---

## 2. 全局错误主题

把 6 块报告的 100 余个局部模式按**根因**合并后，1526 条错例收敛到 13 个可修主题 + 2 个范围外主题。下表按错例数排序：

| 主题 | 错例数 | 主要根因模块 | 修复难度 |
|---|---:|---|---|
| **G3 含氧酸（磷/硫）类名后缀与母体优先级** | 219 | `layer2/parent_select.py`、`layer1/fg_registry.py`、`layer5/assembler.py` | 中 |
| **G4 稠环/多环母体、取向与保留名** | 166 | `layer2/fused_system.py`、`layer4/fused_orientation.py`、`constants.RETAINED_FUSION_ALIASES` | 高 |
| **G9 覆盖率：片段丢失与空输出** | 162 | `layer3/coverage.py`、`layer2/ring_expression_policy.py`、`layer0/salt.py` | 中高 |
| **G7 括号、连字符与词界** | 147 | `layer5/assembler_prefixes.py`、`layer5/assembler.py` | 低 |
| **G1 碳酰/含氮衍生物保留名母体缺失** | 132 | `layer1/fg_registry.py`（FG_SPECS 缺 10+ 类）、`constants.CHAIN_RETAINED` | 中低 |
| **G8 环编号方向与位次 tie-break** | 84 | `layer4/numbering.py`、`layer4/locant_calc.py`、`layer4/fused_numbering.py` | 中 |
| **G5 指示氢与氢化程度** | 81 | `layer4/indicated_hydrogen.py`、`layer2/ring_scaffold.py` | 中 |
| **G6 取代基引用顺序（P-14.5 字母序）** | 72 | `layer5/assembler_prefixes.py::_build_prefix`、`tools/re.py::alkyl_alpha_key` | 低 |
| **G10 盐、抗衡离子与电荷态** | 68 | `layer0/salt.py`、`layer0/charge.py` | 中低 |
| **G2 保留前缀与官能团表达缺失** | 43 | `tools/anchored_table.py`（diazenyl/disulfanyl/carbonimidoyl/…） | 低 |
| **G12 主链与母体骨架选择** | 42 | `layer2/parent_select.py`、`layer2/parent_skeleton.py` | 高 |
| **G11 立体描述符** | 29 | `layer5/stereo.py`、`layer4/locant_calc.py` | 中 |
| **G13 中文侧（EN 通过、ZH 失败）** | 12 | `layer5/assembler*.py` 中文侧词表与括号策略 | 低 |
| 范围外 A：桥环 / 螺环 | ~208* | 需 von Baeyer / spiro 表达通道 | 高 |
| 范围外 B：铵 / 鎓命名 | ~133* | 需鎓母体与电荷归属 | 高 |

\* 范围外计数为 6 块报告的上报值之和；这些分子中有一部分的主因已被计入上表模式，故与 13 个主题相加会超过 1526（13 个主题合计 1255 条，剩余 271 条的主因即落在两个范围外主题）。

**读数提示**：G3、G7、G9、G1 四个主题合计 660 条，占全部错例的 43%，且根因都相对集中（单一排序表 / 单一括号策略 / 单一覆盖校验 / 单一官能团注册表），是性价比最高的四块。

---

## 3. 主题详解

### G3 含氧酸（磷/硫）类名后缀与母体优先级 — 219 条

**症状**：分子同时含磷酸酯（或硫酸酯、膦酸酯、硫代酸酯、多聚磷酸链）与羧酸酯/核苷/糖/甘油骨架时，PRED 一律以**碳骨架**为母体、把含氧酸降级成 `…phosphoryl]oxy` / `…phosphonooxy` 取代前缀；GOLD 反过来，把含氧酸整体作**功能类别后缀**（`… phosphate` / `hydrogen phosphate` / `… sulfate` / `… thioate`）。甘油磷脂类还连带甘油骨架位次 1↔3 反转。

**根因**：
- `layer1/fg_registry.py` 的 p41 优先等级里 `oxoacid`(9) 与 `ester`(9)、`amide`(11) 相对序未体现「羧酸类似酸（-P(O)(OH)₂ / -SO₃H）排在酯与酰胺之前」（P-41/P-42）；`layer2/principal.py::select_principal_group` 直接取最小 priority。
- `layer1/analyzer.py::_OXO_Z_ANCHORED = {phosphate, phosphonate, sulfate}` 只在「P/S 两端都是 P–O–C」时让杂原子自任母体；P 只是二酯一臂、碳骨架另带酯基时该分支不生效。
- `layer5/chain_engine.py` 的 `("phosphate", 2) → dihydrogen phosphate` 词尾钩子在二酯上取错 n；`layer5/assembler.py::_phosphoryl_sub_names` 把多聚磷酸拆成「由内向外」的嵌套前缀，与 GOLD「外层一段做后缀」的切分点相反。
- 硫酯同理：`P_5_Preferred_IUPAC_Names.md` 明确给出 `S-ethyl hexanethioate`、`lithium 4-(ethylsulfanyl)-4-oxobutanethioate` 这类 PIN，PRED 只有取代式出口，写不出 `-thioate` 后缀。

**IUPAC 依据**：P-44.1.2（优先原子序 N > P > … > O > S > C）；P-41/P-43（后缀优先序：羧酸及其类似酸 > 酯 > 酰胺 > …）；P-67（非碳酸含氧酸的酯，`… phosphate` / `… hydrogen phosphate` / `… sulfate`）；P-66.5（硫代酸的 S-酯）。

**全部 id**（219 条）：

chebi-1396 chebi-1705 tiers-107394 tiers-43418 tiers-68805 chebi-2256 tiers-51373 chebi-346
tiers-8543 chebi-3018 chebi-800 chebi-121 chebi-270 tiers-130171 tiers-76395 tiers-10747
chebi-907 chebi-1365 chebi-1812 chebi-2527 chebi-1179 chebi-810 chebi-1427 tiers-155306
chebi-3200 chebi-3131 chebi-1744 chebi-1684 chebi-2100 chebi-2511 chebi-2164 chebi-393
chebi-2957 chebi-571 chebi-227 chebi-333 chebi-369 chebi-1731 chebi-2070 chebi-2721
chebi-2493 chebi-3229 chebi-2025 chebi-2755 chebi-1371 chebi-2204 chebi-925 chebi-2918
chebi-546 chebi-607 chebi-855 chebi-1999 chebi-1028 chebi-1204 chebi-279 chebi-430
chebi-3 chebi-2591 chebi-1401 chebi-295 chebi-1983 chebi-1843 chebi-1860 chebi-1430
chebi-2668 chebi-645 chebi-1058 chebi-2081 chebi-206 chebi-2535 chebi-2606 chebi-1876
chebi-189 chebi-2464 chebi-2898 chebi-2185 chebi-3236 chebi-1997 chebi-2363 chebi-1773
chebi-274 chebi-3182 chebi-1677 chebi-3082 chebi-1169 chebi-2016 chebi-990 chebi-941
chebi-2881 chebi-629 chebi-71 chebi-679 chebi-251 chebi-955 chebi-2289 chebi-2934
chebi-2582 chebi-2680 chebi-3165 chebi-1475 chebi-179 chebi-195 chebi-1083 chebi-1877
chebi-3101 chebi-677 chebi-2389 chebi-1325 chebi-3276 chebi-2350 chebi-3196 chebi-3071
chebi-1075 chebi-2073 chebi-1184 chebi-2501 chebi-1718 chebi-238 chebi-3015 chebi-324
chebi-2630 chebi-2783 chebi-3052 chebi-122 chebi-2968 chebi-2608 chebi-416 chebi-1527
chebi-1133 chebi-2352 chebi-1651 chebi-3039 chebi-2150 chebi-2776 chebi-1324 chebi-1729
chebi-1176 chebi-3238 chebi-1719 chebi-1656 chebi-1950 chebi-203 chebi-3120 chebi-1790
chebi-453 chebi-168 chebi-2726 chebi-2534 chebi-252 chebi-1394 chebi-567 chebi-1062
chebi-2172 chebi-357 chebi-2958 chebi-2012 chebi-2584 chebi-3220 chebi-1045 chebi-1620
chebi-2992 chebi-719 chebi-1936 chebi-331 chebi-2951 chebi-1923 chebi-2056 chebi-2733
chebi-2487 chebi-2279 chebi-1913 chebi-504 chebi-28 chebi-2281 chebi-2836 chebi-887
chebi-2756 chebi-3110 chebi-1630 chebi-128 chebi-1879 chebi-3084 chebi-2844 chebi-3286
chebi-237 chebi-1111 chebi-2054 chebi-2304 chebi-1313 chebi-363 chebi-223 chebi-1280
chebi-1676 chebi-3098 chebi-2383 chebi-2695 chebi-3168 chebi-1044 chebi-904 chebi-2633
chebi-803 chebi-1602 chebi-1390 chebi-3090 chebi-3287 chebi-2624 chebi-762 chebi-1262
chebi-2947 chebi-1637 chebi-2557 chebi-183 chebi-766 chebi-1018 chebi-1144 chebi-2411
chebi-3244 chebi-1645 chebi-2969

---

### G4 稠环/多环母体、取向与保留名 — 166 条

**症状**：同一骨架，GOLD 与 PRED 选了**不同的母体环系或不同的环系名**。四种表现：
1. **缺保留名**：GOLD 用 `picene` / `tetracene` / `porphyrin` / `indolizine` / `1,8-naphthyridine` / `benzo[a]anthracene` / `pyrrolizine` / `benzo[b]thiophene`，PRED 拼出 `cyclohexa[a]chrysene` / `pyrrolo[1,2-a]pyridine` / `pyrido[3,2-b]pyridine` 等并合名。
2. **取向/位标错**：`indolo[4,3-fg]` ↔ `indolo[4,3a-f]`、`naphtho[3,2-h]` ↔ `naphtho[2,3-h]`、`pyrazolo[5,4-d]` ↔ `pyrazolo[3,4-d]`；菲类的固定交点在 GOLD 是 `10a`，PRED 用 `4b`（`layer2/ring_scaffold.py::PHENANTHRENE_LABELS` 没有 `10a` 位），导致取代基位次整体镜像。
3. **缺 von Baeyer 出口**：GOLD 用 `bicyclo[4.2.0]oct-2-ene`（头孢母核）、`7-oxabicyclo[4.1.0]heptene`、`tetracyclo[9.7.0.0²,⁷.0¹²,¹⁶]`，PRED 走稠合名（`azetidino[2,1-b]1,3-thiazine`、`oxepane`）。全库 GOLD 有 47 处 `cyclo[` 骨架而 PRED **一处都没有**——说明 L2 完全没有 von Baeyer 输出通道。β-内酰胺类（8 条）应优先判 von Baeyer 适用性（P-23 优先于 P-25）。
4. **附带机械 bug**：并环前缀拼接时位次串为空，输出 `benzodioxolo[-f]benzodioxol` 这类带空位次的名字（约 15 条，可单独低难度修复）。

**根因**：`constants.RETAINED_FUSION_ALIASES` 白名单过短；`layer2/fused_system.py` 的母体选择未做「等价并合名按环数/最低位次集挑优选名」；`layer2/ring_expression_policy.py` 一律选稠合优先；`layer4/fused_numbering.py` + `fused_orientation.py` 的固定交点与起始边规则不全。

**IUPAC 依据**：P-25.1（保留名优先序）、P-25.2/P-25.3（方位与位标）、P-23（extended von Baeyer，饱和多脂环与桥环）、P-14.4(a)（菲/蒽类固定编号）。

**全部 id**（166 条）：

tiers-6167 tiers-44245 tiers-916 chebi-1919 chebi-155 tiers-52889 tiers-17996 tiers-96844
chebi-1888 chebi-2893 chebi-1820 chebi-167 chebi-2550 chebi-2210 chebi-1098 chebi-1588
chebi-1570 chebi-1802 chebi-1727 tiers-19300 chebi-219 chebi-1386 chebi-3209 chebi-2576
chebi-1639 chebi-1737 chebi-1145 chebi-2173 chebi-1270 chebi-1312 chebi-952 chebi-1894
chebi-2207 chebi-2447 chebi-397 chebi-3153 chebi-1106 chebi-2311 chebi-832 chebi-2148
chebi-465 chebi-626 chebi-604 chebi-289 chebi-241 chebi-188 chebi-1335 chebi-281
chebi-1873 tiers-121680 chebi-667 chebi-1086 chebi-2596 chebi-162 chebi-139 tiers-147852
chebi-1485 chebi-2685 tiers-102877 tiers-133227 chebi-2128 chebi-84 chebi-1824 chebi-3181
tiers-156243 chebi-1087 chebi-2714 chebi-943 tiers-92000 chebi-1775 chebi-3140 tiers-140247
chebi-3193 chebi-1481 chebi-3134 chebi-294 chebi-3225 chebi-2723 chebi-2623 chebi-351
chebi-3122 chebi-632 tiers-96196 chebi-51 chebi-986 chebi-1276 chebi-1990 chebi-2745
chebi-1576 chebi-2953 chebi-1368 chebi-1295 tiers-121896 chebi-2984 chebi-3164 chebi-1026
chebi-1042 chebi-1906 chebi-1985 chebi-1501 chebi-2964 chebi-3262 chebi-2084 chebi-1580
chebi-366 chebi-376 chebi-3112 chebi-406 chebi-3295 chebi-1833 chebi-2126 chebi-788
chebi-2842 chebi-1398 chebi-1121 chebi-840 chebi-2808 chebi-2270 chebi-2495 chebi-1135
chebi-2952 chebi-2237 chebi-2111 chebi-53 chebi-1162 chebi-2863 chebi-758 chebi-1351
chebi-839 chebi-3251 chebi-529 chebi-1916 chebi-2142 chebi-1385 chebi-805 chebi-2425
chebi-1959 chebi-2867 chebi-1358 chebi-2219 chebi-2670 chebi-286 chebi-966 chebi-3147
chebi-541 chebi-1110 chebi-292 chebi-673 chebi-2446 chebi-1816 chebi-2793 chebi-2349
chebi-1340 chebi-2676 chebi-1034 chebi-2580 chebi-3149 chebi-47 chebi-2398 chebi-1188
chebi-3174 chebi-3049 chebi-1292 chebi-2492 chebi-1107 chebi-2149

---

### G9 覆盖率：片段丢失与空输出 — 162 条

**症状**：PRED 只输出分子的一小块（`benzoate`、`propanoic acid`、`3-methylbutanoate`、`methyl dihydrogen phosphate`、`acetate`），或整条为空字符串。空输出共约 60 条，其余是「只命名了一个片段」。

**根因**：
- `layer3/coverage.py::build_coverage_ledger` 的 `complete()` **未被上层强制校验**：未被 principal/取代基认领的原子既不报错也不进名字，`layer5/assembler.py` 照常输出残名（`layer3/substituent_extractor.py` 注释即写明「未命名成功的 claim 静默跳过」）。
- 空输出多因 `layer2/ring_expression_policy.py::_unsupported_typed_ring` 把候选全部拦掉 → 候选数 0；或 `layer2/parent_select.py::_collect_candidates` 对桥/螺/大环内酯/环肽返回空且**不降级到次优候选**。
- `layer0/salt.py::_from_frags` 要求「有机片段恰好 1 个 + 单价碱金属」，Ca²⁺、多有机组分、多价抗衡离子一律返回 None，返回整分子后又只挑到一个片段。

**IUPAC 依据**：无（工程性问题：一个名字必须覆盖分子全部原子，P-14.4 的「名称唯一对应结构」）。这是**唯一一个「加断言就能拿到一部分准确率」**的主题——`complete()` 不通过就拒绝输出，至少不会输出错误答案，但也可能把「残名恰好等于 GOLD」的少数巧合算丢，需实测。

**全部 id**（162 条）：

chebi-2948 chebi-2752 tiers-5188 chebi-1012 chebi-3285 chebi-2815 chebi-3109 chebi-318
chebi-1734 tiers-13388 chebi-3267 tiers-59592 tiers-85289 chebi-896 chebi-1147 chebi-432
chebi-2804 tiers-7355 tiers-80645 tiers-18169 chebi-298 chebi-846 chebi-474 tiers-124039
chebi-22 chebi-2478 tiers-49651 chebi-2517 tiers-83257 chebi-1459 tiers-30123 chebi-2034
chebi-1979 chebi-2136 chebi-1864 chebi-2616 chebi-594 tiers-96746 chebi-1560 chebi-503
chebi-263 chebi-2427 chebi-342 chebi-1356 tiers-146758 chebi-2143 chebi-2409 chebi-1196
chebi-606 chebi-1567 chebi-2589 chebi-2242 chebi-2938 chebi-1248 tiers-104849 chebi-2344
chebi-253 chebi-2351 chebi-661 chebi-916 tiers-108409 tiers-98594 tiers-98600 chebi-782
chebi-3053 chebi-526 chebi-64 chebi-2151 chebi-477 tiers-148964 chebi-968 chebi-3208
tiers-116774 chebi-421 chebi-1017 chebi-2791 tiers-131305 chebi-2831 chebi-1243 tiers-67934
chebi-3226 chebi-17 chebi-665 chebi-2184 chebi-1905 chebi-1636 chebi-2651 chebi-804
chebi-3295 chebi-3064 tiers-121943 chebi-2789 chebi-2944 chebi-2386 chebi-2687 chebi-258
chebi-473 chebi-1687 chebi-1363 chebi-1978 chebi-1896 chebi-1197 chebi-2453 chebi-488
chebi-1187 chebi-520 chebi-690 chebi-1466 chebi-1020 chebi-510 chebi-332 chebi-3119
chebi-104 chebi-1661 chebi-1488 chebi-1146 chebi-2774 chebi-435 chebi-2333 chebi-2963
chebi-1747 chebi-2488 chebi-151 chebi-1082 chebi-2826 chebi-2839 chebi-2923 chebi-1599
chebi-2702 chebi-2665 chebi-1857 chebi-2127 chebi-3275 chebi-1066 chebi-210 chebi-1227
chebi-31 chebi-948 chebi-599 chebi-439 chebi-2241 chebi-68 chebi-3029 chebi-3083
chebi-2586 chebi-2321 chebi-2512 chebi-1077 chebi-1548 chebi-2362 chebi-2382 chebi-1665
chebi-2954 chebi-1918 chebi-3221 chebi-2013 chebi-1854 chebi-2062 chebi-3004 chebi-1763
chebi-3272 chebi-1738

---

### G7 括号、连字符与词界 — 147 条

**症状**：把 GOLD 与 PRED 的**所有括号、连字符、空格删掉后两串完全相同**——原子、位次、词干全对，只差括弧落在哪个词界。子型：
- `(9Z)-octadec-9-enoyloxy` ↔ `(9Z)-octadec-9-enoyl]oxy`（酰氧基的 `-oxy` 被拆出括号外）
- `(…-yl)amino` ↔ `(…-ylamino)`；`[…-yl]sulfanyl` ↔ `[…-ylsulfanyl]`
- 嵌套围栏少一层：`6-[[(3S,…)yl]oxy]-` ↔ `6-[(3S,…)yl]oxy-`
- 首位取代基免括号规则未落地：`methylimino(oxo)methane` ↔ `(methylimino)-oxomethane`
- 数字组方括号丢失：`naphtho[2,1-f][1]benzofuran` ↔ `naphtho[2,1-f]benzofuran`

**根因**：`layer5/assembler_prefixes.py::_front_needs_enclosure()` 的 `if "(" in base:` 分支会把含立体描述符的酰基（`(9Z)-octadec-9-enoyl`）也判为「需要外层围栏」，于是 `-oxy` 被拆出括号；上方 `and not base.endswith("oyl")` 只挡住了 R/S 开头的情形，挡不住 Z/E。`_split_bridge_suffix()` / `_bridge_body()` 的 merge 只对 `sulfanyl/sulfinyl` 生效。整体上缺少一条统一的「前缀是否以位次/数字开头 → 必须括起」判据，也缺少括弧最小化。

**IUPAC 依据**：P-1 §14（Enclosing marks：独立片段的围栏要跳一级、同名围栏不得连写）、P-16.5.1.1~1.4（bis/tris 与嵌套）、P-16.5.1.3.1（首位取代基例外）、P-14.2（倍增前缀）。

**注意（重要）**：本主题中**相当一部分是 GOLD 自身的体例摇摆**，不是 PRED 的错。同一片段 `(9Z)-octadec-9-enoyl` + `oxy`，GOLD 在 chebi-1722 写 `enoyloxy`、在 chebi-122 写 `enoyl]oxy`，同一 chunk 内 `oyloxy` 43 处 vs `oyl]oxy` 17 处。**因此本主题应先按多数体例对齐后再改，否则会引入反向错误。**

**全部 id**（147 条）：

chebi-2995 tiers-60614 chebi-2319 tiers-73893 tiers-137455 tiers-119334 tiers-29650 tiers-8644
chebi-414 chebi-2282 chebi-807 chebi-412 tiers-44863 chebi-2299 chebi-1946 chebi-185
tiers-54007 tiers-83484 chebi-1148 chebi-1650 chebi-2017 chebi-1414 tiers-35976 chebi-377
chebi-2530 chebi-595 chebi-2552 chebi-717 chebi-3055 chebi-464 chebi-3214 chebi-410
chebi-1361 chebi-2345 chebi-2897 chebi-1750 chebi-3000 chebi-3241 chebi-2450 chebi-879
chebi-726 chebi-2573 chebi-173 chebi-860 chebi-215 tiers-154405 tiers-126756 chebi-1952
chebi-265 chebi-2048 tiers-94445 chebi-657 chebi-2137 chebi-2260 chebi-2983 tiers-121294
tiers-110482 tiers-108439 chebi-2818 chebi-689 chebi-828 tiers-91641 chebi-718 chebi-905
chebi-1070 chebi-764 tiers-104153 chebi-1245 tiers-102815 tiers-153455 tiers-132177 tiers-103190
chebi-3106 chebi-3188 chebi-3250 tiers-135315 tiers-96177 chebi-1072 chebi-926 chebi-43
chebi-2773 chebi-3126 tiers-156464 chebi-951 chebi-1473 chebi-1757 chebi-2021 chebi-27
chebi-605 chebi-16 chebi-1724 tiers-156410 chebi-1884 chebi-2885 chebi-1033 tiers-108673
chebi-2700 chebi-454 chebi-2146 chebi-2390 chebi-2254 chebi-353 chebi-598 tiers-84342
chebi-2157 chebi-1668 chebi-429 chebi-99 chebi-1274 chebi-1926 chebi-883 chebi-111
chebi-777 chebi-1722 chebi-1838 chebi-1546 chebi-1716 chebi-1987 chebi-569 chebi-217
chebi-1354 chebi-781 chebi-1005 chebi-348 chebi-861 chebi-819 chebi-1682 chebi-482
chebi-401 chebi-442 chebi-536 chebi-427 chebi-235 chebi-2209 chebi-2579 chebi-2845
chebi-2131 chebi-2082 chebi-697 chebi-3107 chebi-2259 chebi-1425 chebi-1418 chebi-1536
chebi-1565 chebi-349 chebi-2866

---

### G1 碳酰/含氮衍生物保留名母体缺失 — 132 条

**症状**：分子含 `N-C(=O)-N`（尿素）、`O-C(=O)-N`（氨基甲酸酯/ Boc / Cbz）、`O-C(=O)-O`（碳酸酯）、`N-C(=N)-N`（胍/脒）、`N=C-O-`（亚胺酸酯）、`S-C≡N` / `N=C=S`（硫氰酸/异硫氰酸酯）、`N-N`（肼）、`N-OH`（羟胺）时，PRED 把该 C1/N1 单元当开链母体，输出 `formamide` / `formate` / `formonitrile` / `iminomethylamino`；GOLD 用保留名 `urea` / `carbamate` / `carbonate` / `guanidine` / `amidine` / `imidothioate` / `thiocyanate` / `hydrazine` / `hydroxylamine`。作取代基时同理退化成 `amino]-oxomethyl` 而非保留前缀 `carbamoyl`。

**根因**（本报告**置信度最高**的一处，多块独立收敛）：
- `layer1/fg_registry.py::FG_SPECS` 只登记 14 个经典类（radical/acyl/acid/oxoacid/sulfonamide/ester/acyl_halide/amide/nitrile/aldehyde/ketone/alcohol/thiol/amine），**没有** urea / carbamate / carbonate / oxalate / amidine / guanidine / thiocyanate / isothiocyanate / hydroxylamine / hydrazine / imidate / sulfamate 等保留名类。
- 于是 `layer2/principal.py::PRINCIPAL_REGISTRY` 选不中它们，`layer2/principal_expression.py` 只能按 amide/ester/nitrile 表达，最后落到 `constants.CHAIN_RETAINED["amide"][1] = "formamide"` / `["ester"][1] = "formate"` 的 C1 兜底词干。
- `layer1/fg_local_smarts.py:40` 的 amide SMARTS 显式排除了 `N-O-C`（氨基甲酸酯氮）与环内 N，被排除后没有任何类别接管这个 C=O。
- `tools/anchored_table.py` 里其实已有 `carbamoyl` 条目，但只在酰卤路径使用，取代基命名路径不命中。

**修复路径**：这是一处**纯数据驱动**的修复——往 `FG_SPECS` / `CHAIN_RETAINED` / `anchored_table` 补保留名条目即可一次覆盖上百条，与项目「扩展先落数据表、不新加 worker」的原则一致。

**IUPAC 依据**：P-66.1.6.1.1.1（urea 为 PIN，N/N′ 位次）、P-65.2.1.1 / P-66.3.1（carbamic acid / carbamate / carbonate）、P-66.1.1.4（thiocyanate / isothiocyanate）、P-66.1.1.4.1.1（carbamoyl 保留前缀，含 `3-(dimethylcarbamoyl)pentanedioic acid` 等 PIN 例）、P-66.4（amidine / guanidine）、P-66.1.5.1（hydrazine / hydroxylamine）、P-66.5（imidothioate）。

**全部 id**（132 条）：

chebi-1986 chebi-2296 chebi-1152 chebi-2066 chebi-3136 tiers-9441 chebi-1124 chebi-709
tiers-30618 chebi-2007 chebi-3020 chebi-2577 chebi-2798 chebi-475 chebi-2008 chebi-963
chebi-1907 chebi-2011 tiers-74946 chebi-1369 tiers-76866 chebi-2401 tiers-72453 chebi-1766
tiers-105083 chebi-1175 tiers-142794 chebi-1755 tiers-158280 chebi-3067 tiers-28290 chebi-653
chebi-1849 chebi-2042 chebi-2920 chebi-493 tiers-13836 chebi-2063 chebi-90 tiers-76028
tiers-8088 tiers-88220 tiers-67481 chebi-312 tiers-130261 chebi-815 chebi-2002 tiers-14669
tiers-14524 chebi-396 tiers-36233 chebi-1846 chebi-2865 tiers-60725 chebi-119 chebi-834
tiers-20348 tiers-54463 chebi-153 chebi-1326 tiers-19564 chebi-2031 chebi-669 chebi-1180
chebi-538 chebi-1215 tiers-123881 chebi-518 chebi-212 chebi-1680 chebi-965 chebi-136
chebi-2115 chebi-508 chebi-374 tiers-66041 chebi-2113 chebi-2795 chebi-1900 chebi-2276
tiers-100036 chebi-1375 chebi-561 chebi-2230 chebi-192 chebi-2722 chebi-1573 chebi-2385
chebi-218 chebi-1941 tiers-125604 tiers-63826 tiers-143551 chebi-2274 chebi-1272 chebi-2205
chebi-2452 tiers-100116 tiers-118482 chebi-1980 tiers-132818 tiers-57528 chebi-2978 chebi-1195
chebi-1550 chebi-3081 chebi-2267 tiers-11629 tiers-151693 chebi-1016 chebi-2067 chebi-2884
tiers-125022 chebi-1265 chebi-3261 chebi-513 chebi-2093 chebi-3014 tiers-5030 chebi-2807
chebi-729 tiers-148773 chebi-3021 chebi-2499 chebi-1037 chebi-1316 chebi-3265 chebi-56
chebi-1549 chebi-2374 chebi-400 chebi-1429

---

### G8 环编号方向与位次 tie-break — 84 条

**症状**：骨架、母体、取代基集合全对，只是**编号方向取反**，位次集与立体描述符随之整体镜像（糖环 `(2R,3S,4S,5R,6R)` ↔ `(3S,4S,5R)`；肌醇磷酸酯 `(2R,3R,5S,6R)` ↔ `(2R,3S,5R,6R)`；菲类 `4aR,4bS,7S,10aR` ↔ `4aR,4bR,6S,8aS`）。子型还包括唑类 `3-`↔`5-`、`4-`↔`5-`，双键位次 `(3E)`↔`(1E)`，以及前缀引用位次未取最低。

**根因**：`layer4/numbering.py` / `numbering_engine.py` / `fused_numbering.py` / `ring_geometry.py` 在「两条方向位次集相同时」的 tie-break 与 IUPAC 不一致——IUPAC 是两级：先让**后缀特征基**拿最低位次，再按 **P-14.5 字母序**让首个前缀拿最低位次；PRED 疑似按「第一个/最大取代基」取位次。糖环还缺「C1 为异头碳」的保留编号（P-64.1.1.4），改用通用最低位次规则后起点选到 `2-(hydroxymethyl)` 一侧。

**注意**：本主题里有一批是**两个名字描述同一结构**（如 chebi-582 / chebi-2666 的肌醇环，两个方向位次集都是 {2,3,5,6}，最低位次规则判不出方向），属体例差异而非可修 bug。

**IUPAC 依据**：P-14.4(b)(c)(f)(g)(j)（最低位次集与并列时的字母序 tie-break）、P-14.5（字母数字序）、P-64.1.1.4（糖/环醇编号起点）、P-34（糖的立体与位次）。

**全部 id**（84 条）：

chebi-2162 tiers-153 tiers-75 chebi-2814 chebi-1891 chebi-1526 chebi-3019 chebi-485
tiers-9298 chebi-1168 chebi-1614 chebi-496 chebi-2072 tiers-47406 chebi-468 chebi-2335
chebi-208 chebi-341 chebi-2277 chebi-1771 chebi-1858 chebi-287 tiers-148904 chebi-3252
chebi-642 chebi-2667 chebi-935 chebi-2402 chebi-358 chebi-1131 chebi-247 chebi-95
tiers-134235 tiers-111140 chebi-1748 chebi-992 chebi-2337 tiers-153177 chebi-472 chebi-899
chebi-686 chebi-3242 chebi-996 chebi-2771 chebi-70 chebi-2855 chebi-1201 chebi-1212
chebi-1284 chebi-692 chebi-329 chebi-1831 chebi-2516 chebi-1344 chebi-1445 chebi-1413
chebi-408 chebi-562 chebi-3012 chebi-379 chebi-756 chebi-2356 chebi-1320 chebi-2029
chebi-2 chebi-2932 chebi-1704 chebi-2168 chebi-3253 chebi-1339 chebi-2536 chebi-2747
chebi-509 chebi-124 chebi-982 chebi-2396 chebi-559 chebi-953 chebi-1476 chebi-103
chebi-771 chebi-2751 chebi-352 chebi-1544

---

### G5 指示氢与氢化程度 — 81 条

**症状**：双向偏差。
- **漏标**：`4-oxo-1H-pyridine` 写成 `4-oxopyridine`；`indene-1,3-dione` 漏 1H（1 位已是 =O，其实该省）。
- **多标**：`1H-7H-purin-9-yl`（7 位 N 已被取代，不该有 7H）；`9H-7H-purin-3-yl` 双写。
- **位次错**：`5H-` ↔ `4H-`、`2H-tetrazol-5-yl` ↔ `1H-tetrazol-5-yl`。
- **氢化程度错**：`dodecahydro-1H-` ↔ `decahydro-1H,2H,10H`；给已是羰基的 C2/C5 再加 `2,5-dihydro`；`pyrazolidine-3,5-dione`（C4 是 sp³）写成 `pyrazole-3,5-dione`。

**根因**：`layer4/indicated_hydrogen.py::indicated_hydrogen` 只按「环内全单键且带 H」判定（`saturated_ring_atoms`），未与三件事联立：① 该位已被 =O/取代基占据时应省略；② 稠合环系的等效位（嘌呤 7H/9H）；③ 同环多候选时的取低规则。`layer5/assembler.py::_indicated_h_prefix` 直接拼接，多候选时产生 `3H-7H-`。嘌呤的 `7H-` 还叠加了 `layer2/ring_scaffold.py:126` 的 `locant_prefix` 硬编码与 `prefix_nh_conditional` 的重复叠加。

**IUPAC 依据**：P-14.4(b)（未取代化合物的指示氢）、P-14.4(d)（附加指示氢，`3,4-dihydronaphthalen-1(2H)-one`）、P-14.7 / P-25.7.1.3、P-31.1.4 / P-31.2.3（hydro 前缀位次取最低）、P-44.4.1.4（指示氢位次最低者优先）。

**全部 id**（81 条）：

chebi-2314 chebi-2553 chebi-1823 chebi-2433 chebi-1533 chebi-38 chebi-1541 chebi-2231
chebi-2645 chebi-1634 chebi-1514 chebi-1451 chebi-1470 chebi-2737 chebi-2673 chebi-2006
chebi-1468 chebi-1039 chebi-576 chebi-1237 chebi-801 chebi-3223 chebi-3127 chebi-817
chebi-486 chebi-1250 chebi-2601 chebi-2629 chebi-506 chebi-1862 chebi-2744 tiers-53380
chebi-2675 chebi-3169 chebi-1611 chebi-1795 chebi-302 chebi-3280 chebi-1934 chebi-836
chebi-36 chebi-1924 chebi-843 chebi-2903 chebi-3297 tiers-117610 chebi-1808 chebi-975
chebi-446 chebi-37 chebi-462 chebi-3257 chebi-1024 chebi-3144 chebi-1621 chebi-2358
chebi-1595 chebi-2646 tiers-125309 chebi-463 chebi-784 chebi-611 tiers-79323 chebi-1002
chebi-1278 chebi-2908 chebi-350 chebi-3041 chebi-1619 chebi-3183 chebi-1203 chebi-1289
chebi-3194 chebi-2353 chebi-3173 chebi-2545 tiers-101101 chebi-1562 chebi-3061 chebi-1171
chebi-1965

---

### G6 取代基引用顺序（P-14.5 字母序）— 72 条

**症状**：所有原子与基团都对，只是**前缀的引用顺序**不同。GOLD 按 P-14.5 字母数字序排列简单前缀（忽略位次数字、立体描述符与连字符），PRED 按位次或抽取顺序排。例：GOLD `5-hydroxy-3-(hydroxymethyl)-2-[3-hydroxy-4-…phenyl]chromen-4-one`（3 在 2 前），PRED 反过来。

**根因**（单点，修一处即可）：`layer5/assembler_prefixes.py::_build_prefix()` 里 `_mult_rows(..., sort_key=alkyl_alpha_key)` —— **排序键只对烷基生效**（代码注释即写「分组/排序键由调用侧给定，不与英文侧统一」），非烷基前缀退回位次序。且 `tools/re.py::alkyl_alpha_key()` 只剥开头的位次/立体/括号，**不剥词干内部的位次与连字符**，于是 `hydroxy-2-(octadecanoylamino)octadec-4-enoxy` 会以 `hydroxy-` 参与比较，与 GOLD「只比字母」口径不一致。

**顺带发现的关联 bug**：改动排序时若不同步重排位次串，会产生**同一环位次被两个取代基占用**的自相矛盾名字（全库检出 33 处，如 `2-carboxy-…-4,5-dihydroxyoxan-2-yl`）。修排序时必须一起处理。

**IUPAC 依据**：P-14.5.1（简单前缀按字母顺序排列，倍增前缀不改变已建立的顺序）、中国化学会《有机化合物命名原则》§5.5.2。

**全部 id**（72 条）：

chebi-1198 chebi-2627 chebi-852 chebi-2187 chebi-2310 chebi-2434 chebi-2293 chebi-2878
tiers-121960 chebi-123 tiers-121966 chebi-1818 chebi-2810 chebi-753 chebi-1690 chebi-2727
chebi-201 chebi-1109 chebi-737 chebi-987 chebi-1552 chebi-593 chebi-225 chebi-1939
chebi-2611 chebi-902 chebi-2458 chebi-322 chebi-1827 chebi-886 chebi-2835 chebi-1424
chebi-2171 chebi-1792 chebi-1840 chebi-2518 chebi-735 chebi-2418 chebi-2929 chebi-3088
chebi-3113 chebi-2635 chebi-2666 chebi-725 chebi-1582 chebi-2719 chebi-49 chebi-1632
chebi-2988 chebi-1830 chebi-1870 chebi-347 chebi-1247 chebi-838 chebi-1629 chebi-1097
chebi-1720 chebi-154 chebi-101 chebi-649 chebi-76 chebi-979 chebi-1399 chebi-2813
chebi-3065 chebi-650 chebi-54 chebi-1772 chebi-3125 chebi-1977 chebi-388 chebi-1446

---

### G10 盐、抗衡离子与电荷态 — 68 条

**症状**：三类。
1. **抗衡离子丢失**：`dichloride` / `tetrasodium` / `calcium bis(…)` / `hydrobromide` / `perchlorate` 整段不见；GOLD `calcium bis(9-[…]oxynonanoate)`，PRED 只剩有机阴离子。
2. **两性离子电荷不表达**：GOLD `azaniumyl…carboxylate`，PRED 中性 `amino…carboxylic acid`；GOLD `sulfate`，PRED `hydrogen sulfate`。
3. **阴离子后缀 vs 前缀**：分子同时含 COOH 与去质子 O⁻（或酚氧负离子 + 酮）时，GOLD 用 `-olate` 作后缀、COOH 作 `carboxy` 前缀；PRED 反过来用 `-oate` 作后缀、O⁻ 降为 `oxido` 前缀，编号方向随之取反。

**根因**：`layer0/salt.py::dissociate_salt()` 只认碱金属/HCl 片段（`_alkali_en`、`_is_hcl_frag`）且要求 `len(organics) == 1`，Ca²⁺、卤素阴离子之外的抗衡离子、多有机组分一律返回 None；`layer0/charge.py` 只在「酰胺 O⁻ 受体」这一条路径上传递电荷，去质子氧 / 两性离子的电荷不进入特征基清单，于是 `layer2/principal.py` 看不到 `olate` 这类阴离子后缀，只能退成 `oxido-` 前缀。

**注意**：酚盐作后缀、羧酸根作前缀是 **CHEBI 自身体例**，与标准 PIN 取向相反（PIN 更倾向羧酸根作后缀）。修此主题需先与 benchmark 体例对齐，盲修会引入反向错误。

**IUPAC 依据**：P-73（盐与离子化合物的命名）、P-74（阴离子命名 —— 只命名实际存在的电荷）、P-72.2.2.2.2（`-olate` 复合后缀、phenolate 保留名）、P-63.2（`sodium 3-hydroxypropan-1-olate`）。

**全部 id**（68 条）：

chebi-2628 chebi-137 chebi-3160 chebi-792 chebi-3240 chebi-691 chebi-1996 chebi-1564
chebi-2287 chebi-2496 chebi-3161 chebi-384 chebi-10 chebi-1594 chebi-833 chebi-1496
chebi-585 chebi-1508 chebi-1971 chebi-25 chebi-2442 chebi-884 chebi-2001 chebi-622
chebi-1932 chebi-2161 chebi-930 chebi-1765 chebi-2265 chebi-3060 chebi-1861 chebi-280
chebi-438 chebi-2883 chebi-1847 chebi-2782 chebi-2283 chebi-3091 chebi-2213 chebi-1646
chebi-186 chebi-1714 chebi-1078 chebi-1669 chebi-1442 chebi-553 chebi-2503 chebi-117
chebi-2182 chebi-3006 chebi-1330 chebi-2994 chebi-3016 chebi-2315 chebi-365 chebi-445
chebi-2891 chebi-1350 tiers-77209 chebi-549 chebi-2262 chebi-2176 chebi-2841 chebi-1317
chebi-2797 chebi-2278 chebi-958 chebi-2877

---

### G2 保留前缀与官能团表达缺失 — 43 条

**症状**：本应走保留前缀的基团退化成通用拼装式：
- `-N=N-` → `imino` + `amino` 拼成 `…iminoamino`（亚氨氨基），GOLD 用 `diazenyl`
- `-S-S-` → `sulfanyl]sulfanyl`，GOLD 用 `disulfanyl`
- `C=S`（环内或链上）→ `sulfanylidene`（硫烷亚基），GOLD 用 `-thione` / `thioxo`
- `-C(=N-OH)-S-` → `sulfooxyimino…sulfanyl`，GOLD 用 `imidothioate`
- `-SO₂-N<` → `anilino]sulfonyl`，GOLD 用整体前缀 `phenylsulfamoyl`
- `C(=N-OH)-` 作取代基 → `…ylimino…`，GOLD 用 `carbonimidoyl`

**根因**：`layer1/fg_local_smarts.py` 没有 azo（`–N=N–`）规则，N=N 两端被 imine / amine 规则分别吃掉，`tools/anchored_table.py:93` 的 `diazenyl` 锚点（`*N=N`）因此永远命不中；`anchored_table.py:39` 的 `sulfanylidene` 叶子以 `anchored=("*=S",)` 抢占了所有 C=S。缺条目清单：`diazenyl`、`disulfanyl`、`carbonimidoyl`、`benzenesulfonyl`（优于 `phenylsulfonyl`）、`methanesulfonamido`、`carbamothioylamino`。

**IUPAC 依据**：P-68.3.1.3（`diazenyl`，且明确「以 diazenyl 为基础的名字优先于 azo」）、P-66.1.2 / P-64.1.2（`thioxo` / `-thione`）、P-68.3.1.4（`disulfanyl`）、P-66.1.1.4（`carbonimidoyl`）、P-66.4（`sulfamoyl` 整体前缀）、P-66.1.1（保留取代基名）。

**全部 id**（43 条）：

tiers-73757 chebi-409 chebi-2950 chebi-2856 tiers-58291 tiers-148854 chebi-983 tiers-23153
tiers-71605 tiers-25538 chebi-2290 tiers-59838 tiers-50018 chebi-1478 tiers-74062 tiers-880
chebi-3279 chebi-2519 chebi-1982 tiers-92413 chebi-1084 chebi-2973 chebi-3128 chebi-2661
tiers-97781 tiers-152351 chebi-2526 chebi-960 chebi-386 tiers-92497 chebi-1740 chebi-3256
chebi-581 tiers-126705 chebi-2229 chebi-112 chebi-1797 chebi-1643 chebi-1518 chebi-1137
tiers-109289 chebi-2440 chebi-2669

---

### G12 主链与母体骨架选择 — 42 条

**症状**：两条候选骨架都合法，GOLD 与 PRED 各选一边：
- **多羧酸**：羧基 ≥3 时 GOLD 以「含最多羧基的链」为母体并用 `-tricarboxylic acid` 后缀，PRED 选最长链 + `carboxy` 前缀 + `-dioic acid`。
- **链 vs 环 / 环 vs 环**：`benzo[a]quinolizine` 类保留名 vs 系统稠合名；胺的甲基/乙基母体；肽链选哪一端。
- **两个等价环**：如萘醌/酚氧负离子体系里哪个环当母体。

**根因**：`layer2/parent_skeleton.py::select_principal_skeletons` 的骨架评分只看链长/取代基，**不统计「同一主特征基团出现次数」**；`layer2/parent_select.py` 未实现 P-44/P-45 的逐级判据（特征基最多 → 骨架原子最多 → 环数 → 取代基最多 → 位次最低）；`layer2/chain_walk.py` 的链端选择另有独立打分，两处不一致。

**IUPAC 依据**：P-44.3.2（母体氢化物选择顺序）、P-45（链选择）、P-65.2.1（多元羧酸以含最多羧基的链为母体）、P-41.x（环系优先于链）。

**全部 id**（42 条）：

chebi-1721 chebi-652 chebi-2479 chebi-852 chebi-159 chebi-3281 chebi-785 chebi-1585
chebi-984 chebi-2848 chebi-336 chebi-923 chebi-2134 tiers-53012 tiers-12588 chebi-1153
chebi-41 chebi-603 tiers-37791 chebi-2144 chebi-1839 chebi-575 tiers-153130 chebi-2224
chebi-2091 chebi-542 tiers-146666 tiers-127007 tiers-110847 tiers-112315 chebi-1671 chebi-2809
chebi-2693 chebi-2741 tiers-116162 chebi-672 chebi-648 chebi-1189 chebi-1555 chebi-191
chebi-694 chebi-378

---

### G11 立体描述符 — 29 条

**症状**：R/S 缺失、位次错或整体翻转。
- **缺失**：非环骨架上的手性碳（侧链/苄位）不给描述符，如 `(S)-hydroxy(3,4,5-trimethoxyphenyl)methyl` 里的 `(S)`；亚砜硫的 `(S)`。
- **位次错**：`(5R,8R)-2,4-dihydroxynonyl` —— 位次 5,8 与词干「2,4-二羟基壬基」自相矛盾（该用母体链位次 2,4）。
- **C=N 位次**：肟/腙的 E/Z 前缀用了侧链局部编号，未换算回母体位次（`(3E)` ↔ `(1E)`）。
- **整组翻转**：糖环/多环编号方向一变，R↔S 全反（与 G8 同源）。

**根因**：`layer5/stereo.py` 的 `_cip_on_chain` / `_rs_parts` 只覆盖母体链与环上的手性中心，侧链取代基内部的立体中心未收集；`layer4/locant_calc.py` 把立体位次写成分子内原子序号而非母体编号位次。

**IUPAC 依据**：P-91 / P-92（立体描述符须完整且用母体编号的位次）、P-93.4（E/Z）、P-93.5（硫手性）、P-14.4(j)。

**全部 id**（29 条）：

tiers-149669 chebi-818 tiers-29582 chebi-621 chebi-1887 tiers-54380 chebi-712 chebi-134
chebi-994 chebi-1874 chebi-647 chebi-1116 chebi-2018 chebi-2943 chebi-2430 chebi-1592
chebi-2336 chebi-3070 chebi-2451 chebi-582 chebi-2583 chebi-545 chebi-1249 chebi-2199
chebi-875 chebi-2654 chebi-609 chebi-1617 chebi-623

---

### G13 中文侧（EN 通过、ZH 失败）— 12 条

**症状**：英文与 GOLD 完全一致，只有中文不同。子类：① 保留母体中文名表缺条目（`色满`/`异色满`）；② 缺「基」字（`吡啶-4-氧基` ↔ `吡啶-4-基氧基`）；③ 中文前缀次序（N-取代基的位置）；④ 括号/连词策略（`[(4-硝基苯基)硫基]` ↔ `(4-硝基苯硫基)`）；⑤ 词序（`乙酰基乙酸酯` ↔ `乙酸乙酰酯`）。

**根因**：`layer5/assembler.py` / `chain_engine.py` / `assembler_prefixes.py` 的中文侧组装表（后缀词、连接词「基」、括号策略、保留母体中文名）与英文侧不同源，是两套并行逻辑。

**IUPAC 依据**：中国化学会《有机化合物命名原则》第 2/4/5 章（中文后缀与连接词）；`docs/iupac/cn/第5章_命名实施导引.md` §5.5.2。

**全部 id**（12 条）：

tiers-93 tiers-180 tiers-28045 tiers-1819 tiers-48606 tiers-18034 tiers-4449 tiers-67516
tiers-113808 tiers-20508 tiers-119391 tiers-127997

---

## 4. 范围外主题（只计数，本期不修）

按用户要求，以下两个大主题只统计不展开：

### 4.1 桥环 / 螺环

各分块上报：chunk_01 26、chunk_02 26、chunk_03 63、chunk_04 35、chunk_05 22、chunk_06 36，合计约 208 条（含重叠）。表现为 PRED 空输出、或用稠合名替代 von Baeyer 名、或结构解读错误（环氧小环被并进大环）。根因集中在 `layer2/ring_expression_policy.py::_unsupported_typed_ring`（拦截 → 空输出）、`layer2/fused_system.py`、`layer1/ring_systems.py`。

### 4.2 铵 / 鎓命名

各分块上报：chunk_01 40、chunk_02 18、chunk_03 19、chunk_04 15、chunk_05 9、chunk_06 32，合计约 133 条（含重叠）。表现为 GOLD 一律把 N⁺/S⁺ 当母体（`…azanium` / `…sulfonio…` / `…-1-ium`），PRED 一律把胺/酸当母体；或质子化态选择相反（`azaniumyl` ↔ `amino`）。根因集中在 `layer0/charge.py`（质子化/两性离子判定）、`layer5/assembler_prefixes.py`。

---

## 5. benchmark 数据不一致

**这一节不是 PRED 的错，而是评测集自身的问题** —— 它构成准确率的天花板，修代码无法触及。六块分析独立发现的类型如下：

### 5.1 GOLD_EN 与 GOLD_ZH 互不对应

| id | 证据 |
|---|---|
| chebi-1690 | EN 是 C20（`icosa-…`），ZH 写成「二十四碳」（C24） |
| chebi-2041 / chebi-2707 | ZH「二十三碳」↔ EN `henicosa`（C21） |
| chebi-2667 | EN `4-methylsulfinylphenyl`（亚砜）↔ ZH「4-甲磺酰基苯基」（砜） |
| chebi-2110 | ZH「2-硫代苯基」应为「噻吩-2-基」（`thiophen-2-yl`） |
| chebi-186 | ZH「四苯并环丁烯」与 EN `tetracen` 不对应且非标准名 |
| chebi-219 | EN `thiopyrano[4,3-d]pyrimidin-4-one` ↔ ZH「噻吩并[4,3-d]嘧啶」（噻喃并 ≠ 噻吩并） |
| chebi-1459 / chebi-846 | EN `azulene`（薁）↔ ZH「蒽醌」/「氮杂蒽」，中文母体完全错 |
| chebi-3055 | EN `phenazine-2,8-diamine` ↔ ZH 丢二胺、且「菲嗪」≠「吩嗪」 |
| chebi-2221 | EN `hexahydro-1H-azulene` ↔ ZH「六氢-1H-氮杂蒽」 |
| chebi-2296 | EN `sulfanylazaniumylidynemethane` ↔ ZH「硫代异氰酸」，指向不同物种 |
| chebi-792 / 3240 / 384 / 2287 / 2496 / 1996 / 3161 | EN 用阴离子后缀 `-olate`，ZH 写成中性「醇/苯酚」，电荷态不对应 |
| chebi-486 | EN `naphthalene-2-sulfonate`（阴离子）↔ ZH「萘-2-磺酸」（酸） |
| chebi-2984 一组（picene 类） | EN `picene` 在 ZH 里被译成「芴」（fluorene，**另一个环系**）或音译「比森」，chebi-1580 又译成「芘」——同一环系三种互不相容的译名 |
| chebi-1997 / 2363 / 2944 / 679 / 2826 / 1601 | ZH 把 `azaniumyl` 误译为「氮烯/氮烯基」（应为「铵基」） |
| chebi-117 | ZH 把 `oxido` 误译为「氧代」（= `oxo`） |
| chebi-2891 | ZH 明显乱码重复；EN 电荷态与输入 SMILES 不符（结构全中性却命名成氯化物盐） |
| chebi-1442 | ZH 漏译 `-ium` 与位次 |
| chebi-930 | ZH「4-氧代-7-氧代基-2H-1-苯并吡喃」重复字，母体名与 EN 不一致 |
| chebi-694 | ZH 出现「…氧基-5-氧代-2H-色烯-5-酸」——两个 5 位且以「-5-酸」结尾，不可解析 |
| chebi-3106 / chebi-351 | ZH 括号不配对 / 结尾被截断 |
| chebi-173 vs tiers-126756 | 同一 benzodioxole 骨架，两处 ZH 用不同的环名译法 |

### 5.2 GOLD_EN 与给定 SMILES 结构不符 / 非 IUPAC

| id | 证据 |
|---|---|
| chebi-2506 / chebi-554 | GOLD 写 `2,3-dihydroinden-1-one` 省掉 `1H-`，而同批 chebi-1024 / 1621 / 463 又保留指示氢 —— **体例自相矛盾**，PRED 的 `2,3-dihydro-1H-inden-1-one` 更合规 |
| chebi-1909 / chebi-2348 | GOLD `1H-imidazol-4-yl` 与该 SMILES 的互变异构不符（用组氨酸标准 SMILES 对照验证），PRED 的 `5-yl` 才对 |
| chebi-2825 / chebi-1355 | SMILES 是亚胺（`c(=N…)`），GOLD 命名的是对应的胺 —— 命名了另一个互变体 |
| chebi-1359 | SMILES 的 `[nH+]` 在喹啉 N 上，GOLD 把正电荷放在苄位二乙氨基 |
| chebi-2409 | SMILES 是 2 K⁺ + 二价阴离子，GOLD 写 `…carboxylate hydroxide`，结构中不存在 hydroxide 且电荷不平衡 |
| tiers-147852 | GOLD `[1]benzothiol` /「苯并噻咯」均非推荐名（应为 `1-benzothiophene` / 苯并噻吩） |
| chebi-287 | GOLD 用 `1H-pyrazole-5-carboxamide`，按 P-14.4(c) 主特征基应得最低位次（3-），PRED 更合规 —— 需人工裁决 |
| chebi-2814 / chebi-1891 | 咪唑 4/5 位争议，若 GOLD 取自 NH 在另一 N 的互变异构体则自洽 —— 需人工确认 |
| chebi-191 / chebi-694 | GOLD 把色酮当取代基、以单环 phenol/phenolate 为母体，**违反 P-44.2「环数多者为优先母体」**；PRED 的母体选择其实更合规 |
| chebi-2883 / 1847 / 3091 / 186 | 分子内同时有酮/酰胺（优先于酚），GOLD 却把 `-olate` 当后缀母体 —— 与 P-64 特征基次序冲突（**成族**） |

### 5.3 GOLD 排版体例不统一（不改变化学对错，但直接拉低字符串指标）

- **`oyl]oxy` vs `oyloxy`**：同一片段 `(9Z)-octadec-9-enoyl` + `oxy`，chebi-1722 写 `enoyloxy`、chebi-122 写 `enoyl]oxy`、chebi-238 写 `enoyloxy`、chebi-2726 写 `enoyl]oxy`。全库 GOLD 43 处 `oyloxy`、17 处 `oyl]oxy`，而 PRED 一律 `oyl]oxy`。**一个出口不可能同时命中两种写法。**
- **多余/缺失连字符**：chebi-99 GOLD 写 `…-6-(4-methylpiperazin-1-yl)-pyridin-3-yl]propanamide`，同体例在 chebi-2810 与其余 20 余处均无此连字符 —— PRED 与该条的唯一差异就是这个连字符。tiers-151254 / chebi-207 同类。chebi-2137 反之。
- **CAS 式双圆括号**：tiers-100036 / 100116 / 118482 / 132818 用 `1-((1-(…)-5-oxopyrrolidin-3-yl)methyl)-3-…`，同批其它条目用方括号。
- **位次上标三种写法并存**：`五环[9.7.0.0¹,³.0³,⁸.0¹²,¹⁶]`（上标）与 `01,3.03,8`（平文）混用；chebi-1243 的 `tricyclo[22.2.2.211,14]` 上标与主数字粘连、无法解析。
- **N 位标写法**：chebi-1526 / 2072 / 2804 用自造的 `2-N-` / `2-N,4-N-` / `6-N-`，ZH 却写 `N-`；chebi-753 / 3055 用 `2-N`/`4-N`，其余条目用 `N,N'`。
- **指示氢省/留不统一**：chebi-1278 / 1289 写 `furan-2-one` / `pyrrole-2,5-dione`（省），chebi-103 又写 `5-oxo-2H-furan-3-yl`（留）。
- **同位素标记**：chebi-1282 SMILES 含 `[11CH3]`，GOLD 中英均未标 `(¹¹C)`。

### 5.4 建议

修代码前先裁决 5.1 / 5.2 中的条目。可机器判定的部分（同结构不同 gold 名、中英结构不对应）已经核过：**无重复 id、无「同一结构不同 gold EN 名」的硬冲突**（3974 个唯一结构），所以不一致主要发生在**命名体例层**而非结构层。建议在动手改 G7（括号）与 G10（阴离子后缀）之前，先按多数体例把 GOLD 侧统一一遍，否则这两块每改一处都会在另一处产生反向错误。

---

## 6. 候选修复主题（待批准 5 个）

按「预期转正条数 / 改动成本」排序。**加粗**为推荐。

| # | 主题 | 预期转正 | 难度 | 单点程度 | 依据 |
|---|---|---:|---|---|---|
| 1 | **碳酰/含氮衍生物保留名补表**（G1） | ~110–130 | 中低 | 高（扩 `FG_SPECS` + `CHAIN_RETAINED` + `anchored_table` 数据表） | 六块独立收敛到同一处，且属项目「数据驱动优先」范式 |
| 2 | **含氧酸类名后缀与母体优先级**（G3） | ~130–180 | 中 | 中高（p41 等级表 + `_OXO_Z_ANCHORED` + 后缀出口） | 收益最大；P-44.1.2 优先原子序规则明确 |
| 3 | **括号/连字符/词界统一策略**（G7） | ~80–130 | 低 | 高（`assembler_prefixes` 一处策略） | P-1 §14 规则明确。**但必须先统一 GOLD 体例**（见 §5.3） |
| 4 | **取代基引用顺序改字母序**（G6） | ~50–70 | 低 | 极高（`_build_prefix` 排序键 + `alkyl_alpha_key`） | P-14.5 规则唯一、无歧义；顺带修 33 处位次重复占用 |
| 5 | **指示氢规则统一**（G5） | ~50–70 | 中 | 高（`layer4/indicated_hydrogen.py` 一处） | P-14.4(b)(d) + P-44.4.1.4 |
| 6 | 覆盖率强制校验（G9） | ~60–150 | 中高 | 中（`coverage.complete()` 加拦截，但会误伤） | 收益不确定，可能把「残名恰好命中」的少数条目算丢 |
| 7 | 环编号方向 tie-break（G8） | ~50–80 | 中 | 中（`layer4/numbering*` 多处） | P-14.4(f)(g) 两级 tie-break |
| 8 | 保留前缀补表（偶氮/二硫/硫酮/磺酰胺，G2） | ~30–43 | 低 | 极高（纯补 `anchored_table`） | P-68.3.1.3 等，与 #1 同类但更小 |

**推荐组合（覆盖约 420–550 条）**：**#1 + #2 + #4 + #5 + #6**（或把 #6 换成 #3）。

- 选 #1 / #2 的理由：收益最大且根因最集中，多块独立指向同一处。
- 选 #4 / #5 的理由：单点、规则无歧义、低风险，适合作为并行的「稳赢」块。
- #3（括号）虽然条数多，但 §5.3 已证明 GOLD 自身在这一维度上摇摆（`oyl]oxy` 43 vs `oyl]oxy` 17），**建议先做体例统一再做**，否则是负收益。
- #8（43 条、纯补表）可作为 #1 的附赠项一并纳入。

> 以上为**建议**。请指定 5 个主题，之后开 5 个 worktree 并行修复。

---

## 附录 A：全部 1526 个错例 id（按复杂度升序）

序号即复杂度名次（重原子数 → 环数 → 碳数）。每条 id 的 SMILES / GOLD / PRED 全文见 `tmp/analysis/appendix_ids.md`。

   1. chebi-1986 chebi-2296 chebi-1152 chebi-1380 chebi-2066 chebi-2162 chebi-2995 chebi-2055 chebi-86 chebi-3136 tiers-9441 tiers-60614
  13. tiers-73757 tiers-153 chebi-1124 chebi-2948 chebi-709 chebi-2628 tiers-30618 tiers-93 chebi-1282 chebi-2057 chebi-2334 chebi-2319
  25. chebi-2752 chebi-409 tiers-75 chebi-2007 chebi-3020 chebi-1396 chebi-2577 tiers-5188 chebi-2950 tiers-149669 chebi-2798 chebi-475
  37. chebi-789 chebi-1198 chebi-818 chebi-1012 chebi-895 tiers-180 chebi-3285 chebi-2814 chebi-2856 chebi-2008 chebi-1089 chebi-1786
  49. tiers-6167 chebi-971 chebi-2815 chebi-963 tiers-28045 chebi-137 chebi-1705 tiers-107394 chebi-2328 chebi-3160 chebi-320 chebi-624
  61. chebi-3133 chebi-3109 chebi-80 tiers-43418 chebi-318 tiers-58291 chebi-2314 chebi-1907 tiers-1819 tiers-29582 chebi-455 chebi-885
  73. chebi-2415 chebi-826 tiers-66561 chebi-3137 chebi-903 chebi-1686 chebi-730 chebi-1734 chebi-2421 chebi-1528 chebi-2627 chebi-621
  85. chebi-2011 tiers-148854 chebi-1891 chebi-2553 tiers-68805 chebi-792 tiers-74946 chebi-1721 chebi-1823 chebi-3151 tiers-73893 chebi-2433
  97. tiers-44245 chebi-1503 chebi-2147 chebi-2153 tiers-29902 chebi-983 chebi-1758 chebi-1887 chebi-2380 chebi-652 tiers-13388 chebi-2711
 109. tiers-137455 chebi-1526 chebi-1533 chebi-3019 chebi-38 tiers-119334 chebi-1369 chebi-3240 tiers-29650 chebi-3267 chebi-691 tiers-76866
 121. chebi-1841 tiers-54380 tiers-8644 chebi-3086 chebi-1541 chebi-2231 chebi-2645 tiers-6342 tiers-65526 chebi-414 chebi-2141 chebi-2479
 133. chebi-852 chebi-159 chebi-2401 chebi-854 chebi-2585 chebi-485 tiers-72453 chebi-1634 chebi-1514 chebi-1996 tiers-48606 chebi-1766
 145. chebi-3008 tiers-105083 chebi-82 tiers-916 chebi-1451 chebi-1919 chebi-1470 tiers-59592 chebi-155 tiers-23153 tiers-85289 chebi-2737
 157. tiers-52889 chebi-1175 chebi-1796 chebi-2673 tiers-132357 tiers-142794 tiers-17996 chebi-395 chebi-2282 chebi-2256 chebi-1755 chebi-1564
 169. chebi-3281 chebi-785 chebi-829 chebi-807 chebi-698 chebi-1695 chebi-1329 tiers-158280 tiers-51373 chebi-412 chebi-896 chebi-2287
 181. chebi-3067 tiers-28290 chebi-2496 chebi-653 chebi-572 tiers-71605 chebi-3161 tiers-9298 chebi-1147 chebi-2006 chebi-176 chebi-1585
 193. chebi-984 chebi-1849 chebi-2042 chebi-2920 chebi-493 tiers-25538 chebi-2290 chebi-2187 chebi-346 chebi-1914 chebi-384 chebi-2928
 205. chebi-458 tiers-13836 tiers-8543 chebi-10 tiers-18034 tiers-44863 chebi-2299 chebi-1168 chebi-1614 chebi-2221 tiers-96844 chebi-1157
 217. chebi-1256 chebi-2297 chebi-3018 chebi-2848 chebi-336 chebi-923 chebi-1178 chebi-2134 chebi-1946 chebi-2063 chebi-496 chebi-185
 229. chebi-2072 chebi-90 chebi-432 tiers-59838 tiers-76028 chebi-1468 chebi-1039 chebi-576 chebi-1237 chebi-2804 tiers-7355 tiers-8088
 241. chebi-1594 tiers-4449 tiers-50018 tiers-80645 chebi-801 tiers-18169 chebi-194 chebi-800 chebi-833 chebi-1842 chebi-2872 chebi-1348
 253. chebi-121 chebi-270 chebi-3223 tiers-67516 tiers-88220 chebi-1496 chebi-3127 chebi-585 tiers-130171 tiers-67481 tiers-76395 chebi-2310
 265. chebi-817 chebi-1478 chebi-486 tiers-10747 tiers-54007 chebi-312 tiers-68910 tiers-83484 tiers-130261 chebi-1148 chebi-298 chebi-846
 277. chebi-1888 chebi-1250 chebi-815 tiers-53012 chebi-1508 chebi-474 tiers-124039 tiers-12588 chebi-1053 chebi-1296 chebi-1236 chebi-1520
 289. chebi-2893 chebi-1650 chebi-1153 chebi-41 chebi-1436 chebi-907 chebi-1365 chebi-603 chebi-1971 chebi-1812 chebi-2002 tiers-14669
 301. chebi-2017 chebi-2527 chebi-2601 chebi-22 tiers-37200 chebi-2478 tiers-113808 tiers-49651 tiers-74062 chebi-2629 chebi-25 tiers-37791
 313. chebi-1293 chebi-1414 chebi-1820 tiers-14524 tiers-35976 chebi-108 chebi-1703 chebi-2517 chebi-506 chebi-377 chebi-1300 chebi-2530
 325. chebi-396 tiers-36233 chebi-1179 chebi-1846 chebi-2442 chebi-595 chebi-1081 chebi-2552 chebi-810 chebi-1862 chebi-2865 tiers-83257
 337. tiers-880 chebi-717 tiers-60725 chebi-1459 chebi-2138 chebi-3279 chebi-398 chebi-873 chebi-119 chebi-2144 chebi-65 chebi-167
 349. chebi-2550 chebi-2744 chebi-3055 tiers-53380 chebi-587 chebi-1338 chebi-2675 tiers-30123 chebi-2210 tiers-47406 chebi-468 chebi-2434
 361. chebi-1427 chebi-464 tiers-155306 chebi-3214 chebi-3169 chebi-3200 chebi-1125 chebi-834 tiers-20348 tiers-20508 chebi-2034 chebi-410
 373. chebi-1979 chebi-2136 chebi-1098 chebi-1611 chebi-1839 chebi-1361 chebi-884 chebi-1864 chebi-1588 chebi-1795 chebi-1570 chebi-2616
 385. chebi-1802 chebi-2811 chebi-1727 chebi-2477 chebi-594 chebi-2345 chebi-1524 chebi-1868 tiers-54463 chebi-1092 tiers-151254 tiers-19300
 397. chebi-153 chebi-3100 chebi-3054 chebi-1973 tiers-96746 chebi-219 chebi-1326 chebi-1560 chebi-2001 chebi-302 chebi-622 tiers-19564
 409. chebi-1932 chebi-2161 chebi-2897 chebi-3280 chebi-1386 chebi-3209 chebi-930 chebi-2106 chebi-1438 chebi-2576 chebi-1639 chebi-2108
 421. chebi-2594 chebi-1750 chebi-3000 chebi-316 chebi-575 chebi-2870 chebi-2031 chebi-1934 chebi-3241 chebi-669 chebi-1180 chebi-836
 433. chebi-447 chebi-36 chebi-1924 chebi-2519 chebi-538 chebi-843 chebi-1215 tiers-123881 chebi-518 chebi-212 chebi-2335 chebi-503
 445. chebi-207 chebi-2903 chebi-3297 tiers-117610 chebi-1765 chebi-1808 chebi-2450 chebi-2696 chebi-1680 chebi-1737 chebi-1145 chebi-2173
 457. chebi-1270 chebi-1312 chebi-1872 chebi-3131 chebi-1744 chebi-975 chebi-712 chebi-965 chebi-136 chebi-1684 chebi-2100 chebi-1982
 469. chebi-208 chebi-341 chebi-2115 chebi-2265 chebi-508 chebi-374 chebi-446 chebi-263 chebi-3060 chebi-1861 chebi-879 tiers-119391
 481. chebi-2293 chebi-2878 chebi-403 chebi-952 tiers-66041 chebi-1894 chebi-2207 chebi-2447 chebi-397 chebi-2427 chebi-342 chebi-1356
 493. tiers-146758 chebi-2143 chebi-1745 chebi-280 chebi-3153 chebi-438 chebi-2303 chebi-967 chebi-226 chebi-1228 chebi-2277 chebi-2113
 505. chebi-2795 chebi-2532 chebi-1900 chebi-2276 chebi-1771 chebi-87 chebi-726 chebi-37 chebi-1106 chebi-2311 chebi-2409 chebi-1967
 517. tiers-100036 tiers-128583 tiers-153130 chebi-2224 chebi-257 chebi-2573 chebi-3248 chebi-462 chebi-832 chebi-2148 chebi-465 chebi-626
 529. chebi-2327 chebi-1484 chebi-173 chebi-3257 chebi-1908 chebi-2257 chebi-741 chebi-2406 chebi-1099 chebi-2475 chebi-1048 chebi-860
 541. chebi-1375 chebi-561 chebi-2091 chebi-1858 chebi-2230 chebi-2255 chebi-1559 chebi-192 chebi-542 chebi-215 chebi-2722 chebi-1024
 553. chebi-1196 chebi-2883 chebi-604 chebi-287 chebi-289 tiers-146666 tiers-148904 tiers-154405 chebi-241 tiers-92413 chebi-1084 chebi-2825
 565. chebi-1359 chebi-188 chebi-471 chebi-606 chebi-1335 chebi-1567 chebi-1573 chebi-281 tiers-126756 chebi-2589 chebi-1305 chebi-3252
 577. chebi-3266 chebi-434 chebi-1873 chebi-3078 chebi-2485 chebi-557 chebi-1409 chebi-2537 chebi-3144 chebi-2678 chebi-2781 chebi-642
 589. chebi-1847 chebi-1952 chebi-2511 chebi-2385 chebi-218 chebi-1941 chebi-2242 chebi-2931 chebi-2164 chebi-393 chebi-265 chebi-2973
 601. chebi-134 chebi-2269 tiers-125604 tiers-127007 tiers-63826 chebi-3128 tiers-121680 chebi-2782 chebi-2022 chebi-424 chebi-2283 chebi-1621
 613. tiers-143551 chebi-1417 chebi-479 chebi-1791 chebi-360 chebi-667 chebi-2009 chebi-2274 chebi-2999 chebi-1850 chebi-2061 chebi-1698
 625. chebi-1673 chebi-2076 chebi-938 chebi-1362 chebi-3237 chebi-3288 chebi-2041 chebi-2036 chebi-2661 chebi-2871 chebi-2048 chebi-1272
 637. tiers-94445 chebi-657 chebi-2938 chebi-1389 chebi-2957 chebi-571 chebi-227 chebi-1086 chebi-3085 chebi-1248 chebi-2205 chebi-2596
 649. tiers-97781 chebi-162 chebi-2432 tiers-104849 chebi-1968 chebi-2452 chebi-3091 chebi-1540 chebi-2358 tiers-110847 chebi-139 chebi-2344
 661. chebi-3026 tiers-100116 tiers-118482 tiers-107887 tiers-91385 chebi-253 chebi-2667 chebi-3217 chebi-935 tiers-112315 chebi-2137 chebi-2110
 673. chebi-2674 chebi-2603 tiers-147852 chebi-822 chebi-1671 chebi-1980 chebi-1852 chebi-2260 chebi-2402 chebi-2653 chebi-2351 chebi-333
 685. chebi-369 chebi-452 chebi-661 chebi-916 chebi-3206 chebi-1731 chebi-2070 chebi-2213 chebi-2809 tiers-152351 chebi-1355 chebi-1485
 697. chebi-2685 chebi-1692 chebi-1646 tiers-102877 chebi-129 tiers-108409 tiers-133227 chebi-186 tiers-132818 tiers-98594 chebi-358 tiers-57528
 709. chebi-1723 chebi-19 chebi-2506 chebi-2819 chebi-2693 chebi-2707 chebi-2323 chebi-1595 chebi-2646 tiers-95650 chebi-197 chebi-2721
 721. chebi-2493 chebi-3229 chebi-2025 chebi-2755 chebi-2978 chebi-2526 chebi-2983 chebi-994 chebi-1195 chebi-1550 chebi-960 tiers-125309
 733. chebi-386 tiers-121294 chebi-1714 chebi-463 chebi-2128 chebi-784 chebi-2741 chebi-84 chebi-1131 chebi-1824 chebi-247 chebi-95
 745. tiers-110482 tiers-92497 chebi-3181 tiers-116162 tiers-134235 tiers-156243 tiers-98600 chebi-676 chebi-1874 tiers-150280 chebi-554 tiers-9744
 757. chebi-1240 tiers-108439 chebi-611 chebi-1087 chebi-1381 chebi-1962 chebi-2612 chebi-2307 chebi-672 chebi-1371 chebi-1740 chebi-648
 769. chebi-2818 chebi-689 chebi-2204 chebi-925 chebi-2918 chebi-1040 chebi-828 tiers-111140 chebi-782 chebi-546 chebi-2714 tiers-91641
 781. chebi-3053 chebi-526 chebi-64 tiers-79323 chebi-1078 chebi-2175 chebi-718 chebi-943 chebi-2151 chebi-3081 chebi-647 chebi-905
 793. chebi-3278 tiers-92000 chebi-1775 chebi-3140 tiers-140247 chebi-1002 chebi-1278 chebi-271 chebi-3193 chebi-915 chebi-1189 chebi-1555
 805. chebi-3256 chebi-607 chebi-855 chebi-1748 chebi-992 chebi-1999 chebi-2337 chebi-1070 chebi-1028 chebi-1204 chebi-279 chebi-430
 817. chebi-581 chebi-3 chebi-2267 chebi-764 tiers-104153 tiers-126705 chebi-2064 tiers-11629 chebi-1245 chebi-2229 tiers-102815 tiers-153455
 829. chebi-2908 tiers-131559 tiers-132177 tiers-151693 tiers-153177 chebi-1481 tiers-103190 chebi-3106 chebi-3188 chebi-350 chebi-477 chebi-204
 841. tiers-148964 chebi-3134 chebi-294 chebi-3225 chebi-3250 chebi-1837 chebi-849 chebi-921 chebi-472 chebi-899 chebi-686 chebi-1016
 853. chebi-2961 chebi-2591 chebi-112 chebi-1797 chebi-191 chebi-694 tiers-135315 tiers-96177 chebi-2723 chebi-968 chebi-3208 chebi-2067
 865. chebi-125 tiers-116774 chebi-1072 chebi-1787 chebi-2884 chebi-926 tiers-125022 chebi-2623 chebi-351 chebi-1669 chebi-3122 chebi-632
 877. chebi-2959 chebi-816 chebi-1401 chebi-295 chebi-1643 chebi-43 chebi-1265 chebi-2773 chebi-3041 chebi-1983 chebi-3242 chebi-3118
 889. chebi-3261 chebi-513 chebi-421 chebi-1116 chebi-2093 tiers-127997 chebi-1017 chebi-3126 tiers-121960 tiers-156464 tiers-96196 chebi-51
 901. chebi-951 chebi-568 chebi-929 chebi-1473 chebi-1757 chebi-2021 chebi-996 chebi-27 chebi-3032 chebi-605 chebi-3087 chebi-986
 913. chebi-1276 chebi-16 chebi-1990 chebi-2745 chebi-1576 chebi-1619 chebi-2953 chebi-3246 chebi-1115 chebi-1209 chebi-1843 chebi-1860
 925. chebi-378 chebi-1724 chebi-1430 chebi-2668 chebi-2791 chebi-645 chebi-1058 chebi-2081 chebi-2771 chebi-70 tiers-131305 chebi-3014
 937. chebi-123 chebi-1368 chebi-3183 chebi-1295 tiers-156410 chebi-1884 chebi-2885 chebi-3097 tiers-121896 chebi-1033 chebi-2855 tiers-121966
 949. chebi-2984 tiers-108673 chebi-3164 chebi-2700 chebi-1203 chebi-1289 chebi-2831 chebi-1026 chebi-1042 chebi-1906 chebi-3233 chebi-1243
 961. chebi-1149 chebi-795 chebi-1818 chebi-3270 tiers-67934 chebi-1518 tiers-5030 chebi-3226 chebi-206 chebi-2535 chebi-2807 chebi-1201
 973. chebi-2117 chebi-2018 chebi-729 chebi-17 chebi-1985 chebi-1501 chebi-1137 chebi-665 tiers-109289 chebi-2943 chebi-454 chebi-2964
 985. chebi-2184 chebi-3262 chebi-1905 chebi-2606 tiers-112188 chebi-1636 tiers-121931 chebi-2084 tiers-148773 chebi-1580 chebi-366 chebi-376
 997. chebi-3112 chebi-1186 chebi-2651 chebi-791 chebi-1876 chebi-3021 chebi-67 chebi-189 chebi-2146 chebi-2390 chebi-1212 chebi-406
1009. chebi-2464 chebi-2898 chebi-2185 chebi-3194 chebi-3236 chebi-804 chebi-3295 chebi-924 chebi-1512 chebi-327 chebi-1379 chebi-1593
1021. chebi-1833 chebi-1442 chebi-2126 chebi-3064 tiers-121943 chebi-1997 chebi-2363 chebi-2254 chebi-1773 chebi-2353 chebi-274 chebi-3182
1033. chebi-353 chebi-598 chebi-1284 chebi-788 tiers-84342 chebi-692 chebi-553 chebi-2810 chebi-2842 chebi-1398 chebi-1677 chebi-3082
1045. chebi-2503 chebi-1169 chebi-2016 chebi-117 chebi-2182 chebi-1121 chebi-753 chebi-1701 chebi-329 chebi-1831 chebi-840 chebi-3006
1057. chebi-2516 chebi-1690 chebi-1330 chebi-2789 chebi-2994 chebi-990 chebi-2944 chebi-1344 chebi-1445 chebi-941 chebi-2881 chebi-3016
1069. chebi-14 chebi-2808 chebi-2270 chebi-3155 chebi-600 chebi-3173 chebi-2386 chebi-2709 chebi-755 chebi-629 chebi-71 chebi-2925
1081. chebi-2157 chebi-679 chebi-1413 chebi-408 chebi-251 chebi-2545 chebi-2315 chebi-365 chebi-445 chebi-1668 chebi-429 tiers-101101
1093. chebi-2568 chebi-1715 chebi-2727 chebi-562 chebi-2644 chebi-955 chebi-99 chebi-1274 chebi-2430 chebi-3105 chebi-2891 chebi-2687
1105. chebi-248 chebi-2289 chebi-3012 chebi-2495 chebi-1135 chebi-2952 chebi-2426 chebi-2934 chebi-2237 chebi-258 chebi-2111 chebi-1623
1117. chebi-473 chebi-2582 chebi-1350 chebi-1687 tiers-77209 chebi-1926 chebi-883 chebi-201 chebi-111 chebi-1562 chebi-156 chebi-379
1129. chebi-2664 chebi-53 chebi-2203 chebi-2680 chebi-3165 chebi-1475 chebi-179 chebi-195 chebi-1363 chebi-1083 chebi-1978 chebi-1896
1141. chebi-756 chebi-1877 chebi-3101 chebi-1162 chebi-1197 chebi-2453 chebi-488 chebi-2863 chebi-1187 chebi-677 chebi-2389 chebi-2356
1153. chebi-520 chebi-2097 chebi-690 chebi-1325 chebi-1109 chebi-1320 chebi-549 chebi-758 chebi-1351 chebi-1466 chebi-3276 chebi-1020
1165. chebi-839 chebi-2940 chebi-510 chebi-332 chebi-3119 chebi-104 chebi-3061 chebi-3251 chebi-529 chebi-1661 chebi-2350 chebi-3196
1177. chebi-777 chebi-2262 chebi-900 chebi-1488 chebi-1916 chebi-2142 chebi-1722 chebi-3071 chebi-1075 chebi-2073 chebi-1838 chebi-2029
1189. chebi-1146 chebi-2774 chebi-435 chebi-1184 chebi-2501 chebi-2333 chebi-737 chebi-390 chebi-2 chebi-1718 chebi-238 chebi-3015
1201. chebi-324 chebi-2630 chebi-2783 chebi-2932 chebi-2176 chebi-987 chebi-2963 chebi-1385 chebi-1543 chebi-1747 chebi-1546 chebi-1716
1213. chebi-1987 chebi-569 chebi-3052 chebi-217 chebi-2488 chebi-1704 chebi-2168 chebi-3253 chebi-3130 chebi-1461 chebi-805 chebi-122
1225. chebi-1354 chebi-2968 chebi-151 chebi-781 chebi-1339 chebi-1082 chebi-1592 chebi-1552 chebi-1767 chebi-2826 chebi-2608 chebi-416
1237. chebi-1005 chebi-1527 chebi-1133 chebi-2352 chebi-348 chebi-1651 chebi-3039 chebi-249 chebi-2425 chebi-2150 chebi-2776 chebi-1601
1249. chebi-1324 chebi-1729 chebi-1176 chebi-1171 chebi-1965 chebi-1989 chebi-2710 chebi-3238 chebi-2841 chebi-1719 chebi-1656 chebi-1950
1261. chebi-203 chebi-3120 chebi-1959 chebi-1100 chebi-1790 chebi-453 chebi-2867 chebi-1238 chebi-861 chebi-168 chebi-2726 chebi-2839
1273. chebi-2534 chebi-2536 chebi-2747 chebi-1358 chebi-593 chebi-2219 chebi-252 chebi-1394 chebi-2923 chebi-509 chebi-567 chebi-819
1285. chebi-225 chebi-2670 chebi-1599 chebi-1062 chebi-2172 chebi-357 chebi-1682 chebi-2702 chebi-482 chebi-2440 chebi-1939 chebi-2611
1297. chebi-1317 chebi-2336 chebi-286 chebi-2958 chebi-401 chebi-2012 chebi-442 chebi-2669 chebi-2584 chebi-3220 chebi-2665 chebi-966
1309. chebi-1045 chebi-1620 chebi-2992 chebi-719 chebi-1936 chebi-3070 chebi-1857 chebi-2451 chebi-3147 chebi-331 chebi-536 chebi-582
1321. chebi-2951 chebi-124 chebi-541 chebi-902 chebi-2127 chebi-427 chebi-2583 chebi-235 chebi-3275 chebi-2797 chebi-2458 chebi-1923
1333. chebi-322 chebi-2056 chebi-1827 chebi-2733 chebi-1066 chebi-1110 chebi-292 chebi-210 chebi-2487 chebi-886 chebi-2279 chebi-1913
1345. chebi-1227 chebi-2835 chebi-1424 chebi-504 chebi-31 chebi-28 chebi-2278 chebi-673 chebi-2281 chebi-2446 chebi-2171 chebi-2836
1357. chebi-887 chebi-1816 chebi-2793 chebi-948 chebi-2349 chebi-2209 chebi-2756 chebi-3110 chebi-2579 chebi-1340 chebi-599 chebi-1630
1369. chebi-2845 chebi-2676 chebi-1792 chebi-2131 chebi-128 chebi-1879 chebi-3084 chebi-2844 chebi-3286 chebi-237 chebi-958 chebi-545
1381. chebi-1111 chebi-2499 chebi-439 chebi-1034 chebi-2580 chebi-1037 chebi-1840 chebi-3149 chebi-47 chebi-2054 chebi-2304 chebi-1313
1393. chebi-363 chebi-2518 chebi-223 chebi-2241 chebi-68 chebi-1249 chebi-735 chebi-2082 chebi-697 chebi-3029 chebi-2418 chebi-2929
1405. chebi-3083 chebi-2398 chebi-1280 chebi-1676 chebi-3098 chebi-3107 chebi-2586 chebi-2259 chebi-2383 chebi-2695 chebi-2321 chebi-1425
1417. chebi-1418 chebi-3168 chebi-3088 chebi-3113 chebi-1044 chebi-904 chebi-1188 chebi-2635 chebi-2633 chebi-803 chebi-2512 chebi-3174
1429. chebi-1602 chebi-1077 chebi-1390 chebi-1548 chebi-3090 chebi-2362 chebi-2199 chebi-2382 chebi-2666 chebi-725 chebi-875 chebi-3287
1441. chebi-3049 chebi-2624 chebi-1665 chebi-2954 chebi-982 chebi-762 chebi-1918 chebi-1536 chebi-2654 chebi-1292 chebi-1582 chebi-2492
1453. chebi-1262 chebi-2719 chebi-1565 chebi-2947 chebi-49 chebi-1632 chebi-1637 chebi-2988 chebi-1830 chebi-1107 chebi-2557 chebi-1870
1465. chebi-347 chebi-1247 chebi-838 chebi-183 chebi-3221 chebi-2013 chebi-1854 chebi-766 chebi-1316 chebi-1018 chebi-1629 chebi-2149
1477. chebi-2062 chebi-1097 chebi-3265 chebi-609 chebi-56 chebi-1549 chebi-2374 chebi-1720 chebi-400 chebi-349 chebi-2396 chebi-3004
1489. chebi-154 chebi-101 chebi-649 chebi-559 chebi-1763 chebi-3272 chebi-2866 chebi-1617 chebi-76 chebi-979 chebi-1144 chebi-2411
1501. chebi-3244 chebi-1399 chebi-2813 chebi-1738 chebi-3065 chebi-1909 chebi-650 chebi-953 chebi-1645 chebi-1476 chebi-623 chebi-54
1513. chebi-1772 chebi-103 chebi-2877 chebi-3125 chebi-771 chebi-2348 chebi-1977 chebi-388 chebi-2751 chebi-352 chebi-1429 chebi-1544
1525. chebi-1446 chebi-2969

---

## 附录 B：原始材料索引

| 文件 | 内容 |
|---|---|
| `tmp/bench_full.json` | 全量 benchmark 原始 stdout（含进度行） |
| `tmp/fails.json` | 1526 条错例的结构化记录 |
| `tmp/baseline.json` | 基线准确率快照（改动前后对比用） |
| `tmp/chunks/chunk_0N.md` / `.json` | 按复杂度切分后的 6 块错例清单 |
| `tmp/analysis/chunk_0N.md` | 6 个 subagent 的分块根因报告（含每个模式的完整 id、代表例、IUPAC 依据） |
| `tmp/analysis/theme_ids.md` / `.json` | 13 个全局主题的 id 清单（本报告 §3 的来源） |
| `tmp/analysis/appendix_ids.md` | 全部错例的 id / HAC / SMILES / GOLD / PRED 明细表 |
| `tmp/analysis/probe_*.py`、`c0N_probe*.py` | 各 subagent 的 RDKit 核实脚本 |

---

## 7. 五个主题的修复结果（worktree 并行修复，未合并）

用户在 §6 批准了主题 **1 + 2 + 5 + 4 + 8**。5 个 subagent 各在一个独立 git worktree 里修复，
只改 `src/namepredict/**`，未改任何测试。

### 7.1 分主题结果

| # | 主题 | 分支 | 改动文件 | +行 | −行 | acc_en 前→后 | 全库修复 | 回归 | 主题内通过 | 合并后通过 |
|---|---|---|---:|---:|---:|---|---:|---:|---|---|
| 1 | G1 碳酰/含氮衍生物保留名 | `wt/t1-carbonyl-retained` | 5 | 199 | 13 | 61.95 → **62.95%** | 34 | 0 | 33/132 | **40/132** |
| 2 | G3 含氧酸类名后缀与母体优先级 | `wt/t3-oxoacid-parent` | 2 | 32 | 2 | 61.95 → **62.60%** | 26 | 0 | 26/219 | 26/219 |
| 5 | G5 指示氢与氢化程度 | `wt/t5-indicated-h` | 3 | 159 | 10 | 61.95 → **62.45%** | 22 | 2 | 20/81 | 20/81 |
| 4 | G6 取代基引用顺序 | `wt/t6-prefix-order` | 4 | 35 | 11 | 61.95 → **62.50%** | 22 | 0 | 18/72 | 18/72 |
| 8 | G2 保留前缀补表 | `wt/t2-retained-prefix` | 6 | 98 | 13 | 61.95 → **62.50%** | 19 | 0 | 17/43 | 17/43 |

「全库修复 / 回归」由协调者用同一条 benchmark 命令对基线与各分支各跑一次、逐 id 比对得出，
不采信 subagent 自报。所有分支的改动文件 100% 落在 `src/namepredict/` 内。

### 7.2 合并预演结果（`tmp/mergecheck` 分支，**未合入 main**）

| | 基线 main | 合并后 |
|---|---|---|
| acc_en | 61.95%（2463/3976） | **65.32%**（2597/3976，+134 条） |
| acc_zh | 77.41%（579/748） | **79.55%**（595/748，+16 条） |
| 错例总数 | 1526 | 1401 |

- 修复 **127** 条，回归 **2** 条，净 **+125**。
- 5 个分支单独修复的并集是 123 条，合并后 127 条 —— 有 5 条是分支间的协同
  （`chebi-1549 chebi-2374 chebi-2920 chebi-400 chebi-56`：G1 的保留名 + G2 的前缀补表叠加才转正），
  1 条丢失（`tiers-59838`）。
- 合并时有 2 处冲突：`layer5/assembler.py` 与 `layer5/assembler_prefixes.py` 的 import 行
  （取并集），以及 `_mult_rows` 排序键 —— 保留 G6 的 `alpha_order_key`（P-14.5 修复）
  并叠加 G1 的脲 N 位次排序。
- 合并结果上 pytest：**1408 passed / 1 failed**，失败项 `test_amide.py::test_guard_no_false_normalization[Oc1ccccn1]`
  在基线 HEAD 上复跑同样失败，是既有缺陷。

### 7.3 未修复部分与原因

| 子类 | 条数 | 原因 |
|---|---:|---|
| 多聚磷酸 P–O–P 母体 | ~82 | `oxoacid` 的 SMARTS 要求 P 恰 3 个 O 邻居，焦/三/多磷酸不匹配；且 GOLD 有三种互不兼容的酸酐写法 |
| 硫酯 `-thioate` | ~48 | src 中完全没有硫代酸酯官能团类，需整条新链 |
| `-thione` / `thioxo` | 5 | GOLD 判据自相矛盾：饱和的 `1,3-thiazolidine-2-thione` 用 thione，同样饱和的 `1,3-diazinane-4,6-dione` 却用 `sulfanylidene` |
| 脒 / 胍 / 硫代酰胺 / 羟胺 / 肼 | ~10 | 各需新增 FG 检测类 + 新母体 kind，非数据表可解 |
| 5 元内酯环母体（`2H-furan-5-one` vs `2,5-dihydrofuran-2-one`） | ~10 | 属母体/后缀选择而非指示氢 |
| 金标体例冲突类 | 若干 | 见下 |

### 7.4 修复过程中新发现的 4 处金标问题

| 位置 | 问题 |
|---|---|
| `chebi-2982` / `tiers-127` | GOLD 省略了 NH 指示氢（`pyridin-2-one`），但同族 `chebi-2673` 又写 `1H-pyridin-2-one` —— 内部约定冲突，无法同时满足，是本轮唯一的 2 条回归 |
| `chebi-2293` / `chebi-2878` | GOLD 把 `2-oxo`/`2-hydroxy` 写在 `3-(5-benzyloxy-…)` 之前，按 P-14.5 应为 benzyloxy < oxo。PRED 改对后反而与 GOLD 不同 |
| N-取代碳酰围栏 | 三种写法并存：`2-methoxyethylcarbamoyl`（不围栏）/ `(4-methylcyclohexyl)carbamoyl`（围栏）/ `(5-carboxy-2-ethylpentyloxy)carbonyl`（围栏） |
| 磺酰胺 N,N-二取代 | `[2-(4-ethylanilino)-2-oxoethyl]-phenylsulfamoyl` 与 `cyclohexyl(ethyl)sulfamoyl` 两种互斥体例 |

**另修正本报告 §2 / §3 的一处误报**：底层 chunk_06 报告提出的「33 处同一环位次被两个取代基占用」，
经全库插桩核实为**误报** —— 改动前后这类名字数量完全一致（1035 处 / 855 个唯一名），
且绝大多数是合法的**孪位双取代**（`2-(dimethylamino)-2-thiophen-2-yl` 这类季碳形式）。
本报告 §6 中把它列为 G6 的「顺带关联 bug」是不成立的。
