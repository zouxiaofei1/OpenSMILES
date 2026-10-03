# 中英双语命名 (Bilingual Naming)

> **概念层级:** 跨层核心概念 | **涉及层次:** Layer0–Layer5 | **核心数据结构:** `(en, zh)` 元组 | `_Chain` | `_KIND_TABLE` | `NameResult(en, zh)`

---

## 概述

管线同步产出英文与中文两套 IUPAC 名称：语言无关的 `numbered` 字典是唯一数据源，各层分别渲染两种词形，不经「先生成一语再翻译」的中间阶段。双语函数一律返回 `tuple[str, str] | None`——任一语言无法生成即整体返回 `None`，由派发链降级。

## 词表分层

双语文本来源分三层，职责互不重叠：跨层词表与判据在 `constants.py`、碳数词干与盐名派生在 `layer5/stems.py`、kind → 词尾分派在 `chain_engine._KIND_TABLE`。

`constants.py` 存放的词表（中英成对，或本身即双语元组）：

| 词表 | 语义 |
|---|---|
| `MULT_EN` / `MULT_ZH` | 1–9999 倍数数量词，由 `en_num_term`/`zh_numeral` 派生（含百位/千位） |
| `BIS_EN` / `BIS_ZH` | P-16.3.2 复合前缀倍增：`bis`/`tris`/`tetrakis`… ↔ 双/三/四…；`BIS_EN_SET` 供 `mult in …` 常数判据 |
| `BRIDGE_SUFFIX_EN` / `BRIDGE_SUFFIX_ZH` | O/S/N 桥后缀，中英同序：`oxy`/`sulfanyl`/`amino` ↔ 氧基/硫基/氨基 |
| `DIATOMIC_BRIDGE_YL` / `BRIDGE_DIATOMIC_ZH` | 双原子桥合一保留前缀 `diazenyl`/`disulfanyl` 及中文名尾 |
| `BRIDGE_FUSION_YL` | `(中心氢化物词干, 前端名尾)` → `(中文名尾候选, 合一 en/zh)`；含三硫链 `disulfanyl→trisulfanyl` 一行 |
| `BRIDGE_SPLIT_SUFFIX_EN` / `BRIDGE_SPLIT_SUFFIX_ZH` | 可拆桥后缀全集（含双原子桥） |
| `MONONUCLEAR_YL` | `free_en` → `(free_zh, 去氢 yl_en, 组装名中文尾)` |
| `OXO_CENTER_KINDS` / `ESTER_O_SIDE_KINDS` | 中心原子自任母体 / O-侧臂母体 kind 集合 |
| `EXO_RING_SUF` | 环外主基后缀 `(singular, plural\|None)`，消费方 `_exo_ring_spec` |
| `ALKOXY_YLOXY_EN` / `ALKOXY_YLOXY_ZH` | O 锚点 `-yloxy` 收拢为保留烷氧基 |
| `CHAIN_RETAINED` | C1/C2 开链保留名（酸/酰基/醛/酰胺/腈/酯） |
| `AMIDO_RETAINED` / `AZANE_PAREN_SUF` | N-酰基保留式与 azane 内层组括起白名单 |
| `ISO_NUCLIDE_PROP` / `ISO_H_PROPS` | 非氢核素属性名与氘/氚三元组（质量数、重原子属性名、en/zh 词干） |
| `HW_PREFIX_EN` / `HW_PREFIX_ZH` | Hantzsch-Widman 杂单环 'a' 前缀（表 2.4 / 表 3-3 长式） |
| `HW_UNSAT_SIX` / `HW_SAT_SIX` / `HW_UNSAT_TAIL` / `HW_SAT_TAIL` | 六元 A/B/C 组词干与七至十元词尾 |
| `HW_ZH_SHORT` / `HW_ZH_RING` / `HW_MAX_VALENCE` | 中文短式（噁唑/噻唑）、中文环前缀、标准键数参照 |
| `PARENT_HYDRIDE_STEMS` / `HYDRIDE_YL_FORMS` | 原子序数 → 母体氢化物词干（P-21 表 2.1）；去氢 `-yl`/`-ylidene` 六元组 |
| `HYDRIDE_MULT_ZH` / `hydride_chain_stem` | 多核氢化物中文天干倍数（甲乙丙…）与裸词干生成器 |
| `SIMPLE_MOLECULES` / `ELEMENT_METAL_NAMES` | L0 单质与简单分子名（规范 SMILES 或元素符号 → `(en, zh)`） |
| `STANDARD_BONDING_NUMBERS` | 原子序数 → 标准键数（判别 λ 记号，P-14.1.2） |
| `MONONUCLEAR_CATION_SYSTEMATIC` / `CATION_RETAINED_NAMES` / `cation_parent_names` | 单核母体阳离子系统名与保留名（en 取系统、zh 取保留） |
| `RETAINED_DEHYDRO_YL` / `retained_dehydro_yl` | 保留名去氢前缀改写（adamantan-2-yl → 2-adamantyl） |
| `ZH_DIGITS` / `HS_NUMBER` / `zh_bridge_root` | 中文数字词、天干与去「基」函数（酰基「羰基」尾保留） |

`layer5/stems.py` 在词表之上派生：`_en_stem`/`alkane_en`/`alkane_zh`/`zh_num`/`zh_stem` 产碳数词干，`_metal_en_prefix` 取 `MULT_EN` 加金属名、`_metal_zh_suffix` 取 `MULT_ZH` 加中文金属名、`join_anion_names` 把酸转 `-ate`/补「根」、`join_metal_salt_names` 是盐后缀总入口；英文侧无对应 `zh_stem`，C11+ 词干由 `en_num_term(n)` 去尾 `a` 得；`stem_forms` 由「主体名 + 词尾」统一切出（完整名, 裸词干）两形态，供桥环/螺环/生成式共用。

`MONONUCLEAR_HYDRIDES` 每行五元组 `(元素, 中文名, 零价去氢名对, 桥后缀对, 组装名中文尾)`，下游四张表全由它派生、中英永远成对：`MONONUCLEAR_ZERO_YL`（`(free_en, free_zh)` → 零价去氢名对）、`MONONUCLEAR_BRIDGE`（`(free_en, free_zh)` → 桥后缀对，仅 O/N/S 三行）、`MONONUCLEAR_YL`（`free_en` → `(free_zh, 去氢 yl_en, 组装名中文尾)`）、`PHOSPHORYL_STEMS`（P 酰基词干元组，避免 P 被当碳中心）；元素默认 free 名另由 `MONONUCLEAR_BY_ELEMENT` 给出。

键取 `(en, zh)` 二元组而非单独英文名：单核氢化物是同一实体的两个词形，去氢与桥接须同时对两者成立。消费点是 `assembler._mononuclear_radical_names` 与 `free_to_yl`——后者由 `_fg_prefix` 成对约束 `_mononuclear_en`/`_mononuclear_zh`，任一语言未命中即整体返回 `None`。

`PHOSPHORUS_STEM_BY_OXO`、`SULFUR_STEM_BY_OXO`、`NITROGEN_STEM_BY_FREE_DOUBLE` 按 `=O` 数或自由价键级选定锚点母体词干（`phosphanyl`/`phosphoryl`、`sulfane`/`sulfinyl`/`sulfonyl`、`azane`/`imine`），缺一则锚点与饱和类似物同形。

母体氢化物族以 `PARENT_HYDRIDE_STEMS`（原子序数 → 词干）为唯一来源，`HYDRIDE_YL_FORMS` 给出同一行实体的 `-yl`/`-ylidene` 去氢对，`hydride_chain_stem` 用 `en_num_term` 与 `HYDRIDE_MULT_ZH`（天干倍数）拼多核裸词干；`STANDARD_BONDING_NUMBERS` 是 λ 记号的判据来源。单核阳离子侧 `MONONUCLEAR_CATION_SYSTEMATIC` 给系统名、`CATION_RETAINED_NAMES` 给第 15/16/17 族保留名，`cation_parent_names` 合成 `(en, zh)`，`CATION_YL_STEMS` 为这些阳离子派生的去氢前缀（倍数用 `bis(azaniumyl)`）。L0 简单分子与盐反离子名分列 `SIMPLE_MOLECULES`/`ELEMENT_METAL_NAMES` 与 `layer0/salt.py` 的 `POLY_CATION`/`POLY_ANION`。

## kind → 词尾分派

`chain_engine._KIND_TABLE` 每个 entry 是冻结数据类 `_Chain` 实例，`_chain_names` 是唯一渲染入口，`assembler._names_for` 是 L5 派发枢纽。

| kind | `en_suf` | `zh_suf` | `coda` |
|---|---|---|---|
| `alcohol` | `ol` | 醇 | `an` |
| `ketone` | `one` | 酮 | `an` |
| `thione` | `thione` | 硫酮 | `ane` |
| `alkane` | `ane` | 烷 | `""` |
| `acid` | `oic acid` | 酸 | `an` |
| `sulfonic` | `sulfonic acid` | 磺酸 | `ane` |
| `sulfonate` | `sulfonate` | 磺酸 | `ane` |
| `sulfonamide` | `sulfonamide` | 磺酰胺 | `ane` |
| `sulfonyl_chloride` | `sulfonyl chloride` | 磺酰氯 | `ane` |
| `ester` | `oate` | 酸 | `an` |
| `phosphate` | `phosphate` | 磷酸 | `""` |
| `phosphonate` | `phosphonic acid` | 膦酸 | `""` |
| `sulfate` | `sulfate` | 硫酸 | `""` |
| `boronic` | `boronic acid` | 硼酸 | `""` |
| `acyl` | `oyl` | 酰基 | `an` |
| `thiol` | `thiol` | 硫醇 | `ane` |
| `amine` | `amine` | 胺 | `an` |
| `aldehyde` | `al` | 醛 | `an` |
| `nitrile` | `enitrile` | 腈 | `an` |
| `amide` | `amide` | 酰胺 | `an` |
| `acyl_halide` | `oyl {卤}` | 酰{卤} | `an` |
| `radical` | `yl` | 基 | `an` |

要点：

- `acyl_halide` 默认 entry 为 `_ACYL_HALIDE_BY_HAL[Cl]`，`_names_for` 按 `parent["hal_z"]` 覆盖为对应卤素 spec；C1/C2 保留名、苯甲酰与草酰二卤由 `variant` 提供。
- C1/C2 保留名由 `CHAIN_RETAINED` 经 `_retained_plain` 挂进 `variant`；苯环单取代由 `_BENZENE_RETAINED` 在 `_names_for` 按 `scaffold_id` 覆盖，硫代羧酸 S-酯、硫代酰胺、脒、酰肼另有苯环专属保留名（`benzenecarbothioate`/`benzenecarbothioamide`/`benzenecarboximidamide`/`benzohydrazide`）。
- `ester` 在 `parent["thio_side"]` 成立时由 `_names_for` 改写为 `coda="ane"`、`en_suf="thioate"`、`zh_suf="硫"`、`variant=None`，烯/炔基座改 `enethioate`/`烯硫` 与 `ynethioate`/`炔硫`。
- `amide` 在 `parent["amide_z"]` 为 `S` 时改写为 `thioamide`/`硫代酰胺`（coda `ane`、`variant=None`），为 `N` 时改写为 `imidamide`/`亚氨酰胺`（P-43 类 16/17）；`parent["hydrazide_n_idx"]` 存在时经 `hydrazide_chain_spec` 换 `hydrazide`/`酰肼` 尾（coda `ane`，C1/C2 保留名经 `variant` 注入）。
- 多 FG（`mult > 1`）由 `_generated_mult_fields` 把 `MULT_EN`/`MULT_ZH` 注入 `en_suf`/`zh_suf` 并强制 `coda="ane"`、`need=mult`、`no_loc="none"`，故「二醇/二酸」与单 FG 共用同一渲染路径。
- 段式烯/炔段由 `_unsat_seg` 按后缀首字母推导：`en_suf` 以元音（含 `y`）开头取 `en`/`yn`，否则留末端 `e`（`carboxylic acid`、`thiol` 取 `ene`/`yne`，P-16.7.1(a)）；显式 `ene_seg`/`yne_seg` 优先。`radical` 的 `plain_fn=_radical_plain` 在自由价为 1 位时由词干拼 `-yl`/基；母体带 `radical_anchor_element` 时 `_names_for` 转向 `_mononuclear_radical_names`，不查 `_KIND_TABLE`。

## 含氧酸词尾表

`chain_engine._OXO_TAIL` 以 `(oxo_kind, 中心酸式氢数)` 为键出双语功能母体词尾；查表函数 `_oxoacid_tail` 读 `parent["oxo_kind"]`/`parent["n_oh"]`，表外返回 `None`。该函数经 `_Chain.plain_hook` 挂在 `phosphate`/`phosphonate`/`sulfate`/`boronic` 上，`_chain_names` 先走钩子分支再查碳数词表。S 中心两臂皆酸式氧或 O-臂者（含 S–O–S 多硫酸链）判为 `sulfate`，词尾同样由本表给出。十三项：

- `phosphate` 氢数 3/2/1/0 → `phosphoric acid`/磷酸、`dihydrogen phosphate`/磷酸二氢、`hydrogen phosphate`/磷酸氢、`phosphate`/磷酸；
- `phosphonate` 2/1/0 → `phosphonic acid`/膦酸、`hydrogen phosphonate`/膦酸氢、`phosphonate`/膦酸；
- `sulfate` 2/1/0 → `sulfuric acid`/硫酸、`hydrogen sulfate`/硫酸氢、`sulfate`/硫酸；
- `boronic` 2/1/0 → `boronic acid`/硼酸、`hydrogen boronate`/硼酸氢、`boronate`/硼酸。

## 桥环与 C1 保留名

`bridged_namer` 由 L4 裁决后的 `bridged_node` 出名，中英同构：环数词头 `ring_count_prefix` 取 2 环 `bi`/双、≥3 环 `MULT_EN`/`MULT_ZH`，非 `dicyclo`；描述符串的长度数字以 `.` 连接、次级桥位次对直接续在长度数字后；骨架置换前缀按 `P145_SENIOR` 序前置（`4-thia-1-aza` / `4-硫杂-1-氮杂`），数量词尾 `a` 仅在后接元素名以 `a` 开头时省略。`skeleton_replacement.prefix_from_chain` 另可对键数非标准的骨架杂原子在位次后紧贴 `λn`（P-23.6.1，`lambda_ok=False` 时用于稠合组分名不标 λ）。

`assembler._ensure_parent_stem` 按螺环 → 桥环 → 稠环 → 生成式注入词干，同一母体名有两种形态：完整名 `stem_en`/`stem_zh`（FG 分支去末端 `e` 接后缀）与裸词干 `stem_bare_en`/`stem_bare_zh`（链式词干引擎拼 `ane`/`烷` 或不饱和段）。`_c1_retained`/`_c1_amino` 先于系统名出 C1 保留名族：`urea`/脲、`thiourea`/硫脲、`guanidine`/胍、`carbonate`/碳酸、`carbamoyl`/氨基甲酰基、`carbamic acid`/氨基甲酸、`carbamate`/氨基甲酸酯、`carbamothioyl`/氨基硫代羰基。N 位次口径：脲/硫脲/胍（`N_LOCANT_KINDS`）须显式写出，N-取代氨基甲酸（`carbamic_acid`）省略；两者都把 subst `kind` 改写为 `side`，改走数字位次通道。

## 语序：酯 / 硫代酯 / 盐 / 含氧酸中心母体

### 酯

`assembler.join_ester_name` 让 O-侧臂作前缀、酸侧作主体：英文 `{臂} {母体}`，中文 `{母体}{臂}酯`。臂由 `_o_side_arms` 从 `substituents[]` 筛 `o_side` 标记者，交 `_join_o_side_arms(group=False, …)` 拼接：同名臂用 `MULT_EN`/`MULT_ZH` 倍增，任一臂带 `paren` 时英文用 `BIS_EN` 的 `bis(…)`，异名臂平铺。异名臂且能取到母体 O-位次时（`_located_ester_arms`）改走 `{loc}-O-{臂}` 引用序（P-16.6.2 杂原子位次 + P-14.5 字母数字序）。中文臂词形由 `_zh_alkoxy_part` 渲染，英文侧不围栏。

### 硫代酯

`parent["thio_side"]` 成立即硫代羧酸 S-酯（P-65.6.3.3.7.1），语序为 `S-{臂} {母体}` / `S-{臂}{母体}酯`——`S-` 在最前，酯词尾接母体名后。该分支在 `join_ester_name` 内单独构造，S 侧臂中英围栏走 `_fenced_arm_en` 与 `_fenced_arm_zh`；词尾由 `_names_for` 改写为 `thioate`/`硫`，苯环母体取 `_benzene_retained("benzenecarbothioate", "苯硫代甲酸")`，环外母体在 `_exo_ring_spec` 取 `("carbothioate", "硫代甲酸")`。

### 盐

L0 `layer0/salt.py` 先解离反离子：单原子金属（`METAL_ION_EN`/`METAL_ION_ZH`）、卤离子（`HALIDE_EN`/`HALIDE_ZH`）、氢卤酸（`HALIDE_HX_EN`/`HALIDE_HX_ZH`）之外，`POLY_CATION`/`POLY_ANION` 收多原子反离子（`azanium`、`hydroxide`、`nitrate`、`perchlorate`、`azide`、`hexafluorophosphate`），阳离子并入金属通道、阴离子并入卤离子通道；`_from_frags` 输出 `n_metal`/`n_org` 等元数据。

`stems.join_metal_salt_names` 表达语序反转：英文金属名作前缀、仅当 `en.endswith("ate")` 时拼接，数量词取 `MULT_EN`（`_metal_en_prefix`）；中文以金属名替换末尾「根」（`zh[:-1] + suf`），数量词取 `MULT_ZH`（`_metal_zh_suffix`），得「二钠」式。有机阴离子份数 `n_org > 1` 时英文改 `{metal} {BIS_EN}({en})`、中文 `{BIS_ZH}({zh}){metal}`（`calcium bis(…acetate)`）；`azanide` 母体按 `{metal} …azanide` 前置、中文「金属盐」后置。盐酸盐走 `acid_salt`/`acid_salt_zh`：英文空格相接，中文直接相接。盐后缀由 `namer._apply_salt_suffix` 在组装之外施加，母体 `parent_kind` 属 `OXO_CENTER_KINDS` 时跳过（盐已在整名内组装）。

### 含氧酸中心母体

`kind in OXO_CENTER_KINDS`（`phosphate`/`phosphonate`/`sulfate`/`boronic`）时 `join_kind_name` 交给 `assembler.join_oxoacid_name`，整名按「取代前缀 + O-侧臂 + 功能母体词尾 + 金属盐」构造，盐元取自 `parent["salt_meta"]`：

- 无臂有金属：酸式盐，`{metal} {tail}` / `{tail}{metal_zh}`，不出现「酯」字；
- 无臂且 `n_om > 0`：游离酸式根，中文补「根」；
- 有臂有金属：`{metal} {alk} {tail}` / `{tail}{alk}酯 {metal}盐`；
- 有臂无金属：`{alk} {tail}` / `{tail}{alk}酯`。

臂围栏分两路：`phosphate` 走 `_fenced_arms_phosphate`（缩合磷酸 P-O-P 二酯、全部臂名带前导立体描述符时经 `_all_arms_stereo_lead` 整体围栏，否则按 `oxo_arm_fence` 或「多个位次段」判定），其余中心母体走 `_fenced_arms`。

### 官能母体 (硫)氰酸酯

腈碳唯一重邻居为 S/O 时（`parent.kind == "nitrile"` 且单碳单取代、取代基名尾即桥元素），`_chalcogen_nitrile_alias` 把分子式名改写为官能母体名：`_CHALCOGEN_NITRILE_TAIL` 以桥尾为键给出 `sulfanyl → ("thiocyanate", "硫氰酸", "硫基")`、`oxy → ("cyanate", "氰酸", "氧基")`，得英文 `{R} {tail}`、中文 `{tail}{R}酯`（`methyl thiocyanate` / 硫氰酸甲酯，P-65.6.3.1）。

## 取代基词形的成对产出

保留式取代基以「英文取保留式、中文取对应基名」成对产出（`tools/anchored_table.py` 的 `RetainedSubstituent` 登记锚定 SMARTS 与 `paren`）：

| 结构 | EN | ZH | 依据 |
|---|---|---|---|
| N-酰基单取代（方法 1） | `acetamido` / `formamido` / `benzamido` | 乙酰氨基 / 甲酰胺基 / 苯甲酰胺基 | P-66.1.1.4.3，`AMIDO_RETAINED` |
| 长链/烯酰/杂环羰酰（方法 2） | `acetylamino` 等 | …酰氨基 | 不入表，走 `free_to_yl` |
| 质子化伯胺 | `azaniumyl` | 铵基 | P-62.4.1，锚定保留叶 |
| 苯胺去氢 | `anilino` / `4-chloroanilino` | 苯胺基 / 4-氯苯胺基 | P-62.1.1，`_mononuclear_en` / `_mononuclear_zh` |
| N-取代磺酰胺基 | `butyl(methyl)sulfamoyl` | …磺酰基 | P-66.1.1.4.2，与 `sulfamoyl` 融合 |
| O 锚点 `-yloxy` | `ethyloxy` → `ethoxy` | 乙基氧基 → 乙氧基 | P-66.5.2.1.2，`ALKOXY_YLOXY_EN` / `ALKOXY_YLOXY_ZH` |
| 氰氧基 / 硫氰基 | `cyanato` / `thiocyanato` | 氰氧基 / 硫氰基 | P-67.1.4.2，锚 `*OC#N` / `*SC#N` |
| 重氮鎓 | `diazonio` | 重氮鎓基 | P-65.3，锚 `*[N+]#N` |
| 甲硒基 / 三甲基甲硅烷基 | `methylselanyl` / `trimethylsilyl` | 甲硒基 / 三甲基甲硅烷基 | P-63.2 / P-67.1.4.2 |
| 硼酸 / 二氟硼烷基 / 二苯基硼烷基 | `borono` / `difluoroboranyl` / `diphenylboranyl` | 硼酸 / 二氟硼烷基 / 二苯基硼烷基 | P-68，`*B(O)O` 等锚 |

成对性靠「英文词干承担判据、中文侧同名同形」：`AMIDO_RETAINED_EN` 供 L3 判免括号（两侧均为「乙酰氨基」）；`AZANE_PAREN_SUF` 决定方法 2 中哪些尾缀须把内层组整体括起再加 `amino`；`free_to_yl` 第三位返回值 `need_paren` 只由英文名形判定（`methylamino`、带环取代基的 `anilino`），中英共用该布尔。

## 中文侧特有处理

- **词干路径**：`_chain_zh_base` 用 `zh_stem` 剥去 `alkane_zh(n)` 的「烷」；`_chain_stem_pair` 在 `unsat_polyol` 且 `zh_full`/`cyclic` 时改取完整中文名。
- **S 链后缀保「烷」**：`_ZH_CODA_KINDS`（`sulfonic`/`sulfonate`/`sulfonamide`/`sulfonyl_chloride`）在 FG 分支取 `alkane_zh(n)`，保留「烷」（丙烷-1-磺酰氯）。
- **酯 O-侧基名**（`_zh_alkoxy_part`）：简单烃基省「基」（乙酸乙酯）；带位次者保留「基」（乙酸丙-2-基酯）；其余加围栏，前端已含圆/方括号时升方括号。
- **含氧酸 O-侧臂**（`_oxoacid_arm_zh`）：简单基去「基」；整串为中文数字（`ZH_DIGITS`+「十」）时补「烷基」；含 `-` 或 `(` 开头的复合名原样保留。
- **「基」的回补**（`_zh_front_ji`）：中文复合前端的「基」被上游切掉时在拆分点补回；酰基前端走融合式，或 `_front_needs_enclosure` 为假时不补。桥后缀前的中文烃基名去尾「基」由 `zh_bridge_root` 承担，以「羰基」结尾的酰基名（羰基/甲氧羰基）保尾不削（金标作羰基氨基）。
- **指示氢前缀的取代式落位**：`_with_indicated_h` 把 L4 定好的位次（`_indicated_h_prefix`，同位次去重）落到词干前并取代词干自带前缀，以避 `1H-7H-` 双写；`_nh_only`/`_ring_n_free` 决定未取代 NH 环是否补指示氢；`_stem_prefix_stale` 在词干自带位次已成羰基碳时丢弃它；`_ensure_fused_stem` 为未注册稠环注入词干时同样带前缀。
- **环阳离子与保留名**：`join_ring_cation_suffix` 英文缀 `-{位次}-ium`，中文沿用母体名（`_zh_ring_cation` 插「-{位次}-正离子」，但鎓类保留名沿用「鎓」）；`zh_1h_parent` 在环被取代时给中文保留名补 `1H-`。

## 围栏与排序

- `_enclose(s)`：已带完整方括号者不二次围栏；含圆括号（排除核素描述符 `(14C)`）则升方括号，否则加圆括号；双语共用的围栏原语。
- `oxo_arm_fence(name, sub, mol)` 判含氧酸中心母体臂名是否整体括起：前导 `(` 真（仅立体描述符）、前导 `[` 假（已是完整围栏）、自带两个以上位次段真、带 `paren` 的复合臂挂环上真。
- `_arm_fence_needed(name)`：臂名内含括号/方括号却无前导围栏（前导方括号者假），或前导立体描述符、酰基尾（`oyl`/`carbonyl`）→ 须整体围栏。
- `_fenced_arm_en`/`_fenced_arm_zh`：酯 O-侧臂围栏，前导位次或自带括号时整体括起。
- `_flat_bridge_stem`：简单前端（烃基/桥前端/`SIMPLE_ALKOXY_NO_PAREN`）+ `R-imino` 融合式（`methylimino`/`hydroxyimino`）整体平铺不围栏（P-66.4.1.2.1）；`_BRIDGE_SPLIT_EN/ZH` 在既有桥后缀外另收 `imino`/`亚氨基` 供拆分。
- `_bridge_body(base, suf, merge)`：平铺式 `_enclose(base) + suf`，桥后缀留括号外；双原子桥整段再括一层；`merge=True` 时前端为「方括号组 + 直链 -yl」者整段同括，否则外包方括号使括界自闭合。
- `_sbridge_flat_stem`：磺酰/亚磺酰桥 + 无取代直链 `-yl` 前端，英文侧平铺不围栏。
- **倍数前缀的围栏**：`mult in BIS_EN_SET`（`bis`/`tris`…）时操作数须整体围栏（`bis(carboxymethyl)`），桥后缀留围栏外；N-型取代基混编且带位次者经 `_n_prefix` 整体加括号。
- **中英同步判据**：`_split_bridge_suffix_zh` 的唯一判据是英文 stem 是否拆分（先调 `_split_bridge_suffix`），故中英围栏必然同形；`bridge_self_enclosed` 标记经 `meta` 回传，使 L5 对已自含围栏的桥前端跳过二次加括号。
- **排序**：`tools.re.alpha_order_key` 是 P-14.5 字母数字序的唯一实现，先比非斜体字母序列、再比首字母前位次；中英共用同一键（`_mult_rows`/`_build_prefix`/`_join_o_side_arms`/`_phosphoryl_sub_names` 均以它作 sort key），故取代基两语排列一致。`alkyl_alpha_key` 是其循环剥离步骤；`_alpha_key` 是 `assembler` 内只留字母的简化键，用于 N-芳基-N-某基胺的 `anilino` 引用序判定。N-型取代基的撇号按引用序查 `_n_prime_map`（同 N 用 `N,N-`、跨不同 N 用 `N,N'-`），酰肼两端由 `_hydrazide_primes` 固定（羰基侧 N、肼远端 N′）。`re.py` 的 `normalize_en`/`normalize_zh`/`nospace` 是判分口径，命名管线不调用。

## 递归命名：双语缓存与宿主根传播

`tools/common_names.CommonNameCache` 以 SMILES 为键存 `NameResult`（含 `en`/`zh`），有 `get`/`put`（容量满抛 `ValueError`）/`clear`。`SMILESNNamer.__init__` 默认构造 20000 条缓存；`name` 先查缓存，未命中才走 `_pipeline` 并 `_cache_put`（容量满的 `ValueError` 静默忽略），写前经 `_canonical_result` 用 `CanonicalRankAtoms` 归一 `meta.parent_chain`。

`root_ctx` 是 `(root_mol, to_root)`——根分子与「当前块原子 → 根原子」映射。`_name_mol` 在 `root_ctx is None`（顶层整分子）时取解盐后的 `organic` 与恒等映射，并另开 2000 条容量的 `run_cache`，使内部 `*` 片段名只在本次运行内共享；非空时沿用上游传入的根分子、映射与缓存。该上下文写入 `info["root_ctx"]`。

L3 的 `as_substituent._radical_yl_from_sub` 解包 `root_ctx`、构造锚定块的 `anchored_to_root` 映射后递归调用 `_name_mol(..., root_ctx=...)`；`_fix_rs_with_real` 用原始根分子 CIP 校正被 `*` 顶替后翻转的异头位 R/S，fresh 与缓存命中均按当前根分子校正一次。桥前端自含围栏的标记同样经 `meta` 回传，保证中英围栏形态一致。

相关页面：[[architecture/overview]]、[[architecture/layer1-analyzer]]、[[architecture/layer5-name-assembly]]、[[concepts/functional-group-priority]]。
