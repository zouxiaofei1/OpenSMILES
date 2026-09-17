# 中英双语命名 (Bilingual Naming)

> **概念层级:** 跨层核心概念 | **涉及层次:** Layer0–Layer5 | **核心数据结构:** `(en, zh)` 元组 | `_Chain` | `_KIND_TABLE` | `NameResult(en, zh)`

---

## 概述

管线从起点到终点同步产出英文与中文两套 IUPAC 名称：语言无关的 `numbered` 字典是唯一数据源，每层各自渲染两种词形，不经"先生成一语再翻译"的中间阶段。所有双语函数的约定是返回 `tuple[str, str] | None`——任一语言无法生成即整体返回 `None`，由派发链降级到下一策略。

## 词表分层

双语的文本来源切成三层，职责互不重叠：**跨层词表与判据集中在 `constants.py`**、**碳数词干与盐名派生在 `layer5/stems.py`**、**kind → 词尾分派在 `chain_engine._KIND_TABLE`**。

`constants.py` 存放的词表（中英成对，或本身即双语元组）：

| 词表 | 语义 |
|---|---|
| `MULT_EN` / `MULT_ZH` | 倍数数量词 1–99，分别由 `en_num_term` / `zh_numeral` 派生 |
| `BIS_EN` / `BIS_ZH` | P-16.3.2 复合前缀倍增：`bis`/`tris` ↔ 双/三 |
| `BRIDGE_SUFFIX_EN` / `BRIDGE_SUFFIX_ZH` | O/S/N 桥后缀，同序同位：`oxy`/`sulfanyl`/`amino` ↔ 氧基/硫基/氨基 |
| `DIATOMIC_BRIDGE_YL` / `BRIDGE_DIATOMIC_ZH` | 双原子桥合一保留前缀 `diazenyl`/`disulfanyl` 及其中文名尾 |
| `BRIDGE_FUSION_YL` | `(中心氢化物词干, 前端名尾)` → `(中文名尾候选, 合一 en, 合一 zh)` |
| `BRIDGE_SPLIT_SUFFIX_EN` / `BRIDGE_SPLIT_SUFFIX_ZH` | 可拆桥后缀全集（含双原子桥） |
| `MONONUCLEAR_YL` | `free_en` → `(free_zh, 去氢 yl_en, 组装名中文尾)` |
| `OXO_CENTER_KINDS` / `ESTER_O_SIDE_KINDS` | 中心原子自任母体 / O-侧臂母体的 kind 集合 |
| `EXO_RING_SUF` | 环外主基后缀 `(singular, plural\|None)`，消费方 `chain_engine._exo_ring_spec` |
| `ALKOXY_YLOXY_EN` / `ALKOXY_YLOXY_ZH` | O 锚点 `-yloxy` 收拢为保留烷氧基 |
| `CHAIN_RETAINED` | C1/C2 开链保留名（酸/酰基/醛/酰胺/腈/酯） |
| `AMIDO_RETAINED` / `AZANE_PAREN_SUF` | N-酰基保留式与 azane 内层组括起白名单 |
| `ZH_DIGITS` / `HS_NUMBER` / `zh_bridge_root` | 中文数字词与桥后缀前的去「基」函数 |

`layer5/stems.py` 在词表之上派生：`_en_stem`/`alkane_en`/`alkane_zh`/`zh_num`/`zh_stem` 产碳数词干，`_metal_en_prefix` 取 `MULT_EN` 加金属名、`_metal_zh_suffix` 取 `MULT_ZH` 加中文金属名、`join_anion_names` 把酸转 `-ate`/补「根」、`join_metal_salt_names` 是盐后缀总入口。英文侧无 `zh_stem` 的对应函数：C11+ 词干由 `en_num_term(n)` 去尾 `a` 得到。

`MONONUCLEAR_HYDRIDES` 每行是五元组 `(元素, 中文名, 零价去氢名对, 桥后缀对, 组装名中文尾)`，其下游四张表全由它派生，中英永远成对：

| 派生表 | 键 → 值 |
|---|---|
| `MONONUCLEAR_BY_ELEMENT` | 原子序数 → 该元素的默认 free 名 |
| `MONONUCLEAR_ZERO_YL` | `(free_en, free_zh)` → 零价去氢名对 |
| `MONONUCLEAR_BRIDGE` | `(free_en, free_zh)` → 桥后缀对（仅 O/N/S 三行） |
| `MONONUCLEAR_YL` | `free_en` → `(free_zh, 去氢 yl_en, 组装名中文尾)` |
| `PHOSPHORYL_STEMS` | P 酰基词干元组（避免 P 被当碳中心） |

键取 `(en, zh)` 二元组而非单独英文名：单核氢化物是同一化学实体的两个词形，去氢与桥接必须同时对两者成立。消费点是 `assembler._mononuclear_radical_names` 与 `free_to_yl`——后者由 `_fg_prefix` 成对约束 `_mononuclear_en` 与 `_mononuclear_zh`，任一语言未命中即整体返回 `None`。

`PHOSPHORUS_STEM_BY_OXO`、`SULFUR_STEM_BY_OXO`、`NITROGEN_STEM_BY_FREE_DOUBLE` 按 `=O` 数或自由价键级选定锚点母体词干（`phosphanyl`/`phosphoryl`、`sulfane`/`sulfinyl`/`sulfonyl`、`azane`/`imine`），缺一则锚点与饱和类似物同形。

## kind → 词尾分派

`chain_engine._KIND_TABLE` 每个 entry 是冻结数据类 `_Chain` 实例，`_chain_names` 是唯一渲染入口，`assembler._names_for` 是 L5 派发枢纽。

| kind | `en_suf` | `zh_suf` | `coda` |
|---|---|---|---|
| `alcohol` | `ol` | 醇 | `an` |
| `ketone` | `one` | 酮 | `an` |
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
| `acyl` | `oyl` | 酰基 | `an` |
| `thiol` | `thiol` | 硫醇 | `ane` |
| `amine` | `amine` | 胺 | `an` |
| `aldehyde` | `al` | 醛 | `an` |
| `nitrile` | `enitrile` | 腈 | `an` |
| `amide` | `amide` | 酰胺 | `an` |
| `acyl_halide` | `oyl {卤}` | 酰{卤} | `an` |
| `radical` | `yl` | 基 | `an` |

要点：

- `acyl_halide` 的默认 entry 是 `_ACYL_HALIDE_BY_HAL[Cl]`，`_names_for` 按 `parent["hal_z"]` 覆盖为对应卤素 spec；C1/C2 保留名与苯甲酰由 `variant` 提供。
- C1/C2 保留名由 `constants.CHAIN_RETAINED` 经 `_retained_plain` 挂进 `variant`；苯环单取代由 `chain_engine._BENZENE_RETAINED` 在 `_names_for` 中按 `scaffold_id` 覆盖。
- `ester` 在 `parent["thio_side"]` 成立时由 `_names_for` 改写为 `coda="ane"`、`en_suf="thioate"`、`zh_suf="硫"`、`variant=None`，烯/炔基座同步改为 `enethioate`/`烯硫` 与 `ynethioate`/`炔硫`。
- 多 FG（`mult > 1`）由 `_generated_mult_fields` 把 `MULT_EN`/`MULT_ZH` 注入 `en_suf`/`zh_suf` 并强制 `coda="ane"`、`need=mult`、`no_loc="none"`，故「二醇/二酸」与单 FG 共用同一条渲染路径。
- `thiol` 的两个不饱和段取 `("ene","烯")`/`("yne","炔")`（P-57 保留 e）；`radical` 的 `plain_fn=_radical_plain` 在自由价为 1 位时直接由词干拼 `-yl`/基。
- `_names_for` 在 `kind == "radical"` 且母体带 `radical_anchor_element` 时转向 `_mononuclear_radical_names`，不查 `_KIND_TABLE`。

## 含氧酸词尾表

`chain_engine._OXO_TAIL` 以 `(oxo_kind, 中心酸式氢数)` 为键给出双语功能母体词尾；查表函数 `_oxoacid_tail` 读 `parent["oxo_kind"]` 与 `parent["n_oh"]`，表外返回 `None`。该函数经 `_Chain.plain_hook` 挂在 `phosphate`/`phosphonate`/`sulfate` 三个 kind 上，`_chain_names` 先走钩子分支，再查碳数词表。

| `oxo_kind` | 酸式氢 | EN | ZH |
|---|---|---|---|
| `phosphate` | 3 | `phosphoric acid` | 磷酸 |
| `phosphate` | 2 | `dihydrogen phosphate` | 磷酸二氢 |
| `phosphate` | 1 | `hydrogen phosphate` | 磷酸氢 |
| `phosphate` | 0 | `phosphate` | 磷酸 |
| `phosphonate` | 2 | `phosphonic acid` | 膦酸 |
| `phosphonate` | 1 | `hydrogen phosphonate` | 膦酸氢 |
| `phosphonate` | 0 | `phosphonate` | 膦酸 |
| `sulfate` | 2 | `sulfuric acid` | 硫酸 |
| `sulfate` | 1 | `hydrogen sulfate` | 硫酸氢 |
| `sulfate` | 0 | `sulfate` | 硫酸 |

## 语序：酯 / 硫代酯 / 盐 / 含氧酸中心母体

### 酯

`assembler.join_ester_name` 让 O-侧臂作前缀、酸侧作主体：英文 `{臂} {母体}`，中文 `{母体}{臂}酯`。臂由 `_o_side_arms` 从 `substituents[]` 中筛 `o_side` 标记者得到，交 `_join_o_side_arms(group=False, …)` 拼接：同名臂用 `MULT_EN`/`MULT_ZH` 倍增，任一臂带 `paren` 时英文改用 `BIS_EN` 的 `bis(…)`，异名臂依次平铺。中文臂词形由 `_zh_alkoxy_part` 渲染，英文侧不围栏。

### 硫代酯

`parent["thio_side"]` 成立即硫代羧酸 S-酯（P-65.6.3.3.7.1），语序为 `S-{臂} {母体}` / `S-{臂}{母体}酯`——位次符号 `S-` 在最前，酯词尾接在母体名后。该分支在 `join_ester_name` 内单独构造，S 侧臂的中英围栏分别走 `_fenced_arm_en` 与 `_fenced_arm_zh`。词尾由 `_names_for` 改写为 `thioate`/`硫`；苯环母体取 `_benzene_retained("benzenecarbothioate", "苯硫代甲酸")`，环外母体在 `_exo_ring_spec` 中取 `("carbothioate", "硫代甲酸")`。

### 盐

`stems.join_metal_salt_names` 表达语序反转：英文金属名作前缀、且仅当 `en.endswith("ate")` 时拼接，数量词取 `MULT_EN`（`_metal_en_prefix`）；中文以金属名替换末尾的「根」（`zh[:-1] + suf`），数量词取 `MULT_ZH`（`_metal_zh_suffix`），得「二钠」式。盐酸盐走 `acid_salt`/`acid_salt_zh`：英文空格相接，中文直接相接。盐后缀由 `namer._apply_salt_suffix` 在组装之外施加。

### 含氧酸中心母体

`kind in OXO_CENTER_KINDS` 时 `join_kind_name` 把名称交给 `assembler.join_oxoacid_name`，整名按「取代前缀 + O-侧臂 + 功能母体词尾 + 金属盐」构造：

- 无臂且有金属：酸式盐，英文 `{metal} {tail}`、中文 `{tail}{metal_zh}`，不出现「酯」字；
- 无臂且 `parent["n_om"] > 0`：游离酸式根，中文补「根」；
- 有臂且有金属：英文 `{metal} {alk} {tail}`、中文 `{tail}{alk}酯 {metal}盐`；
- 有臂无金属：英文 `{alk} {tail}`、中文 `{tail}{alk}酯`。

臂围栏分两路：`phosphate` 走 `_fenced_arms_phosphate`（若全部臂名均带前导立体描述符则由 `_all_arms_stereo_lead` 整体围栏，否则按 `oxo_arm_fence` 或「多个位次段」判定），其余中心母体走 `_fenced_arms`。

## 取代基词形的成对产出

保留式取代基以「英文取保留式、中文取对应基名」成对产出：

| 结构 | EN | ZH | 依据 |
|---|---|---|---|
| N-酰基单取代（方法 1） | `acetamido` / `formamido` / `benzamido` | 乙酰氨基 / 甲酰胺基 / 苯甲酰胺基 | P-66.1.1.4.3，查 `constants.AMIDO_RETAINED` |
| 长链/烯酰/杂环羰酰（方法 2） | `acetylamino` 等 | …酰氨基 | 不入表，走 `free_to_yl` |
| 质子化伯胺 | `azaniumyl` | 铵基 | P-62.4.1，锚定保留叶 |
| 苯胺去氢 | `anilino` / `4-chloroanilino` | 苯胺基 / 4-氯苯胺基 | P-62.1.1，`_mononuclear_en` / `_mononuclear_zh` |
| N-取代磺酰胺基 | `butyl(methyl)sulfamoyl` | …磺酰基 | P-66.1.1.4.2，与 `sulfamoyl` 融合 |
| O 锚点 `-yloxy` | `ethyloxy` → `ethoxy` | 乙基氧基 → 乙氧基 | P-66.5.2.1.2，`ALKOXY_YLOXY_EN` / `ALKOXY_YLOXY_ZH` |

成对性靠「英文词干承担判据、中文侧同名同形」维持：`AMIDO_RETAINED_EN` 供 L3 判定免括号，中文两侧词形一致（均为「乙酰氨基」）；`AZANE_PAREN_SUF` 决定方法 2 中哪些尾缀须把内层组整体括起再加 `amino`；`free_to_yl` 的第三位返回值 `need_paren` 只由英文名形判定（`methylamino`、带环取代基的 `anilino`），中英共用该布尔。

## 中文侧特有处理

- **词干路径**：`_chain_zh_base` 用 `zh_stem` 剥去 `alkane_zh(n)` 的「烷」得词干；`_chain_stem_pair` 在 `unsat_polyol` 且 `zh_full`/`cyclic` 时改用完整中文全名。
- **酯 O-侧基名**（`_zh_alkoxy_part`）：简单烃基省「基」（乙酸乙酯）；仅带位次者保留「基」（乙酸丙-2-基酯）；其余加围栏，前端已含圆/方括号时升级为方括号。
- **含氧酸 O-侧臂**（`_oxoacid_arm_zh`）：简单基去「基」；整串为中文数字（`ZH_DIGITS` +「十」）时补「烷基」；含 `-` 或以 `(` 开头的复合名原样保留。
- **「基」的回补**（`_zh_front_ji`）：中文复合前端的「基」被上游切掉时，在拆分点补回；酰基前端走融合式，或 `_front_needs_enclosure` 为假时不补。
- **指示氢前缀的取代式落位**：`_with_indicated_h` 把 L4 定好的位次（`_indicated_h_prefix`，同位次去重）落到词干前，并**取代**词干自带前缀，避免 `1H-7H-` 双写；`_nh_only` 与 `_ring_n_free` 共同决定未取代的 NH 环是否补指示氢；`_stem_prefix_stale` 在词干自带位次已成羰基碳时丢弃该前缀。`_ensure_fused_stem` 为未注册稠环注入词干时同样带上前缀。
- **环阳离子与保留名**：`join_ring_cation_suffix` 英文缀 `-{位次}-ium`，中文沿用母体名；`zh_1h_parent` 在环被取代时给中文保留名补 `1H-`。

## 围栏与排序

- `_enclose(s)`：已含圆括号则升为方括号，否则加圆括号。它是双语共用的围栏原语。
- `oxo_arm_fence(name, sub, mol)` 判含氧酸中心母体的臂名是否整体括起：前导 `(` 返回真（仅为立体描述符），前导 `[` 返回假（已是完整围栏），自带两个以上位次段返回真，带 `paren` 的复合臂挂环上返回真。
- `_fenced_arm_en` / `_fenced_arm_zh`：酯 O-侧臂围栏，前导位次或自带括号时整体括起。
- `_bridge_body(base, suf, merge)`：平铺式为 `_enclose(base) + suf`，桥后缀留在括号外；双原子桥整段再括一层；`merge=True` 时前端已含方括号则同括。
- `_sbridge_flat_stem`：磺酰/亚磺酰桥 + 无取代直链 `-yl` 前端，英文侧平铺不加围栏。
- **中英同步判据**：`_split_bridge_suffix_zh` 的唯一判据是英文 stem 是否拆分（先调 `_split_bridge_suffix`），故中英围栏必然同形；`_bridge_self_enclosed` 标记经 `meta` 回传，使 L5 对自含围栏的桥前端跳过二次加括号。
- **排序**：`tools/re.alpha_order_key` 是 P-14.5 字母数字序的唯一实现，先比非斜体字母序列、再比首字母前位次；中英两侧共用同一键（`_mult_rows`、`_build_prefix`、`_join_o_side_arms`、`_phosphoryl_sub_names` 均以它作 sort key），故取代基在两语中排列一致。`alkyl_alpha_key` 是其循环剥离步骤；`_alpha_key` 是 `assembler` 内只保留字母的简化键，用于 N-芳基-N-某基胺的 `anilino` 引用序判定。`re.py` 中的 `normalize_en`/`normalize_zh`/`nospace` 是判分口径，命名管线不调用。

## 递归命名：双语缓存与宿主根传播

`tools/common_names.CommonNameCache` 以 SMILES 为键存 `NameResult`（含 `en`/`zh`），提供 `get`、`put`（容量满抛 `ValueError`）、`clear`。`SMILESNNamer.__init__` 默认构造 20000 条容量的缓存；`name` 先查缓存，未命中才走 `_pipeline` 并 `_cache_put`（容量满的 `ValueError` 静默忽略），写前经 `_canonical_result` 用 `CanonicalRankAtoms` 归一 `meta.parent_chain`。

`root_ctx` 是 `(root_mol, to_root)`——根分子对象与「当前块原子序号 → 根分子原子序号」映射。`_name_mol` 在 `root_ctx is None`（顶层整分子）时取解离盐后的 `organic` 与恒等映射，并另开 2000 条容量的 `run_cache`，使内部 `*` 片段名只在本次分子运行内共享；`root_ctx` 非空时沿用上游传入的根分子、映射与缓存。该上下文写入 `info["root_ctx"]`。

L3 的 `as_substituent._radical_yl_from_sub` 解包 `root_ctx`、构造锚定块的 `anchored_to_root` 映射后递归调用 `_name_mol(..., root_ctx=...)`；`_fix_rs_with_real` 用原始根分子 CIP 校正被 `*` 顶替后翻转的异头位 R/S，fresh 计算与缓存命中都按当前根分子校正一次。桥前端自含围栏的标记同样经 `meta` 回传，保证中英两侧的围栏形态一致。
