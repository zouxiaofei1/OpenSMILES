# 中英双语命名约定 (Bilingual Naming Convention)

> **概念页面** | **最后更新:** 2026-09-13

---

## 概述

NamePredict 的核心设计决策之一是从管线起点到终点**同步生成英文与中文两套 IUPAC 名称**。不存在"先生成英文再翻译为中文"的二次处理阶段——每一层、每一个函数，都同时产出或传递双语名称。这一约定贯穿全部 6 层架构，是 NamePredict 区别于仅支持单语言的命名引擎的根本特征。

**管线中的双语数据流：**

```
SMILES 输入 (语言无关)
    ↓
Layer0: salt_meta = {metal: "sodium", metal_zh: "钠"}   ← 双语元数据
    ↓
Layer1–4: parent.stem_en / parent.stem_zh, substituent.en / substituent.zh   ← 双语中间表示
    ↓
Layer5: (en, zh) tuples → NameResult(en="ethanol", zh="乙醇")   ← 双语输出
```

## 核心约定：`(en, zh)` 元组

整个 Layer5 组装管线中，所有名称生成函数遵循一个统一的返回值约定：**返回 `tuple[str, str] | None`，其中第一个元素为英文名、第二个元素为中文名**。这一约定确保了函数的可组合性——上游的输出可以直接作为下游的输入，无需额外转换。少数辅助函数的形状不同：`free_to_yl`（`src/namepredict/layer5/assembler.py:280`）返回 `(en, zh, need_paren)` 三元组，第三位是中英两侧共用的围栏标记；`_anilino_en`（`assembler.py:245`）只加工英文侧。

### 管线的双语流水线

以 `assemble()` 函数（`src/namepredict/layer5/assembler.py:562`）为例，每一步都操作 `(en, zh)` 对：

```
1. _ensure_fused_stem(numbered)                 → bool              未注册稠环词干注入（assembler.py:352）
2. _names_for(kind, n, numbered)                → (en, zh) | None   母体名称（assembler.py:372）
3. join_hydro_prefix(names, numbered)           → (en, zh)          hydro + 指示氢前缀（assembler.py:510）
4. join_ring_cation_suffix(numbered, names)     → (en, zh)          环内 N+/O+ 母体名缀 -{位次}-ium（assembler.py:523）
5. _prefix_for(numbered, kind, n)               → (pre_en, pre_zh)  取代基前缀（assembler_prefixes.py:299）
6. join_kind_name(kind, pre, names, numbered)   → (en, zh) | None   前缀+母体拼接（assembler.py:491；酯走 join_ester_name）
7. join_anion_names(numbered, en, zh)           → (en, zh)          羧酸根转换（stems.py:74）
8. join_rs_prefix(numbered, en, zh)             → (en, zh)          R/S 立体化学（stereo.py:234）
9. join_metal_salt_names(numbered, en, zh)      → (en, zh)          盐后缀（stems.py:129）
```

每一步都是纯文本转换：输入 `(en, zh)`，输出 `(en, zh)`。这种设计具有三个优点：

1. **可测试性**：每个步骤可以独立验证其中英双语输出的正确性。
2. **可追踪性**：当名称出现错误时，可以精确定位到出错的步骤（例如中文盐后缀错位但英文正确）。
3. **可扩展性**：新增官能团命名时，只需实现一个新的 `(kind, n, numbered) -> (en, zh) | None` 函数并注册到 `_names_for` 的派发链中。

### 母体名称派发中的双语约定

`_names_for` 函数（`src/namepredict/layer5/assembler.py:372`）先处理两条特例：`kind == "phosphate"` 直接调 `layer5/phosphate.py` 的 `phosphate_names`（`assembler.py:375`，定义于 `phosphate.py:111`）组装磷酸整名；环外主基（`kind` 落在 `constants.EXO_RING_SUF`，`constants.py:179`）交 `_exocyclic_ring_names`（`assembler.py:78`），覆盖 acid/aldehyde/ester/amide/nitrile/acyl 六类，苯单取代时不接手、回落 chain_engine 保留名。之后按 `kind` 查 `chain_engine._KIND_TABLE`（14 个 `_Chain` spec：alcohol/ketone/alkane/acid/sulfonic/ester/acyl/thiol/amine/aldehyde/nitrile/amide/acyl_halide/radical，按 `scaffold_id` 运行时注入环前缀/稠环词干），命中即返回；杂原子锚点自由基走 `_mononuclear_radical_names`（`assembler.py:292`）的单核氢化物去氢管线（`assembler.free_to_yl`，`assembler.py:280`），苯单取代的 `phenyl` 保留名由 chain_engine 的 `variant` 提供；都不命中时由 `_parent_stem_names`（`assembler.py:410`）回退。每个 worker 都是返回 `(en, zh)` 或 `None` 的双语函数——返回 `None` 时派发器降级到下一策略。这种"尝试-失败-降级"模式在双语的上下文中尤其重要——如果一个命名策略返回了英文名但无法生成中文名（或反之），则整个结果不被接受。

## 双语词表：constants.py（跨层）与 stems.py（词干生成）

双语词表分置两处：`src/namepredict/constants.py` 存放**跨层共享的词表常量**（原子序数、倍数前缀、单核氢化物、保留名、桥后缀等，被 L0–L5 与 `tools/` 共同 import），`src/namepredict/layer5/stems.py` 存放**碳数词干生成器**（C1–C99 烷烃词干与阴离子/盐转换）。`stems.py` 本身也从 `constants` 取词（`stems.py:5`），两者构成"表在 constants、派生在 stems"的分工。

`constants.py` 中与双语命名直接相关的词表：

| 词表 | 位置 | 内容 |
|---|---|---|
| `en_num_term(n)` / `zh_numeral(n)` | `constants.py:56` / `constants.py:76` | 英文/中文数值词干（1–99）；英文末带 'a'（undeca/icosa…），中文 1 返空串 |
| `MULT_EN` / `MULT_ZH` | `constants.py:87` / `constants.py:88` | 由上面两函数派生的 1–99 倍数前缀字典（di/二、tri/三…） |
| `HALO_ZH` / `HALIDE_EN` | `constants.py:37` / `constants.py:38` | 卤素中文名（氟/氯/溴/碘）与卤离子英文名（fluoride…） |
| `ALKALI_EN` / `METAL_ZH` | `constants.py:110` / `constants.py:111` | 碱金属英文名 → 中文名（锂/钠/钾） |
| `AMIDO_RETAINED` / `AMIDO_RETAINED_EN` | `constants.py:94` / `constants.py:99` | N-酰基保留式 amido 的 (en, zh) 与英文名集合（供 L3 判定免括号） |
| `MONONUCLEAR_HYDRIDES` 及派生表 | `constants.py:126`–`constants.py:143` | P-15.4.1 表 2.1 单核母体氢化物：free 名、(free_en, free_zh) 键的零价去氢名 `MONONUCLEAR_ZERO_YL`、桥后缀 `MONONUCLEAR_BRIDGE`、去氢名 `MONONUCLEAR_YL` |
| `zh_bridge_root(name)` | `constants.py:90` | 桥后缀前的中文烃基去尾「基」（甲基→甲、叔丁基→叔丁、环己基→环己） |
| `SIMPLE_ALKOXY_NO_PAREN` | `constants.py:146` | 简单保留烷氧基作前缀免括号白名单 |
| `BIS_EN` / `BIS_ZH` | `constants.py:173` / `constants.py:174` | P-16.3.2 复合前缀倍增（bis/tris、双/三/四） |
| `BRIDGE_SUFFIX_EN` / `BRIDGE_SUFFIX_ZH` 及 `-yl` 型派生 | `constants.py:175`–`constants.py:178` | O/S/N 桥后缀（oxy/sulfanyl/amino ↔ 氧基/硫基/氨基）与 `-yl` 型桥基对 |
| `EXO_RING_SUF` | `constants.py:179` | 环外主基系统名后缀表（acid/aldehyde/ester/amide/nitrile/acyl 六类的 singular/plural/位次要求） |
| `AZANE_PAREN_SUF` | `constants.py:187` | azane 单取代基内层组须加括号的尾缀白名单 |
| `ALKOXY_YLOXY_EN` / `ALKOXY_YLOXY_ZH` | `constants.py:188` / `constants.py:192` | O 锚点 `-yloxy` 收拢为 IUPAC 保留烷氧基（ethyloxy→ethoxy ↔ 乙基氧基→乙氧基） |
| `CHAIN_RETAINED` | `constants.py:196` | C1/C2 开链英文保留名（formic/acetic、formamide/acetamide…）及其中文对 |
| `RETAINED_FUSION_ALIASES` | `constants.py:204` | 稠合组装名 → 保留名整名替换（benzo[c]furan→2-benzofuran/2-苯并呋喃…） |
| `HS_NUMBER` / `ZH_DIGITS` | `constants.py:213` / `constants.py:53` | 天干串「甲乙丙丁戊己庚辛壬癸」与数字汉字「一二三…九」 |

### 基表结构

**C1-C10 保留名基表**（`src/namepredict/layer5/stems.py:8-15`）：

```python
_ALKANE_EN_BASE = {1: "methane", 2: "ethane", ..., 10: "decane"}
_ALKANE_ZH_BASE = {1: "甲烷",  2: "乙烷",   ..., 10: "癸烷"}
```

英文和中文基表严格一一对应，使用相同的整数 key。这种 `{n: str}` 的字典结构确保了查找效率和数据的不可变性。

**C11+ 英文词干复用 `constants.en_num_term`**（数量词与母链碳数同源）：`_en_stem`（`stems.py:39`）对 n≥11 直接 `en_num_term(n)[:-1]` 去尾 'a'（undec、icos、docos…），`constants.MULT_EN/MULT_ZH`（`constants.py:87-88`）是同一组数值词干在倍数前缀上的复用。中文则由 `zh_num(n)` 函数动态生成 `十一烷`、`十二烷`…。

### 烷烃词干与 FG 名称的生成

`stems.py` 提供烷烃英中词干生成器：`alkane_en(n)`（`stems.py:47`）/`alkane_zh(n)`（`stems.py:53`）生成 C1-C99 烷烃全名，`_en_stem(n)`（`stems.py:39`）/`zh_stem(n)`（`stems.py:31`）取去后缀词干。**FG 名称不由 stems 逐官能团派生**（无 `alcohol_en`/`acid_en`/`amide_en` 等函数）——由 chain_engine 用 `_en_stem` + `_Chain.en_suf`/`zh_suf` 直接拼接（`chain_engine.py:276` 的 `_chain_plain`），C1/C2 保留名（formic/acetic、甲酸/乙酸）由 `constants.CHAIN_RETAINED`（`constants.py:196`）经 `_retained_plain`（`chain_engine.py:34`）提供：

| 官能团 | 生成方式 | 示例 (C2) | 示例 (C2 中文) |
|--------|---------|-----------|---------------|
| 烷烃 | `alkane_en(n)` / `alkane_zh(n)` | `ethane` | `乙烷` |
| 醇 | chain_engine `alcohol` entry（`_en_stem` + `-ol`） | `ethanol` | `乙醇` |
| 酸 | chain_engine `acid` entry（C1/C2 走 `CHAIN_RETAINED` 保留名） | `acetic acid` | `乙酸` |
| 醛 | chain_engine `aldehyde` entry | `acetaldehyde` | `乙醛` |
| 酰胺 | chain_engine `amide` entry | `acetamide` | `乙酰胺` |
| 腈 | chain_engine `nitrile` entry | `acetonitrile` | `乙腈` |
| 酯 | chain_engine `ester` entry | `acetate` | `乙酸` |

公开词干表 `ALKANE_EN`/`ALKANE_ZH` 通过 `_fill(fn, lo=1, hi=99)`（`stems.py:138`）自动填充为完整字典（`stems.py:147-148`），从 C1 覆盖到 C99，确保长链名称无需手动维护。中文的 C1/C2 酸、醛、酰胺、腈等使用保留名（如 `甲酸`/`乙酸`、`甲醛`/`乙醛`），与英文保留名（`formic acid`/`acetic acid`、`formaldehyde`/`acetaldehyde`）保持对应。

盐与阴离子转换同在 `stems.py`：`acid_to_anion_en`（`stems.py:60`）/`acid_to_anion_zh`（`stems.py:69`）、`join_anion_names`（`stems.py:74`）、`join_metal_salt_names`（`stems.py:129`）及其内部 `_metal_en_prefix`（`stems.py:91`）/`_metal_zh_suffix`（`stems.py:96`）。

### amido 保留式的双语词形

N-酰基取代基有两条表达路径，双语词形由 `constants.AMIDO_RETAINED`（`constants.py:94`）区分：

| 路径 | 英文 | 中文 | 说明 |
|---|---|---|---|
| retained amido（P-66.1.1.4.3 方法 1） | `acetamido` / `formamido` / `benzamido` | `乙酰氨基` / `甲酰胺基` / `苯甲酰胺基` | 仅 acetyl/formyl/benzoyl 三词入表，由 `_mononuclear_radical_names` 查表（`assembler.py:313`）；中文按 P-66 CN 译本取「酰氨基」式 |
| acylamino（方法 2） | `acetylamino` 等 | `…酰氨基` | 长链/烯酰/被取代苯甲酰/杂环羰酰，不入表，走 `free_to_yl` 的 `_mononuclear_en`（`assembler.py:250`） |

amido 与 acylamino 的区分由**英文词干**承担（`AMIDO_RETAINED_EN` 供 L3 判定免括号，`as_substituent.py:95`），中文两者词形一致（`乙酰氨基`）。formyl/benzoyl 的中文（`甲酰胺基`/`苯甲酰胺基`）待与 gold 统一。

### 其他保留式取代基的双语词形

| 结构 | 英文 | 中文 | 依据（实现位置） |
|---|---|---|---|
| 质子化伯胺 NH&#8323;&#8314; | `azaniumyl` | `铵基` | P-62.4.1（锚定保留叶 `*[NH3+]`，`tools/anchored_table.py:72`；同类 `methylazaniumyl`/`dimethylazaniumyl` 见 `:73-74`） |
| 苯胺去氢（N-苯基氨基） | `anilino`；带环取代基 `4-chloroanilino` | `苯胺基`；`4-氯苯胺基` | P-62.2.1.1（`assembler._anilino_en`，`assembler.py:245`） |
| N-取代磺酰胺基 | 与 `sulfamoyl` 融合：`butyl(methyl)sulfamoyl` | `…磺酰基` | P-66.1.1.4.2（`assembler._mononuclear_radical_names`，`assembler.py:305`；裸 `sulfamoyl`/氨磺酰基 为锚定保留叶 `anchored_table.py:78`） |

三者都是"英文取保留式、中文取对应基名"的成对词形：`azaniumyl`/`铵基` 由 L3 的锚定保留表直接产出，`anilino`/`苯胺基` 由 `free_to_yl` 的去氢转换产出，`sulfamoyl` 融合在 L5 单核自由基母体管线中定形。

## 中英文命名差异

虽然英文和中文 IUPAC 名称共享相同的位次号和化学结构语义，但两者的语序和词法规则存在系统性的差异。以下表格总结了主要差异：

| 特征 | 英文 | 中文 |
|------|------|------|
| **取代基前缀** | `2,4-dichloro-3-methyl-` | `2,4-二氯-3-甲基` |
| **取代基命名** | `methyl`, `ethyl`, `chloro` | `甲基`, `乙基`, `氯` (带 `基` 介词) |
| **多重度前缀** | `di`, `tri`, `tetra` | `二`, `三`, `四` |
| **位次号位置** | 在取代基名前：`2-chloro-` | 在取代基名前：`2-氯`（与英文一致） |
| **字母序** | 英文名按字母序排列前缀 | 中文仍用 `alkyl_alpha_key`（英文排序键）以保持一致 |
| **母体名称** | `butanamide` | `丁酰胺` |
| **醇后缀** | `-ol` | `醇` |
| **醛后缀** | `-al` | `醛` |
| **酸后缀** | `-oic acid` / `-ic acid` | `酸` |
| **酮后缀** | `-one` | `酮` |
| **酯命名** | `methyl benzoate` (醇在前、酸在后) | `苯甲酸甲酯` (酸在前、醇在后) |
| **磷酸酯命名** | `trimethyl phosphate` (烷基在前) | `磷酸三甲酯` (酸在前、醇在后) |
| **盐命名** | `sodium acetate` (阳离子在前) | `乙酸钠` (阳离子在末尾) |
| **盐酸盐** | `… hydrochloride`（空格连接） | `…盐酸盐`（直接相接） |
| **立体化学** | `(E,2S)-` | `(E,2S)-` (与英文一致，不翻译) |
| **环前缀** | `cyclohexane` | `环己烷` |
| **编号前缀** | `1H-pyrrole-` | `1H-吡咯`（`1H-` 保留，母体翻译） |
| **环阳离子** | `-{位次}-ium`（`pyridin-1-ium`） | 沿用母体名（`吡啶`） |

编号前缀（`1H-`、`1,2-`、`1,3,5-` 等）由 `ring_scaffold._TEMPLATES` 的 `locant_prefix` 提供（`layer2/ring_scaffold.py:35-36` 定义字段、`ring_scaffold.py:71` 起为模板注册表）：无条件注入者（`1,2-` 异噁唑、`1,3,5-` 三嗪）恒带前缀；`prefix_nh_conditional` 为真者（吡咯型 `1H-` 四唑）仅当环含未取代 NH 时注入，N 被取代时省略。

### 酯命名的语序反转

酯命名是中英文差异最显著的例子。英文的酯命名为 `{alkyl} {acyl}ate` 格式，烷基（醇部分）在酰基（酸部分）之前：

- EN: `methyl benzoate`（甲基 + 苯甲酸酯）

中文则遵循"酸在前、醇在后"的规则，格式为 `{酸}{醇}酯`：

- ZH: `苯甲酸甲酯`（苯甲酸 + 甲基 + 酯）

这一差异在代码中由 `assembler.py` 的 `join_ester_name` 函数（`assembler.py:460`）处理，分别构造英文和中文的语序。对于不饱和酯和带取代基的酯，中英文的插入位置也有不同的处理逻辑。

### 磷酸/磷酸酯的双语词形

磷酸/磷酸酯（`kind == "phosphate"`）由 L5 `layer5/phosphate.py` 的 `phosphate_names`（`phosphate.py:111`）整体组装，不经 `join_ester_name`：

- **中性酯**：EN `{臂} {tail}`（`trimethyl phosphate`），ZH `磷酸{臂}酯`（`磷酸三甲酯`）——同样是"酸在前、醇在后"的中文语序（`phosphate.py:131-150`）
- **酸式词尾**：`_tail_en`（`phosphate.py:8`）/`_hyd_zh`（`phosphate.py:17`）按剩余酸式 H 数（`n_oh`）给出 `phosphate`/`hydrogen phosphate`/`dihydrogen phosphate` ↔ `磷酸`/`氢`/`二氢`
- **带酯臂的碱金属盐**：EN 金属名前置（`disodium tridecyl phosphate`），ZH 金属名置末（`磷酸十三烷基酯 二钠盐`）——与羧酸盐相同的语序反转，由 `_ester_salt_names`（`phosphate.py:58`）构造
- **游离阴离子**（无抗衡金属）：有臂时 `{臂} {tail}`/`磷酸{氢}{臂}酯`，无臂时 `phosphate`/`磷酸根`、`dihydrogen phosphate`/`磷酸二氢根`；负电荷不标注（`_free_anion_names`，`phosphate.py:86`）

中文臂词形分两档：`_arm_zh_root`（`phosphate.py:26`）对单字根（甲基→甲、苯基→苯）去「基」用于中性酯；`_arm_ester_zh`（`phosphate.py:49`）对多位纯中文数字根（十三基→十三烷基）补「烷」以对齐 gold。相同英文的臂由 `_group_arms`（`phosphate.py:34`）合并计数后按英文序排列。

### P 酰基前缀的双语词干

P 锚点自由基母体的词干按氧化态选定（`constants.PHOSPHORUS_STEM_BY_OXO`，`constants.py:138`，由 `principal_expression.py:414` 消费）：0 个 =O 取 `phosphanyl`/`磷烷基`，1 个 =O 取 `phosphoryl`/`磷酰`（P-67.1.4.1.1.2/6）。取代基拼接中英同步（`_phosphoryl_sub_names`，`assembler.py:148`）：

- 全为简单基：首基平铺、其余括起，中英同形（`hydroxy(methyl)phosphoryl` / `羟基(甲基)磷酰基`）
- 同基倍增：用 `di-`/`二` 而非逐基括号（`dimethoxyphosphoryl` / `二甲氧基磷酰基`）
- 含复合组分：逐组分以连字符分隔、需围栏者加方括号（P-16.5.2 嵌套标记；-yl 型桥基把桥后缀挪到方括号外）
- P 酰基经 O/N/S 桥连母体（P-67.1.4.1.3 复合前缀）：整体加方括号后接桥后缀（`[hydroxy(methoxy)phosphoryl]oxy` / `[羟基(甲氧基)磷酰]氧基`）

### 盐命名的后缀位置反转

盐的命名也体现出语序差异。英文将金属阳离子放在羧酸根名称之前作为前缀（`sodium dodecanoate`），而中文将金属名放在名称末尾（`十二酸钠`）。这一转换由 `stems.py` 中的 `join_metal_salt_names` 函数（`src/namepredict/layer5/stems.py:129`）完成：

- `_salt_en`（`stems.py:101`）：将金属名作为前缀拼接到英文名（仅当英文名以 `ate` 结尾时），数量词取 `MULT_EN`
- `_salt_zh`（`stems.py:107`）：将金属中文名替换中文名的末尾 `酸根` → `酸{金属}`，数量词取 `MULT_ZH`（`二钠`）

盐酸盐（HCl salt）采用不同的机制：通过 `acid_salt` / `acid_salt_zh` 字段追加后缀（空格式 `…… hydrochloride` / 直接相接 `……盐酸盐`），由 `_with_acid_salt`（`stems.py:121`）拼装。

磷酸/磷酸酯是例外：盐形态（含 `disodium … phosphate` / `磷酸…酯 二钠盐` 语序）由 L5 `phosphate_names` 在整名内组装，`namer._apply_salt_suffix`（`namer.py:192`）对 `parent_kind == "phosphate"` 直接返回，通用盐后缀路径对该 kind 不生效。

### 取代基前缀的排序一致性

虽然英文前缀按取代基英文名的字母序排列（符合 IUPAC P-14.5），但中文前缀的排列仍使用 `alkyl_alpha_key`——一个基于英文名的字母数字序排序键，定义在 `tools/re.py:61`（剥掉 `sec-`/`tert-`/`N-`/括号/前导位次/前导立体组后比较），由 L5 的 `assembler_prefixes.py:7`、`assembler.py:150` 与 L4 的 `numbering_engine.py:232`、`_chain_orient.py:4` 共享。这确保了在双语输出中，取代基的排列顺序始终保持一致：中英文的前缀中的基团顺序总是相同的，避免了因语种不同导致取代基排列顺序不一致的混淆。

### 桥平铺式围栏的中英同形

O/S/N 桥后缀（`-oxy`/`-sulfanyl`/`-amino` ↔ `氧基`/`硫基`/`氨基`，表在 `constants.BRIDGE_SUFFIX_EN`/`BRIDGE_SUFFIX_ZH`，`constants.py:175-176`）采取"平铺式"：括号闭在前端基的 `-yl`/`基` 之后，桥后缀留在括号外——`(4-methoxyphenyl)sulfonyl` ↔ `(4-甲氧基苯基)磺酰基`。中文侧的拆分由 `_split_bridge_suffix_zh`（`assembler_prefixes.py:173`）完成，其**唯一判据是英文 stem 是否拆**（先调 `_split_bridge_suffix`，`assembler_prefixes.py:120`），故中英围栏必然同形。

前端是否拆分由 `_front_needs_enclosure`（`assembler_prefixes.py:99`）判定：

- 整括不拆：简单保留基（`methyl`/`benzyl`）、无取代直链 -yl（`propan-2-yloxy` 平铺）、酰基前端（`acetyloxy`/`benzoyloxy`）、自由价在端碳的链基（`…phenyl)methyl`、`5-(X)pentyl`）、以及磺酰/亚磺酰基前端（其围栏由 L3 定形）
- 须围栏：自由价碳带手性描述符（`(2R)-2-amino-2-carboxyethyl`）、环型内嵌位次（`…oxan-2-yl`）以及已自带方括号/括号的复合前端；`amino` 桥在括号后仍接位次前缀时再叠一层方括号（P-63.2.2.1.2：`[[X]amino]propanoyl`）
- 磺酰/亚磺酰桥 + 无取代直链 -yl 前端：英文侧平铺不加围栏（`propan-2-ylsulfonyl`，`_sbridge_flat_stem`），中文侧仍括注

前导立体描述符须整体围栏：取代基名以 `(1Z)-`/`(2R,4R)-` 开头时中英两侧同时置位（`_STEREO_LEAD_RE`，`assembler_prefixes.py:87`），如 `5-[(1Z)-prop-1-enyl]…` ↔ `5-[(1Z)-丙-1-烯基]…`。

同基倍增（N 桥）的中文词根由 `zh_bridge_root`（`constants.py:90`）去尾「基」后再接桥后缀：甲基→甲、叔丁基→叔丁、环己基→环己、丙-2-基→丙-2-，故英文 `dimethylamino` 对应 `二甲氨基`（而非 `二甲基氨基`）。英文侧以 `di-`/`tri-` 表达同一倍增（`assembler.py:327-330`）。

### 环阳离子：英文 `-{位次}-ium` 与中文母体名

净正电荷分子的环内 N&#8314;/O&#8314; 由 `join_ring_cation_suffix`（`assembler.py:523`）缀在母体名上（P-62.4.1）：英文取 `-{位次}-ium`（`pyridin-1-ium`），色烯型氧鎓保留名整词干替换（`chromene` → `chromenylium`）；中文侧沿用母体名（`吡啶`、`色烯`），不译出 `-ium` 的对应词形。位次取 `parent.numbering_scaffold.labels`，缺失时回落到链内序号。

### 位次省略的双语同步

不饱和位次的省略由 L4 的 `omit_ene_locant`/`omit_yne_locant` 与 L5 `_Chain` 的 `ene_loc_omit`/`yne_loc_omit` 共同决定，中英两侧共用同一判定、同一字段：

- 开链烃二核烯与二/三核炔（`ethene`/`乙烯`、`ethyne`/`乙炔`、`propyne`/`丙炔`）位次 1 隐含省略（P-14.3.4.2(d)）；C3 烯仍保留（`prop-1-ene`/`丙-1-烯`）
- 带 FG 后缀的炔恒保留炔位次（`spec.yne_loc_omit` 为假时不给省略，P-14.3.4 例外：`prop-2-ynoic acid`）
- 环单烯在双键起点与 FG/自由价同为 1 时省去冗余的 `1` 并融合（P-31.1.2：`cyclohexen-1-yl`），FG 位次仍显式保留

## `zh_stem` 转换

`zh_stem` 函数（`src/namepredict/layer5/stems.py:31-36`）是一个关键的中文词干提取工具。它从中文全名中剥离末端的官能团/母体后缀，返回"裸词干"供后续拼接。

**后缀剥离表**（`src/namepredict/layer5/stems.py:17`）：

```python
_ZH_SUFFIXES = ("酰胺", "酰氯", "硫醇", "烷", "醇", "酸", "醛", "腈", "胺", "酮", "烯", "炔")
```

**转换示例：**

| 输入 (`zh_full`) | 输出 (词干) | 用途 |
|-----------------|-----------|------|
| `十一烷` | `十一` | 构造 FG 位次名称，如 `十一烷-2-醇` |
| `乙醇` | `乙` | 构造带位次的醇名，如 `乙-1,2-二醇` |
| `丁酰胺` | `丁` | 构造取代酰胺名 |
| `环己烷` | `环己`（剥「烷」） | 环母体的词干 |

`zh_stem` 的工作原理是按 `_ZH_SUFFIXES` 元组中的顺序依次尝试匹配——一旦某个后缀匹配成功，就切除该后缀并返回剩余部分。匹配顺序从长后缀（`酰胺`）到短后缀（`烯`、`炔`），避免短后缀误匹配长后缀的尾部。如果没有任何后缀匹配，函数返回原始输入——这意味着它可以安全地用于任何中文名称字符串。

在 `chain_engine.py` 的环系命名中广泛使用了 `zh_stem`——`_chain_zh_base`（`chain_engine.py:255`）用它对 `alkane_zh(n)` 去「烷」得词干；`assembler._names_for` 对 `scaffold_id=="carbocycle"` 的母体注入 `cyclic=True`（`assembler.py:398-401`，恒加 `环`/`cyclo` 前缀），`_chain_names`（`chain_engine.py:332`）再在其上追加烯/醇等后缀。

## `zh_num`：中文数字生成

`zh_num(n)`（`src/namepredict/layer5/stems.py:20-28`）将整数 n（1-99）转换为中文链数字：

- 1-10：取天干串 `constants.HS_NUMBER`（`constants.py:213`，「甲乙丙丁戊己庚辛壬癸」），即该档用甲/乙/丙…而非一/二/三…
- 11-19：`十` + 个位（`十一` `十二`…而非 `一十` `二十`…）
- 20-99：十位 + `十` + 个位（`二十五`、`九十九`），十位/个位汉字取自 `_DIGIT_ZH`（`stems.py:16`，即 `"零" + constants.ZH_DIGITS`）

消费方是 `alkane_zh(n)`（`stems.py:53`）：C1–C10 直接查 `_ALKANE_ZH_BASE`（`stems.py:12`），C11+ 用 `zh_num(n) + "烷"` 生成，故该函数是中文长链烷烃名称（`十一烷`、`三十五烷` 等）的数据源。`zh_num` 不存在对应的英文函数，因为英文 C11+ 词干由 `constants.en_num_term`（数量词同源）去尾 'a' 派生（undec/icos/docos…）；倍数前缀侧的中文数值词另走 `constants.zh_numeral`（`constants.py:76`，1 返空串，与 `zh_num` 的天干口径不同）。

## 递归命名的 `(en, zh)` 缓存与宿主根传播

管线只产出一种命名风格（通用名），`(en, zh)` 对不带风格开关；双语结果通过**结果缓存**与**宿主根上下文**在递归调用之间传播。

**结果缓存**（`src/namepredict/tools/common_names.py`）：`CommonNameCache`（`common_names.py:7`）以 SMILES 为键存 `NameResult`（含 `en`/`zh` 两字段），提供 `get`（`common_names.py:15`）、`put`（`common_names.py:19`，容量满抛 `ValueError`）、`clear`（`common_names.py:28`）。`SMILESNNamer.__init__`（`namer.py:277`）默认建 20000 条容量的缓存；`name`（`namer.py:281`）先查缓存，未命中才走 `_name_uncached`（`namer.py:269`）并写回。顶层整分子另行开一个 2000 条容量的 `run_cache`（`namer.py:222`）：内部 `*` 片段名只在本次分子运行内共享，跨根缓存会带入别的宿主的异头立体。

**宿主根上下文**：`_name_mol`（`namer.py:208`）把 `root_ctx = (根分子, 原子→根索引映射)` 写入 `info["root_ctx"]`（`namer.py:227`），取代基命名端据此把 `(en, zh)` 里的 R/S 描述符按真实宿主根重算（`as_substituent.py:87` 的 `_fix_rs_with_real`），因此缓存命中的片段名不会带跨根立体污染。写缓存前还会经 `_canonical_result`（`namer.py:252`）把 `meta.parent_chain` 重写为规范原子序——缓存键是子结构 SMILES，链序不能依赖随母体变化的切分点。

## 盐元数据的双语结构

Layer0 的 `dissociate_salt` 函数（`src/namepredict/layer0/salt.py:81`）在检测到盐结构时，返回的元数据字典包含完整的中英双语字段（碱金属盐由 `_meta_metal`（`salt.py:52`）产出，盐酸盐由 `_meta_hcl`（`salt.py:60`）产出）：

| 字段 | 示例值 | 说明 |
|------|--------|------|
| `metal` | `"sodium"` | 英文金属名 |
| `metal_zh` | `"钠"` | 中文金属名 |
| `n_metal` | `1` | 金属阳离子数量 |
| `acid_salt` | `"hydrochloride"` | 英文酸盐后缀 |
| `acid_salt_zh` | `"盐酸盐"` | 中文酸盐后缀 |

这些元数据在 `_name_mol` 中先注入 `info["salt"]`（`src/namepredict/namer.py:228`，磷酸母体的盐门控与 `salt_meta` 来源），命名成功后再注入 `NameResult.meta["salt"]`（`src/namepredict/namer.py:231-232`），供 Layer5 组装器在生成中英双语名称时追加盐后缀（`_apply_salt_suffix`，`namer.py:192`）。由于盐解离在 Layer1 之前完成，Layer1-5 处理的始终是解离后的纯有机片段；除磷酸需按 `n_om` 与金属数配平外，各层代码无需关心盐的存在。

## 设计原则总结

NamePredict 的双语命名设计遵循以下原则：

1. **同步而非翻译**：英文和中文名称从同一结构化数据源（`numbered` 字典）同时生成，而非先生成一种语言再翻译为另一种。这避免了翻译过程中引入的二次歧义。

2. **跨层词表集中在 `constants.py`、词干派生集中在 `layer5/stems.py`**：前者是数值词/倍数前缀、单核氢化物、保留名、桥后缀、环外主基后缀等双语词表的唯一存放处；后者在其上派生烷烃词干与碳数词形（`alkane_en(n)` / `alkane_zh(n)`、`_en_stem(n)` / `zh_stem()`、`zh_num()`），C11+ 数量词取自 `constants.en_num_term`。官能团名称由 chain_engine 用这些词干 + 后缀拼出，杜绝了不同模块各造一套翻译的可能性。

3. **成对函数约定**：英中双语始终以成对形式出现——表（`_ALKANE_EN_BASE` / `_ALKANE_ZH_BASE`）、函数（`alkane_en(n)` / `alkane_zh(n)`、`acid_to_anion_en` / `acid_to_anion_zh`）、元组（`(en, zh)`）。不存在"只有英文"或"只有中文"的命名代码路径。

4. **差异显式化**：中英文命名规则的差异（如酯的语序、盐的后缀位置）在对应模块中显式处理，而非隐藏在通用逻辑中。这使得特定语种的规则变更不会意外影响另一种语言。

5. **空值安全**：所有双语返回函数在任一语言无法生成时返回 `None`（而非部分成功）。客户端代码（如 `_names_for` 的派发链）在返回 `None` 时自动降级到备选策略，确保管线始终产出完整的中英双语名称。

---

## 相关页面

- [[architecture/overview]] — 6 层架构总览，展示双语数据在各层之间的流转
- [[architecture/layer5-name-assembly]] — Layer5 名称组装层的完整文档，包含 `(en, zh)` 流水线的实现细节
- [[architecture/layer0-preprocessor]] — Layer0 预处理层，盐元数据的双语字段定义
- [[reference/core-data-contracts]] — `NameResult` 及 `numbered` 字典的双语字段结构
- [[index]] — Wiki 首页
