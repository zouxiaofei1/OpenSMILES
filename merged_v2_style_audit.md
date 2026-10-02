# merged_v2 数据集命名风格核查报告

对象：`benchmarks/merged_v2_benchmark.json`（6679 条，6 个来源：`chebi` / `chebi_truefail` / `chebi_tautomer` / `pubchem_truefail` / `pubchem_tautomer` / `smiles_tiers`）。
上次运行 `dual_ok` 4088 条（61.2%）。本次核查识别出 **25 组数据风格冲突（涉及 314 个条目）** 与 **23 组可由 IUPAC 规则唯一解释的差异（涉及 145 个条目）**，合计涉及 451 个条目。

裁决基准：风格统一优先，默认向 **PIN**（Preferred IUPAC Name）对齐；个别维度按实测改动面择优（见 1.2）。规则出处为 `docs/iupac/cn_translated/`（IUPAC Bluebook 中文译本）与 `docs/iupac/cn/`。

---

## 一、数据风格冲突（同一条目特征两种写法并存，namer 无法同时通过）

### 1.1 指示氢 `1H-` / `7H-` 缺写（最大宗，94 个条目）

PIN 规则：P-58.2.1.1「在优选 IUPAC 名中，当相应结构中存在指示氢时，必须**始终标出**」；P-25 表 2.8 明列 `indole` PIN = `1H-indole`、`purine` PIN = `7H-purine`；P-12 表 1.2 例 9 给出 `(1H-indol-1-yl)acetic acid` (PIN) —— 即使自由价落在带氢的 N 上，`1H-` 仍须写出。P-64 例 `di(1H-imidazol-1-yl)methanethione` (PIN) 同理。

**判定口径**：只有当该环 N 上仍带氢时才须写指示氢；N 已被取代或作为稠合组分时不须写（这部分属二、规则性差异，不是冲突）。以下条目经 rdkit 复核，环内确有游离 N-H 却省略了指示氢。

- **carbazole / pyrrole 型（N-H 存在但未标）**
  `pubchem-20760233`、`pubchem-22846941`、`pubchem-1829444`
  对照正确写法：`truefail-0317`「5,11-dimethyl-6H-pyrido[4,3-b]carbazole」、`chebi-3216`「9H-pyrido[3,4-b]indole」、`chebi-919`「1H-pyrrole-2-carboxylic acid」。
- **嘌呤 / 吲哚 / 咪唑 / 苯并咪唑型**
  省略侧：`chebi-207`、`chebi-2574`、`chebi-587`、`truefail-0068`、`truefail-0229`、`truefail-0261`、`truefail-0608`
  写出侧：`chebi-2629`、`chebi-1204`、`chebi-1790`、`chebi-2934`、`chebi-2969`、`chebi-3041`、`chebi-1132`、`chebi-1751`、`chebi-484`、`truefail-0003`、`truefail-0189`、`truefail-0209`、`truefail-0323`
  另：`pubchem-45078109` 写 `3H-benzimidazol-5-yl`，非最低位次（应为 `1H-`）。
- **吲哚 / 唑类（片 3）**
  省略侧：`chebi-876`、`chebi-2829`、`tiers-119391`、`tiers-93716`、`pubchem-143388932`、`tiers-154405`
  写出侧：`chebi-1180`、`pubchem-49423148`、`tiers-155716`、`tiers-65164`、`chebi-404`
- **嘌呤 9-位取代型（48 条，ChEBI 内部并存）**
  `chebi-1028`、`chebi-1280`、`chebi-1299`、`chebi-1390`、`chebi-1416`、`chebi-1425`、`chebi-1453`、`chebi-1505`、`chebi-1530`、`chebi-1598`、`chebi-1624`、`chebi-1630`、`chebi-1712`、`chebi-1877`、`chebi-1917`、`chebi-2164`、`chebi-2169`、`chebi-2209`、`chebi-227`、`chebi-2279`、`chebi-2324`、`chebi-2339`、`chebi-2501`、`chebi-252`、`chebi-2534`、`chebi-2621`、`chebi-264`、`chebi-2733`、`chebi-279`、`chebi-2951`、`chebi-3061`、`chebi-3084`、`chebi-3107`、`chebi-37`、`chebi-453`、`chebi-546`、`chebi-819`、`pubchem-135618584`、`pubchem-139488604`、`truefail-0015`、`truefail-0146`、`truefail-0263`、`truefail-0268`、`truefail-0278`、`truefail-0293`、`truefail-0298`、`truefail-0590`、`truefail-0610`
- **四唑（C5-连接）**
  `chebi-2908` 写 `2H-tetrazol-5-yl`，`pubchem-1156852` 裸写 `tetrazol-5-yl`。四唑不在 P-25.7.1.3.1 允许省略清单内。

**⚠ 需人工裁定的子项**：`chebi-37`、`chebi-453`、`chebi-801`、`chebi-1299`、`chebi-1468`、`chebi-1470`、`truefail-0004`、`truefail-0125`、`truefail-0240`、`truefail-0287`、`truefail-0356`、`truefail-0512`、`truefail-0629` 用了 `1H-purin-9-yl` / `3H-purin-9-yl` / `1H-purin-6-one`。P-25 表 2.8 记 purine PIN 为 `7H-purine`，但 6-氨基（腺嘌呤型）与 6-氧代（次黄嘌呤型）的游离氢分别在 N7 与 N1，**不能一律改成 `1H-`**，需按结构逐个定指示氢位次。

---

### 1.2 N-取代基定位符：`N1/N2/N4` vs `N/N′` vs `4-N-`（三式并存）

**统一目标：撇号式 `N` / `N′`**（字符形态取 ASCII 撇号 `'`，U+0027）。

选择依据（实测）：全库三式计数 —— 撇号式 214 条、N 数字式 34 条、位次-N- 式 14 条；namer 现有输出已以撇号式为主（撇号式组 26 条中 24 条 namer 输出撇号式、该组 13 条已通过；数字式组 17 条 namer 输出无一条匹配、0 条通过）。取撇号式改动面最小、与 namer 现状一致。字符形态：数据集 151 条、namer 输出 106 条均用 ASCII `'`，全库无 `′`(U+2032)，故统一为 `'`。

IUPAC 说明：PIN 用数字式（P-16 译本 895-896 行、P-62.2.3 例 `N1-(2-aminoethyl)-N1,N2,N2-trimethylethane-1,2-diamine`），撇号式属**一般 IUPAC 名**（同处列出 `N-ethyl-N′-methylethane-1,2-diamine`）。本次按"风格统一优先"取撇号式，属有意选用的一般名风格，不违反规则。

**改法**：位次按从小到大排序，最低位次的 N 不带撇号，其后依次 `′`、`″`；同一 N 上的多个取代基仍写 `N,N-`。
例：`N1-ethyl-N2-methylethane-1,2-diamine` → `N-ethyl-N'-methylethane-1,2-diamine`；`N2,N2-dimethylpropane-1,2-diamine` → `N,N-dimethylpropane-1,2-diamine`；pyrimidine-2,4-diamine 的 `N4-(...)` + `N2-(...)` → `N-(...)` + `N'-(...)`。

- **保留（已是撇号式，26 条）**
  `pubchem-112885672`、`pubchem-117086128`、`pubchem-114944686`、`pubchem-65504564`、`pubchem-94518320`、`pubchem-62083986`、`pubchem-21816602`、`pubchem-75786591`、`pubchem-155414522`、`pubchem-109090906`、`pubchem-90957759`、`pubchem-46831782`、`pubchem-136606542`、`tiers-91385`、`pubchem-110214999`、`pubchem-16782759`、`pubchem-38282670`、`pubchem-114133123`、`pubchem-64830644`、`pubchem-119615229`、`pubchem-112948434`、`pubchem-90361091`、`pubchem-171597532`、`chebi-2072`、`chebi-992`、`chebi-3021`
- **改：N 数字式 → 撇号式（34 条）**
  `pubchem-80403558`、`pubchem-105857869`、`pubchem-103968155`、`pubchem-80769779`、`pubchem-29067523`、`pubchem-54776072`、`pubchem-108422691`、`pubchem-115346335`、`pubchem-130470015`、`pubchem-82735202`、`pubchem-103088699`、`pubchem-103360109`、`pubchem-114785757`、`pubchem-65342627`、`pubchem-67105263`、`pubchem-79860090`、`pubchem-80982295`、`pubchem-91071327`、`pubchem-64301719`、`pubchem-80829418`、`pubchem-13801875`、`pubchem-113334397`、`pubchem-54868854`、`pubchem-101039887`、`pubchem-62984033`、`pubchem-140664268`、`pubchem-69014081`、`pubchem-47071973`、`pubchem-112944029`、`pubchem-91270484`、`pubchem-19139821`、`pubchem-112912028`、`pubchem-80849416`、`pubchem-28921016`
- **改：倒序位次-N- 式 → 撇号式（14 条）**
  `truefail-0382`、`pubchem-118421203`、`pubchem-15961883`、`pubchem-106747792`、`pubchem-163390445`、`pubchem-172624190`、`truefail-0116`、`truefail-0127`、`truefail-0447`、`truefail-0230`、`truefail-0283`、`truefail-0404`、`pubchem-104868178`、`pubchem-175950225`

同源铁证：`pubchem_truefail` 内 `pubchem-80403558`「6-chloro-N4-(thiolan-2-ylmethyl)pyrimidine-2,4-diamine」与 `pubchem-112885672`「N-(2-methoxyethyl)-N'-[(4-methylphenyl)methyl]pyrimidine-2,4-diamine」同为 pyrimidine-2,4-diamine 母体，一用 `N4-` 一用 `N'`。

**例外（不属本条）**：`N,N′-dimethylurea`、`N,N′-dimethylguanidine` 这类**无位次**官能母体用撇号本就是 PIN。

---

### 1.3 并列取代基的圆括号（`X(Y)Z` vs `X-YZ`；`(X)Y` vs `XY`）

P-16.5.1.1/.2：圆括号用于分隔涉及不同结构要素的同类位次、并括起简单取代基前缀；P-16 例 `ethyl(methyl)(propyl)phosphane`、`(chloromethyl)silane`；P-62.2.3 例 `3-[methyl(phenyl)amino]phenol`；P-29 明列优选前缀 `(furan-2-yl)methyl`、`(thiophen-2-yl)methyl`，并把 `2-thienylmethyl` 列为非优选。

**(a) 同一碳/氮上两个取代基**
- 括号式（PIN）：`pubchem-92212504`、`pubchem-93313275`、`pubchem-94197071`、`pubchem-39594423`、`pubchem-34014695`、`pubchem-109989873`、`pubchem-47460037`、`pubchem-62984315`、`pubchem-111235761`
- 连字符式（应改）：`chebi-1070`、`chebi-788`、`chebi-661`、`pubchem-10034971`、`pubchem-8453702`、`pubchem-92810812`、`pubchem-12557327`、`pubchem-153141168`、`pubchem-29442476`、`pubchem-78818982`

**(b) 含位次的取代基构成复合前缀（`(furan-2-yl)ethyl` vs `furan-2-ylethyl`）**
- 括号式（PIN）：`pubchem-111355188`、`pubchem-172015788`、`pubchem-38212454`、`pubchem-69324583`
- 裸写（应改）：`pubchem-111008773`、`pubchem-111304885`、`pubchem-16907213`、`pubchem-124295005`、`pubchem-169248170`、`pubchem-28122539`、`pubchem-39628531`、`pubchem-41442795`、`pubchem-52787841`、`pubchem-55754891`、`pubchem-72374086`、`pubchem-92584871`、`pubchem-100095121`、`pubchem-111804822`、`pubchem-112155309`、`tiers-102344`、`tiers-18034`、`truefail-0038`、`truefail-0516`

**(c) 同一 P 上并列两个简单含氧前缀**
P-67 酰基前缀表对并列在 P 上的简单前缀一律加圆括号：`hydroxy(sulfanyl)phosphoryl`、`dimethoxy(selenophosphoryl)`、`(hydroperoxy)phosphoryl`；PIN 例 `4,4′-(hydroxyphosphoryl)dibenzoic acid`。
- 括号式（PIN）：`chebi-2256`、`chebi-1005`、`truefail-0524`
- 连字符式（应改）：`chebi-546`、`chebi-122`、`chebi-279`、`chebi-819`、`chebi-1416`、`chebi-1877`、`chebi-1917`、`chebi-1983`、`truefail-0146`、`truefail-0278`、`truefail-0298`、`truefail-0330`、`truefail-0675` 等（该组共 31 条，另含 `chebi-264`、`chebi-679`、`chebi-1028`、`chebi-1144`、`chebi-1425`、`chebi-1505`、`chebi-1712`、`chebi-1965`、`chebi-1997`、`chebi-2169`、`chebi-2209`、`chebi-2324`、`chebi-2339`、`chebi-2501`、`chebi-2994`、`chebi-3107`）

**(d) 乘数前缀后多余圆括号**
`chebi-1282`「2-hydroxyethyl**tri(methyl)**azanium」。methyl 是简单前缀，P-16.5.1.1 的括号只对复合/复杂前缀；同源 64 条写 `trimethylazanium`（如 `chebi-1062`、`chebi-1677`）。

---

### 1.4 电中性 N-H 氨基被写成阳离子式的 `azaniumyl`

P-103.2.4.1 / P-74.1.3：氨基酸类按中性形式命名（`2-aminopropanoic acid`，而非 `2-azaniumylpropanoate`）；`azaniumyl` 仅对应 `-NH3+`（P-64）。

rdkit 实测（N 形式电荷 = 0、分子净电荷 = −1、N 上 H 数 = 2）：

| 条目 | 名称用词 | 实测 |
|---|---|---|
| `chebi-122` | `azaniumyl` | N⁺=0, net=−1 |
| `chebi-2781` | `azaniumyl` | N⁺=0, net=−1 |
| `chebi-2421` | `azaniumyl` | N⁺=0, net=−1 |
| `chebi-182` | `amino` | N⁺=0, net=−1 |
| `chebi-1844` | `amino` | N⁺=0, net=−1 |

同一电子结构两种写法并存：
- 写 `azaniumyl` 侧：`chebi-112`、`chebi-122`、`chebi-1997`、`chebi-2380`、`chebi-2421`、`chebi-2653`、`chebi-2871`、`chebi-2940`、`chebi-316`、`chebi-67`、`chebi-679`、`chebi-730`、`chebi-2781`、`chebi-829`、`chebi-1178`、`chebi-445`、`chebi-1048`、`chebi-1718`、`chebi-2141`、`chebi-2877`、`truefail-0025`、`truefail-0045`、`truefail-0055`、`truefail-0065`、`truefail-0096`、`truefail-0146`、`truefail-0278`、`truefail-0308`、`truefail-0330`、`truefail-0415`、`truefail-0416`、`truefail-0427`、`truefail-0620`、`truefail-0675`、`truefail-0707`、`truefail-0712`
- 写 `amino` 侧：`chebi-182`、`chebi-1844`、`chebi-2413`、`chebi-616`、`chebi-2462`

→ 统一为 `amino`。

---

### 1.5 外层嵌套封闭符号用 `()` 而非 `[]`

P-16.5.2.4：方括号括起**其中已用圆括号**的取代基前缀。例 `4-[(hydroxyselanyl)methyl]benzoic acid` (PIN)。

- 错（圆括号套圆括号）：`tiers-100116`「1-((1-(1,3-benzodioxol-5-yl)-5-oxopyrrolidin-3-yl)methyl)-…」、`tiers-40911`、`tiers-125461`、`tiers-98594`
- 对：`tiers-147895`「1-[(4-fluorophenyl)methyl]-2-methylsulfonylbenzimidazole」、`tiers-54438`「1-chloro-3-[chloro(phenyl)methyl]benzene」

同为 `smiles_tiers` 来源内部并存；全库嵌套圆括号仅 34 条。

---

### 1.6 取代基闭括号与后续片段之间多写一个连字符

P-16.2.4：连字符只用于位次与字母之间、以及后缀/母体连接；`(X)` 闭括号后直接接母体，不加连字符。

- 错：`chebi-657`「5-(trifluoromethyl)**-**pyridin-2-yl」、`chebi-661`、`chebi-265`「(3,4-dimethoxy**-**pyridin-2-yl)」、`chebi-1896`、`truefail-0157`、`truefail-0691`、`pubchem-32141430`、`pubchem-54992377`、`pubchem-27677160`
- 对：同源大量条目写 `(3,4-dimethoxypyridin-2-yl)`、`(trifluoromethyl)pyridin-2-yl`

---

### 1.7 杂单环骨架位次被完全省略

P-14.3.4.4：只有当移动或互换位次不产生异构体时才允许省略。这些环存在 1,2-/1,3-/1,4- 异构体，位次不可省。

- 裸写（应补位次）：`pubchem-154825743`「diazinan-4-yl」、`pubchem-163926779`「oxathiane」、`pubchem-70082595`、`pubchem-137258909`「oxadiazol-3-ium」
- 对照正确写法：`pubchem-86051096`「1,3-diazinan-1-yl」、`truefail-0129`「1,3-oxathiane」、`tiers-131305`「1,3-oxazinan-2-yl」

同属 `pubchem_truefail` 来源内部并存。

---

### 1.8 环烯基自由价碳上 `en` 的位次被省略

省去后 `cyclohexen-1-yl` 与 `cyclohex-2-en-1-yl` 无法区分，不满足 P-14.3.4.4；P-14.3.3 要求 PIN 引用必需位次。

- 省略：`chebi-1285`「(5Z)-7-[…-5-oxo**cyclopenten-1-yl**]hept-5-enoate」、`chebi-2489`、`pubchem-11131210`、`pubchem-154232309`、`truefail-0020`
- 全位次（PIN）：`chebi-823`「…-5-oxo**cyclopent-2-en-1-yl**…」（与 `chebi-1285` 同骨架）、`chebi-2376`、`pubchem-45179479`、`pubchem-9635845`、`truefail-0141`
- 另有 `chebi-912`「2,6,6-trimethylcyclohexen-1-yl」

ChEBI 同一来源内部并存。

---

### 1.9 酰氧基中 `oxy` 的归属：整体式 vs 断开式（34 条）

- A 整体式 `[(9Z)-octadec-9-enoyloxy]`（oxy 在方括号内），全库 151 条
- B 断开式 `[(9Z)-octadec-9-enoyl]oxy`（方括号在 `oyl` 后闭合），全库 64 条

Blue Book 两式皆有例（P-65 例 `2-(acetyloxy)-3-(hexadecanoyloxy)propyl` 为整体式；P-67 例 `3-[(chlorosulfonyl)oxy]propanoic acid` 为断开式），故属**需统一**的风格问题，建议取多数派整体式 A。

**同一名称内两式并现的铁证**：`chebi-1544`（`…hexadecanoyloxytetradecanoyl…` 与 `…tetradecanoyloxytetradecanoyl]oxy` 同现）、`chebi-2751`、`chebi-1876`（`2-acetyloxy-3-[(9Z)-hexadec-9-enoyl]oxypropyl`）。

全部 34 条：`chebi-21`、`chebi-49`、`chebi-215`、`chebi-238`、`chebi-429`、`chebi-600`、`chebi-719`、`chebi-831`、`chebi-1013`、`chebi-1133`、`chebi-1475`、`chebi-1544`、`chebi-1716`、`chebi-1852`、`chebi-1876`、`chebi-1926`、`chebi-1987`、`chebi-2013`、`chebi-2117`、`chebi-2172`、`chebi-2492`、`chebi-2696`、`chebi-2751`、`chebi-2793`、`chebi-2907`、`chebi-3038`、`chebi-3088`、`chebi-3196`、`chebi-3270`、`truefail-0094`、`truefail-0501`、`truefail-0679`、`pubchem-150594676`、`tiers-10669`

---

### 1.10 前缀次序与低频孤例

| 子项 | 条目 | 说明 |
|---|---|---|
| CoA 型连接臂前缀次序（P-14.5.1 字母序 e<o） | `chebi-85`、`chebi-781`、`chebi-3000`（`3-oxo-3-[…乙氨基]propyl`，应改）vs `chebi-499`、`chebi-1091`、`chebi-1214`、`chebi-1370`、`chebi-1519`、`chebi-1597`、`chebi-1950`、`chebi-2052`、`chebi-2121`、`chebi-2288`、`chebi-2529`、`chebi-2682`、`chebi-2861`、`chebi-3058`、`chebi-3117`（`3-[…]-3-oxopropyl`） | 同一 chebi 来源内同骨架两式并存 |
| 甲苯磺酰基用词 | `pubchem-277242`「2-methylbenzenesulfonyl fluoride」vs `tiers-116637`「1-(4-methylphenyl)sulfonylindole-3-carboxylate」 | 建议统一为括号式；全库未取代情形已统一（`phenylsulfonyl` 0 条 / `benzenesulfonyl` 12 条） |
| 苯并噻吩稠合位次 | `tiers-19300`「…5,7-dihydro-2-benzothiophene-1-carboxylate」，全库孤例；对照 `tiers-35463`、`pubchem-20175335`、`pubchem-40114806` 均为 `1-benzothiophene` | 需人工确认 S 位置 |
| 英文名尾随空格 | `chebi-963`：`'(carbamoylcarbamoyl)carbamic acid '` | 仅此一条 |

---

## 二、可由 IUPAC 规则唯一解释的差异（数据自洽，属 namer 侧偏差）

这些不是数据冲突，**不要改数据**；但可作为 namer 待修清单。

| 维度 | 特征 | 规则 | 代表条目 |
|---|---|---|---|
| 位次 | 二羧酸/二酰胺按 P-14.3.4.1 省末端位次；二胺/二醇两位次不同须保留 | P-14.3.4.1/.2 | 省：`chebi-1774`、`chebi-1230`、`truefail-0632`；留：`chebi-2568`、`chebi-617`、`chebi-2378`、`tiers-173`、`tiers-55`、`tiers-98`、`tiers-174`、`chebi-2913`、`chebi-1970`、`chebi-1303`、`chebi-1660`、`pubchem-115126140`、`pubchem-28255181`、`pubchem-107058380` |
| 用词 | 母体末尾 `e` 省略（元音/`y` 前省 e，辅音后缀前留 e） | P-14.4 | `chebi-1949`、`chebi-950`、`chebi-1067`、`chebi-3152`、`chebi-919`、`tiers-138602`、`chebi-1126`、`chebi-1971`、`chebi-1550`、`chebi-988`、`tiers-35976`、`tiers-7653`、`chebi-1331`、`chebi-2413` |
| 用词 | `carboxylate`（后缀/阴离子）vs `carboxylato-`（取代基前缀）vs `carboxymethyl`（中性） | P-65 / P-66 / P-29.1 | `chebi-176`、`chebi-1980`、`chebi-2434`、`chebi-2845`、`chebi-2351`、`chebi-3056`、`chebi-3275`、`chebi-2781`、`chebi-2877`、`chebi-1066`、`chebi-1857`、`chebi-869`、`chebi-1376`、`truefail-0062`、`truefail-0133`、`truefail-0305`、`truefail-0617`、`truefail-0210`、`truefail-0165`、`truefail-0421`、`truefail-0666` |
| 用词 | `amino`(-NH₂) vs `azaniumyl`(-NH₃⁺) | P-62.2 / P-62.5 / P-66.4 | `chebi-103`、`chebi-1436`、`chebi-1686`、`chebi-1380`、`chebi-121`、`chebi-1842`、`chebi-1852`、`chebi-1629`、`truefail-0517` |
| 用词 | `thione` 后缀 vs `sulfanylidene` 前缀（C=O 优先于 C=S） | P-64.6.1/.2 / P-66.5 | `tiers-19564`、`truefail-0392`、`tiers-73757`、`tiers-148854`、`tiers-25538`、`tiers-60725`、`truefail-0625`、`chebi-796`、`chebi-2427`、`chebi-3200`、`chebi-409`、`chebi-2950`、`truefail-0104`、`chebi-972`、`chebi-2856`、`chebi-3074`、`pubchem-139488604`、`pubchem-14450841`、`pubchem-165776602`、`pubchem-167248404`、`pubchem-375375` |
| 用词 | 磺酸系 `sulfo`(-SO₃H) / `sulfonato`(-SO₃⁻) / `sulfonyl`(-SO₂-) / `sulfooxy` | P-67.3 / P-64 | `chebi-106`、`chebi-804`、`chebi-1982`、`chebi-247`、`chebi-687`、`chebi-139`、`tiers-92413`、`pubchem-72520845`、`tiers-152730`、`pubchem-106001362` |
| 用词 | `sulfanyl` 为优选前缀；`thio` 仅存于保留名 | P-63.1.5 | `chebi-2156`、`chebi-2072`、`chebi-531`、`chebi-2006`、`pubchem-106001362`、`truefail-0591`、`truefail-0596` |
| 用词 | 五元氧杂环统一 `oxolane`/`oxolan-2-yl`；全库无 `tetrahydrofuran`/`thien-`/`furyl` | P-22.2.4 / P-25.2.1 | `chebi-592`、`chebi-1180`、`chebi-1067`、`chebi-3152`、`chebi-30`、`chebi-2152`、`tiers-13374`、`truefail-0024`、`truefail-0069`、`truefail-0549`、`truefail-0684` |
| 用词 | 高价态 `lambda5-phosphane` / `lambda4-sulfanyl` 写法（数据一致，namer 该类 32 条全败） | P-14.3.4.2 | `chebi-1179`、`chebi-1812`、`tiers-59592`、`truefail-0188`、`truefail-0198`、`pubchem-142954085`、`pubchem-155553813` |
| 用词 | 阴离子母体/后缀（`-olate`、`-ate`、`oxido`、`sulfonate`）按 SMILES 电荷态命名 | P-66.3 | `chebi-25`、`chebi-2001`、`chebi-581`、`chebi-186` |
| 括号 | 稠环内杂原子位次用方括号 `[1,3]benzothiazol-`；独立母体写 `1,3-benzothiazol-` | P-16.5.2.2 / P-25.3.5 | `chebi-2115`、`chebi-25`、`tiers-152873`、`tiers-146758`、`tiers-140247`、`tiers-96177`、`pubchem-130428842`、`tiers-108673`、`tiers-94714`、`truefail-0309`、`truefail-0546`、`pubchem-18399665` |
| 指示氢 | `1H-` 有无由取代/稠合上下文决定 | P-58.2.1 / P-14.7.1 | `chebi-190`、`chebi-2048`、`chebi-2151`、`tiers-25570`、`pubchem-128446264`、`pubchem-129768153` |

---

## 三、统一方向（按影响面排序）

1. **指示氢一律按结构补写**（约 94 条）——对齐 PIN，`1H-indole` / `7H-purine` 为准；嘌呤的指示氢位次需按 6-取代类型分别判定。
2. **N 定位符三式统一为撇号式 `N/N′`**（48 条已改：N 数字式 34 + 位次-N- 式 14，另 26 条已是撇号式保留）——字符统一为 ASCII `'`；转换时最低位次 N 不加撇号。此维度按改动面择优，采用一般 IUPAC 名风格而非 PIN 数字式（见 1.2 与第四节）。
3. **并列取代基一律加圆括号**（约 70 条）——包括 `X(Y)Z`、`(furan-2-yl)ethyl`、P 上的 `hydroxy(oxido)phosphoryl`。
4. **电中性 N-H 一律写 `amino`**，`azaniumyl` 只留给 `-NH₃⁺`（约 40 条）。
5. **外层嵌套封闭符号用 `[]`**（4 条）、**删去闭括号后多余连字符**（9 条）。
6. **补全异构体必需的骨架位次**：杂单环骨架（4 条）、环烯基 `en` 位次（11 条）。
7. **酰氧基 `oxy` 归属取整体式**（34 条，多数派 151:64）。

统一后的风格基线：除第 2 项按实测择优取撇号式外，其余向 PIN 对齐。`pubchem_truefail` 是第 2、3、4 类偏差的主要来源，`chebi` 是第 1、9 类的主要来源，`smiles_tiers` 是第 5 类的主要来源。

---

## 四、本次已应用的改动（1.2 N 定位符）

改动范围：仅 `benchmarks/merged_v2_benchmark.json` 的 `english_name` 字段，共 **48 条**；`smiles`、`chinese_name` 及其余字段未动。备份在 `tmp/v2audit/merged_v2_benchmark.orig.json`。

改法：N 数字式 34 条 + 位次-N- 式 14 条 → 撇号式，按 P-16.9.2「最低位次 N 不加撇号，其后依次 ′、″」赋值并保留原有前缀书写位置。
全库计数：撇号式 214 → 237，N 数字式 34 → 0，位次-N- 式 14 → 0。

**结果：`dual_ok` 4088/6679（61.21%）→ 4112/6679（61.57%），+24，0 条回退。** 48 条中 24 条转为通过。

### 剩余 24 条：属 namer 侧问题，数据现为 IUPAC 正解

**（a）撇号归属相反——namer 按"名称中先出现的 N"给无撇号，IUPAC（P-16.9.2）按最低位次**
P-16 译本 895-896 行成对给出 `N1-ethyl-N2-methylethane-1,2-diamine` (PIN) 与 `N-ethyl-N′-methylethane-1,2-diamine` —— 无撇号者对应 **N1（最低位次）**。数据已按此写入，namer 未实现：

| 条目 | 数据（IUPAC） | namer 输出 |
|---|---|---|
| `pubchem-105857869` | `N'-[(1-methoxycyclobutyl)methyl]-…-N-propyl…` | `N-[(1-methoxycyclobutyl)methyl]-…-N'-propyl…` |
| `pubchem-80769779` | `N'-cyclopentyl-N-ethyl…` | `N-cyclopentyl-N'-ethyl…` |
| `pubchem-115346335` | `N,N,N'-trimethyl…` | `N,N',N'-trimethyl…` |
| `pubchem-130470015` | `N'-ethyl-N,N-dimethyl…` | `N-ethyl-N',N'-dimethyl…` |
| `pubchem-82735202` | `N'-ethyl-N,3-dimethyl-N-…` | `N-ethyl-N',3-dimethyl-N'-…` |
| `pubchem-79860090` | `N',3,3-trimethyl-N-…-N-…` | `N,3,3-trimethyl-N'-…-N'-…` |
| `truefail-0382` | `N'-(2-fluorophenyl)-N-(4-methoxyphenyl)…` | `N-(2-fluorophenyl)-N'-(4-methoxyphenyl)…` |
| `truefail-0127` | `N'-(2-cyanopropan-2-yl)-N-[…]…` | `N-(2-cyanopropan-2-yl)-N'-[…]…` |
| `pubchem-80829418` | `N'-ethyl-N-(2-methylcyclopentyl)…` | `N-ethyl-N'-(2-methylcyclopentyl)…` |
| `pubchem-54868854` | `N',N',2-trimethyl-N-[…]…` | `N,N,2-trimethyl-N'-[…]…` |
| `pubchem-62984033` | `N'-[1-(4-chlorophenyl)ethyl]-N-ethyl-N'-methyl…` | `N-[1-(4-chlorophenyl)ethyl]-N'-ethyl-N-methyl…` |
| `pubchem-112944029` | `N'-(3,4-difluorophenyl)-N-(3-methoxypropyl)…` | 反向 |
| `pubchem-91270484` | `N'-(4-ethylphenyl)-N-hydroxy…` | 反向 |
| `pubchem-112912028` | `N'-(2,3-dihydro-1,4-benzodioxin-6-yl)-N-(3-methoxypropyl)…` | 反向 |
| `pubchem-80849416` | `N',6-dimethyl-N-propyl-N'-(2-thiophen-2-ylethyl)…` | `N,6-dimethyl-N'-propyl-N-(2-thiophen-2-ylethyl)…` |
| `pubchem-19139821` | `N'-benzhydryl-N-butyl-N-methyl…` | `N-butyl-N'-(diphenylmethyl)-N-methyl…`（撇号 + `benzhydryl`/`diphenylmethyl` 用词） |

**（b）其余与 N 定位符无关的独立差异（数据本身另有其它维度问题，见第一节）**

| 条目 | 差异 | 归属维度 |
|---|---|---|
| `pubchem-118421203` | namer 未写出 `N-oxo`（对应 SMILES `C(=O)N=O`） | namer 缺口 |
| `pubchem-15961883` | `[3-[[…]oxy]propyl]` vs `[3-[…]oxypropyl]` | 括号 P-16.5.2.4 |
| `pubchem-106747792` | `N-(bicyclo[2.2.1]heptan-2-yl)` vs `N-bicyclo[2.2.1]heptan-2-yl`（namer 缺圆括号） | 括号 P-16.5.1.2 |
| `pubchem-47071973` | `bis(1H-pyrazol-4-yl)` vs `di(1H-pyrazol-4-yl)`（数据用 `bis` 正确） | P-16.3 |
| `pubchem-140664268` | `cyclohex-4-ene-1,2-dicarboxamide` vs namer `cyclohexane-…`（namer 丢双键） | namer 缺口 |
| `pubchem-13801875` | 数据取代式 vs namer 官能类式 `3,5-bis(tert-butylamino)-…` | 命名类型选择 |
| `truefail-0447` | 稠环名不同（furo[3,2-d][1,3,2]dioxaphosphinine vs phosphabicyclo[4.3.0]nonane） | 稠环命名 |
| `truefail-0230` | 整名构造完全不同（母体选择不同） | 母体选择 |
