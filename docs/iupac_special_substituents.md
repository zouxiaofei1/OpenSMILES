# IUPAC 2013 Blue Book: "地位特殊的取代基" 全面名单

> 依据：IUPAC Blue Book (2013) P-29, P-31, P-32, P-57, P-65；
> 项目现状：`src/namepredict/layer2/`–`layer3/` retained naming。
>
> **重要**：2013 Blue Book 大幅缩减了保留名。以下用三种级别标注：
> - ✅ **PIN** = Preferred IUPAC Name（首选名）
> - ⚠️ **General only** = 仅限普通命名，不能出现在 PIN 中
> - ❌ **Not recommended** = 2013 已废弃（P-57.1.4）

---

## 1. 支链烷基取代基 (P-29.3.1 / P-29.6 / P-57.1)

### 1a. PIN 级别的保留名（仅一个！）

| # | 保留名 | 系统名 | 中文 | 项目状态 |
|---|--------|--------|------|---------|
| 1 | **tert-butyl** | 1,1-dimethylethyl | 叔丁基 | ✅ 已实现 |

> Note: `tert-` italic 前缀不参与字母排序（P-14.5）。

### 1b. 普通命名级别（非 PIN，但 IUPAC 允许在非 PIN 中使用）

| # | 保留名 | PIN (必须用这个) | 中文 | 项目状态 |
|---|--------|-----------------|------|---------|
| 2 | **isopropyl** | propan-2-yl | 异丙基 | ✅ 已实现 |
| 3 | **isopropylidene** | propan-2-ylidene | 异亚丙基 | ❌ 未实现（二价） |

> 仅 `isopropyl` 是唯一幸存的 iso- 前缀。其余 iso- 全部废除。

### 1c. 2013 年已废弃（P-57.1.4 Not Recommended）

| # | 废弃名 | PIN（必须用这个） | 中文 | 项目现状 |
|---|--------|-----------------|------|---------|
| 4 | **isobutyl** | 2-methylpropyl | 异丁基 | ⚠️ 项目仍用保留名 |
| 5 | **sec-butyl** | butan-2-yl | 仲丁基 | ⚠️ 项目仍用保留名 |
| 6 | **isopentyl** | 3-methylbutyl | 异戊基 | ⚠️ 项目仍用保留名 |
| 7 | **neopentyl** | 2,2-dimethylpropyl | 新戊基 | ⚠️ 项目仍用保留名 |
| 8 | **tert-pentyl** | 2-methylbutan-2-yl | 叔戊基 | ✅ 项目用系统名 |
| 9 | **phenethyl** | 2-phenylethyl | 2-苯乙基 | ❌ 未实现 |
| 10 | **benzhydryl** | diphenylmethyl | 二苯甲基 | ❌ 未实现 |
| 11 | **trityl** | triphenylmethyl | 三苯甲基 | ❌ 未实现（注：trityl 归入普通命名） |

---

## 2. 不饱和无环烃基取代基 (P-29.3.2 / P-29.6.1)

### 2a. 普通命名级别（非 PIN）

| # | 保留名 | PIN (必须用这个) | 中文 | 项目状态 |
|---|--------|-----------------|------|---------|
| 12 | **vinyl** | ethenyl | 乙烯基 | ❌ **缺失** |
| 13 | **allyl** | prop-2-en-1-yl | 烯丙基 | ❌ **缺失** |
| 14 | **isopropenyl** | prop-1-en-2-yl | 异丙烯基 | ❌ **缺失** |
| 15 | **vinylidene** | ethenylidene | 亚乙烯基 | ❌ 未实现 |
| 16 | **allylidene** | prop-2-en-1-ylidene | 亚烯丙基 | ❌ 未实现 |
| 17 | **allylidyne** | prop-2-en-1-ylidyne | 次烯丙基 | ❌ 未实现 |

### 2b. 不再保留（2013 确认非 PIN）

| # | 非PIN名 | PIN (必须用这个) | 中文 | 项目状态 |
|---|---------|-----------------|------|---------|
| 18 | propargyl | prop-2-yn-1-yl | 炔丙基 | ❌ **缺失**（propargyl 非 PIN） |
| 19 | prenyl | 3-methylbut-2-enyl | 异戊烯基 | ✅ 项目用系统名 |
| 20 | crotyl | (E)-but-2-enyl | 巴豆基 | ❌ 未实现 |
| 21 | cinnamyl | 3-phenylprop-2-en-1-yl | 肉桂基 | ❌ 未实现 |

**关键发现**：IUPAC 2013 将 vinyl/allyl/isopropenyl 降级为"普通命名"(P-29.6.2.2 P-57.1.3)，**不允许在 PIN 中使用**。但它们在普通命名中仍极为常用。`ethenyl`、`prop-2-en-1-yl`、`prop-1-en-2-yl` 才是 PIN。

---

## 3. 环状取代基 (P-29.6 / P-52.2.8)

### 3a. 环烷基 — 系统名，非保留名

环烷基（cyclopropyl ~ cyclooctyl 等）**并非保留名**，是从母体氢化物按标准后缀规则系统派生的。它们完全属于 PIN。

| # | 名称 | 中文 | 项目状态 |
|---|------|------|---------|
| 22 | cyclopropyl | 环丙基 | ✅ C3–C8 已实现 |
| 23 | cyclobutyl | 环丁基 | ✅ |
| 24 | cyclopentyl | 环戊基 | ✅ |
| 25 | cyclohexyl | 环己基 | ✅ |
| 26 | cycloheptyl | 环庚基 | ✅ |
| 27 | cyclooctyl | 环辛基 | ✅ |
| 28 | cyclononyl | 环壬基 | ❌ (C9+) |
| 29 | cyclodecyl | 环癸基 | ❌ (C10+) |

### 3b. 饱和杂环侧链

| # | 名称 | 中文 | 项目状态 |
|---|------|------|---------|
| 30 | pyrrolidin-n-yl | 吡咯烷-n-基 | ✅ |
| 31 | piperidin-n-yl | 哌啶-n-基 | ✅ |
| 32 | oxan-n-yl | 噁烷-n-基 | ✅ |
| 33 | oxolan-n-yl | 氧杂环戊烷-n-基 | ✅ |
| 34 | morpholin-n-yl | 吗啉-n-基 | ✅ |

### 3c. 1-环烷基乙基

| # | 名称 | 项目状态 |
|---|------|---------|
| 35 | 1-cyclohexylethyl | ✅ 已实现（C6；C3–C8 可扩展） |

---

## 4. 芳基与芳氧基取代基 (P-29.3.2 / P-29.6 / P-57.1.2)

### 4a. PIN 级别保留名

| # | 保留名 | 系统名 | 中文 | 项目状态 |
|---|--------|--------|------|---------|
| 36 | **phenyl** | — | 苯基 | ✅ 已实现（含取代变体） |
| 37 | **benzyl** | phenylmethyl | 苄基 | ✅ 已实现（P-57.1.2：不可被取代） |
| 38 | **benzylidene** | phenylmethylidene | 亚苄基 | ❌ 未实现 |
| 39 | **benzylidyne** | phenylmethylidyne | 次苄基 | ❌ 未实现 |

### 4b. 系统名（非"保留"但常用为特殊形态）

| # | 名称 | 中文 | 项目状态 |
|---|------|------|---------|
| 40 | phenoxy | 苯氧基 | ✅ 已实现 |
| 41 | benzyloxy | 苄氧基 | ✅ 已实现 |
| 42 | naphthalen-1-yl | 萘-1-基 | ✅ 已实现 |
| 43 | naphthalen-2-yl | 萘-2-基 | ✅ 已实现 |
| 44 | pyridin-n-yl | 吡啶-n-基 | ✅ 已实现 |
| 45 | anthracen-n-yl | 蒽-n-基 | ⚠️ 母体已实现 |
| 46 | phenanthren-n-yl | 菲-n-基 | ❌ |
| 47 | biphenyl-n-yl | 联苯-n-基 | ❌ |

### 4c. 杂芳基 (P-57.1.5.3 — 全部为 General Only，PIN 必须用系统名)

| # | General 名 | PIN (必须用这个) | 中文 |
|---|-----------|-----------------|------|
| 47a | 2-furyl / 3-furyl | furan-2-yl / furan-3-yl | 呋喃基 |
| 47b | 2-thienyl / 3-thienyl | thiophen-2-yl / thiophen-3-yl | 噻吩基 |
| 47c | 2-pyridyl / 3-pyridyl / 4-pyridyl | pyridin-2-yl / pyridin-3-yl / pyridin-4-yl | 吡啶基 |
| 47d | 2-quinolyl (及异构体) | quinolin-n-yl | 喹啉基 |
| 47e | 1-/3-/4-/5-/6-/7-/8-isoquinolyl | isoquinolin-n-yl | 异喹啉基 |
| 47f | 2-anthryl / 9-anthryl | anthracen-2-yl / anthracen-9-yl | 蒽基 |
| 47g | 9-phenanthryl | phenanthren-9-yl (及1,2,3,4) | 菲基 |
| 47h | 2-adamantyl | adamantan-2-yl | 金刚烷基 |
| 47i | 2-/3-/4-piperidyl | piperidin-2/3/4-yl | 哌啶基 |

### 4d. 已废弃的芳基/杂芳基名 (P-57.1.5.4)

| # | 废弃名 | PIN | 中文 |
|---|--------|-----|------|
| 47j | furfuryl | furan-2-ylmethyl | 糠基 |
| 47k | thenyl | thiophen-2-ylmethyl | 噻吩甲基 |
| 47l | o/m/p-tolyl | 2-/3-/4-methylphenyl (不可被取代) | 甲苯基 |

---

## 5. 二价/多价取代基 (P-29.2/P-29.5/P-57.1.1.2)

### 5a. PIN 级别保留名

| # | 保留名 | 系统名 | 中文 | 项目状态 |
|---|--------|--------|------|---------|
| 48 | **methylidene** (=CH2, 同原子) | methylene (非PIN) | 亚甲基 | ❌ P-71.2.2.1 |
| 49 | **methanediyl** (–CH2–, 桥连) | methylene (非PIN) | 亚甲基 | ❌ 二价桥连 PIN |
| 50 | **ethylene** | ethane-1,2-diyl | 1,2-亚乙基 | ❌ General only |
| 51 | **1,2-phenylene** | — | 1,2-亚苯基 | ❌ P-57.1.5.2 PIN |
| 52 | **1,3-phenylene** | — | 1,3-亚苯基 | ❌ P-57.1.5.2 PIN |
| 53 | **1,4-phenylene** | — | 1,4-亚苯基 | ❌ P-57.1.5.2 PIN |

> ⚠️ 2013 重要变化：`methylene` 不再推荐作为 =CH2 的 PIN，必须用 `methylidene` (P-71.2.2.1)。桥连 –CH2– 的 PIN 是 `methanediyl`。

### 5b. 普通命名级别

| # | 保留名 | PIN | 中文 | 项目状态 |
|---|--------|-----|------|---------|
| 53 | vinylene | ethene-1,2-diyl | 1,2-亚乙烯基 | ❌ |
| 54 | trimethylene | propane-1,3-diyl | 1,3-亚丙基 | ❌ |
| 55 | ethynylene | ethyne-1,2-diyl | 1,2-亚乙炔基 | ❌ |

### 5c. 二价功能基（系统命名，非保留）

| # | 名称 | 中文 | 项目状态 |
|---|------|------|---------|
| 56 | carbonyl | 羰基 (—CO—) | ❌ |
| 57 | sulfonyl | 磺酰基 (—SO2—) | ❌ |
| 58 | sulfinyl | 亚磺酰基 (—SO—) | ❌ |
| 59 | oxy | 氧基 (—O—) | ✅ 醚母体 |
| 60 | thio | 硫基 (—S—) | ✅ 硫醚母体 |
| 61 | dioxy | 二氧基 (—O—O—) | ❌ |
| 62 | azo | 偶氮基 (—N=N—) | ❌ |
| 63 | hydrazo | 1,2-亚肼基 | ❌ |
| 64 | carbonyldioxy | 羰基二氧基 | ❌ |
| 65 | oxydicarbonyl | 氧基二羰基 | ❌ |

---

## 6. 含杂原子取代基前缀 (P-31.1 / P-32 / P-63–P-66)

### 6a. 含氧前缀

| # | 名称 | 中文 | 类别 | 项目状态 |
|---|------|------|------|---------|
| 66 | hydroxy | 羟基 | PIN | ✅ |
| 67 | methoxy | 甲氧基 | PIN | ✅ |
| 68 | ethoxy | 乙氧基 | PIN | ✅ |
| 69 | propoxy | 丙氧基 | PIN | ✅ |
| 70 | isopropoxy | 异丙氧基 | PIN | ✅ |
| 71 | butoxy | 丁氧基 | PIN | ✅ |
| 72 | isobutoxy | 异丁氧基 | PIN | ✅ |
| 73 | hydroperoxy | 氢过氧基 | PIN | ❌ |
| 74 | peroxy (bridge) | 过氧基 | PIN | ❌ |

### 6b. 含硫前缀

| # | 名称 | 中文 | 类别 | 项目状态 |
|---|------|------|------|---------|
| 75 | sulfanyl | 硫烷基 | PIN | ⚠️ **mercapto 已废除** |
| 76 | methylsulfanyl | 甲硫基 | PIN | ✅ 已实现 |
| 77 | ethylsulfanyl | 乙硫基 | 系统 | ❌ |
| 78 | methylsulfinyl | 甲亚磺酰基 | 系统 | ❌ |
| 79 | methylsulfonyl | 甲磺酰基 | PIN | ❌ |
| 80 | mesyl | 甲磺酰基 | ❌非PIN | ❌ |
| 81 | tosyl | 对甲苯磺酰基 | ❌非PIN | ❌ |
| 82 | triflyl | 三氟甲磺酰基 | ❌非PIN | ❌ |

### 6c. 含氮前缀

| # | 名称 | 中文 | 类别 | 项目状态 |
|---|------|------|------|---------|
| 83 | amino | 氨基 | PIN | ✅ |
| 84 | methylamino | 甲氨基 | 系统 | ❌ |
| 85 | dimethylamino | 二甲氨基 | 系统 | ❌ |
| 86 | nitro | 硝基 | PIN | ✅ |
| 87 | nitroso | 亚硝基 | PIN | ❌ |
| 88 | cyano | 氰基 | PIN | ✅ 芳环叶 |
| 89 | isocyano | 异氰基 | PIN | ❌ |
| 90 | isocyanato | 异氰酸根合 | PIN | ✅ |
| 91 | isothiocyanato | 异硫氰酸根合 | PIN | ✅ |
| 92 | azido | 叠氮基 | PIN | ❌ |
| 93 | diazenyl | 二氮烯基 (HN=N–) | PIN | ❌ (azo 旧式非PIN) |
| 94 | diazo | 重氮基 (=N2) | PIN | ❌ |
| 95 | **hydrazinyl** | 肼基 (NH2–NH–) | ✅ PIN (P-29 Type 2a) | ❌ 注意: 不是 hydrazino |
| 95a | **phosphanyl** | 膦基 (H2P–) | PIN | ❌ phosphino 已废除 |
| 95b | **boranyl** | 硼烷基 (H2B–) | PIN | ❌ boryl 已废除 |
| 95c | **sulfanyl** (重复) | 硫烷基 (HS–) | PIN | ⚠️ mercapto 已废除 |
| 95d | **selanyl** | 硒烷基 (HSe–) | PIN | ❌ selenyl 已废除 |
| 96 | amidino | 脒基 | ❌ **已废除 (P-2013)** | ❌ PIN=carbaminidoyl |
| 97 | guanidino | 胍基 | ❌ **已废除 (P-2013)** | ❌ PIN=carbamimidamido |
| 98 | ureido | 脲基 | ❌ **已废除 (P-2013)** | ❌ PIN=carbamoylamino |
| 98a | **anilino** | 苯胺基 | ✅ PIN (P-62.2.1.1.1) | ❌ |

### 6d. 卤素前缀（全部为PIN系统前缀）

| # | 名称 | 中文 | 项目状态 |
|---|------|------|---------|
| 98 | fluoro | 氟 | ✅ |
| 99 | chloro | 氯 | ✅ |
| 100 | bromo | 溴 | ✅ |
| 101 | iodo | 碘 | ✅ |
| 102 | trifluoromethyl | 三氟甲基 | ✅ |
| 103 | trichloromethyl | 三氯甲基 | ❌ |
| 104 | tribromomethyl | 三溴甲基 | ❌ |
| 105 | difluoromethyl | 二氟甲基 | ❌ |
| 106 | pentafluoroethyl | 五氟乙基 | ❌ |

### 6e. 酰基/羰基前缀（P-65/P-66）

| # | PIN 名 | 非PIN保留名 | 中文 | 项目状态 |
|---|--------|------------|------|---------|
| 114 | formyl | — | 甲酰基 | ✅ 五大保留酰基PIN之一 | ⚠️ 作醛母体，非前缀 |
| 115 | acetyl | — | 乙酰基 | ✅ 五大保留酰基PIN之一 | ⚠️ acetyl chloride/bromide 保留 |
| 116 | benzoyl | — | 苯甲酰基 | ✅ 五大保留酰基PIN之一 | ✅ benzoyl Cl/Br |
| 117 | oxalyl (二价) | — | 草酰基/乙二酰基 | ✅ 五大保留酰基PIN之一 | ❌ |
| 118 | oxamoyl | — | 草氨酰基 | ✅ 五大保留酰基PIN之一 | ❌ |
| 119 | propanoyl | propionyl (非PIN) | 丙酰基 | ✅ propanoyl 系统 | ✅ propanoyl 系统 |
| 120 | butanoyl | butyryl (非PIN) | 丁酰基 | ✅ butanoyl 系统 | ✅ butanoyl 系统 |
| 111 | 2-methylpropanoyl | isobutyryl (非PIN) | 异丁酰基 | ✅ isobutyryl bromide |
| 112 | pentanoyl | valeryl (非PIN) | 戊酰基 | ✅ pentanoyl 系统 |
| 113 | benzoyl | — | 苯甲酰基 | ✅ benzoyl Cl/Br |
| 114 | oxalyl (二价) | — | 草酰基/乙二酰基 | ❌ |
| 115 | malonyl (二价) | — | 丙二酰基 | ❌ |
| 116 | succinyl (二价) | — | 丁二酰基 | ❌ |
| 117 | glutaryl (二价) | — | 戊二酰基 | ❌ |
| 118 | adipoyl (二价) | — | 己二酰基 | ❌ |
| 119 | phthaloyl (二价) | — | 邻苯二甲酰基 | ❌ |
| 120 | isophthaloyl (二价) | — | 间苯二甲酰基 | ❌ |
| 121 | terephthaloyl (二价) | — | 对苯二甲酰基 | ❌ |
| 122 | carboxy | — | 羧基 | ✅ 酸母体 |
| 123 | methoxycarbonyl | — | 甲氧羰基 | ❌ |
| 124 | ethoxycarbonyl | — | 乙氧羰基 | ❌ |
| 125 | carbamoyl | — | 氨基甲酰基 | ❌ |
| 126 | sulfamoyl | — | 氨基磺酰基 | ❌ |

---

## 7. 项目当前明确缺失的高优先级取代基

### Tier 1 — IUPAC 承认（PIN 或 General），化工实践中不可或缺

| 优先级 | 取代基 | PIN/General | 说明 |
|--------|--------|------------|------|
| 🔴 | **vinyl** | General only | 太常用，ethenyl 在实践中极少见 |
| 🔴 | **allyl** | General only | 同上 |
| 🔴 | **isopropenyl** | General only | 同上 |
| 🟡 | **methylene** | PIN | 二价，桥环必需 |
| 🟡 | **1,2/1,3/1,4-phenylene** | PIN | 二价，聚合物命名必需 |
| 🟡 | **oxalyl / malonyl / succinyl** | PIN（二价酰基） | 二酸衍生对称取代基 |

### Tier 2 — 普通命名常用但非 PIN

| 优先级 | 取代基 | 说明 |
|--------|--------|------|
| 🟡 | **formyl** (作为前缀) | 当前仅作醛母体 |
| 🟡 | **acetyl** (作为前缀) | 仅 acyl chloride/bromide 保留 |
| 🟡 | **benzoyl** (作为通用前缀) | 泛化至 benzoyl-X 类型 |
| 🟢 | **nitroso** | 亚硝基基本前缀 |
| 🟢 | **azido** | 叠氮基基本前缀 |
| 🟢 | **carbamoyl / sulfamoyl** | 酰胺/磺酰胺语境 |
| 🟢 | **cyano** (自由前缀) | 当前仅在芳环叶上作为叶 |

### Tier 3 — 保护基/工业常用名，非 PIN

| 优先级 | 取代基 | 说明 |
|--------|--------|------|
| 🟢 | tosyl, mesyl, triflyl | 保护基语境，非 PIN 但极常用 |
| 🟢 | isobutyryl, pivaloyl | 工业名，非 PIN |
| 🟢 | ethylene, trimethylene (二价) | 普通命名保留 |

---

## 8. IUPAC 2013 vs 1993 关键变化

| 取代基 | 1993 Guide | 2013 Blue Book | 项目当前行为 |
|--------|-----------|----------------|------------|
| isopropyl | Retained | **General only** (PIN=propan-2-yl) | ⚠️ 使用保留名 |
| isobutyl | Retained | **Not recommended** (PIN=2-methylpropyl) | ⚠️ 使用保留名 |
| sec-butyl | Retained | **Not recommended** (PIN=butan-2-yl) | ⚠️ 使用保留名 |
| neopentyl | Retained | **Not recommended** (PIN=2,2-dimethylpropyl) | ⚠️ 使用保留名 |
| isopentyl | Retained | **Not recommended** (PIN=3-methylbutyl) | ⚠️ 使用保留名 |
| tert-pentyl | Retained | **Not recommended** (PIN=2-methylbutan-2-yl) | ✅ 使用系统名 |
| vinyl | Retained | **General only** (PIN=ethenyl) | ❌ 缺失 |
| allyl | Retained | **General only** (PIN=prop-2-en-1-yl) | ❌ 缺失 |
| isopropenyl | Retained | **General only** (PIN=prop-1-en-2-yl) | ❌ 缺失 |
| tert-butyl | Retained | **PIN** (唯一支链烷基 PIN！) | ✅ 使用保留名 |
| phenyl | Retained | **PIN** | ✅ |
| benzyl | Retained | **PIN** | ✅ |

---

## 9. 决策要点

项目的核心问题是：**IUPAC 2013 大幅缩减了保留名，但化工界实际使用并未跟进。**

需要考虑两个维度：

### 维度 A：PIN 严格合规
- 支链烷基几乎全部使用系统名：propan-2-yl, 2-methylpropyl, butan-2-yl...
- vinyl/allyl → ethenyl/prop-2-en-1-yl
- 仅保留 tert-butyl, phenyl, benzyl, methylene, phenylene

### 维度 B：化工实践（大多数期刊/软件的做法）
- 保留 isopropyl, isobutyl, sec-butyl, neopentyl, isopentyl
- 保留 vinyl, allyl, isopropenyl
- propargyl 可保留（实践中远多于 prop-2-yn-1-yl）
- 用系统名处理 tert-pentyl → 2-methylbutan-2-yl（项目目前已这样做）

### 建议：采用"双轨制"
类似于项目已有的 `rooted_tree` 系统名 + `retained` 保留名的 backend 架构，可以明确定义：
- **PIN mode**：所有取代基严格使用 IUPAC 2013 PIN
- **General mode**（默认）：使用化工界公认的保留名，即便不是 PIN

---

## 10. 架构建议：集中化 Retained 注册表

当前"地位特殊取代基"分散在 ~8 个文件中。推荐一个统一注册表：

```python
# src/namepredict/layer3/retained_substituents.py
from dataclasses import dataclass
from enum import Enum

class IupacLevel(Enum):
    PIN = "pin"                    # IUPAC Preferred
    GENERAL = "general"            # Allowed in general nomenclature
    NOT_RECOMMENDED = "not_rec"    # 2013 deprecated
    SYSTEMATIC = "systematic"      # Not a retained name

@dataclass(frozen=True)
class RetainedSubstituent:
    en: str           # retained/common name
    zh: str           # Chinese name
    systematic_en: str  # systematic PIN equivalent
    systematic_zh: str
    level: IupacLevel
    rule_ref: str     # e.g. "P-29.6.1", "P-57.1.2"

RETAINED_MONOVALENT: dict[str, RetainedSubstituent] = {
    # PIN
    "tert-butyl": RetainedSubstituent("tert-butyl", "叔丁基",
        "1,1-dimethylethyl", "1,1-二甲基乙基", IupacLevel.PIN, "P-57.1.2"),
    "phenyl": RetainedSubstituent("phenyl", "苯基",
        "", "", IupacLevel.PIN, "P-29.6"),
    "benzyl": RetainedSubstituent("benzyl", "苄基",
        "phenylmethyl", "苯甲基", IupacLevel.PIN, "P-57.1.2"),
    # General only
    "vinyl": RetainedSubstituent("vinyl", "乙烯基",
        "ethenyl", "乙烯基", IupacLevel.GENERAL, "P-29.6.2.2"),
    "allyl": RetainedSubstituent("allyl", "烯丙基",
        "prop-2-en-1-yl", "丙-2-烯-1-基", IupacLevel.GENERAL, "P-29.6.2.2"),
    "isopropenyl": RetainedSubstituent("isopropenyl", "异丙烯基",
        "prop-1-en-2-yl", "丙-1-烯-2-基", IupacLevel.GENERAL, "P-29.6.2.2"),
    "isopropyl": RetainedSubstituent("isopropyl", "异丙基",
        "propan-2-yl", "丙-2-基", IupacLevel.GENERAL, "P-29.6.2.2"),
    # Not recommended (but common)
    "isobutyl": RetainedSubstituent("isobutyl", "异丁基",
        "2-methylpropyl", "2-甲基丙基", IupacLevel.NOT_RECOMMENDED, "P-57.1.4"),
    "sec-butyl": RetainedSubstituent("sec-butyl", "仲丁基",
        "butan-2-yl", "丁-2-基", IupacLevel.NOT_RECOMMENDED, "P-57.1.4"),
    "neopentyl": RetainedSubstituent("neopentyl", "新戊基",
        "2,2-dimethylpropyl", "2,2-二甲基丙基", IupacLevel.NOT_RECOMMENDED, "P-57.1.4"),
    "isopentyl": RetainedSubstituent("isopentyl", "异戊基",
        "3-methylbutyl", "3-甲基丁基", IupacLevel.NOT_RECOMMENDED, "P-57.1.4"),
    # ... 100+ more
}
```

这样每次新增一个"地位特殊的取代基"只需注册一个条目，无需改动多个 backend 和提取器。
