# IUPAC 2013 Blue Book: "地位特殊的母体" 全面名单

> 依据：IUPAC Blue Book (2013) P-21, P-22, P-23, P-24, P-25, P-31, P-44, P-52, P-62–P-66；
> 项目现状：`src/namepredict/layer2/` scaffold, kind_registry, ring/fg/unsat_producers。
>
> 母体 (parent hydride / parent structure) 是系统命名的骨架基础，取代基和功能基以此为根进行命名。
> IUPAC 2013 保留了相当数量的母体惯用名。

---

## 1. 无环饱和烃母体 (P-21.2.1: Alkanes)

### 1a. 保留母体名（前 4 个为保留 PIN）

| n | 保留名 (PIN) | 中文 | 系统名 | 项目状态 |
|---|------------|------|--------|---------|
| 1 | **methane** | 甲烷 | carbane (从未采用) | ✅ |
| 2 | **ethane** | 乙烷 | — | ✅ |
| 3 | **propane** | 丙烷 | — | ✅ |
| 4 | **butane** | 丁烷 | — | ✅ |
| 5+ | pentane, hexane, heptane... | 戊烷… | 系统构造名 (PIN) | ✅ C5–C35 |

> C1–C4 为 IUPAC 保留名；C5+ 为系统构造名但也是 PIN。IUPAC 2013 推荐 C20+ 使用 `icos-` 拼写（非 `eicos-`）。

### 1b. 无环不饱和烃母体 (P-31.1)

| 类型 | 命名规则 | 保留名 | 项目状态 |
|------|---------|--------|---------|
| alkene | -ane→-ene, 带位次 | — | ✅ alkene/polyene |
| alkyne | -ane→-yne, 带位次 | acetylene (C2), propyne (C3) 可接受 | ✅ alkyne |

---

## 2. 单环烃母体 (P-22.1: Monocyclic Hydrocarbons)

### 2a. 饱和单环 (Cycloalkanes) — 系统 PIN

| n | 名称 (PIN) | 中文 | 项目状态 |
|---|-----------|------|---------|
| 3 | **cyclopropane** | 环丙烷 | ✅ |
| 4 | **cyclobutane** | 环丁烷 | ✅ |
| 5 | **cyclopentane** | 环戊烷 | ✅ |
| 6 | **cyclohexane** | 环己烷 | ✅ |
| 7 | **cycloheptane** | 环庚烷 | ✅ |
| 8 | **cyclooctane** | 环辛烷 | ✅ |
| 9 | **cyclononane** | 环壬烷 | ✅ (Spec) |
| 10 | **cyclodecane** | 环癸烷 | ✅ (Spec) |
| 11+ | cycloundecane... | 环十一烷… | ❌ |

### 2b. 苯系芳烃 — 保留母体名

| # | 保留名 | IUPAC 2013 级别 | 中文 | 系统名 | 项目状态 |
|---|--------|----------------|------|--------|---------|
| 1 | **benzene** | ✅ PIN (P-22.1.3) | 苯 | cyclohexatriene (非PIN) | ✅ |
| 2 | **toluene** | ✅ PIN (unsubstituted) | 甲苯 | 不可取代(PIN)；普通命名仅限卤/硝基/亚硝基/–OR/–SO-R/–SO₂-R (P-15.1.8.2.2) | ✅ |
| 3 | **xylene** (o/m/p) | ✅ PIN (unsubstituted) | 二甲苯 | **禁止任何取代** (PIN或普通命名均不可) | ✅ (en xylene) |
| 4 | **mesitylene** | ❌ General only | 均三甲苯 | 1,3,5-trimethylbenzene (PIN) | ❌ |
| 5 | **styrene** | ❌ General only | 苯乙烯 | ethenylbenzene (PIN) | ❌ |
| 6 | **cumene** | ❌ **2013年未保留** | 异丙苯 | isopropylbenzene (PIN) | ❌ |
| 7 | **cymene** (o/m/p) | ❌ **2013年未保留** | 伞花烃 | methylisopropylbenzene (PIN) | ❌ |

### 2c. 环烯烃与环多烯 (P-22.1.1 / P-31.1)

| 名称 (PIN) | 中文 | 保留? | 项目状态 |
|-----------|------|-------|---------|
| cyclohexene | 环己烯 | 系统 PIN | ✅ cycloalkene |
| cyclohexa-1,3-diene | 环己-1,3-二烯 | 系统 PIN | ✅ cyclopolyene |
| cyclohexa-1,4-diene | 环己-1,4-二烯 | 系统 PIN | ✅ cyclopolyene |
| cyclopentadiene | 环戊二烯 | 系统 PIN | ⚠️ |

> 注：无环多烯母体名（buta-1,3-diene, hexa-1,3,5-triene 等）均为系统构造 PIN。

---

## 3. 稠环芳烃母体 (P-25: Fused Polycyclic Systems)

### 3a. Table 2.7 — 保留 PIN 碳环稠环（按 seniority 降序）

> P-25.1.1 Table 2.7 共 19 个保留母体，seniority 越高越优先选为稠合母体组分。

| Seniority | 保留名 (PIN) | 中文 | 分子式 | 环数 | 编号 | 项目状态 |
|-----------|------------|------|--------|------|------|---------|
| 1 | **ovalene** | 卵苯 | C32H14 | 10 | 系统 | ❌ |
| 2 | **pyranthrene** | 皮蒽 | C30H16 | 8 | 系统 | ❌ |
| 3 | **coronene** | 晕苯 | C24H12 | 7 | 系统 | ❌ |
| 4 | **rubicene** | 茹比省 | C26H14 | 7 | 系统 | ❌ |
| 5 | **perylene** | 苝 | C20H12 | 5 | 系统 | ❌ |
| 7 | **pleiadene** | 七曜烯 | C18H12 | 4 | 系统 | ❌ |
| 8 | **chrysene** | 屈 | C18H12 | 4 | 系统 | ❌ |
| 9 | **pyrene** | 芘 | C16H10 | 4 | 系统 | ❌ |
| 10 | **fluoranthene** | 荧蒽 | C16H10 | 4 | 系统 | ❌ |
| 11 | **anthracene** | 蒽 | C14H10 | 3 | **传统编号例外** | ✅ |
| 12 | **phenanthrene** | 菲 | C14H10 | 3 | **传统编号例外** | ❌ |
| 13 | **phenalene** (1H) | 菲那烯 | C13H10 | 3 | 系统 | ❌ |
| 14 | **fluorene** (9H) | 芴 | C13H10 | 3 | 系统 | ❌ |
| 15 | **s-indacene** | 对称引达省 | C12H8 | 3 | 系统 | ❌ |
| 16 | **as-indacene** | 不对称引达省 | C12H8 | 3 | 系统 | ❌ |
| 17 | **azulene** | 薁 | C10H8 | 2 | 系统 | ❌ |
| 18 | **naphthalene** | 萘 | C10H8 | 2 | 系统 (4a,8a) | ✅ |
| 19 | **indene** (1H) | 茚 | C9H8 | 2 | 系统 | ❌ |

### 3b. 其他碳环保留/系统名

| 名称 | 级别 | 中文 | 项目状态 |
|------|------|------|---------|
| **tetracene** (naphthacene) | Polyacene 系统 PIN | 并四苯 | ❌ |
| **pentacene** | Polyacene 系统 PIN | 并五苯 | ❌ |
| **biphenylene** | PIN (P-25.1.2.4) | 亚联苯 | ❌ |
| **acenaphthylene** | PIN (半系统) | 苊 | ❌ |
| **indane** | ⚠️ General only (PIN = 2,3-dihydro-1H-indene) | 茚满 | ❌ |
| **[60]fullerene** | PIN (P-27) | 富勒烯-C60 | ❌ |

### 3c. 7 个传统编号例外 (P-25.3.3.3)

| 体系 | 类型 |
|------|------|
| **anthracene** | 保留 PIN |
| **phenanthrene** | 保留 PIN |
| **acridine** | 保留 PIN (N=10) |
| **carbazole** (9H) | 保留 PIN (N=9) |
| **xanthene** (9H) + S/Se/Te 类似物 | 保留 PIN |
| **purine** (7H) | 保留 PIN |
| **cyclopenta[a]phenanthrene** | 系统融合名（唯一一个非保留名的例外） |

### 3d. 9 个保留融合前缀 PIN (P-25.3.2.2.3)

> 仅以下收缩前缀可作为 PIN 用于融合命名：**anthra, naphtho, benzo, phenanthro, furo, imidazo, pyrido, pyrimido, thieno**
> 
> ⚠️ `acenaphtho, perylo, isoquino, quino` 仅限普通命名。

| 保留名 | 中文 | PIN? | 项目状态 |
|--------|------|------|---------|
| 1,4-dihydronaphthalene | 1,4-二氢萘 | 系统 | ❌ |
| 9,10-dihydroanthracene | 9,10-二氢蒽 | 系统 | ❌ |

---

## 4. 杂环母体 (P-22.2: Heterocyclic Parent Hydrides)

### 4a. 保留名的单杂环芳烃 (Table 2.2, P-22.2.1)

**五元单杂:**

| # | 保留名 | 中文 | 项目状态 |
|---|--------|------|---------|
| 1 | **furan** | 呋喃 | ✅ |
| 2 | **thiophene** | 噻吩 | ✅ |
| 3 | **selenophene** | 硒吩 | ❌ |
| 4 | **tellurophene** | 碲吩 | ❌ |
| 5 | **1H-pyrrole** | 吡咯 | ✅ |

**五元 1,3-二杂:**

| # | 保留名 | 中文 | 项目状态 |
|---|--------|------|---------|
| 6 | **1H-imidazole** | 咪唑 | ✅ |
| 7 | **1H-pyrazole** | 吡唑 | ✅ |
| 8 | **1,3-oxazole** | 噁唑 | ✅ |
| 9 | **1,3-thiazole** | 噻唑 | ✅ |
| 10 | **1,3-selenazole** | 硒唑 | ❌ |
| 11 | **1,3-tellurazole** | 碲唑 | ❌ |

**五元 1,2-二杂:**

| # | 保留名 | 中文 | 项目状态 |
|---|--------|------|---------|
| 12 | **isoxazole** (1,2-oxazole) | 异噁唑 | ❌ |
| 13 | **isothiazole** (1,2-thiazole) | 异噻唑 | ❌ |

**六元单杂:**

| # | 保留名 | 中文 | 项目状态 |
|---|--------|------|---------|
| 14 | **pyridine** | 吡啶 | ✅ |
| 15 | **2H-pyran** | 吡喃 | ❌ |
| 16 | **2H-thiopyran** | 噻喃 | ❌ |

**六元二氮:**

| # | 保留名 | 中文 | 项目状态 |
|---|--------|------|---------|
| 17 | **pyrimidine** | 嘧啶 | ✅ |
| 18 | **pyrazine** | 吡嗪 | ✅ |
| 19 | **pyridazine** | 哒嗪 | ✅ |

> 注: triazole, tetrazole, oxadiazole, thiadiazole, triazine 等更多 N 杂环保留名也存在于 Table 2.2。

### 4b. 保留名的饱和单杂环 (P-22.2.2)

> **关键区分**：仅 **6 个**含氮饱和杂环是 IUPAC 保留 PIN（Table 2.3），其余饱和杂环使用 Hantzsch-Widman 系统名作为 PIN。

#### 保留 PIN（Table 2.3 六母体 + 硫代类似物）

| PIN 保留名 | HW 系统名 | 中文 | 环大小 | 杂原子 | 项目状态 |
|-----------|----------|------|--------|--------|---------|
| **pyrrolidine** | azolidine | 吡咯烷 | 5 | 1N | ✅ |
| **piperidine** | azinane | 哌啶 | 6 | 1N | ✅ |
| **morpholine** | 1,4-oxazinane | 吗啉 | 6 | 1O+1N | ✅ |
| **piperazine** | 1,4-diazinane | 哌嗪 | 6 | 2N | ✅ |
| **imidazolidine** | 1,3-diazolidine | 咪唑烷 | 5 | 2N | ❌ |
| **pyrazolidine** | 1,2-diazolidine | 吡唑烷 | 5 | 2N | ❌ |
| **thiomorpholine** | — | 硫吗啉 | 6 | 1S+1N | ❌ |

#### Hantzsch-Widman 系统 PIN（非"保留"，但是 PIN）

| PIN 名 | 中文 | 杂原子 | 环大小 | 项目状态 |
|--------|------|--------|--------|---------|
| **aziridine** | 氮杂环丙烷 | N | 3 | ✅ |
| **azetidine** | 氮杂环丁烷 | N | 4 | ❌ |
| **oxirane** | 环氧乙烷 | O | 3 | ✅ |
| **oxetane** | 氧杂环丁烷 | O | 4 | ❌ |
| **oxolane** | 氧杂环戊烷 | O | 5 | ✅ (**THF 非 PIN**) |
| **oxane** | 氧杂环己烷 | O | 6 | ✅ (**THP 非 PIN**) |
| **oxepane** | 氧杂环庚烷 | O | 7 | ❌ |
| **thiolane** | 硫杂环戊烷 | S | 5 | ✅ |
| **thiane** | 硫杂环己烷 | S | 6 | ❌ |
| **1,4-dioxane** | 1,4-二氧六环 | 2O | 6 | ✅ |
| **1,3-dioxolane** | 1,3-二氧戊环 | 2O | 5 | ✅ |
| **1,3-dioxane** | 1,3-二氧六环 | 2O | 6 | ❌ |
| **thiazolidine** | 噻唑烷 | S,N | 5 | ❌ |
| **oxazolidine** | 噁唑烷 | O,N | 5 | ❌ |

> 注：aziridine/azetidine 虽非 Table 2.3 保留名，但因含氮使用传统词干 (-iridine/-etidine)，在 PIN 中优先于替换命名法。

---

## 5. 稠杂环母体 (P-25: Fused Heterocycles)

### 5a. 5+6 稠合体系 — 非保留 PIN（融合系统名才是 PIN）

> **关键澄清**：benzofuran/benzothiophene/benzoxazole/benzothiazole/benzimidazole **不是 IUPAC 保留 PIN**！
> 它们的 PIN 是融合系统名 (P-25.2.2)：`1-benzofuran`, `1-benzothiophene`, `1,3-benzoxazole`, `1,3-benzothiazole`, `1H-benzimidazole`。
> 项目 ScaffoldSpec 将其标为 `retained=True`，但这应理解为"可直接使用该名作为母体名"，而非 IUPAC 严格意义上的保留 PIN。

| # | 保留名 | 中文 | 编号 | 项目状态 |
|---|--------|------|------|---------|
| 1 | **1H-indole** | 吲哚 | N=1…7a | ✅ |
| 2 | **isoindole** | 异吲哚 | | ❌ |
| 3 | **1H-indazole** | 吲唑 | N=1,N=2…7a | ✅ |
| 4 | **1H-benzimidazole** | 苯并咪唑 | N=1,N=3…7a | ✅ |
| 5 | **benzofuran** (1-benzofuran) | 苯并呋喃 | O=1…7a | ✅ |
| 6 | **isobenzofuran** | 异苯并呋喃 | | ❌ |
| 7 | **1-benzothiophene** | 苯并[b]噻吩 | S=1…7a | ✅ |
| 8 | **1,3-benzothiazole** | 苯并噻唑 | S=1,N=3…7a | ✅ |
| 9 | **1,3-benzoxazole** | 苯并噁唑 | O=1,N=3…7a | ✅ |
| 10 | **1,2-benzisothiazole** | 1,2-苯并异噻唑 | | ❌ |
| 11 | **1,2-benzisoxazole** | 1,2-苯并异噁唑 | | ❌ |

### 5b. 6+6 稠合体系

| # | 保留名 | 中文 | 编号 | 项目状态 |
|---|--------|------|------|---------|
| 12 | **quinoline** | 喹啉 | N=1…8a | ✅ |
| 13 | **isoquinoline** | 异喹啉 | N=2…8a | ✅ |
| 14 | **quinazoline** | 喹唑啉 | 1,3-N₂ | ✅ |
| 15 | **quinoxaline** | 喹喔啉 | 1,4-N₂ | ✅ |
| 16 | **cinnoline** | 噌啉 | 1,2-N₂ | ❌ |
| 17 | **phthalazine** | 酞嗪 | 2,3-N₂ | ❌ |
| 18 | **chromene** (2H/4H) | 色烯 | O=1…8a | ❌ |
| 19 | **coumarin** | 香豆素 | ⚠️ non-PIN (PIN = 2H-1-benzopyran-2-one) | ✅ (chromen-2-one) |
| 20 | **isocoumarin** | 异香豆素 | | ❌ |

### 5c. 三环及以上稠杂环 — 保留 PIN

| # | 保留名 | 中文 | 项目状态 |
|---|--------|------|---------|
| 21 | **carbazole** | 咔唑 (9H) | ❌ |
| 22 | **acridine** | 吖啶 | ❌ |
| 23 | **phenazine** | 吩嗪 | ❌ |
| 24 | **phenoxazine** | 吩噁嗪 | ❌ |
| 25 | **phenothiazine** | 吩噻嗪 | ❌ |
| 26 | **xanthene** | 呫吨 (9H) | ❌ |
| 27 | **purine** (7H/9H) | 嘌呤 | ❌ |
| 28 | **pteridine** | 蝶啶 | ❌ |
| 29 | **thianthrene** | 噻蒽 | ❌ |
| 30 | **phenanthridine** | 菲啶 | ❌ |
| 31 | **phenanthroline** (o/m/p) | 菲咯啉 (1,7/1,10/4,7) | ❌ |

---

## 6. 桥环与螺环母体 (P-23 / P-24)

### 6a. 保留名的桥环烃

| # | 保留名 | 中文 | 系统名 | 项目状态 |
|---|--------|------|--------|---------|
| 1 | **adamantane** | 金刚烷 | tricyclo[3.3.1.1³·⁷]decane | ✅ PIN | ✅ (作为取代基) |
| 2 | **cubane** | 立方烷 | pentacyclo[4.2.0.0²·⁵.0³·⁸.0⁴·⁷]octane | ✅ PIN | ❌ |
| 3 | **norbornane** | 降冰片烷 | bicyclo[2.2.1]heptane | ⚠️ General only | ❌ |
| 4 | **quinuclidine** | 奎宁环 | 1-azabicyclo[2.2.2]octane | ⚠️ General only | ❌ |
| 5 | **basketane** | 篮烷 | pentacyclo[4.4.0.0²·⁵.0³·⁸.0⁴·⁷]decane | ⚠️ General only | ❌ |
| 6 | **prismane** | 棱柱烷 | tetracyclo[2.2.0.0²·⁶.0³·⁵]hexane | ❌ **2013废除** | ❌ |

### 6b. 螺环体系 — 无保留名

> 螺环体系全部使用系统命名 (`spiro[x.y]alkane`)，无 IUPAC 保留名。

---

## 7. Hantzsch-Widman 系统单杂环命名 (P-22.2.3)

当杂单环**没有**保留名时，使用 Hantzsch-Widman 系统名。

### 7a. 词干表 (Stem Table)

| 环大小 | 含N (不饱和) | 含N (饱和) | 无N (不饱和) | 无N (饱和) |
|--------|------------|-----------|------------|-----------|
| 3 | **-irine** | **-iridine** | **-irene** | **-irane** |
| 4 | **-ete** | **-etidine** | **-ete** | **-etane** |
| 5 | **-ole** | **-olidine** | **-ole** | **-olane** |
| 6A (O,S,Se,Te,Bi) | **-ine** | — | **-ine** | **-ane** |
| 6B (N,Si,Ge,Sn,Pb) | **-ine** | **-inane** | — | — |
| 6C (P,As,Sb; B,Al,Ga,In,Tl; 卤素) | **-inine** | **-inane** | — | — |
| 7 | **-epine** | — | **-epine** | **-epane** |
| 8 | **-ocine** | — | **-ocine** | **-ocane** |
| 9 | **-onine** | — | **-onine** | **-onane** |
| 10 | **-ecine** | — | **-ecine** | **-ecane** |

### 7b. 杂原子前缀表（优先级降序）

| 前缀 | 杂原子 | 价态 |
|------|--------|------|
| **oxa** | O | II |
| **thia** | S | II |
| **selena** | Se | II |
| **tellura** | Te | II |
| **aza** | N | III |
| **phospha** | P | III |
| **arsa** | As | III |
| **stiba** | Sb | III |
| **bisma** | Bi | III |
| **sila** | Si | IV |
| **germa** | Ge | IV |
| **stanna** | Sn | IV |
| **plumba** | Pb | IV |
| **bora** | B | III |

> elision 规则：两个 `a` 相邻时，去掉一个（如 oxa + azete → oxazete）。

### 7c. Hantzsch-Widman 系统名示例

| 系统名 | 中文 | 可替代的保留名 |
|--------|------|-------------|
| 1,4-oxazine | 噁嗪(1,4) | — (无保留名) |
| 1,3-dioxole | 二氧杂环戊烯(1,3) | — |
| oxaziridine | 氧杂氮丙啶 | — |
| thiaziridine | 硫杂氮丙啶 | — |
| 1,4,7-trioxonane | 1,4,7-三氧杂环壬烷 | — (冠醚) |

---

## 8. 含功能基的保留母体名 (Functional Parents)

这类母体**不是纯碳氢骨架**，而是含特征基团（羧基、羰基等）的保留完整结构名。

### 8a. 羧酸保留母体 (P-65.1)

**仅 5 个羧酸是 PIN (P-65.1.1.1)：**

| 保留名 | 中文 | 级别 | 项目状态 |
|--------|------|------|---------|
| **formic acid** | 甲酸 | ✅ PIN | ✅ |
| **acetic acid** | 乙酸 | ✅ PIN | ✅ |
| **benzoic acid** | 苯甲酸 | ✅ PIN | ✅ |
| **oxalic acid** | 草酸 | ✅ PIN | ✅ |
| **oxamic acid** | 草氨酸 | ✅ PIN | ❌ |

**⚠️ General only (非 PIN，不可在 PIN 中用作母体)：**

| 保留名 | 中文 | PIN 替代 | 项目状态 |
|--------|------|---------|---------|
| propionic acid | 丙酸 | propanoic acid | ❌ |
| butyric acid | 丁酸 | butanoic acid | ❌ |
| acrylic acid | 丙烯酸 | prop-2-enoic acid | ❌ |
| malonic acid | 丙二酸 | propanedioic acid | ❌ |
| succinic acid | 丁二酸 | butanedioic acid | ❌ |
| glutaric acid | 戊二酸 | pentanedioic acid | ❌ |
| adipic acid | 己二酸 | hexanedioic acid | ❌ |
| phthalic acid | 邻苯二甲酸 | benzene-1,2-dicarboxylic acid | ❌ |
| isophthalic acid | 间苯二甲酸 | benzene-1,3-dicarboxylic acid | ❌ |
| terephthalic acid | 对苯二甲酸 | benzene-1,4-dicarboxylic acid | ❌ |
| cinnamic acid | 肉桂酸 | (2E)-3-phenylprop-2-enoic acid | ❌ |

> ⚠️ salicylic acid, gallic acid 等传统名在 2013 年**不在保留功能母体列表中**。

### 8b. 醛保留母体 (P-66.6)

| 保留名 | 中文 | 级别 | 项目状态 |
|--------|------|------|---------|
| **benzaldehyde** | 苯甲醛 | ✅ PIN（唯一醛 PIN！） | ✅ |
| formaldehyde | 甲醛 | ⚠️ General only (PIN=methanal) | ✅ |
| acetaldehyde | 乙醛 | ⚠️ General only (PIN=ethanal) | ✅ |

> formaldehyde/acetaldehyde 均非 PIN。**仅有 benzaldehyde 是醛类保留 PIN**。

### 8c. 酮保留母体 (P-64)

| 保留名 | 中文 | 级别 | 注意 | 项目状态 |
|--------|------|------|------|---------|
| **chalcone** | 查尔酮 | ✅ PIN（唯一酮 PIN！） | 仅反式异构体 | ❌ |
| acetone | 丙酮 | ⚠️ General only | PIN = propan-2-one | ❌ |
| acetophenone | 苯乙酮 | ⚠️ General only | **不可被取代** | ✅ (作为母体) |
| benzophenone | 二苯甲酮 | ⚠️ General only | **不可被取代** | ❌ |

> **benzil, biacetyl, propiophenone** 在 2013 年**明确废除** (P-64.2.1.3)。

### 8d. 醇/酚保留母体 (P-63)

| 保留名 | 中文 | 级别 | 项目状态 |
|--------|------|------|---------|
| **phenol** | 苯酚 | ✅ PIN（唯一酚 PIN！） | ✅ |
| cresol (o/m/p) | 甲酚 | ⚠️ General only (OH 必须位次=1) | ❌ |
| 1-naphthol / 2-naphthol | 萘酚 | ⚠️ General only（不可被取代） | ❌ |
| pyrocatechol / resorcinol / hydroquinone | 苯二酚 | ⚠️ General only（不可被取代） | ❌ |
| ethylene glycol / glycerol | 乙二醇/丙三醇 | ⚠️ General only | ❌ |

> **phenanthrol** — 2013 年明确不再保留。

### 8e. 胺保留母体 (P-62)

| 保留名 | 中文 | 级别 | 项目状态 |
|--------|------|------|---------|
| **aniline** | 苯胺 | ✅ PIN（唯一胺 PIN！） | ✅ |
| benzidine | 联苯胺 | ⚠️ General only (仅4,4'-异构体) | ❌ |

> **toluidine, anisidine, phenetidine, xylidine — 全部被废除** (P-62.2.1.1.2)，必须用系统名。

### 8f. 酰胺保留母体 (P-66.1)

| # | 保留名 | 中文 | 级别 | 项目状态 |
|---|--------|------|------|---------|
| 1 | **formamide** | 甲酰胺 | ✅ PIN | ✅ |
| 2 | **acetamide** | 乙酰胺 | ✅ PIN | ✅ |
| 3 | **benzamide** | 苯甲酰胺 | ✅ PIN | ✅ |
| 4 | **urea** | 脲/尿素 | ✅ PIN | ✅ |
| 5 | **oxamide** | 草酰胺 | ⚠️ General only | ❌ |

### 8g. 腈保留母体 (P-66.5)

| # | 保留名 | 中文 | 级别 | 项目状态 |
|---|--------|------|------|---------|
| 1 | **formonitrile** | 甲腈 | ✅ PIN | ✅ |
| 2 | **acetonitrile** | 乙腈 | ✅ PIN | ✅ |
| 3 | **benzonitrile** | 苯甲腈 | ✅ PIN | ✅ |

### 8h. 其他功能母体

| # | 保留名 | 中文 | 类型 | 项目状态 |
|---|--------|------|------|---------|
| 1 | **guanidine** | 胍 | PIN | ✅ |
| 2 | **hydrazine** | 肼 | PIN | ✅ |
| 3 | **sulfanilic acid** | 对氨基苯磺酸 | ⚠️ General | ❌ |

---

## 9. 母体选择层级规则 (P-44 / P-52)

IUPAC 2013 的母体选择有严格优先级（用于 PIN）：

### 9a. 特征基团优先 (P-44.1.1)
首先选含最优先特征基团（后缀）的结构为母体。
优先级: `acid > anhydride > ester > acyl halide > amide > nitrile > aldehyde > ketone > alcohol > thiol > amine > ether > ...`

### 9b. 环优先于链 (P-44.1.2.2 / P-52.2.8)
当特征基团相同时，**环 > 链**是 PIN 的铁律。
- PIN: `cyclopropylpentane` (环优先) → 实际 PIN 是 pentylcyclopropane
- ⚠️ General: `1-cyclopropylpentane` (链优先，大者优先)

### 9c. 杂原子数 (P-44.2)
环体系之间：杂原子多的优先。

### 9d. 环数 (P-44.2.2.2)
杂原子数相同时，环数多的优先。

### 9e. 环大小 (P-44.2.2.3)
环数相同时，大环优先。

### 9f. 不饱和度 (P-44.2.2.4)
不饱和度高的优先。

### 9g. 杂原子优先级 (P-25.6)
O > S > N > P > ...

---

## 10. 项目当前母体覆盖总结

### 已实现的保留母体（约 50 个）

| 类别 | 已实现 | 主要缺失 |
|------|--------|---------|
| 无环烷烃 | C1–C35 (✅) | — |
| 无环不饱和 | alkene/alkyne/polyene (✅) | cumulene |
| 单环碳环 | cycloalkane C3–C10, cycloalkene, cyclopolyene (✅) | C11+ |
| 苯系 | benzene, toluene, xylene (✅) | mesitylene, styrene |
| 稠环碳环 | naphthalene, anthracene (✅) | **phenanthrene**, pyrene, azulene, fluorene, indene |
| 单杂芳环 | pyridine, furan, thiophene, pyrrole, imidazole, pyrazole, oxazole, thiazole, pyrimidine, pyrazine, pyridazine (✅) | **triazole, tetrazole**, oxadiazole, thiadiazole, triazine |
| 饱和杂环 | aziridine, oxirane, oxolane, oxane, pyrrolidine, piperidine, morpholine, piperazine, dioxolane, dioxane, thiolane (✅) | **azetidine**, azepane, oxetane, thiane, thiazolidine, imidazolidine |
| 5+6稠杂环 | indole, indazole, benzimidazole, benzofuran, benzothiophene, benzothiazole, benzoxazole (✅) | **isoindole**, isobenzofuran |
| 6+6稠杂环 | quinoline, isoquinoline, quinazoline, quinoxaline (✅) | **cinnoline, phthalazine**, chromene |
| 三环+杂环 | — (❌) | **carbazole, acridine, purine, pteridine**, phenazine |
| 功能母体 | formic/acetic/oxalic acid, formaldehyde, acetaldehyde, benzaldehyde, acetophenone, phenol, aniline, benzamide, formamide, acetamide, benzonitrile, formonitrile, acetonitrile, urea, guanidine, hydrazine (✅) | propionic/butyric acid, acetone, naphthol, cresol |
| 桥环 | adamantane (作为取代基) (⚠️) | norbornane |

### 最关键缺失 (按优先级)

1. **phenanthrene** — 三大萘/蒽/菲仅缺此
2. **purine** — 生物化学核心杂环
3. **carbazole** — 材料化学重要母体
4. **triazole / tetrazole** — 点击化学核心
5. **azetidine** — 四元含氮杂环
6. **chromene / xanthene** — 染料核心
7. **cinnoline / phthalazine** — 完备6+6二嗪系列
8. **norbornane** — 经典桥环

---

## 11. 架构建议

参考 substituents 的集中注册表方案，为母体也建立类似结构：

```python
# src/namepredict/layer2/retained_parents.py  
@dataclass(frozen=True)
class RetainedParent:
    id: str
    en: str
    zh: str
    naming_class: str       # carbocycle | fused56 | naph_family | monohetero | ...
    n_rings: int
    ring_type: str          # carbo | hetero
    retained_level: str     # PIN | general | deprecated
    numbering_labels: tuple[str, ...]  # 1…8a etc.
    rule_ref: str           # P-22.1.3 etc.
```

当前项目的 `ScaffoldSpec` 已经接近这个设计，主要差距是：
- 三环及以上（anthracene 之外）尚未有 Spec
- 功能母体（formic acid, phenol 等）尚未纳入 Spec 体系
- 编号策略（NumberingPolicy）对三环+体系需要扩展
