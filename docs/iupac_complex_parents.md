# IUPAC 2013 复杂母体（complex parent structures）全面盘点

> 研究方式：deep-research 工作流（5 个搜索角度 → 20 个来源 → 65 条声明 → 25 条经 3 票对抗核验 → 23 通过 / 2 否决）。
> 相关文档：`iupac_parent_hydrides.md`（§8c 酮保留母体）、`docs/iupac/`（P-2/P-4/P-44 等章节）。
>
> **问题界定**：本文盘点的"复杂母体"指——
> 1. 不是一个单一官能团（如羧酸、醇）；
> 2. 也不是一个裸的环/环系（如苯、萘）；
> 3. 可以是"环 + 环外挂接部分"（如苯醌 = 环 + 环外 =O），也可以是 IUPAC 的特殊规定（retained names / trivial retained names）。

---

## 0. 核心概念界定（判定边界）

**IUPAC 对"functional parent"（功能母体）有明确排除条款**（Gold Book F02556，出处 R-0.2.1.2 / 蓝书 p.13）：

> "A parent hydride bearing a characteristic group denoted by a suffix, for example, **cyclohexanol**, is not considered to be a functional parent, but may be described as a **functionalized parent hydride**."

即 **"环 + 单一后缀官能团"（cyclohexanol = 环 + -ol）不算复杂母体**，只能叫 functionalized parent hydride。这正是问题界定 (1) 的判定线。

**保留名（retained names）正式定义**（P-15.1.8.2；2004 草案编号 P-55）：
> 为有机命名而保留的俗名/半系统名（trivial or semisystematic names retained for naming organic compounds），按可取代性分三类：
> - **Type 1**：无限取代
> - **Type 2a**：仅限前缀取代，以承认 functional parent 明示/暗示的官能团（如醌的 =O）
> - **Type 2b**：强制前缀；**Type 2c**：其他；**Type 3**：不允许取代

复杂母体的合法性来源就是这套保留名框架：名称内蕴的 =O / 内酯被当作 functional parent 的明示/暗示官能团，允许前缀取代。

**结论：复杂母体 = "环/环系 + 环外悬挂 =O"或"名称内蕴官能团的保留名"，不是"环 + 单个后缀官能团"。**

---

## 1. 五类来源总览

| # | 类型 | 代表母体 | 章节 | 保留名? | 后缀/命名机制 |
|---|------|---------|------|---------|--------------|
| 2 | 醌类 | benzoquinone, naphthoquinone, anthraquinone | P-64.2.2.2.3（+ P-22.2） | ✅ Type 2a，**非 PIN** | PIN 一律退回 `-dione` |
| 3 | 环酮类 | cyclohexanone, camphor | P-64.2（P-26.2.2.2） | 环烷酮非保留（modified parent hydride） | 系统 `-one` 后缀；camphor 保留名 |
| 4 | "杂环 + 环外 =O" | pyran-4-one, 2H-chromen-2-one, isochromen-1-one, xanthen-9-one | P-64.2 / 1993 表 23(a) | coumarin 等为 trivial（保留）类名 | `-one` 后缀 + 保留名 |
| 5 | 内酯类 | oxolan-2-one, oxiran-2-one, oxepan-2-one | P-65.6.3.5 + Gold Book L03439 | 非保留（系统名） | 1-oxacycloalkan-2-one 核心 |
| 6 | 黄酮 functional parents | flavone, isoflavone, neoflavone, flavanone | 2017 IUPAC/IUBMB 黄酮建议 | 保留 functional parent | 名称内蕴一个羰基；多羰基退回 `-dione/-trione` |
| 7 | 天然产物保留母体 | androstane, cholestane, morphinan 等 | 附录 3、表 10.1、P-101.2.2/7 | ✅ 整体保留 | 叠 `-one/-ol/-oic acid` 后缀 |

> 章节定位修正：酮类全部集中在 **P-64**（"Ketones, pseudoketones and heterones, and chalcogen analogues"，位于 P6.pdf 即 P-60~P-65），**不在 P-66**。内酯（环酯）在 **P-65.6.3.5**。P-66（酰胺/腙/腈/醛）、P-67（oxoacids）、P-68、P-69 均不覆盖这类母体。
> 另：2004 草案编号（P-55.4.2 等）在 2013 终稿中重编号（P-64.2.2.2.3 等），系位置迁移而非内容变化。

---

## 2. 醌类：环 + 两个环外 =O（P-64.2.2.2.3）

**核心规则**：苯醌/萘醌/蒽醌三个保留名是 **Type 2a 保留名，但一律不是 PIN**。P-64.2.2.2.3 明文：
> "No retained quinone names are used as preferred IUPAC names" —— 仅"供一般命名（general nomenclature）用"。

**系统 PIN 用 `-dione` 后缀加在母环上**：

| 保留名（general） | 系统 PIN | CAS |
|-------------------|----------|-----|
| 1,4-benzoquinone | **cyclohexa-2,5-diene-1,4-dione** | 106-51-4 |
| 1,4-naphthoquinone | **naphthalene-1,4-dione** | 130-15-4 |
| 9,10-anthraquinone | **anthracene-9,10-dione** | 84-65-1 |
| — | quinoline-5,8-dione | — |
| — | acenaphthylene-1,2-dione | — |

- 无位次的裸词 **`benzoquinone` / `anthraquinone` 被明确禁用**（"not benzoquinone"）。
- 带取代时的 PIN 反例：`2-methylanthracene-9,10-dione`，**而非** `2-methyl-9,10-anthraquinone`。
- 历史注记：2004 草案曾允许 quinone 型作 PIN，终稿收窄。
- 1979 规则（C-311.1 / C-317.1）另设专属 `-quinone` / `-diquinone` 后缀（区别于普通 `-one`），是 benzoquinone、camphorquinone 等俗名的来源。

> 与苯/蒽等"裸环"母体（PIN 即保留名）形成对比：裸环 PIN=保留名；带环外 =O 的醌 PIN=系统 `-dione`。

---

## 3. 环酮类（P-64.2）

### 3a. 环烷酮 — 非保留，系统 `-one`
- 环烷酮按 `母环 + -one` 后缀命名，**cyclohexanone 本身就是 PIN**，归类为 P-26.2.2.2 的 **"modified parent hydride names"**（修正母氢化物名），不是裸保留名。
- 开链酮同规则：C-312.1 的 `-one/-dione` 加在母烃名上；**acetone 是保留名**（P-66.5.1.1 仍保留 propan-2-one/acetone）；酸名派生俗名 butyrone/valerone/stearone 被明确劝阻。

### 3b. 带环组分的保留酮名（1979 C-313.2(d) 装置）
专属保留名装置：把对应酸名的 `-ic acid`/`-oic acid` 词尾改为 `-ophenone`（苯环）或 `-onaphthone`（萘环），得 **acetophenone** 等。C-313 列出的保留俗名：

| 保留名 | 中文 | 2013 级别 | 说明 |
|--------|------|----------|------|
| **acetophenone** | 苯乙酮 | ⚠️ General only | PIN = `1-phenylethan-1-one`；1993 Guide R-9.1 表 27(a) 亦保留 |
| propiophenone | 苯丙酮 | （拼写非 'propionophenone'） | — |
| **chalcone** | 查尔酮 | ✅ **唯一酮 PIN** | 优先于 cinnamophenone / benzylideneacetophenone / 3-phenylacrylophenone；仅反式 |
| deoxybenzoin | 脱氧苯偶姻 | — | — |

> ⚠️ **未决细节**：acetophenone / benzophenone 的确切 substitutability 类型（是否"不可取代"）在核验中被否决（1-2），未能确认。现有 `iupac_parent_hydrides.md` §8c 表格写 acetophenone"不可被取代"，与 `benzil / biacetyl / propiophenone 在 2013 明确废除（P-64.2.1.3）"一致，但该点请以 P6.pdf 原文复核为准。

### 3c. 樟脑（camphor）
camphor 为保留俗名；系统名 **`1,7,7-trimethylbicyclo[2.2.1]heptan-2-one`**（用 `-one` 后缀），属"稠合环 + 环外 =O + 甲基取代"的复杂母体。

---

## 4. "杂环 + 环外 =O"母体（pyran-4-one 家族）

1993 IUPAC Guide 表 23(a) 给出杂环父氢化物优先级顺序（高位含 Pyran、Isobenzofuran、Isochromene、Chromene、Xanthene 等 O-杂环），多杂环化合物选最高优先级环为母体。这些 O-杂环母体叠 `-one` 后缀即得"杂环 + 环外 =O"复杂母体：

| 母体 | 结构 | 旧称 | CAS | 备注 |
|------|------|------|-----|------|
| **pyran-4-one**（= 4H-pyran-4-one） | 吡喃-4-酮 | pyrone | — | PIN |
| **2H-chromen-2-one** | 苯并吡喃-2-酮（δ-内酯稠苯） | **coumarin** / 1,2-benzopyrone | 91-64-5 | trivial（保留）类名 |
| **isochromen-1-one** | 异香豆素 | isocoumarin | — | PIN 名 1H-isochromen-1-one |
| **xanthen-9-one** | 占吨酮 | xanthone | — | PIN 名 9H-xanthen-9-one |

**coumarin 的定义**（Gold Book C01369，来源 PAC 1995, 67, 1307）：C₉H₆O₂，δ-内酯稠合苯环，C2 环碳带环外 C=O 的双环 benzopyran-2-one，系 trivial（保留）类名。InChI 编码确认 C=O 在 c10-9、O11 桥连 C9-C8、苯环稠合吡喃酮——正是"环 + 环外 =O"复杂母体的教科书实例。

> 注意：表 23(a) 的优先级标准在 2013 蓝书中被 **P-44 取代**，但环母体名本身仍有效。

---

## 5. 内酯类（P-65.6.3.5 + Gold Book L03439）

**定义**（Gold Book L03439）：内酯（lactone）= 羟基羧酸的环状酯，核心结构为 **`1-oxacycloalkan-2-one`**（杂环 2 位环碳带环外 =O）。

- 明示允许环内含双键（不饱和）或杂原子替代环碳的类似物仍属内酯。
- 覆盖范围很宽：oxepan-2-one（七元）、oxolan-2-one（五元）、**oxiran-2-one（3 元 α-内酯）**、氮杂双环内酯等。
- **内酯不是保留名**，全部为系统 `-one` 命名。
- 常见俗名对照：butyrolactone → oxolan-2-one；valerolactone → oxan-2-one。

> 内酯与醌/香豆素的区别：内酯的环外 =O 来自酯基（环内已有酯 O），醌/香豆素是"环 + 独立环外 =O"或"δ-内酯稠环"。判定边界仍以 Gold Book 定义为准。

---

## 6. 黄酮 functional parents（2017 IUPAC/IUBMB 建议，PAC 90(9):1429-1486）

**flavone / isoflavone / neoflavone 是 functional parents**：其名**已内蕴一个羰基**（"ring + exocyclic =O"），无需显式叠 `-one`。当结构含**多个**羰基时，退回以母氢化物（flavan / isoflavan / neoflavan）+ `-dione/-trione` 系统命名。

| 名称 | 结构 | CAS |
|------|------|-----|
| **flavone** | 2-phenyl-4H-1-benzopyran-4-one（= 2-phenyl-4H-chromen-4-one） | 525-82-6 |
| **flavanone** | flavan-4-one = 2-phenyl-2,3-dihydro-4H-1-benzopyran-4-one（= 2-phenylchroman-4-one） | 487-26-3 |
| isoflavone | 3-苯基位置异构 | — |
| neoflavone | 4-苯基位置异构 | — |

> 这是"名称内蕴官能团 + 环外 =O"的典型：flavone 的名字本身就是"母体 + 一个隐含羰基"，比 `-one` 后缀更"高阶"——多羰基才退回系统 `-dione`。来源为 2017 黄酮建议（独立现行文件），非 2013 蓝书正文。

---

## 7. 天然产物保留母体（附录 3、表 10.1、P-101.2.2/7）

蓝书**附录 3 "Structures for Alkaloids, Steroids, Terpenoids, and Similar Compounds"** 提供甾体/生物碱/萜类母体环系。保留理由（P-101.2.2 / 2004 草案 P-55.4.6.1）：**隐含立体构型、广泛使用、结构复杂**。

### 7a. 甾体母体
androstane、cholestane、estrane、gonane、**pregnane**（非 'pregnaane'，核验发现的抄录笔误）、bufanolide、cardanolide、spirostan、stigmastane（另含 cholane、ergostane、lanostane 等）。以 5α/5β 立体化学定义（IUPAC-IUB 甾体命名 1989）。

### 7b. 生物碱母体
morphinan 等约 70 个母体。

### 7c. 后缀用法（"复杂母体 + 官能化后缀"典型）
- **androstan-3-one**（甾酮）
- **cholest-5-en-3β-ol**（胆固醇）
- **estra-1,3,5(10)-triene**（雌三烯）

> 这类母体本质是"带隐含立体构型的多环骨架 + 叠官能后缀"，与醌/内酯的"环 + 环外 =O"不同——复杂在骨架本身是保留名，而非悬挂部分。

---

## 8. 判定边界总结

| 结构 | 分类 | 原因 |
|------|------|------|
| cyclohexanol（环 + -ol） | functionalized parent hydride，**非** complex parent | F02556 明文排除 |
| cyclohexanone（环 + -one） | modified parent hydride（PIN） | 系统 `-one` 后缀，非保留 |
| 1,4-benzoquinone | **complex parent**（保留名 Type 2a，general） | 名称内蕴环外 =O |
| coumarin（2H-chromen-2-one） | **complex parent**（trivial 保留类名） | δ-内酯 + 环外 =O |
| flavone | **complex parent**（functional parent） | 名称内蕴一个羰基 |
| androstane | **complex parent**（整体保留） | 隐含立体构型 + 广泛使用 |

---

## 9. 与 chem 项目的对应关系

| 复杂母体 | 项目现状 | 实现文件 |
|----------|---------|---------|
| 1,4-benzoquinone / 1,2-benzoquinone / 9,10-anthraquinone | ✅ 已实现（专用检测器） | `layer2/scaffold/benzoquinone.py` |
| chromenone（coumarin） | ✅ 已实现（专用检测器） | `layer2/scaffold/chromenone.py` |
| 模板表 `_TEMPLATES`（24 个纯母环） | ✅ | `layer2/scaffold/retained_templates.py` |

**本项目设计佐证了研究结论**：模板表刻意只匹配环系统原子集合，不向外延申——因为真实母环普遍带取代基。碳羰母环（苯醌/香豆素）含环外 =O 是超集，进不了纯母环模板表，所以必须走专用检测器（`benzoquinone.py` / `chromenone.py`）。这与 IUPAC 的定位完全一致：这些就是需要单独处理的"复杂母体"。

**研究发现的扩展候选**（尚未覆盖，按 IUPAC 归类）：
- 醌系扩展：1,4-naphthoquinone、phenanthrenequinone、acenaphthylene-1,2-dione（PIN）
- "杂环 + =O"系：isocoumarin（isochromen-1-one）、xanthone（xanthen-9-one）、flavone/flavanone
- 内酯系：oxolan-2-one（butyrolactone）等
- 保留酮：camphor（1,7,7-trimethylbicyclo[2.2.1]heptan-2-one）

---

## 10. 局限与未决问题

1. **网络限制**：核验阶段 iupac.qmul.ac.uk / goldbook.iupac.org / media.iupac.org / acdlabs.com 的直接抓取被沙箱网络策略阻断，部分引用依赖搜索引擎索引快照；但每条声明均经 2-3 个独立来源交叉印证。
2. **P-64.2.1.2 保留酮名精确清单**：声称清单"仅含 acetone 与带位次苯醌/萘醌/蒽醌，acetophenone/benzophenone 无取代保留"的细化断言被 3 票否决（1-2），各词 substitutability 类型未澄清，需直接核对 P6.pdf 原文。
3. **"杂环 + 环外 =O"完整枚举**：除已覆盖的 pyran/chromene/isochromene/xanthene/oxiran-2-one 外，quinolone、chromen-4-one 本身、氮杂内酯、五/七元内酯、fur-2(3H)-one 等在 P-64/P-22 下的完整枚举未系统建立。
4. **2017 黄酮建议与 2013 蓝书 P-44/P-22.2 的逐条对应**未系统建立。
5. **NCI/CAS/药物骨架全面表格**：本次仅确认约 7 个化合物的 CAS 与数据库名（coumarin 91-64-5、flavanone 487-26-3、flavone 525-82-6、1,4-benzoquinone 106-51-4、1,4-naphthoquinone 130-15-4、9,10-anthraquinone 84-65-1、acetophenone 98-86-2），缺一张覆盖全部母体的系统表。
6. **'pregnaane' 笔误**：正确保留名为 'pregnane'（无第二个 a）。

---

## 引用来源

**一级来源（primary）**：
- IUPAC Blue Book 2013（QMUL 官方镜像）：https://iupac.qmul.ac.uk/BlueBook/
  - P-6 章（酮/醌）：https://iupac.qmul.ac.uk/BlueBook/PDF/P6.pdf （P-64, P-64.2.2.2.3）
  - 附录 3（甾体/生物碱/萜类）：https://iupac.qmul.ac.uk/BlueBook/Papp3.html
- 2004 草案 Chapter 5（P-55 保留名分类）：http://media.iupac.org/reports/provisional/abstract04/BB-prs310305/Chapter5.pdf
- IUPAC Guide 1993 表 23(a)（杂环优先级）：https://www.acdlabs.com/iupac/nomenclature/93/r93_673.htm
- IUPAC 1979 C-311.1 / C-312.1 / C-313.2(d)（-one/-quinone/酮保留名装置）：https://www.acdlabs.com/iupac/nomenclature/79/r79_270.htm
- IUPAC/IUBMB 黄酮命名（2017, PAC 90(9):1429-1486）：https://iupac.qmul.ac.uk/flavonoid/index.html
- Gold Book F02556（functional parent）：https://goldbook.iupac.org/terms/view/F02556
- Gold Book C01369（coumarin）：https://www.dev.goldbook.iupac.org/terms/view/C01369
- Gold Book L03439（lactone）：https://goldbook.iupac.org/terms/view/L03439

**次级来源（secondary）**：Peeref review（pyran 结构单元综述）、IDRBLab INTEDE（9,10-anthraquinone 数据库条目）。
