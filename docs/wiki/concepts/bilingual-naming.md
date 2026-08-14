# 中英双语命名约定 (Bilingual Naming Convention)

> **概念页面** | **最后更新:** 2026-07-20

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

整个 Layer5 组装管线中，所有名称生成函数遵循一个统一的返回值约定：**返回 `tuple[str, str] | None`，其中第一个元素为英文名、第二个元素为中文名**。这一定约定确保了函数的可组合性——上游的输出可以直接作为下游的输入，无需额外转换。

### 管线的 5 步双语流水线

以 `assemble()` 函数（`src/namepredict/layer5/assembler.py:573-583`）为例，每一步都操作 `(en, zh)` 对：

```
1. _names_for(kind, n, numbered)     → (en, zh) | None   母体名称
2. _prefix_for(numbered, kind, n)    → (pre_en, pre_zh)   取代基前缀
3. join_kind_name(kind, pre, names)  → (en, zh)           前缀+母体拼接
4. maybe_anion_names(numbered, en, zh) → (en, zh)         羧酸根转换
5. apply_rs_prefix(numbered, en, zh) → (en, zh)           R/S 立体化学
6. maybe_metal_salt_names(numbered, en, zh) → (en, zh)    盐后缀
```

每一步都是纯文本转换：输入 `(en, zh)`，输出 `(en, zh)`。这种设计具有三个优点：

1. **可测试性**：每个步骤可以独立验证其中英双语输出的正确性。
2. **可追踪性**：当名称出现错误时，可以精确定位到出错的步骤（例如中文盐后缀错位但英文正确）。
3. **可扩展性**：新增官能团命名时，只需实现一个新的 `(kind, n, numbered) -> (en, zh) | None` 函数并注册到 `_names_for` 的派发链中。

### 母体名称派发中的双语约定

`_names_for` 函数（`src/namepredict/layer5/assembler.py:79-112`）先查 `chain_engine._KIND_TABLE`（15 个 `_Chain` spec，按 `scaffold_id` 运行时注入环前缀/稠环词干），命中即返回；否则落到特殊 worker（`_exocyclic_acid_names`/`phenyl`/`benzene`/`_ARENE_RETAINED_NAMES`/`benzenediol`）或 `_parent_stem_names` 回退。每个 worker 都是返回 `(en, zh)` 或 `None` 的双语函数——返回 `None` 时派发器降级到下一策略。这种"尝试-失败-降级"模式在双语的上下文中尤其重要——如果一个命名策略返回了英文名但无法生成中文名（或反之），则整个结果不被接受。

## 词干表：双语单词的单一权威来源

`src/namepredict/layer5/stems.py` 是整个 NamePredict 双语命名的数据基础。它定义了从 C1 到 C35+ 碳数范围的英中双语词干生成器，确保整个系统中同一碳数的烷烃词干、官能团名称保持一致性。

### 基表结构

**C1-C10 保留名基表**（`src/namepredict/layer5/stems.py:9-17`）：

```python
_ALKANE_EN_BASE = {1: "methane", 2: "ethane", ..., 10: "decane"}
_ALKANE_ZH_BASE = {1: "甲烷",  2: "乙烷",   ..., 10: "癸烷"}
```

英文和中文基表严格一一对应，使用相同的整数 key。这种 `{n: str}` 的字典结构确保了查找效率和数据的不可变性。

**C11-C19 半系统命名**（`src/namepredict/layer5/stems.py:18-21`）：

英文使用复合词干 `_SEMI_EN`（`undec`, `dodec`, `tridec`...），中文则由 `zh_num(n)` 函数动态生成 `十一烷`、`十二烷`...。

**C20+ 乘性组合**（`src/namepredict/layer5/stems.py:53-63`）：

英文词干通过 `_compose_en_stem(n)` 按十位（`icos`, `triacont`, `tetracont`...）+ 个位（`hen`, `do`, `tri`...）组合生成，如 C22 = `docos`（do + cos，elide icos 的 i）。中文仍由 `zh_num(n)` 生成数字+烷。

### 官能团词干的派生

每种官能团都有独立的一对英/中函数，从烷烃词干派生：

| 官能团 | 英文函数 | 中文函数 | 示例 (C2) | 示例 (C2 中文) |
|--------|---------|---------|-----------|---------------|
| 烷烃 | `alkane_en(n)` | `alkane_zh(n)` | `ethane` | `乙烷` |
| 醇 | `alcohol_en(n)` | `alcohol_zh(n)` | `ethanol` | `乙醇` |
| 酸 | `acid_en(n)` | `acid_zh(n)` | `acetic acid` | `乙酸` |
| 醛 | `aldehyde_en(n)` | `aldehyde_zh(n)` | `acetaldehyde` | `乙醛` |
| 酰胺 | `amide_en(n)` | `amide_zh(n)` | `acetamide` | `乙酰胺` |
| 腈 | `nitrile_en(n)` | `nitrile_zh(n)` | `acetonitrile` | `乙腈` |
| 酯酰基 | `ester_acyl_en(n)` | — | `acetate` | — |

所有公开词干表通过 `_fill(fn, lo=1, hi=35)` 自动填充为完整字典（`src/namepredict/layer5/stems.py:333-352`），从 C1 覆盖到 C35，确保长链名称无需手动维护。中文的 C1/C2 酸、醛、酰胺、腈等使用保留名（如 `甲酸`/`乙酸`、`甲醛`/`乙醛`），与英文保留名（`formic acid`/`acetic acid`、`formaldehyde`/`acetaldehyde`）保持对应。

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
| **盐命名** | `sodium acetate` (阳离子在前) | `乙酸钠` (阳离子在末尾) |
| **盐酸盐** | `;hydrochloride` | `;盐酸盐` |
| **立体化学** | `(E,2S)-` | `(E,2S)-` (与英文一致，不翻译) |
| **环前缀** | `cyclohexane` | `环己烷` |
| **编号前缀** | `1H-pyrrole-` | `1H-吡咯`（`1H-` 保留，母体翻译） |

### 酯命名的语序反转

酯命名是中英文差异最显著的例子。英文的酯命名为 `{alkyl} {acyl}ate` 格式，烷基（醇部分）在酰基（酸部分）之前：

- EN: `methyl benzoate`（甲基 + 苯甲酸酯）

中文则遵循"酸在前、醇在后"的规则，格式为 `{酸}{醇}酯`：

- ZH: `苯甲酸甲酯`（苯甲酸 + 甲基 + 酯）

这一差异在代码中由 `benzene_names.py` 的 `benzoate_parent_names` 和 `join_ester_name` 函数处理，分别构造英文和中文的语序。对于不饱和酯和带取代基的酯，中英文的插入位置也有不同的处理逻辑。

### 盐命名的后缀位置反转

盐的命名也体现出语序差异。英文将金属阳离子放在羧酸根名称之前作为前缀（`sodium dodecanoate`），而中文将金属名放在名称末尾（`十二酸钠`）。这一转换由 `stems.py` 中的 `maybe_metal_salt_names` 函数（`src/namepredict/layer5/stems.py:285-290`）完成：

- `_salt_en`：将金属名作为前缀拼接到英文名（仅当英文名以 `ate` 结尾时）
- `_salt_zh`：将金属中文名替换中文名的末尾 `酸根` → `酸{金属}`

盐酸盐（HCl salt）采用不同的机制：通过 `acid_salt` 字段追加 `;hydrochloride` / `;盐酸盐` 后缀。

### 取代基前缀的排序一致性

虽然英文前缀按取代基英文名的字母序排列（符合 IUPAC P-14.5），但中文前缀的排列仍使用 `alkyl_alpha_key`——一个基于英文名的排序键，定义在 `assembler_prefixes.py` 中并由 Layer3 的取代基提取器共享。这确保了在双语输出中，取代基的排列顺序始终保持一致：中英文的前缀中的基团顺序总是相同的，避免了因语种不同导致取代基排列顺序不一致的混淆。

## `zh_stem` 转换

`zh_stem` 函数（`src/namepredict/layer5/stems.py:45-50`）是一个关键的中文词干提取工具。它从中文全名中剥离末端的官能团/母体后缀，返回"裸词干"供后续拼接。

**后缀剥离表**（`src/namepredict/layer5/stems.py:31`）：

```python
_ZH_SUFFIXES = ("酰胺", "酰氯", "硫醇", "烷", "醇", "酸", "醛", "腈", "胺", "酮", "烯", "炔")
```

**转换示例：**

| 输入 (`zh_full`) | 输出 (词干) | 用途 |
|-----------------|-----------|------|
| `十一烷` | `十一` | 构造 FG 位次名称，如 `十一烷-2-醇` |
| `乙醇` | `乙` | 构造带位次的醇名，如 `乙-1,2-二醇` |
| `丁酰胺` | `丁` | 构造取代酰胺名 |
| `环己烷` | `环己` (无匹配后缀，原样返回) | 环母体的词干 |

`zh_stem` 的工作原理是按 `_ZH_SUFFIXES` 元组中的顺序依次尝试匹配——一旦某个后缀匹配成功，就切除该后缀并返回剩余部分。匹配顺序从长后缀（`酰胺`）到短后缀（`烯`、`炔`），避免短后缀误匹配长后缀的尾部。如果没有任何后缀匹配，函数返回原始输入——这意味着它可以安全地用于任何中文名称字符串。

在 `chain_engine.py` 的环系命名中广泛使用了 `zh_stem`——`_chain_names` 对 `scaffold_id=="carbocycle"` 的母体注入 `cyclic=True`（恒加 `环`/`cyclo` 前缀），并用 `zh_stem` 获取环烷烃词干后再追加烯/醇等后缀。

## `zh_num`：中文数字生成

`zh_num(n)`（`src/namepredict/layer5/stems.py:34-42`）将整数 n（1-99）转换为中文基数词：

- 1-9：直接使用 `_DIGIT_ZH` 字符串索引（`一` `二` `三`...`九`）
- 10-19：`十` + 个位（`十一` `十二`...而非 `一十` `二十`...）
- 20-99：十位 + `十` + 个位（`二十五`、`九十九`）

该函数是中文 C11+ 烷烃名称（`十一烷`、`三十五烷` 等）的数据源，同时也是 C20+ 乘性组合名称中中文部分的唯一来源。`zh_num` 不存在对应的英文函数，因为英文 C11+ 使用半系统词干表 `_SEMI_EN` 而非数字前缀。

## `name_mode` 的双语传播

`name_mode` 参数控制命名风格的选择（`"general"` 通用名 vs `"pin"` 中国药品通用名称），它从 `SMILESNNamer` 构造函数开始，经过整个管线传播到 Layer3 的取代基命名和 Layer5 的 FG 命名模块。

**传播路径：**

```
SMILESNNamer(name_mode="general")         # namer.py:222
  → _pipeline(smiles, t0, name_mode)       # namer.py:203
    → _name_mol(mol, ..., name_mode)       # namer.py:187
      → _run_candidates(info, ..., name_mode)  # namer.py:145
        → _prepare_candidate(info, parent, name_mode)  # namer.py:112
          → extract_substituents(info, parent, name_mode)  # Layer3
        → _assemble_candidate(parent, subst, ..., name_mode)  # namer.py:104
          → _ok_result(numbered, ..., name_mode)  # namer.py:68
            → numbered["name_mode"] = name_mode
```

在 Layer5 中，命名模块可以通过 `numbered.get("name_mode")` 或 `parent.get("name_mode")` 获取当前模式，从而在通用名和 PIN 名之间切换。当前实现以通用名为主要支持目标，`name_mode` 机制为未来的命名风格扩展预留了接口。

## 盐元数据的双语结构

Layer0 的 `dissociate_salt` 函数（`src/namepredict/layer0/salt.py:90-96`）在检测到盐结构时，返回的元数据字典包含完整的中英双语字段：

| 字段 | 示例值 | 说明 |
|------|--------|------|
| `metal` | `"sodium"` | 英文金属名 |
| `metal_zh` | `"钠"` | 中文金属名 |
| `n_metal` | `1` | 金属阳离子数量 |
| `acid_salt` | `"hydrochloride"` | 英文酸盐后缀 |
| `acid_salt_zh` | `"盐酸盐"` | 中文酸盐后缀 |

这些元数据在 `_name_mol` 中注入 `NameResult.meta["salt"]`（`src/namepredict/namer.py:198-199`），供 Layer5 组装器在生成中英双语名称时追加盐后缀。由于盐解离在 Layer1 之前完成，Layer1-5 处理的始终是解离后的纯有机片段，各层代码无需关心盐的存在。

## 设计原则总结

NamePredict 的双语命名设计遵循以下原则：

1. **同步而非翻译**：英文和中文名称从同一结构化数据源（`numbered` 字典）同时生成，而非先生成一种语言再翻译为另一种。这避免了翻译过程中引入的二次歧义。

2. **词干表为单一权威来源**：`stems.py` 是双语单词的唯一生成点。任何碳数的烷烃和官能团名称都从此处派生，杜绝了不同模块使用不同翻译的可能性。

3. **成对函数约定**：英中双语始终以成对形式出现——表（`_ALKANE_EN_BASE` / `_ALKANE_ZH_BASE`）、函数（`alcohol_en(n)` / `alcohol_zh(n)`）、元组（`(en, zh)`）。不存在"只有英文"或"只有中文"的命名代码路径。

4. **差异显式化**：中英文命名规则的差异（如酯的语序、盐的后缀位置）在对应模块中显式处理，而非隐藏在通用逻辑中。这使得特定语种的规则变更不会意外影响另一种语言。

5. **空值安全**：所有双语返回函数在任一语言无法生成时返回 `None`（而非部分成功）。客户端代码（如 `_names_for` 的派发链）在返回 `None` 时自动降级到备选策略，确保管线始终产出完整的中英双语名称。

---

## 相关页面

- [[architecture/overview]] — 6 层架构总览，展示双语数据在各层之间的流转
- [[architecture/layer5-name-assembly]] — Layer5 名称组装层的完整文档，包含 `(en, zh)` 流水线的实现细节
- [[architecture/layer0-preprocessor]] — Layer0 预处理层，盐元数据的双语字段定义
- [[concepts/iupac-rules]] — IUPAC 蓝皮书规则映射（P-14, P-21, P-63, P-65）
- [[reference/core-data-contracts]] — `NameResult` 及 `numbered` 字典的双语字段结构
- [[reference/api]] — `SMILESNNamer.name()` 中 `name_mode` 参数的公开接口
- [[reference/molecule-model]] — `NameResult(en, zh, ...)` 的类型定义
- [[index]] — Wiki 首页
