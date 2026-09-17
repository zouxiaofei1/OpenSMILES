# Layer5: 名称组装 (Name Assembly)

> 源文件 7 个（约 2146 行） · 对外接口 `assemble`

## 概述

Layer5 是命名流水线的末层，把 Layer4 输出的编号结果（numbered dict）渲染为双语 IUPAC 名称。本层不做化学推断，只做名称拼接：官能团词尾、取代基前缀、立体描述符、盐后缀都在此定形。

### 输入 numbered dict

| 键 | 含义 | 本层用途 |
| --- | --- | --- |
| `parent.kind` | 母体类别 | 选 `_KIND_TABLE` 条目 |
| `parent.n_carbons` | 母体碳数 | 词干生成 |
| `parent.stem_en` / `parent.stem_zh` | 环式/稠环词干 | 覆盖碳数词干 |
| `parent.scaffold_id` | 骨架 ID | `benzene` / `carbocycle` / 稠环 |
| `parent.chain` | 母体原子序列 | 指示氢位次映射 |
| `parent.mol` | RDKit 分子 | 臂根成环判定与手性 |
| `parent.numbering_scaffold.labels` | 原子 → 位次标签 | 位次串渲染 |
| `parent.oxo_kind` / `parent.n_oh` | 含氧酸中心元素与酸式氢数 | 词尾查表 |
| `parent.thio_side` | 硫代羧酸 S-酯标记 | 走 thioate 词尾与 `S-` 斜体 |
| `parent.radical_anchor_element` | 杂原子锚点自由基标记 | 走单核自由基路径 |
| `parent.hydro_prefix` | 加氢前缀 | 拼到母体名前 |
| `parent.indicated_h_locants` | 指示氢位次 | 拼 `nH-` 前缀 |
| `parent.salt_meta` / `parent.n_om` | 金属盐与游离根计数 | 盐式与「根」后缀 |
| `substituents` | 取代基记录 | `o_side`、`paren`、`attach_idx`、`atoms`、`arm_fenced` |
| `fg_locants` | 主官能团位次记录 | `kind` / `locants` / `omit` |

### 输出 NameResult

`NameResult(en, zh, success, source, time_ms, meta)`。失败时经 `_fail` 构造，`success=False` 且 `meta` 携带 `reason`、`n_carbons`、`kind`；成功时经 `_ok` 构造。`_unsupported` 是带上述三项诊断的 `_fail` 快捷封装。

### assemble 步骤

`assemble(numbered, *, time_ms=0.0, source="iupac")` 按固定顺序执行：

1. `_ensure_fused_stem` — 稠环词干缺失时由 `fused_namer.fused_parent_names` 组装注入，并前置 L4 指示氢位次；失败即返回 unsupported。
2. `_names_for(kind, n, numbered)` — 取母体双语名。
3. `join_hydro_prefix` — 拼接 hydro 前缀与指示氢前缀。
4. `join_ring_cation_suffix` — 环内 N⁺/O⁺ 缀 `-{位次}-ium`（P-62.4.1）。
5. `join_kind_name(kind, _prefix_for(numbered, kind, n), names, numbered)` — 取代基前缀与母体拼接。
6. `join_anion_names` — 羧酸阴离子转 `-ate` / 酸根。
7. `join_ez_prefix` — 母体外挂双键的 E/Z 部件。
8. `join_rs_prefix` — 母体手性中心的 CIP R/S 部件。

`_parent_stem_names` 是词干直取通道，仅返回 `parent.stem_en` / `parent.stem_zh`，不做任何后缀拼接。

## 链引擎 `chain_engine.py`

链引擎以数据驱动方式生成链状母体名：一张 `_KIND_TABLE` 描述每个 kind 的词尾、位次规则与特例，由 `_chain_names` 统一消费。

### `_Chain` 规格

`_Chain` 是 `frozen dataclass`，字段分四类。

- 词尾：`kind`、`en_suf`、`zh_suf`、`coda`（饱和词干后接段，默认 `"an"`，部分 kind 用 `""` 或 `"ane"`）。
- 位次：`no_loc`（`"plain"` 无位次回落普通名 / `"none"` 直接返回 None）、`omit_rule`（`(n, loc, omit) -> bool`，True 表示省略位次）、`fg`、`need`、`ene_loc_omit`、`yne_loc_omit`、`yl_loc_omit`、`zh_loc_omit`。
- 不饱和段：`ene_seg`、`yne_seg`（段式单烯/炔段）、`ene_base`、`yne_suf`（融合式基座）、`unsat_polyol`、`zh_full`、`ez_ene`、`ez_ene_multi`。
- 覆盖与派生：`plain_maps`（俗名表）、`plain_fn`（按碳数派生命名）、`plain_hook`（按分子计数直接给词尾）、`variant`（`{scaffold 或 multiplicity: 字段覆盖 dict}`）、`stem`、`aromatic`、`cyclic`、`cyclic_unsat`、`wrap`。

三个位次规则常量：`_NO_OMIT`（恒不省）、`_omit_all`（恒省）、`_omit_term_locant`（`omit` 或无位次，或 `n <= 2` 且位次为 1 时省）。

### 主流程 `_chain_names`

1. `plain_hook` 非空时先调用；返回非 None 即直接作为结果（无碳链母体，如磷酸等 P 中心）。
2. 芳香环醇统一改 `zh_suf="酚"`。
3. `_exo_ring_spec` 改写环外主基的 spec。
4. 按 `_parent_multiplicity`（优先 `principal_expression_facts.multiplicity`，回落 `principal_group_count`）走多 FG 或单 FG 分支：多 FG 经 `_generated_mult_fields` 生成数量后缀字段并叠加 `variant[mult]` 覆盖；单 FG 关闭 `unsat_polyol` / `zh_full` 后叠加 `variant[1]`。
5. `_chain_enyne` 尝试不饱和段渲染；`unsat_polyol` 时再补 FG 位次串，且要求 `fg_locants` 记录位数等于 `need`。
6. 不饱和失败则走饱和词干 + FG 位次 + 后缀路径，最后按需加 `cyclo`/`环` 前缀与 `_ylidene_form`。

### 不饱和段 `_chain_enyne`

按 `ene` / `yne` 位次记录收集段，段数决定形态：

- 两段共存时须共用同一形态（`fused` / `polyol` / `seg`），否则回落。
- `fused`：`ene_base` / `yne_suf` 提供融合基座，末段产出后缀。
- `polyol`：词干内嵌烯/炔段，FG 后缀后拼。
- `seg`：段式后缀（`ene_seg` / `yne_seg`），FG 位次经 `_fg_yl_tail` 拼接，C-1 自由价经 `_yl_loc_omitted` 省略位次（P-29.2）。

位次省略由 `_unsat_loc_omit` 统一判定，仅单段且形态为 `fused` / `polyol` 时（P-14.3.4.2(d)）生效。

### 多 FG 字段派生 `_generated_mult_fields`

多官能团时按倍数生成后缀：英文数量词以元音开头的场合去尾 `a`（P-14.3.2），`coda="ane"`、`need=mult`、`no_loc="none"`、`omit_rule=_NO_OMIT`。

- `acid`：`ene_base = ("ene{n}oic acid", "烯{n}酸")`，`yne_suf=None`。
- `ester`：`ene_base = ("ene{n}{en_suf}", "烯{n}{zh_suf}")`，`yne_suf = ("yne{n}{en_suf}", "炔{n}{zh_suf}")`，均由 `spec.en_suf` / `spec.zh_suf` 派生，硫代酯因此自动取 `ene{n}thioate` / `烯{n}硫` 一类词尾。

### 含氧酸词尾 `_OXO_TAIL`

以 `(oxo_kind, 中心酸式氢数)` 为键给出双语功能母体词尾。查表函数 `_oxoacid_tail(numbered)` 读 `parent.oxo_kind` 与 `parent.n_oh`，表外返回 None。

| oxo_kind | 氢数 | 英文词尾 | 中文词尾 |
| --- | --- | --- | --- |
| phosphate | 3 | phosphoric acid | 磷酸 |
| phosphate | 2 | dihydrogen phosphate | 磷酸二氢 |
| phosphate | 1 | hydrogen phosphate | 磷酸氢 |
| phosphate | 0 | phosphate | 磷酸 |
| phosphonate | 2 | phosphonic acid | 膦酸 |
| phosphonate | 1 | hydrogen phosphonate | 膦酸氢 |
| phosphonate | 0 | phosphonate | 膦酸 |
| sulfate | 2 | sulfuric acid | 硫酸 |
| sulfate | 1 | hydrogen sulfate | 硫酸氢 |
| sulfate | 0 | sulfate | 硫酸 |

`_oxoacid_tail` 经 `_Chain.plain_hook` 挂在 `phosphate` / `phosphonate` / `sulfate` 三条目上，使这三类中心母体的词尾由分子计数直接决定，不查碳数词表。

### `_KIND_TABLE` 条目

| kind | 词尾 (en / zh) | 关键字段 |
| --- | --- | --- |
| `alcohol` | `ol` / `醇` | `fg="alcohol"`、`need=1`、`omit_rule=_omit_term_locant`、`zh_full`、`unsat_polyol` |
| `ketone` | `one` / `酮` | `fg="ketone"`、`need=1`、`no_loc="none"`、`omit_rule`：`n<=2 and loc==1` |
| `alkane` | `ane` / `烷` | `coda=""`、`omit_rule`：`omit or n<=3`、`ene_base`/`yne_suf`、`ene_loc_omit`、`yne_loc_omit`、`wrap=_with_ez` |
| `acid` | `oic acid` / `酸` | `ene_base=("enoic acid","烯酸")`、`yne_suf=("ynoic acid","炔酸")`、`variant` C1/C2 保留名与草酸 |
| `sulfonic` | `sulfonic acid` / `磺酸` | `coda="ane"`、`fg="oxoacid"`、`need=1`、`omit_rule=_omit_term_locant`（P-65.3.1 取代式：链/环母体 + 磺酸后缀） |
| `sulfonate` | `sulfonate` / `磺酸` | `coda="ane"`、`fg="oxoacid"`、`need=1`、`omit_rule=_omit_term_locant`（磺酸酯：O 侧臂 + 磺酸酯后缀） |
| `sulfonamide` | `sulfonamide` / `磺酰胺` | `coda="ane"`、`fg="sulfonamide"`、`need=1`、`omit_rule=_omit_term_locant` |
| `sulfonyl_chloride` | `sulfonyl chloride` / `磺酰氯` | `coda="ane"`、`fg="oxoacid"`、`need=1`、`omit_rule=_omit_term_locant` |
| `ester` | `oate` / `酸` | `ene_base`、`yne_suf`、`variant` C1/C2 保留名与 oxalate |
| `phosphate` | `phosphate` / `磷酸` | `coda=""`、`fg="oxoacid"`、`plain_hook=_oxoacid_tail` |
| `phosphonate` | `phosphonic acid` / `膦酸` | `coda=""`、`fg="oxoacid"`、`plain_hook=_oxoacid_tail` |
| `sulfate` | `sulfate` / `硫酸` | `coda=""`、`fg="oxoacid"`、`plain_hook=_oxoacid_tail` |
| `acyl` | `oyl` / `酰基` | `ene_base`、`yne_suf`、`variant` C1/C2 保留名 |
| `thiol` | `thiol` / `硫醇` | `coda="ane"`、`fg="thiol"`、`need=1`、`ene_seg=("ene","烯")`、`yne_seg=("yne","炔")` |
| `amine` | `amine` / `胺` | `fg="amine"`、`need=1`、`zh_full`、`unsat_polyol` |
| `aldehyde` | `al` / `醛` | `ene_base=("enal","烯醛")`、`yne_suf`、`variant` C1/C2 |
| `nitrile` | `enitrile` / `腈` | `ene_base=("enenitrile","烯腈")`、`yne_suf`、`variant` C1/C2 |
| `amide` | `amide` / `酰胺` | `ene_base`、`yne_suf`、`variant` C1/C2 与草酰胺 |
| `acyl_halide` | `oyl {halide}` / `酰{卤}` | 由 `_ACYL_HALIDE_BY_HAL` 按卤素原子序数索引 |
| `radical` | `yl` / `基` | `coda="an"`、`no_loc="none"`、`yl_loc_omit`、`plain_fn=_radical_plain`、`wrap=_with_ez` |

`_ac_hal_chain(hal_z)` 按卤素构造酰卤 spec：词尾 `oyl {halide}` / `酰{卤}`，融合式烯/炔基座同步带卤素词；`variant` 给出 C1/C2 保留名（`formyl`/`acetyl`）与苯环外保留名 `benzoyl {halide}`。

### 苯环单取代保留名 `_BENZENE_RETAINED`

`_benzene_retained(en, zh, *, fg_drop=False, omit_all=False)` 构造 `variant = {1: {...}}`：`plain_fn` 恒返回给定双语对；`fg_drop` 清空 `fg`；`omit_all` 把 `omit_rule` 置为 `_omit_all`。

表内含 `alcohol`（phenol/苯酚，`fg_drop`）、`amine`（aniline/苯胺，`fg_drop`）、`radical`（phenyl/苯基，`omit_all`）、`acid`（benzoic acid/苯甲酸）、`sulfonic`（benzenesulfonic acid/苯磺酸，`omit_all`）、`sulfonate`（benzenesulfonate/苯磺酸，`omit_all`）、`sulfonamide`（benzenesulfonamide/苯磺酰胺，`omit_all`）、`sulfonyl_chloride`（benzenesulfonyl chloride/苯磺酰氯，`omit_all`）、`ester`（benzoate/苯甲酸）、`acyl`（benzoyl/苯甲酰基）、`aldehyde`（benzaldehyde/苯甲醛）、`nitrile`（benzonitrile/苯甲腈）、`amide`（benzamide/苯甲酰胺）。

## 母体拼接入口

### `join_kind_name` 三条路径

```
kind ∈ ESTER_O_SIDE_KINDS 且 kind ∈ OXO_CENTER_KINDS  →  join_oxoacid_name
kind ∈ ESTER_O_SIDE_KINDS                              →  join_ester_name
其余                                                   →  join_parent_name
```

- `ESTER_O_SIDE_KINDS = {ester, phosphate, phosphonate, sulfate, sulfonate}`：O 侧臂由本层整名消费，取代基前缀中已滤除 `o_side` 记录。
- `OXO_CENTER_KINDS = {phosphate, phosphonate, sulfate}`：中心原子自任母体，母体链位次无意义，臂直接挂在中心原子上。
- 普通路径对中英文分别调用 `join_parent_name`；中文侧经 `zh_1h_parent` 补 `1H-`（英文母体以 `1H-` 开头而中文尚未带时）。

`join_parent_name(prefix, parent)` 拼前缀与母体：前导为数字、`[` 或 `1H-` 时以连字符相接，并对母体名做 `_stereo_lead` 切分，使立体块保持在最前。

### `_names_for` 派发

`_names_for(kind, n, numbered)` 先取 `_KIND_TABLE[kind]` 作为 `entry`，再按分子上下文逐层改写：

1. **杂原子锚点**：`kind == "radical"` 且 `parent.radical_anchor_element` 存在时转 `_mononuclear_radical_names`。
2. **酰卤**：`kind == "acyl_halide"` 时用 `_ACYL_HALIDE_BY_HAL.get(parent.hal_z)` 覆盖默认的氯化物 spec。
3. **硫代酯条目改写**：`kind == "ester"` 且 `parent.thio_side` 时经 `dataclasses.replace` 把词尾换成 thioate 族并清空 `variant`（详见「酯与硫代酯」）。
4. **稠环词干直取**：`kind == "alkane"` 且 `parent.fused_tree` 存在、`scaffold_id != "benzene"` 时，返回 `_parent_stem_names`，由 `_ensure_fused_stem` 注入的稠合 base 名直接作为母体名。
5. **苯 base**：`scaffold_id == "benzene"` 且 `kind == "alkane"` 时返回 `("benzene", "苯")`。
6. **环式 FG 词干注入**：`parent.stem_en` / `parent.stem_zh` 同时存在时覆盖 `entry.stem` 并把 `coda` 置空（位次 omit 由 L4 决定），`aromatic` 置为 `scaffold_id == "benzene"`。
7. **单环饱和烃自由基**：`scaffold_id == "carbocycle"` 时置 `cyclic=True`、`ene_loc_omit=True`，`omit_rule` 对自由基为「`loc == 1`」，其余沿用 L4 的 `omit`。
8. **scaffold 专属保留名**：`variant` 先取苯环硫代酯专用的 `benzenecarbothioate`，再取 `_BENZENE_RETAINED[kind]`（`benzene` scaffold）或 `entry.variant.get(sid)`。

改写完成后交 `_chain_names(entry, n, numbered)` 生成双语名。无对应 kind 条目时回落 `_parent_stem_names`。

### 含氧酸中心母体 `join_oxoacid_name`

`join_oxoacid_name(pre, names, numbered)` 产出整名，结构为「取代前缀 + O-侧臂 + 功能母体词尾 + 金属盐」：

1. `tail_en = join_parent_name(pre[0], names[0])`；中文侧以 `zh_1h_parent` 处理后同法拼接。
2. 取 `_o_side_arms(numbered)`（`substituents` 中 `o_side` 为真的记录）。`oxo_kind == "phosphate"` 时走 `_fenced_arms_phosphate`，其余中心母体走 `_fenced_arms`。
3. `_join_o_side_arms(arms, group=True, arm_zh_fn=_oxoacid_arm_zh)` 拼接臂名（同基按 `_mult_rows` 归并）。
4. 末段分四种情形：无臂且有金属 → 酸式盐（英文金属前置、中文金属直接缀于酸式词后，如磷酸二氢钾）；无臂且 `n_om > 0` → 游离根（负电荷不标注，中文补「根」）；无臂 → 纯酸；有臂且有金属 → 酯盐（英文 `{metal} {arms} {tail}`，中文 `{tail}{arms}酯 {metal}盐`）；有臂无金属 → 英文 `{arms} {tail}`，中文 `{tail}{arms}酯`。

臂的中文词由 `_oxoacid_arm_zh` 渲染：简单基去「基」（`甲基` → `甲`），全数字词干补「烷基」，含连字符或括号的复合名原样保留。

## 酯与硫代酯

### `join_ester_name`

`join_ester_name(pre_en, pre_zh, names, numbered)` 以酸侧为母体、O 侧臂作前缀。中文 O 侧渲染由 `_zh_alkoxy_part` 决定：简单烃基省「基」（乙酸乙酯、十八酸苄酯），仅带位次的基保留「基」不加围栏（乙酸噻吩-3-基酯），其余按需加圆括号，已含围栏时升为方括号。异名臂依次平铺，同名臂用 `di`/`tri` 或复合前缀 `bis`/`tris`（`_join_o_side_arms`，P-16.3.2）。

### 硫代羧酸 S-酯

`_is_thio_side(numbered)` 读 `parent.thio_side`。为真时：

- `_names_for` 把 `ester` 条目改写为 `coda="ane"`、`en_suf="thioate"`、`zh_suf="硫"`、`ene_base=("enethioate","烯硫")`、`yne_suf=("ynethioate","炔硫")`、`variant=None`（P-65.6.3.3.7.1 硫代羧酸 S-酯无 C1/C2 保留名）。
- scaffold 为 `benzene` 时改取 `_benzene_retained("benzenecarbothioate", "苯硫代甲酸")`。
- 臂名经 `_fenced_arm_en` / `_fenced_arm_zh` 围栏：前导位次或自带括号的复合名须括起，内含圆括号时升方括号。
- 整名写作 `S-<臂> <母体>` / `S-<臂><母体>酯`，位次符号在最前；无臂时退化为 `S-<母体>`。

### O 侧臂拼接 `_join_o_side_arms`

`_o_side_arms(numbered)` 从 `substituents` 中挑出 `o_side` 为真的记录。`_join_o_side_arms(arms, *, group, arm_zh_fn)` 分三种情形：

- 单臂：英文名与经 `arm_zh_fn` 渲染的中文名直接返回。
- `group=True`（含氧酸中心母体）：按 `_mult_rows` 以 `alpha_order_key` 分组计数，同名臂前拼 `MULT_EN` / `MULT_ZH`，异名臂以空格分隔、中文直接相连。
- `group=False`（酯）：同名臂若任一记录带 `paren`，英文用 `BIS_EN` 复合前缀并括起（`bis(...)`），否则用 `MULT_EN` 平铺；异名臂依次平铺（`methyl ethyl oxalate`）。

### 单核自由基与桥前缀

`_mononuclear_radical_names(numbered)` 生成杂原子锚点自由基名，是 `radical` kind 在 `radical_anchor_element` 存在时的出口。流程按取代基数分派：

- 无取代基：取 `MONONUCLEAR_ZERO_YL[(stem_en, stem_zh)]`。
- P 酰基词干（`stem_en in PHOSPHORYL_STEMS`）：转 `_phosphoryl_sub_names`，按取代基拼接（P-67.1.4.1.1.5）。
- 单取代基：先试 sulfamoyl 融合（`sulfonyl` + `amino` → `…sulfamoyl`；`sulfonyl` + N-芳基按苯基并入，见下）；再试 `_fused_bridge_name` 双原子桥合一；再试 `_bridge_enclosed_names` 桥围栏；`azane` 中心另走 `AMIDO_RETAINED` 保留式与内层括起加 `amino` 两条支路；最后对 `f"{sub}-{stem}"` 调 `free_to_yl`，`oxidane` 中心再经 `_retained_alkoxy` 收拢为保留烷氧基。
- 多取代基：仅 `azane` 中心的双取代基成立（O/S 双烷基非标准自由基），按英文名排序后分「同基倍增」「N-芳基-N-某基取 anilino（P-62.2.1.1）」「首基平铺、其余各基加括号（P-62.2.2.1）」三路。

**苯磺酰/苯亚磺酰保留名**：`_mononuclear_en` 在基为 `phenyl` 且后缀为 `sulfonyl`/`sulfinyl` 时取 `benzenesulfonyl`（P-66.4.1）；`_mononuclear_zh` 对应 `base + 「基」` 走 `zh_bridge_root`。**N-芳基并入 sulfamoyl**：`stem_en == "sulfonyl"` 且取代基名以 `anilino` / `苯胺基` 结尾时，把 `anilino` 还原为 `phenyl` / `苯基` 后并入 `sulfamoyl` / `氨基磺酰基`，裸 `anilino` 时前端免括、带取代基时经 `_enclose` 括起（P-66.1.1.4.2）。

**双原子桥合一**：`_fused_bridge_name(stem_en, a)` 按 `BRIDGE_FUSION_YL` 把前端尾 `imino`/`sulfanyl` 与中心单核氢化物 `azane`/`sulfane` 合一为 `diazenyl` / `disulfanyl`（中文 二氮烯基 / 二硫代基，P-68.3.1.3/.4），该路径优先于桥围栏路径。

**桥围栏**：`_bridge_enclosed_names(a, stem_en, stem_zh)` 处理两类情形——P 酰基经 O/N/S 桥连母体时由 `MONONUCLEAR_BRIDGE` 给方括号；`sulfinyl`/`sulfonyl` 且取代基为复合名（`a["paren"]`）时经 `_enclose` 升级方括号，但直链 `-yl` 前端（`_SIMPLE_CHAIN_YL_RE`）与桥融合平铺而不拆。这些围栏在取代基词条定形时一次产出中英文双语，L5 的前缀渲染不重复拆分。命中时置 `numbered["bridge_self_enclosed"] = True`（`azane` 中心除外），标记该组已整体围栏。

**自由基的 free→yl 转换**：`free_to_yl(en, zh, attach_locant, *, paren=True)` 是 P-63.2.2 醇/胺的公共出口，`_fg_prefix` 要求中英文同时命中单核氢化物名。带环取代基的 `anilino` 与复合 `amino` 前缀（如 `methylamino`）须加括号，裸 `amino` / `anilino` 免括（P-29.3.6）。

## 环外主基命名

`_exo_ring_spec(spec, n, numbered)` 依据 `principal_expression_facts.relation == "exocyclic"` 改写 spec，后缀取自 `EXO_RING_SUF[kind]`（singular 单取代、plural 多取代基底、None 表示该主基无多取代系统名）。

- 环外硫代羧酸 S-酯：`ester` 且 `thio_side` 时后缀改用 `("carbothioate", "硫代甲酸")`。
- 环外酰卤：卤素词按 `parent.hal_z` 由 `HALIDE_EN` / `HALO_ZH` 动态拼接（如 `carbonyl chloride` / `甲酰氯`）。
- 多取代且 plural 为 None 时放弃；苯环单取代由 `variant` 保留名承担，不改写 spec。
- 稠环/杂环：词干已注入，位次恒显式。
- 纯碳环：环内有不饱和键时用段式后缀（`coda="ane"`、`cyclic=True`、`zh_loc_omit=False`）；饱和时取完整氢化物词干 `cyclo{alkane}` / `环{烷}` 并关闭 `cyclic`。`omit_rule` 为「`omit` 或本环无位次取代基」，故干净环省略主基位次（`cyclohexanecarboxylic acid`），环上另带被编号前缀时保留（P-66.6.1）。

## 前缀组装与围栏

`assembler_prefixes.py` 负责取代基前缀的分组、围栏与双语渲染。入口 `_prefix_for` 从 numbered 取上下文后委托 `_build_prefix`。

### `_build_prefix` 流程

1. `_fence_o_side_arms` 就地围栏 O 侧臂。
2. 滤除 `o_side` 记录（O 侧臂由 `join_kind_name` 消费）。
3. `_omit_sub_locants` 判定取代基位次可否省略。
4. `flat = kind in OXO_CENTER_KINDS`：为真时用 `oxo_arm_fence` 重算每个取代基的 `paren` 标记，使臂名围栏由环基判据决定，并去掉多余位次围栏。
5. `_mult_rows` 按英文名分组计数，排序键统一为 `alpha_order_key`（P-14.5 字母数字序）。
6. 计算 `bracket = omit and len(groups) >= 2 and _groups_simple(groups)`，据以选择 `sep`（`""` 或 `"-"`）并调用 `_collect_parts`。

### 位次省略 `_omit_sub_locants`

按母体情形分派：自由基母体仅单碳链省位次；`n_carbons <= 1` 省位次，但 N- 型与 C- 型取代基共存时 C 侧须带位次；纯烃环/苯环单取代位次隐含（环烯除外）；酰胺仅在取代基全为 N- 型时省；`acyl` / `ketone` / `acid` 恒不省；`kind in OXO_CENTER_KINDS` 恒返回 True（臂挂在中心原子上，母体链位次无意义）；含显式括号或自带位次的复合取代基不省；`n_carbons == 2` 的单取代仅在端碳无 H 的母体成立，`alcohol` / `amine` / `thiol` / `sulfonic` / `sulfonate` / `sulfonamide` / `sulfonyl_chloride` 除外。

### 通用臂围栏 `oxo_arm_fence`

`oxo_arm_fence(name, sub, mol)` 是含氧酸中心母体臂名的通用判据，命中任一条件即须围栏：

- 臂名已自带 `[...]` 完整围栏 → 不加。
- 前导 `(` 仅为立体描述符（`(2S)-…`），非完整围栏 → 加。
- 去掉立体描述符后仍有 ≥2 个位次段（`\d+(,\d+)*[a-z]*-`，`2,3-` 只算一段）的复合臂名 → 加（P-16.5.1.3.1）。
- 复合基（`sub["paren"]`）且臂根原子在环上 → 加；臂根取 claim 原子中与 `attach_idx` 成键者。

`_fenced_arm` 命中时调用 `_enclose`（含圆括号则升为方括号），`_fenced_arms` 对整组臂中英文同步改写。这两者服务于非磷酸的其余中心母体。

### 磷酸酯专用围栏 `_fenced_arm_phosphate`

- 臂名已自带括号或方括号 → 不加；简单基（`sub` 无 `paren`）→ 平铺；复合臂名再按 `oxo_arm_fence` 或位次段数 > 1 判定是否整体括起。
- `_fenced_arms_phosphate` 先检查 `_all_arms_stereo_lead`：臂数 ≥2 且全部臂名匹配 `_STEREO_LEAD_RE`（`(1Z)-`、`(2R,4R)-` 等前导立体描述符）时整体围栏（P-16.5.1.3.1 第二臂规则仅对多臂生效，单臂不受支配）。

### O 侧臂就地围栏 `_fence_o_side_arms`

`_O_SIDE_ARM_FENCE_KINDS = {"sulfonate", "phosphate"}`，仅这两类中心母体走此路径（assembler 侧对这两类不判臂围栏）。`_o_side_arm_fence(name, sub)` 判据：`sub` 非复合、臂名不带前导 `[`/`(` 且不含任何括号、且位次段数 ≥2 → 加围栏。处理时置 `arm_fenced` 标记，同一 numbered 被多次组装时不二次围栏；附着原子为硫（`GetAtomicNum() == 16`）的硫代酯 S-侧臂跳过，其围栏已在 assembler 侧完成。

### `flat` 参数

`flat` 贯穿 `_build_prefix` → `_collect_parts` → `_parts_for_stem` → `_prefix_one_en` / `_prefix_one_zh`：

- 英文 `_stem_needs_paren(stem, subs, omit, flat)`：词干以前导数字开头时返回 `not flat`，即中心母体的臂名前导位次不须围栏消歧。
- 中文 `_prefix_one_zh` 的围栏判据 `en[:1].isdigit() and not flat` 同步放宽。

### 括号式（P-16.5.1.3.1）

适用条件为 `omit and len(groups) >= 2 and _groups_simple(groups)`。`_groups_simple` 要求全部词干无前导位次且非 N- 类；`_collect_parts` 中首个引用的取代基从不加围栏（`bracket and not _LOCANT_RE.search(stem)` 时把该组 `paren` 全部清零），自带位次者除外；第 2 个起按 `MULT(词干)` 形式括起。

### 桥后缀拆分

`_split_bridge_suffix(stem)` 对 `BRIDGE_SPLIT_SUFFIX_EN`（`oxy` / `sulfanyl` / `amino` + 双原子桥 `diazenyl` / `disulfanyl`）逐个尝试，要求前端以 `_FRONT_TAILS`（`yl` / `sulfanyl` / `amino`）结尾且 `_front_needs_enclosure` 为真。中文侧 `_split_bridge_suffix_zh` 与英文侧同步判定，并在「基」被上游切掉时经 `_zh_front_ji` 在拆分点补回（酰基前端走融合式除外）。

`_bridge_body(base, suf, merge)` 把前端包上围栏、桥后缀留在括号外（P-63.2.2.1.2）；双原子桥因前端已围栏，整段再括一层；`merge` 用于前端已含方括号的场合，使桥后缀并入同一围栏。

相关正则：`_STEREO_LEAD_RE`、`_SIMPLE_CHAIN_YL_RE`、`_SUBST_CHAIN_YL_RE`、`_TERMINAL_CHAIN_YL_RE`、`_BENZYL_TAIL_RE`、`_LOCANT_RE`、`_ACYL_FRONT_RE`、`_LOCANT_SUBST_TAIL_RE`、`_MULT_WRAP_RE`。

### 倍数与 N- 前缀

`_mult_of(lang, stem, subs, n)`：词干含 `carboxy`/`羧` 哨兵或该组任一取代基带 `paren` 时用 `bis`/`tris`（`_COMPLEX_MULT_LANG`），否则用 `di`/`tri`。`_MULT_WRAP_RE`（`-\d+-yl$`、`oyloxy$`）命中时强制加括号并把倍数词换为 `bis`。

`_parts_for_stem` 在整组取代基均为 N- 型（`N_PREFIX_KINDS = {n_alkyl, n_block}`）时走 N- 计数通道：位次以 `N` 标注，`_n_prime_map` 按引用序给不同 N 原子分配撇号个数（`N,N-` 同氮、`N,N'-` 跨氮），由 `_n_prefix` 拼接。

## 指示氢与 hydro 前缀

`join_hydro_prefix(names, numbered)` 依次处理三类前缀。

1. **L4 指示氢位次**：`_indicated_h_prefix` 把 `parent.indicated_h_locants` 去重后拼成 `nH,` 串。当该串非空且（存在 `hydro_prefix`、或 `parent.indicated_h_forced`、或 `_nh_only` 且 `_ring_n_free`）时，用 `_with_indicated_h` 把它落到词干前：先由 `_split_stem_h_prefix` 经 `_STEM_H_PREFIX_RE`（`^(?:\d+[a-z]?H[-,])+`）剥掉词干自带前缀，再前置 L4 位次，避免 `1H-7H-` 双写。
2. **未取代 NH 环补前缀**：`_nh_only` 经 `_h_prefix_atoms` 把前缀位次映射回母体链原子，要求全部为氮；`_ring_n_free` 要求环内每个氮都不带环外邻居。二者同时成立即补写，N-取代基环（如 N-糖基尿嘧啶）省略 NH 前缀。
3. **词干自带前缀失效**：L4 无位次时，`_stem_prefix_stale` 检查词干自带前缀的位次原子是否已成为羰基碳（无 H、且有指向环外非碳原子的非单键），命中即丢弃该前缀（如茚-1,3-二酮的 C1）。

最后把 `parent.hydro_prefix` 经 `join_parent_name` 拼到最前。`_ensure_fused_stem` 注入稠环词干时同样前置一次指示氢位次。

## 环内阳离子与阴离子后缀

`join_ring_cation_suffix(numbered, names)` 处理环内 N⁺ / O⁺（P-62.4.1）：

- 净电荷为负的分子直接跳过（盐与两性离子都可能含环阳离子，净负时不处理）。
- 名称已含 `ium` 时不重复；在 `chain` 内找出带 +1 形式电荷的 N / O，取其位次标签。
- 词干以 `ene` 结尾时按色烯型氧鎓处理（`chromene` → `chromenylium`），以 `e` 结尾时缀 `-{位次}-ium`，否则直接在词干后接 `-{位次}-ium`。
- 位次插入点分三种：完整母体名出现词干时整词干替换；出现 `base + "e"` 时在该 token 后插入；否则在 `base` 后插入。

`join_anion_names` 在 `parent.anion` 为真时把羧酸后缀转为阴离子式：`acid_to_anion_en` 把 `-oic acid` → `-oate`、`-ic acid` → `-ate`，中文补「根」（已带「根」则不重复）。`stems.py` 另提供 `join_metal_salt_names`，按 `salt.metal` 加金属数量前缀（英文前置、中文把「酸根」替换为金属后缀），或在 `salt.acid_salt` 存在时改缀酸式盐词。

## 稠合名组装 `fused_namer.py`

`fused_parent_names(mol, node)` 由 L2 打包的 `FusedNode` 树生成未注册稠环的 base 名，编号仍走 L4。

- `_component_prefix(node)`：附加组分前缀取 `node.fused_prefix`，否则「去尾 `e` 加 `o`」/「词干 + 并」（P-25.3.2.2.2）。
- `_component_numbering`：委托 L4 的 `fused_component_numbering` 按共享原子取向给子环编号。
- `_fusion_letter(parent_chain, shared)`：共享边须落在母体外周边界上，返回侧字母 `a`/`b`/…。
- `_fusion_numbers(child_chain, child_labels, parent_chain, shared)`：附加组分共享原子的位次，顺序沿母体低位次端到高位次端。
- `_inner_atoms(node, rings)`：出现在 ≥3 个环的 perifused 中心原子，从母体外周链中剔除。
- `_fused_one`：一级单环烃附加组分省略数字位次（`fused_omit_numbers`，P-25.3.8.1）时输出 `前缀[字母]`；否则输出 `前缀[数字-字母]`（P-25.3.2）。
- `_collect_attached`：递归收集附加组分，二级组分前缀排在附着的一级组分前。
- 最终把附加前缀串到根词干前，并查 `RETAINED_FUSION_ALIASES` 换取保留名；单节点（无附加组分）交由 `_parent_stem_names` 处理。

## 词干表 `stems.py`

`stems.py` 提供碳数词干与盐后缀的底层词表：`_ALKANE_EN_BASE` / `_ALKANE_ZH_BASE` 覆盖 C1–C10 保留形，`_en_stem(n)` 去尾 `ane`（≥11 去尾 `a`，如 undec/icos），`alkane_en` / `alkane_zh` 产出全名，`zh_num(n)` 对 1–10 用天干（甲…癸）、11 以上复用 `zh_numeral`，`zh_stem(zh_full)` 剥掉末端官能团后缀（十一烷 → 十一）。`_ZH_SUFFIXES` 列出参与剥离的中文后缀。盐侧 `_metal_prefix(metal, n, mult)` 对 `n == 1` 不加数量词、`n > 1` 取 `MULT_EN` / `MULT_ZH`。

## 立体化学 `stereo.py`

本模块产出名称最前的 `(…)-` 立体块，对 E/Z 与 R/S 共用一套解析与排序设施。

- **切分与解析**：`_split_stereo_lead` 切出前导 `(…)-` 块；`_parse_stereo` / `_parse_token` 解析为 `(位次|None, 字母)` 部件列表；`_format_stereo` 以 `_part_key` 排序（无位次者在前，余按 `locant_key`）后拼回前缀。
- **E/Z**：`_bond_stereo` 读 RDKit 双键 `GetStereo()`；`_ez_prefix` 对单双键取名并补母体链最小位次（`_bond_min_loc`）；`_ez_multi_prefix` 对 `double_bonds` 收集所有已定义立体的键；`ez_for_parent` 按有无 `double_bonds` 分派。三者经 `_Chain.ez_ene` / `ez_ene_multi` 挂到链引擎。
- **外挂双键**：`_exo_ez_parts` 只收集一端在母体内的立体双键，位次取母体侧；`join_ez_prefix` 是 E/Z 对外入口，先剔除母体名已带位次再补写，`ester` 时经 `_ester_en_rs` 把部件插到烷基词与酰基词之间。
- **R/S**：`_cip_on_chain` 调 `assign_cip` 后逐链原子读 `_CIPCode`；`_chain_locant` 把链序号映射为位次标签；`_collapsed_parent` 在母体 C1 由更大或环状分子折叠而来时跳过；`_rs_parts` 汇总部件。`join_rs_prefix` 为对外入口，`_with_rs` 把 R/S 部件并入名称已有立体前缀（保留非 R/S 部件并按位次重排）。

## 文件与接口

| 文件 | 职责 |
| --- | --- |
| `__init__.py` | 包声明，指向入口 `assemble` |
| `assembler.py` | 组装入口、母体拼接分派、酯/含氧酸整名、自由基名、指示氢 |
| `assembler_prefixes.py` | 取代基前缀分组、围栏判据与双语渲染 |
| `chain_engine.py` | `_Chain` 规格、`_KIND_TABLE`、词尾表、链式词干引擎 |
| `stems.py` | 碳数词干表、酸转阴离子、金属盐后缀 |
| `stereo.py` | E/Z 与 CIP R/S 立体前缀 |
| `fused_namer.py` | 未注册稠环的稠合 base 名组装 |

### 对外接口

```python
assemble(numbered: dict, *, time_ms: float = 0.0, source: str = "iupac") -> NameResult
```

上游只需提供 numbered dict；失败时返回 `success=False` 的 `NameResult`，`meta` 中给出 `reason` / `n_carbons` / `kind`。

### 跨层依赖

- 自 L1：`mol`（RDKit）、单核母体氢化物词表、官能团类别。
- 自 L2：`kind`、`n_carbons`、`scaffold_id`、`oxo_kind`、`n_oh`、`salt_meta`、`fused_tree`。
- 自 L4：`chain`、`numbering_scaffold.labels`、`fg_locants`、`indicated_h_locants`、`hydro_prefix`、`double_bond(s)`。
- 本层回注：`fused_namer` 依赖 L4 的 `fused_component_numbering`；`stereo` 依赖 L4 的 `assign_cip` 与 `locant_key` / `locant_str_sort`。
