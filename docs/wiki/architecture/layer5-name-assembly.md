# Layer5: 名称组装 (Name Assembly)

> 源文件 10 个（约 2890 行） · 对外接口 `assemble`

## 概述

Layer5 是命名流水线的末层：把 Layer4 输出的编号结果（numbered dict）渲染为双语 IUPAC 名称。本层不做化学推断，只做名称定形——母体词干注入、官能团词尾、取代基前缀、立体描述符、盐后缀都在此拼接。产物为 `NameResult(en, zh, success, source, time_ms, meta)`：失败经 `_fail` 构造（`success=False`，`meta` 携带 `reason`、`n_carbons`、`kind`），成功经 `_ok` 构造，`_unsupported(n, kind)` 是上述三项诊断的 `_fail` 快捷封装。

### 输入 numbered dict

| 键 | 本层用途 |
| --- | --- |
| `parent.kind` / `parent.n_carbons` | 选 `_KIND_TABLE` 条目、生成碳数词干 |
| `parent.stem_en` / `parent.stem_zh`（含 `stem_bare_*` 裸词干） | 覆盖碳数词干；裸词干供链引擎自拼 `ane` / 「烷」与不饱和段 |
| `parent.scaffold_id` | 骨架 ID（`benzene` / `carbocycle` / `bridged` / `mono_spiro` / `fused_bridged_spiro` / 稠环）：分派保留名与环前缀 |
| `parent.spiro_node` / `parent.fbs_node` / `parent.bridged_node` / `parent.fused_tree` | L2 打包、L4 选中的环系节点：螺环 / 桥环 / 稠合词干注入 |
| `parent.chain` / `parent.mol` | 母体原子序列与 RDKit 分子：指示氢位次映射、臂根成环判定、手性中心、环阳离子、C1 杂原子环境 |
| `parent.numbering_scaffold.labels` | 原子 → 位次标签（组分式螺环带撇号位次 `1'`）：位次串与环阳离子位次 |
| `parent.locant_kind` / `parent.subs_consumed` / `parent.stem_generated` / `parent.bridge_self_enclosed` | 位次省略 kind、取代基已并入保留名、生成式词干、桥前缀已自含围栏：`_omit_sub_locants` 判据、前缀关闭、自由基取裸词干、臂前缀不二次围栏 |
| `parent.oxo_kind` / `parent.n_oh` / `parent.thio_side` / `parent.hal_z` | 含氧酸中心元素与酸式氢数、硫代羧酸 S-酯标记、卤素原子序数：`_OXO_TAIL` 查表、thioate 词尾与 `S-` 位次符号、酰卤词尾 |
| `parent.radical_anchor_element` / `parent.free_valence_order` | 杂原子锚点自由基标记 / 自由价键级：单核自由基路径、`-ylidene` / `-ylidyne` |
| `parent.indicated_h_locants` / `parent.indicated_h_forced` / `parent.hydro_fallback` / `parent.hydro_prefix` | 指示氢位次与加氢前缀：拼 `nH-` 前缀（`hydro_fallback` 为真时强制写出指示氢） |
| `parent.anion` / `parent.salt_meta` / `parent.n_om` | 酸根标记、金属盐记录、游离根计数：阴离子式与盐式 |
| `substituents` | 取代基记录（`o_side` / `paren` / `attach_idx` / `atoms` / `kind` / `en` / `zh` / `arm_fenced` / `bare`）：前缀分组、围栏与位次省略 |
| `fg_locants` | 主官能团位次记录（`kind` / `locants` / `omit`）：词尾位次串 |
| `ene_locants` / `yne_locants` / `omit_ene_locant` / `omit_yne_locant` / `double_bond` / `double_bonds` | 不饱和键位次与省略标记、母体内立体双键：不饱和段与 E/Z 前缀 |

### `assemble` 步骤

`assemble(numbered, *, time_ms=0.0, source="iupac")` 按固定顺序执行：

1. `_ensure_parent_stem` — 母体词干注入（螺环 / 桥环 / 稠环 / 大环杂单环生成式），任一环节失败即返回 unsupported。
2. `_names_for(kind, n, numbered)` — 取母体双语名。
3. `join_hydro_prefix` — 拼接 hydro 前缀与指示氢前缀。
4. `join_ring_cation_suffix` — 环内 N⁺ / O⁺ 缀 `-{位次}-ium`（P-62.4.1）。
5. `join_kind_name(kind, _prefix_for(numbered, kind, n), names, numbered)` — 取代基前缀与母体拼接。
6. `join_anion_names` → `join_ez_prefix` → `join_rs_prefix` — 依次补羧酸阴离子式（`-ate` / 酸根）、母体外挂双键的 E/Z 部件与母体手性中心的 CIP R/S 部件（E/Z 先补，R/S 再归并排序）。

### 职责边界

本层不选母体、不编号、不推断官能团：`kind`、`scaffold_id`、`chain`、`fg_locants` 与全部位次均由上游给定；也不 import L2，环系节点只按鸭子类型读字段。失败一律以 `success=False` 返回，不抛异常。

## 核心逻辑

### 母体词干注入分派 `_ensure_parent_stem`

`_ensure_parent_stem(numbered)` 是 `assemble` 入口的词干注入分派，各子注入器共用同一幂等守卫：`parent.stem_en` 与 `parent.stem_zh` 已在即直接放行。分派顺序为 `fbs_node` 非空走 `_ensure_fbs_stem`、`spiro_node` 非空走 `_ensure_spiro_stem`、`bridged_node` 非空走 `_ensure_bridged_stem`，其余先 `_ensure_fused_stem` 再 `_ensure_generated_stem`。

- `_ensure_ring_stem(numbered, node_key, namer)` 是 P-23 / P-24 的公共注入器，由 `_ensure_spiro_stem`（`spiro_node` + `spiro_parent_names`）与 `_ensure_bridged_stem`（`bridged_node` + `bridged_parent_names`）调用。注入器返回 `((完整英, 完整中), (裸词干英, 裸词干中))`：裸词干写入 `stem_bare_en` / `stem_bare_zh`；环内有 `ene_locants` / `yne_locants` 时 `stem_en` / `stem_zh` 取裸词干（由链引擎自拼 `ane` / 「烷」与不饱和段），否则取完整名。分子、节点或链缺失、注入器返回 None 均判失败。
- `_ensure_fbs_stem` 由 `spiro_namer.fbs_parent_name` 取单一双名对写入 `stem_en` / `stem_zh`，不设裸词干形态（组分名自带不饱和，母体层不出 `ene` / `yne` 位次）。
- `_ensure_fused_stem` 在 `fused_tree` 缺失时放行，否则取 `fused_namer.fused_parent_names` 的稠合 base 名、前置一次 `_indicated_h_prefix`（L4 已按整体编号定好指示氢位次，P-58.2.1）后写入 `stem_en` / `stem_zh`。
- `_ensure_generated_stem` 处理 >10 元无保留名的环（`scaffold_id == "carbocycle"` 且 `chain` 长度 > 10）：`prefix_from_chain` 取 'a' 前缀，`stems.stem_forms` 拼 `{前缀}cyclo{烷词干}` / `{前缀}环{烷}` 两形态并置 `parent.stem_generated`；无杂原子（空串）放行，词表外元素（None）判失败，≤10 元环归 Hantzsch-Widman 保留名直接放行。

### 母体名派发 `_names_for`

`_names_for(kind, n, numbered)` 先试 `_c1_retained`（单碳母体保留名，命中即直接返回），再做一次 radical 锚点短路，随后取 `_KIND_TABLE` 条目并按分子上下文逐层改写：

1. **杂原子锚点**：`kind == "radical"` 且 `parent.radical_anchor_element` 存在时转 `_mononuclear_radical_names`。
2. **酰卤**：`kind == "acyl_halide"` 时用 `_ACYL_HALIDE_BY_HAL.get(parent.hal_z)` 覆盖默认的氯化物 spec。
3. **裸词干环系直取**：`parent.stem_bare_en` 存在且 `kind` 属 `("bridged", "alkane", *_SPIRO_KINDS)` 时，取 `_KIND_TABLE["alkane"]` 并把 `stem` 换成裸词干、`coda` 置空后交 `_chain_names`，使环内 `ene` / `yne` 位次照常渲染（`_SPIRO_KINDS = ("mono_spiro", "fused_bridged_spiro")` 是 P-24 螺环 scaffold 直取的 kind，不在 `_KIND_TABLE` 内）。
4. **硫代酯条目改写**：`kind == "ester"` 且 `parent.thio_side` 时经 `dataclasses.replace` 把词尾换成 thioate 族并清空 `variant`。
5. **稠环词干直取**：`kind == "alkane"` 且 `parent.fused_tree` 存在、`scaffold_id != "benzene"` 时返回 `_parent_stem_names`；`scaffold_id == "benzene"` 且 `kind == "alkane"` 时直接返回 `("benzene", "苯")`。
6. **环式 FG 词干注入**：`parent.stem_en` / `parent.stem_zh` 同时存在时覆盖 `entry.stem`、把 `coda` 置空、`omit_rule` 收敛为「沿用 L4 的 `omit`」、`aromatic` 置为 `scaffold_id == "benzene"`；词干再经 `_dihydro_stem_fix` 校正。生成式词干上的自由基（`stem_generated`）改取裸词干（P-29.2），使名形为 `…tetrazacyclododec-1-yl`。
7. **单环碳环**：`scaffold_id == "carbocycle"` 时置 `cyclic=True`、`ene_loc_omit=True`，`omit_rule` 对自由基为「`loc == 1`」，其余沿用 L4 的 `omit`。
8. **scaffold 专属保留名**：`variant` 先取苯环硫代酯专用的 `benzenecarbothioate`，再取 `_BENZENE_RETAINED[kind]`（`benzene` scaffold）或 `entry.variant.get(sid)`。

改写完成后交 `_chain_names(entry, n, numbered)` 生成双语名，结果再经 `_c1_amino` 做氨基甲酸类词尾替换；`kind` 无表条目（含螺环 scaffold kind）时回落 `_parent_stem_names`（仅返回 `parent.stem_en` / `parent.stem_zh`，不做任何后缀拼接）。

### 母体拼接与 C1 保留名

`join_kind_name` 三条路径：`kind` 同时属 `ESTER_O_SIDE_KINDS` 与 `OXO_CENTER_KINDS` 时走 `join_oxoacid_name`；仅属 `ESTER_O_SIDE_KINDS` 时走 `join_ester_name`；其余走 `join_parent_name`。

- `ESTER_O_SIDE_KINDS = {ester, phosphate, phosphonate, sulfate, sulfonate}`：O 侧臂由本层整名消费（取代基前缀中已滤除 `o_side` 记录）；`OXO_CENTER_KINDS = {phosphate, phosphonate, sulfate, boronic}` 的中心元素自任母体、母体链位次无意义、臂直接挂在中心原子上（`boronic` 不在 `ESTER_O_SIDE_KINDS` 内，其碳基作普通取代基前缀）。
- 普通路径对中英文分别调用 `join_parent_name`（经 `_stereo_lead` 切出立体块后，余下词干前导为数字、`[` 或 `1H-` 时以连字符相接，否则直连，最后把立体块复原到最前）；中文侧先经 `zh_1h_parent` 补 `1H-`（英文母体以 `1H-` 开头而中文尚未带时）。

**`_c1_retained`** 识别单碳母体官能团碳上的两个杂原子，命中即返回保留名双语对。`_c1_two_hetero(parent)` 按 `parent.chain` 长度为 1 定位该碳，统计单键 N、单键 O 与双键杂原子（氢与锚定哑原子不计，非单 / 双键键型即不适用）：`radical` + 一个单键 N + `=O` 且无双键 O → 伯酰胺取 `carbamoyl` / `氨基甲酰基`（N 带取代基时经 `carbamoyl_prefix_name` 产出 `<N-取代基>carbamoyl`，取不到名即回落系统名），命中后置 `parent.subs_consumed` 关闭取代基前缀；`amide` + 两个单键 N + `=O` → `urea` / `脲`，`amine` + 两个单键 N + `=S` → 先撤下 `=S` 侧取代基块再返回 `thiourea` / `硫脲`，`=N` 同法得 `guanidine` / `胍`；`ester` + 两个单键 O + `=O` 且无单键 N → `carbonate` / `碳酸`。

脲 / 硫脲 / 胍经 `_urea_subs` 把 N-取代基改走数字位次通道：全部取代基 `kind` 置 `"side"`，按所连 N 分组后以组内首个取代基名的 `alpha_order_key` 排序，依次取 `("1","3")`（胍为 `("2","3")` 且强制编号）写入 `s["locant"]`；全分子只有一个取代基且非强制时不写位次（`(4-甲基苯基)脲`），成功时另置 `parent.locant_kind`。

**`_c1_amino`** 判官能团碳是否直连 N（`chain` 长度 1、该碳为碳、邻居含氮），命中时替换 `_chain_names` 结果的词尾：`ester` 的 `formate` / `甲酸` → `carbamate` / `氨基甲酸`（P-66.3.2）；`acid` 的 `formic acid` / `甲酸` → `carbamic acid` / `氨基甲酸`，并置 `parent.locant_kind = "carbamic_acid"`、把全部取代基 `kind` 改为 `"side"`（N-取代基因而不带位次，`dimethylcarbamic acid`）；`acyl_halide` 的 `formyl` / `甲酰` → `carbamoyl` / `氨基甲酰`（P-66.1.1.4.1）。

### 链引擎与 `_KIND_TABLE`

链引擎以数据驱动方式生成链状母体名：一张 `_KIND_TABLE` 描述每个 kind 的词尾、位次规则与特例，由 `_chain_names` 统一消费。

**`_Chain` 规格**（`frozen dataclass`）：词尾 `kind` / `en_suf` / `zh_suf` / `coda`（饱和词干后接段，默认 `"an"`，部分 kind 用 `""` 或 `"ane"`）；位次 `no_loc`（`"plain"` 无位次回落普通名 / `"none"` 返回 None）、`omit_rule`（`(n, loc, omit) -> bool`，True 表示省略位次）、`fg`、`need`、`ene_loc_omit`、`yne_loc_omit`、`yl_loc_omit`、`zh_loc_omit`；不饱和段 `ene_seg` / `yne_seg`（段式单烯 / 炔段，默认 None 表示由后缀首字母推导）、`ene_base` / `yne_suf`（融合式基座）、`unsat_polyol`、`zh_full`、`ez_ene` / `ez_ene_multi`；覆盖与派生 `plain_maps` / `plain_fn` / `plain_hook`（俗名表 / 按碳数派生命名 / 按分子计数直接给词尾）、`variant`（`{scaffold 或 multiplicity: 字段覆盖 dict}`）、`stem`、`aromatic`、`cyclic`、`cyclic_unsat`、`wrap`。位次规则常量 `_NO_OMIT`（恒不省）、`_omit_all`（恒省）、`_omit_term_locant`（`omit` 或无位次，或 `n <= 2` 且位次为 1 时省）由 `omit_rule` 引用。

**主流程 `_chain_names`**：先调 `plain_hook`（返回非 None 即直接作为结果，无碳链母体如磷酸等 P / B 中心）；芳香环醇统一改 `zh_suf="酚"`；`_exo_ring_spec` 改写环外主基 spec；按 `_parent_multiplicity`（优先 `principal_expression_facts.multiplicity`，回落 `principal_group_count`）走多 FG（`_generated_mult_fields` 生成数量后缀字段并叠加 `variant[mult]`）或单 FG（关闭 `unsat_polyol` / `zh_full` 后叠加 `variant[1]`）分支；`parent.free_valence_order` 仅在 `kind == "radical"` 时生效，>1 时经 `_free_valence_form` 出 `-ylidene` / 亚基（双键）或 `-ylidyne` / 次基（三键，P-31.2.3）；再交 `_chain_enyne` 尝试不饱和段渲染（`unsat_polyol` 时补 FG 位次串，位数须等于 `need`），失败则走饱和词干 + FG 位次 + 后缀路径，最后按需加 `cyclo` / 环 前缀与 `wrap`。

**不饱和段 `_chain_enyne`**：按 `ene` / `yne` 位次记录收集段，两段共存时须共用同一形态（`fused` / `polyol` / `seg`），否则回落。`fused` 由 `ene_base` / `yne_suf` 提供融合基座、末段产出后缀；`polyol` 把烯 / 炔段嵌入词干、FG 后缀后拼；`seg` 走段式后缀，FG 位次经 `_fg_yl_tail` 拼接，C-1 自由价经 `_yl_loc_omitted` 省略位次（P-29.2）。注入的完整母体名（`spec.stem`）以 `ane` 结尾时先去尾 `ane` 再接不饱和段（`azacyclododecane` → `azacyclododec-9-en`），两段以上插连接元音 `a`。段式段的词形由 `_unsat_seg(spec, b)` 定：`ene_seg` / `yne_seg` 显式给出时直接采用，否则看 `spec.en_suf` 首字母——非元音（`aeiouy` 之外）保留末端 `e`（`carboxylic acid`、`thiol` 一类，P-16.7.1(a)、P-59.1.9），元音开头用 `en` / `yn`。段位次省略由 `_unsat_loc_omit` 判定，仅单段且形态为 `fused` / `polyol` 时生效（P-14.3.4.2(d)）：烯段看 `spec.ene_loc_omit` 或形态为 `polyol`，炔段看 `spec.yne_loc_omit` 或形态为 `polyol`，两者都另要求 L4 的 `omit_ene_locant` / `omit_yne_locant`。

**多 FG 字段派生 `_generated_mult_fields`**：多官能团时按倍数生成后缀（`coda="ane"`、`need=mult`、`no_loc="none"`、`omit_rule=_NO_OMIT`、清空 `plain_maps`，英文数量词尾 `a` 在后缀元音开头时省略，P-14.3.2）。`acid` 在环外羧酸（`fg` 存在且非 `cyclic`）时关闭 `ene_base` / `yne_suf` 并开 `unsat_polyol`，不饱和段独立保留 `e`（`prop-1-ene-1,2,3-tricarboxylic acid`），其余场合取 `ene_base = ("ene{n}oic acid", "烯{n}酸")`、`yne_suf = None`；`ester` 取 `ene_base = ("ene{n}{en_suf}", "烯{n}{zh_suf}")`、`yne_suf = ("yne{n}{en_suf}", "炔{n}{zh_suf}")`，由 `spec.en_suf` / `spec.zh_suf` 派生，硫代酯因此自动取 `ene{n}thioate` / `烯{n}硫` 一类词尾。

**含氧酸词尾 `_OXO_TAIL`**：以 `(oxo_kind, 中心酸式氢数)` 为键给出双语功能母体词尾，由 `_oxoacid_tail(numbered)` 读 `parent.oxo_kind` 与 `parent.n_oh` 查表，表外返回 None。

| oxo_kind | 氢数 | 英文词尾 | 中文词尾 |
| --- | --- | --- | --- |
| phosphate | 3 / 2 / 1 / 0 | phosphoric acid / dihydrogen phosphate / hydrogen phosphate / phosphate | 磷酸 / 磷酸二氢 / 磷酸氢 / 磷酸 |
| phosphonate | 2 / 1 / 0 | phosphonic acid / hydrogen phosphonate / phosphonate | 膦酸 / 膦酸氢 / 膦酸 |
| sulfate | 2 / 1 / 0 | sulfuric acid / hydrogen sulfate / sulfate | 硫酸 / 硫酸氢 / 硫酸 |
| boronic | 2 / 1 / 0 | boronic acid / hydrogen boronate / boronate | 硼酸 / 硼酸氢 / 硼酸 |

`_oxoacid_tail` 经 `_Chain.plain_hook` 挂在 `phosphate` / `phosphonate` / `sulfate` / `boronic` 四条目上，故这几类中心母体的词尾由分子计数直接决定，不查碳数词表。

**`_KIND_TABLE` 条目**：

| kind | 词尾 (en / zh) | 关键字段 |
| --- | --- | --- |
| `alcohol` / `amine` / `thiol` | `ol`/`醇`、`amine`/`胺`、`thiol`/`硫醇` | `need=1`、`zh_full`、`unsat_polyol`；`thiol` 用 `coda="ane"` |
| `ketone` | `one` / `酮` | `no_loc="none"`、`unsat_polyol`、`omit_rule`：`n<=2 and loc==1` |
| `alkane` | `ane` / `烷` | `coda=""`、`omit_rule`：`omit or n<=3`、`ene_base`/`yne_suf`、`ene_loc_omit`、`yne_loc_omit`、`wrap=_with_ez` |
| `acid` / `ester` | `oic acid`/`酸`、`oate`/`酸` | `ene_base`/`yne_suf`、C1/C2 保留名（草酸 / oxalate） |
| `sulfonic` / `sulfonate` / `sulfonamide` / `sulfonyl_chloride` | `sulfonic acid`/`磺酸`、`sulfonate`/`磺酸`、`sulfonamide`/`磺酰胺`、`sulfonyl chloride`/`磺酰氯` | `coda="ane"`；`sulfonamide` 的 `fg="sulfonamide"`，余为 `fg="oxoacid"` |
| `phosphate` / `phosphonate` / `sulfate` / `boronic` | `phosphate`/`磷酸`、`phosphonic acid`/`膦酸`、`sulfate`/`硫酸`、`boronic acid`/`硼酸` | `coda=""`、`plain_hook=_oxoacid_tail` |
| `acyl` / `aldehyde` / `nitrile` / `amide` | `oyl`/`酰基`、`al`/`醛`、`enitrile`/`腈`、`amide`/`酰胺` | `ene_base`/`yne_suf`、C1/C2 保留名（`amide` 另带草酰胺） |
| `acyl_halide` | `oyl {halide}` / `酰{卤}` | 由 `_ACYL_HALIDE_BY_HAL` 索引 |
| `radical` | `yl` / `基` | `coda="an"`、`no_loc="none"`、`yl_loc_omit`、`plain_fn=_radical_plain`、`wrap=_with_ez` |

`_ac_hal_chain(hal_z)` 按卤素构造酰卤 spec：词尾 `oyl {halide}` / `酰{卤}`，融合式烯 / 炔基座同步带卤素词；`variant` 给出 C1/C2 保留名（`formyl` / `acetyl`）、C2 的 `oxalyl di{halide}` / `草酰二{卤}`（P-65.1.1）与苯环外保留名 `benzoyl {halide}`。

**苯环单取代保留名 `_BENZENE_RETAINED`**：`_benzene_retained(en, zh, *, fg_drop=False, omit_all=False)` 构造 `variant = {1: {...}}`，`plain_fn` 恒返回给定双语对；`fg_drop` 清空 `fg`，`omit_all` 把 `omit_rule` 置为 `_omit_all`。表内含 `alcohol`（phenol / 苯酚）、`amine`（aniline / 苯胺）两条 `fg_drop`；`radical`（phenyl / 苯基）与磺酸族四条（benzenesulfonic acid / benzenesulfonate / benzenesulfonamide / benzenesulfonyl chloride）共五条 `omit_all`；其余 `acid` / `ester` / `acyl` / `aldehyde` / `nitrile` / `amide` 取苯甲酸族系统保留名（benzoic acid / benzoate / benzoyl / benzaldehyde / benzonitrile / benzamide）。

**指示氢与 hydro 前缀**：`join_hydro_prefix(names, numbered)` 依次处理三类前缀。

1. **L4 指示氢位次**：`_indicated_h_prefix` 把 `parent.indicated_h_locants` 去重后拼成 `nH,` 串。该串非空且（存在 `hydro_prefix`、或 `parent.indicated_h_forced`、或 `parent.hydro_fallback`、或 `_nh_only` 且 `_ring_n_free`）时用 `_with_indicated_h` 落到词干前：先由 `_split_stem_h_prefix` 经 `_STEM_H_PREFIX_RE`（`^(?:\d+[a-z]?H[-,])+`）剥掉词干自带前缀，再前置 L4 位次，避免 `1H-7H-` 双写。
2. **未取代 NH 环补前缀**：`_nh_only` 经 `_h_prefix_atoms` 把前缀位次映射回母体链原子，要求全部为氮；`_ring_n_free` 要求环内每个氮都不带环外邻居；二者同时成立即补写，N-取代基环（如 N-糖基尿嘧啶）省略 NH 前缀。
3. **词干自带前缀失效**：L4 无位次时，`_stem_prefix_stale` 检查词干自带前缀的位次原子是否已成为羰基碳（无 H、且有指向环外非碳原子的非单键），命中即丢弃该前缀（如茚-1,3-二酮的 C1）。

最后把 `parent.hydro_prefix` 经 `join_parent_name` 拼到最前。

**二氢词干校正 `_dihydro_stem_fix`**：词干自带的 `N,M-二氢`（`_DIHYDRO_STEM_RE` / `_DIHYDRO_STEM_ZH_RE`）中，已带多重键的位次不能再算二氢（同一位不得既是酮又是 CH2，P-58.2.1），`_dihydro_stem_kept` 逐位次用 `_h_prefix_atoms` 映射回原子并检查其全部键，带非单键者剔除；剩余位次改写成 `nH-` 形式（`5-oxo-2,5-dihydrofuran-3-yl` → `5-oxo-2H-furan-3-yl`），全部剔除时只留母体名。

### 环系母体名：稠合 / 桥环 / 螺旋

**稠合名 `fused_namer.fused_parent_names(mol, node)`** 由 L2 打包的 `FusedNode` 树生成未注册稠环的 base 名，编号仍走 L4。

`_component_prefix(node)` 取 `node.fused_prefix`，否则「去尾 `e` 加 `o`」/「词干 + 并」（P-25.3.2.2.2）；`_component_numbering` 委托 L4 的 `fused_component_numbering` 按共享原子取向给子环编号；`_fusion_letter(parent_chain, shared)` 要求共享边落在母体外周边界上并返回侧字母 `a`/`b`/…（`_inner_atoms` 先把出现在 ≥3 环的 perifused 中心原子从母体外周链中剔除）；`_fusion_numbers(child_chain, child_labels, parent_chain, shared)` 给出附加组分共享原子的位次，顺序沿母体低位次端到高位次端。`_fused_one` 在一级单环烃附加组分省略数字位次（`fused_omit_numbers`，P-25.3.8.1）时输出 `前缀[字母]`，否则输出 `前缀[数字-字母]`（P-25.3.2）；`_collect_attached` 递归收集附加组分，二级组分前缀排在附着的一级组分前。入口把附加前缀串到根词干前并查 `RETAINED_FUSION_ALIASES` 换取保留名；单节点（无附加组分）返回 None，交由 `_parent_stem_names` 处理。

**桥环 `bridged_namer.py`**：`bridged_parent_names(mol, node, chain)` 产出 `((完整英, 完整中), (裸词干英, 裸词干中))`，不可组装返回 None。`ring_count_prefix(n_rings)` 对 1 环无前缀、2 环给 `bi` / `双`、≥3 环取 `MULT_EN` / `MULT_ZH`（不用 `dicyclo`）。`bridged_body_names(node)` 取环数 = `len(node.descriptor) - 1`、骨架碳数 = `sum(node.descriptor) + 2`、烷词干 `stems.alkane_en` / `alkane_zh`，主体名即 `{环数词头}cyclo[{描述符}]` / `{环数词头}环[{描述符}]`，经 `stem_forms` 出两形态；该主体名不含 'a' 前缀，供组分式螺环复用（P-24.5.2：'a' 前缀须挂在 spiro 之前而非组分名内）。`descriptor_str(descriptor, locant_pairs)` 把次级桥的上标位次对直接续在对应长度数字之后（`2.2.1`、`2.2.2.0,2` 一类），各长度数字以 `.` 分隔。整名即 `{a}{主体名}`。

**螺旋 `spiro_namer.py`** 处理两类 P-24 母体。

`spiro_parent_names(mol, node, chain)` 处理单环组分螺环（`mono_spiro`，注入 `parent.spiro_node`）：`spiro_multiplier(n_spiro)` 对 1 个螺原子用 `spiro` / `螺`，2 起用 `MULT_EN` / `MULT_ZH` 计数词加 `spiro` / `螺`（`dispiro` / `二螺`，不是 `bispiro`），数量超出词表返回 None，螺原子数取 `len(node.free_spiro_atoms)`；`spiro_descriptor_str(node.descriptor, node.descriptor_superscripts)` 把段长与重访螺原子的上标位次紧接写出、段间以 `.` 分隔，判分口径不带上标花括号；整名经 `stem_forms` 出 `{a}{螺词头}[{描述符}]{烷}` 两形态，'a' 前缀由 `prefix_from_chain(mol, chain)` 给出，烷词干取 `alkane_en(len(chain))` / `alkane_zh(len(chain))`。

`fbs_parent_name(mol, node)` 处理组分式螺环（`fused_bridged_spiro`，注入 `parent.fbs_node`），按 P-24.5.1 拼 `{a 前缀}{螺词头}[组分名-位次对-组分名…]`：

- 前置校验：`node.components` 至少两项、`node.cites` 非空、`node.links` 条数恰为组分数 - 1；螺词头由 `spiro_multiplier(len(comps) - 1)` 给出。每个组分按本组分自身编号渲染：`_numbering(comp)` 取 `comp.numberings[0]`（L4 已选定的组分编号），`_component_name` 出组分全名，`_a_prefix` 出该组分的杂原子 'a' 前缀。
- `_component_name`：`comp.bare_en` 为 None（保留名 / 稠合名）时原样返回 `base_en` / `base_zh`；否则由 `_unsat_locants` 用 `kekulized(mol)` 取组分环内双 / 三键的**本位次**（沿编号方向取较低端），有不饱和时以 `_KIND_TABLE["alkane"]` 的裸词干交 `_chain_names` 生成，无可接裸词干时回落 `base_*`。`_a_prefix` 只对 `comp.kind == "bridged_ring"` 的组分出前缀（保留名 / 稠合名已含杂原子），位次取 `num.locants`，经 `skeleton_replacement_prefix(by_z, "'" * comp.index)` 带上该组分的撇号；各组分前缀以 `-` 相连后置于螺词头之前。
- `_cite_names(cite, names)`：单项直接取组分名，多项用 `BIS_EN` / `BIS_ZH` 出 `bis(...)` / `双(...)`（P-24.7.1 / 24.7.2），倍数超出词表返回 None。
- 相邻引用项之间写螺位次对：由 `links` 建 `槽位 → (父组分, 螺原子)` 映射，每对写成 `{父组分位次}{父组撇号},{本组分位次}{本组撇号}`，同一引用项内多对以 `:` 连接、项间以 `-` 连接；父槽位须已在更早的引用项中出现，否则失败。

### 骨架置换前缀 `skeleton_replacement.py`

`prefix_from_chain(mol, chain)` 由骨架原子序表产出 P-23.3.1 的 'a' 前缀（位次 = 原子在 `chain` 中的序号 + 1）：非碳原子按元素归组，无杂原子返回两个空串、词表外元素返回 `(None, None)`。`skeleton_replacement_prefix(locants_by_z, marks="")` 按元素分组渲染：组内位次升序、逗号连接并加数量前缀，组间按 `P145_SENIOR` 的引用顺序（F > Cl > Br > I > O > S > Se > Te > N > P > …）用连字符连接，得到 `4-thia-1-aza` / `4-硫杂-1-氮杂`；`marks` 为位次撇号后缀（组分式螺环第 k 个组分传 k 个撇号，P-24.6），英文数量词尾 `a` 仅在后接元素名以 `a` 开头时省略（`tetraza` 对 `tetraoxa`）。`A_PREFIX_EN` / `A_PREFIX_ZH` 是唯一的元素词表，同时服务桥环 / 螺环母体名与生成式词干。

### 酯、含氧酸与硫代酯

**`join_oxoacid_name(pre, names, numbered)`** 产出「取代前缀 + O-侧臂 + 功能母体词尾 + 金属盐」整名：`tail_en = join_parent_name(pre[0], names[0])`，中文侧以 `zh_1h_parent` 处理后同法拼接；取 `_o_side_arms(numbered)`（`substituents` 中 `o_side` 为真的记录），`oxo_kind == "phosphate"` 时走 `_fenced_arms_phosphate`、其余中心母体走 `_fenced_arms`；交 `_join_o_side_arms(arms, group=True, arm_zh_fn=_oxoacid_arm_zh)` 拼接臂名（同基按 `_mult_rows` 归并）。末段分四种情形：无臂且有金属 → 酸式盐（英文金属前置、中文金属直接缀于酸式词后，如磷酸二氢钾）；无臂且 `n_om > 0` → 游离根（负电荷不标注，中文补「根」）；无臂 → 纯酸；有臂且有金属 → 酯盐（英文 `{metal} {arms} {tail}`，中文 `{tail}{arms}酯 {metal}盐`）；有臂无金属 → 英文 `{arms} {tail}`，中文 `{tail}{arms}酯`。臂的中文词由 `_oxoacid_arm_zh` 渲染：简单基去「基」（`甲基` → `甲`），全数字词干补「烷基」，含连字符或括号的复合名原样保留。

**`join_ester_name(pre_en, pre_zh, names, numbered)`** 以酸侧为母体、O 侧臂作前缀。中文 O 侧渲染由 `_zh_alkoxy_part` 决定：简单烃基省「基」（乙酸乙酯、十八酸苄酯），仅带位次的基保留「基」不加围栏（乙酸噻吩-3-基酯），其余按需加圆括号、已含围栏时升为方括号。异名臂依次平铺，同名臂用 `di` / `tri` 或复合前缀 `bis` / `tris`（`_join_o_side_arms`，P-16.3.2）。

**硫代羧酸 S-酯**：`_is_thio_side(numbered)` 读 `parent.thio_side`。为真时 `_names_for` 把 `ester` 条目改写为 `coda="ane"`、`en_suf="thioate"`、`zh_suf="硫"`、`ene_base=("enethioate", "烯硫")`、`yne_suf=("ynethioate", "炔硫")`、`variant=None`（P-65.6.3.3.7.1 无 C1/C2 保留名），scaffold 为 `benzene` 时改取 `_benzene_retained("benzenecarbothioate", "苯硫代甲酸")`。臂名经 `_fenced_arm_en` / `_fenced_arm_zh` 围栏（前导位次或自带括号的复合名须括起，内含圆括号时升方括号），整名写作 `S-<臂> <母体>` / `S-<臂><母体>酯`，位次符号在最前；无臂时退化为 `S-<母体>`。

**O 侧臂拼接 `_join_o_side_arms(arms, *, group, arm_zh_fn)`**：单臂时英文名与经 `arm_zh_fn` 渲染的中文名直接返回；`group=True`（含氧酸中心母体）时按 `_mult_rows` 以 `alpha_order_key` 分组计数，同名臂前拼 `MULT_EN` / `MULT_ZH`，异名臂以空格分隔、中文直接相连；`group=False`（酯）时同名臂若任一记录带 `paren` 则英文用 `BIS_EN` 复合前缀并括起（`bis(...)`）、否则用 `MULT_EN` 平铺，异名臂依次平铺（`methyl ethyl oxalate`）。

**环外主基 `_exo_ring_spec(spec, n, numbered)`**：依据 `principal_expression_facts.relation == "exocyclic"` 改写 spec，后缀取自 `EXO_RING_SUF[kind]`（singular 单取代、plural 多取代基底、None 表示该主基无多取代系统名）。环外硫代羧酸 S-酯（`ester` 且 `thio_side`）后缀改用 `("carbothioate", "硫代甲酸")`；环外酰卤的卤素词按 `parent.hal_z` 由 `HALIDE_EN` / `HALO_ZH` 动态拼接（如 `carbonyl chloride` / `甲酰氯`）；多取代而 plural 为 None 时放弃，苯环单取代由 `variant` 保留名承担、不改写 spec；`scaffold_id` 非 `carbocycle`、或已注入 `fused_tree`、或 `stem_generated` 为真时词干已注入，位次恒显式。纯碳环另按环内不饱和分派：有不饱和键时用段式后缀（`coda="ane"`、`cyclic=True`、`zh_loc_omit=False`，段形由 `_unsat_seg` 按后缀推导），饱和时取完整氢化物词干 `cyclo{alkane}` / `环{烷}` 并关闭 `cyclic`；`omit_rule` 为「`omit` 或本环无位次取代基」（`_ring_prefix_located` 判环上是否有非 O-侧取代基），故干净环省略主基位次（`cyclohexanecarboxylic acid`），环上另带被编号前缀时保留（P-66.6.1）。

### 杂原子锚点自由基 `_mononuclear_radical_names`

`kind == "radical"` 且 `parent.radical_anchor_element` 存在时由本函数出口，`stem_en` / `stem_zh` 取 `parent.stem_*`（未注入即失败），按取代基数分派。

- 无取代基：取 `MONONUCLEAR_ZERO_YL[(stem_en, stem_zh)]`。
- P 酰基词干（`stem_en in PHOSPHORYL_STEMS`）：转 `_phosphoryl_sub_names`（P-67.1.4.1.1.5）。已被夺去酸式氢的 O⁻ 臂经 `_oxido_arm` 改写为 `oxido` / `氧化`；同基倍增走 `MULT` 前缀；含自身带括号的复合组分时逐组分以连字符相接、`_bracket_bridge_suffix` 出方括号（`-yl` 型桥基把桥后缀留在括号外），中文补「基」。
- 单取代基：依次试 sulfamoyl 融合、`_fused_bridge_name` 双原子桥合一、`_bridge_enclosed_names` 桥围栏；`azane` 中心另走 `AMIDO_RETAINED` 保留式、内层括起加 `amino`、复合前端括起加 `amino` 三条支路；最后对 `f"{sub}-{stem}"` 调 `free_to_yl`，`oxidane` 中心再经 `_retained_alkoxy` 收拢为保留烷氧基。
- 多取代基：仅 `azane` 中心的双取代基与单核阳离子词干（`stem_en in CATION_STEMS`）的 2–3 个烃基臂成立（O/S 双烷基非标准自由基，明确失败），按英文名排序后分三路：同基倍增；N-芳基-N-某基取 anilino（P-62.2.1.1，阳离子无此保留名）；首基平铺、其余各基加括号（P-62.2.2.1）。

**azane 支路**：`_azane_acyl_stereo_lead` 判单 N-酰基残基名是否带立体前导且词干为 `oyl` / `carbonyl`，`_azane_sub_needs_paren` 判单取代基是否带 `paren` 且尾缀属 `AZANE_PAREN_SUF`（`benzoyl` / `carbonyl` / `acetyl`），命中即经 `_enclose` 内层围栏后缀 `amino`；`_azane_front_needs_paren` 处理前端自带多个位次段（`SUB_LOCANT_RE` 命中 ≥2）的复合取代基（P-63.2.2.1.1），命中即前端圆括后再缀 `amino`——酰基前端（`_ACYL_FRONT_RE`）、苯基前端（走 `anilino` 保留式，P-62.2.1.1）与已自带括号的前端不括。

**保留名与桥合一**：`_mononuclear_en` 在基为 `phenyl` 且后缀为 `sulfonyl` / `sulfinyl` 时取 `benzenesulfonyl`（P-66.4.1），`_mononuclear_zh` 对应走 `zh_bridge_root`；`stem_en == "sulfonyl"` 且取代基名以 `anilino` / `苯胺基` 结尾时把 `anilino` 还原为 `phenyl` / `苯基` 并入 `sulfamoyl` / `氨基磺酰基`（裸 `anilino` 免括，带取代基者经 `_enclose` 括起，P-66.1.1.4.2）。`_fused_bridge_name` 按 `BRIDGE_FUSION_YL` 把前端尾 `imino` / `sulfanyl` 与中心 `azane` / `sulfane` 合一为 `diazenyl` / `disulfanyl`（二氮烯基 / 二硫代基，P-68.3.1.3/.4），优先于桥围栏路径。

**桥围栏**：`_bridge_enclosed_names(a, stem_en, stem_zh)` 处理两类情形——P 酰基经 O/N/S 桥连母体时由 `MONONUCLEAR_BRIDGE` 给方括号；`sulfinyl` / `sulfonyl` 且取代基为复合名（`a["paren"]`）时经 `_enclose` 升级方括号，直链 `-yl` 前端（`_SIMPLE_CHAIN_YL_RE`）与桥融合平铺而不拆。围栏在取代基词条定形时一次产出双语，L5 前缀渲染不重复拆分，命中时置 `numbered["bridge_self_enclosed"] = True`（`azane` 中心除外）。

**free→yl 转换**：`free_to_yl(en, zh, attach_locant, *, paren=True)` 是 P-63.2.2 醇 / 胺的公共出口，`_fg_prefix` 要求中英文同时命中单核氢化物名（词表取 `MONONUCLEAR_HYDRIDES`，含氧化烷 / 氮烷 / 硫烷与 `azanium` / `oxidanium` / `phosphanium` / `sulfanium` 阳离子词干——词表本体 `MONONUCLEAR_CATION_SYSTEMATIC` / `CATION_RETAINED_NAMES` 与 `cation_parent_names` 在 L1/L2 侧，L5 判据是 `CATION_STEMS` / `CATION_YL_STEMS`——以及亚磺酰 / 磺酰 / 亚胺 / 磷酰 / 磷烷基，`PHOSPHORYL_STEMS` 显式排除 `phosphanium`）。带环取代基的 `anilino` 与复合 `amino` 前缀（如 `methylamino`）须加括号，裸 `amino` / `anilino` 免括（P-29.3.6）；阳离子词干的烃基尾「基」保留（`dimethylsulfonio`，P-73.1.1）。

### 前缀组装与围栏 `assembler_prefixes.py`

入口 `_prefix_for(numbered, kind, n)` 从 numbered 取上下文（`scaffold_id`、`double_bond(s)`、`mol`、`locant_kind`）后委托 `_build_prefix`。两种情况直接返回空前缀：`parent.subs_consumed` 为真（取代基已并入 C1 保留名），或 `kind == "radical"` 且 `radical_anchor_element` 存在（烷基取代基已并入组装名）。

**`_build_prefix` 流程**：`_fence_o_side_arms` 就地围栏 O 侧臂后滤除 `o_side` 记录（O 侧臂由 `join_kind_name` 消费），无剩余即返回空串对；`_sub_root_in_ring` 给每个取代基打 `ring_yl` 标记（根原子是否在环上，供 `_stem_needs_paren` 免围栏）；`_omit_sub_locants` 判定位次可否省略（母体 kind 取 `locant_kind or kind`）；`bare = kind == "cation"`（单核母体阳离子，臂名直接前置于阳离子名，无母体位次可混，经 `cation_arm_bare` 打 `bare` 标记）；`flat = kind in OXO_CENTER_KINDS` 时用 `oxo_arm_fence` 重算每个取代基的 `paren` 标记，使臂名围栏由环基判据决定并去掉多余位次围栏；`_mult_rows` 按英文名分组计数，排序键统一为 `alpha_order_key`（P-14.5 字母数字序）；最后算 `bracket = omit and len(groups) >= 2 and (bare or _groups_simple(groups))`（`bare` 时再算 `_cation_arm_hyphen` 决定复合臂间是否以连字符分段），据此选择 `sep`（`""` 或 `"-"`）并调用 `_collect_parts`。

**位次省略 `_omit_sub_locants(n_carbons, substituents, kind, scaffold, has_ene)`** 前几条优先判据：`kind` 属 `N_LOCANT_KINDS`（`urea` / `thiourea` / `guanidine`）时恒不省（`1,3-二甲基脲`）；`kind == "carbamic_acid"` 时恒省（`dimethylcarbamic acid`）；自由基母体仅单碳链省位次；`kind == "cation"` 恒省（取代基全挂在同一原子上，P-73.1.1）。其后按母体情形分派：`n_carbons <= 1` 省位次，但 N- 型与 C- 型取代基共存时 C 侧须带位次；`alkane` 母体在 `carbocycle` / `benzene` scaffold 上单取代且环无不饱和键时位次隐含（环烯除外）；酰胺仅在取代基全为 N- 型时省；`acyl` / `ketone` / `acid` 恒不省；`kind in OXO_CENTER_KINDS` 恒返回 True；含显式括号或自带位次（名内出现数字）的复合取代基不省；`n_carbons == 2` 的单取代仅在端碳无 H 的母体成立（`alcohol` / `amine` / `thiol` / 磺酸族除外）。

**通用臂围栏 `oxo_arm_fence(name, sub, mol)`** 是含氧酸中心母体臂名的判据，命中任一条件即须围栏：臂名已自带 `[...]` 完整围栏 → 不加；前导 `(` 仅为立体描述符（`(2S)-…`）→ 加；去掉立体描述符后仍有 ≥2 个位次段（`_LOCANT_RUN_RE` = `\d+(?:,\d+)*[a-z]*-`，`2,3-` 只算一段）的复合臂名 → 加（P-16.5.1.3.1）；复合基（`sub["paren"]`）且臂根原子在环上 → 加（臂根取 claim 原子中与 `attach_idx` 成键者）。`_fenced_arm` 命中时调用 `_enclose`（含圆括号则升为方括号），`_fenced_arms` 对整组臂中英文同步改写，服务非磷酸的其余中心母体。

**磷酸酯专用围栏 `_fenced_arm_phosphate`**：臂名已自带括号或方括号 → 不加；简单基（`sub` 无 `paren`）→ 平铺；复合臂名再按 `oxo_arm_fence` 或位次段数 > 1 判定是否整体括起。`_fenced_arms_phosphate` 先检查 `_all_arms_stereo_lead`：臂数 ≥2 且全部臂名匹配 `_STEREO_LEAD_ENCLOSE_RE`（`(1Z)-`、`(2R,4R)-` 等前导立体描述符）时整体围栏（P-16.5.1.3.1 第二臂规则仅对多臂生效，单臂不受支配）。

**O 侧臂就地围栏 `_fence_o_side_arms`**：`_O_SIDE_ARM_FENCE_KINDS = {"sulfonate", "phosphate"}`，仅这两类中心母体走此路径（assembler 侧对这两类不判臂围栏）。`_o_side_arm_fence(name, sub)` 判据：`sub` 非复合、臂名不带前导 `[` / `(` 且不含任何括号、且位次段数 ≥2 → 加围栏。处理时置 `arm_fenced` 标记防二次围栏；附着原子为硫（原子序数 16）的硫代酯 S-侧臂跳过，其围栏已在 assembler 侧完成。

**`flat` 与 `_stem_needs_paren`**：`flat` 贯穿 `_build_prefix` → `_collect_parts` → `_parts_for_stem` → `_prefix_one_en` / `_prefix_one_zh`，在 `kind in OXO_CENTER_KINDS` 时为真。英文 `_stem_needs_paren(stem, subs, omit, flat)` 判词干是否须围栏：取代基带 `paren` 时，词干尾为单核阳离子去氢名（`_CATION_YL_TAIL_RE` = `(onio|onium|inium)$`）则免括（`dimethylsulfonio` 整体一词），否则加；`_ring_fused_lead` 判「环基 + 连接组分」式复合前缀（前导指示氢 `_INDICATED_H_LEAD_RE`，或数字开头且含 `-yl`），母体位次省略时不加；词干前导数字时返回 `not flat`，若母体位次已省略且该组带 `ring_yl` 标记也不加（P-16.5.1.2）；余下情形在母体位次未省略且词干为 `trifluoromethyl` 或匹配 `_STEREO_LEAD_ENCLOSE_RE` 时加。中文 `_prefix_one_zh` 的判据同步放宽。

**括号式（P-16.5.1.3.1）**：适用条件为 `omit and len(groups) >= 2 and (bare or _groups_simple(groups))`，`_groups_simple` 要求全部词干无前导位次且非 N- 类。`_collect_parts` 中首个引用的取代基从不加围栏（`bracket and not _LOCANT_RE.search(stem)` 时把该组 `paren` 清零），自带位次者除外；第 2 个起交 `_arm_tail`，普通场合出 `倍数(词干)`，阳离子臂 `bare` 时按 P-73.1.1 改直连式。

**桥后缀拆分**：`_split_bridge_suffix(stem)` 对 `BRIDGE_SPLIT_SUFFIX_EN`（`oxy` / `sulfanyl` / `amino` + 双原子桥 `diazenyl` / `disulfanyl`）逐个尝试，要求前端以 `_FRONT_TAILS`（`yl` / `sulfanyl` / `amino`）结尾且 `_front_needs_enclosure` 为真。该判据逐层排除不拆的场合：`…sulfonyl` / `…sulfinyl`（围栏由 L3 定形）、苄基型前端（`_BENZYL_TAIL_RE`）、`amino` 桥的酰基前端（`_ACYL_FRONT_RE`，按取代式融合）与 `amino` 桥的一般场合；前端含 `[` / `(` 时按端基链型（`_TERMINAL_CHAIN_YL_RE`）与位次取代基尾（`_LOCANT_SUBST_TAIL_RE`）细分；其余以 `-\d+-yl$` 且非 `_SUBST_CHAIN_YL_RE` 判拆。中文侧 `_split_bridge_suffix_zh` 与英文侧同步判定，并在「基」被上游切掉时经 `_zh_front_ji` 在拆分点补回。`_bridge_body(base, suf, merge)` 把前端包上围栏、桥后缀留在括号外（P-63.2.2.1.2），双原子桥整段再括一层，`merge` 用于前端已含方括号且后续还接别的前缀的场合；`_sbridge_flat_stem` 另把「磺酰基 / 亚磺酰基 + 直链 `-yl` 前端」判为英文平铺（P-63.2.1）。

相关正则：`_STEREO_LEAD_ENCLOSE_RE`、`_SIMPLE_CHAIN_YL_RE`、`_SUBST_CHAIN_YL_RE`、`_TERMINAL_CHAIN_YL_RE`、`_BENZYL_TAIL_RE`、`_LOCANT_RE`、`_LOCANT_RUN_RE`、`_ACYL_FRONT_RE`、`_LOCANT_SUBST_TAIL_RE`、`_MULT_WRAP_RE`、`_CATION_YL_TAIL_RE`、`_INDICATED_H_LEAD_RE`。其中 `_LOCANT_RUN_RE` 是 `tools/re.SUB_LOCANT_RE`（本层 `assembler` 与 `tools/anchored_table.carbamoyl_prefix_name` 共用）的无段首锚点变体，`_STEREO_LEAD_ENCLOSE_RE` 与 `tools/re._STEREO_LEAD_STRIP_RE` 分别服务围栏与排序，两者语义不同、不可互换。

**倍数与 N- 前缀**：`_mult_of(lang, stem, subs, n)` 在词干含 `carboxy` / `羧` 哨兵或该组任一取代基带 `paren` 时用 `bis` / `tris`（`_COMPLEX_MULT_LANG`），否则用 `di` / `tri`；`_MULT_WRAP_RE`（`-\d+-yl$`、`oyloxy$`）命中或阳离子去氢前缀（`CATION_YL_STEMS`）组内多于一个时改用 `BIS_EN` / `BIS_ZH` 并加括号。`_place(mult, s, subs, omit)` 在省略位次时直接相接（倍数前缀不得与位次数字直连，必要时加括号），否则以位次串加连字符前置（`_locant_str` 把 N- 型取代基渲染为字母位次 `N`，其余按 `locant_str_sort` 排序）。`_parts_for_stem` 在整组取代基均为 N- 型（`N_PREFIX_KINDS = {n_alkyl, n_block}`）时改走 N- 计数通道：`_n_prime_map` 按字母序给不同 N 原子分配撇号个数（`N,N-` 同氮、`N,N'-` 跨氮），由 `_n_prefix` 拼接。

### 立体化学 `stereo.py`

本模块产出名称最前的 `(…)-` 立体块，对 E/Z 与 R/S 共用一套解析与排序设施。

- **切分与解析**：`_split_stereo_lead` 切出前导 `(…)-` 块；`_parse_stereo` / `_parse_token` 解析为 `(位次|None, 字母)` 部件列表（位次可为 `None`：单原子链只写字母，取代基侧由 L3 按同一格式约定产出）；`_format_stereo` 以 `_part_key` 排序（无位次者在前，余按 `locant_key`）后拼回前缀。
- **E/Z**：`_bond_stereo` 读 RDKit 双键 `GetStereo()`；`_ez_prefix` 对单双键取名并补母体链最小位次（`_bond_min_loc`）；`_ez_multi_prefix` 对 `double_bonds` 收集所有已定义立体的键；`ez_for_parent` 按有无 `double_bonds` 分派。三者经 `_Chain.ez_ene` / `ez_ene_multi` 挂到链引擎。外挂双键由 `_exo_ez_parts` 收集（只取一端在母体内者，位次取母体侧），`join_ez_prefix` 为对外入口，先剔除母体名已带位次再补写，`ester` 时经 `_ester_en_rs` 把部件插到烷基词与酰基词之间。
- **R/S**：`_cip_on_chain` 调 `assign_cip` 后逐链原子读 `_CIPCode`；`_chain_locant` 把链序号映射为位次标签；`_collapsed_parent` 在母体 C1 由更大或环状分子折叠而来时跳过；`_rs_parts` 汇总部件（母体 kind 须属 `_RS_KINDS`，或本身带 `scaffold_id` 的环母体）。`join_rs_prefix` 为对外入口，`_with_rs` 把 R/S 部件并入名称已有立体前缀（保留非 R/S 部件并按位次重排），`ester` 时经 `_ester_en_rs` 插入。

### 盐后缀与阴离子

`join_ring_cation_suffix(numbered, names)` 处理环内 N⁺ / O⁺（P-62.4.1）：净电荷为负的分子直接跳过（盐与两性离子都可能含环阳离子，净负时不处理），名称已含 `ium` 时不重复；在 `chain` 内找出带 +1 形式电荷的 N / O，取其位次标签（`numbering_scaffold.labels`，长度不符时回落链序号）。词干以 `ene` 结尾时按色烯型氧鎓处理（`chromene` → `chromenylium`）；其余经 `_cation_insert` 在英文母体名中定位词干，跳过连接元音 `a`，词干末为 `e` 时吞掉 `e`，再插入 `-{位次}-ium`。中文侧 `_zh_ring_cation` 在母体词干后插 `-{位次}-鎓`（词干自带指示氢前缀时按去前缀词干定位），名中已含「鎓」则不重复。

`join_anion_names(numbered, en, zh)` 在 `parent.anion` 为真时把羧酸后缀转为阴离子式：`acid_to_anion_en` 把 `-oic acid` → `-oate`、`-ic acid` → `-ate`，中文补「根」（已带「根」则不重复）。

`stems.py` 另提供 `join_metal_salt_names(numbered, en, zh)`，由顶层 `namer` 在含氧酸中心母体之外调用（P-71.2 / P-71.3），消费 L0 依 `METAL_ION_EN` / `METAL_ION_ZH`（含 Li / Na / K / Mg / Ca）与 `HALIDE_ZH` / `HALIDE_HX_EN` / `HALIDE_HX_ZH` 产出的 `salt` 记录：`salt.metal` 存在时经 `_metal_en_prefix` / `_metal_zh_suffix` 加金属数量前缀（`_metal_prefix` 对 `n == 1` 不加数量词、`n > 1` 取 `MULT_EN` / `MULT_ZH`），`salt.n_org > 1` 时英文用 `bis(...)`、中文用 `双(...)`；`salt.acid_salt` 与 `salt.halide` 分别缀酸式盐词与卤化物盐 / 氢卤酸盐词。`HALIDE_EN` / `HALO_Z` 是另一套词表，供 `chain_engine` 的酰卤 spec 按 `parent.hal_z` 取卤素词。

**词干表**：`_ALKANE_EN_BASE` / `_ALKANE_ZH_BASE` 覆盖 C1–C10 保留形，`_en_stem(n)` 去尾 `ane`（≥11 去尾 `a`，如 undec / icos），`alkane_en` / `alkane_zh` 产出全名，`zh_num(n)` 对 1–10 用天干（甲…癸）、11 以上复用 `zh_numeral`，`zh_stem` 剥掉末端官能团后缀（十一烷 → 十一）。

## 与其他层的契约

- 自 L1：`mol`、`constants.MONONUCLEAR_HYDRIDES` 等单核母体词表、官能团类别枚举（`stereo._RS_KINDS` 由它派生）、`ring_systems.kekulized` / `build_ring_systems` / `sssr_rings`。
- 自 L2：`kind_registry.pack_parent_stem` 注入 `stem_en` / `stem_zh`；螺环 scaffold 的 `kind` 直接取 scaffold id，由 `_SPIRO_KINDS` 承接。
- 自 L3：`meta.bridge_self_enclosed` 标记 S 桥复合前端名已自含围栏（L5 不整体加括号）。
- 自 L4：`orient_numbering` 对 `fbs_nodes` / `spiro_nodes` / `bridged_nodes` 分别短路到 `layer4.spiro_numbering.fbs_numbering` / `layer4.spiro_numbering.spiro_numbering` / `layer4.bridged_numbering.bridged_numbering`，无论成败都不下落（故 P-25 稠合编号不会误吞螺环与桥环），并把选中的节点写回 `fbs_node` / `spiro_node` / `bridged_node`。
- 本层不 import L2：环系节点只按鸭子类型读字段——`bridged_node` 读 `descriptor` 与 `locant_pairs`；`spiro_node` 读 `descriptor`、`descriptor_superscripts`、`free_spiro_atoms`；`fbs_node` 读 `components`（`index` / `kind` / `atom_ids` / `base_en` / `base_zh` / `bare_en` / `bare_zh` / `numberings`）、`links`（`(螺原子, 父槽位, 本槽位)`）、`cites`。候选由 L2 给出、L4 选中写回，须与本次编号自洽。
- 本层回注：`fused_namer` 依赖 L4 的 `fused_component_numbering`；`stereo` 依赖 L4 的 `assign_cip`、`locant_key`、`locant_str_sort`。
- 对外：`namer` 只调用 `layer5.assembler.assemble`；含氧酸中心母体之外的结果另由 `namer` 调 `join_metal_salt_names` 追加盐后缀（`assemble` 的 `meta.parent_kind` 属 `OXO_CENTER_KINDS` 时跳过，盐已在整名内组装）。

相关页面：[[architecture/overview]]、[[architecture/layer2-parent-selector]]、[[architecture/layer4-numbering]]、[[concepts/bilingual-naming]]、[[reference/core-data-contracts]]。
