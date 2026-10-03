# trim_code 参考：历史考古证据

`git log --numstat -- src` 聚合出的 src/*.py 净减少提交（按净值排序）：

| 提交 | 日期 | 净值 | 文件 | 标题 | 家族 |
|---|---|---|---|---|---|
| 7acd176 | 08-07 | -1712 | 28 | remove useless code | B 函数级 |
| 3712138 | 08-06 | -1442 | 14 | remove files | A 整文件（合并式） |
| a0128de | 08-06 | -1299 | 13 | remove outdated files | A 整文件 |
| 37b9fdd | 08-06 | -1120 | 11 | remove files | A 整文件 |
| 967edd0 | 08-08 | -1117 | 15 | remove files | A 整文件（先禁用后删） |
| ac65db1 | 08-10 | -942 | 13 | remove codes | A/B |
| 97aafbc | 08-14 | -902 | 13 | remove dead code | A/B |
| 74b5fb1 | 08-11 | -842 | 6 | remove files | A 整文件（级联孤儿） |
| a0862aa | 08-06 | -775 | 23 | remove some dead code | B + 删 17 测试文件 |
| 34b4940 | 08-08 | -754 | 7 | remove dead code | A 整文件 |
| 1d4260f | 08-12 | -700 | 14 | move files | D 移动（净值幻觉） |
| 2aac96b | 08-08 | -589 | 5 | remove files | A 整文件（producer 退役） |
| 2e0c5ba | 08-11 | -422 | **1** | remove unused defs | B 死链传播 |
| 7fac90c | 09-13 | -370 | 20 | remove some codes | B + 落地 dead-code-sweep |
| 02cd17f | 10-03 | -84 | 4 | 删除三个「注入 return None 后预测不变」的死函数 | B 变异法 |
| 43627cb | 09-24 | — | — | 恢复被 tail 清理误删的 purine 模板行 | 反面案例 |

## 判定手法的演进时间线

| 时期 | 判定依据 | 证据 |
|---|---|---|
| 08-06 ~ 08-11 | 静态零引用 / cProfile `ncalls`，**无 evidence 落库** | 提交正文只有一句 `remove files`，无分数无回归数 |
| 08-16 | 首次写进正文：benchmark 四项 + pytest passed 全等 | e601ccd `chore(dead-code): remove unused helpers in layer2` |
| 09-05 | 首次点名方法：AST ref-graph + 全量 pytest/benchmark | 7988502（pytest 1701 passed，benchmark 不变） |
| 09-13 | 固化成 skill：静态 + 变异双法 | 7fac90c 提交 `.claude/skills/dead-code-sweep/` 全套脚本 |
| 10-03 | 变异扫描法大规模落地 | 02cd17f |

08-11 那批是**按清单执行**的：498cd26 新增 `_dead_candidates.txt`（318 行，
格式 `相对路径.py::函数名`），a375f4c 用完删除。该格式与
`dead-code-sweep/scripts/ast_scan.py` 的 `key = f"{rel}::{qual}"` 一致。

## 逐案证据

### A：合并式重构（08-06 晚，两段式提交）
删除与新增**跨两个 commit**，间隔 3~10 秒：
- 20:45:10 `b823cc7 remove files -445`（纯 D 10 个文件）
- 20:45:20 `d30f689 test +401`（新增 `cyclo_fg.py`/`fg_helpers.py`/`principal.py`/
  `polycyclic_parent.py`/`fused56_mono.py`，并删掉 `parent_selector.py` 里对应的 import）

函数名在新旧模块都出现，证明是「接管」而非丢弃：
- `spiro_parent.py` + `bridged_parent.py` → `scaffold/polycyclic_parent.py`
  （`_is_simple_spiro`/`_try_spiro_parent`/`_is_simple_bridged`/`_try_bridged_parent` 原样存在）
- `principal_registry.py` + `principal_selection.py` → `principal.py`
- `aliph_fg.py` + `parent_selector_common.py` → `fg_helpers.py`

### A：级联孤儿化
- `0dde978`（08-11 11:54）删 dispatcher `layer5/special_fg_names.py` + `phosphate_names.py`
- `74b5fb1`（08-11 17:41）删它下面的 `carbamate_names.py`/`carbonate_names.py`/
  `cyclo_exo_fg_names.py`/`diester_names.py`/`sulfur_names.py` —— 这 5 个在父提交上
  已经无 importer，是被上一笔孤儿化的。
- 注意 74b5fb1 的 `unsat_acid.py -244` 是 **M（截断）不是 D**，`--stat` 里看不出来。

### A：先禁用后删除
`967edd0` 删的 `leaves/protocol.py` 文件里自带说明：
`"""Aryr leaf-kind enum (leaf matching retired in 446ce69; enum kept for typing)."""`
前一笔 `446ce69` 已把调用点剪断（`_recurse_leaf_kind(...)` 替换为 `return None`），
33 分钟后才删文件。

### B：变异扫描法的完整 diff（02cd17f）
三处删除 + 一处测试同步：
1. `layer3/claimable_block.py`：删 `_canonical_edge`（零调用）、删 `_attach_parents_of`
   并**内联回** `claim_block`，保留原短路顺序 `if not atoms or len(...) != 1`。
2. `layer2/ring_expression_policy.py`：整模块 63 行删除（`RingExpressionPolicy` +
   `_POLICIES` 24 条 + `_POLICY_INDEX` + `supports_ring_expression`），连带去掉
   `principal_expression._scaffold_fields` 的 `typed_ring_expression_supported` 字段写入。
3. `tests/unit/test_architecture_contracts.py`：删 3 条断言 + 改名 1 个测试。

提交正文给出的验证：「pytest FAILED 集合与基线逐条一致（14 条）；merged_v2
en=63.4%(4232/6679) zh=90.2%(675/748) fails=2452 不变」。

### C：结构收敛的三个完整案例
1. `a03bb07`（09-15）analyzer.py **-270**：一族手写谓词 `_is_carboxyl_carbon`/
   `_is_amide_carbon`/`_is_ketone_carbon`/`_is_ester_carbon`/`_is_lactone_carbon`/
   `_is_aldehyde_carbon`/`_is_hydroxyl_oxygen`/`_is_neutral_ester_alkoxy_o`
   → `fg_local_smarts.py` 的 SMARTS 表 + 一个 `match_local_fg(mol)`。
   非局部判据（酰基 heads 联动、磷酸臂回接、整分子纯度）留在 analyzer 后置。
2. `d26c336`（08-06）+212/-541：`ring_core.py` 里成排的
   `_benzene_core`/`_pyridine_core`/`_diazine_core`/`_hetero5_core`/`_azole13_core`/
   `_diazole_core`/`_fused56_mono`/`_fused56_di13` → 单个 `_template_mother(info)`，
   内部只做一次 SMILES 子图同构。新建 `retained_templates.py`（162 行），核心
   `_TEMPLATES = {"benzene": "c1ccccc1", "pyridine": "n1ccccc1", ...}`。
   被删的 `pyridine.py`(40)/`quinoline.py`(35)/`heteroarene5.py`(62)/
   `sat_hetero_one.py`(148) 各自塌缩成表里一行 SMILES。
3. `d783658`（08-06）+53/-768：删手维护的 `_H5_COOH`/`_SAT_COOH`/`_SAT_ONE` 元组，
   注释写 "Spec is authority when present; this table only covers residual"；
   `benzene_names.py` 同步删掉 `benzenediamine_names`/`pyridinecarboxylic_names` 等
   （-230），这些 kind 名在 HEAD 已彻底不存在。
4. `a5df744`（09-13）**-131**：`anchored_table.py` 由 `_build_registry()` 函数体里
   构造大 dict 改为模块级 `_REGISTRY = {...}` 字面量，顺带删掉 `IupacLevel` 枚举
   及 `systematic_en/zh`、`level`、`kind` 字段。

### D：move 的净值幻觉（1d4260f）
按字节统计 src/*.py：
```
94c5d3c（62381d5 前）  348,151 B  85 文件   ← 基线
62381d5                370,977 B  93 文件   ← 复制，+22,826 B
1d4260f                347,150 B  84 文件   ← 删原件
```
两提交合起来真实净减 ≈ 1,001 B ≈ 40 行 = 被删的 `scaffold/__init__.py`
（40 行 re-export shim）。原文件与新文件逐字节只差 import 行。
`git show -M` / `--find-copies-harder` 都无法配对成 `R`，因为目标文件在
**上一提交**就已存在。

## 反面案例

| 提交 | 事件 | 根因 | 教训 |
|---|---|---|---|
| 43627cb | tail 清理连带删掉 `purine` 模板行 → merged 回归 159 条 → 按原行恢复 | 数据表按尾部批量删时串行 | 表增删一行 = 改一整类命名，必须全量逐行验证 |
| 9fe48d2 | 删 `layer1/urea.py` 等后在调用点留 `return not None`/`and not None`/`if None or ...` | 删调用点没有保语义内联 | 内联要保原短路顺序，不塞 `None` 占位 |
| 9039c75 | 同类，但直接把 entries 换成字面量 `None` 再删模块 | 同上 | 同上 |
| 6467515 | 更早的重构删掉 `anchored_table.resolve_name` 的 general/pin 分派 → `anilino` 恒输出系统名 | 删「一个分派分支」，函数还在还被调用 | AST 引用图与「注入 return None」都抓不到，**只有 benchmark 逐行比对能发现** |
| b46cf40 | 标题 `remove dead code`，实际改 22 行金标 | 删除与金标改动混编 | 归因不清 = 放弃验证 |
| 4f98f15 | 一次删 675 个失败用例，「Full suite: 1795 passed, 0 failed (was 675 failed)」 | 删测试过门禁 | 基线从明红变暗红，后续长期漂移到 14 failed |
| a65e2c6 | 新增 `thiochromene` 模板 → `test_acyl` 2 条既有用例失败 → 43627cb 撤除 | 表增行砸既有测试 | 同上，表改动跑全量 |
| 75e27f9 | `_collapsed_parent` 判据收窄 → 挡没了单碳母体 R/S 描述符，ef98f4a 才补回 | 守卫收紧也是一种删除 | 启发式收紧走同样门禁（全量 3883→3970，fixed 87 / broke 0） |
| 6032f66 | 分支 `_mergecheck` 上实现的 `layer5/isotope.py`(229 行) 未进 main | 多 worktree 下代码会丢 | trim 前确认动的是主干 |

### 侥幸没出事（未经验证的大删波）
2026-08-06 ~ 08-14 的大删波（a0862aa/3712138/a0128de/ac65db1/7acd176/97aafbc）
提交信息只有一句 `remove files`，**无 benchmark / pytest 证据**，紧邻后继只有
`test` / `docs`。功能最终没丢，靠的是**后续系统性重构**把 ether、多元羧酸等
改走通用路径——不是删除当时验证过。

同一时期 `a0862aa` 曾向仓库加入 `conftest.py` + `known_failures.txt`(571 行)，
用 `pytest.mark.skip` 跳过已知失败用例；次日 `ff394d0` 又把这两个文件删掉，
机制只活了一天。这条路线与 4f98f15 一样属于「用跳过/删除让门禁变绿」，不要重走。

### 提交标题不可信
`69085d1` 标题 `fix bugs`，实际删 228 行 L5 模块；`d26c336`/`d783658`/`9039c75`/
`7d63894`/`62381d5` 标题全是 `test`，实为抽象收敛 / 去重 / 删死码 / 纯删 / 复制。
**按标题分类必然误判，必须看 diff。**

## 门禁命令输出参考

```
pytest  : 14 failed, 1757 passed, 1 xfailed          （约 8s，-n auto）
benchmark: en=65.8% (4397/6679) zh=91.0% (681/748)
           dual=65.7% (4390/6679) fails=2289         （约 10s，--time）
```
数值随工作区未提交改动浮动，**开工前自己重测一遍**；要固化的是
「FAILED 集合逐条同名」与「逐行预测全等」这两条不变量。

14 条基线失败（与删除区域多有重叠，按名字集合比对才能区分「我删红的」和「原本就红」）：
test_acyl、test_amide::test_guard_no_false_normalization、test_bridged_naming、
test_bridged_numbering::test_numbering_is_deterministic、test_bridged_system、
test_carbonyl::test_formyl_prefix×3、test_prefixes::test_alkyl_alpha_key_tert_pentyl×2、
test_ring_engine×2、test_registry::test_retained_prefix_whole_molecule、
test_stereo_rs::test_non_stereo_leading_parens_unaffected。

## 其它遗留项
`src/namepredict.egg-info/SOURCES.txt` 仍引用已删文件（`parent_candidate.py`、
`ring_expression_policy.py`、`ring_parent.py`、`scoring.py`、`tools/free_to_yl.py`）
——删除未同步打包元数据，重装包时会暴露。
