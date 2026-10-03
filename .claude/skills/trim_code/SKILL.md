---
name: trim_code
description: namepredict 代码瘦身总纲。把 src/ 行数降下来的四路历史手法（整文件退役 / 函数级删除 / 结构收敛 / 移动合并）+ 不可协商的删除门禁 + 反面案例。函数级扫描转 dead-code-sweep。触发词：/trim_code、瘦身、减行数、清理代码、收缩代码、删模块。
---

# 代码瘦身（trim_code）

本仓 src/ 历史上减行只靠四类手法，做法与风险各不相同。先按下表定位，
再照该路执行。**任何一路都必须先过「门禁」。**

| 想删什么 | 走哪路 | 判据 | 历史提交 |
|---|---|---|---|
| 整个模块 / 文件 | 路 A 整文件退役 | 无 importer / 功能已被接管 / 上游退役 | 3712138、967edd0、74b5fb1 |
| 单个函数、分支、常量表 | 路 B 函数级删除 | 零引用；或注入 `return None` 后逐行预测不变 | 2e0c5ba、7acd176、02cd17f |
| 一堆同构小函数、重复清单 | 路 C 结构收敛 | 差异只是常量 → 抽成数据表 | d26c336、d783658、a03bb07 |
| 文件挪位置 / 合并 | 路 D 移动合并 | 注意净值是记账幻觉 | 1d4260f |

判据交叉时按「破坏面」选路：**路 C 最安全也最有效**（改结构不改行为），
路 B 最依赖工具，路 A 连带面最大。

## 门禁（四路共用，不可协商）

1. **先测基线**，改动全程用它比对：
   ```
   .venv/Scripts/python -m pytest -n auto -q --tb=no
   .venv/Scripts/python -m benchmarks.benchmark_parallel --data benchmarks/merged_v2_benchmark.json --time --timeout 1
   ```
   **基线不是全绿**（实测 14 failed / 1757 passed，2026-10-03）。判据是
   **pytest FAILED 集合逐条同名**，不是 0 failed。要求全绿等于逼人删测试。
2. **benchmark 逐行预测 + 总分全等**，只看总分不够。删一个分派分支
   （6467515：`anilino` 退化成 `phenylamino`）或收紧一个守卫
   （75e27f9→ef98f4a：挡没了单碳母体 R/S）时，总分可能微动甚至不动。
3. **删除、改测试、改金标（`benchmarks/*.json`）分三个提交**。混在一起
   就无法归因（b46cf40 标题 `remove dead code` 却夹带 22 行金标改动）；
   也无法区分金标改动是纠错还是「适配我当前的输出」。
4. **禁止删或跳过测试来通过门禁**。4f98f15 一次删掉 675 个失败用例，
   基线从「明红」变「暗红」，随后长期漂移到 14 failed。
5. **「只被测试引用」≠ 可删**。先查那条测试是活测试还是已知坏测试；
   删了不同步测试，FAILED 集合就变。
6. **改测试要改内容，不是只删断言**。
7. **确认动的是主干**，不是某个 worktree / 分支快照。

## 路 A：整文件退役

三种判据，任选其一：
- **零 importer**：全模块路径 grep 无命中，且全库无动态 import
  （`git grep -n "importlib\|__import__"` 为空）。
- **功能已被新模块接管**：函数名在两处都存在。先做对照表
  （`git show <新模块> | grep -E "^def|^class"` vs 被删文件），逐个核对。
- **上游 dispatcher 退役导致的孤儿化**：下游模块的 importer 被删后成片变死。
  这是**级联**，分两轮走——先删 dispatcher，下一轮再删被孤儿化的下游
  （0dde978 删 `special_fg_names.py` → 74b5fb1 删它下面 5 个 `*_names.py`）。

步骤：
1. 判据是「接管」时，**先把实现搬进新模块**，别先删。
2. 同一次提交里 `git rm` 旧模块 **+ 改完所有 import 点**。
3. 模块有专属测试（`test_<module>.py`）就一起删。

硬坑：
- **grep 必须用全模块路径**（`namepredict.layer2.urea`）。basename 会被跨层
  同名误报——`layer1/urea.py`、`layer2/urea.py`、`layer5/urea_names.py` 同名。
- **删文件与改 import 若拆成两次提交，中间提交 import-broken**。历史犯过 4 次，
  最长断链 37 分钟（b823cc7 删文件留 10 处悬空 import，d30f689 才修）。
- 扫**硬编码文件名的测试断言**：`_TARGET_L5` 那类参数化清单配
  `assert path.is_file()`，删文件必红。
- 用 `--name-status` 区分 `D` 与 `M`，`--stat` 的柱状图会骗人：
  74b5fb1 的 `unsat_acid.py -244` 是 `M` 不是 `D`。
- 拿不准就用**先禁用后删除**：把调用点换成中性值 + 注释写明退役提交号
  （967edd0 的 `leaves/protocol.py` 留 "retired in 446ce69"），跑一轮确认无变化，再删文件。

## 路 B：函数级删除 → 转 dead-code-sweep

扫描与判定全部走 `/dead-code-sweep`（静态 `ast_scan.py` / 变异 `sweep_*.py` 两选一，
或两法都跑，差集最有信息量）。本路专属要点：

- **死链传播**：删一个函数常让它唯一的 callee 也变死，需**迭代重扫**。
  2e0c5ba 在单文件里连删 422 行就是滚出来的。
- **删调用点必须保语义内联**，不要塞 `None` 占位。反例：9fe48d2 删掉
  `layer1/urea.py` 后在调用点留下 `return not None`、`and not None`、
  `if None or ...`，隔天 7c73e5b 才修干净。正确姿势见 02cd17f：删
  `_attach_parents_of` 时把它**按原短路顺序**内联回 `claim_block`。
- 判定口径：零引用可直接删；**有引用但 benchmark 不变 → 出报告，不直接删**。

## 路 C：结构收敛（优先考虑）

三招，按收益排序：

1. **代码 → 数据表外移**（本仓最推崇，符合「数据驱动 > 代码驱动」）。
   - if/elif 谓词长链 → SMARTS 表 + 一个匹配器。a03bb07 把 `_is_carboxyl_carbon`
     /`_is_amide_carbon`/`_is_ketone_carbon`… 一族手写谓词换成
     `fg_local_smarts.py`，analyzer.py **-270 行**。表头注释即方法论：
     「新增官能团只加表项，不改检测代码」。
   - 函数体里构造大 dict → 模块级字面量。a5df744 把 `_build_registry()` 改成
     `_REGISTRY = {...}`，**-131 行**，顺带删掉已无人查的列。
2. **抽象收敛**：N 个「只差一个常量」的专用函数 → 1 个通用匹配器 + 表。
   d26c336 把 `_benzene_core`/`_pyridine_core`/`_diazine_core`/`_hetero5_core`…
   换成单个 `_template_mother` + `_TEMPLATES` dict（值就是 SMILES），
   被删的 4 个专用模块的定位逻辑各自塌缩成表里一行。
3. **去重**：多份手写清单 → 单一 registry，并在注释里声明唯一事实来源。
   d783658 删掉手维护的 `_H5_COOH`/`_SAT_COOH` 等元组，注释写明
   "Spec is authority when present; this table only covers residual"；
   8ab9f93 删掉 `fg_producers.py` 这个 bootstrap 壳 + 4 个同构 `_try_*_parent` 分支。

风险：**数据表增删一行 = 改一整类命名，静态引用图完全看不出**。
43627cb 记录 tail 清理时连带删掉 `purine` 模板行，导致嘌呤系母体退化为
生成式名（merged 回归 159 条）。表改动**必须跑全量逐行**。

## 路 D：移动合并

- **净值是记账幻觉**。`git show -M --stat` 里只有 `D` 没有 `R` 的「move」，
  去**前后各一个提交**找对应的 `A`。1d4260f `move files -723`：真实净减
  ≈40 行（只是不再需要 `scaffold/__init__.py` 那层 re-export），另外 660 行的
  「增」落在 5 秒前的 62381d5——copy 在前一提交、delete 在本提交，
  git 无法配对成重命名。
- copy-then-delete 会留下**两份几乎相同的代码并存窗口**（那 5 秒里
  `layer2/` 与 `layer2/scaffold/` 同时存在）。复制与删除**合成一笔提交**。

## 命令速查

```bash
# 找删除类提交与规模
git log --all --date=short --format='%h|%ad|%s' --shortstat | grep -iE 'remove|delete|删'
# 单笔粒度：先看 D/M 分型，再看内容
git show --name-status <sha>
git show <sha> -- src/namepredict/layer2/xxx.py
# 是否有配套脚本 / 报告落库
git show --stat <sha> | grep -viE 'src/|tests/'
# 反复检测（文件级：ADD 计数 >1 即复活）
git log --all --diff-filter=A --oneline -- <path> | wc -l
# 反复检测（符号级）
git log --all --oneline -S'<symbol>' -- src/
# 提交信息里的危险词
git log --all --pretty='%h|%ad|%s%n%b' --date=short | grep -iE '误删|连带|补回|修回|撤除|restore|regress|broke'
```

历史手法演进与逐案证据见 `reference.md`。
