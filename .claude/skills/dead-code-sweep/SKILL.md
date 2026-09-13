---
name: dead-code-sweep
description: 清理 namepredict 死代码。先让用户二选一：①常规 AST 引用图分析（秒级，抓「没人调用」）②变异扫描法（慢，逐个函数注入 return None 跑 benchmark，抓「调用但没效果」）。触发词：/dead-code-sweep、扫死代码、删无用函数、死代码清理。
---

# 死代码清扫

## 第一步：让用户选方法（必须先问）

用 `AskUserQuestion` 问，不要替用户决定：

**选项 1 — 常规静态分析**（秒级）
以 AST 引用图 + 调用图闭包找不可达函数。快、可反复跑，但**漏动态分派**：注册表 / probe 元组 / 装饰器注册 / 契约测试里的成员会被误判死。
适合：先摸清盘子、或做小范围快速排查。

**选项 2 — 变异扫描法**（耗时较长，几百个函数约 20–60 分钟）
给每个函数开头插 `return None`，跑全量 benchmark，**逐行预测不变**即候选无用。逐个函数跑，带早期退出（任一行变了就停）。
优点：抓静态法抓不到的「函数被调用、但结果没人用/被下游覆盖」。
代价：要跑几百到上千次 benchmark；结论依赖 benchmark 覆盖度。

两个都做也行——**两者的差集最有信息量**：静态说死、变异说活 ⇒ 有动态引用链；静态说活、变异说不动 ⇒ 白跑的计算。

## 通用前置（两个选项都要）

1. **建隔离工作树**，全程只在副本上改写，别碰主仓库（项目有自动 commit 进程）：
   ```bash
   python .claude/skills/dead-code-sweep/scripts/setup_tree.py --dest tmp/deadcode
   # 从归档副本建树：--from-archive src.7z
   ```
2. 之后每条命令都带 `SWEEP_ROOT=<工作树绝对路径>`。
3. **确认基线**：先原样跑一次 benchmark 记下分数；变异法还要跑 `sweep_base.py`（跑两遍校验逐行确定性）。
4. **任何 src 相关测量都不能与变异扫描并发**。工作树被反复改写/还原，并发读源码的进程会采到掺杂状态——本次实测出现过「920 vs 干净 991 个函数」的假矛盾。

## 选项 1：常规静态分析

```bash
SWEEP_ROOT=<树> python .claude/skills/dead-code-sweep/scripts/ast_scan.py --json ast.json
```

输出不可达函数清单。**不可达 ≠ 可以删**，必须过下面「删除前门禁」。

## 选项 2：变异扫描法

```bash
S=<仓库>/.claude/skills/dead-code-sweep/scripts      # 脚本目录
T=<工作树绝对路径>                                    # 上面 setup_tree 建的
export SWEEP_ROOT=$T PYTHONDONTWRITEBYTECODE=1

python $S/sweep_base.py    # 逐行基线（跑两遍校验确定性）
python $S/sweep_cov.py     # 覆盖率（一次全量，~1 分钟）→ 给出「从未执行」集合
python $S/sweep_base.py    # 覆盖序列变了，重排探针序（前 10 行即覆盖 96% 被执行函数）
python $S/sweep_drive.py   # 全量扫描，断点续跑，中途可 Ctrl-C
python $S/sweep_verify.py  # 终局确认：候选一次性全注入，跑无早退全量 benchmark
python $S/sweep_analyze.py # 分档 + 扫下游引用 + 生成 $T/dead-code-report.md
```

先 `--limit 5` 冒烟一遍再放开全量。

判读 `_mut_results.json` 的 status：

| status | 含义 | 处置 |
|---|---|---|
| `diff` | 注入后预测改变 | 有用 |
| `same` | 全量扫完预测不变 | 候选无用 |
| `import_error` | 注入后包加载不了 | **导入期必需，有用**，别删 |

## 判定标准（用户口径）

- **零引用的可删** —— 候选里在 `src/` + `tests/` + `benchmarks/` + `server/` + `tools/` 全库零引用。
- **有引用但 benchmark 不变的，不直接删，输出 .md** —— `sweep_analyze.py` 生成的 `dead-code-report.md` 分两栏：①可直接清理（零引用）②需同步修改下游才能删（列出 `tests/` / `server/` / `benchmarks/` 的引用位置）。
- 判定前**必须**跑 `sweep_verify.py`：候选全注入后，全量 benchmark 的**逐行预测与总分都要与基线一字不差**。单函数跑分只管「这一行变了没」，这一步管「全删了总分动不动」。

## 删除前门禁（漏一个就可能删错）

1. `sweep_verify.py` 逐行预测全等 + 总分全等，源码逐字节还原。
2. `pytest` FAILED 集合全等（不是全绿——本项目基线本来就有失败）。**只被测试引用 ≠ 可删**：删了 fail 测试的，按本口径不算死代码，要么同步改测试，要么保留。
3. 整模块未被导入的，先确认是整套旧实现被取代，而不是被 `__init__.py` 或动态 import 延迟加载。
4. 删完**重跑一遍**：删一个函数可能让它的唯一 callee 变死（死链传播），需迭代重扫。

## 硬性约束

- **禁止**改 `benchmarks/benchmark.py` 的计分逻辑，或删/跳过测试来「通过门禁」。
- **禁止**未经允许删 `print`（管线里有遗留 debug print，是既定事实）。
- 临时产物只放工作树（`tmp/`，已 gitignore），不落主仓库被跟踪处。
- 结论必须写清口径边界：**「从未执行」是在这份 benchmark 的输入分布下的结论，不等于绝对死代码**。

细节坑与实现要点见 `reference.md`。
