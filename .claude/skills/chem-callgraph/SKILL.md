---
name: chem-callgraph
description: Use when analyzing the namepredict function call graph — finding hotspots, tracing call paths, hunting dead/uncovered functions, profiling performance, or generating call-chain visualizations. Triggers on 调用链/调用图/热点/性能剖析/找死函数/cProfile/call graph/hotspot, and on questions about which function calls what, cum time vs self time, or layer-level call structure. Commands: tools/callchain_cli.py (CLI) and the web "调用链" tab.
---

# chem-callgraph

分析 **namepredict 函数调用链**：动态 cProfile 采样 N 个 benchmark 分子 → 调用图（节点=函数，带累计/自耗时与调用次数）+ 调用次数/耗时排行 + graphviz 分层 SVG。

**两套入口，同一底层**：
- **命令行** `python tools/callchain_cli.py`（脚本化/批量/无界面）
- **Web** `http://127.0.0.1:8766/` →「调用链」tab（交互探索：缩放/拖拽/阈值滑块/排行展开排序）

底层是 `server/profile_sampler.py`（多进程 cProfile）+ `server/routes_call_graph.py`（nodes/edges JSON + dot→SVG）。CLI 与 web 结果完全一致。

## 何时用 / 不用

**用**：想知道「某函数被谁调/调了谁」；定位耗时热点；查为什么某些分子慢；找死函数（静态可达但采样未覆盖的分支）；按 layer 看内部调用结构；性能对比改动前后。

**不用**：纯 AST 静态调用关系（那是 `tools/call_graph.py`，layer2 专用、无耗时权重）；改代码时（这是分析工具，不是 lint）。

## 快速开始

```bash
# 采样并打印调用次数/累计耗时 Top 20
python tools/callchain_cli.py --n 200

# 生成分层 SVG
python tools/callchain_cli.py --n 200 --svg callgraph.svg
```

## 完整参数

| 参数 | 默认 | 说明 |
|---|---|---|
| `--n` | 200 | 采样前 N 个 benchmark 分子（`data/merged_benchmark.json`，全量 4070） |
| `--pct` | — | 按全量百分比采样（`--pct 100` = 全量，等价 `--n 4070`） |
| `--workers` | 自动 | 多进程并行数；`n < 300` 自动单进程（免 spawn 开销） |
| `--module` | `namepredict` | 只保留含该子串路径的函数节点 |
| `--top` | 20 | 排行条数 |
| `--json` | — | 导出完整 nodes/edges JSON |
| `--svg` | — | 生成 graphviz 分层 SVG |
| `--layer` | — | 仅某 pipeline layer（0-5；`-1`=核心调度/namer 入口）；跨层调用聚合成灰虚线节点 |
| `--floor` | 1.0 | SVG 的累计耗时阈值%（0=全量；调低会增大 SVG，dot 渲染 1.6k 节点约 1.5s） |
| `--agg` | off | 聚合外部依赖（rdkit/stdlib）为命名空间节点 |

## 输出解读

```
采样: 100 分子 | 墙钟 0.3s | CPU 0.2s | 节点 1427 | 边 1744

调用次数 (Top 20):
    1       8,281×  analyzer:35      ← <genexpr> 等高频微调用
累计耗时 (Top 20):
    1    100.6%  (自耗 0.06%)  namer:227   ← name 入口
```

- **调用次数** `ncalls`：函数被调的总次数（高频 = 可 memoize/内联的微优化点，如 `_dbl_o_on`、`_has_double_bonded_o`）
- **累计耗时** `cum%`：含子调用的总耗时占比；**自耗** `self%` = 函数自身纯开销。找瓶颈看 cum，找纯 CPU 段看 self
- **墙钟 vs CPU**：多进程下 `CPU` 是各 worker 累计，`墙钟` 才是实际等待时间

## 常见任务

**找热点（哪个函数最耗时）**：
```bash
python tools/callchain_cli.py --n 200 --top 30
```
历史上 layer1 `analyzer._collect_fgs`（~40%）与 layer2 候选主链 `_run_candidates`（~55%）是两大头。

**找死函数 / 未覆盖分支**：把阈值调到 0 看全量；反向后排序看低频函数；用 `--pct 100` 跑全量覆盖更多稀有分子触发的分支（多进程下约 11s）。
```bash
python tools/callchain_cli.py --pct 100 --svg all.svg   # 全量 2331 节点
python tools/callchain_cli.py --n 200 --floor 0 --top 50
```

**按 layer 看内部调用结构**：
```bash
python tools/callchain_cli.py --n 200 --layer 2 --svg layer2.svg
```
layer 规模参考：L1 官能团 ~128 节点、L2 骨架 ~282、L3 取代基 ~37、L4 编号 ~32、L5 组装 ~21、L0 预处理 2（太小无意义）。

**追踪某条具体调用路径**：
```bash
python tools/callchain_cli.py --n 200 --json out.json
# 然后 jq/python 查某函数的 callers/callees
```

## 概念澄清（易混淆）

- **「分层」= graphviz dot 的布局层**（按调用深度上下排布），**不是** pipeline Layer 0-5。图里自上而下的层是布局算法的产物。
- **layer 参数**（`--layer`）过滤的是**项目流水线阶段**（L0 预处理→L5 组装），跨 layer 的调用会聚合成灰色虚线「外部节点」并标注来源 layer。
- **累计耗时** cum_s 是含子调用的总时间；**自耗时** self_s 是函数自身。`cum_pct = cum_s / total_s`（分子分母同尺度，多进程并行下占比自洽）。

## Web 前端速查

浏览器打开 `http://127.0.0.1:8766/` →「调用链」tab：
- 「力导向图 / 分层 SVG」切换视图；「生成分层 SVG 图」手动出图
- **采样比例滑块**（1-100%，100% = 全量 4070）、**累计耗时阈值滑块**（0-50%，0 = 全量节点）
- 图下方**排行**：调用次数 / 累计耗时，各 Top 20，可展开全部、正/反向排序，随阈值联动
- 滚轮缩放、鼠标拖拽、双击复位（SVG 视图）；节点悬停看详情

## 坑 / 注意

- **全量 `--pct 100` 首次采样约 11s**（多进程），之后后端缓存命中秒回；日常用 `--n 200` 快速迭代
- `--svg` 节点超过 2500 会报错（全量 < 2500，通常不会触发）
- 采样结果按 `src/namepredict/**/*.py` 的 mtime 签名缓存，**改代码后自动失效重采样**
- 输出含中文，Windows 终端若乱码属编码显示问题，实际为 UTF-8
- CLI/Web 运行前需 server 或直接跑命令即可（CLI 不依赖 server 在运行）
