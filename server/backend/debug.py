"""SMILES → 命名管线的 print 调试入口。

三种跑法: 单个/多个 SMILES 作命令行参数; 无参数时从 stdin 逐行读; --dataset 跟
数据集 JSON 时整表跑一遍(随机顺序, 只展示管线里 print 有输出的分子, 单分子异常
不中断整批)。
"""

from __future__ import annotations

import io
import json
import os
import random
import sys
import time
import traceback

_HERE = os.path.dirname(os.path.abspath(__file__))  # E:\dev\chem\server\backend
_ROOT = os.path.dirname(os.path.dirname(_HERE))     # E:\dev\chem
sys.path = [p for p in sys.path if p != _HERE]
for _p in (_ROOT, os.path.join(_ROOT, "src")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from namepredict.namer import SMILESNNamer  # noqa: E402

# 批量模式: 数据集仍整个跑完(某条 print 可能只在很靠后的分子上触发), 这里的上限
# 只截"展示"——几千个分子全打出来没人看得完, 超过上限的照跑但输出丢掉。
BATCH_PRINT_LIMIT = 200
BATCH_PROGRESS_EVERY = 250


def describe(s: str, r) -> str:
    """把一次命名排成展示用的几行文本(SMILES / EN / ZH, 失败再补 meta)。"""
    lines = [
        f"SMILES : {s}",
        f"  EN   : {r.en or '<fail>'}",
        f"  ZH   : {r.zh or '<fail>'}",
    ]
    if not r.success:
        lines.append(f"  meta : {r.meta}")
    return "\n".join(lines)


def show(namer: SMILESNNamer, s: str) -> None:
    """单个 SMILES: 打印命名结果(管线内部的 print 会先于这里输出)。"""
    print(describe(s, namer.name(s)))
    print()


def run_quiet(namer: SMILESNNamer, s: str) -> tuple[str, object, str | None]:
    """跑一个分子, 把管线里的 print 收进字符串: 返回 (输出, 命名结果, 异常回溯)。"""
    buf = io.StringIO()
    real_stdout = sys.stdout
    sys.stdout = buf
    result = None
    exc = None
    try:
        result = namer.name(s)
    except Exception:
        exc = traceback.format_exc()
    finally:
        sys.stdout = real_stdout
    return buf.getvalue(), result, exc


def dataset_smiles(path: str) -> list[str]:
    """从 benchmark 数据集 JSON 读出全部 SMILES(支持 [{"smiles":...}] 与 ["..."] 两种形态)。"""
    with open(path, encoding="utf-8") as f:
        rows = json.load(f)
    out = []
    for row in rows if isinstance(rows, list) else []:
        s = row.get("smiles") if isinstance(row, dict) else row
        if isinstance(s, str) and s.strip():
            out.append(s.strip())
    return out


def batch(namer: SMILESNNamer, smiles_list: list[str]) -> None:
    """整表跑一遍: 只展示管线里有 print 冒泡(或抛异常)的分子, 最多 200 个。

    跑之前先打乱顺序: 数据集开头是一串甲烷/乙烷, 固定顺序等于只看最简单的那批,
    某个函数的 print 到底要什么样的分子才触发, 随机撒才碰得到。
    """
    total = len(smiles_list)
    order = list(smiles_list)
    random.shuffle(order)
    t0 = time.perf_counter()
    shown = 0
    over_limit = 0
    silent = 0
    fails = 0

    for i, s in enumerate(order, 1):
        out, result, exc = run_quiet(namer, s)
        if exc is not None:
            fails += 1
        # 只看管线里 print 冒泡的分子: 几千个分子的命名结果(含命名失败)全打出来会把
        # 真正想看的 print 淹掉。异常照留——那是崩溃, 藏起来比淹没更糟。
        if not (out.strip() or exc is not None):
            silent += 1
            continue
        if shown >= BATCH_PRINT_LIMIT:
            over_limit += 1
        else:
            shown += 1
            print(f"\n===== [{i}/{total}] {s} =====", flush=True)
            if out:
                print(out, end="", flush=True)
            if exc is not None:
                print(f"----- 异常 -----\n{exc}", end="", flush=True)
            else:
                print(describe(s, result) + "\n", flush=True)

        # 进度单独打一行: 展示是稀疏的, 分组头的序号不足以当进度用
        if i % BATCH_PROGRESS_EVERY == 0:
            print(f"\n##### 进度 {i}/{total} (已展示 {shown}) #####\n", flush=True)

    print(f"\n===== 完成: {total} 个分子, 展示 {shown} 个(超上限略过 {over_limit} 个), "
          f"无 print 跳过 {silent} 个, 异常 {fails} 个, {time.perf_counter() - t0:.1f}s =====",
          flush=True)


def main() -> None:
    namer = SMILESNNamer()
    args = [a for a in sys.argv[1:] if a.strip()]
    if args and args[0] in ("--dataset", "-d"):
        if len(args) < 2:
            print("用法: debug.py --dataset <数据集.json>", file=sys.stderr)
            raise SystemExit(2)
        batch(namer, dataset_smiles(args[1]))
        return
    if args:
        for s in args:
            show(namer, s)
        return
    print("输入 SMILES：")
    for line in sys.stdin:
        s = line.strip()
        if s:
            show(namer, s)


if __name__ == "__main__":
    main()
