"""命名过程中的一次性记忆：同一次命名内同一个分子对象只算一次确定性中间结果。

命名同一个分子时，氢化骨架、CIP 标签、环感知、锚定子分子这些量会被不同层反复
问到（实测同一分子 `AtomRings()` 被取十几次、`_hydrogenated` 被重建几十次），
但它们只依赖分子自身，在一次命名内不会变。这里按分子对象身份记忆，纯粹消除
重复计算，不改变任何返回值。

不跨分子共享：跨分子共享会带入宿主相关的立体上下文（同 `namer.run_cache` 的
理由）。每次顶层命名开始时由 `begin_run` 清空。存储按线程隔离，避免并发请求
互相清空对方的记忆。
"""

from __future__ import annotations

import threading

_store = threading.local()


def begin_run() -> None:
    """开始一次顶层命名：清空本次记忆。"""
    _store.data = {}


def _bucket(name: str) -> dict:
    """取本次命名中某个记忆表的 dict，不存在则建。"""
    data = getattr(_store, "data", None)
    if data is None:
        data = _store.data = {}
    return data.setdefault(name, {})


def by_mol(name: str, fn, mol):
    """按 (记忆表名, 分子对象) 记忆 `fn(mol)` 的结果；结果可为 None。

    键用 `id(mol)`，同时把 mol 存进记忆值里保活：本次命名内对象不会被回收，
    id 也就不会复用，因此不会串到别的分子上。
    """
    bucket = _bucket(name)
    entry = bucket.get(id(mol))
    if entry is not None:
        return entry[1]
    value = fn(mol)
    bucket[id(mol)] = (mol, value)  # 首项保活 mol
    return value


def by_key(name: str, key: tuple, fn, *keepalive):
    """按自定义键记忆 `fn()` 的结果；结果可为 None，`keepalive` 里的对象随条目保活。"""
    bucket = _bucket(name)
    entry = bucket.get(key)
    if entry is not None:
        return entry[1]
    value = fn()
    bucket[key] = (keepalive, value)
    return value
