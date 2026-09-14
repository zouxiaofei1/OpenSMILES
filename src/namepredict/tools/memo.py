"""命名过程内的一次性记忆：同一分子对象的确定性中间结果只算一次。
按对象身份或自定义键缓存消除重复计算；不跨分子共享，begin_run 清空。
存储按线程隔离，避免并发请求互相清空记忆。"""

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
    """按 (记忆表名, 分子对象) 记忆 fn(mol) 的结果。"""
    bucket = _bucket(name)
    entry = bucket.get(id(mol))
    if entry is not None:
        return entry[1]
    value = fn(mol)
    bucket[id(mol)] = (mol, value)  # 首项保活 mol
    return value


def by_key(name: str, key: tuple, fn, *keepalive):
    """按自定义键记忆 fn() 的结果，keepalive 对象随条目保活。"""
    bucket = _bucket(name)
    entry = bucket.get(key)
    if entry is not None:
        return entry[1]
    value = fn()
    bucket[key] = (keepalive, value)
    return value
