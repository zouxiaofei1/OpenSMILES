"""历史命名引擎(tools/namepredict-v2/v3)别名加载与命名封装。

目录名带连字符(namepredict-v2/-v3)不能直接 import, 用 importlib spec 以别名注册进
sys.modules, 使包内相对导入按别名包名正常解析。v2 内部依赖旧版顶层
namepredict.core.capitalization(现役 src 无 core 子包), 用 shim 挂到已加载的
opensmiles 包上指向 v2 自带实现, 使 v2 与 src 能同进程共存。命名结果规整成与 src
NameResult 同构的 dict(en/zh/success/source/time_ms/meta), 供 routes_name 统一再
追加 engine / gold 字段。
"""

from __future__ import annotations

import importlib.util
import sys
import threading
import types
from pathlib import Path
from typing import Any

_ROOT = Path(__file__).resolve().parents[2]

# 引擎 → (工具目录名, sys.modules 别名, 是否需要 capitalization shim)
ENGINES: dict[str, dict[str, Any]] = {
    "v2": {"dirname": "namepredict-v2", "alias": "namepredict_v2", "need_shim": True},
    "v3": {"dirname": "namepredict-v3", "alias": "namepredict_v3", "need_shim": False},
}

_lock = threading.Lock()


def _install_v2_shim(alias: str) -> None:
    """把 namepredict.core.capitalization shim 挂到顶层 opensmiles 包。

    v2 的 cache/common_names.py 里 `from namepredict.core.capitalization import
    capitalize_chemical_name` 依赖旧版结构。现役 src/opensmiles 已无 core 子包;
    此处把函数指向 v2 自带实现, 避免改 v2 源码。real opensmiles 若已导入(srv 启动时
    即导入), 挂在其上即可共存; 否则建空壳父包承载。
    """
    cap_mod = sys.modules.get(f"{alias}.core.capitalization")
    if cap_mod is None:
        cap_mod = importlib.import_module(f"{alias}.core.capitalization")
    np = sys.modules.get("opensmiles")
    if np is None:
        np = types.ModuleType("opensmiles")
        np.__path__ = []
        sys.modules["opensmiles"] = np
    core = sys.modules.get("namepredict.core")
    if core is None:
        core = types.ModuleType("namepredict.core")
        core.__path__ = []
        sys.modules["namepredict.core"] = core
        np.core = core
    cap = sys.modules.get("namepredict.core.capitalization")
    if cap is None:
        cap = types.ModuleType("namepredict.core.capitalization")
        sys.modules["namepredict.core.capitalization"] = cap
        core.capitalization = cap
    cap.capitalize_chemical_name = cap_mod.capitalize_chemical_name


def _package(engine: str) -> Any:
    """线程安全地 alias 加载指定引擎包一次; 之后从 sys.modules 直取。"""
    spec = ENGINES[engine]
    alias = spec["alias"]
    cached = sys.modules.get(alias)
    if cached is not None:
        return cached
    with _lock:
        cached = sys.modules.get(alias)
        if cached is not None:
            return cached
        init_py = _ROOT / "tools" / spec["dirname"] / "__init__.py"
        pkg_spec = importlib.util.spec_from_file_location(
            alias, str(init_py),
            submodule_search_locations=[str(_ROOT / "tools" / spec["dirname"])],
        )
        mod = importlib.util.module_from_spec(pkg_spec)
        sys.modules[alias] = mod
        pkg_spec.loader.exec_module(mod)
        if spec.get("need_shim"):
            _install_v2_shim(alias)
        return mod


def get_namer(engine: str) -> Any:
    """返回该引擎的一个全新命名器实例(与 src get_namer 一致, 无进程级共享单例)。"""
    return _package(engine).SMILESNNamerV2()


def name_result(smiles: str, engine: str) -> dict[str, Any]:
    """用历史引擎命名, 返回 NameResult 同构 dict(en/zh/success/source/time_ms/meta)。

    失败(无法解析 / 引擎异常 / error 标记)以 success=False + 空名返回; 引擎的
    name_type 写入 source, 便于与 src 的 source(iupac) 区分。
    """
    if engine not in ENGINES:
        return {"en": "", "zh": "", "success": False, "source": "error",
                "time_ms": 0.0, "meta": {"error": f"unknown engine: {engine}"}}
    try:
        raw = get_namer(engine).name(smiles)
    except Exception as exc:  # 引擎内部兜底之外的异常(双保险)
        raw = {"error": f"{type(exc).__name__}: {exc}"}
    if not isinstance(raw, dict):
        raw = {"error": f"unexpected {engine} result: {type(raw).__name__}"}
    en = (raw.get("en") or "").strip()
    zh = (raw.get("zh") or "").strip()
    error = raw.get("error")
    failed = bool(error)
    meta = {}
    if failed:
        meta["error"] = str(error)
    elif raw.get("name_type"):
        meta["name_type"] = raw.get("name_type")
    return {
        "en": en,
        "zh": zh,
        "success": bool(en or zh) and not failed,
        "source": ("error" if failed else (raw.get("name_type") or "iupac")),
        "time_ms": float(raw.get("time_ms") or 0.0),
        "meta": meta,
    }
