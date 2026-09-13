"""FastAPI app entry: Namer only.

Run:
  uvicorn server.backend.app:app --host 127.0.0.1 --port 8765
"""

from __future__ import annotations

import importlib
import logging
import time
import traceback
from pathlib import Path

from fastapi import FastAPI, Request
from fastapi.staticfiles import StaticFiles
from starlette.responses import JSONResponse, Response
from starlette.types import Scope

from server.backend.routes_name import router as name_router
from server.backend.routes_benchmark import router as benchmark_router
from server.backend.routes_debug import router as debug_router
from server.backend.routes_layer_benchmark import router as layer_benchmark_router
from server.backend.routes_code_analysis import router as code_analysis_router
from server.backend.routes_call_graph import router as call_graph_router
from server.backend.routes_history import router as history_router

app = FastAPI(title="ChemAgent Namer", version="0.1.0")
app.include_router(name_router)
app.include_router(benchmark_router)
app.include_router(debug_router)
app.include_router(layer_benchmark_router)
app.include_router(code_analysis_router)
app.include_router(call_graph_router)
app.include_router(history_router)


# src/ 健康检查: 导入成功即缓存(进程内已加载的模块不会因磁盘改动失效, 与请求
# 实际取到的代码一致); 导入失败不缓存, 修好后下一次轮询立刻转好。
_SRC_CHECK: dict[str, object] = {}


ROOT = Path(__file__).resolve().parents[2]


def _error_location(exc: BaseException) -> str:
    """取异常里最靠 namepredict 的那一帧, 返回 "src/namepredict/...:行号"(取不到则空串)。"""
    frames = traceback.extract_tb(exc.__traceback__)
    where = next((f for f in reversed(frames) if "namepredict" in f.filename), None)
    where = where or (frames[-1] if frames else None)
    if where is None:
        return ""
    try:
        return f"{Path(where.filename).relative_to(ROOT)}:{where.lineno}"
    except ValueError:  # 帧不在项目内(如 site-packages), 退回绝对路径
        return f"{where.filename}:{where.lineno}"


def src_error() -> str | None:
    """返回 namepredict 的导入失败原因(含出错文件行号); 能正常导入时返回 None。"""
    if _SRC_CHECK.get("ok"):
        return None
    try:
        importlib.invalidate_caches()
        importlib.import_module("namepredict.namer")
    except Exception as exc:
        where = _error_location(exc)
        return f"{type(exc).__name__}: {exc}" + (f"  ({where})" if where else "")
    _SRC_CHECK["ok"] = True
    return None


@app.get("/health")
def health() -> dict[str, str]:
    """Liveness + src 可用性: 前端据此区分「后端没起」和「src/ 编译失败」。"""
    err = src_error()
    return {"status": "degraded" if err else "ok", "src": "broken" if err else "ok",
            "src_error": err or ""}


@app.exception_handler(ImportError)
async def import_error_handler(request: Request, exc: ImportError) -> JSONResponse:
    """src/ 编译失败只影响本请求: 回 500 + 出错文件行号, 而不是裸 traceback。"""
    where = _error_location(exc)
    logging.getLogger("chemagent").error("import error on %s: %s", request.url.path, exc)
    return JSONResponse(
        status_code=500,
        content={
            "detail": f"导入 namepredict 失败, src/ 可能编译不过: {type(exc).__name__}: {exc}"
                      + (f"  ({where})" if where else ""),
            "error": f"{type(exc).__name__}: {exc}",
            "where": where,
        },
    )


class NoCacheHTMLJSStaticFiles(StaticFiles):
    """Disable heuristic caching for shell HTML/JS so fixes take effect on refresh.

    Large hashed vendor assets under /vendor/ketcher/assets/ keep default caching.
    """

    async def get_response(self, path: str, scope: Scope) -> Response:
        response = await super().get_response(path, scope)
        lower = path.replace("\\", "/").lower().lstrip("/")
        base = lower.rsplit("/", 1)[-1] if lower else ""
        if (
            lower.endswith(".html")
            or lower in ("", "index.html")
            or lower.startswith("js/")
            or base in ("app.js", "namer-ketcher.js", "app.css")
            or lower.endswith("/app.js")
            or lower.endswith("/namer-ketcher.js")
            or lower.endswith("/app.css")
        ):
            response.headers["Cache-Control"] = "no-cache, must-revalidate"
            response.headers["Pragma"] = "no-cache"
        return response


# Static UI last so /api/* and /health win over StaticFiles.
WEB_DIR = ROOT / "server" / "web"
app.mount("/", NoCacheHTMLJSStaticFiles(directory=str(WEB_DIR), html=True), name="web")
