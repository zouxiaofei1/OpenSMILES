"""FastAPI app entry: Namer only.

Run:
  uvicorn server.app:app --host 127.0.0.1 --port 8765
"""

from __future__ import annotations

from fastapi import FastAPI
from fastapi.staticfiles import StaticFiles
from starlette.responses import Response
from starlette.types import Scope

from server.routes_name import router as name_router
from server.routes_benchmark import router as benchmark_router
from server.routes_benchmark_run import router as benchmark_run_router
from server.routes_debug import router as debug_router

app = FastAPI(title="ChemAgent Namer", version="0.1.0")
app.include_router(name_router)
app.include_router(benchmark_router)
app.include_router(benchmark_run_router)
app.include_router(debug_router)


@app.get("/health")
def health() -> dict[str, str]:
    """Liveness probe."""
    return {"status": "ok"}


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
            or base in ("app.js", "namer-ketcher.js", "app.css")
            or lower.endswith("/app.js")
            or lower.endswith("/namer-ketcher.js")
            or lower.endswith("/app.css")
        ):
            response.headers["Cache-Control"] = "no-cache, must-revalidate"
            response.headers["Pragma"] = "no-cache"
        return response


# Static UI last so /api/* and /health win over StaticFiles.
app.mount("/", NoCacheHTMLJSStaticFiles(directory="web", html=True), name="web")
