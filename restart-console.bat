@echo off
setlocal EnableExtensions
chcp 65001 >nul

REM ChemAgent Console — restart http://127.0.0.1:8666/  (auto-restarts on src/server changes)
cd /d "%~dp0"

echo [restart-console] Stopping any process listening on port 8666...
powershell -NoProfile -ExecutionPolicy Bypass -File "server\backend\stop-port.ps1" -Port 8666

timeout /t 1 /nobreak >nul

if not exist ".venv\Scripts\uvicorn.exe" (
  echo [restart-console] ERROR: .venv\Scripts\uvicorn.exe not found.
  echo   Run: python -m venv .venv ^&^& .venv\Scripts\pip install -e ".[dev]"
  pause
  exit /b 1
)

echo [restart-console] Starting dev server on http://127.0.0.1:8666/
echo   Changes under src\ and server\ restart uvicorn automatically.
echo   Close this window or press Ctrl+C to stop the server.
echo.

".venv\Scripts\python.exe" server\backend\dev_server.py

echo.
echo [restart-console] Server exited.
pause
