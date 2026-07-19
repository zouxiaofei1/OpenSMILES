@echo off
setlocal EnableExtensions
chcp 65001 >nul

REM ChemAgent Console — double-click to restart http://127.0.0.1:8765/
cd /d "%~dp0"

echo [restart-console] Stopping any process listening on port 8765...
for /f "tokens=5" %%P in ('netstat -ano ^| findstr ":8765" ^| findstr "LISTENING"') do (
  echo   killing PID %%P
  taskkill /F /PID %%P >nul 2>&1
)

REM Brief wait so the port is released
timeout /t 1 /nobreak >nul

if not exist ".venv\Scripts\uvicorn.exe" (
  echo [restart-console] ERROR: .venv\Scripts\uvicorn.exe not found.
  echo   Run: python -m venv .venv ^&^& .venv\Scripts\pip install -e ".[dev]"
  pause
  exit /b 1
)

echo [restart-console] Starting uvicorn on http://127.0.0.1:8766/
echo   Close this window or press Ctrl+C to stop the server.
echo.

".venv\Scripts\uvicorn.exe" server.app:app --host 127.0.0.1 --port 8766

echo.
echo [restart-console] Server exited.
pause
