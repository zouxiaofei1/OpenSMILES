@echo off
rem ============================================================
rem 启动 3 个 Claude Code 会话，放进同一个 Windows Terminal 窗口
rem （窗口命名为 chem，内含 3 个标签页 claude-1/2/3）
rem 重复运行：若 chem 窗口已存在，会在其中追加 3 个新标签页
rem ============================================================
setlocal

set "WORKDIR=E:\dev\chem"

wt.exe -w chem new-tab -d "%WORKDIR%" --title "claude-1" --suppressApplicationTitle cmd /k claude ; new-tab -d "%WORKDIR%" --title "claude-2" --suppressApplicationTitle cmd /k claude ; new-tab -d "%WORKDIR%" --title "claude-3" --suppressApplicationTitle cmd /k claude

endlocal
