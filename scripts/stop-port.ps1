param(
    [Parameter(Mandatory = $true)]
    [int]$Port
)

# Kill the uvicorn service tree listening on $Port.
#
# Works for both plain uvicorn (single worker process) and `--reload` mode
# (reloader parent + spawned worker child). For --reload, the listening socket
# is owned by the reloader while the accept loop runs in the spawned worker, so
# killing only the listener would leave an orphaned worker alive and the port
# half-held. We therefore collect the listener, all its descendants, and its
# uvicorn/python ancestors, and kill the whole set.
#
# Ancestor climb stops at the first process whose command line does not contain
# "uvicorn" or "python" (e.g. cmd / explorer / bash), so the launching shell is
# never killed.

$ErrorActionPreference = 'SilentlyContinue'
$all = @()

$conn = Get-NetTCPConnection -LocalPort $Port -State Listen
if (-not $conn) {
    Write-Host "  port $Port not in use"
    exit 0
}
$listen = ($conn | Select-Object -First 1).OwningProcess
Write-Host "  listener PID: $listen"
$all += $listen

function Get-Descendants([int]$id) {
    $kids = Get-CimInstance Win32_Process -Filter "ParentProcessId=$id"
    foreach ($k in $kids) {
        $script:all += $k.ProcessId
        Get-Descendants $k.ProcessId
    }
}
Get-Descendants $listen

# Climb ancestors, stopping at the first non-uvicorn/python process.
$p = Get-CimInstance Win32_Process -Filter "ProcessId=$listen"
while ($p) {
    $parent = Get-CimInstance Win32_Process -Filter "ProcessId=$($p.ParentProcessId)"
    if (-not $parent) { break }
    if ($parent.CommandLine -match 'uvicorn|python') {
        $script:all += $parent.ProcessId
        $p = $parent
    } else {
        break
    }
}

foreach ($id in ($all | Sort-Object -Unique)) {
    if (Get-Process -Id $id) {
        Stop-Process -Id $id -Force
        Write-Host "  killed PID $id"
    }
}
