"""Loop configuration for the self-improving agent."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass
class LoopConfig:
    """Stop/gate knobs and runtime paths for AgentLoop."""

    data_path: Path
    cwd: Path
    k: int = 5
    target_dual: float = 0.99
    max_iters: int = 100
    drop_tol: float = 0.005
    n_rules: int = 3
    mock_pi: bool = False
    bench_limit: int | None = None
