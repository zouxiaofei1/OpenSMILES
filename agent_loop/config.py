"""Loop configuration for the self-improving agent."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class LoopConfig:
    """Stop/gate knobs and runtime paths for AgentLoop."""

    data_path: Path = field(default_factory=lambda: Path("data/merged_benchmark.json"))
    cwd: Path = field(default_factory=lambda: Path(".").resolve())
    k: int = 5
    target_dual: float = 0.99
    max_iters: int = 100
    drop_tol: float = 0.005
    n_rules: int = 3
    mock_pi: bool = False
    bench_limit: int | None = None
