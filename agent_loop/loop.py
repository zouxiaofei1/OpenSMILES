"""Agent loop state machine: bench → cluster → pi → lint → pytest → gate."""

from __future__ import annotations

import json
import subprocess
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from agent_loop.config import LoopConfig
from agent_loop.events import EventBus
from agent_loop.git_gate import DEFAULT_ALLOWED, GitGate
from agent_loop.pi_runner import MockPiRunner, PiResult, PiRunner
from agent_loop.state import StateStore
from benchmarks.benchmark import run_benchmark
from tools.fail_cluster import Cluster, cluster_failures
from tools.structure_lint import lint_tree

NamerFactory = Callable[[], Any]
BenchFn = Callable[[], dict[str, Any]]
BoolFn = Callable[[], bool]
LintFn = Callable[[], list[str]]

_MEMORY = Path("agent_loop/memory")
_PROMPT_TMPL = Path("agent_loop/prompts/cycle.md")
_PROMPT_OUT = Path("agent_loop/prompts/.cycle_current.md")
_SKILL = "skills/chem-tdd-skill/"
_NAMEPREDICT = Path("src/namepredict")
_CONSTRAINTS = (
    "no ML/LLM naming; no SMILES special-case; no gold/scoring edits; "
    "file<=500 lines; func body<=10; allowlist paths only; S3 N=3"
)


@dataclass
class CycleResult:
    """One cycle outcome for callers and tests."""

    iter: int
    dual_before: float
    dual_after: float
    decision: str
    session_path: Path


def _default_namer() -> Any:
    from namepredict.namer import SMILESNNamer

    return SMILESNNamer()


def _default_pi(cfg: LoopConfig) -> PiRunner | MockPiRunner:
    return MockPiRunner() if cfg.mock_pi else PiRunner()


def _ensure_paths(cfg: LoopConfig) -> LoopConfig:
    cfg.data_path = Path(getattr(cfg, "data_path", None) or "data/merged_benchmark.json")
    cfg.cwd = Path(getattr(cfg, "cwd", None) or Path(".").resolve())
    return cfg


def _dual(report: dict[str, Any]) -> float:
    return float(report.get("acc_dual") or 0.0)


def _n_fails(report: dict[str, Any]) -> int:
    return len(report.get("fails") or [])


def _is_improve(before: dict[str, Any], after: dict[str, Any]) -> bool:
    d0, d1 = _dual(before), _dual(after)
    if d1 > d0:
        return True
    return d1 == d0 and _n_fails(after) < _n_fails(before)


def _one_cluster_line(c: Cluster) -> str:
    return (
        f"- key={c['key']} size={c['size']} layers={c['suggest_layers']} "
        f"feats={c['features']} samples={c['sample_smiles'][:3]}"
    )


def _cluster_text(clusters: list[Cluster]) -> str:
    if not clusters:
        return "(no failures)"
    return "\n".join(_one_cluster_line(c) for c in clusters)


def _bench_cmd(cfg: LoopConfig) -> str:
    limit = f" --limit {cfg.bench_limit}" if cfg.bench_limit is not None else ""
    return f"python -m benchmarks.benchmark --data {cfg.data_path}{limit}"


def _apply_placeholders(text: str, mapping: dict[str, str]) -> str:
    for k, v in mapping.items():
        text = text.replace(k, v)
    return text


def _prompt_map(
    iter_n: int, dual: float, cluster: str, skill: str, constraints: str, cmd: str
) -> dict[str, str]:
    return {
        "{{iter}}": str(iter_n),
        "{{dual}}": f"{dual:.6f}",
        "{{cluster}}": cluster,
        "{{skill_path}}": skill,
        "{{constraints}}": constraints,
        "{{bench_cmd}}": cmd,
    }


def _render_prompt(tmpl: Path, out: Path, mapping: dict[str, str]) -> Path:
    text = _apply_placeholders(tmpl.read_text(encoding="utf-8"), mapping)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(text, encoding="utf-8")
    return out


def _run_pytest(cwd: Path) -> bool:
    r = subprocess.run(
        ["pytest", "tests/unit", "-q"],
        cwd=cwd,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
    )
    return r.returncode == 0


def _parse_status_line(line: str) -> str | None:
    if len(line) < 4:
        return None
    rel = line[3:].strip().replace("\\", "/")
    return rel.split(" -> ", 1)[1] if " -> " in rel else rel


def _changed_paths(cwd: Path) -> list[str]:
    r = subprocess.run(
        ["git", "status", "--porcelain", "-uall"],
        cwd=cwd,
        capture_output=True,
        text=True,
        check=False,
    )
    return [p for p in (_parse_status_line(ln) for ln in (r.stdout or "").splitlines()) if p]


def _is_allowed(path: str) -> bool:
    posix = path.replace("\\", "/")
    return any(posix == a.rstrip("/") or posix.startswith(a) for a in DEFAULT_ALLOWED)


def _filter_allowed(paths: list[str]) -> list[str]:
    return [p.replace("\\", "/") for p in paths if _is_allowed(p)]


def _session_path(root: Path, iter_n: int) -> Path:
    d = root / "sessions"
    d.mkdir(parents=True, exist_ok=True)
    return d / f"cycle-{iter_n:04d}.json"


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )


def _progress_line(
    iter_n: int, dual0: float, dual1: float, decision: str, cluster_key: str
) -> str:
    return (
        f"[iter {iter_n}] dual {dual0:.4f}→{dual1:.4f} "
        f"decision={decision} cluster={cluster_key}"
    )


def _should_stop(state: dict[str, Any], cfg: LoopConfig, memory_root: Path) -> bool:
    """Internal stop policy; prefer AgentLoop.should_stop for public callers."""
    if (memory_root / "STOP").is_file():
        return True
    if int(state.get("no_improve") or 0) >= int(cfg.k):
        return True
    if float(state.get("last_dual") or 0.0) >= float(cfg.target_dual):
        return True
    return int(state.get("iter") or 0) >= int(cfg.max_iters)


def _bench_side(report: dict[str, Any]) -> dict[str, Any]:
    return {
        "dual": _dual(report),
        "en": float(report.get("acc_en") or 0.0),
        "zh": float(report.get("acc_zh") or 0.0),
        "fails": _n_fails(report),
    }


def _session_head(
    iter_n: int, decision: str, base_sha: str, commit_sha: str | None
) -> dict[str, Any]:
    status = "gate_pass" if decision == "commit" else "gate_fail"
    return {
        "iter": iter_n,
        "status": status,
        "base_sha": base_sha,
        "commit_sha": commit_sha,
        "rules": [],
        "decision": decision,
    }


def _session_cluster_bits(clusters: list[Cluster]) -> dict[str, Any]:
    layers = clusters[0]["suggest_layers"] if clusters else []
    return {
        "layers": layers,
        "cluster_summary": _cluster_text(clusters),
        "prompt_path": str(_PROMPT_OUT).replace("\\", "/"),
    }


def _session_bench_bits(ctx: dict[str, Any]) -> dict[str, Any]:
    return {
        "bench_before": _bench_side(ctx["before"]),
        "bench_after": _bench_side(ctx["after"]),
        "tdd": {"passed": ctx["pytest_ok"]},
        "lint_ok": ctx["lint_ok"],
        "log_excerpt": (ctx["pi_res"].log or "")[:2000],
    }


def _session_payload(iter_n: int, ctx: dict[str, Any]) -> dict[str, Any]:
    p = _session_head(iter_n, ctx["decision"], ctx["base_sha"], ctx["commit_sha"])
    p.update(_session_cluster_bits(ctx["clusters"]))
    p.update(_session_bench_bits(ctx))
    return p


def _capability_dual(ctx: dict[str, Any]) -> float:
    """Restored dual after revert/noop; observed after only on real commit."""
    if ctx["decision"] == "commit":
        return _dual(ctx["after"])
    return _dual(ctx["before"])


def _pack_gate(
    decision: str, base_sha: str, commit_sha: str | None, lint_ok: bool, pytest_ok: bool
) -> dict[str, Any]:
    return {
        "decision": decision,
        "base_sha": base_sha,
        "commit_sha": commit_sha,
        "lint_ok": lint_ok,
        "pytest_ok": pytest_ok,
    }


def _pack_ctx(
    before: dict[str, Any],
    after: dict[str, Any],
    clusters: list[Cluster],
    pi_res: PiResult,
    gate: dict[str, Any],
) -> dict[str, Any]:
    out = {"before": before, "after": after, "clusters": clusters, "pi_res": pi_res}
    out.update(gate)
    return out


class AgentLoop:
    """Orchestrate one or many self-improve cycles with git gate."""

    def __init__(
        self,
        config: LoopConfig,
        namer_factory: NamerFactory | None = None,
        pi_runner: PiRunner | MockPiRunner | None = None,
        bus: EventBus | None = None,
        *,
        bench_fn: BenchFn | None = None,
        pytest_fn: BoolFn | None = None,
        lint_fn: LintFn | None = None,
        git: GitGate | None = None,
        store: StateStore | None = None,
    ) -> None:
        self._wire_core(config, namer_factory, pi_runner, bus)
        self._wire_deps(bench_fn, pytest_fn, lint_fn, git, store)

    def _wire_core(
        self,
        config: LoopConfig,
        namer_factory: NamerFactory | None,
        pi_runner: PiRunner | MockPiRunner | None,
        bus: EventBus | None,
    ) -> None:
        self.config = _ensure_paths(config)
        self.namer_factory = namer_factory or _default_namer
        self.pi_runner = pi_runner or _default_pi(self.config)
        self.bus = bus or EventBus()
        self.cwd = Path(self.config.cwd)

    def _wire_deps(
        self,
        bench_fn: BenchFn | None,
        pytest_fn: BoolFn | None,
        lint_fn: LintFn | None,
        git: GitGate | None,
        store: StateStore | None,
    ) -> None:
        self.memory = store.root if store else (self.cwd / _MEMORY)
        self.store = store or StateStore(self.memory)
        self.git = git or GitGate(self.cwd)
        self.bench_fn = bench_fn
        self.pytest_fn = pytest_fn
        self.lint_fn = lint_fn

    @staticmethod
    def should_stop(
        state: dict[str, Any], config: LoopConfig, memory: Path
    ) -> bool:
        """True if STOP file, no_improve>=K, dual target, or max iters."""
        return _should_stop(state, config, memory)

    def run(self) -> None:
        """Run cycles until a stop condition holds."""
        while True:
            if self.should_stop(self.store.load(), self.config, self.memory):
                return
            self.run_once()
            if self.should_stop(self.store.load(), self.config, self.memory):
                return

    def run_once(self) -> CycleResult:
        """Execute one full state-machine cycle; return CycleResult."""
        state = self.store.load()
        iter_n = int(state.get("iter") or 0) + 1
        self._publish("cycle_start", {"iter": iter_n})
        ctx = self._cycle_body(iter_n, self.git.snapshot())
        return self._finish(state, iter_n, ctx)

    def _cycle_body(self, iter_n: int, base_sha: str) -> dict[str, Any]:
        before = self._bench()
        clusters = cluster_failures(list(before.get("fails") or []), top_k=5)
        self._write_prompt(iter_n, _dual(before), clusters)
        return self._after_pi(before, clusters, base_sha)

    def _after_pi(
        self, before: dict[str, Any], clusters: list[Cluster], base_sha: str
    ) -> dict[str, Any]:
        pi_res = self._run_pi()
        lint_ok, pytest_ok = self._checks()
        after = self._bench()
        decision, commit_sha = self._gate(
            before, after, lint_ok, pytest_ok, pi_res, base_sha
        )
        gate = _pack_gate(decision, base_sha, commit_sha, lint_ok, pytest_ok)
        return _pack_ctx(before, after, clusters, pi_res, gate)

    def _finish(
        self, state: dict[str, Any], iter_n: int, ctx: dict[str, Any]
    ) -> CycleResult:
        cluster_key = ctx["clusters"][0]["key"] if ctx["clusters"] else None
        self._persist(state, iter_n, ctx, cluster_key)
        session = self._write_session(iter_n, ctx)
        result = CycleResult(
            iter_n, _dual(ctx["before"]), _dual(ctx["after"]), ctx["decision"], session
        )
        self._publish("cycle_end", {"iter": iter_n, "decision": ctx["decision"]})
        return result

    def _persist(
        self,
        state: dict[str, Any],
        iter_n: int,
        ctx: dict[str, Any],
        cluster_key: str | None,
    ) -> None:
        state = self._update_state(state, iter_n, ctx, cluster_key)
        self.store.save(state)
        self._log_progress(iter_n, ctx, cluster_key)

    def _log_progress(
        self, iter_n: int, ctx: dict[str, Any], cluster_key: str | None
    ) -> None:
        line = _progress_line(
            iter_n, _dual(ctx["before"]), _dual(ctx["after"]), ctx["decision"], str(cluster_key)
        )
        self.store.append_progress(line)

    def _bench(self) -> dict[str, Any]:
        if self.bench_fn is not None:
            return self.bench_fn()
        namer = self.namer_factory()
        return run_benchmark(namer, self.config.data_path, limit=self.config.bench_limit)

    def _resolve_tmpl(self) -> Path:
        tmpl = self.cwd / _PROMPT_TMPL
        return tmpl if tmpl.is_file() else Path(_PROMPT_TMPL)

    def _write_prompt(self, iter_n: int, dual: float, clusters: list[Cluster]) -> Path:
        mapping = _prompt_map(
            iter_n,
            dual,
            _cluster_text(clusters),
            _SKILL,
            _CONSTRAINTS,
            _bench_cmd(self.config),
        )
        return _render_prompt(self._resolve_tmpl(), self.cwd / _PROMPT_OUT, mapping)

    def _run_pi(self) -> PiResult:
        return self.pi_runner.run(self.cwd / _PROMPT_OUT, self.cwd, timeout_s=3600)

    def _checks(self) -> tuple[bool, bool]:
        issues = self.lint_fn() if self.lint_fn else lint_tree(self.cwd / _NAMEPREDICT)
        pytest_ok = self.pytest_fn() if self.pytest_fn else _run_pytest(self.cwd)
        return len(issues) == 0, pytest_ok

    def _gate(
        self,
        before: dict[str, Any],
        after: dict[str, Any],
        lint_ok: bool,
        pytest_ok: bool,
        pi_res: PiResult,
        base_sha: str,
    ) -> tuple[str, str | None]:
        if not pi_res.success or not lint_ok or not pytest_ok:
            return self._do_revert(base_sha)
        if not _is_improve(before, after):
            return self._do_revert(base_sha)
        return self._commit_if_any(base_sha, before, after)

    def _do_revert(self, base_sha: str) -> tuple[str, str | None]:
        self.git.revert(base_sha)
        return "revert", None

    def _commit_if_any(
        self, base_sha: str, before: dict[str, Any], after: dict[str, Any]
    ) -> tuple[str, str | None]:
        paths = _filter_allowed(_changed_paths(self.cwd))
        if not paths:
            return "noop", None
        msg = (
            f"agent cycle dual {_dual(before):.4f}->{_dual(after):.4f} "
            f"fails {_n_fails(before)}->{_n_fails(after)}"
        )
        return "commit", self.git.commit(msg, paths)

    def _update_state(
        self,
        state: dict[str, Any],
        iter_n: int,
        ctx: dict[str, Any],
        cluster_key: str | None,
    ) -> dict[str, Any]:
        dual = _capability_dual(ctx)
        self._fill_state_common(state, iter_n, dual, cluster_key)
        self._fill_state_gate(state, dual, ctx["decision"])
        return state

    def _fill_state_common(
        self, state: dict[str, Any], iter_n: int, dual: float, cluster_key: str | None
    ) -> None:
        state["iter"] = iter_n
        state["last_dual"] = dual
        state["last_cluster"] = cluster_key
        state["target"] = self.config.target_dual
        state["K"] = self.config.k

    def _fill_state_gate(
        self, state: dict[str, Any], dual: float, decision: str
    ) -> None:
        if decision == "commit":
            state["no_improve"] = 0
            state["best_dual"] = max(float(state.get("best_dual") or 0.0), dual)
            return
        state["no_improve"] = int(state.get("no_improve") or 0) + 1

    def _write_session(self, iter_n: int, ctx: dict[str, Any]) -> Path:
        path = _session_path(self.memory, iter_n)
        _write_json(path, _session_payload(iter_n, ctx))
        return path

    def _publish(self, event_type: str, payload: dict[str, Any]) -> None:
        self.bus.publish(event_type, payload)
