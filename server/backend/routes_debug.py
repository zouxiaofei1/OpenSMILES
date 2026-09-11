"""Debug API: SMILES -> layer-by-layer pipeline trace."""

from __future__ import annotations

import os
import subprocess
import sys
import threading
import time
from pathlib import Path
from typing import Any

from fastapi import APIRouter
from fastapi.responses import StreamingResponse
from pydantic import BaseModel, Field
from rdkit import Chem

from namepredict.layer0.preprocessor import preprocess
from namepredict.layer1.analyzer import analyze
from namepredict.layer2.parent_selector import select_parent
from namepredict.layer2.parent_ownership import finalize_parent_ownership
from namepredict.layer3.substituent_extractor import extract_substituents
from namepredict.layer3.coverage import build_coverage_ledger
from namepredict.layer3.substituent_namer import SubstituentName
from namepredict.layer3.claimable_block import ClaimedBlock, SideSlot
from namepredict.layer4.numbering import number
from namepredict.layer5.assembler import assemble

router = APIRouter(prefix="/api/v1", tags=["debug"])

_ROOT = Path(__file__).resolve().parents[2]


class DebugBody(BaseModel):
    smiles: str = Field(..., min_length=1)


def _mol_info(mol) -> dict:
    """Serializable subset of Mol, incl. raw atom/bond tables for the debug UI.

    The `atoms` / `bonds` / `smiles` fields mirror the full RDKit graph so the
    frontend can expand the raw molecule without the Mol object ever crossing
    the API boundary (JSON has no object references).
    """
    from collections import Counter
    elem = Counter()
    atoms = []
    for a in mol.GetAtoms():
        elem[a.GetSymbol()] += 1
        atoms.append({
            "idx": a.GetIdx(),
            "symbol": a.GetSymbol(),
            "charge": a.GetFormalCharge(),
            "h": a.GetTotalNumHs(),
            "aromatic": a.GetIsAromatic(),
            "degree": a.GetDegree(),
            "nbrs": [n.GetIdx() for n in a.GetNeighbors()],
        })
    bonds = []
    for b in mol.GetBonds():
        bonds.append({
            "a": b.GetBeginAtomIdx(),
            "b": b.GetEndAtomIdx(),
            "order": int(b.GetBondTypeAsDouble()),
            "aromatic": b.GetIsAromatic(),
        })
    return {
        "num_atoms": mol.GetNumAtoms(),
        "num_bonds": mol.GetNumBonds(),
        "composition": " ".join(f"{e}={n}" for e, n in sorted(elem.items())),
        "smiles": Chem.MolToSmiles(mol),
        "atoms": atoms,
        "bonds": bonds,
    }


def _info_serializable(info: dict) -> dict:
    """Clone info dict, drop mol and make JSON-safe."""
    out = {}
    for k, v in info.items():
        if k == "mol":
            out["mol"] = _mol_info(v)
        elif isinstance(v, set):
            out[k] = sorted(v)
        elif isinstance(v, list):
            out[k] = [_clean_item(x) for x in v]
        elif isinstance(v, dict):
            out[k] = _clean_item(v)
        else:
            out[k] = v
    return out


def _clean_item(x):
    if isinstance(x, dict):
        return {kk: list(vv) if isinstance(vv, (set, tuple)) else vv for kk, vv in x.items()}
    if isinstance(x, (set, tuple)):
        return list(x)
    return x


def _parent_serializable(parent: dict) -> dict:
    """Make parent dict JSON-safe."""
    out = {}
    for k, v in parent.items():
        if k in ("mol",):
            continue
        if isinstance(v, set):
            out[k] = sorted(v)
        elif isinstance(v, list) and k in ("chain", "owned_atoms"):
            out[k] = list(v)
        elif isinstance(v, list):
            out[k] = [_clean_item(x) for x in v]
        elif isinstance(v, dict):
            out[k] = _clean_item(v)
        else:
            out[k] = v
    return out


def _subst_serializable(subst: list[dict]) -> list[dict]:
    out = []
    for s in subst:
        d = {}
        for k, v in s.items():
            if isinstance(v, set):
                d[k] = sorted(v)
            elif isinstance(v, frozenset):
                d[k] = sorted(v)
            elif isinstance(v, list):
                d[k] = list(v)
            elif isinstance(v, dict):
                d[k] = _clean_item(v)
            else:
                d[k] = v
        out.append(d)
    return out


@router.post("/debug")
def debug_smiles(body: DebugBody) -> dict[str, Any]:
    t0 = time.perf_counter()

    # ── L0: preprocess ──
    l0_start = time.perf_counter()
    mol = preprocess(body.smiles)
    l0_ms = (time.perf_counter() - l0_start) * 1000
    if mol is None:
        return {"success": False, "error": "SMILES parse failed", "layer": "L0"}

    l0_out = {
        "smiles": body.smiles,
        ** _mol_info(mol),
        "time_ms": round(l0_ms, 2),
    }

    # ── L1: analyze ──
    l1_start = time.perf_counter()
    info = analyze(mol)
    l1_ms = (time.perf_counter() - l1_start) * 1000
    l1_out = {
        **_info_serializable(info),
        "time_ms": round(l1_ms, 2),
    }

    # ── L2: parent candidates / selection ──
    l2_start = time.perf_counter()
    candidates_raw = select_parent(info, all_candidates=True)
    selected = select_parent(info)
    candidates = []
    for c in candidates_raw:
        c_ser = _parent_serializable(c)
        c_ser["is_selected"] = (c.get("kind") == selected.get("kind"))
        candidates.append(c_ser)
    selected_ser = _parent_serializable(selected)
    l2_ms = (time.perf_counter() - l2_start) * 1000
    l2_out = {
        "n_candidates": len(candidates),
        "candidates": candidates,
        "selected": selected_ser,
        "time_ms": round(l2_ms, 2),
    }

    # ── L3: substituents ──
    l3_start = time.perf_counter()
    parent = finalize_parent_ownership(selected, mol)
    owned = parent.get("owned_atoms", set())
    subst = extract_substituents(info, parent)
    # coverage
    names_for_ledger = []
    for s in subst:
        atoms = frozenset(s.get("atoms") or [])
        if atoms:
            attach = s.get("attach_idx")
            names_for_ledger.append(
                SubstituentName(
                    claim=ClaimedBlock(
                        slot=SideSlot.OTHER,
                        attach_parent=int(attach) if attach is not None else -1,
                        root=min(atoms),
                        atoms=atoms,
                    ),
                    en=s.get("en") or "x",
                    zh=s.get("zh") or "x",
                    requires_parentheses=bool(s.get("paren")),
                    backend=s.get("backend") or "extract",
                )
            )
    ledger = build_coverage_ledger(mol, owned_atoms=owned, names=names_for_ledger)
    l3_ms = (time.perf_counter() - l3_start) * 1000
    l3_out = {
        "n_substituents": len(subst),
        "substituents": _subst_serializable(subst),
        "coverage": {
            "complete": ledger.complete,
            "gap": sorted(ledger.gap) if ledger.gap else [],
            "overlap": sorted(ledger.overlap) if ledger.overlap else [],
        },
        "owned_atoms": sorted(owned),
        "time_ms": round(l3_ms, 2),
    }

    # ── L4: numbering ──
    l4_start = time.perf_counter()
    numbered = None
    l4_error = None
    try:
        # filter subs to those on chain
        chain_set = set(parent.get("chain") or [])
        subs_numbered = []
        for s in subst:
            attach = s.get("attach_idx")
            if attach in chain_set:
                subs_numbered.append(s)
        numbered = number(parent, subs_numbered)
    except Exception as e:
        l4_error = str(e)
    l4_ms = (time.perf_counter() - l4_start) * 1000
    if numbered is not None:
        l4_out = {
            "parent_kind": numbered.get("parent", {}).get("kind"),
            "chain": list(numbered.get("parent", {}).get("chain") or []),
            "substituents": _subst_serializable(numbered.get("substituents", [])),
            "time_ms": round(l4_ms, 2),
        }
        # extract key locants
        p = numbered.get("parent", {})
        locants = {}
        for k in p:
            if k.endswith("_locant") or k.endswith("_locants"):
                v = p[k]
                locants[k] = list(v) if isinstance(v, (list, tuple)) else v
        l4_out["locants"] = locants
    else:
        l4_out = {"error": l4_error, "time_ms": round(l4_ms, 2)}

    # ── L5: assembly ──
    l5_start = time.perf_counter()
    l5_out = {}
    if numbered is not None:
        try:
            result = assemble(numbered, time_ms=0)
            l5_out = {
                "en": result.en,
                "zh": result.zh,
                "success": result.success,
                "source": result.source,
                "meta": {k: v for k, v in (result.meta or {}).items()
                         if k not in ("parent_chain",)},
            }
        except Exception as e:
            l5_out = {"error": str(e)}
    else:
        l5_out = {"error": "L4 failed, no numbered dict"}
    l5_ms = (time.perf_counter() - l5_start) * 1000
    l5_out["time_ms"] = round(l5_ms, 2)

    total_ms = (time.perf_counter() - t0) * 1000

    return {
        "success": True,
        "total_time_ms": round(total_ms, 2),
        "layers": {
            "L0_preprocess": l0_out,
            "L1_analyze": l1_out,
            "L2_parent": l2_out,
            "L3_substituents": l3_out,
            "L4_numbering": l4_out,
            "L5_assemble": l5_out,
        },
    }


@router.post("/debug-print")
def debug_print(body: DebugBody) -> dict[str, Any]:
    """Run server/backend/debug.py on the SMILES and return its captured stdout.

    Executes in a fresh subprocess so every call reflects the latest src/
    code, and any print() executed anywhere in the naming path ends up in
    `stdout` — exactly what the "打印调试" panel in debug.html shows.
    """
    script = _ROOT / "server" / "backend" / "debug.py"
    env = {**os.environ, "PYTHONIOENCODING": "utf-8"}
    try:
        proc = subprocess.run(
            [sys.executable, str(script), body.smiles],
            capture_output=True,
            encoding="utf-8",
            errors="replace",
            timeout=60,
            env=env,
            cwd=str(_ROOT),
        )
    except subprocess.TimeoutExpired:
        return {"success": False, "error": "debug.py timed out after 60s"}
    return {
        "success": proc.returncode == 0,
        "returncode": proc.returncode,
        "stdout": proc.stdout,
        "stderr": proc.stderr,
        "smiles": body.smiles,
    }


# 命名路径中 print() 产物需能在产生时立刻看到，/debug-print 是跑完才一次性
# 返回，故这里用 Popen + 流式响应：stdout 逐行即时下发(实时滚出)，stderr 由
# 守护线程并行收集(避免管道写满阻塞子进程)，结束时把 stderr/退出码补在尾部。
_DEBUG_STREAM_TIMEOUT = 60  # 与旧 /debug-print 的 timeout 保持一致
_ERR_KEEP_LINES = 400


def _debug_print_stream(smiles: str):
    """Yield debug.py stdout lines live, then stderr tail / exit code at the end."""
    script = _ROOT / "server" / "backend" / "debug.py"
    # PYTHONUNBUFFERED + -u: 子进程写管道默认块缓冲，会把 src 命名管线的 print
    # 攒到缓冲满才出现；去掉缓冲才能逐行实时到达。
    env = {**os.environ, "PYTHONIOENCODING": "utf-8", "PYTHONUNBUFFERED": "1"}
    try:
        proc = subprocess.Popen(
            [sys.executable, "-u", str(script), smiles],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            encoding="utf-8",
            errors="replace",
            env=env,
            cwd=str(_ROOT),
        )
    except Exception as exc:
        yield f"[error] 启动 debug.py 失败: {exc}\n"
        return

    err_lines: list[str] = []

    def _drain_err() -> None:
        if proc.stderr is None:
            return
        try:
            for line in proc.stderr:
                err_lines.append(line)
                if len(err_lines) > _ERR_KEEP_LINES + 100:
                    del err_lines[: len(err_lines) - _ERR_KEEP_LINES]
        except Exception:
            pass

    drain_thread = threading.Thread(target=_drain_err, daemon=True)
    drain_thread.start()
    guard = threading.Timer(_DEBUG_STREAM_TIMEOUT, proc.kill)  # 兜底, 避免永久挂起
    guard.start()

    try:
        if proc.stdout is not None:
            for line in proc.stdout:
                yield line
    except GeneratorExit:
        # 前端取消/断开: 立即终止子进程, 丢弃尾部
        guard.cancel()
        try:
            proc.kill()
        except Exception:
            pass
        try:
            proc.wait(timeout=10)
        except Exception:
            pass
        drain_thread.join(timeout=3)
        raise
    finally:
        guard.cancel()
        try:
            proc.wait(timeout=15)
        except Exception:
            try:
                proc.kill()
            except Exception:
                pass
        drain_thread.join(timeout=3)

    # 正常跑完后把退出码 / stderr 追到末尾, 与 /debug-print 的展示对齐
    if proc.returncode != 0 or err_lines:
        yield "\n"
        if proc.returncode != 0:
            yield f"--- 退出码 {proc.returncode} ---\n"
        if err_lines:
            yield "--- stderr ---\n"
            for line in err_lines:
                yield line


@router.post("/debug-print-stream")
def debug_print_stream(body: DebugBody) -> StreamingResponse:
    """Stream server/backend/debug.py output live (text/plain), line by line."""
    return StreamingResponse(
        _debug_print_stream(body.smiles),
        media_type="text/plain; charset=utf-8",
    )
