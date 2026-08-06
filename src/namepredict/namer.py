from __future__ import annotations

import time

from namepredict.cache.common_names import CommonNameCache
from namepredict.layer0.preprocessor import preprocess
from namepredict.layer0.salt import dissociate_salt
from namepredict.layer1.analyzer import analyze
from namepredict.layer3.claimable_block import ClaimedBlock, SideSlot
from namepredict.layer2.parent_ownership import finalize_parent_ownership
from namepredict.layer2.parent_selector import iter_parent_candidates, select_parent
from namepredict.layer2.parent_candidate import principal_phases
from namepredict.layer3.coverage import build_coverage_ledger
from namepredict.layer3.substituent_extractor import extract_substituents
from namepredict.layer3.substituent_namer import SubstituentName
from namepredict.layer4.numbering import number
from namepredict.layer5.assembler import assemble
from namepredict.types import NameResult


def _fail(time_ms: float = 0.0, reason: str = "parse", **meta) -> NameResult:
    return NameResult(
        en="", zh="", success=False, source="iupac",
        time_ms=time_ms, meta={"reason": reason, **meta},
    )


def _elapsed_ms(t0: float) -> float:
    return (time.perf_counter() - t0) * 1000.0


def _chain_meta(numbered: dict) -> dict:
    parent = numbered.get("parent") or {}
    return {"parent_chain": list(parent.get("chain") or []), "parent_kind": parent.get("kind")}


def _claim_from_sub(s: dict, atoms: frozenset[int]) -> ClaimedBlock:
    attach = s.get("attach_idx")
    return ClaimedBlock(
        slot=SideSlot.OTHER,
        attach_parent=int(attach) if attach is not None else -1,
        root=min(atoms),
        atoms=atoms,
    )


def _one_name_from_sub(s: dict) -> SubstituentName | None:
    atoms = frozenset(s.get("atoms") or [])
    if not atoms:
        return None
    return SubstituentName(
        claim=_claim_from_sub(s, atoms),
        en=s.get("en") or "x",
        zh=s.get("zh") or "x",
        requires_parentheses=bool(s.get("paren")),
        backend=s.get("backend") or "extract",
    )


def _names_from_subs(subs: list[dict]) -> list[SubstituentName]:
    return [n for s in subs if (n := _one_name_from_sub(s)) is not None]


def _ledger_complete(mol, owned, subst: list[dict]) -> bool:
    names = _names_from_subs(subst)
    return build_coverage_ledger(mol, owned_atoms=owned, names=names).complete


def _ok_result(numbered: dict, *, depth: int, t0: float, name_mode: str = "general") -> NameResult | None:
    numbered["name_mode"] = name_mode
    result = assemble(numbered, time_ms=_elapsed_ms(t0))
    if not result.success or not result.en:
        return None
    result.meta = {**(result.meta or {}), **_chain_meta(numbered), "depth": depth, "coverage_complete": True}
    return result


def _chain_set(parent: dict) -> set[int]:
    return set(parent.get("chain") or [])


def _remap_attach(parent: dict, s: dict) -> dict:
    """Ensure attach_idx is on parent chain for L4 orient (ring FG attach)."""
    chain = _chain_set(parent)
    attach = s.get("attach_idx")
    if attach in chain:
        return s
    for key in ("ring_attach_idx", "amide_c_idx", "ketone_c_idx", "oh_c_idx"):
        alt = parent.get(key)
        if alt in chain:
            return {**s, "attach_idx": alt}
    return s


def _subs_for_numbering(parent: dict, subst: list[dict]) -> list[dict]:
    chain = _chain_set(parent)
    out: list[dict] = []
    for s in subst:
        s2 = _remap_attach(parent, s)
        if s2.get("attach_idx") in chain:
            out.append(s2)
    return out


def _assemble_candidate(parent, subst, *, depth: int, t0: float, name_mode: str = "general") -> NameResult | None:
    try:
        numbered = number(parent, _subs_for_numbering(parent, subst))
    except (ValueError, KeyError, TypeError):
        return None
    return _ok_result(numbered, depth=depth, t0=t0, name_mode=name_mode)


def _prepare_candidate(
    info: dict, parent: dict, *, name_mode: str = "general", cache: CommonNameCache | None = None,
) -> tuple[dict, list[dict], bool]:
    mol = info["mol"]
    parent = finalize_parent_ownership(parent, mol)
    if not parent.get("owned_atoms"):
        return parent, [], False
    if not parent.get("chain") and not info.get("has_ring"):
        return parent, [], False
    subst = extract_substituents(info, parent, name_mode=name_mode, cache=cache)
    complete = _ledger_complete(mol, parent["owned_atoms"], subst)
    return parent, subst, complete


def try_candidate(
    info: dict,
    parent: dict,
    *,
    depth: int = 0,
    t0: float | None = None,
    require_complete: bool = True,
    name_mode: str = "general",
) -> NameResult | None:
    """Finalize ownership, extract, optionally require complete ledger, then assemble."""
    t0 = t0 if t0 is not None else time.perf_counter()
    parent, subst, complete = _prepare_candidate(info, parent, name_mode=name_mode)
    if require_complete and not complete:
        return None
    hit = _assemble_candidate(parent, subst, depth=depth, t0=t0)
    if hit is None:
        return None
    hit.meta = {**(hit.meta or {}), "coverage_complete": complete}
    return hit


def _complete_hit(prepared, *, depth, t0, name_mode):
    for parent, subst, complete in prepared:
        if complete and (hit := _assemble_candidate(parent, subst, depth=depth, t0=t0, name_mode=name_mode)):
            hit.meta = {**(hit.meta or {}), "coverage_complete": True}
            return hit
    return None


def _partial_hit(prepared, *, depth, t0, name_mode, attempts):
    for parent, subst, complete in prepared:
        hit = _assemble_candidate(parent, subst, depth=depth, t0=t0, name_mode=name_mode)
        if hit is not None:
            hit.meta = {**(hit.meta or {}), "coverage_complete": complete, "fallback": "no_coverage_gate", "attempts": attempts}
            return hit
    return None


def _try_phase(prepared, *, depth, t0, name_mode, attempts):
    hit = _complete_hit(prepared, depth=depth, t0=t0, name_mode=name_mode)
    return hit or _partial_hit(
        prepared, depth=depth, t0=t0, name_mode=name_mode, attempts=attempts,
    )


def _candidate_phases(info: dict, depth: int) -> list[list[dict]]:
    if depth > 0:
        return [[select_parent(info)]]
    return  [[select_parent(info)]]


def _run_candidates(
    info: dict, *, depth: int, t0: float, name_mode: str = "general", cache: CommonNameCache | None = None,
) -> NameResult:
    """Try only the P-44.1.1-senior phase; never capability-downgrade."""
    attempts: list[dict] = []
    phase = _candidate_phases(info, depth)[0]
    prepared = [_prepare_candidate(info, cand, name_mode=name_mode, cache=cache) for cand in phase]
    hit = _try_phase(prepared, depth=depth, t0=t0, name_mode=name_mode, attempts=attempts)
    return hit or _fail(_elapsed_ms(t0), "no_assemblable_candidate", attempts=attempts)


def _name_mol(
    mol,
    *,
    depth: int = 0,
    cache: CommonNameCache | None = None,
    t0: float | None = None,
    name_mode: str = "general",
) -> NameResult:
    """Run L1–L5 from mol with coverage-gated candidate retry."""
    t0 = t0 if t0 is not None else time.perf_counter()
    if mol is None:
        return _fail(_elapsed_ms(t0), "parse")
    organic, salt = dissociate_salt(mol)
    result = _run_candidates(analyze(organic), depth=depth, t0=t0, name_mode=name_mode, cache=cache)
    if salt and result.success:
        result.meta = {**(result.meta or {}), "salt": salt}
    return result


def _pipeline(smiles: str, t0: float, *, name_mode: str = "general", cache: CommonNameCache | None = None) -> NameResult:
    mol = preprocess(smiles)
    if mol is None:
        return _fail(_elapsed_ms(t0), "parse")
    return _name_mol(mol, depth=0, t0=t0, name_mode=name_mode, cache=cache)


def _cache_put(cache: CommonNameCache, smiles: str, result: NameResult) -> None:
    try:
        cache.put(smiles, result)
    except ValueError:
        pass


def _name_uncached(smiles: str, t0: float, *, name_mode: str = "general", cache: CommonNameCache | None = None) -> NameResult:
    return _pipeline(smiles, t0, name_mode=name_mode, cache=cache)


class SMILESNNamer:
    def __init__(self, cache: CommonNameCache | None = None, *, name_mode: str = "general") -> None:
        # 容量留足给主分子 + 递归子结构命名（全量去重后约 8.7k 条）
        self.cache = cache if cache is not None else CommonNameCache(max_entries=20000)
        self._name_mode = name_mode

    def name(self, smiles: str) -> NameResult:
        t0 = time.perf_counter()
        hit = self.cache.get(smiles)
        if hit is not None:
            return hit
        result = _name_uncached(smiles, t0, name_mode=self._name_mode, cache=self.cache)
        if result.success:
            _cache_put(self.cache, smiles, result)
        return result
