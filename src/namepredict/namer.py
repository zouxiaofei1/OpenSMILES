from __future__ import annotations

import time

from namepredict.cache.common_names import CommonNameCache
from namepredict.layer0.preprocessor import preprocess
from namepredict.layer0.salt import dissociate_salt
from namepredict.layer1.analyzer import analyze
from namepredict.layer2.claimable_block import ClaimedBlock, SideSlot
from namepredict.layer2.parent_ownership import finalize_parent_ownership
from namepredict.layer2.parent_selector import iter_parent_candidates, select_parent
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


def _ok_result(numbered: dict, *, depth: int, t0: float) -> NameResult | None:
    result = assemble(numbered, time_ms=_elapsed_ms(t0))
    if not result.success or not result.en:
        return None
    result.meta = {**(result.meta or {}), **_chain_meta(numbered), "depth": depth, "coverage_complete": True}
    return result


def _assemble_candidate(parent, subst, *, depth: int, t0: float) -> NameResult | None:
    return _ok_result(number(parent, subst), depth=depth, t0=t0)


def try_candidate(
    info: dict, parent: dict, *, depth: int = 0, t0: float | None = None,
) -> NameResult | None:
    """Finalize ownership, extract, require complete ledger, then number/assemble."""
    t0 = t0 if t0 is not None else time.perf_counter()
    mol = info["mol"]
    parent = finalize_parent_ownership(parent, mol)
    subst = extract_substituents(info, parent)
    if not _ledger_complete(mol, parent["owned_atoms"], subst):
        return None
    return _assemble_candidate(parent, subst, depth=depth, t0=t0)


def _run_candidates(info: dict, *, depth: int, t0: float) -> NameResult:
    attempts: list[dict] = []
    cands = list(iter_parent_candidates(info)) or [select_parent(info)]
    for parent in cands:
        hit = try_candidate(info, parent, depth=depth, t0=t0)
        if hit is not None and hit.success:
            return hit
        attempts.append({"kind": parent.get("kind"), "reason": "incomplete_or_unnamed"})
    return _fail(_elapsed_ms(t0), "no_complete_candidate", attempts=attempts)


def _name_mol(
    mol,
    *,
    depth: int = 0,
    cache: CommonNameCache | None = None,
    t0: float | None = None,
) -> NameResult:
    """Run L1–L5 from mol with coverage-gated candidate retry."""
    t0 = t0 if t0 is not None else time.perf_counter()
    if mol is None:
        return _fail(_elapsed_ms(t0), "parse")
    organic, salt = dissociate_salt(mol)
    result = _run_candidates(analyze(organic), depth=depth, t0=t0)
    if salt and result.success:
        result.meta = {**(result.meta or {}), "salt": salt}
    return result


def _pipeline(smiles: str, t0: float) -> NameResult:
    mol = preprocess(smiles)
    if mol is None:
        return _fail(_elapsed_ms(t0), "parse")
    return _name_mol(mol, depth=0, t0=t0)


def _cache_put(cache: CommonNameCache, smiles: str, result: NameResult) -> None:
    try:
        cache.put(smiles, result)
    except ValueError:
        pass


def _name_uncached(smiles: str, t0: float) -> NameResult:
    return _pipeline(smiles, t0)


class SMILESNNamer:
    def __init__(self, cache: CommonNameCache | None = None) -> None:
        self.cache = cache if cache is not None else CommonNameCache()

    def name(self, smiles: str) -> NameResult:
        t0 = time.perf_counter()
        hit = self.cache.get(smiles)
        if hit is not None:
            return hit
        result = _name_uncached(smiles, t0)
        if result.success:
            _cache_put(self.cache, smiles, result)
        return result
