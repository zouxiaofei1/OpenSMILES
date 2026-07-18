from __future__ import annotations

import time

from namepredict.cache.common_names import CommonNameCache
from namepredict.layer0.preprocessor import preprocess
from namepredict.layer0.salt import dissociate_salt
from namepredict.layer1.analyzer import analyze
from namepredict.layer2.parent_selector import select_parent
from namepredict.layer3.substituent_extractor import extract_substituents
from namepredict.layer4.numbering import number
from namepredict.layer5.assembler import assemble
from namepredict.types import NameResult


def _fail(time_ms: float = 0.0, reason: str = "parse") -> NameResult:
    return NameResult(
        en="", zh="", success=False, source="iupac",
        time_ms=time_ms, meta={"reason": reason},
    )


def _elapsed_ms(t0: float) -> float:
    return (time.perf_counter() - t0) * 1000.0


def _run_layers(mol) -> dict:
    info = analyze(mol)
    parent = select_parent(info)
    subst = extract_substituents(info, parent)
    return number(parent, subst)


def _with_salt_meta(numbered: dict, salt: dict) -> dict:
    return {**numbered, "salt": salt} if salt else numbered


def _chain_meta(numbered: dict) -> dict:
    parent = numbered.get("parent") or {}
    return {"parent_chain": list(parent.get("chain") or []), "parent_kind": parent.get("kind")}


def _name_mol(mol, *, depth: int = 0, cache: CommonNameCache | None = None) -> NameResult:
    """Run L1–L5 from mol. depth>0 disables cache write (caller must not put)."""
    t0 = time.perf_counter()
    if mol is None:
        return _fail(_elapsed_ms(t0), "parse")
    organic, salt = dissociate_salt(mol)
    numbered = _with_salt_meta(_run_layers(organic), salt)
    result = assemble(numbered, time_ms=_elapsed_ms(t0))
    result.meta = {**(result.meta or {}), **_chain_meta(numbered), "depth": depth}
    return result


def _pipeline(smiles: str, t0: float) -> NameResult:
    mol = preprocess(smiles)
    if mol is None:
        return _fail(_elapsed_ms(t0), "parse")
    return _name_mol(mol, depth=0)


def _cache_put(cache: CommonNameCache, smiles: str, result: NameResult) -> None:
    try:
        cache.put(smiles, result)
    except ValueError:
        pass


def _name_uncached(smiles: str, t0: float) -> NameResult:
    result = _pipeline(smiles, t0)
    return result


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
