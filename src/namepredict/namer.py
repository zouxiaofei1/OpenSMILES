from __future__ import annotations

import time

from namepredict.cache.common_names import CommonNameCache
from namepredict.layer0.preprocessor import preprocess
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


def _pipeline(smiles: str, t0: float) -> NameResult:
    mol = preprocess(smiles)
    if mol is None:
        return _fail(_elapsed_ms(t0), "parse")
    info = analyze(mol)
    parent = select_parent(info)
    subst = extract_substituents(info, parent)
    numbered = number(parent, subst)
    return assemble(numbered, time_ms=_elapsed_ms(t0))


def _cache_put(cache: CommonNameCache, smiles: str, result: NameResult) -> None:
    try:
        cache.put(smiles, result)
    except ValueError:
        pass


class SMILESNNamer:
    def __init__(self, cache: CommonNameCache | None = None) -> None:
        self.cache = cache if cache is not None else CommonNameCache()

    def name(self, smiles: str) -> NameResult:
        t0 = time.perf_counter()
        hit = self.cache.get(smiles)
        if hit is not None:
            return hit
        result = _pipeline(smiles, t0)
        if result.success:
            _cache_put(self.cache, smiles, result)
        return result
