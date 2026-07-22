"""Coverage-gated candidate retry: first incomplete candidate yields to next."""
from __future__ import annotations

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer, try_candidate
from namepredict.layer1.analyzer import analyze
from rdkit import Chem


def test_methylbutylbenzene_prefers_benzene_parent():
    r = SMILESNNamer().name("c1ccc(cc1)CC(C)CC")
    assert r.success
    assert normalize_en(r.en) == normalize_en("(2-methylbutyl)benzene")
    assert normalize_zh(r.zh) == normalize_zh("(2-甲基丁基)苯")
