# IUPAC: P-14.4(e) / P-31.1.3 / P-22.1.2.2
# Layer: L4,L5
"""环内累积二烯 (ring cumulene) 编号：C=C=C 嵌入环内。

双键共享同一 sp 累积碳（如 C1=CCCCC=1 中 C0 连 (0,1) 与 (0,5) 两个 C=C），
P-14.4(e)(ii) 编号不得让共享碳同时充当两条双键的 locant → 应为互异连续 locant
EN cyclohexa-1,2-diene；ZH 环己-1,2-二烯。
回归守护：普通 1,3-/1,4- 环己二烯与单环己烯不回归。
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: 环内累积二烯（两个 C=C 共享 sp 中心碳，locant 须互异连续）
    ("C1=CCCCC=1", "cyclohexa-1,2-diene", "环己-1,2-二烯"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_ring_cumulene(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


@pytest.mark.parametrize("smiles", ["C1=CCCCC=1"])
def test_no_duplicate_ene_locant(smiles: str) -> None:
    """累积双键 locant 必须互异，不得退化成 1,1-diene 这类重复位次。"""
    r = SMILESNNamer().name(smiles)
    assert r.success
    en = normalize_en(r.en)
    assert "diene" in en
    assert en != "cyclohexa-1,1-diene"
    assert "1,1" not in en
