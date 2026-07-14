# IUPAC: P-22.2.2
# Layer: L2,L3,L4,L5
"""Retained / Hantzsch–Widman saturated monocyclic hetero parents (P-22.2.2 / P-22.2.3 / P-25).

Unsubstituted (and optional mono-methyl / mono-halo on ring C) monoheterocycles:
aziridine, oxirane, oxolane, oxane, pyrrolidine, piperidine, morpholine,
piperazine, 1,3-dioxolane, 1,4-dioxane; optional thiolane.

Ring heteroatoms (O/S/N) belong to the heterocyclic parent (Hantzsch–Widman /
retained), NOT open-chain ether / dialkyl sulfide / sec-amine parents. Open-
chain ether/sulfide/sec|tert amine generators must require the bridging
heteroatom is not in a ring (atom.IsInRing()).
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: unsubstituted mono-O / mono-N
    ("C1CCOC1", "oxolane", "氧杂环戊烷"),
    ("C1CN1", "aziridine", "氮杂环丙烷"),
    ("C1CO1", "oxirane", "环氧乙烷"),
    ("C1CCOCC1", "oxane", "氧杂环己烷"),
    ("C1CCNC1", "pyrrolidine", "吡咯烷"),
    ("C1CCNCC1", "piperidine", "哌啶"),
    # positive: di-hetero retained
    ("C1COCO1", "1,3-dioxolane", "1,3-二氧戊环"),
    ("C1COCCO1", "1,4-dioxane", "1,4-二氧六环"),
    ("C1COCCN1", "morpholine", "吗啉"),
    ("C1CNCCN1", "piperazine", "哌嗪"),
    # positive: optional thiolane
    ("C1CCSC1", "thiolane", "硫杂环戊烷"),
    # positive: mono-methyl (hetero fixed = 1 → lowest set)
    ("CC1CCCO1", "2-methyloxolane", "2-甲基氧杂环戊烷"),
    ("CC1CCNCC1", "4-methylpiperidine", "4-甲基哌啶"),
    # negative: carbocycle / arene / open chain must not break
    ("C1CCCCC1", "cyclohexane", "环己烷"),
    ("c1ccncc1", "pyridine", "吡啶"),
    ("CCCCC", "pentane", "戊烷"),
    ("c1ccoc1", "furan", "呋喃"),
    # negative: open-chain ether / sec-amine / sulfide must keep open-chain names
    ("CCOCC", "diethyl ether", "二乙基醚"),
    ("CCNCC", "N-ethylethanamine", "N-乙基乙胺"),
    ("CCSCC", "diethyl sulfide", "二乙硫醚"),
    ("COC", "dimethyl ether", "二甲基醚"),
    ("CNC", "N-methylmethanamine", "N-甲基甲胺"),
    ("CSC", "dimethyl sulfide", "二甲硫醚"),
]


# Ring hetero must not be named as open-chain ether / sec-amine / sulfide.
# ("smiles", "expected_en", forbidden substrings in en)
NOT_OPEN_CHAIN = [
    ("C1CCOC1", "oxolane", ("ether", "dimethyl", "butane")),
    ("C1CO1", "oxirane", ("ether", "dimethyl")),
    ("C1CCOCC1", "oxane", ("ether", "dimethyl")),
    ("C1CN1", "aziridine", ("methanamine", "ethanamine", "N-methyl")),
    ("C1CCNC1", "pyrrolidine", ("methanamine", "ethanamine", "N-methyl")),
    ("C1CCNCC1", "piperidine", ("methanamine", "ethanamine", "N-methyl")),
    ("C1COCCN1", "morpholine", ("methanamine", "ether", "N-methyl")),
    ("C1CCSC1", "thiolane", ("sulfide", "dimethyl")),
    ("CC1CCCO1", "2-methyloxolane", ("ether",)),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_sat_hetero(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


@pytest.mark.parametrize("smiles,en,forbidden", NOT_OPEN_CHAIN)
def test_ring_hetero_not_open_chain(
    smiles: str, en: str, forbidden: tuple[str, ...],
) -> None:
    """Endocyclic O/S/N must not become dialkyl ether / sulfide / sec-amine."""
    r = SMILESNNamer().name(smiles)
    assert r.success
    got = normalize_en(r.en)
    assert got == normalize_en(en)
    for frag in forbidden:
        assert frag not in got


def test_oxolane_not_butane() -> None:
    """Current failure: tetrahydrofuran core must not fall to butane."""
    r = SMILESNNamer().name("C1CCOC1")
    assert r.success
    assert normalize_en(r.en) == "oxolane"
    assert "butane" not in normalize_en(r.en)
    assert "ether" not in normalize_en(r.en)


def test_n_methylpiperidine_not_parent() -> None:
    """N-alkyl on ring N is out of scope (negative: must not claim piperidine)."""
    r = SMILESNNamer().name("CN1CCCCC1")
    assert r.success
    assert "piperidine" not in normalize_en(r.en)


def test_piperidine_carboxylic_not_oxolane() -> None:
    """Principal acid FG: must not invent sat-hetero parent name."""
    r = SMILESNNamer().name("O=C(O)C1CCCCN1")
    assert r.success
    en = normalize_en(r.en)
    assert "oxolane" not in en
    assert "oxane" not in en


def test_open_chain_ether_not_sat_hetero() -> None:
    """Open diethyl ether must remain ether (not invent a sat-hetero parent)."""
    r = SMILESNNamer().name("CCOCC")
    assert r.success
    en = normalize_en(r.en)
    assert en == "diethyl ether"
    assert "oxolane" not in en
    assert "oxane" not in en


def test_open_chain_sec_amine_not_sat_hetero() -> None:
    """Open sec-amine must remain N-alkylethanamine (not aziridine/pyrrolidine)."""
    r = SMILESNNamer().name("CCNCC")
    assert r.success
    en = normalize_en(r.en)
    assert en == "n-ethylethanamine"
    assert "aziridine" not in en
    assert "pyrrolidine" not in en
    assert "piperidine" not in en


def test_open_chain_sulfide_not_sat_hetero() -> None:
    """Open dialkyl sulfide must remain sulfide (not thiolane)."""
    r = SMILESNNamer().name("CCSCC")
    assert r.success
    en = normalize_en(r.en)
    assert en == "diethyl sulfide"
    assert "thiolane" not in en
