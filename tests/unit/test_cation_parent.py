# 单核母体阳离子的 L1 登记与 L2/L5 命名（P-73）。
"""
test_mono_cation.py: 非环阳离子在 L1 被登记为 cation FG，L2 收敛为单原子母体，
其余原子全作取代基；自由价双键/三键的臂出 -ylidene/-ylidyne。
"""
from __future__ import annotations

import pytest

from namepredict.layer1.analyzer import analyze
from namepredict.layer1.functional_group_inventory import FunctionalGroupClass as FG, inventory_from_info
from namepredict.namer import SMILESNNamer
from namepredict.tools.re import normalize_en, normalize_zh
from rdkit import Chem

# ==========================================================================
# IUPAC: P-73.1.1 / P-41 表 4.1 类 6
# Layer: L1,L2,L3,L4,L5
#
# 单核母体阳离子的保留名/系统名（P-21 表 2.1 去 ne）作母体名，
# 其余原子全作取代基；母体只有 1 个原子，故取代基位次恒省。
# ==========================================================================
mono_cation__CASES = [
    # ── 母体只有阳离子本身 ──
    ("[NH4+]", "azanium", "铵"),
    ("C[NH3+]", "methylazanium", "甲基铵"),
    ("C[NH+](C)C", "trimethylazanium", "三甲基铵"),
    # ── 臂为碳链/官能团取代基 ──
    ("[NH3+]CCc1ccc(O)cc1", "2-(4-hydroxyphenyl)ethylazanium", "2-(4-羟基苯基)乙基铵"),
    ("CC(=O)C[NH3+]", "2-oxopropylazanium", "2-氧代丙基铵"),
    ("C[N+](C)(C)CCP(=O)(O)O", "trimethyl(2-phosphonoethyl)azanium", "三甲基(2-膦酸乙基)铵"),
    ("[NH3+][C@@H](Cc1c[nH]c2ccccc12)C(=O)O", "[(1S)-1-carboxy-2-(1H-indol-3-yl)ethyl]azanium",
     "[(1S)-1-羧基-2-(1H-吲哚-3-基)乙基]铵"),
    # ── 自由价非单键的臂：双键 -ylidene，三键 -ylidyne（P-31.2.3） ──
    ("C[NH+]=C", "methyl(methylidene)azanium", "甲基(亚甲基)铵"),
    ("CC(=[NH2+])C", "propan-2-ylideneazanium", "亚丙-2-基铵"),
    ("C#[N+]C", "methyl(methylidyne)azanium", "甲基(次甲基)铵"),
    ("CC#[N+]C", "ethylidyne(methyl)azanium", "次乙基(甲基)铵"),
]


@pytest.mark.parametrize("smiles,en,zh", mono_cation__CASES)
def test_mono_cation(smiles: str, en: str, zh: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# IUPAC: P-73.1.1
# Layer: L1
#
# 阳离子中心须非环，且 radius=1 内无负形式电荷原子：
# 硝基/N-氧化物/异氰等"阳离子寄居在别的基团里"的结构不登记 cation。
# ==========================================================================
def test_cation_only_for_unclaimed_non_ring_centers() -> None:
    claimed = inventory_from_info(analyze(Chem.MolFromSmiles("C[NH+]=C")))
    assert len(claimed.occurrences(FG.CATION)) == 1
    assert not claimed.occurrences(FG.AMINE)

    for smiles in ("[O-][N+](=O)c1ccccc1", "[C-]#[N+]C1CCCCC1", "C[n+]1ccccc1"):
        inv = inventory_from_info(analyze(Chem.MolFromSmiles(smiles)))
        assert not inv.occurrences(FG.CATION)


# ==========================================================================
# IUPAC: P-41 表 4.1 类 4 > 类 6
# Layer: L1,L2
#
# 分子内存在负形式电荷原子时阳离子让位：酸根/阴离子作母体，阳离子退为前缀。
# ==========================================================================
def test_anion_outranks_cation() -> None:
    inv = inventory_from_info(analyze(Chem.MolFromSmiles("[NH3+]CC(=O)[O-]")))
    assert inv.has_anion
    assert inv.occurrences(FG.CATION)


# ==========================================================================
# IUPAC: P-73.1.1
# Layer: L2
#
# 反例：环内 N+/O+ 仍走环阳离子后缀，不归本类的非环母体阳离子。
# ==========================================================================
anti_cation__CASES = [
    ("[O-][N+](=O)c1ccccc1", "nitrobenzene", "硝基苯"),
    ("[C-]#[N+]C1CCCCC1", "isocyanocyclohexane", "异氰基环己烷"),
    ("C[n+]1ccccc1", "1-methylpyridin-1-ium", "1-甲基吡啶-1-鎓"),
]


@pytest.mark.parametrize("smiles,en,zh", anti_cation__CASES)
def test_ring_or_inner_balanced_cation_unchanged(smiles: str, en: str, zh: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)
