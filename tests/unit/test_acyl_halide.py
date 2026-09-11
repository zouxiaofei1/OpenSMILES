# IUPAC: P-65.5
# Layer: L1,L2,L5
"""Acyl halides: alkanoyl/benzoyl/enoyl F·Cl·Br·I (P-65.5).

主 FG 命名后缀随卤素（hal_z）变：-oyl halide / 酰氟氯溴碘（保留名 acetyl/
benzoyl 亦同）。覆盖开链 C2–C4、2-甲基丙酰、苯甲酰与烯酰。负数：酸、
溴代烷、醛不得判成酰卤。
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # ── 正例：C1/C2 保留名与开链后缀随卤素变 ──
    ("CC(=O)F", "acetyl fluoride", "乙酰氟"),
    ("CC(=O)Cl", "acetyl chloride", "乙酰氯"),
    ("CC(=O)Br", "acetyl bromide", "乙酰溴"),
    ("CC(=O)I", "acetyl iodide", "乙酰碘"),
    ("CCC(=O)Cl", "propanoyl chloride", "丙酰氯"),
    ("CCC(=O)F", "propanoyl fluoride", "丙酰氟"),
    ("CCC(=O)I", "propanoyl iodide", "丙酰碘"),
    # 2-甲基丙酰（基准金标：Br 必须 bromide 而非 chloride）
    ("CC(C)C(=O)Br", "2-methylpropanoyl bromide", "2-甲基丙酰溴"),
    ("CC(C)C(=O)I", "2-methylpropanoyl iodide", "2-甲基丙酰碘"),
    # ── 苯甲酰卤（环外酰卤主 FG 同样带卤素）──
    ("O=C(F)c1ccccc1", "benzoyl fluoride", "苯甲酰氟"),
    ("O=C(Br)c1ccccc1", "benzoyl bromide", "苯甲酰溴"),
    ("O=C(I)c1ccccc1", "benzoyl iodide", "苯甲酰碘"),
    # ── 烯酰卤 ──
    ("O=C(Cl)C=C", "prop-2-enoyl chloride", "丙-2-烯酰氯"),
    ("O=C(F)C=C", "prop-2-enoyl fluoride", "丙-2-烯酰氟"),
    ("O=C(Br)C=C", "prop-2-enoyl bromide", "丙-2-烯酰溴"),
    # ── 炔酰卤（P-14.3.4：三碳炔酰保留炔位次，与烯酰卤同形）──
    ("O=C(Cl)C#C", "prop-2-ynoyl chloride", "丙-2-炔酰氯"),
    ("O=C(Br)C#C", "prop-2-ynoyl bromide", "丙-2-炔酰溴"),
    ("O=C(Cl)CC#C", "but-3-ynoyl chloride", "丁-3-炔酰氯"),
    # ── 负数：必须不判成酰卤 ──
    ("CC(=O)O", "acetic acid", "乙酸"),
    ("CCBr", "bromoethane", "溴乙烷"),
    ("CC=O", "acetaldehyde", "乙醛"),
    # F 在 α 碳而非羰基：仍为醛（不误判为酰卤）
    ("FCC=O", "fluoroacetaldehyde", "氟乙醛"),
    # 羰基碳带 F 但有酯烷氧基侧：是氟甲酸酯，不是酰氟
    ("O=C(F)OC", "methyl fluoroformate", None),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_acyl_halide(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
