# IUPAC: P-44.1.2 / P-44.3 / P-44.4
# Layer: L2
#
# 验证 rule_driven_parent_candidates 在【无主官能团】时独立产出纯烃候选，
# 不依赖经典生产者通道（ring/unsat/benzene/alkane fallback）。
# 通过 _principal_name 直接走 principal 管线 + L3-L5 装配，避免经典通道假绿。
from rdkit import Chem

from namepredict.constants import normalize_en, normalize_zh
from namepredict.layer0.preprocessor import preprocess
from namepredict.layer1.analyzer import analyze
from namepredict.layer2.principal_parent import rule_driven_parent_candidates
from namepredict.namer import try_candidate

# 正例：无主官能团的纯烃（开链饱和/烯/炔/多烯、单环、苯、稠环、带烷基侧链）
CASES = [
    ("CC", "ethane", "乙烷"),
    ("CCCCC", "pentane", "戊烷"),
    ("C=C", "ethene", "乙烯"),
    ("C#CC", "propyne", "丙炔"),
    ("C=CC=C", "buta-1,3-diene", "丁-1,3-二烯"),
    ("C1CCCCC1", "cyclohexane", "环己烷"),
    ("C1CC=CCC1", "cyclohexene", "环己烯"),
    ("c1ccccc1", "benzene", "苯"),
    ("c1ccc2ccccc2c1", "naphthalene", "萘"),
    ("CC1CCCCC1", "methylcyclohexane", "甲基环己烷"),
    ("c1ccccc1CC", "ethylbenzene", "乙基苯"),
]

# 负例：含主官能团，必须继续走主官能团管线，不被纯烃逻辑误伤
NEGATIVE = [
    ("CC(C)O", "propan-2-ol", "丙-2-醇"),
    ("OC1CCCCC1", "cyclohexanol", "环己醇"),
]


def _principal_name(smiles: str):
    """Only the rule-driven principal pipeline + L3-L5, no legacy producers."""
    from namepredict.layer2.kind_registry import pack_parent_stem

    mol = preprocess(smiles)
    if mol is None:
        return None
    info = analyze(mol)
    for parent in rule_driven_parent_candidates(info):
        # 模拟完整管线的词干填充（iter_parent_candidates 内部会 pack）。
        packed = pack_parent_stem(parent, info["mol"])
        hit = try_candidate(info, packed, depth=0)
        if hit is not None and hit.success:
            return hit
    return None


def test_principal_hydrocarbon_cases():
    for smiles, en, zh in CASES:
        r = _principal_name(smiles)
        assert r is not None, f"principal 管线未能产出候选: {smiles}"
        assert r.success, f"装配失败 {smiles}: {r.meta}"
        assert normalize_en(r.en) == normalize_en(en), f"EN {smiles}: got {r.en!r} want {en!r}"
        assert normalize_zh(r.zh) == normalize_zh(zh), f"ZH {smiles}: got {r.zh!r} want {zh!r}"


def test_principal_hydrocarbon_negative_fg_untouched():
    for smiles, en, zh in NEGATIVE:
        r = _principal_name(smiles)
        assert r is not None, f"principal 管线未产出: {smiles}"
        assert normalize_en(r.en) == normalize_en(en), f"EN {smiles}: got {r.en!r} want {en!r}"
        assert normalize_zh(r.zh) == normalize_zh(zh), f"ZH {smiles}: got {r.zh!r} want {zh!r}"


