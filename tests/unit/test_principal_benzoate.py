# IUPAC: P-65.6.1 / P-22.1.3
# Layer: L2
#
# 验证 principal 规则管线对「苯环上直接连酯基」产出 benzoate 保留母体，
# 不依赖经典生产者通道（_try_arene_other_fg / _benzoate_parent）。
from namepredict.tools.re import normalize_en, normalize_zh
from namepredict.layer0.preprocessor import preprocess
from namepredict.layer1.analyzer import analyze
from namepredict.layer2.kind_registry import pack_parent_stem
from namepredict.layer2.principal_parent import rule_driven_parent_candidates
from namepredict.namer import try_candidate

# 正例：alkyl benzoate（苯甲酸酯），侧链 C1–C16
CASES = [
    ("COC(=O)c1ccccc1", "methyl benzoate", "苯甲酸甲酯"),
    ("CCOC(=O)c1ccccc1", "ethyl benzoate", "苯甲酸乙酯"),
    ("CCCCCCOC(=O)c1ccccc1", "hexyl benzoate", "苯甲酸己酯"),
    ("CCCCCCCCOC(=O)c1ccccc1", "octyl benzoate", "苯甲酸辛酯"),
    ("CCCCCCCCCCCCCCCCOC(=O)c1ccccc1", "hexadecyl benzoate", "苯甲酸十六酯"),
]

# 负例：近邻但不应被 benzoate 规则误伤
NEGATIVE = [
    ("OC(=O)c1ccccc1", "benzoic acid", "苯甲酸"),  # 苯甲酸走 ACID 保留名
    ("CC(=O)OC", "methyl acetate", "乙酸甲酯"),    # 开链酯走开链母体
]


def _principal_name(smiles: str):
    """Only the rule-driven principal pipeline + L3-L5, no legacy producers."""
    mol = preprocess(smiles)
    if mol is None:
        return None
    info = analyze(mol)
    for parent in rule_driven_parent_candidates(info):
        packed = pack_parent_stem(parent, info["mol"])
        hit = try_candidate(info, packed, depth=0)
        if hit is not None and hit.success:
            return hit
    return None


def test_principal_benzoate_cases():
    for smiles, en, zh in CASES:
        r = _principal_name(smiles)
        assert r is not None, f"principal 管线未能产出候选: {smiles}"
        assert r.success, f"装配失败 {smiles}: {r.meta}"
        assert normalize_en(r.en) == normalize_en(en), f"EN {smiles}: got {r.en!r} want {en!r}"
        assert normalize_zh(r.zh) == normalize_zh(zh), f"ZH {smiles}: got {r.zh!r} want {zh!r}"


def test_principal_benzoate_negative_untouched():
    for smiles, en, zh in NEGATIVE:
        r = _principal_name(smiles)
        assert r is not None, f"principal 管线未产出: {smiles}"
        assert normalize_en(r.en) == normalize_en(en), f"EN {smiles}: got {r.en!r} want {en!r}"
        assert normalize_zh(r.zh) == normalize_zh(zh), f"ZH {smiles}: got {r.zh!r} want {zh!r}"
