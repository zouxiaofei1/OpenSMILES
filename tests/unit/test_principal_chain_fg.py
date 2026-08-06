# IUPAC: P-44 / P-65.6 / P-66
# Layer: L2,L3,L4,L5
#
# 验证 rule-driven principal 管线对无环单官能团 ester/amide/aldehyde/nitrile
# 也产出信息完整的 chain parent dict（kind + 单数 *_c_idx + 不饱和 + 酯烷氧基），
# 与 builder 路径 _open_chain_expression 经 _dedupe_parents 按 (kind, chain)
# 合并后命名行为不变。
from __future__ import annotations

from namepredict.constants import normalize_en, normalize_zh
from namepredict.layer0.preprocessor import preprocess
from namepredict.layer1.analyzer import analyze
from namepredict.layer2.kind_registry import pack_parent_stem
from namepredict.layer2.parent_skeleton import SkeletonTopology
from namepredict.layer2.principal_expression import express_chain_principal
from namepredict.layer2.principal_parent import (
    rule_driven_parent_candidates,
    select_principal_parent_skeletons,
)
from namepredict.namer import try_candidate

# 正例：无环单官能团（饱和 + C=C / C≡C）
CASES = [
    # esters
    ("CC(=O)OCC", "ethyl acetate", "乙酸乙酯"),
    ("CCCC(=O)OC", "methyl butanoate", "丁酸甲酯"),
    ("C=CC(=O)OC", "methyl prop-2-enoate", "丙-2-烯酸甲酯"),
    ("CC=CC(=O)OC", "methyl but-2-enoate", "丁-2-烯酸甲酯"),
    ("C#CCC(=O)OC", "methyl but-3-ynoate", "丁-3-炔酸甲酯"),
    # amides
    ("CC(=O)N", "acetamide", "乙酰胺"),
    ("CCCC(=O)N", "butanamide", "丁酰胺"),
    ("C=CC(=O)N", "acrylamide", "丙烯酰胺"),
    ("CC=CC(=O)N", "but-2-enamide", "丁-2-烯酰胺"),
    # aldehydes
    ("CCCCC=O", "pentanal", "戊醛"),
    ("CC(C)C=O", "2-methylpropanal", "2-甲基丙醛"),
    ("CC=CC=O", "but-2-enal", "丁-2-烯醛"),
    ("C=CCC=O", "but-3-enal", "丁-3-烯醛"),
    # nitriles
    ("CCCC#N", "butanenitrile", "丁腈"),
    ("CCC#N", "propanenitrile", "丙腈"),
    ("C=CC#N", "prop-2-enenitrile", "丙-2-烯腈"),
    ("CC=CC#N", "but-2-enenitrile", "丁-2-烯腈"),
]

# 负例：有环 / 二酯 / 高优先级官能团，必须走既有路径不被误伤
NEGATIVE = [
    ("CC(=O)O", "acetic acid", "乙酸"),                        # 羧酸 > 酯
    ("COC(=O)c1ccccc1", "methyl benzoate", "苯甲酸甲酯"),       # 苯环酯 → benzoate
    ("COC(=O)CCC(=O)OC", "dimethyl butanedioate", "丁二酸二甲酯"),  # 二酯
    ("O=C(O)C1CCCCC1", "cyclohexanecarboxylic acid", "环己烷甲酸"),  # 环酸
    ("N#Cc1ccccc1", "benzonitrile", "苯甲腈"),                 # 环腈
    ("C1CCCCC1", "cyclohexane", "环己烷"),                     # 纯烃环
]


def _principal_name(smiles: str):
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


def test_chain_fg_cases():
    for smiles, en, zh in CASES:
        r = _principal_name(smiles)
        assert r is not None, f"principal 管线未能产出候选: {smiles}"
        assert r.success, f"装配失败 {smiles}: {r.meta}"
        assert normalize_en(r.en) == normalize_en(en), f"EN {smiles}: got {r.en!r} want {en!r}"
        assert normalize_zh(r.zh) == normalize_zh(zh), f"ZH {smiles}: got {r.zh!r} want {zh!r}"


def test_chain_fg_negative_untouched():
    for smiles, en, zh in NEGATIVE:
        r = _principal_name(smiles)
        assert r is not None, f"principal 管线未产出: {smiles}"
        assert normalize_en(r.en) == normalize_en(en), f"EN {smiles}: got {r.en!r} want {en!r}"
        assert normalize_zh(r.zh) == normalize_zh(zh), f"ZH {smiles}: got {r.zh!r} want {zh!r}"


# ── 表达层字段完整性：express_chain_principal 直接产出完整 parent dict ──

def _chain_parent(smiles: str, kind: str):
    mol = preprocess(smiles)
    if mol is None:
        return None
    info = analyze(mol)
    selection = select_principal_parent_skeletons(info)
    if selection.principal is None or selection.skeletons is None:
        return None
    for skeleton in selection.skeletons.candidates:
        if skeleton.topology is not SkeletonTopology.ACYCLIC:
            continue
        parent = express_chain_principal(info, selection.principal, skeleton)
        if parent is not None and parent["kind"] == kind:
            return parent
    return None


def test_chain_ester_expression_fields():
    p = _chain_parent("CC=CC(=O)OC", "ester")
    assert p is not None
    assert p["ester_c_idx"] is not None          # 单数（L4 _orient_ester 依赖）
    assert p["ester_c_idxs"] == [p["ester_c_idx"]]
    assert p["n_carbons"] == 4
    assert isinstance(p["double_bond"], tuple) and len(p["double_bond"]) == 2
    assert p["alkoxy_n"] == 1                    # 甲酯


def test_chain_ester_saturated_no_unsat():
    p = _chain_parent("CCCC(=O)OC", "ester")
    assert p is not None
    assert p.get("double_bond") is None and p.get("triple_bond") is None
    assert p.get("alkoxy_n") == 1


def test_chain_amide_aldehyde_nitrile_fields():
    for smiles, kind, key in [
        ("CC=CC(=O)N", "amide", "amide_c_idx"),
        ("CC=CC=O", "aldehyde", "aldehyde_c_idx"),
        ("CC=CC#N", "nitrile", "nitrile_c_idx"),
    ]:
        p = _chain_parent(smiles, kind)
        assert p is not None, f"{smiles} 无 {kind} 表达候选"
        assert p[key] is not None, f"{smiles} 缺单数 {key}"
        assert p["n_carbons"] == 4
        assert isinstance(p["double_bond"], tuple), f"{smiles} 缺 double_bond"


def test_chain_nitrile_yne_field():
    p = _chain_parent("CC#CC#N", "nitrile")
    assert p is not None
    assert isinstance(p["triple_bond"], tuple) and len(p["triple_bond"]) == 2
