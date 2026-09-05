# IUPAC: P-14.4(j)
# Layer: L4 numbering
"""P-14.4(j)：镜像/反向等价编号平局用 CIP 立体描述符破局——较低位次赋予 R/M/r（而非 S/P/s），
编号方向不随输入 SMILES 写法漂移。chebi-2118 内消旋 cis-1,7-dimethyl-4-(propan-2-yl)cyclodecane
两条镜像编号均给最低位次集 {1,4,7}，须确定性地让 R 中心取 locant 1 → (1R,7S)。"""

from __future__ import annotations

from rdkit import Chem

from namepredict.layer0.preprocessor import preprocess
from namepredict.layer1.analyzer import analyze
from namepredict.layer2.parent_ownership import finalize_parent_ownership
from namepredict.layer2.parent_selector import select_parent
from namepredict.layer3.substituent_extractor import extract_substituents
from namepredict.layer4.numbering_engine import orient_numbering
from namepredict.namer import SMILESNNamer

# (smiles, expected_en) —— 同一内消旋环的两种等价写法都必须指向同一方向
_CASES = [
    (
        "CC(C)C1CC[C@H](C)CCC[C@H](C)CC1",  # canonical：R 中心位于环书写起点附近
        "(1R,7S)-1,7-dimethyl-4-propan-2-ylcyclodecane",
    ),
    (
        "[C@H]1(C)CCC[C@@H](C)CCC(C(C)C)CC1",  # 重排原子序：S 中心位于分子开头
        "(1R,7S)-1,7-dimethyl-4-propan-2-ylcyclodecane",
    ),
]


def test_p144j_chebi_2118_meso_ring_rs_locant_stable():
    """同一分子不同输入写法，P-14.4(j) 保证都给出 R 中心在 locant 1 的 (1R,7S)。"""
    namer = SMILESNNamer()
    outputs = [namer.name(s).en for s, _ in _CASES]
    assert outputs[0] == outputs[1]  # 不再随原子序漂移
    for s, expected in _CASES:
        assert namer.name(s).en == expected


def test_p144j_orient_numbering_direction_invariant():
    """numbering 引擎层面：无论初始 chain 从哪个镜像方向给入，定向后 locant 1 恒为 CIP-R 中心。"""
    smi = "CC(C)C1CC[C@H](C)CCC[C@H](C)CC1"
    mol = preprocess(smi)
    info = analyze(mol)
    parent = select_parent(info)
    parent = finalize_parent_ownership(parent, mol)
    subs = extract_substituents(info, parent)
    base = parent["chain"]
    assert base[0] == 3  # 初始链从异丙基碳起，验证下面的收窄确实改写了起点
    # 环上两个手性中心的 CIP（与输入无关的绝对构型）
    from rdkit.Chem import rdCIPLabeler

    m2 = Chem.MolFromSmiles(Chem.MolToSmiles(mol, isomericSmiles=True))
    Chem.AssignStereochemistry(m2, force=True, cleanIt=True)
    rdCIPLabeler.AssignCIPLabels(m2)
    cip = {a.GetIdx(): a.GetProp("_CIPCode") for a in m2.GetAtoms() if a.HasProp("_CIPCode")}
    # 多种初始编号方向
    start_from_s = [11, 13, 14, 3, 4, 5, 6, 8, 9, 10]  # S 中心做 locant 1 的镜像方向
    orderings = {
        "base": base,
        "rev": list(reversed(base)),
        "start_from_s": start_from_s,
    }
    for name, chain0 in orderings.items():
        oriented = orient_numbering(
            {"mol": mol, "kind": "alkane", "scaffold_id": "carbocycle", "chain": chain0},
            subs,
        )
        locant1 = oriented[0]
        assert cip[locant1] == "R", f"{name}: locant 1 ({locant1}) 应为 R，实际 {cip[locant1]}"
