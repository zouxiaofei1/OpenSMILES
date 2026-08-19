# IUPAC: P-25.3.2.4
# Layer: L2
"""稠环拆解 fused_system: 稠环系统 → 保留母体组分树（P-25.3.2.4(a)-(f)）。"""
from __future__ import annotations

from namepredict.layer0.preprocessor import preprocess
from namepredict.layer1.analyzer import analyze
from namepredict.layer2.fused_system import FusedNode, decompose_fused_system, fused_info_for_mol, fused_node_dict


def _info(smiles: str):
    mol = preprocess(smiles)
    assert mol is not None
    return analyze(mol)


def _root(smiles: str) -> FusedNode:
    """首个环系的拆解根节点。"""
    info = _info(smiles)
    assert info["ring_systems"], smiles
    return decompose_fused_system(info, info["ring_systems"][0])


def test_chrysene_decomposes_to_phenanthrene_plus_benzene():
    """未注册稠环 chrysene → 根 phenanthrene + 附加 benzene(P-25.3.2.4(j) 稠合碳位次低)。"""
    root = _root("c1cc2c3ccccc3ccc2c2ccccc12")
    assert root.scaffold_id == "phenanthrene"
    assert root.ring_indices == frozenset({0, 2, 3})
    assert len(root.atom_ids) == 14
    assert root.fusion_shared == ()
    assert len(root.attached) == 1
    child = root.attached[0]
    assert child.scaffold_id == "benzene"
    assert child.ring_indices == frozenset({1})
    assert child.fusion_shared == (frozenset({3, 8}),)
    assert child.attached == ()


def test_tetracene_decomposes_to_anthracene_plus_benzene():
    """并四苯(未注册) → 根 anthracene + 附加 benzene。"""
    root = _root("c1ccc2cc3cc4ccccc4cc3cc2c1")
    assert root.scaffold_id == "anthracene"
    assert root.ring_indices == frozenset({0, 1, 2})
    assert len(root.attached) == 1
    assert root.attached[0].scaffold_id == "benzene"


def test_registered_phenanthrene_pyrene_are_single_node():
    """注册保留名后, phenanthrene/pyrene 整环匹配 → 单节点。"""
    assert _root("c1ccc2c(c1)ccc1ccccc12").scaffold_id == "phenanthrene"
    assert _root("c1cc2ccc3cccc4ccc(c1)c2c34").scaffold_id == "pyrene"


def test_registered_fused_systems_are_single_node():
    for smiles in (
        "c1ccc2ccccc2c1",     # naphthalene
        "c1ccc2cc3ccccc3cc2c1",  # anthracene
        "c1ccc2[nH]c3ccccc3c2c1",  # carbazole
        "c1ccc2nc3ccccc3cc2c1",  # acridine
        "c1ccc2Sc3ccccc3Nc2c1",  # phenothiazine
        "c1ccc2ncccc2c1",     # quinoline
        "c1ccc2[nH]ccc2c1",   # indole
        "c1ccc2c(c1)ccc1ccccc12",  # phenanthrene
        "c1cc2ccc3cccc4ccc(c1)c2c34",  # pyrene
    ):
        root = _root(smiles)
        assert root.attached == (), smiles
        assert root.scaffold_id in ("naphthalene", "anthracene", "carbazole",
                                    "acridine", "phenothiazine", "quinoline", "indole",
                                    "phenanthrene", "pyrene")


def test_mono_ring_is_single_node():
    for smiles, sid in (("c1ccccc1", "benzene"), ("c1ccncc1", "pyridine")):
        root = _root(smiles)
        assert root.scaffold_id == sid
        assert root.attached == ()


def test_biphenyl_yields_two_independent_systems():
    info = _info("c1ccc(-c2ccccc2)cc1")
    nodes = fused_info_for_mol(info)
    assert len(nodes) == 2
    assert all(n.scaffold_id == "benzene" and n.attached == () for n in nodes)


def test_saturated_ring_has_no_retained_candidate():
    info = _info("C1CCCCC1")
    assert fused_info_for_mol(info) == []


def test_fused_node_dict_is_json_safe():
    import json
    root = _root("c1cc2c3ccccc3ccc2c2ccccc12")
    d = fused_node_dict(root)
    json.dumps(d)  # 无异常即 JSON 安全
    assert d["scaffold"] == "phenanthrene"
    assert d["attached_components"][0]["scaffold"] == "benzene"
    assert d["attached_components"][0]["fusion_shared"] == [[3, 8]]
