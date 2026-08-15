# IUPAC: P-22 / P-25 template-based retained ring identification
# Layer: L2
"""SMILES 模板子图同构识别（原 scaffold/retained_templates.py 移植）。

五元组 _TOPOLOGY 无法区分位置异构体（azulene/naphthalene、
isoindole/indole 等字段全等），模板方案修复这些误配。
"""
from __future__ import annotations

from namepredict.layer0.preprocessor import preprocess
from namepredict.layer1.analyzer import analyze
from namepredict.layer2.parent_skeleton import ParentSkeleton, SkeletonTopology
from namepredict.layer2.ring_scaffold import match_retained, resolve_ring_scaffold


def _resolve(smiles: str) -> str | None:
    mol = preprocess(smiles)
    assert mol is not None
    info = analyze(mol)
    system = info["ring_systems"][0]
    skeleton = ParentSkeleton(
        SkeletonTopology.RING_SYSTEM, tuple(system["atom_ids"]), frozenset(),
    )
    identity = resolve_ring_scaffold(info, skeleton)
    return identity.id if identity else None


def _template_id(smiles: str) -> str | None:
    mol = preprocess(smiles)
    assert mol is not None
    info = analyze(mol)
    return match_retained(info, info["ring_systems"][0]["atom_ids"])


def test_registered_templates_resolve():
    assert _resolve("c1ccccc1") == "benzene"
    assert _resolve("c1ccncc1") == "pyridine"
    assert _resolve("c1ccc2ccccc2c1") == "naphthalene"
    assert _resolve("c1ccc2[nH]ccc2c1") == "indole"


def test_azulene_no_longer_misidentified_as_naphthalene():
    # 5+7 稠合全碳环：五元组与 naphthalene 字段全等曾误配为 naphthalene。
    assert _resolve("C1=CC2=CC=CC=CC2=CC1") == "carbocycle"


def test_isoindole_no_longer_misidentified_as_indole():
    # 苯并[c]吡咯：五元组与 indole 字段全等曾误配为 indole；含 N 回落 None。
    assert _resolve("c1ccc2c(c1)c[nH]c2") is None


def test_substituted_ring_still_resolves():
    assert _resolve("Oc1ccccc1") == "benzene"
    assert _resolve("Oc1cccc2ccccc12") == "naphthalene"


def test_all_templates_resolve_to_own_scaffold():
    # 唯一事实来源：模板全部派生 spec，命中即解析为自身。
    assert _resolve("c1ccoc1") == "furan"
    assert _resolve("c1ccc2cc3ccccc3cc2c1") == "anthracene"


def test_positional_isomers_hit_own_template_only():
    # 元素标注的子图同构区分位置异构体，各自单命中。
    cases = {
        "quinoline": "c1ccc2ncccc2c1",
        "isoquinoline": "c1nccc2ccccc21",
        "pyridazine": "c1ccnnc1",
        "pyrimidine": "c1cncnc1",
        "pyrazine": "c1cnccn1",
        "imidazole": "c1cnc[nH]1",
        "pyrazole": "c1ccn[nH]1",
        "benzimidazole": "c1ccc2[nH]cnc2c1",
        "indazole": "c1ccc2cn[nH]c2c1",
        "quinazoline": "c1ccc2ncncc2c1",
        "quinoxaline": "c1ccc2nccnc2c1",
    }
    for sid, smi in cases.items():
        assert _template_id(smi) == sid
