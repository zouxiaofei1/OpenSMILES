# 合并自 6 个原测试文件（按主题分组，内容与断言未改动）。
"""
test_kind_registry.py: ParentKind registry: single source for scoring meta + L5 parent stems.
test_retained_registry.py: Retained registry matches ring_systems for benzene/pyridine/naphthalene/indole.
test_leaf_registry_extend.py: Registry leaves: n-alkyl C2–C4, methylthio, cyano; nested benzyl skeleton.
test_principal_registry.py: Principal characteristic-group metadata is the sole priority authority.
test_scaffold_registry_single_source.py: Single source of truth: ScaffoldSpec drives kind_registry stems + retained ids.
test_retained_names.py: 保留名（retained names）层：取代基保留名 anilino 与官能团衍生物保留前缀 carbamoyl/carbamoylamino/
"""
from __future__ import annotations

import pytest

from namepredict.layer0.preprocessor import preprocess
from namepredict.layer1.analyzer import analyze
from namepredict.layer1.functional_group_inventory import FunctionalGroupClass as FG, FunctionalGroupInventory, FunctionalGroupOccurrence
from namepredict.layer2 import kind_registry as kr
from namepredict.layer2.parent_selector import select_parent
from namepredict.layer2.principal import PRINCIPAL_REGISTRY, PrincipalExpression, PrincipalFeatureSpec, PrincipalPriority, select_principal_group
from namepredict.layer2.ring_scaffold import all_specs, get_spec, match_retained
from namepredict.namer import SMILESNNamer
from namepredict.tools.re import normalize_en, normalize_zh

# ==========================================================================
# 合并自 test_kind_registry.py
# IUPAC: P-44 / architecture
# Layer: L2,L5
#
# ParentKind registry: single source for scoring meta + L5 parent stems.
#
# 正交化后 kind 只表达 FG 类别 / 保留 scaffold；环系维度由 scaffold_id 承载，
# 命名 kind（cycloalcohol/benzoic/…）由 L5 typed_kinds 决定，不注册在此。
# ==========================================================================
kind_registry___RETAINED_KINDS = {
    "benzene": ("benzene", "苯"),
    "pyridine": ("pyridine", "吡啶"),
    "naphthalene": ("naphthalene", "萘"),
    "indole": ("1H-indole", "吲哚"),
}

@pytest.mark.parametrize("kind,names", kind_registry___RETAINED_KINDS.items())
def test_retained_kind_stems(kind: str, names: tuple[str, str]) -> None:
    m = kr.get(kind)
    assert m is not None
    assert kr.parent_names(kind) == names
    assert m.retained is True


def test_ether_not_registered() -> None:
    # ether 非主官能团（compat 0），不作为 scaffold 词干注册。
    assert kr.get("ether") is None


def test_unknown_kind_negative() -> None:
    assert kr.get("no_such") is None
    assert kr.parent_names("no_such") is None


# 正交化契约：环 + 主 FG → FG 类别 kind + scaffold_id；无主 FG 保留 scaffold 结构 kind。
kind_registry___ORTHO_CASES = [
    ("Oc1ccccc1", "alcohol", "benzene"),
    ("Nc1ccccc1", "amine", "benzene"),
    ("O=C(O)c1ccccc1", "acid", "benzene"),
    ("n1ccccc1", "pyridine", "pyridine"),
    ("c1ccc2ccccc2c1", "naphthalene", "naphthalene"),
    ("c1ccc2[nH]ccc2c1", "indole", "indole"),
]


@pytest.mark.parametrize("smiles,kind,scaffold", kind_registry___ORTHO_CASES)
def test_select_parent_orthogonalized_kind(smiles: str, kind: str, scaffold: str) -> None:
    mol = preprocess(smiles)
    assert mol is not None
    parent = select_parent(analyze(mol))[0]
    assert parent["kind"] == kind
    assert parent.get("scaffold_id") == scaffold


def test_select_parent_preloads_registry_stem() -> None:
    """有 spec stem 的无 FG 环 scaffold：select_parent 预填 kind_registry 词干。"""
    for smiles, kind in (("n1ccccc1", "pyridine"), ("c1ccc2ccccc2c1", "naphthalene")):
        mol = preprocess(smiles)
        assert mol is not None
        parent = select_parent(analyze(mol))[0]
        assert parent["kind"] == kind
        assert (parent.get("stem_en"), parent.get("stem_zh")) == kr.parent_names(kind)


def test_select_parent_preserves_existing_stem(monkeypatch: pytest.MonkeyPatch) -> None:
    from namepredict.layer2 import candidates

    parent = {"kind": "phenol", "stem_en": "custom", "stem_zh": "自定义",
              "principal_group_count": 1, "mol": object()}
    monkeypatch.setattr(candidates, "_collect_candidates", lambda _: [parent])
    selected = select_parent({"mol": parent["mol"]})[0]
    assert selected["stem_en"] == "custom"
    assert selected["stem_zh"] == "自定义"


def kind_registry___p(kind: str, **kw) -> dict:
    return {"kind": kind, **kw}


kind_registry__REGRESSION = [
    ("c1ccccc1", "benzene", "苯"),
    ("c1ccc2ccccc2c1", "naphthalene", "萘"),
    ("CCO", "ethanol", "乙醇"),
]


@pytest.mark.parametrize("smiles,en,zh", kind_registry__REGRESSION)
def test_behavior_regression(smiles: str, en: str, zh: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_retained_registry.py
# IUPAC: P-22 / P-25
# Layer: L2
#
# Retained registry matches ring_systems for benzene/pyridine/naphthalene/indole.
# ==========================================================================
def retained_registry___ids(smiles: str) -> list[str]:
    mol = preprocess(smiles)
    assert mol is not None
    info = analyze(mol)
    return [sid for s in info["ring_systems"] if (sid := match_retained(info, s["atom_ids"]))]


def test_benzene_registry():
    assert "benzene" in retained_registry___ids("c1ccccc1")


def test_pyridine_registry():
    assert "pyridine" in retained_registry___ids("c1ccncc1")


def test_naphthalene_registry():
    assert "naphthalene" in retained_registry___ids("c1ccc2ccccc2c1")


def test_indole_registry():
    assert "indole" in retained_registry___ids("c1ccc2[nH]ccc2c1")


def test_biphenyl_two_benzene():
    ids = retained_registry___ids("c1ccc(-c2ccccc2)cc1")
    assert ids.count("benzene") == 2


def test_open_chain_no_match():
    assert retained_registry___ids("CCCC") == []


def test_match_systems_has_entry():
    mol = preprocess("c1ccc2ccccc2c1")
    info = analyze(mol)
    system = info["ring_systems"][0]
    sid = match_retained(info, system["atom_ids"])
    assert sid == "naphthalene"
    assert get_spec(sid).stem_en == "naphthalene"
    assert system["n_rings"] == 2


# ==========================================================================
# 合并自 test_leaf_registry_extend.py
# IUPAC: P-29.3
# Layer: L2
#
# Registry leaves: n-alkyl C2–C4, methylthio, cyano; nested benzyl skeleton.
#
# Ph arm on chain alcohol/amine/acid; Ar-CN must not become formonitrile parent.
# Prior halo/MeO/nested Ph and bare arene alkyl parents stay correct.
# ==========================================================================
leaf_registry_extend__CASES = [
    # n-alkyl leaves
    ("CCCCc1ccc(CCO)cc1", "2-(4-butylphenyl)ethanol", "2-(4-丁基苯基)乙醇"),
    # methylthio
    (
        "CSc1ccc(CCO)cc1",
        "2-(4-methylsulfanylphenyl)ethanol",
        "2-(4-甲硫基苯基)乙醇",
    ),
    # cyano (alcohol parent, not formonitrile)
    ("N#Cc1ccc(CCO)cc1", "2-(4-cyanophenyl)ethanol", "2-(4-氰基苯基)乙醇"),
    # nested benzyl
    (
        "c1ccc(Cc2ccc(CCO)cc2)cc1",
        "2-(4-benzylphenyl)ethanol",
        "2-(4-苄基苯基)乙醇",
    ),
    (
        "Clc1ccc(Cc2ccc(CCO)cc2)cc1",
        "2-[4-(4-chlorobenzyl)phenyl]ethanol",
        "2-[4-(4-氯苄基)苯基]乙醇",
    ),
    # negatives
    ("CCc1ccccc1", "ethylbenzene", "乙基苯"),
    ("N#Cc1ccccc1", "benzonitrile", "苯甲腈"),
    ("Clc1ccc(CCO)cc1", "2-(4-chlorophenyl)ethanol", "2-(4-氯苯基)乙醇"),
    ("c1ccc(-c2ccc(CCO)cc2)cc1", "2-(4-phenylphenyl)ethanol", "2-(4-苯基苯基)乙醇"),
    ("COc1ccc(CCO)cc1", "2-(4-methoxyphenyl)ethanol", "2-(4-甲氧基苯基)乙醇"),
]


# ==========================================================================
# 合并自 test_principal_registry.py
# IUPAC: P-41/P-43
# Layer: L2
#
# Principal characteristic-group metadata is the sole priority authority.
# ==========================================================================
principal_registry__CASES = [
    ("CC(=O)OCC(=O)O", FG.ACID, FG.ESTER),
    ("N#CCC=O", FG.NITRILE, FG.ALDEHYDE),
]


def principal_registry___occ(group_class: FG) -> FunctionalGroupOccurrence:
    return FunctionalGroupOccurrence(group_class.value, group_class, frozenset(), frozenset(), {})


def principal_registry___inventory(*classes: FG) -> FunctionalGroupInventory:
    return FunctionalGroupInventory(tuple(principal_registry___occ(group_class) for group_class in classes))


@pytest.mark.parametrize("smiles,expected,other", principal_registry__CASES)
def test_principal_boundaries(smiles: str, expected: FG | None, other: FG) -> None:
    del smiles
    selected = select_principal_group(principal_registry___inventory(*(() if expected is None else (expected,)), other))
    assert (selected.group_class if selected else None) == expected


def test_registry_supports_new_p41_class_and_p43_subpath() -> None:
    """默认注册表未收录的类别可由外部注册表注入，并按优先级参与选择。"""
    classes = principal_registry___inventory(FG.NONE, FG.AMIDE)
    assert select_principal_group(classes).group_class == FG.AMIDE  # 未注册时不参与
    custom = dict(PRINCIPAL_REGISTRY)
    custom[FG.NONE] = PrincipalFeatureSpec(
        PrincipalPriority(7, (3, 2)), PrincipalExpression.SUFFIX,
    )
    selected = select_principal_group(classes, custom)
    assert selected is not None and selected.group_class == FG.NONE


def test_unregistered_kind_is_absent_not_guessed() -> None:
    """未注册 kind 一律返回 None：kind 只按注册表精确查表，不靠子串猜条目。"""
    for kind in ("carboxylate", "formamide_like", "amine_oxide", "unknown_one"):
        assert kr.get(kind) is None
    assert kr.parent_names("unknown_one") is None


# ==========================================================================
# 合并自 test_scaffold_registry_single_source.py
# IUPAC: P-22.2.1 / P-25
# Layer: L2
#
# Single source of truth: ScaffoldSpec drives kind_registry stems + retained ids.
#
# 唯一事实来源 _TEMPLATES 派生全部 ScaffoldSpec（25 个保留母体），编号标签
# standard_path 暂为空（L4 编号消费端为 stub）。
# ==========================================================================
def scaffold_registry_single_source___stem_specs():
    return [
        s for s in all_specs()
        if s.stem_en is not None and s.stem_zh is not None
    ]


@pytest.mark.parametrize("sp", scaffold_registry_single_source___stem_specs(), ids=lambda s: s.id)
def test_spec_stem_matches_kind_registry(sp) -> None:
    assert kr.parent_names(sp.id) == (sp.stem_en, sp.stem_zh)


def test_required_scaffolds_in_spec() -> None:
    for sid in ("benzene", "pyridine", "naphthalene", "indole"):
        assert get_spec(sid) is not None, sid


@pytest.mark.parametrize("sid", ["benzene", "pyridine", "naphthalene", "indole"])
def test_spec_meta(sid: str) -> None:
    sp = get_spec(sid)
    assert sp is not None
    assert sp.retained is True
    assert sp.ring in ("carbo", "hetero")
    assert sp.n_rings in (1, 2)


def scaffold_registry_single_source___ids_by_class() -> dict[str, frozenset[str]]:
    """按 naming_class 分组的 scaffold id 集合。"""
    groups: dict[str, set[str]] = {}
    for sp in all_specs():
        groups.setdefault(sp.naming_class, set()).add(sp.id)
    return {nc: frozenset(ids) for nc, ids in groups.items()}


def test_kind_ids_derive_from_specs() -> None:
    """spec 按 naming_class 分组；各组只钉锚点 id，registry 新增 spec 时断言自动跟随。"""
    by_class = scaffold_registry_single_source___ids_by_class()
    assert by_class["mono_carbo"] == {"benzene"}
    assert by_class["monohetero"] >= {"furan", "thiophene", "pyrrole", "pyridine",
                                      "pyridazine", "pyrimidine", "pyrazine",
                                      "imidazole", "pyrazole", "oxazole", "thiazole",
                                      "pyrrolidine", "piperidine", "morpholine",
                                      "piperazine", "oxolane", "oxane"}
    assert by_class["naph_family"] >= {"naphthalene", "quinoline", "isoquinoline",
                                       "quinazoline", "quinoxaline"}
    assert by_class["fused56"] >= {"indole", "indazole", "benzimidazole",
                                   "benzofuran", "benzothiophene",
                                   "benzothiazole", "benzoxazole"}
    assert by_class["anthra"] >= {"anthracene"}


def test_fused56_ids_retained_in_kind_registry() -> None:
    ids = scaffold_registry_single_source___ids_by_class()["fused56"]
    assert ids
    for sid in ids:
        m = kr.get(sid)
        assert m is not None, sid
        assert m.retained is True
        assert m.n_rings == 2


def test_monohetero_ids_retained_in_kind_registry() -> None:
    ids = scaffold_registry_single_source___ids_by_class()["monohetero"]
    assert ids
    for sid in ids:
        m = kr.get(sid)
        assert m is not None, sid
        assert m.retained is True
        assert m.ring == "hetero"
        assert m.n_rings == 1


def test_no_bidirectional_drift_stem_specs() -> None:
    """Every retained Spec with stem is in kind_registry with same stems."""
    for sp in all_specs():
        if sp.stem_en is None or sp.stem_zh is None:
            continue
        assert kr.parent_names(sp.id) == (sp.stem_en, sp.stem_zh)


def test_indole_parent_names_stable() -> None:
    assert kr.parent_names("indole") == ("1H-indole", "吲哚")


# ==========================================================================
# 合并自 test_retained_names.py
# IUPAC: P-62.2.1.1,P-66.1.1.4,P-29.2
# Layer: L3,L5
#
# 保留名（retained names）层：取代基保留名 anilino 与官能团衍生物保留前缀 carbamoyl/carbamoylamino/
# carbamoyloxy/carbamothioylamino/sulfamoyl 必须真正进入输出，不得退化为系统名或逐原子拼装。
#
# - P-62.2.1.1：phenylamino = anilino*（Glossary 818）；N-苯环带取代基时取代基前置于 anilino
#   （2-methylanilino* = (2-methylphenyl)amino，Glossary 447；4-[(4-hydroxyanilino)methyl]phenol 为 PIN）。
# - P-66.1.1.4.1：氨基甲酸（carbamic acid）的酰基保留前缀 carbamoyl；P-66.1.1.6 明确 ureido 不再
#   使用，优选 carbamoylamino（P_1 附录）。
# - P-66.1.1.4.2：磺酰胺保留前缀 sulfamoyl；(phenylamino)sulfonyl = phenylsulfamoyl*
#   （Glossary 820），N-取代基与 sulfamoyl 融合（丁基(甲基)sulfamoyl）。
#
# 说明：gold/ChEBI 全量对 anilino 取 41:0、carbamoyl 取 54 处、sulfamoyl 取 20 处，故这些保留名
# 是 benchmark 与 IUPAC 一致的口径；而 vinyl/isobutyl/tosyl 等 general/not_rec 级保留名 gold 一律
# 取系统名（见 anchored_table.resolve_name），不在此列。
# ==========================================================================
retained_names__ANILINO = [
    ("O=C(O)c1ccccc1Nc1ccccc1", "2-anilinobenzoic acid", "2-苯胺基苯甲酸"),
    ("O=C(O)c1ccccc1Nc1ccc(Cl)cc1", "2-(4-chloroanilino)benzoic acid", "2-(4-氯苯胺基)苯甲酸"),
    ("Cc1ccc(Cl)c(Nc2ccccc2C(=O)O)c1Cl",
     "2-(2,6-dichloro-3-methylanilino)benzoic acid", "2-(2,6-二氯-3-甲基苯胺基)苯甲酸"),
    ("O=C(O)c1ccccc1N(C(=O)C)c1ccccc1",
     "2-(N-acetylanilino)benzoic acid", "2-(N-乙酰基苯胺基)苯甲酸"),
]


@pytest.mark.parametrize("smiles,en,zh", retained_names__ANILINO)
def test_anilino_retained(smiles: str, en: str, zh: str) -> None:
    """anilino 取代基（裸/环取代/N-取代）输出保留名而非 phenylamino。"""
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)


# 锚定表保留前缀：片段级精确断言（碳/杂原子锚点）
retained_names__RETAINED_FRAG = [
    ("*C(=O)N", "carbamoyl", "氨基甲酰基"),
    ("*NC(=O)N", "carbamoylamino", "氨基甲酰氨基"),
    ("*OC(=O)N", "carbamoyloxy", "氨基甲酰氧基"),
    ("*NC(N)=S", "carbamothioylamino", "氨基硫代羰基氨基"),
    ("*S(=O)(=O)N", "sulfamoyl", "氨磺酰基"),
]


@pytest.mark.parametrize("smiles,en,zh", retained_names__RETAINED_FRAG)
def test_retained_fragment_prefixes(smiles: str, en: str, zh: str) -> None:
    """官能团衍生物保留前缀片段名精确匹配（不退化为 amino(oxo)methyl 等拼装式）。"""
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)


# 整分子：保留前缀嵌入组装名
retained_names__RETAINED_WHOLE = [
    ("O=C(N)c1ccccc1C(=O)O", "2-carbamoylbenzoic acid", "2-氨基甲酰基苯甲酸"),
    ("N=C(O)N[C@@H](CS)C(=O)O",
     "(2R)-2-(carbamoylamino)-3-sulfanylpropanoic acid", "(2R)-2-(氨基甲酰氨基)-3-巯基丙酸"),
]


@pytest.mark.parametrize("smiles,en,zh", retained_names__RETAINED_WHOLE)
def test_retained_prefix_whole_molecule(smiles: str, en: str, zh: str) -> None:
    """整分子组装名使用 carbamoyl/carbamoylamino 保留前缀。"""
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)


def test_sulfamoyl_n_substituted_fuses() -> None:
    """N-取代磺酰胺取 sulfamoyl 且 N-取代基与之融合（P-66.1.1.4.2）。"""
    r = SMILESNNamer().name("COC1=CC=C(C(=O)NCCS(NCC=2C=NC=CC2)(=O)=O)C=C1")
    assert r.success
    assert normalize_en(r.en) == normalize_en(
        "4-methoxy-N-[2-(pyridin-3-ylmethylsulfamoyl)ethyl]benzamide")
    assert normalize_zh(r.zh) == normalize_zh(
        "4-甲氧基-N-[2-(吡啶-3-基甲氨基磺酰基)乙基]苯甲酰胺")


def test_general_level_retained_stays_systematic() -> None:
    """general/not_rec 级保留名（vinyl/isobutyl/tosyl）在 general 模式下仍取系统名，避免回退。"""
    cases = [
        ("C=Cc1ccccc1", "ethenylbenzene"),          # vinyl → gold 全量取 ethenyl（非 vinyl）
        ("CC(C)Cc1ccccc1", "2-methylpropyl"),       # isobutyl → 2-methylpropyl（非 isobutyl）
    ]
    for smiles, token in cases:
        r = SMILESNNamer().name(smiles)
        assert r.success
        assert token in normalize_en(r.en)
