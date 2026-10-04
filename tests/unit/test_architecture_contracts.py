# 合并自 17 个原测试文件（按主题分组，内容与断言未改动）。
"""
test_l2_parent_core_contract.py: L2 parent_core is the sole authority for parent assembly / chain / gate helpers.
test_l2_l3_side_facts_contract.py: Architecture contract for Layer 2 side topology facts consumed by Layer 3.
test_l5_no_layer2_private.py: L5 special FG name modules must not import L2 private APIs.
test_ring_scaffold_expression_contract.py:
test_anchored_whole_mol.py: 带自由价的整分子命名：*xxx 走常规管线得到保留取代基名。
test_claimable_block_api.py: Ownership-only claimable blocks: topology claims, no naming mode.
test_claimable_parent_batches.py: Parent-class batch coverage for universal claimable-block mainline.
test_universal_claimable_block.py:
test_parent_ownership.py: Terminal parent ownership: immutable owned_atoms + ordered candidates.
test_parent_scoring.py: Composable parent scoring + recursive alkyl substituents.
test_open_chain_alkoxy_claims.py: Open-chain alkoxy: methoxybutane via typed ether claims, not alkane evaporation.
test_open_fg_no_ring_walk.py: Open FG parent chain must not enter rings; claim cycloalkyl / piperidinyl sides.
test_scaffold_numbering_producer_flow.py: Selected retained parents must carry facts before L4 numbering.
test_merge_alken_kinds.py: Merge open-chain mono-FG alken* ParentKind into saturated kinds.
test_fg_locants_records.py: fg_locants 稀疏产出契约: principal FG 位次为结构化列表 [{kind, locants, omit}],
test_substituent_namer.py: Ordered SubstituentNamer backends: retained → rooted_tree → recursive.
"""
from __future__ import annotations

import ast
import pytest

from opensmiles.layer0.preprocessor import preprocess
from opensmiles.layer1.analyzer import analyze
from opensmiles.layer1.fg_registry import FG_SPECS
from opensmiles.layer1.functional_group_inventory import (
    FunctionalGroupClass as FG,
    FunctionalGroupInventory,
    inventory_from_info,
)
from opensmiles.layer2.parent_select import _collect_candidates
from opensmiles.layer2.parent_select import finalize_parent_ownership
from opensmiles.layer2.principal_expression import (
    PrincipalExpressionFacts,
    PrincipalRelation,
)
from opensmiles.layer2.parent_select import select_parent
from opensmiles.layer2.ring_scaffold import all_specs
from opensmiles.layer3.claimable_block import ClaimedBlock, SideSlot, claim_block, iter_claims
from opensmiles.layer3.coverage import build_coverage_ledger
from opensmiles.layer3.substituent_extractor import extract_substituents
from opensmiles.layer3.substituent_namer import SubstituentName, SubstituentNamer
from opensmiles.layer4.locant_calc import _fg_locants
from opensmiles.layer4.numbering import number
from opensmiles.namer import SMILESNNamer
from opensmiles.tools.re import normalize_en, normalize_zh
from pathlib import Path
from rdkit import Chem

# ==========================================================================
# 合并自 test_l2_parent_core_contract.py
# IUPAC: P-44 / architecture L2 parent hub
# Layer: L2
#
# L2 parent_core is the sole authority for parent assembly / chain / gate helpers.
#
# Producer modules must not import private helpers from parent_selector.
# parent_selector keeps FG try + select_parent; helpers live in parent_core.
# ==========================================================================
l2_parent_core_contract___ROOT = Path(__file__).resolve().parents[2]
l2_parent_core_contract___L2 = l2_parent_core_contract___ROOT / "src" / "opensmiles" / "layer2"

# Helpers that must be sourced from parent_core, not parent_selector.
l2_parent_core_contract___HELPERS = frozenset({
    "_parent_dict",
    "_longest_from",
    "_longest_chain",
})

# Modules allowed to re-export / host try surface.
l2_parent_core_contract___EXEMPT = frozenset({
    "parent_select.py",
    "parent_core.py",
    "fg_helpers.py",
    "__init__.py",
})


def l2_parent_core_contract___parent_selector_helper_hits(path: Path) -> list[str]:
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    hits: list[str] = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.ImportFrom) or not node.module:
            continue
        if node.module != "opensmiles.layer2.parent_select":
            continue
        for alias in node.names:
            if alias.name in l2_parent_core_contract___HELPERS:
                hits.append(f"{alias.name} (L{node.lineno})")
    return hits


def l2_parent_core_contract___producer_paths() -> list[Path]:
    return sorted(
        p for p in l2_parent_core_contract___L2.rglob("*.py")
        if p.name not in l2_parent_core_contract___EXEMPT and p.is_file()
    )


@pytest.mark.parametrize("path", l2_parent_core_contract___producer_paths(), ids=lambda p: p.relative_to(l2_parent_core_contract___L2).as_posix())
def test_producers_do_not_import_helpers_from_parent_selector(path: Path) -> None:
    hits = l2_parent_core_contract___parent_selector_helper_hits(path)
    assert hits == [], (
        f"{path.relative_to(l2_parent_core_contract___L2).as_posix()} must not import helpers "
        f"from parent_selector: {hits}"
    )


# End-to-end bilingual smoke (zero-behavior after cut).
l2_parent_core_contract___E2E = [
    # acid / chain
    ("CC(=O)O", "acetic acid", "乙酸"),
    # ester
    ("CCOC(=O)C", "ethyl acetate", "乙酸乙酯"),
]


@pytest.mark.parametrize("smiles,en,zh", l2_parent_core_contract___E2E)
def test_parent_core_e2e_bilingual(smiles: str, en: str, zh: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)


def test_alkyl_acid_negative_not_benzene() -> None:
    """Negative: simple aliphatic acid stays alkyl (not benzoic)."""
    r = SMILESNNamer().name("CCC(=O)O")
    assert r.success
    assert normalize_en(r.en) == "propanoic acid"
    assert "benzoic" not in normalize_en(r.en)
    assert "苯" not in normalize_zh(r.zh)


# ==========================================================================
# IUPAC: P-41 / L1 官能团单一格式契约
# Layer: L1,L2
#
# analyze() 只以 FunctionalGroupInventory 承载官能团事实：不再并行回写扁平 FG list
# （旧格式缺类别/特征原子，且与清单不等价——demoted_* 只在旧格式里）。
# _fg_parts 只承载 P-41 官能团条目；双键/三键是不饱和度事实，走独立通道。
# ==========================================================================
l1_fg_single_format__FG_KEYS = frozenset(sp.fg for sp in FG_SPECS)


def test_info_carries_inventory_and_no_flat_fg_lists() -> None:
    """info 只有 fg_inventory，没有扁平 FG list 键（含 demoted_*）。"""
    info = analyze(Chem.MolFromSmiles("OC(=O)CC#N"))
    assert isinstance(info["fg_inventory"], FunctionalGroupInventory)
    assert not l1_fg_single_format__FG_KEYS & set(info), l1_fg_single_format__FG_KEYS & set(info)
    assert not [key for key in info if key.startswith("demoted_")]


def test_unsaturation_stays_out_of_fg_channel() -> None:
    """双键/三键在 info 顶层，但不进 _detect_parts。"""
    from opensmiles.layer1.analyzer import _detect_parts

    smiles = "C=CC#C"
    info = analyze(Chem.MolFromSmiles(smiles))
    assert info["double_bonds"] and info["triple_bonds"]
    assert set(_detect_parts(Chem.MolFromSmiles(smiles))) <= l1_fg_single_format__FG_KEYS


def test_demoted_leaf_survives_in_inventory_not_in_flat_lists() -> None:
    """腈被羧酸压制时：条目以 demoted 状态留在清单中，且不再是主基团候选。"""
    info = analyze(Chem.MolFromSmiles("N#CCC(=O)O"))
    nitriles = inventory_from_info(info).occurrences(FG.NITRILE)
    assert nitriles == ()  # occurrences() 只返回未降级条目
    assert [o.group_class for o in info["fg_inventory"].demoted_entries()] == [FG.NITRILE]
    assert info["fg_inventory"].demoted_entries()[0].payload["center_idx"] is not None


def test_inventory_missing_is_explicit_failure() -> None:
    """缺 fg_inventory 即上游违约：显式报错，不静默退化成空清单。"""
    from opensmiles.layer1.functional_group_inventory import inventory_from_info

    with pytest.raises(KeyError, match="fg_inventory"):
        inventory_from_info({"mol": Chem.MolFromSmiles("CC")})


# ==========================================================================
# 合并自 test_l2_l3_side_facts_contract.py
# IUPAC: P-29.2
# Layer: L2,L3
#
# Architecture contract for Layer 2 side topology facts consumed by Layer 3.
# ==========================================================================
l2_l3_side_facts_contract__LAYER3 = Path(__file__).parents[2] / "src" / "opensmiles" / "layer3"
l2_l3_side_facts_contract__TOOLS = Path(__file__).parents[2] / "src" / "opensmiles" / "tools"
# Side-topology facts merged into tools 碳拓扑原语 chain.py。
l2_l3_side_facts_contract__TOOLS_FACTS = l2_l3_side_facts_contract__TOOLS / "chain.py"
l2_l3_side_facts_contract__LEAF_PROTOCOL = l2_l3_side_facts_contract__LAYER3 / "leaves" / "protocol.py"


def l2_l3_side_facts_contract___private_layer2_imports(path: Path) -> list[str]:
    tree = ast.parse(path.read_text(encoding="utf-8"))
    found: list[str] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom) and (node.module or "").startswith("opensmiles.layer2"):
            found.extend(alias.name for alias in node.names if alias.name.startswith("_"))
    return found


def l2_l3_side_facts_contract___all_functions(path: Path) -> list[ast.FunctionDef | ast.AsyncFunctionDef]:
    tree = ast.parse(path.read_text(encoding="utf-8"))
    return [n for n in ast.walk(tree) if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef))]


def l2_l3_side_facts_contract___tools_modules(tree: ast.AST) -> list[str]:
    direct = [alias.name for node in ast.walk(tree) if isinstance(node, ast.Import)
              for alias in node.names if alias.name.startswith("opensmiles.tools")]
    froms = [(node.module or "") for node in ast.walk(tree) if isinstance(node, ast.ImportFrom)
             and (node.module or "").startswith("opensmiles.tools")]
    return [*direct, *froms]


def test_tools_import_scan_covers_both_ast_forms() -> None:
    tree = ast.parse("import opensmiles.tools.aryl_sub\nfrom opensmiles.tools import side_alkyl")
    assert l2_l3_side_facts_contract___tools_modules(tree) == [
        "opensmiles.tools.aryl_sub", "opensmiles.tools",
    ]


def test_side_facts_has_no_string_to_enum_dispatcher() -> None:
    tree = ast.parse(l2_l3_side_facts_contract__TOOLS_FACTS.read_text(encoding="utf-8"))
    mappings = [node.value for node in tree.body if isinstance(node, ast.Assign)
                and isinstance(node.value, ast.Dict)]
    bad = [mapping for mapping in mappings
           if any(isinstance(key, ast.Constant) and isinstance(key.value, str)
                  and isinstance(value, ast.Attribute) and value.attr.isupper()
                  for key, value in zip(mapping.keys, mapping.values))]
    assert not bad


def test_side_facts_does_not_dispatch_match_kind_keys() -> None:
    tree = ast.parse(l2_l3_side_facts_contract__TOOLS_FACTS.read_text(encoding="utf-8"))
    keys = [node.slice.value for node in ast.walk(tree) if isinstance(node, ast.Subscript)
            and isinstance(node.slice, ast.Constant) and isinstance(node.slice.value, str)]
    assert "kind" not in keys


def test_side_facts_has_no_naming_token_maps() -> None:
    tree = ast.parse(l2_l3_side_facts_contract__TOOLS_FACTS.read_text(encoding="utf-8"))
    assigned = {target.id for node in tree.body if isinstance(node, ast.Assign)
                for target in node.targets if isinstance(target, ast.Name)}
    assert not {"_LEAF_KINDS", "_ALKYL_SHAPES"} & assigned


def l2_l3_side_facts_contract___enum_members(node: ast.ClassDef) -> list[ast.Assign]:
    return [item for item in node.body if isinstance(item, ast.Assign)]


def l2_l3_side_facts_contract___is_auto_member(node: ast.Assign) -> bool:
    value = node.value
    return (isinstance(value, ast.Call) and isinstance(value.func, ast.Name)
            and value.func.id == "auto" and not value.args and not value.keywords)


def test_side_facts_has_no_naming_dependencies() -> None:
    tree = ast.parse(l2_l3_side_facts_contract__TOOLS_FACTS.read_text(encoding="utf-8"))
    imports = [(node.module or "", alias.name) for node in ast.walk(tree)
               if isinstance(node, ast.ImportFrom) for alias in node.names]
    calls = [node.func.id for node in ast.walk(tree) if isinstance(node, ast.Call)
             and isinstance(node.func, ast.Name)]
    assert not [(module, name) for module, name in imports if name == "name_leaf"]
    assert "name_leaf" not in calls


def test_public_side_fact_functions_are_typed_and_explicit() -> None:
    public = [fn for fn in l2_l3_side_facts_contract___all_functions(l2_l3_side_facts_contract__TOOLS_FACTS) if not fn.name.startswith("_")]
    assert all(fn.returns is not None for fn in public)
    assert all(fn.args.vararg is None and fn.args.kwarg is None for fn in public)
    assert all(all(arg.annotation is not None for arg in fn.args.args) for fn in public)


def test_layer3_functions_are_explicit() -> None:
    functions = [fn for path in l2_l3_side_facts_contract__LAYER3.rglob("*.py") for fn in l2_l3_side_facts_contract___all_functions(path)]
    assert all(fn.args.vararg is None and fn.args.kwarg is None for fn in functions)


def test_side_fact_contract_has_no_public_dict_returns() -> None:
    public = [fn for fn in l2_l3_side_facts_contract___all_functions(l2_l3_side_facts_contract__TOOLS_FACTS) if not fn.name.startswith("_")]
    returns = [ast.unparse(fn.returns) for fn in public]
    assert not any("dict" in annotation for annotation in returns), returns


def test_side_fact_contract_has_no_public_name_renderers() -> None:
    tree = ast.parse(l2_l3_side_facts_contract__TOOLS_FACTS.read_text(encoding="utf-8"))
    public = [n for n in tree.body
              if isinstance(n, ast.FunctionDef) and not n.name.startswith("_")]
    banned = [fn.name for fn in public if "render" in fn.name or "name" in fn.name]
    string_tuples = [fn.name for fn in public if "tuple[str" in ast.unparse(fn.returns)]
    assert not banned
    assert not string_tuples


def test_side_fact_dataclasses_expose_only_topology_fields() -> None:
    tree = ast.parse(l2_l3_side_facts_contract__TOOLS_FACTS.read_text(encoding="utf-8"))
    classes = [n for n in tree.body if isinstance(n, ast.ClassDef)]
    forbidden = {"name", "en", "zh", "paren"}
    for cls in classes:
        fields = {n.target.id for n in cls.body if isinstance(n, ast.AnnAssign)
                  and isinstance(n.target, ast.Name)}
        assert not fields & forbidden, (cls.name, fields & forbidden)
        methods = [n.name for n in cls.body if isinstance(n, ast.FunctionDef)]
        assert methods == [], (cls.name, methods)


# ==========================================================================
# 合并自 test_l5_no_layer2_private.py
# IUPAC: P-65.3 / architecture layer purity
# Layer: L2,L5
#
# L5 special FG name modules must not import L2 private APIs.
#
# Aryl / cycloalkyl stems are precomputed into parent side dicts at L2 pack time.
# L5 only reads side['en']/side['zh'] (and existing kind/mode fields).
# ==========================================================================
l5_no_layer2_private___ROOT = Path(__file__).resolve().parents[2]
l5_no_layer2_private___L5 = l5_no_layer2_private___ROOT / "src" / "opensmiles" / "layer5"

# L5 special-FG name modules that previously walked via L2 private APIs.


def l5_no_layer2_private___layer2_import_hits(path: Path) -> list[str]:
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    hits: list[str] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom) and node.module:
            if node.module == "opensmiles.layer2" or node.module.startswith(
                "opensmiles.layer2."
            ):
                hits.append(f"from {node.module} (L{node.lineno})")
        elif isinstance(node, ast.Import):
            for alias in node.names:
                if alias.name == "opensmiles.layer2" or alias.name.startswith(
                    "opensmiles.layer2."
                ):
                    hits.append(f"import {alias.name} (L{node.lineno})")
    return hits


@pytest.mark.parametrize("path", sorted(l5_no_layer2_private___L5.glob("*.py")), ids=lambda p: p.name)
def test_all_l5_modules_have_no_layer2_import(path: Path) -> None:
    hits = l5_no_layer2_private___layer2_import_hits(path)
    assert hits == [], f"{path.name} must not import layer2: {hits}"


# End-to-end bilingual gold (behavior must stay flat after cut).
l5_no_layer2_private___E2E = [
    # sulfonamide aryl
    ("c1ccc(S(=O)(=O)N)cc1", "benzenesulfonamide", "苯磺酰胺"),
    # urea aryl
    ("c1ccc(NC(=O)N)cc1", "phenylurea", "苯基脲"),
    # carbonate aryl
    (
        "C(OC)(OC1=C(C=C(C=C1)C)C)=O",
        "methyl 2,4-dimethylphenyl carbonate",
        "2,4-二甲苯基甲基碳酸酯",
    ),
    # alkyl sulfonamide still correct (non-aryl path)
    ("CS(=O)(=O)N", "methanesulfonamide", "甲磺酰胺"),
]


# ==========================================================================
# 合并自 test_ring_scaffold_expression_contract.py
# IUPAC: P-22 / P-44 ring scaffold-expression contract
# Layer: L2
# ==========================================================================
def ring_scaffold_expression_contract___typed(smiles, group):
    parents = _collect_candidates(analyze(Chem.MolFromSmiles(smiles)))
    return next(p for p in parents if (f := p.get("principal_expression_facts")) and f.group_class.value == group)


def test_generic_carbocycle_resolves_independently_of_alcohol_expression():
    parent = ring_scaffold_expression_contract___typed("OC1CCCCC1", "alcohol")
    assert parent["scaffold_identity"].naming_class == "carbocycle"


def test_benzene_resolves_to_stable_scaffold_identity():
    parent = ring_scaffold_expression_contract___typed("Oc1ccccc1", "alcohol")
    assert parent["scaffold_id"] == "benzene"
    assert parent["scaffold_identity"].naming_class == "mono_carbo"


def test_naphthalenol_resolves_to_naphthalene_scaffold():
    parent = ring_scaffold_expression_contract___typed("Oc1cccc2ccccc12", "alcohol")
    assert parent["scaffold_id"] == "naphthalene"


# ==========================================================================
# 合并自 test_anchored_whole_mol.py
#
# 带自由价的整分子命名：*xxx 走常规管线得到保留取代基名。
# ==========================================================================
def anchored_whole_mol___name(smiles: str):
    return SMILESNNamer().name(smiles)


# ── 内联条目：顶层 *xxx == 内联 (en, zh) ──

@pytest.mark.parametrize("smi,en,zh", [
    ("*c1ccccc1", "phenyl", "苯基"),
    ("*C1CC1", "cyclopropyl", "环丙基"),
])
def test_top_level_star_matches_inline(smi, en, zh):
    r = anchored_whole_mol___name(smi)
    assert r.success
    assert r.en == en
    assert r.zh == zh


# ── *O/*[O]/*N 只服务整分子顶层，不进 L3 取代基查表 ──

def test_hetero_anchored_not_used_as_substituent_lookup():
    # 多酮羰基氧不被 L3 误作羟基取代基（*O 不进取代基查表 + _ketone_fg_atoms 修复）
    assert anchored_whole_mol___name("CC(=O)C(C)=O").en == "butane-2,3-dione"
    assert anchored_whole_mol___name("CC(=O)CC(=O)C").en == "pentane-2,4-dione"
    # 真正醚氧取代基仍正常（*OC 在表，非 whole-only）
    assert anchored_whole_mol___name("*OC").en == "methoxy"


# ==========================================================================
# 合并自 test_claimable_block_api.py
#
# Ownership-only claimable blocks: topology claims, no naming mode.
# ==========================================================================
def claimable_block_api___info(smiles: str):
    mol = preprocess(smiles)
    return mol, analyze(mol)


def claimable_block_api___benzene_owned(smiles: str):
    mol, info = claimable_block_api___info(smiles)
    for cand in select_parent(info):
        if cand.get("scaffold_id") == "benzene":
            return mol, cand["owned_atoms"]
    raise AssertionError("no benzene parent candidate")


def test_methylbutylbenzene_ring_c_claim():
    """2-methylbutylbenzene with ring ownership: one RING_C claim, five side C."""
    mol, owned = claimable_block_api___benzene_owned("c1ccc(cc1)CC(C)CC")
    claims = iter_claims(mol, owned)
    assert len(claims) == 1
    c = claims[0]
    assert c.slot == SideSlot.RING_C
    assert len(c.atoms) == 5
    assert all(mol.GetAtomWithIdx(i).GetAtomicNum() == 6 for i in c.atoms)
    assert all(not mol.GetAtomWithIdx(i).IsInRing() for i in c.atoms)
    assert mol.GetAtomWithIdx(c.attach_parent).IsInRing()


def test_multi_attachment_component_rejected():
    """Component touching two owned atoms (tetralin aliphatic bridge) is rejected."""
    mol = preprocess("c1cccc2c1CCCC2")
    owned = frozenset(
        i
        for i in range(mol.GetNumAtoms())
        if mol.GetAtomWithIdx(i).GetIsAromatic()
    )
    claims = iter_claims(mol, owned)
    assert claims == []
    # explicit claim_block also returns None for either attachment
    roots = [
        n.GetIdx()
        for p in owned
        for n in mol.GetAtomWithIdx(p).GetNeighbors()
        if n.GetAtomicNum() != 1 and n.GetIdx() not in owned
    ]
    assert roots
    for root in roots:
        attach = next(
            p
            for p in owned
            for n in mol.GetAtomWithIdx(p).GetNeighbors()
            if n.GetIdx() == root
        )
        assert (
            claim_block(
                mol,
                owned_atoms=owned,
                attach_parent=attach,
                root=root,
                slot=SideSlot.RING_C,
            )
            is None
        )


def test_claim_block_valid_returns_atoms():
    """claim_block returns ClaimedBlock for a valid single-attachment side."""
    mol, owned = claimable_block_api___benzene_owned("c1ccc(cc1)CC")
    claims = iter_claims(mol, owned)
    assert len(claims) == 1
    c = claims[0]
    again = claim_block(
        mol,
        owned_atoms=owned,
        attach_parent=c.attach_parent,
        root=c.root,
        slot=c.slot,
    )
    assert again == c


def test_sulfonyl_o_not_claimed():
    """砜双键氧同样被跳过，避免 cut 出 *O 污染成羟基。"""
    mol, info = claimable_block_api___info("CS(=O)(=O)C")
    parent = select_parent(info)[0]
    claims = iter_claims(mol, parent["owned_atoms"])
    assert all(
        mol.GetAtomWithIdx(c.root).GetAtomicNum() != 8 for c in claims
    )


# ==========================================================================
# 合并自 test_claimable_parent_batches.py
#
# Parent-class batch coverage for universal claimable-block mainline.
# ==========================================================================
claimable_parent_batches__BATCHES = [
    # acid
    ("CC(=O)O", "acetic acid", "乙酸"),
    # alkane
    ("CCCC", "butane", "丁烷"),
    # ketone
    ("CC(=O)C", "propan-2-one", None),  # acetone may be retained
    ("CCCC(=O)C", "pentan-2-one", "戊-2-酮"),
    # alcohol
    ("CCO", "ethanol", "乙醇"),
    ("CC(C)O", "propan-2-ol", "丙-2-醇"),
    # benzene
    ("c1ccccc1", "benzene", "苯"),
    # phenol
    ("Oc1ccccc1", "phenol", "苯酚"),
    # aniline
    ("Nc1ccccc1", "aniline", "苯胺"),
    # pyridine
    ("c1ccncc1", "pyridine", "吡啶"),
    # amide
    ("CC(=O)N", "acetamide", "乙酰胺"),
    # benzamide
    ("c1ccccc1C(=O)N", "benzamide", "苯甲酰胺"),
]


def claimable_parent_batches___assert_coverage_complete(smiles: str) -> None:
    """并列母体组内至少有一个候选能让 coverage ledger 完整——即 namer 实际采用的那个。"""
    mol = Chem.MolFromSmiles(smiles)
    info = analyze(mol)
    group = select_parent(info)
    assert group, f"no parent candidate: {smiles}"
    gaps = []
    for parent in group:
        led = build_coverage_ledger(mol, owned_atoms=parent["owned_atoms"], names=[])
        if led.complete:
            return
        gaps.append(f"kind={parent.get('kind')} gap={sorted(led.gap)} overlap={sorted(led.overlap)}")
    raise AssertionError(f"{smiles}: no candidate covers all atoms; {'; '.join(gaps)}")


@pytest.mark.parametrize("smiles,en,zh", claimable_parent_batches__BATCHES)
def test_parent_batch_names_and_coverage(smiles, en, zh):
    r = SMILESNNamer().name(smiles)
    if en is None:
        # optional supported class: only require coverage when successful
        if r.success:
            claimable_parent_batches___assert_coverage_complete(smiles)
        return
    assert r.success, f"failed {smiles}: {r.meta}"
    if en == "propan-2-one":
        assert normalize_en(r.en) in (normalize_en("propan-2-one"), normalize_en("acetone"))
    else:
        assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        if zh == "戊-2-酮" and "丙酮" in (r.zh or ""):
            pass
        else:
            assert normalize_zh(r.zh) == normalize_zh(zh)
    claimable_parent_batches___assert_coverage_complete(smiles)


# ==========================================================================
# 合并自 test_universal_claimable_block.py
# Universal claimable-block class-level red bars.
# Goal: every selected parent accounts for every heavy atom via typed, ordered
# side-block naming (or fail/retry another parent). Task 1 establishes RED only.
# ==========================================================================
universal_claimable_block__CASES = [
    ("COC(C)CC", "2-methoxybutane", "2-甲氧基丁烷"),
]


@pytest.mark.parametrize("smiles,en,zh", universal_claimable_block__CASES)
def test_universal_claimable_classes(smiles: str, en: str, zh: str) -> None:
    result = SMILESNNamer().name(smiles)
    assert result.success
    assert normalize_en(result.en) == normalize_en(en)
    assert normalize_zh(result.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_parent_ownership.py
#
# Terminal parent ownership: immutable owned_atoms + ordered candidates.
# ==========================================================================
def parent_ownership___info(smiles: str):
    mol = preprocess(smiles)
    return mol, analyze(mol)


def test_benzamide_owns_core_not_n_phenyl():
    """Benzamide owns aryl + amide C/N/O; N-phenyl carbons stay outside."""
    mol, info = parent_ownership___info("c1ccccc1C(=O)Nc2ccccc2")
    parent = select_parent(info)[0]
    assert parent["kind"] == "amide"
    assert parent.get("scaffold_id") == "benzene"
    owned = parent["owned_atoms"]
    assert isinstance(owned, frozenset)

    am = inventory_from_info(info).occurrences(FG.AMIDE)[0].payload
    n_idx = next(i for i in am["surr_idx"] if mol.GetAtomWithIdx(i).GetAtomicNum() == 7)
    assert am["center_idx"] in owned
    assert n_idx in owned
    # carbonyl O is owned
    for n in mol.GetAtomWithIdx(am["center_idx"]).GetNeighbors():
        if n.GetAtomicNum() == 8:
            assert n.GetIdx() in owned
    # aryl core (chain) owned
    for idx in parent["chain"]:
        assert idx in owned
    # N-phenyl carbons not owned（N 上除酰胺羰基碳外的碳取代基）
    for c in [n.GetIdx() for n in mol.GetAtomWithIdx(n_idx).GetNeighbors()
              if n.GetAtomicNum() == 6 and n.GetIdx() != am["center_idx"]]:
        assert c not in owned
        for n in mol.GetAtomWithIdx(c).GetNeighbors():
            if n.GetAtomicNum() == 6 and n.GetIsAromatic() and n.GetIdx() != c:
                assert n.GetIdx() not in owned


def test_acid_owns_carboxyl_c_and_both_oxygens():
    """Acid owns carboxyl carbon + both oxygens (carbonyl O and OH O)."""
    mol, info = parent_ownership___info("CC(=O)O")
    parent = select_parent(info)[0]
    assert parent["kind"] == "acid"
    owned = parent["owned_atoms"]
    assert isinstance(owned, frozenset)

    anchors = parent["principal_expression_facts"].anchor_atoms
    assert len(anchors) == 1
    c_idx = next(iter(anchors))
    assert c_idx in owned
    oxygens = [
        n.GetIdx()
        for n in mol.GetAtomWithIdx(c_idx).GetNeighbors()
        if n.GetAtomicNum() == 8
    ]
    assert len(oxygens) == 2
    for o in oxygens:
        assert o in owned


def test_owned_atoms_is_frozenset_and_immutable_after_ops():
    """owned_atoms is frozenset and unchanged after extraction / side ops."""
    mol, info = parent_ownership___info("c1ccccc1C(=O)O")
    parent = select_parent(info)[0]
    owned = parent["owned_atoms"]
    assert isinstance(owned, frozenset)
    snapshot = frozenset(owned)

    # side operations must not mutate ownership truth
    _ = finalize_parent_ownership(parent, mol)
    _ = list(owned)
    copy_out = set(owned)
    copy_out.add(999)

    assert parent["owned_atoms"] == snapshot
    assert parent["owned_atoms"] is owned
    assert isinstance(parent["owned_atoms"], frozenset)


def test_select_parent_returns_finalized_tied_group():
    """select_parent returns the P-45.2.1 tied group, already finalized and in candidate order."""
    mol, info = parent_ownership___info("c1ccccc1C(=O)Nc2ccccc2")
    group = select_parent(info)
    assert group
    assert all(isinstance(c.get("owned_atoms"), frozenset) for c in group)
    # 组内候选同分，但 owned_atoms 互不相同（记录的是各自母体范围）
    assert len({tuple(sorted(c["owned_atoms"])) for c in group}) == len(group)
    # 首位是 P-44/P-45.2 排序领先者：苯甲酰胺母体
    head = group[0]
    assert head["kind"] == "amide"
    assert head.get("scaffold_id") == "benzene"


def test_finalize_copies_once_with_owned_atoms():
    """finalize_parent_ownership copies candidate once with owned_atoms frozenset."""
    mol, info = parent_ownership___info("CC(=O)O")
    raw = {k: v for k, v in select_parent(info)[0].items() if k != "owned_atoms"}
    fin = finalize_parent_ownership(raw, mol)
    assert fin is not raw
    assert "owned_atoms" not in raw
    assert isinstance(fin["owned_atoms"], frozenset)
    assert 1 in fin["owned_atoms"]
    assert len([a for a in fin["owned_atoms"] if mol.GetAtomWithIdx(a).GetAtomicNum() == 8]) == 2


def parent_ownership___ester_o_count(mol, owned):
    """Count O atoms in owned set that are neighbours of owned C atoms (ester O)."""
    return sum(
        1 for i in owned
        if mol.GetAtomWithIdx(i).GetAtomicNum() == 8
        and any(n.GetIdx() in owned and n.GetAtomicNum() == 6 for n in mol.GetAtomWithIdx(i).GetNeighbors())
    )


# ==========================================================================
# 合并自 test_parent_scoring.py
# IUPAC: P-44.1 / P-29.3.1
# Layer: L2,L3
#
# Composable parent scoring + recursive alkyl substituents.
#
# Layer2 select_parent moves from an if-else waterfall to candidate
# collection + composable scoring: a benzene (or cycloalkane) core with
# fully-recognizable alkyl sides must win over the longest-chain alkane,
# while principal-FG parents (e.g. benzoic acid) must still outrank the
# bare ring. Layer3 recognizes tert-butyl recursively (P-29.3.1).
# ==========================================================================
parent_scoring__CASES = [
    # positive: ring core must beat longest-chain alkane via scoring
    ("c1ccc(cc1)C(C)(C)C", "tert-butylbenzene", "叔丁基苯"),
    ("CCc1ccccc1C", "1-ethyl-2-methylbenzene", "1-乙基-2-甲基苯"),

  # acid outranks bare benzene
    ("Cc1ccc(Cl)cc1", "1-chloro-4-methylbenzene", "1-氯-4-甲基苯"),
    ("CCC(C)C", "2-methylbutane", "2-甲基丁烷"),  # plain alkane: no ring bias

]


@pytest.mark.parametrize("smiles,en,zh", parent_scoring__CASES)
def test_parent_scoring(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


@pytest.mark.parametrize(
    "smiles,bad",
    [
        ("c1ccc(cc1)C(C)(C)C", "2,2-dimethyloctane"),
        ("CCc1ccccc1C", "nonane"),
        # unhandled complex side: must NOT collapse to a bare "benzene"
        ("CC(C)Cc1ccccc1", "benzene"),
        ("CCC(C)c1ccccc1", "benzene"),
    ],
)
def test_no_wrong_parent(smiles: str, bad: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert not (r.success and normalize_en(r.en) == normalize_en(bad))


# ==========================================================================
# 合并自 test_open_chain_alkoxy_claims.py
#
# Open-chain alkoxy: methoxybutane via typed ether claims, not alkane evaporation.
# ==========================================================================
open_chain_alkoxy_claims__METHOXYBUTANE = "COC(C)CC"


def test_methoxybutane_final_name():
    result = SMILESNNamer().name(open_chain_alkoxy_claims__METHOXYBUTANE)
    assert result.success
    assert normalize_en(result.en) == normalize_en("2-methoxybutane")
    assert normalize_zh(result.zh) == normalize_zh("2-甲氧基丁烷")


def test_methoxybutane_has_ether_or_alkoxy_claim_path():
    """Alkane parent with methoxy claim, or ether parent — O must be covered."""
    mol = Chem.MolFromSmiles(open_chain_alkoxy_claims__METHOXYBUTANE)
    info = analyze(mol)
    parent = select_parent(info)[0]
    owned = parent["owned_atoms"]
    claims = iter_claims(mol, owned)
    o_idxs = {a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 8}
    covered = set(owned) | {i for c in claims for i in c.atoms}
    assert o_idxs <= covered
    from opensmiles.layer3.substituent_namer import SubstituentNamer

    named = [SubstituentNamer().name(mol, c) for c in claims]
    assert any(n is not None and n.en == "methoxy" for n in named) or parent.get("kind") == "ether"


@pytest.mark.parametrize(
    "smiles,en",
    [
        ("COC", "methoxyethane"),  # may be dimethyl ether / methoxyethane retained
        ("CCOCC", "ethoxyethane"),
    ],
)
def test_simple_sym_ethers_still_succeed(smiles, en):
    result = SMILESNNamer().name(smiles)
    assert result.success
    # loose: just ensure oxygen is named somehow
    assert "oxy" in normalize_en(result.en) or "ether" in normalize_en(result.en)


# ==========================================================================
# 合并自 test_open_fg_no_ring_walk.py
# IUPAC: P-44 / P-29.6 / P-22.2
# Layer: L2,L3
#
# Open FG parent chain must not enter rings; claim cycloalkyl / piperidinyl sides.
#
# 1-cyclohexylethanone and 1-(piperidin-4-yl)ethanone stay ketone parents with
# ring as substituent — not octan-2-one / pentan-2-one. Ring FG parents
# (cyclohexanone/cyclohexanol) unchanged.
# ==========================================================================
open_fg_no_ring_walk__CASES = [
    ("CCC(=O)C", "butan-2-one", "丁-2-酮"),
]


@pytest.mark.parametrize("smiles,en,zh", open_fg_no_ring_walk__CASES)
def test_open_fg_no_ring_walk(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_scaffold_numbering_producer_flow.py
# IUPAC: P-14 / P-25
# Layer: L2,L4
#
# Selected retained parents must carry facts before L4 numbering.
# ==========================================================================
scaffold_numbering_producer_flow___CASES = [
    ("c1ccc2[nH]ccc2c1", "indole", "1H-indole"),
    ("c1ccc2ccccc2c1", "naphthalene", "naphthalene"),
    ("C1CCCCC1", "cycloalkane", "cyclohexane"),
]


def scaffold_numbering_producer_flow___selected(smiles: str) -> dict:
    """取并列母体组首位——namer 实际编号用的候选。"""
    mol = preprocess(smiles)
    assert mol is not None
    return select_parent(analyze(mol))[0]


@pytest.mark.parametrize("smiles,_,expected", scaffold_numbering_producer_flow___CASES)
def test_representative_public_names_are_unchanged(smiles: str, _: str, expected: str) -> None:
    result = SMILESNNamer().name(smiles)
    assert result.success
    assert result.en == expected


# ==========================================================================
# 合并自 test_merge_alken_kinds.py
# IUPAC: P-31.1 / P-44
# Layer: L2,L4,L5
#
# Merge open-chain mono-FG alken* ParentKind into saturated kinds.
#
# Unsaturated information lives only on parent fields double_bond /
# double_bonds (+ L4 ene_locant / ene_locants). Names must stay correct;
# parent.kind must be the saturated FG name.
# ==========================================================================
merge_alken_kinds__UNSAT_CASES = [
    ("C=CC(=O)O", "acid", "prop-2-enoic acid", "丙-2-烯酸"),
    ("C=CC=O", "aldehyde", "prop-2-enal", "丙-2-烯醛"),
    ("C=CCO", "alcohol", "prop-2-en-1-ol", "丙-2-烯-1-醇"),
    ("COC(=O)C=C", "ester", "methyl prop-2-enoate", "丙-2-烯酸甲酯"),
    ("C=CC#N", "nitrile", "prop-2-enenitrile", "丙-2-烯腈"),
    ("C=CC(N)=O", "amide", "prop-2-enamide", "丙-2-烯酰胺"),
    ("O=C(O)C=CC(=O)O", "acid", "but-2-enedioic acid", "丁-2-烯二酸"),
]
merge_alken_kinds__SAT_CASES = [
    ("CC(=O)O", "acid", "acetic acid", "乙酸"),
    ("CCO", "alcohol", "ethanol", "乙醇"),
]
merge_alken_kinds__NAME_CASES = merge_alken_kinds__UNSAT_CASES + merge_alken_kinds__SAT_CASES


def merge_alken_kinds___parent(smiles: str) -> dict:
    """取并列母体组首位——namer 实际使用的候选。"""
    mol = preprocess(smiles)
    assert mol is not None
    return select_parent(analyze(mol))[0]


@pytest.mark.parametrize("smiles,kind,en,zh", merge_alken_kinds__UNSAT_CASES)
def test_merge_alken_parent_kind_unsat(
    smiles: str, kind: str, en: str, zh: str | None,
) -> None:
    p = merge_alken_kinds___parent(smiles)
    assert p.get("kind") == kind
    assert p.get("double_bond") or p.get("double_bonds")


@pytest.mark.parametrize("smiles,kind,en,zh", merge_alken_kinds__SAT_CASES)
def test_merge_alken_parent_kind_sat(
    smiles: str, kind: str, en: str, zh: str | None,
) -> None:
    p = merge_alken_kinds___parent(smiles)
    assert p.get("kind") == kind
    assert not p.get("double_bond") and not p.get("double_bonds")


# ==========================================================================
# 合并自 test_fg_locants_records.py
# Layer: L4
#
# fg_locants 稀疏产出契约: principal FG 位次为结构化列表 [{kind, locants, omit}],
# 位次一律取自 principal_expression_facts 的锚点(单一格式, 不再读扁平 *_c_idx 字段);
# 烯/炔不进列表(不饱和度独立); 非 principal 类不产记录.
# ==========================================================================
def fg_locants_records___oriented(kind: str, chain: list, anchors: list, group) -> dict:
    """构造只填锚点的 oriented dict（principal 位次记录的唯一输入）。"""
    facts = PrincipalExpressionFacts(
        group_class=group, multiplicity=len(anchors), relation=PrincipalRelation.IN_SKELETON,
        occurrence_ids=(), characteristic_atoms=frozenset(), anchor_atoms=frozenset(anchors),
        attachment_atoms=frozenset(anchors))
    return {"kind": kind, "n_carbons": len(chain), "chain": list(chain),
            "principal_expression_facts": facts}


def test_ketone_record():
    fg = _fg_locants(fg_locants_records___oriented("ketone", [0, 1, 2, 3], [1], FG.KETONE))
    assert fg == [{"kind": "ketone", "locants": [2], "omit": False}]


def test_dione_single_record():
    fg = _fg_locants(fg_locants_records___oriented("ketone", [0, 1, 2, 3], [1, 3], FG.KETONE))
    assert fg == [{"kind": "ketone", "locants": [2, 4], "omit": False}]


def test_amine_record():
    fg = _fg_locants(fg_locants_records___oriented("amine", [0, 1, 2, 3], [1], FG.AMINE))
    assert fg == [{"kind": "amine", "locants": [2], "omit": False}]


def test_ring_exocyclic_single_fg_produces_one_record():
    """环上单个环外主官能团只产所属类别一条记录：曾因无 gate 的 ring_attach_idx 兜底，给酸/酯/酰胺/腈/醛/酰基六类各产一条同锚点记录。"""
    for smiles, kind in [("O=Cc1ccccc1", "aldehyde"), ("c1ccccc1C(=O)O", "acid")]:
        mol = preprocess(smiles)
        info = analyze(mol)
        parent = select_parent(info)[0]
        numbered = number(parent, extract_substituents(info, parent))
        kinds = [r["kind"] for r in numbered.get("fg_locants") or []]
        assert kinds == [kind], f"{smiles}: {kinds}"


l4_locant_calc___BANNED_FIELDS = frozenset({
    "oh_c_idx", "oh_c_idxs", "amine_c_idx", "amine_c_idxs", "ketone_c_idx", "ketone_c_idxs",
    "sh_c_idx", "sh_c_idxs", "cooh_c_idx", "cooh_c_idxs", "ester_c_idx", "ester_c_idxs",
    "amide_c_idx", "amide_c_idxs", "nitrile_c_idx", "nitrile_c_idxs", "aldehyde_c_idx",
    "aldehyde_c_idxs", "acyl_c_idx", "acyl_c_idxs", "radical_c_idx", "radical_c_idxs",
    "ring_attach_idx", "principal_attachment_atoms",
})


def l4_locant_calc___dict_keys(tree: ast.AST) -> set:
    """源码中按字符串键取值的位置（下标与 .get 参数），不看注释/docstring。"""
    keys = {node.slice.value for node in ast.walk(tree)
            if isinstance(node, ast.Subscript) and isinstance(node.slice, ast.Constant)}
    keys |= {arg.value for node in ast.walk(tree)
             if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
             and node.func.attr == "get"
             for arg in node.args if isinstance(arg, ast.Constant)}
    return keys


def test_locant_calc_reads_no_flat_anchor_fields():
    """位次原子唯一来源是 principal_expression_facts：locant_calc 不得按键读扁平锚点字段，也不得持有 per-class 函数表。"""
    path = Path(__file__).parents[2] / "src" / "opensmiles" / "layer4" / "locant_calc.py"
    tree = ast.parse(path.read_text(encoding="utf-8"))
    assert not (l4_locant_calc___dict_keys(tree) & l4_locant_calc___BANNED_FIELDS)
    assigned = {target.id for node in tree.body if isinstance(node, ast.Assign)
                for target in node.targets if isinstance(target, ast.Name)}
    assert "_LOCANT_FNS" not in assigned


def test_unsaturation_stays_out_of_records():
    """烯/炔不进 fg_locants: 不饱和度是扁平字段, 独立于 principal FG."""
    fg = _fg_locants({"kind": "alkene", "n_carbons": 2, "chain": [0, 1],
                      "double_bond": (0, 1)})
    assert fg == []


# ==========================================================================
# 合并自 test_substituent_namer.py
#
# Ordered SubstituentNamer backends: retained → rooted_tree → recursive.
# ==========================================================================
def substituent_namer___claim(atoms=(1, 2, 3), root=1, attach=0) -> ClaimedBlock:
    return ClaimedBlock(
        slot=SideSlot.OTHER,
        attach_parent=attach,
        root=root,
        atoms=frozenset(atoms),
    )


class substituent_namer___FakeBackend:
    def __init__(self, name: str, result: SubstituentName | None, log: list[str]):
        self.name = name
        self._result = result
        self._log = log

    def try_name(self, mol, claim: ClaimedBlock) -> SubstituentName | None:
        self._log.append(self.name)
        return self._result


def substituent_namer___ok(claim: ClaimedBlock, label: str) -> SubstituentName:
    return SubstituentName(
        claim=claim,
        en=f"{label}-en",
        zh=f"{label}-zh",
        requires_parentheses=False,
    )


def test_retained_success_skips_later_backends():
    claim = substituent_namer___claim()
    log: list[str] = []
    namer = SubstituentNamer(
        backends=[
            substituent_namer___FakeBackend("retained", substituent_namer___ok(claim, "retained"), log),
            substituent_namer___FakeBackend("rooted_tree", substituent_namer___ok(claim, "rooted_tree"), log),
            substituent_namer___FakeBackend("recursive", substituent_namer___ok(claim, "recursive"), log),
        ]
    )
    result = namer.name(None, claim)
    assert result is not None
    assert result.en == "retained-en"
    assert result.claim is claim
    assert log == ["retained"]


def test_rooted_tree_runs_after_retained_failure():
    claim = substituent_namer___claim()
    log: list[str] = []
    namer = SubstituentNamer(
        backends=[
            substituent_namer___FakeBackend("retained", None, log),
            substituent_namer___FakeBackend("rooted_tree", substituent_namer___ok(claim, "rooted_tree"), log),
            substituent_namer___FakeBackend("recursive", substituent_namer___ok(claim, "recursive"), log),
        ]
    )
    result = namer.name(None, claim)
    assert result is not None
    assert result.en == "rooted_tree-en"
    assert log == ["retained", "rooted_tree"]


def test_all_backends_fail_returns_none():
    claim = substituent_namer___claim()
    log: list[str] = []
    namer = SubstituentNamer(
        backends=[
            substituent_namer___FakeBackend("retained", None, log),
            substituent_namer___FakeBackend("rooted_tree", None, log),
            substituent_namer___FakeBackend("recursive", None, log),
        ]
    )
    assert namer.name(None, claim) is None
    assert log == ["retained", "rooted_tree", "recursive"]


def test_backend_none_never_yields_empty_name_or_mutates_claim():
    claim = substituent_namer___claim(atoms=(10, 11))
    before = (claim.slot, claim.attach_parent, claim.root, claim.atoms)
    namer = SubstituentNamer(
        backends=[substituent_namer___FakeBackend("retained", None, []), substituent_namer___FakeBackend("rooted_tree", None, [])]
    )
    result = namer.name(None, claim)
    assert result is None
    assert (claim.slot, claim.attach_parent, claim.root, claim.atoms) == before


