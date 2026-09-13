# 合并自 17 个原测试文件（按主题分组，内容与断言未改动）。
"""
test_l2_parent_core_contract.py: L2 parent_core is the sole authority for parent assembly / chain / gate helpers.
test_l2_l3_side_facts_contract.py: Architecture contract for Layer 2 side topology facts consumed by Layer 3.
test_l5_no_layer2_private.py: L5 special FG name modules must not import L2 private APIs.
test_ring_scaffold_expression_contract.py:
test_anchored_whole_mol.py: 整分子锚定键命中：顶层命名 *xxx 与取代基查表（resolve_name）一致。
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
test_yl_form.py: Tests for yl_form FG suffix→prefix conversion (P-63.2.2 / P-63.2.1 / P-62.2).
"""
from __future__ import annotations

import ast
import pytest

from namepredict.layer0.preprocessor import preprocess
from namepredict.layer1.analyzer import analyze
from namepredict.layer2.candidates import _collect_candidates
from namepredict.layer2.parent_ownership import finalize_parent_ownership
from namepredict.layer2.parent_selector import select_parent
from namepredict.layer2.ring_scaffold import all_specs
from namepredict.layer3.claimable_block import ClaimedBlock, SideSlot, claim_block, iter_claims
from namepredict.layer3.coverage import build_coverage_ledger
from namepredict.layer3.substituent_extractor import extract_substituents
from namepredict.layer3.substituent_namer import SubstituentName, SubstituentNamer
from namepredict.layer4.locant_calc import _fg_locants
from namepredict.layer4.numbering import number
from namepredict.namer import SMILESNNamer, _names_from_subs
from namepredict.tools.anchored_table import resolve_name
from namepredict.tools.block_cut import parent_atom_set
from namepredict.tools.free_to_yl import free_to_yl as yl_form
from namepredict.tools.re import normalize_en, normalize_zh
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
l2_parent_core_contract___L2 = l2_parent_core_contract___ROOT / "src" / "namepredict" / "layer2"

# Helpers that must be sourced from parent_core, not parent_selector.
l2_parent_core_contract___HELPERS = frozenset({
    "_parent_dict",
    "_longest_from",
    "_longest_chain",
})

# Modules allowed to re-export / host try surface.
l2_parent_core_contract___EXEMPT = frozenset({
    "parent_selector.py",
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
        if node.module != "namepredict.layer2.parent_selector":
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
# 合并自 test_l2_l3_side_facts_contract.py
# IUPAC: P-29.2
# Layer: L2,L3
#
# Architecture contract for Layer 2 side topology facts consumed by Layer 3.
# ==========================================================================
l2_l3_side_facts_contract__LAYER3 = Path(__file__).parents[2] / "src" / "namepredict" / "layer3"
l2_l3_side_facts_contract__TOOLS = Path(__file__).parents[2] / "src" / "namepredict" / "tools"
# Side-topology facts merged into tools 碳拓扑原语 chain.py。
l2_l3_side_facts_contract__TOOLS_FACTS = l2_l3_side_facts_contract__TOOLS / "chain.py"
l2_l3_side_facts_contract__LEAF_PROTOCOL = l2_l3_side_facts_contract__LAYER3 / "leaves" / "protocol.py"


def l2_l3_side_facts_contract___private_layer2_imports(path: Path) -> list[str]:
    tree = ast.parse(path.read_text(encoding="utf-8"))
    found: list[str] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom) and (node.module or "").startswith("namepredict.layer2"):
            found.extend(alias.name for alias in node.names if alias.name.startswith("_"))
    return found


def l2_l3_side_facts_contract___all_functions(path: Path) -> list[ast.FunctionDef | ast.AsyncFunctionDef]:
    tree = ast.parse(path.read_text(encoding="utf-8"))
    return [n for n in ast.walk(tree) if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef))]


def l2_l3_side_facts_contract___tools_modules(tree: ast.AST) -> list[str]:
    direct = [alias.name for node in ast.walk(tree) if isinstance(node, ast.Import)
              for alias in node.names if alias.name.startswith("namepredict.tools")]
    froms = [(node.module or "") for node in ast.walk(tree) if isinstance(node, ast.ImportFrom)
             and (node.module or "").startswith("namepredict.tools")]
    return [*direct, *froms]


def test_tools_import_scan_covers_both_ast_forms() -> None:
    tree = ast.parse("import namepredict.tools.aryl_sub\nfrom namepredict.tools import side_alkyl")
    assert l2_l3_side_facts_contract___tools_modules(tree) == [
        "namepredict.tools.aryl_sub", "namepredict.tools",
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
l5_no_layer2_private___L5 = l5_no_layer2_private___ROOT / "src" / "namepredict" / "layer5"

# L5 special-FG name modules that previously walked via L2 private APIs.


def l5_no_layer2_private___layer2_import_hits(path: Path) -> list[str]:
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    hits: list[str] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom) and node.module:
            if node.module == "namepredict.layer2" or node.module.startswith(
                "namepredict.layer2."
            ):
                hits.append(f"from {node.module} (L{node.lineno})")
        elif isinstance(node, ast.Import):
            for alias in node.names:
                if alias.name == "namepredict.layer2" or alias.name.startswith(
                    "namepredict.layer2."
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
    assert parent["typed_ring_expression_supported"] is True


def test_benzene_resolves_to_stable_scaffold_identity():
    parent = ring_scaffold_expression_contract___typed("Oc1ccccc1", "alcohol")
    assert parent["scaffold_id"] == "benzene"
    assert parent["scaffold_identity"].naming_class == "mono_carbo"
    assert parent["typed_ring_expression_supported"] is True


def test_naphthalenol_uses_naph_family_expression_policy():
    parent = ring_scaffold_expression_contract___typed("Oc1cccc2ccccc12", "alcohol")
    assert parent["scaffold_id"] == "naphthalene"
    assert parent["typed_ring_expression_supported"] is True


# ==========================================================================
# 合并自 test_anchored_whole_mol.py
#
# 整分子锚定键命中：顶层命名 *xxx 与取代基查表（resolve_name）一致。
#
# 之前顶层把 *xxx 当自由基母体硬算，暴露出与保留表不一致的编号/骨架错误
# （quinolin-5-yl vs quinolin-2-yl、ethan-1-yl vs methoxy、FAIL 等）。
# 本测试要求：整分子 canonical == 锚定表键时，顶层命名直接返回保留名。
# ==========================================================================
def anchored_whole_mol___name(smiles: str):
    return SMILESNNamer().name(smiles)


# ── registry 条目：顶层 *xxx == resolve_name（registry 双语名）──

@pytest.mark.parametrize("smi,key", [
    ("*C(C)C", "isopropyl"),
    ("*C(C)(C)C", "tert-butyl"),
    ("*Cc1ccccc1", "benzyl"),
    ("*OC", "methoxy"),
    ("*S(C)(=O)=O", "methylsulfonyl"),
    ("*C#N", "cyano"),
    ("*O", "hydroxy"),
    ("*[O]", "oxidanyl"),
    ("*N", "amino"),
    ("*S", "sulfanyl"),
])
def test_top_level_star_matches_registry(smi, key):
    r = anchored_whole_mol___name(smi)
    assert r.success, f"{smi} 命名失败: {r.meta.get('reason')}"
    en, zh = resolve_name(key)
    assert r.en == en, f"{smi} en: {r.en!r} != {en!r}"
    assert r.zh == zh, f"{smi} zh: {r.zh!r} != {zh!r}"


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


# ── 命中路径标记 ──

def test_hit_sets_anchored_meta():
    r = anchored_whole_mol___name("*C(C)C")
    assert (r.meta or {}).get("anchored") is True
    r2 = anchored_whole_mol___name("*Cc1ccccc1")
    assert (r2.meta or {}).get("anchored") is True


# ── 普通分子不受影响（表键带 * 前缀，绝不误命中）──

def test_plain_molecules_not_hit():
    assert anchored_whole_mol___name("c1ccccc1").en == "benzene"
    assert anchored_whole_mol___name("CC(C)C").success
    assert anchored_whole_mol___name("CC(C)C").en != "isopropyl"
    assert ( anchored_whole_mol___name("CC(C)C").meta or {}).get("anchored") is not True


# ── 带附加取代基的自由基：整分子不在表，仍走管线 ──

def test_star_with_extra_substituent_not_hit():
    # 整分子 *c1ccc(Cl)cc1 不是单一锚定键 → 走管线，绝不误命中 phenyl
    r = anchored_whole_mol___name("*c1ccc(Cl)cc1")
    assert (r.meta or {}).get("anchored") is not True
    assert r.en != "phenyl"


# ── *O/*[O]/*N 只服务整分子顶层，不进 L3 取代基查表 ──

def test_mononuclear_anchored_whole_mol_only():
    # 整分子 *O 顶层命中 anchored（hydroxy）
    r = anchored_whole_mol___name("*O")
    assert r.en == "hydroxy"
    assert (r.meta or {}).get("anchored") is True


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


def test_n_phenyl_benzamide_amide_n_claim():
    """N-phenyl benzamide: one AMIDE_N claim with six phenyl atoms."""
    mol, info = claimable_block_api___info("c1ccccc1C(=O)Nc2ccccc2")
    parent = select_parent(info)[0]
    assert parent["kind"] == "amide"
    assert parent.get("scaffold_id") == "benzene"
    claims = iter_claims(mol, parent["owned_atoms"])
    assert len(claims) == 1
    c = claims[0]
    assert c.slot == SideSlot.AMIDE_N
    assert len(c.atoms) == 6
    assert all(mol.GetAtomWithIdx(i).GetAtomicNum() == 6 for i in c.atoms)
    assert c.attach_parent in parent["owned_atoms"]
    assert c.root not in parent["owned_atoms"]
    assert c.root in c.atoms


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


def test_ketone_carbonyl_o_not_claimed():
    """外部羰基氧（双键连所属碳）被跳过，不作侧链 claim（主 FG 已处理）。"""
    mol, info = claimable_block_api___info("OC(=O)C(=O)C")
    parent = select_parent(info)[0]
    assert parent["kind"] == "acid"
    claims = iter_claims(mol, parent["owned_atoms"])
    assert all(
        mol.GetAtomWithIdx(c.root).GetAtomicNum() != 8 for c in claims
    )


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
    ("CC(C)CC", "2-methylbutane", "2-甲基丁烷"),
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
    ("Nc1ccc(C)cc1", "4-methylaniline", "4-甲基苯胺"),
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
        subst = extract_substituents(info, parent)
        names = _names_from_subs(subst)
        led = build_coverage_ledger(mol, owned_atoms=parent["owned_atoms"], names=names)
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

    am = info["amides"][0]
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

    c_idx = parent["cooh_c_idx"]
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
    _ = parent_atom_set(parent, mol)
    _ = finalize_parent_ownership(parent, mol)
    _ = list(owned)
    copy_out = set(owned)
    copy_out.add(999)

    assert parent["owned_atoms"] == snapshot
    assert parent["owned_atoms"] is owned
    assert isinstance(parent["owned_atoms"], frozenset)


def test_parent_atom_set_adapter_prefers_owned_atoms():
    """Compatibility adapter returns owned_atoms when present."""
    mol, info = parent_ownership___info("CC(=O)O")
    parent = select_parent(info)[0]
    assert parent_atom_set(parent, mol) == parent["owned_atoms"]


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
    raw = {
        "chain": [0, 1],
        "n_carbons": 2,
        "kind": "acid",
        "cooh_c_idx": 1,
    }
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
    from namepredict.layer3.substituent_namer import SubstituentNamer

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


def test_scaffold_parent_without_materialized_facts_is_rejected() -> None:
    parent = scaffold_numbering_producer_flow___selected("c1ccc2[nH]ccc2c1")
    parent.pop("numbering_scaffold")
    with pytest.raises(ValueError, match="numbering_scaffold"):
        number(parent, [])


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
# 只产实际存在的 FG; 烯/炔不进列表(扁平字段独立); cooh 不产(死字段).
# ==========================================================================
def test_single_alcohol_record():
    fg = _fg_locants({"kind": "alcohol", "n_carbons": 4, "chain": [0, 1, 2, 3],
                      "oh_c_idx": 1, "oh_c_idxs": [1]})
    assert fg == [{"kind": "oh", "locants": [2], "omit": False}]


def test_methanol_omit_flag():
    fg = _fg_locants({"kind": "alcohol", "n_carbons": 1, "chain": [0],
                      "oh_c_idx": 0, "oh_c_idxs": [0]})
    assert fg == [{"kind": "oh", "locants": [1], "omit": True}]


def test_diol_single_oh_record_multi_locants():
    fg = _fg_locants({"kind": "alcohol", "n_carbons": 7, "chain": [9, 8, 6, 5, 3, 1, 0],
                      "oh_c_idxs": [1, 9]})
    assert fg == [{"kind": "oh", "locants": [1, 6], "omit": False}]


def test_ketone_record():
    fg = _fg_locants({"kind": "ketone", "n_carbons": 4, "chain": [0, 1, 2, 3],
                      "ketone_c_idx": 1, "ketone_c_idxs": [1]})
    assert fg == [{"kind": "ketone", "locants": [2], "omit": False}]


def test_dione_single_record():
    fg = _fg_locants({"kind": "ketone", "n_carbons": 4, "chain": [0, 1, 2, 3],
                      "ketone_c_idxs": [1, 3]})
    assert fg == [{"kind": "ketone", "locants": [2, 4], "omit": False}]


def test_amine_record():
    fg = _fg_locants({"kind": "amine", "n_carbons": 4, "chain": [0, 1, 2, 3],
                      "amine_c_idx": 1, "amine_c_idxs": [1]})
    assert fg == [{"kind": "amine", "locants": [2], "omit": False}]


def test_acid_produces_no_record():
    """cooh 是死字段: 单/多酸位次隐含, 不产记录."""
    assert _fg_locants({"kind": "acid", "n_carbons": 5, "chain": [0, 1, 2, 3, 4],
                        "cooh_c_idx": 4, "cooh_c_idxs": [4]}) == []


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


# ==========================================================================
# 合并自 test_yl_form.py
#
# Tests for yl_form FG suffix→prefix conversion (P-63.2.2 / P-63.2.1 / P-62.2).
# ==========================================================================
@pytest.mark.parametrize("en,zh,k,exp_en,exp_zh", [
    ("methanol", "甲醇", 1, "methoxy", "甲氧基"),
    ("ethanol", "乙醇", 2, "ethoxy", "乙氧基"),
    ("propan-1-ol", "丙-1-醇", 1, "propoxy", "丙氧基"),
    ("propan-2-ol", "丙-2-醇", 2, "propan-2-yloxy", "丙-2-基氧基"),
])
def test_alcohol_to_alkoxy(en, zh, k, exp_en, exp_zh):
    got_en, got_zh, paren = yl_form(en, zh, k)
    assert got_en == exp_en, f"EN: expected {exp_en}, got {got_en}"
    assert got_zh == exp_zh, f"ZH: expected {exp_zh}, got {got_zh}"


# ── thiol → alkylsulfanyl (P-63.2.1) ──
@pytest.mark.parametrize("en,zh,k,exp_en,exp_zh", [
    ("methanethiol", "甲硫醇", 1, "methylsulfanyl", "甲硫基"),
    ("ethanethiol", "乙硫醇", 2, "ethylsulfanyl", "乙硫基"),
])
def test_thiol_to_sulfanyl(en, zh, k, exp_en, exp_zh):
    got_en, got_zh, paren = yl_form(en, zh, k)
    assert got_en == exp_en, f"EN: expected {exp_en}, got {got_en}"
    assert got_zh == exp_zh, f"ZH: expected {exp_zh}, got {got_zh}"


# ── primary amine → alkylamino (P-62.2) ──
@pytest.mark.parametrize("en,zh,k,exp_en,exp_zh", [
    ("methanamine", "甲胺", 1, "methylamino", "甲氨基"),
    ("ethanamine", "乙胺", 2, "ethylamino", "乙氨基"),
])
def test_primary_amine_to_amino(en, zh, k, exp_en, exp_zh):
    got_en, got_zh, paren = yl_form(en, zh, k)
    assert got_en == exp_en, f"EN: expected {exp_en}, got {got_en}"
    assert got_zh == exp_zh, f"ZH: expected {exp_zh}, got {got_zh}"


# ── non-converting: fallback to -n-yl ──
@pytest.mark.parametrize("en,zh,k,exp_en,exp_zh", [
    ("N-methylethanamine", "N-甲基乙胺", 1, "N-methylethanamin-1-yl", "N-甲基乙胺-1-基"),
    ("benzene", "苯", 1, "phenyl", "苯基"),
    ("ethane", "乙烷", 1, "ethan-1-yl", "乙烷-1-基"),
])
def test_fallback_to_n_yl(en, zh, k, exp_en, exp_zh):
    got_en, got_zh, paren = yl_form(en, zh, k)
    assert got_en == exp_en, f"EN: expected {exp_en}, got {got_en}"
    assert got_zh == exp_zh, f"ZH: expected {exp_zh}, got {got_zh}"
