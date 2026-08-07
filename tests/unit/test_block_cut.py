from namepredict.layer0.preprocessor import preprocess
from namepredict.layer1.analyzer import analyze
from namepredict.layer2.parent_selector import select_parent
from namepredict.tools.block_cut import parent_atom_set, cut_block, side_roots
from namepredict.layer3.submol_build import build_cut_submol


def test_cut_n_phenyl_from_benzamide():
    mol = preprocess("c1ccccc1C(=O)Nc2ccccc2")
    parent = select_parent(analyze(mol))
    assert parent["kind"] == "benzamide"
    patoms = parent_atom_set(parent, mol)
    # amide N is in parent; N-phenyl ring is outside
    am = analyze(mol)["amides"][0]
    n_idx = am["n_idx"]
    roots = [r for r in side_roots(mol, patoms) if r in set(am.get("n_c_idxs") or [])]
    assert len(roots) == 1
    block = cut_block(mol, roots[0], patoms)
    assert block is not None and len(block) == 6
    sub = build_cut_submol(mol, block, roots[0])
    assert sub is not None
    assert sub.mol.GetNumAtoms() >= 6
    assert sub.attach_new in sub.atom_map
