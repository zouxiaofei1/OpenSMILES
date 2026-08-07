from namepredict.layer0.preprocessor import preprocess
from namepredict.layer1.analyzer import analyze
from namepredict.layer2.parent_selector import select_parent
from namepredict.tools.block_cut import parent_atom_set, cut_block, side_roots
from namepredict.layer3.as_substituent import name_as_substituent


def test_name_thiazol_2_yl_from_benzamide_context():
    mol = preprocess("c1ccccc1C(=O)Nc2nccs2")
    parent = select_parent(analyze(mol))
    patoms = parent_atom_set(parent, mol)
    am = analyze(mol)["amides"][0]
    roots = [r for r in side_roots(mol, patoms) if r in set(am.get("n_c_idxs") or [])]
    assert len(roots) == 1
    atoms = cut_block(mol, roots[0], patoms)
    assert atoms is not None
    attach_old = roots[0]
    got = name_as_substituent(mol, attach_old, atoms, depth=0)
    assert got is not None
    en, zh, paren = got
    assert en == "1,3-thiazol-2-yl"
