# IUPAC: P-22 / P-44 ring scaffold-expression contract
# Layer: L2
from rdkit import Chem

from namepredict.layer1.analyzer import analyze
from namepredict.layer2.candidates import _collect_candidates


def _typed(smiles, group):
    parents = _collect_candidates(analyze(Chem.MolFromSmiles(smiles)))
    return next(p for p in parents if (f := p.get("principal_expression_facts")) and f.group_class.value == group)


def test_generic_carbocycle_resolves_independently_of_alcohol_expression():
    parent = _typed("OC1CCCCC1", "alcohol")
    assert parent["scaffold_identity"].naming_class == "carbocycle"
    assert parent["typed_ring_expression_supported"] is True


def test_benzene_resolves_to_stable_scaffold_identity():
    parent = _typed("Oc1ccccc1", "alcohol")
    assert parent["scaffold_id"] == "benzene"
    assert parent["scaffold_identity"].naming_class == "mono_carbo"
    assert parent["typed_ring_expression_supported"] is True


def test_naphthalenol_uses_naph_family_expression_policy():
    parent = _typed("Oc1cccc2ccccc12", "alcohol")
    assert parent["scaffold_id"] == "naphthalene"
    assert parent["typed_ring_expression_supported"] is True


def test_unmigrated_naphthalene_polyol_is_explicitly_not_supported():
    parent = _typed("Oc1ccc2ccccc2c1O", "alcohol")
    assert parent["scaffold_id"] == "naphthalene"
    assert parent["typed_ring_expression_supported"] is False
