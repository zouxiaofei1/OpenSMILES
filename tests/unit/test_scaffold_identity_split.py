# IUPAC: P-22 / P-25 scaffold identity contract
# Layer: L2
from namepredict.layer2.scaffold.identity import ScaffoldIdentity
from namepredict.layer2.scaffold.specs import all_identities, all_specs, get_identity, get_spec


def test_scaffold_identity_is_unique_and_stable():
    identities = all_identities()
    assert identities
    assert len({identity.id for identity in identities}) == len(identities)
    assert all(isinstance(identity, ScaffoldIdentity) for identity in identities)


def test_scaffold_spec_compatibility_facade_matches_identity():
    for spec in all_specs():
        identity = get_identity(spec.id)
        assert identity == spec.identity
        assert get_spec(identity.id) is spec
        assert (identity.naming_class, identity.n_rings, identity.ring) == (
            spec.naming_class, spec.n_rings, spec.ring,
        )
