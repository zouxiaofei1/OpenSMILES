"""Fused56 scaffold builder: engine Spec constants + ScaffoldSpec id lookup.

Thin ring modules import Mono/Di13 engine specs from here.
Chain orientation stays in layer2.fused56; labels authority is specs.FUSED56_LABELS.
"""
from __future__ import annotations

from namepredict.layer2.scaffold.fused56 import Fused56Di13Spec, Fused56MonoSpec
from namepredict.layer2.scaffold.specs import FUSED56_LABELS, get_spec

__all__ = [
    "FUSED56_LABELS",
    "FUSED56_AZA_KINDS",
    "BF_MONO",
    "BF_AMINE",
    "BT_MONO",
    "BT_OL",
    "BTZ_DI13",
    "BOX_DI13",
    "scaffold_id_for_kind",
]

# Aza fused56 kinds (indole / indazole / bim); O/S engine specs below.
FUSED56_AZA_KINDS = frozenset({
    "indole", "indolecarboxylic", "indazole", "indazolecarbonitrile",
    "indazolecarbaldehyde", "benzimidazole", "benzimidazolamine",
})

# Engine specs aligned with ScaffoldSpec ids in specs.FUSED56_SPECS.
BF_MONO = Fused56MonoSpec(
    kind="benzofuran", hetero_z=8, hetero_key="o_idx", sub_cap=1,
)
BF_AMINE = Fused56MonoSpec(
    kind="benzofuranamine", hetero_z=8, hetero_key="o_idx", sub_cap=1,
)
BT_MONO = Fused56MonoSpec(
    kind="benzothiophene", hetero_z=16, hetero_key="s_idx", sub_cap=1,
)
BT_OL = Fused56MonoSpec(
    kind="benzothiophenol", hetero_z=16, hetero_key="s_idx", sub_cap=1,
)
BTZ_DI13 = Fused56Di13Spec(
    kind="benzothiazole", hetero_z=16, hetero_key="s_idx",
    amine_kind="benzothiazolamine",
)
BOX_DI13 = Fused56Di13Spec(
    kind="benzoxazole", hetero_z=8, hetero_key="o_idx",
    amine_kind="benzoxazolamine",
)


def scaffold_id_for_kind(kind: str | None) -> str | None:
    """Return registered fused56 scaffold id for kind, else None."""
    if not kind:
        return None
    sp = get_spec(kind)
    return sp.id if sp is not None and sp.naming_class == "fused56" else None
