from __future__ import annotations
from namepredict.layer4.locants.adapt import plan_from_chain
from namepredict.layer4.numbering_engine import orient_numbering
from namepredict.layer4.orienters import (
    _kind_orienters, _orient_benzoic, _orient_by_kind, _orient_chain, _orient_polycarboxylic,
)  # noqa: F401  (re-export for tests)
from namepredict.layer4.locant_calc import _amine_locant, _oh_locant, _pack, _sub_locant, _with_locants  # noqa: F401  (re-export for tests)
def number(parent: dict, substituents: list) -> dict:
    chain = orient_numbering(parent, substituents)
    kind = parent.get("kind")
    if parent.get("numbering_scaffold_required") and not parent.get("numbering_scaffold"):
        raise ValueError("numbering_scaffold facts required for selected scaffold")
    oriented = {**parent, "chain": chain}
    plan = plan_from_chain(chain, oriented.get("scaffold_id") or kind, oriented.get("numbering_scaffold"), required=oriented.get("numbering_scaffold_required", False))
    if plan is not None: oriented["numbering"] = plan
    return _pack(oriented, _with_locants(chain, substituents, kind, oriented.get("numbering_scaffold")))
