"""Backward-compatible re-export: recursive aryl naming via LeafRegistry.

Implementation lives in layer2/leaves/ (protocol + registry + handlers).
"""
from __future__ import annotations

from namepredict.layer2.leaves.registry import match_leaf_kind as _recurse_leaf_kind
from namepredict.layer2.leaves.ring_namer import recursive_ph_name as _recursive_ph_name

__all__ = ["_recurse_leaf_kind", "_recursive_ph_name"]
