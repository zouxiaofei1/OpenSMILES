"""Cluster benchmark failures by features (or source+error_langs fallback)."""

from __future__ import annotations

from collections import defaultdict
from typing import Any, TypedDict

# Feature → suggested naming pipeline layers (L0–L5).
_FEATURE_LAYERS: dict[str, tuple[str, ...]] = {
    "alkane": ("L2", "L5"),
    "alkene": ("L1", "L2", "L5"),
    "alkyne": ("L1", "L2", "L5"),
    "alcohol": ("L1", "L3", "L5"),
    "phenol": ("L1", "L3", "L5"),
    "ether": ("L1", "L3", "L5"),
    "amine": ("L1", "L3", "L5"),
    "primary_amine": ("L1", "L3", "L5"),
    "secondary_amine": ("L1", "L3", "L5"),
    "ketone": ("L1", "L3", "L5"),
    "aldehyde": ("L1", "L3", "L5"),
    "carboxylic_acid": ("L1", "L3", "L5"),
    "ester": ("L1", "L3", "L5"),
    "amide": ("L1", "L3", "L5"),
    "nitrile": ("L1", "L3", "L5"),
    "nitro": ("L1", "L3", "L5"),
    "halide": ("L1", "L3", "L5"),
    "halogen": ("L1", "L3", "L5"),
    "chloride": ("L1", "L3", "L5"),
    "bromide": ("L1", "L3", "L5"),
    "fluorine": ("L1", "L3", "L5"),
    "chloro": ("L1", "L3", "L5"),
    "bromo": ("L1", "L3", "L5"),
    "aromatic": ("L1", "L2", "L4"),
    "benzene": ("L1", "L2", "L4"),
    "phenyl": ("L1", "L3", "L4"),
    "heterocycle": ("L1", "L2", "L4"),
    "pyridine": ("L1", "L2", "L4"),
    "fused_ring": ("L1", "L2", "L4"),
    "locant": ("L4", "L5"),
    "retained_name": ("L5",),
    "single_carbon": ("L2", "L5"),
}

_DEFAULT_LAYERS: tuple[str, ...] = ("L1", "L5")
_SAMPLE_CAP = 8


class Cluster(TypedDict):
    key: str
    size: int
    features: list[str]
    sample_smiles: list[str]
    suggest_layers: list[str]
    error_langs: list[str]


def _error_langs(fail: dict[str, Any]) -> list[str]:
    langs: list[str] = []
    if fail.get("en_ok") is False:
        langs.append("en")
    if fail.get("zh_ok") is False:
        langs.append("zh")
    return langs or ["unknown"]


def _norm_features(fail: dict[str, Any]) -> tuple[str, ...]:
    feats = fail.get("features") or []
    return tuple(sorted({str(f) for f in feats if f}))


def _cluster_key(fail: dict[str, Any]) -> tuple:
    feats = _norm_features(fail)
    if feats:
        return ("features", feats)
    langs = tuple(_error_langs(fail))
    source = str(fail.get("source") or "unknown")
    return ("fallback", source, langs)


def _suggest_layers(features: list[str]) -> list[str]:
    if not features:
        return list(_DEFAULT_LAYERS)
    seen: list[str] = []
    for feat in features:
        for layer in _FEATURE_LAYERS.get(feat, _DEFAULT_LAYERS):
            if layer not in seen:
                seen.append(layer)
    return seen


def _key_str(key: tuple) -> str:
    kind = key[0]
    if kind == "features":
        return "features:" + "+".join(key[1])
    return f"fallback:{key[1]}:{'+'.join(key[2])}"


def _group_features(key: tuple) -> list[str]:
    return list(key[1]) if key[0] == "features" else []


def _group_error_langs(group: list[dict[str, Any]]) -> list[str]:
    langs: list[str] = []
    for f in group:
        for lang in _error_langs(f):
            if lang not in langs and lang != "unknown":
                langs.append(lang)
    return langs or ["unknown"]


def _sample_smiles(group: list[dict[str, Any]]) -> list[str]:
    return [str(f.get("smiles") or "") for f in group[:_SAMPLE_CAP]]


def _build_cluster(key: tuple, group: list[dict[str, Any]]) -> Cluster:
    feats = _group_features(key)
    return {
        "key": _key_str(key),
        "size": len(group),
        "features": feats,
        "sample_smiles": _sample_smiles(group),
        "suggest_layers": _suggest_layers(feats),
        "error_langs": _group_error_langs(group),
    }


def cluster_failures(fails: list[dict], top_k: int = 5) -> list[Cluster]:
    """Group fail records by features (or source+error_langs); return top_k by size."""
    groups: dict[tuple, list[dict[str, Any]]] = defaultdict(list)
    for fail in fails:
        groups[_cluster_key(fail)].append(fail)
    ordered = sorted(groups.items(), key=lambda kv: (-len(kv[1]), _key_str(kv[0])))
    return [_build_cluster(k, g) for k, g in ordered[:top_k]]
