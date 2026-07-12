"""Unit tests for tools.fail_cluster — cluster by feature aggregation."""

from __future__ import annotations

from tools.fail_cluster import cluster_failures


def _fail(
    smiles: str,
    features: list[str] | None = None,
    *,
    source: str = "smiles_tiers",
    en_ok: bool | None = False,
    zh_ok: bool | None = False,
) -> dict:
    return {
        "smiles": smiles,
        "features": list(features or []),
        "source": source,
        "en_ok": en_ok,
        "zh_ok": zh_ok,
        "dual_ok": False,
    }


def test_cluster_by_features_aggregation():
    fails = [
        _fail("CCO", ["alcohol"], en_ok=False, zh_ok=True),
        _fail("CCCO", ["alcohol"], en_ok=False, zh_ok=True),
        _fail("C", ["alkane"], en_ok=False, zh_ok=False),
    ]
    clusters = cluster_failures(fails, top_k=5)
    assert len(clusters) == 2
    top = clusters[0]
    assert top["size"] == 2
    assert top["features"] == ["alcohol"]
    assert set(top["sample_smiles"]) == {"CCO", "CCCO"}
    assert "en" in top["error_langs"]
    assert "zh" not in top["error_langs"] or top["error_langs"] == ["en"]


def test_empty_features_uses_source_and_error_langs():
    fails = [
        _fail("CCO", [], source="chebi", en_ok=False, zh_ok=None),
        _fail("CCN", [], source="chebi", en_ok=False, zh_ok=None),
        _fail("CCC", [], source="smiles_tiers", en_ok=True, zh_ok=False),
    ]
    clusters = cluster_failures(fails, top_k=5)
    assert len(clusters) == 2
    sizes = sorted(c["size"] for c in clusters)
    assert sizes == [1, 2]
    big = next(c for c in clusters if c["size"] == 2)
    assert big["features"] == []
    assert "chebi" in big["key"]
    assert big["error_langs"] == ["en"]


def test_top_k_limits_and_orders_by_size():
    fails = [
        _fail("A1", ["alcohol"]),
        _fail("A2", ["alcohol"]),
        _fail("A3", ["alcohol"]),
        _fail("K1", ["ketone"]),
        _fail("K2", ["ketone"]),
        _fail("N1", ["nitrile"]),
    ]
    clusters = cluster_failures(fails, top_k=2)
    assert len(clusters) == 2
    assert clusters[0]["size"] == 3
    assert clusters[0]["features"] == ["alcohol"]
    assert clusters[1]["size"] == 2
    assert clusters[1]["features"] == ["ketone"]


def test_features_order_normalized_same_cluster():
    fails = [
        _fail("X", ["alcohol", "alkene"]),
        _fail("Y", ["alkene", "alcohol"]),
    ]
    clusters = cluster_failures(fails, top_k=5)
    assert len(clusters) == 1
    assert clusters[0]["size"] == 2
    assert clusters[0]["features"] == ["alcohol", "alkene"]


def test_cluster_fields_present():
    fails = [_fail("CCO", ["alcohol"], en_ok=False, zh_ok=False)]
    c = cluster_failures(fails)[0]
    for k in ("key", "size", "features", "sample_smiles", "suggest_layers", "error_langs"):
        assert k in c
    assert c["size"] == 1
    assert isinstance(c["suggest_layers"], list)
    assert c["suggest_layers"]
    assert set(c["error_langs"]) == {"en", "zh"}


def test_empty_fails_returns_empty():
    assert cluster_failures([]) == []
