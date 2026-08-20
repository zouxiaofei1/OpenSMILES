"""Tests for L4 locant structure image (server.locants_svg + API endpoint)."""

from __future__ import annotations

import re

from fastapi.testclient import TestClient

from server.app import app
from server.locants_svg import build_locants_svg

_TEXT_RE = re.compile(r"<text[^>]*>([^<]+)</text>")


def _texts(smiles: str) -> list[str]:
    res = build_locants_svg(smiles)
    assert res is not None
    return _TEXT_RE.findall(res["svg"])


def test_chain_locants_ethanol():
    assert _texts("CCO") == ["1", "2"]


def test_ring_locants_cyclohexane():
    assert _texts("C1CCCCC1") == ["1", "2", "3", "4", "5", "6"]


def test_branched_chain_locants():
    assert _texts("CCC(CC)CC") == ["1", "2", "3", "4", "5"]


def test_fused_letter_locants_quinoline():
    texts = _texts("c1ccc2ncccc2c1")
    assert "4a" in texts and "8a" in texts


def test_invalid_smiles_returns_none():
    assert build_locants_svg("not-a-smiles!!") is None


def test_api_endpoint_ok():
    c = TestClient(app)
    r = c.post("/api/v1/name/locants-svg", json={"smiles": "c1ccc2ncccc2c1"})
    assert r.status_code == 200
    body = r.json()
    assert body["ok"] is True
    assert "4a" in body["svg"]


def test_api_endpoint_fail():
    c = TestClient(app)
    r = c.post("/api/v1/name/locants-svg", json={"smiles": "not-a-smiles!!"})
    assert r.status_code == 200
    assert r.json()["ok"] is False
