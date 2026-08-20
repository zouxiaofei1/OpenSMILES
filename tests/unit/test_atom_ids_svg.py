"""Tests for atom/SSSR index structure image (server.atom_ids_svg + API)."""

from __future__ import annotations

import re

from fastapi.testclient import TestClient

from server.app import app
from server.atom_ids_svg import build_atom_ids_svg

_TEXT_RE = re.compile(r"<text[^>]*>([^<]+)</text>")


def _labels(smiles: str) -> tuple[list[str], list[str]]:
    res = build_atom_ids_svg(smiles)
    assert res is not None
    texts = _TEXT_RE.findall(res["svg"])
    atom_ids = [t for t in texts if t.isdigit()]
    ring_ids = [t for t in texts if t.startswith("R")]
    return atom_ids, ring_ids


def test_chain_atom_ids_ethanol():
    atom_ids, ring_ids = _labels("CCO")
    assert atom_ids == ["0", "1", "2"]
    assert ring_ids == []


def test_ring_atom_ids_and_ring_id_cyclohexane():
    atom_ids, ring_ids = _labels("C1CCCCC1")
    assert atom_ids == ["0", "1", "2", "3", "4", "5"]
    assert ring_ids == ["R1"]


def test_fused_two_rings_naphthalene():
    atom_ids, ring_ids = _labels("c1ccc2ccccc2c1")
    assert atom_ids == ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9"]
    assert ring_ids == ["R1", "R2"]


def test_invalid_smiles_returns_none():
    assert build_atom_ids_svg("not-a-smiles!!") is None


def test_api_endpoint_ok():
    c = TestClient(app)
    r = c.post("/api/v1/name/atom-ids-svg", json={"smiles": "c1ccc2ccccc2c1"})
    assert r.status_code == 200
    body = r.json()
    assert body["ok"] is True
    assert "R1" in body["svg"] and "R2" in body["svg"]


def test_api_endpoint_fail():
    c = TestClient(app)
    r = c.post("/api/v1/name/atom-ids-svg", json={"smiles": "not-a-smiles!!"})
    assert r.status_code == 200
    assert r.json()["ok"] is False
