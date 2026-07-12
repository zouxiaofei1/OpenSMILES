"""API tests for POST /api/v1/name and light loop state checks."""

from __future__ import annotations

from fastapi.testclient import TestClient

from server.app import app


def test_name_ethanol():
    c = TestClient(app)
    r = c.post("/api/v1/name", json={"smiles": "CCO"})
    assert r.status_code == 200
    body = r.json()
    assert body["success"] is True
    assert "ethanol" in body["en"].lower()


def test_name_invalid_still_200():
    c = TestClient(app)
    r = c.post("/api/v1/name", json={"smiles": "not-a-smiles!!!"})
    assert r.status_code == 200
    body = r.json()
    assert body["success"] is False


def test_loop_state_200():
    c = TestClient(app)
    r = c.get("/api/v1/loop/state")
    assert r.status_code == 200
    body = r.json()
    assert "status" in body
    assert "iter" in body
