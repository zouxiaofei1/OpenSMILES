# tests/unit/test_settings_api.py
from pathlib import Path
from fastapi.testclient import TestClient
import server.deps as deps
from agent_loop.secrets import SecretsStore
from server.app import app

def test_api_key_roundtrip(tmp_path: Path, monkeypatch):
    store = SecretsStore(tmp_path / "secrets.json")
    monkeypatch.setattr(deps, "get_secrets", lambda: store)
    c = TestClient(app)
    r = c.get("/api/v1/settings/api-key")
    assert r.status_code == 200
    assert r.json()["configured"] is False
    r2 = c.post("/api/v1/settings/api-key", json={"api_key": "sk-ant-secret-9999"})
    assert r2.status_code == 200
    body = r2.json()
    assert body["configured"] is True
    assert body["hint"] == "9999"
    assert "secret" not in body.get("hint", "")
    # GET must not leak full key
    g = c.get("/api/v1/settings/api-key").json()
    assert "sk-ant-secret" not in str(g)

def test_start_without_key_400(tmp_path: Path, monkeypatch):
    store = SecretsStore(tmp_path / "secrets.json")
    monkeypatch.setattr(deps, "get_secrets", lambda: store)
    # ensure controller uses same store gate
    c = TestClient(app)
    r = c.post("/api/v1/loop/start")
    assert r.status_code == 400
    assert r.json()["error"] == "missing_api_key"
