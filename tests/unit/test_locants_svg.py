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


# ── preferred_orientation 水平行开关 ──────────────────────────────

def test_orient_toggle_changes_fused_layout():
    """chrysene(rdkit 默认非水平) 开关应改变布局, 证明 orient 参数生效。"""
    smi = "c1ccc2c(c1)ccc3c4ccccc4ccc23"
    on = build_locants_svg(smi, orient=True)
    off = build_locants_svg(smi, orient=False)
    assert on is not None and off is not None
    assert on["svg"] != off["svg"]


def test_orient_defaults_on():
    """默认 orient=True: 稠环按优选取向摆放。"""
    res = build_locants_svg("c1ccc2c(c1)ccc3c4ccccc4ccc23")
    assert res is not None
    assert "svg" in res


def test_orient_ignored_for_single_ring():
    """单环无优选取向, 开关不影响输出。"""
    smi = "C1CCCCC1"
    on = build_locants_svg(smi, orient=True)
    off = build_locants_svg(smi, orient=False)
    assert on is not None and off is not None
    assert on["svg"] == off["svg"]


def test_api_orient_param():
    c = TestClient(app)
    r = c.post("/api/v1/name/locants-svg",
               json={"smiles": "c1ccc2c(c1)ccc3c4ccccc4ccc23", "orient": False})
    assert r.status_code == 200
    assert r.json()["ok"] is True
