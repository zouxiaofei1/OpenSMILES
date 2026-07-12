# tests/unit/test_secrets.py
from pathlib import Path
from agent_loop.secrets import SecretsStore

def test_set_get_configured(tmp_path: Path):
    store = SecretsStore(tmp_path / "secrets.json")
    assert store.is_configured() is False
    assert store.public_status() == {"configured": False, "hint": None}
    store.set_key("sk-ant-test-1234")
    assert store.is_configured() is True
    assert store.get_key() == "sk-ant-test-1234"
    st = store.public_status()
    assert st["configured"] is True
    assert st["hint"] == "1234"
    assert "sk-ant" not in (st["hint"] or "")

def test_clear_and_empty_set(tmp_path: Path):
    store = SecretsStore(tmp_path / "secrets.json")
    store.set_key("abcd")
    store.clear_key()
    assert store.get_key() is None
    store.set_key("  ")
    assert store.is_configured() is False

def test_corrupt_file_unconfigured(tmp_path: Path):
    path = tmp_path / "secrets.json"
    path.write_text("{not json", encoding="utf-8")
    store = SecretsStore(path)
    assert store.is_configured() is False
