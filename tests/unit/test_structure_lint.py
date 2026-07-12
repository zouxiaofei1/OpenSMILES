"""Unit tests for tools.structure_lint using temp violation trees."""

from __future__ import annotations

from pathlib import Path

from tools.structure_lint import lint_file, lint_tree, main


def _write(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path


def test_file_over_500_lines(tmp_path: Path):
    root = tmp_path / "pkg"
    f = _write(root / "big.py", "\n".join(f"x{i} = {i}" for i in range(501)) + "\n")
    issues = lint_file(f, root)
    assert any("501 lines" in m for m in issues)


def test_function_body_over_10(tmp_path: Path):
    root = tmp_path / "pkg"
    body = "\n".join(f"    a{i} = {i}" for i in range(11))
    src = f"def too_long():\n{body}\n"
    f = _write(root / "long_fn.py", src)
    issues = lint_file(f, root)
    assert any("too_long" in m and "11" in m for m in issues)


def test_function_body_exactly_10_ok(tmp_path: Path):
    root = tmp_path / "pkg"
    body = "\n".join(f"    a{i} = {i}" for i in range(10))
    src = f"def ok_fn():\n{body}\n"
    f = _write(root / "ok_fn.py", src)
    issues = lint_file(f, root)
    assert not any("ok_fn" in m for m in issues)


def test_bans_if_smiles_eq(tmp_path: Path):
    root = tmp_path / "pkg"
    src = 'def f(smiles):\n    if smiles == "CCO":\n        return "ethanol"\n'
    f = _write(root / "special.py", src)
    issues = lint_file(f, root)
    assert any("smiles equality" in m for m in issues)


def test_layer_cannot_import_higher(tmp_path: Path):
    root = tmp_path / "pkg"
    src = "from namepredict.layer2.parent_selector import select_parent\n"
    f = _write(root / "layer0" / "bad.py", src)
    issues = lint_file(f, root)
    assert any("layer0 imports higher layer2" in m for m in issues)


def test_layer_can_import_lower_or_same(tmp_path: Path):
    root = tmp_path / "pkg"
    src = (
        "from namepredict.layer0.preprocessor import preprocess\n"
        "from namepredict.layer2.x import y\n"
    )
    f = _write(root / "layer2" / "ok.py", src)
    issues = lint_file(f, root)
    assert not any("imports higher" in m for m in issues)


def test_common_names_preload_over_100(tmp_path: Path):
    root = tmp_path / "pkg"
    entries = ", ".join(f'"C{i}": "n{i}"' for i in range(101))
    src = f"PRELOAD = {{{entries}}}\n"
    f = _write(root / "cache" / "common_names.py", src)
    issues = lint_file(f, root)
    assert any("preloaded entries 101" in m for m in issues)


def test_common_names_preload_100_ok(tmp_path: Path):
    root = tmp_path / "pkg"
    entries = ", ".join(f'"C{i}": "n{i}"' for i in range(100))
    src = f"PRELOAD = {{{entries}}}\n"
    f = _write(root / "cache" / "common_names.py", src)
    issues = lint_file(f, root)
    assert not any("preloaded" in m for m in issues)


def test_lint_tree_aggregates(tmp_path: Path):
    root = tmp_path / "pkg"
    _write(root / "a.py", "x = 1\n")
    _write(
        root / "layer1" / "up.py",
        "from namepredict.layer5.assembler import assemble\n",
    )
    issues = lint_tree(root)
    assert any("layer1 imports higher layer5" in m for m in issues)


def test_cli_exit_codes(tmp_path: Path, capsys):
    good = tmp_path / "good"
    _write(good / "ok.py", "def f():\n    return 1\n")
    assert main(["--root", str(good)]) == 0

    bad = tmp_path / "bad"
    _write(bad / "x.py", 'def f(smiles):\n    if smiles == "x":\n        pass\n')
    assert main(["--root", str(bad)]) == 1
