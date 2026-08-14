# IUPAC: P-29.2
# Layer: L2,L3
"""Architecture contract for Layer 2 side topology facts consumed by Layer 3."""
from __future__ import annotations

import ast
from pathlib import Path


LAYER3 = Path(__file__).parents[2] / "src" / "namepredict" / "layer3"
TOOLS = Path(__file__).parents[2] / "src" / "namepredict" / "tools"
# Side-topology facts merged into tools 碳拓扑原语 chain.py。
TOOLS_FACTS = TOOLS / "chain.py"
LEAF_PROTOCOL = LAYER3 / "leaves" / "protocol.py"


def _private_layer2_imports(path: Path) -> list[str]:
    tree = ast.parse(path.read_text(encoding="utf-8"))
    found: list[str] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom) and (node.module or "").startswith("namepredict.layer2"):
            found.extend(alias.name for alias in node.names if alias.name.startswith("_"))
    return found


def _all_functions(path: Path) -> list[ast.FunctionDef | ast.AsyncFunctionDef]:
    tree = ast.parse(path.read_text(encoding="utf-8"))
    return [n for n in ast.walk(tree) if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef))]


def _tools_modules(tree: ast.AST) -> list[str]:
    direct = [alias.name for node in ast.walk(tree) if isinstance(node, ast.Import)
              for alias in node.names if alias.name.startswith("namepredict.tools")]
    froms = [(node.module or "") for node in ast.walk(tree) if isinstance(node, ast.ImportFrom)
             and (node.module or "").startswith("namepredict.tools")]
    return [*direct, *froms]


def test_tools_import_scan_covers_both_ast_forms() -> None:
    tree = ast.parse("import namepredict.tools.aryl_sub\nfrom namepredict.tools import side_alkyl")
    assert _tools_modules(tree) == [
        "namepredict.tools.aryl_sub", "namepredict.tools",
    ]


def test_side_facts_has_no_string_to_enum_dispatcher() -> None:
    tree = ast.parse(TOOLS_FACTS.read_text(encoding="utf-8"))
    mappings = [node.value for node in tree.body if isinstance(node, ast.Assign)
                and isinstance(node.value, ast.Dict)]
    bad = [mapping for mapping in mappings
           if any(isinstance(key, ast.Constant) and isinstance(key.value, str)
                  and isinstance(value, ast.Attribute) and value.attr.isupper()
                  for key, value in zip(mapping.keys, mapping.values))]
    assert not bad


def test_side_facts_does_not_dispatch_match_kind_keys() -> None:
    tree = ast.parse(TOOLS_FACTS.read_text(encoding="utf-8"))
    keys = [node.slice.value for node in ast.walk(tree) if isinstance(node, ast.Subscript)
            and isinstance(node.slice, ast.Constant) and isinstance(node.slice.value, str)]
    assert "kind" not in keys


def test_side_facts_has_no_naming_token_maps() -> None:
    tree = ast.parse(TOOLS_FACTS.read_text(encoding="utf-8"))
    assigned = {target.id for node in tree.body if isinstance(node, ast.Assign)
                for target in node.targets if isinstance(target, ast.Name)}
    assert not {"_LEAF_KINDS", "_ALKYL_SHAPES"} & assigned


def _enum_members(node: ast.ClassDef) -> list[ast.Assign]:
    return [item for item in node.body if isinstance(item, ast.Assign)]


def _is_auto_member(node: ast.Assign) -> bool:
    value = node.value
    return (isinstance(value, ast.Call) and isinstance(value.func, ast.Name)
            and value.func.id == "auto" and not value.args and not value.keywords)


def test_side_facts_has_no_naming_dependencies() -> None:
    tree = ast.parse(TOOLS_FACTS.read_text(encoding="utf-8"))
    imports = [(node.module or "", alias.name) for node in ast.walk(tree)
               if isinstance(node, ast.ImportFrom) for alias in node.names]
    calls = [node.func.id for node in ast.walk(tree) if isinstance(node, ast.Call)
             and isinstance(node.func, ast.Name)]
    assert not [(module, name) for module, name in imports if name == "name_leaf"]
    assert "name_leaf" not in calls


def test_public_side_fact_functions_are_typed_and_explicit() -> None:
    public = [fn for fn in _all_functions(TOOLS_FACTS) if not fn.name.startswith("_")]
    assert all(fn.returns is not None for fn in public)
    assert all(fn.args.vararg is None and fn.args.kwarg is None for fn in public)
    assert all(all(arg.annotation is not None for arg in fn.args.args) for fn in public)


def test_layer3_functions_are_explicit() -> None:
    functions = [fn for path in LAYER3.rglob("*.py") for fn in _all_functions(path)]
    assert all(fn.args.vararg is None and fn.args.kwarg is None for fn in functions)


def test_side_fact_contract_has_no_public_dict_returns() -> None:
    public = [fn for fn in _all_functions(TOOLS_FACTS) if not fn.name.startswith("_")]
    returns = [ast.unparse(fn.returns) for fn in public]
    assert not any("dict" in annotation for annotation in returns), returns


def test_side_fact_contract_has_no_public_name_renderers() -> None:
    tree = ast.parse(TOOLS_FACTS.read_text(encoding="utf-8"))
    public = [n for n in tree.body
              if isinstance(n, ast.FunctionDef) and not n.name.startswith("_")]
    banned = [fn.name for fn in public if "render" in fn.name or "name" in fn.name]
    string_tuples = [fn.name for fn in public if "tuple[str" in ast.unparse(fn.returns)]
    assert not banned
    assert not string_tuples


def test_side_fact_dataclasses_expose_only_topology_fields() -> None:
    tree = ast.parse(TOOLS_FACTS.read_text(encoding="utf-8"))
    classes = [n for n in tree.body if isinstance(n, ast.ClassDef)]
    forbidden = {"name", "en", "zh", "paren"}
    for cls in classes:
        fields = {n.target.id for n in cls.body if isinstance(n, ast.AnnAssign)
                  and isinstance(n.target, ast.Name)}
        assert not fields & forbidden, (cls.name, fields & forbidden)
        methods = [n.name for n in cls.body if isinstance(n, ast.FunctionDef)]
        assert methods == [], (cls.name, methods)
