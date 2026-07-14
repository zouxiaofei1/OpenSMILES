"""One-shot patch: multi-ring + aryl for diazine/hetero5 in heteroarene5.py."""
from __future__ import annotations

from pathlib import Path

p = Path("src/namepredict/layer2/heteroarene5.py")
t = p.read_text(encoding="utf-8")

if "from namepredict.layer2.aryl_sub import" not in t:
    t = t.replace(
        "from namepredict.layer2.side_alkyl import _linear_n_alkyl_sides_ok\n",
        "from namepredict.layer2.side_alkyl import _linear_n_alkyl_sides_ok\n"
        "from namepredict.layer2.aryl_sub import _aryl_atoms, _aryl_exclude, _aryl_sub_n\n",
    )
    print("import added")

needle = (
    "def _ring_atoms_if_mono(info: dict) -> list[int] | None:\n"
    "    rings = info.get(\"rings\") or []\n"
    "    if len(rings) != 1:\n"
    "        return None\n"
    "    return list(rings[0][\"atom_ids\"])\n"
    "\n"
    "\n"
    "def _hetero5_zs"
)
insert = (
    "def _ring_atoms_if_mono(info: dict) -> list[int] | None:\n"
    "    rings = info.get(\"rings\") or []\n"
    "    if len(rings) != 1:\n"
    "        return None\n"
    "    return list(rings[0][\"atom_ids\"])\n"
    "\n"
    "\n"
    "def _unfused_ring(mol: Mol, ring: set[int]) -> bool:\n"
    "    for i in ring:\n"
    "        if sum(1 for r in mol.GetRingInfo().AtomRings() if i in r) != 1:\n"
    "            return False\n"
    "    return True\n"
    "\n"
    "\n"
    "def _arom_ring_lists(mol: Mol, n: int, pred) -> list[list[int]]:\n"
    "    out: list[list[int]] = []\n"
    "    for r in mol.GetRingInfo().AtomRings():\n"
    "        ids = list(r)\n"
    "        if len(ids) != n or not _unfused_ring(mol, set(ids)):\n"
    "            continue\n"
    "        if not all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ids):\n"
    "            continue\n"
    "        if pred(mol, ids):\n"
    "            out.append(ids)\n"
    "    return out\n"
    "\n"
    "\n"
    "def _is_c4n2(mol: Mol, atoms: list[int]) -> bool:\n"
    "    zs = _hetero5_zs(mol, atoms)\n"
    "    return zs.count(7) == 2 and zs.count(6) == 4\n"
    "\n"
    "\n"
    "def _is_c4x_mono(mol: Mol, atoms: list[int]) -> bool:\n"
    "    zs = _hetero5_zs(mol, atoms)\n"
    "    return zs.count(6) == 4 and sum(1 for z in zs if z in _HETERO5_KIND) == 1\n"
    "\n"
    "\n"
    "def _diazine_rings(mol: Mol) -> list[list[int]]:\n"
    "    return _arom_ring_lists(mol, 6, _is_c4n2)\n"
    "\n"
    "\n"
    "def _hetero5_rings(mol: Mol) -> list[list[int]]:\n"
    "    return _arom_ring_lists(mol, 5, _is_c4x_mono)\n"
    "\n"
    "\n"
    "def _hetero5_zs"
)
if needle not in t:
    if "def _diazine_rings" in t:
        print("helpers already present")
    else:
        raise SystemExit("needle1 missing")
else:
    t = t.replace(needle, insert)
    print("helpers inserted")

old = (
    "def _is_hetero5_core(info: dict) -> bool:\n"
    "    atom_ids = _ring_atoms_if_mono(info)\n"
    "    if atom_ids is None or len(atom_ids) != 5:\n"
    "        return False\n"
    "    mol: Mol = info[\"mol\"]\n"
    "    if not all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in atom_ids):\n"
    "        return False\n"
    "    zs = _hetero5_zs(mol, atom_ids)\n"
    "    return zs.count(6) == 4 and sum(1 for z in zs if z in _HETERO5_KIND) == 1\n"
    "\n"
    "\n"
    "def _hetero5_idx(info: dict) -> int | None:\n"
    "    if not _is_hetero5_core(info):\n"
    "        return None\n"
    "    mol: Mol = info[\"mol\"]\n"
    "    for i in info[\"rings\"][0][\"atom_ids\"]:\n"
    "        if mol.GetAtomWithIdx(i).GetAtomicNum() in _HETERO5_KIND:\n"
    "            return i\n"
    "    return None\n"
    "\n"
    "\n"
    "def _hetero5_kind(info: dict) -> str | None:\n"
    "    idx = _hetero5_idx(info)\n"
    "    if idx is None:\n"
    "        return None\n"
    "    z = info[\"mol\"].GetAtomWithIdx(idx).GetAtomicNum()\n"
    "    if z == 7 and info[\"mol\"].GetAtomWithIdx(idx).GetTotalNumHs() < 1:\n"
    "        return None\n"
    "    return _HETERO5_KIND.get(z)\n"
)
new = (
    "def _is_hetero5_core(info: dict) -> bool:\n"
    "    return bool(_hetero5_rings(info[\"mol\"]))\n"
    "\n"
    "\n"
    "def _hetero5_idx_in(mol: Mol, ring: list[int]) -> int | None:\n"
    "    for i in ring:\n"
    "        if mol.GetAtomWithIdx(i).GetAtomicNum() in _HETERO5_KIND:\n"
    "            return i\n"
    "    return None\n"
    "\n"
    "\n"
    "def _hetero5_kind_of(mol: Mol, ring: list[int]) -> str | None:\n"
    "    idx = _hetero5_idx_in(mol, ring)\n"
    "    if idx is None:\n"
    "        return None\n"
    "    z = mol.GetAtomWithIdx(idx).GetAtomicNum()\n"
    "    if z == 7 and mol.GetAtomWithIdx(idx).GetTotalNumHs() < 1:\n"
    "        return None\n"
    "    return _HETERO5_KIND.get(z)\n"
    "\n"
    "\n"
    "def _hetero5_idx(info: dict) -> int | None:\n"
    "    rings = _hetero5_rings(info[\"mol\"])\n"
    "    return None if not rings else _hetero5_idx_in(info[\"mol\"], rings[0])\n"
    "\n"
    "\n"
    "def _hetero5_kind(info: dict) -> str | None:\n"
    "    rings = _hetero5_rings(info[\"mol\"])\n"
    "    return None if not rings else _hetero5_kind_of(info[\"mol\"], rings[0])\n"
)
if old not in t:
    if "def _hetero5_kind_of" in t:
        print("h5 core already patched")
    else:
        raise SystemExit("h5 core missing")
else:
    t = t.replace(old, new)
    print("h5 core ok")

old = (
    "def _hetero5_subs_ok(mol: Mol, ring: set[int]) -> bool:\n"
    "    \"\"\"Allow ≤2 ring simple subs: monohalo and/or n-alkyl C1–C2 (P-14.3.4).\"\"\"\n"
    "    h, starts = _ring_halo_n(mol, ring), _ring_side_starts(mol, ring)\n"
    "    if h + len(starts) > 2:\n"
    "        return False\n"
    "    return _linear_n_alkyl_sides_ok(mol, ring, starts, 2)\n"
    "\n"
    "\n"
    "def _is_simple_hetero5(info: dict) -> bool:\n"
    "    if _hetero5_kind(info) is None or _hetero5_fg_block(info):\n"
    "        return False\n"
    "    mol, ring = info[\"mol\"], set(info[\"rings\"][0][\"atom_ids\"])\n"
    "    if not _outside_ok(mol, ring):\n"
    "        return False\n"
    "    return _hetero5_subs_ok(mol, ring)\n"
    "\n"
    "\n"
    "def _hetero5_parent(info: dict) -> dict:\n"
    "    ring = list(info[\"rings\"][0][\"atom_ids\"])\n"
    "    kind = _hetero5_kind(info)\n"
    "    return {\n"
    "        \"chain\": ring, \"n_carbons\": 5, \"kind\": kind,\n"
    "        \"hetero_idx\": _hetero5_idx(info),\n"
    "    }\n"
)
new = (
    "def _hetero5_subs_ok(info: dict, mol: Mol, ring: set[int]) -> bool:\n"
    "    \"\"\"Allow ≤2 ring subs: halo / methyl / aryl (P-14.3.4 / P-29.3).\"\"\"\n"
    "    n_aryl = _aryl_sub_n(info, ring)\n"
    "    excl = _aryl_exclude(info, ring)\n"
    "    h = _ring_halo_n(mol, ring)\n"
    "    starts = _ring_side_starts(mol, ring, excl)\n"
    "    if h + len(starts) + n_aryl > 2:\n"
    "        return False\n"
    "    return True if not starts else _linear_n_alkyl_sides_ok(mol, ring, starts, 2)\n"
    "\n"
    "\n"
    "def _h5_ring_ok(info: dict, ring: list[int]) -> bool:\n"
    "    mol, rs = info[\"mol\"], set(ring)\n"
    "    if not _outside_ok(mol, rs, _aryl_atoms(info, rs)):\n"
    "        return False\n"
    "    return _hetero5_subs_ok(info, mol, rs)\n"
    "\n"
    "\n"
    "def _pick_hetero5_ring(info: dict) -> list[int] | None:\n"
    "    cands = [r for r in _hetero5_rings(info[\"mol\"]) if _h5_ring_ok(info, r)]\n"
    "    return min(cands, key=min) if cands else None\n"
    "\n"
    "\n"
    "def _is_simple_hetero5(info: dict) -> bool:\n"
    "    if _hetero5_fg_block(info):\n"
    "        return False\n"
    "    return _pick_hetero5_ring(info) is not None\n"
    "\n"
    "\n"
    "def _hetero5_parent(info: dict) -> dict:\n"
    "    ring = _pick_hetero5_ring(info) or _hetero5_rings(info[\"mol\"])[0]\n"
    "    mol = info[\"mol\"]\n"
    "    return {\n"
    "        \"chain\": ring, \"n_carbons\": 5, \"kind\": _hetero5_kind_of(mol, ring),\n"
    "        \"hetero_idx\": _hetero5_idx_in(mol, ring),\n"
    "    }\n"
)
if old not in t:
    if "def _pick_hetero5_ring" in t:
        print("h5 simple already")
    else:
        raise SystemExit("h5 simple missing")
else:
    t = t.replace(old, new)
    print("h5 simple ok")

old = (
    "def _is_diazine_core(info: dict) -> bool:\n"
    "    atom_ids = _ring_atoms_if_mono(info)\n"
    "    if atom_ids is None or len(atom_ids) != 6:\n"
    "        return False\n"
    "    mol: Mol = info[\"mol\"]\n"
    "    if not all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in atom_ids):\n"
    "        return False\n"
    "    zs = _hetero5_zs(mol, atom_ids)\n"
    "    return zs.count(7) == 2 and zs.count(6) == 4\n"
    "\n"
    "\n"
    "def _diazine_n_idxs(info: dict) -> list[int] | None:\n"
    "    if not _is_diazine_core(info):\n"
    "        return None\n"
    "    mol: Mol = info[\"mol\"]\n"
    "    return [i for i in info[\"rings\"][0][\"atom_ids\"] if mol.GetAtomWithIdx(i).GetAtomicNum() == 7]\n"
    "\n"
    "\n"
    "def _ring_nn_dist(atom_ids: list[int], n_idxs: list[int]) -> int:\n"
    "    ia, ib = atom_ids.index(n_idxs[0]), atom_ids.index(n_idxs[1])\n"
    "    d = abs(ia - ib)\n"
    "    return min(d, len(atom_ids) - d)\n"
    "\n"
    "\n"
    "def _diazine_kind(info: dict) -> str | None:\n"
    "    n_idxs = _diazine_n_idxs(info)\n"
    "    if n_idxs is None or len(n_idxs) != 2:\n"
    "        return None\n"
    "    ring = list(info[\"rings\"][0][\"atom_ids\"])\n"
    "    return _DIAZINE_KIND.get(_ring_nn_dist(ring, n_idxs))\n"
)
new = (
    "def _is_diazine_core(info: dict) -> bool:\n"
    "    return bool(_diazine_rings(info[\"mol\"]))\n"
    "\n"
    "\n"
    "def _n_idxs_in(mol: Mol, ring: list[int]) -> list[int]:\n"
    "    return [i for i in ring if mol.GetAtomWithIdx(i).GetAtomicNum() == 7]\n"
    "\n"
    "\n"
    "def _diazine_n_idxs(info: dict, ring: list[int] | None = None) -> list[int] | None:\n"
    "    rings = _diazine_rings(info[\"mol\"])\n"
    "    if ring is None:\n"
    "        ring = rings[0] if rings else None\n"
    "    if ring is None:\n"
    "        return None\n"
    "    ns = _n_idxs_in(info[\"mol\"], ring)\n"
    "    return ns if len(ns) == 2 else None\n"
    "\n"
    "\n"
    "def _ring_nn_dist(atom_ids: list[int], n_idxs: list[int]) -> int:\n"
    "    ia, ib = atom_ids.index(n_idxs[0]), atom_ids.index(n_idxs[1])\n"
    "    d = abs(ia - ib)\n"
    "    return min(d, len(atom_ids) - d)\n"
    "\n"
    "\n"
    "def _diazine_kind_of(ring: list[int], n_idxs: list[int]) -> str | None:\n"
    "    return _DIAZINE_KIND.get(_ring_nn_dist(ring, n_idxs))\n"
    "\n"
    "\n"
    "def _diazine_kind(info: dict, ring: list[int] | None = None) -> str | None:\n"
    "    rings = _diazine_rings(info[\"mol\"])\n"
    "    ring = ring or (rings[0] if rings else None)\n"
    "    ns = _diazine_n_idxs(info, ring)\n"
    "    if ring is None or ns is None:\n"
    "        return None\n"
    "    return _diazine_kind_of(ring, ns)\n"
)
if old not in t:
    if "def _diazine_kind_of" in t:
        print("dz core already")
    else:
        raise SystemExit("dz core missing")
else:
    t = t.replace(old, new)
    print("dz core ok")

old = (
    "def _diazine_subs_ok(info: dict, mol: Mol, ring: set[int]) -> bool:\n"
    "    \"\"\"Allow ≤4 simple ring subs: halo / methyl / methoxy / nitro (P-14.3.4).\"\"\"\n"
    "    alk, n_alk = _arene_alkoxy(info, ring)\n"
    "    h = _ring_halo_n(mol, ring)\n"
    "    starts = _ring_side_starts(mol, ring, alk)\n"
    "    n_sub = h + len(starts) + _ring_nitro_n(info, ring) + n_alk\n"
    "    if n_sub > 4:\n"
    "        return False\n"
    "    return _diazine_side_methyl_ok(mol, ring, starts, alk)\n"
    "\n"
    "\n"
    "def _diazine_outside_ok(info: dict, mol: Mol, ring: set[int]) -> bool:\n"
    "    alk = _arene_alkoxy(info, ring)[0]\n"
    "    allowed = _ring_nitro_atoms(info, ring) | alk\n"
    "    return _outside_ok(mol, ring, allowed)\n"
    "\n"
    "\n"
    "def _is_simple_diazine(info: dict) -> bool:\n"
    "    if _diazine_kind(info) is None or _diazine_fg_block(info):\n"
    "        return False\n"
    "    mol, ring = info[\"mol\"], set(info[\"rings\"][0][\"atom_ids\"])\n"
    "    if not _diazine_outside_ok(info, mol, ring):\n"
    "        return False\n"
    "    return _diazine_subs_ok(info, mol, ring)\n"
    "\n"
    "\n"
    "def _diazine_parent(info: dict) -> dict:\n"
    "    ring = list(info[\"rings\"][0][\"atom_ids\"])\n"
    "    return {\n"
    "        \"chain\": ring, \"n_carbons\": 6, \"kind\": _diazine_kind(info),\n"
    "        \"n_idxs\": _diazine_n_idxs(info),\n"
    "    }\n"
)
new = (
    "def _diazine_subs_ok(info: dict, mol: Mol, ring: set[int]) -> bool:\n"
    "    \"\"\"Allow ≤4 ring subs: halo / methyl / alkoxy / nitro / aryl.\"\"\"\n"
    "    alk, n_alk = _arene_alkoxy(info, ring)\n"
    "    excl = alk | _aryl_exclude(info, ring)\n"
    "    h = _ring_halo_n(mol, ring)\n"
    "    starts = _ring_side_starts(mol, ring, excl)\n"
    "    n_sub = h + len(starts) + _ring_nitro_n(info, ring) + n_alk + _aryl_sub_n(info, ring)\n"
    "    if n_sub > 4:\n"
    "        return False\n"
    "    return True if not starts else _diazine_side_methyl_ok(mol, ring, starts, excl)\n"
    "\n"
    "\n"
    "def _diazine_outside_ok(info: dict, mol: Mol, ring: set[int]) -> bool:\n"
    "    alk = _arene_alkoxy(info, ring)[0]\n"
    "    allowed = _ring_nitro_atoms(info, ring) | alk | _aryl_atoms(info, ring)\n"
    "    return _outside_ok(mol, ring, allowed)\n"
    "\n"
    "\n"
    "def _dz_ring_ok(info: dict, ring: list[int]) -> bool:\n"
    "    mol, rs = info[\"mol\"], set(ring)\n"
    "    return _diazine_outside_ok(info, mol, rs) and _diazine_subs_ok(info, mol, rs)\n"
    "\n"
    "\n"
    "def _pick_diazine_ring(info: dict) -> list[int] | None:\n"
    "    cands = [r for r in _diazine_rings(info[\"mol\"]) if _dz_ring_ok(info, r)]\n"
    "    return min(cands, key=min) if cands else None\n"
    "\n"
    "\n"
    "def _is_simple_diazine(info: dict) -> bool:\n"
    "    if _diazine_fg_block(info):\n"
    "        return False\n"
    "    return _pick_diazine_ring(info) is not None\n"
    "\n"
    "\n"
    "def _diazine_parent(info: dict) -> dict:\n"
    "    ring = _pick_diazine_ring(info) or _diazine_rings(info[\"mol\"])[0]\n"
    "    ns = _diazine_n_idxs(info, ring) or []\n"
    "    return {\n"
    "        \"chain\": ring, \"n_carbons\": 6, \"kind\": _diazine_kind_of(ring, ns),\n"
    "        \"n_idxs\": ns,\n"
    "    }\n"
)
if old not in t:
    if "def _pick_diazine_ring" in t:
        print("dz simple already")
    else:
        raise SystemExit("dz simple missing")
else:
    t = t.replace(old, new)
    print("dz simple ok")

# pyrimidinamine multi-ring
old = (
    "def _is_simple_pyrimidinamine(info: dict) -> bool:\n"
    "    if _diazine_kind(info) != \"pyrimidine\" or _pyrimidinamine_fg_block(info):\n"
    "        return False\n"
    "    mol, ring = info[\"mol\"], set(info[\"rings\"][0][\"atom_ids\"])\n"
    "    am = _mono_amine_on_ring(info, ring)\n"
    "    if am is None or am.get(\"degree\") != 1:\n"
    "        return False\n"
    "    if not _pyrimidinamine_outside(info, mol, ring, am[\"n_idx\"]):\n"
    "        return False\n"
    "    return _pyrimidinamine_subs_ok(info, mol, ring, am[\"n_idx\"])\n"
    "\n"
    "\n"
    "def _pyrimidinamine_parent(info: dict) -> dict:\n"
    "    ring = list(info[\"rings\"][0][\"atom_ids\"])\n"
    "    return {\n"
    "        \"chain\": ring, \"n_carbons\": 6, \"kind\": \"pyrimidinamine\",\n"
    "        \"n_idxs\": _diazine_n_idxs(info),\n"
    "        \"amine_c_idx\": info[\"amines\"][0][\"c_idx\"],\n"
    "    }\n"
)
new = (
    "def _is_simple_pyrimidinamine(info: dict) -> bool:\n"
    "    ring = _pick_diazine_ring(info)\n"
    "    if ring is None or _diazine_kind(info, ring) != \"pyrimidine\":\n"
    "        return False\n"
    "    if _pyrimidinamine_fg_block(info):\n"
    "        return False\n"
    "    mol, rs = info[\"mol\"], set(ring)\n"
    "    am = _mono_amine_on_ring(info, rs)\n"
    "    if am is None or am.get(\"degree\") != 1:\n"
    "        return False\n"
    "    if not _pyrimidinamine_outside(info, mol, rs, am[\"n_idx\"]):\n"
    "        return False\n"
    "    return _pyrimidinamine_subs_ok(info, mol, rs, am[\"n_idx\"])\n"
    "\n"
    "\n"
    "def _pyrimidinamine_parent(info: dict) -> dict:\n"
    "    ring = _pick_diazine_ring(info) or _diazine_rings(info[\"mol\"])[0]\n"
    "    return {\n"
    "        \"chain\": ring, \"n_carbons\": 6, \"kind\": \"pyrimidinamine\",\n"
    "        \"n_idxs\": _diazine_n_idxs(info, ring),\n"
    "        \"amine_c_idx\": info[\"amines\"][0][\"c_idx\"],\n"
    "    }\n"
)
if old in t:
    t = t.replace(old, new)
    print("pyrimidinamine ok")
else:
    print("pyrimidinamine skip")

p.write_text(t, encoding="utf-8")
import ast

ast.parse(t)
print("lines", len(t.splitlines()), "syntax ok")
