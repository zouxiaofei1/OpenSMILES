from __future__ import annotations

from rdkit.Chem import Mol

import namepredict.layer2.side_facts as side_facts
from namepredict.layer3.aryl_names import aryl_arm_name, heteroaryl_name
from namepredict.layer3.alkoxy_names import _extract_alkoxys
from namepredict.layer3.cycloalkyl_names import _one_cycloalkyl_side
from namepredict.layer3.amino_side import _extract_aminos as _extract_aminos_impl


def _strip_ital_prefix(stem: str) -> str:
    if stem.startswith("tert-") or stem.startswith("sec-"):
        return stem[stem.index("-") + 1 :]
    return stem


def _strip_n_prefix(stem: str) -> str:
    if stem.startswith("N,"):
        return stem.split("-")[-1] if "-" in stem else stem
    return stem[2:] if stem.startswith("N-") else stem


def _strip_lead_locant(stem: str) -> str:
    """Strip one leading locant set: '4-', '1,3-', '1,1,1-', '1H-' (P-14.5)."""
    i = 0
    n = len(stem)
    while i < n and stem[i].isdigit():
        i += 1
        while i < n and stem[i] == ",":
            i += 1
            while i < n and stem[i].isdigit():
                i += 1
    # Indicated-hydrogen prefix: 1H-, 2H-, 3H- (P-14.5 / P-65.3.2.5)
    if i and i + 1 < n and stem[i] == "H" and stem[i + 1] == "-":
        i += 1
    return stem[i + 1 :] if i and i < n and stem[i] == "-" else stem


def _strip_outer_parens(stem: str) -> str:
    if len(stem) >= 2 and stem[0] == "(" and stem[-1] == ")":
        return stem[1:-1]
    return stem


def alkyl_alpha_key(stem: str) -> str:
    """Alphanumerical-order key: sec-/tert-/N-/parens/leading locants ignored (P-14.5)."""
    s = _strip_n_prefix(_strip_ital_prefix(stem))
    s = _strip_outer_parens(s)
    return _strip_lead_locant(s)


# n-alkyl C1–C12 (P-29.3 / BQ P-64.2 Round B); C11+ ZH uses …烷基
ALKYL_EN = {
    1: "methyl", 2: "ethyl", 3: "propyl", 4: "butyl", 5: "pentyl", 6: "hexyl",
    7: "heptyl", 8: "octyl", 9: "nonyl", 10: "decyl", 11: "undecyl", 12: "dodecyl",
}
ALKYL_ZH = {
    1: "甲基", 2: "乙基", 3: "丙基", 4: "丁基", 5: "戊基", 6: "己基",
    7: "庚基", 8: "辛基", 9: "壬基", 10: "癸基", 11: "十一烷基", 12: "十二烷基",
}
HALO_EN = {9: "fluoro", 17: "chloro", 35: "bromo", 53: "iodo"}
HALO_ZH = {9: "氟", 17: "氯", 35: "溴", 53: "碘"}
_HALOALKYL_STEM_EN = {1: "methyl", 2: "ethyl", 3: "propyl", 4: "butyl"}
_HALOALKYL_STEM_ZH = {1: "甲基", 2: "乙基", 3: "丙基", 4: "丁基"}


def _side_starts(mol: Mol, chain: list[int]) -> list[tuple[int, int]]:
    cs = set(chain)
    return [(c, n) for c in chain for n in side_facts.carbon_neighbors(mol, c) if n not in cs]


def _make_alkyl(attach: int, path: list[int]) -> dict:
    n = len(path)
    return {
        "kind": "alkyl", "n_carbons": n, "attach_idx": attach,
        "atoms": path, "en": ALKYL_EN[n], "zh": ALKYL_ZH[n],
    }


def _make_branch(attach: int, atoms: list[int], n: int, en: str, zh: str) -> dict:
    return {
        "kind": "alkyl", "n_carbons": n, "attach_idx": attach,
        "atoms": atoms, "en": en, "zh": zh,
    }


def _make_cf3(attach: int, atoms: list[int]) -> dict:
    return {
        "kind": "trifluoromethyl", "n_carbons": 1, "attach_idx": attach,
        "atoms": atoms, "en": "trifluoromethyl", "zh": "三氟甲基",
    }


def _haloalkyl_names(n: int, z: int) -> tuple[str, str]:
    se, sz = _HALOALKYL_STEM_EN[n], _HALOALKYL_STEM_ZH[n]
    if n == 1:
        return f"{HALO_EN[z]}{se}", f"{HALO_ZH[z]}{sz}"
    return f"{n}-{HALO_EN[z]}{se}", f"{n}-{HALO_ZH[z]}{sz}"


def _halo_atom(mol: Mol, carbon: int) -> int | None:
    atom = mol.GetAtomWithIdx(carbon)
    hits = [nb.GetIdx() for nb in atom.GetNeighbors() if nb.GetAtomicNum() in HALO_EN]
    return hits[0] if len(hits) == 1 else None


def _make_haloalkyl(attach: int, path: list[int], z: int, halo: int) -> dict:
    n = len(path)
    en, zh = _haloalkyl_names(n, z)
    return {
        "kind": "haloalkyl", "n_carbons": n, "attach_idx": attach,
        "atoms": [*path, halo], "en": en, "zh": zh, "paren": True,
    }


def _one_haloalkyl(mol: Mol, attach: int, start: int, chain_set: set[int]) -> dict | None:
    fact = side_facts.omega_halo_alkyl(mol, start, chain_set)
    path = list(fact.atoms) if fact else None
    if not path:
        return None
    z = side_facts.terminal_halogen(mol, path[-1])
    halo = _halo_atom(mol, path[-1])
    return _make_haloalkyl(attach, path, z, halo) if z and halo is not None else None


def _one_alkyl(mol: Mol, attach: int, start: int, chain_set: set[int], *, name_mode: str = "general") -> dict | None:
    fact = side_facts.linear_alkyl(mol, start, chain_set, 12)
    path = list(fact.atoms) if fact else None
    if path and len(path) in ALKYL_EN:
        return _make_alkyl(attach, path)
    ha = _one_haloalkyl(mol, attach, start, chain_set)
    return ha if ha is not None else _one_branched(mol, attach, start, chain_set, name_mode=name_mode)


_BRANCH_CHECKS = (
    (side_facts.AlkylShape.C2_VINYL, 2, "vinyl"),
    (side_facts.AlkylShape.C3_ALLYL, 3, "allyl"),
    (side_facts.AlkylShape.C3_ISOPROPENYL, 3, "isopropenyl"),
    (side_facts.AlkylShape.C3_BRANCH_AT_ROOT, 3, "isopropyl"),
    (side_facts.AlkylShape.C4_TRIPLE_BRANCH_AT_ROOT, 4, "tert-butyl"),
    (side_facts.AlkylShape.C5_ASYMMETRIC_ROOT_BRANCH, 5, "2-methylbutan-2-yl"),
    (side_facts.AlkylShape.C4_BRANCH_AFTER_ROOT, 4, "isobutyl"),
    (side_facts.AlkylShape.C4_BRANCH_AT_SECOND, 4, "sec-butyl"),
    (side_facts.AlkylShape.C5_DOUBLE_BRANCH_AFTER_ROOT, 5, "neopentyl"),
    (side_facts.AlkylShape.C5_PRENYL, 5, "3-methylbut-2-enyl"),
    (side_facts.AlkylShape.C5_BRANCH_NEAR_LEAF, 5, "isopentyl"),
)


def _retained_branch(mol, attach, start, chain_set, name_mode):
    from namepredict.layer3.retained_substituents import resolve_name

    for shape, n, key in _BRANCH_CHECKS:
        if fact := side_facts.alkyl_shape(mol, start, chain_set, shape):
            en, zh = resolve_name(key, name_mode=name_mode)
            return _make_branch(attach, list(fact.atoms), n, en, zh)
    return None


def _one_branched(mol: Mol, attach: int, start: int, chain_set: set[int], *, name_mode: str = "general") -> dict | None:
    if branch := _retained_branch(mol, attach, start, chain_set, name_mode):
        return branch
    if cyc := _one_cycloalkyl_side(mol, attach, start, chain_set):
        return cyc
    shape = side_facts.AlkylShape.C1_THREE_HALOGEN_LEAVES
    cf3 = side_facts.alkyl_shape(mol, start, chain_set, shape)
    return _make_cf3(attach, list(cf3.atoms)) if cf3 else None


def _make_halo(attach: int, halo_idx: int, z: int) -> dict:
    return {
        "kind": "halo", "attach_idx": attach, "atoms": [halo_idx],
        "en": HALO_EN[z], "zh": HALO_ZH[z],
    }


def _halo_on_carbon(mol: Mol, c_idx: int) -> list[dict]:
    return [
        _make_halo(c_idx, n.GetIdx(), n.GetAtomicNum())
        for n in mol.GetAtomWithIdx(c_idx).GetNeighbors()
        if n.GetAtomicNum() in HALO_EN
    ]


def _extract_halos(mol: Mol, chain: list[int]) -> list[dict]:
    return [h for c in chain for h in _halo_on_carbon(mol, c)]


def _extract_alkyls(mol: Mol, chain: list[int]) -> list[dict]:
    cs, out = set(chain), []
    for attach, start in _side_starts(mol, chain):
        one = _one_alkyl(mol, attach, start, cs)
        if one is not None:
            out.append(one)
    return out


def _filter_fg_halos(halos: list, parent: dict) -> list:
    # Functional-class ether arms already encode F (e.g. HFIP); do not re-prefix.
    if parent.get("kind") == "ether" and parent.get("ether_arms"):
        return []
    if parent.get("kind") not in ("acyl_chloride", "acyl_bromide"):
        return halos
    cl = parent.get("cl_idx") or parent.get("hal_idx")
    return [h for h in halos if cl not in (h.get("atoms") or [])]


_PARENT_OH_KINDS = frozenset({
    "alcohol", "diol", "triol", "cycloalcohol", "cycloalkanediol",
    "phenol", "benzenediol", "pyridinol", "benzothiophenol", "quinolinol",
    "naphthalenol", "naphthalenediol", "quinolinediol",
})
_PARENT_NH2_KINDS = frozenset({
    "amine", "diamine", "triamine", "tetraamine", "cycloamine", "sec_amine", "tert_amine",
    "aniline", "pyridinamine", "pyrimidinamine", "benzofuranamine",
    "benzothiazolamine", "benzoxazolamine", "benzimidazolamine", "benzenediamine",
    "naphthalenamine", "pyrazolamine", "thiazolamine", "quinazolinamine",
})
_PARENT_OXO_KINDS = frozenset(
    {
        "ketone", "dione", "cycloketone", "cycloalkanedione", "anthraquinone",
        "benzoquinone", "ortho_benzoquinone",
    },
)


def _make_hydroxy(attach: int, o_idx: int) -> dict:
    return {
        "kind": "hydroxy", "attach_idx": attach, "atoms": [o_idx],
        "en": "hydroxy", "zh": "羟基",
    }


def _make_oxo(attach: int) -> dict:
    return {
        "kind": "oxo", "attach_idx": attach, "atoms": [attach],
        "en": "oxo", "zh": "氧代",
    }


def _make_nitro(attach: int, n_idx: int, o_idxs: list[int]) -> dict:
    return {
        "kind": "nitro", "attach_idx": attach, "atoms": [n_idx] + list(o_idxs),
        "en": "nitro", "zh": "硝基",
    }


def _extract_nitros(info: dict, parent: dict) -> list[dict]:
    chain = set(parent.get("chain") or [])
    return [
        _make_nitro(n["c_idx"], n["n_idx"], n.get("o_idxs") or [])
        for n in info.get("nitros") or [] if n["c_idx"] in chain
    ]


def _make_iso(kind: str, attach: int, e: dict, en: str, zh: str) -> dict:
    return {
        "kind": kind, "attach_idx": attach,
        "atoms": [e["n_idx"], e["c_idx"], e["x_idx"]], "en": en, "zh": zh,
    }


def _extract_iso_kind(
    info: dict, parent: dict, key: str, kind: str, en: str, zh: str,
) -> list[dict]:
    if parent.get("kind") in ("isocyanate", "isothiocyanate"):
        return []
    chain = set(parent.get("chain") or [])
    return [
        _make_iso(kind, e["r_c_idx"], e, en, zh)
        for e in info.get(key) or [] if e["r_c_idx"] in chain
    ]


def _extract_isocyanates(info: dict, parent: dict) -> list[dict]:
    a = _extract_iso_kind(info, parent, "isocyanates", "isocyanato", "isocyanato", "异氰酸根合")
    b = _extract_iso_kind(
        info, parent, "isothiocyanates", "isothiocyanato", "isothiocyanato", "异硫氰酸根合",
    )
    return a + b


def _principal_attachments(parent: dict, group: str) -> frozenset[int]:
    facts = parent.get("principal_expression_facts")
    return facts.attachment_atoms if facts and facts.group_class.value == group else frozenset()


def _extract_hydroxys(info: dict, parent: dict) -> list[dict]:
    principal = _principal_attachments(parent, "alcohol")
    if parent.get("kind") in _PARENT_OH_KINDS and not principal:
        return []
    chain = set(parent.get("chain") or [])
    return [
        _make_hydroxy(h["c_idx"], h["o_idx"])
        for h in info.get("hydroxyls") or []
        if h["c_idx"] in chain and h["c_idx"] not in principal
    ]


def _extract_aminos(info: dict, parent: dict) -> list[dict]:
    return _extract_aminos_impl(info, parent, _PARENT_NH2_KINDS)


def _extract_oxos(info: dict, parent: dict) -> list[dict]:
    if parent.get("kind") in _PARENT_OXO_KINDS:
        return []
    chain = set(parent.get("chain") or [])
    return [
        _make_oxo(k["c_idx"]) for k in info.get("ketones") or [] if k["c_idx"] in chain
    ]


_AMIDE_KINDS = frozenset({"amide", "benzamide"})


def _n_benzyl_en_zh(en: str, zh: str, paren: bool) -> tuple[str, str]:
    if paren or en != "benzyl":
        return f"N-({en})", f"N-({zh})"
    return f"N-{en}", f"N-{zh}"


def _n_benzyl_sub(attach: int, en: str, zh: str, paren: bool) -> dict:
    ne, nz = _n_benzyl_en_zh(en, zh, paren)
    return {
        "kind": "n_benzyl", "n_carbons": 7, "attach_idx": attach,
        "atoms": [], "en": ne, "zh": nz,
    }


def _n_benzyl_named(info: dict, ch2: int) -> tuple[str, str, bool] | None:
    ams = info.get("amides") or []
    n_idx = ams[0]["n_idx"] if ams else -1
    ring = side_facts.benzyl_ring(info["mol"], ch2, n_idx)
    if ring is None:
        return None
    fact = side_facts.ArylArmFact(
        side_facts.ArylArmKind.METHYLENE_C, n_idx, ch2, ch2, ring, (ch2, *ring),
    )
    return aryl_arm_name(info["mol"], fact)


def _extract_n_benzyl(info: dict, parent: dict) -> list[dict]:
    if parent.get("kind") not in _AMIDE_KINDS or not parent.get("n_benzyl"):
        return []
    attach, ch2 = parent.get("amide_c_idx"), parent.get("n_benzyl_ch2")
    if attach is None or ch2 is None:
        return []
    named = _n_benzyl_named(info, ch2)
    return [_n_benzyl_sub(attach, *named)] if named else []


def _make_aryl(
    kind: str, attach: int, atoms: list[int], en: str, zh: str, paren: bool,
) -> dict:
    n = 7 if kind in ("benzyl", "benzyloxy") else 6
    return {
        "kind": kind, "attach_idx": attach, "atoms": atoms,
        "n_carbons": n, "en": en, "zh": zh, "paren": paren,
    }


_ARYL_KINDS = {
    side_facts.ArylArmKind.DIRECT_C: "phenyl",
    side_facts.ArylArmKind.METHYLENE_C: "benzyl",
    side_facts.ArylArmKind.DIRECT_O: "phenoxy",
    side_facts.ArylArmKind.O_METHYLENE_C: "benzyloxy",
}


def _one_aryl(mol: Mol, fact: side_facts.ArylArmFact) -> dict:
    en, zh, paren = aryl_arm_name(mol, fact)
    kind = _ARYL_KINDS[fact.kind]
    return _make_aryl(kind, fact.attachment, list(fact.atoms), en, zh, paren)


def _aryl_facts(info: dict, parent: dict) -> list[side_facts.ArylArmFact]:
    return side_facts.aryl_arms(info, set(parent.get("chain") or []))


def _one_heteroaryl(fact: side_facts.HeteroarylFact) -> dict:
    en, zh = heteroaryl_name(fact)
    pyridinyl = fact.kind == side_facts.HeteroarylKind.SIX_MEMBER_ONE_N
    return {
        "kind": "pyridinyl" if pyridinyl else "naphthyl",
        "attach_idx": fact.attachment, "atoms": list(fact.atoms),
        "n_carbons": 5 if pyridinyl else 10, "en": en, "zh": zh, "paren": True,
    }


def _aryl_outer_starts(info: dict, parent: dict) -> set[int]:
    mol, chain = info["mol"], set(parent.get("chain") or [])
    aromatic = {fact.outer for fact in side_facts.aryl_arms(info, chain)}
    return aromatic | side_facts.heteroaryl_outers(mol, chain)


def _extract_alkyls_no_aryl(mol: Mol, chain: list[int], skip: set[int], *, name_mode: str = "general") -> list[dict]:
    cs, out = set(chain), []
    for attach, start in _side_starts(mol, chain):
        if start in skip:
            continue
        one = _one_alkyl(mol, attach, start, cs, name_mode=name_mode)
        if one is not None:
            out.append(one)
    return out


def _extract_core_subs(info: dict, parent: dict) -> list:
    mol, chain = info["mol"], parent.get("chain") or []
    halo = _filter_fg_halos(_extract_halos(mol, chain), parent)
    return (
        halo + _extract_hydroxys(info, parent) + _extract_aminos(info, parent)
        + _extract_oxos(info, parent) + _extract_nitros(info, parent)
        + _extract_isocyanates(info, parent)
    )


def _extract_aryls(info: dict, parent: dict) -> list[dict]:
    mol, chain = info["mol"], set(parent.get("chain") or [])
    aryl = [_one_aryl(mol, fact) for fact in _aryl_facts(info, parent)]
    hetero = side_facts.pyridinyl_facts(mol, chain) + side_facts.naphthyl_facts(mol, chain)
    return aryl + [_one_heteroaryl(fact) for fact in hetero]


def _extract_n_subs(info: dict, parent: dict) -> list[dict]:
    from namepredict.layer3.n_block_extract import extract_n_blocks
    from namepredict.layer3.n_side_extract import extract_n_alkyl, extract_n_phenyl

    # Mutual exclusion: n_block claims complex N; skip simple n_alkyl/phenyl/benzyl.
    n_blocks = extract_n_blocks(info, parent)
    if n_blocks:
        return n_blocks
    mol = info.get("mol")
    return (
        extract_n_alkyl(parent, mol) + extract_n_phenyl(parent, mol)
        + _extract_n_benzyl(info, parent)
    )


def _extract_carboxymethyls(parent: dict) -> list[dict]:
    facts = parent.get("carboxymethyl_arms") or ()
    return [_carboxyalkyl_sub(fact) for fact in facts]


def _carboxyalkyl_sub(fact: side_facts.CarboxyalkylArm) -> dict:
    n = len(fact.atoms)
    names = {1: ("carboxymethyl", "羧甲基"), 2: ("2-carboxyethyl", "2-羧乙基"), 3: ("3-carboxypropyl", "3-羧丙基"), 4: ("4-carboxybutyl", "4-羧丁基")}
    en, zh = names.get(n, (f"{n}-carboxy{ALKYL_EN[n]}", f"{n}-羧{ALKYL_ZH[n]}"))
    return {"kind": "carboxyalkyl", "attach_idx": fact.attachment,
            "atoms": [*fact.atoms, fact.carboxyl], "en": en, "zh": zh, "paren": True}


def extract_substituents(info: dict, parent: dict, *, name_mode: str = "general") -> list:
    from namepredict.layer3.claim_extract import extract_claimed_sides

    mol, chain = info["mol"], parent.get("chain") or []
    alkyl = _extract_alkyls_no_aryl(mol, chain, _aryl_outer_starts(info, parent), name_mode=name_mode)
    base = (
        alkyl + _extract_carboxymethyls(parent) + _extract_core_subs(info, parent)
        + _extract_alkoxys(info, parent) + _extract_aryls(info, parent)
        + _extract_n_subs(info, parent)
    )
    return base + extract_claimed_sides(info, parent, base, name_mode=name_mode)
