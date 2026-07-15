from __future__ import annotations
from namepredict.layer3.substituent_extractor import alkyl_alpha_key
from namepredict.layer4.locants.adapt import (
    FUSED56_KINDS as _FUSED56_KINDS,
    INDOLE_LOCANTS as _INDOLE_LOCANTS,
    INDOLE_ORIENT_KINDS as _INDOLE_ORIENT_KINDS,
    NAPH_KINDS as _NAPH_KINDS,
    NAPH_LOCANTS as _NAPH_LOCANTS,
    Q_KINDS as _Q_KINDS,
    effective_sub_locant,
    plan_from_chain,
)
from namepredict.layer4.polyene import (
    ene_locants, orient_alkenol, orient_alkenedioic, orient_cyclopolyene, orient_polyene,
)
from namepredict.layer4.sat_hetero_orient import orient_sat_hetero_repl as _orient_sat_hetero_repl
def _pos_on(chain: list[int], c: int | None) -> int | None:
    if c is None or c not in chain:
        return None
    return chain.index(c) + 1
def _maybe_reverse(chain: list[int], pos: int) -> list[int]:
    if pos > (len(chain) + 1) // 2:
        return list(reversed(chain))
    return chain
def _locants_on(chain: list[int], substituents: list) -> list[int]:
    return sorted(chain.index(s["attach_idx"]) + 1 for s in substituents)
def _locant_key(locs: list[int]) -> tuple:
    return (locs, len(locs))
def _stem_loc_pairs(chain: list[int], substituents: list) -> list[tuple]:
    return sorted(
        (alkyl_alpha_key(s.get("en") or ""), chain.index(s["attach_idx"]) + 1)
        for s in substituents
    )
def _orient_key(chain: list[int], substituents: list) -> tuple:
    return (_locant_key(_locants_on(chain, substituents)), _stem_loc_pairs(chain, substituents))
def _prefer_chain(a: list[int], b: list[int], substituents: list) -> list[int]:
    return a if _orient_key(a, substituents) <= _orient_key(b, substituents) else b
def _orient_alkane(chain: list[int], substituents: list) -> list[int]:
    if not substituents or not chain:
        return chain
    return _prefer_chain(chain, list(reversed(chain)), substituents)
def _orient_by_single_fg(
    chain: list[int], parent: dict, substituents: list, key: str
) -> list[int]:
    pos = _pos_on(chain, parent.get(key))
    if pos is None:
        return chain
    base = _maybe_reverse(chain, pos)
    rev = list(reversed(base))
    b, r = _pos_on(base, parent.get(key)), _pos_on(rev, parent.get(key))
    if b is not None and r is not None and b == r:
        return _prefer_chain(base, rev, substituents)
    return base
def _orient_alcohol(chain: list[int], parent: dict, substituents: list) -> list[int]:
    if parent.get("double_bonds"): return orient_alkenol(chain, parent, substituents, _prefer_chain)
    return _orient_by_single_fg(chain, parent, substituents, "oh_c_idx")
def _orient_thiol(chain: list[int], parent: dict, substituents: list) -> list[int]:
    return _orient_by_single_fg(chain, parent, substituents, "sh_c_idx")
def _orient_amine(chain: list[int], parent: dict, substituents: list) -> list[int]:
    return _orient_by_single_fg(chain, parent, substituents, "amine_c_idx")
def _orient_ketone(chain: list[int], parent: dict, substituents: list) -> list[int]:
    return _orient_by_single_fg(chain, parent, substituents, "ketone_c_idx")
def _pair_locs_on(chain: list[int], cs) -> tuple[int, ...] | None:
    if not cs: return None
    locs = sorted(chain.index(c) + 1 for c in cs if c in chain)
    return tuple(locs) if len(locs) == len(cs) else None
def _better_pair_orient(a: list[int], b: list[int], cs, subs: list) -> list[int]:
    la, lb = _pair_locs_on(a, cs), _pair_locs_on(b, cs)
    if la is None: return b
    if lb is None or la < lb: return a
    if lb < la: return b
    return _prefer_chain(a, b, subs)
def _orient_pair(chain: list[int], parent: dict, key: str, subs: list) -> list[int]:
    cs = parent.get(key)
    return chain if not cs else _better_pair_orient(chain, list(reversed(chain)), cs, subs)
def _orient_polyol(chain: list[int], parent: dict, substituents: list) -> list[int]:
    return _orient_pair(chain, parent, "oh_c_idxs", substituents)
def _orient_diamine(chain: list[int], parent: dict, substituents: list) -> list[int]:
    return _orient_pair(chain, parent, "amine_c_idxs", substituents)
def _orient_dione(chain: list[int], parent: dict, substituents: list) -> list[int]:
    return _orient_pair(chain, parent, "ketone_c_idxs", substituents)
def _orient_to_terminal(chain: list[int], c_idx: int | None) -> list[int]:
    pos = _pos_on(chain, c_idx)
    return chain if pos is None or pos == 1 else list(reversed(chain))
def _orient_acid(chain: list[int], parent: dict) -> list[int]:
    return _orient_to_terminal(chain, parent.get("cooh_c_idx"))
def _orient_aldehyde(chain: list[int], parent: dict) -> list[int]:
    return _orient_to_terminal(chain, parent.get("aldehyde_c_idx"))
def _orient_ester(chain: list[int], parent: dict) -> list[int]:
    return _orient_to_terminal(chain, parent.get("ester_c_idx"))
def _orient_amide(chain: list[int], parent: dict) -> list[int]:
    return _orient_to_terminal(chain, parent.get("amide_c_idx"))
def _orient_nitrile(chain: list[int], parent: dict) -> list[int]:
    return _orient_to_terminal(chain, parent.get("nitrile_c_idx"))
def _orient_acyl_chloride(chain: list[int], parent: dict) -> list[int]:
    return _orient_to_terminal(chain, parent.get("acyl_c_idx"))
def _ene_ends_on(chain: list[int], ends: tuple[int, int] | None) -> tuple[int, int] | None:
    if not ends or ends[0] not in chain or ends[1] not in chain: return None
    return chain.index(ends[0]) + 1, chain.index(ends[1]) + 1
def _ene_locant_of(ends: tuple[int, int] | None) -> int | None:
    return None if ends is None else min(ends)
def _prefer_ene_orient(a, b, ends0, subs):
    """Pick chain direction with lower unsaturation locant (P-31.1)."""
    la = _ene_locant_of(_ene_ends_on(a, ends0))
    lb = _ene_locant_of(_ene_ends_on(b, ends0))
    if la is None: return b
    if lb is None or la < lb: return a
    return b if lb < la else _prefer_chain(a, b, subs)
def _orient_by_bond(chain, parent, substituents, key):
    ends0 = parent.get(key)
    if _ene_ends_on(chain, ends0) is None: return chain
    return _prefer_ene_orient(chain, list(reversed(chain)), ends0, substituents)
def _orient_alkene(chain: list[int], parent: dict, substituents: list) -> list[int]:
    return _orient_by_bond(chain, parent, substituents, "double_bond")
def _orient_cycloalkene(chain: list[int], parent: dict, substituents: list) -> list[int]:
    ends = parent.get("double_bond")
    if not ends or ends[0] not in chain: return chain
    base = _rotate_to_front(chain, ends[0])
    if len(base) > 1 and base[1] != ends[1]:
        base = _rotate_to_front(list(reversed(chain)), ends[0])
    return base
def _orient_alkyne(chain: list[int], parent: dict, substituents: list) -> list[int]:
    return _orient_by_bond(chain, parent, substituents, "triple_bond")
def _term_fn(fn):
    return lambda c, p, s: fn(c, p)
def _terminal_orienters() -> dict:
    return {
        "acid": _term_fn(_orient_acid), "aldehyde": _term_fn(_orient_aldehyde),
        "ester": _term_fn(_orient_ester), "amide": _term_fn(_orient_amide),
        "nitrile": _term_fn(_orient_nitrile),
        "acyl_chloride": _term_fn(_orient_acyl_chloride),
    }
def _carbonyl_orienters() -> dict:
    return {**_terminal_orienters(), "ketone": _orient_ketone, "dione": _orient_dione}
def _rotate_to_front(chain: list[int], atom: int) -> list[int]:
    if atom not in chain: return chain
    i = chain.index(atom)
    return chain[i:] + chain[:i]
def _rotations(chain: list[int]) -> list[list[int]]:
    return [chain[i:] + chain[:i] for i in range(len(chain))]
def _ring_candidates(chain: list[int]) -> list[list[int]]:
    out: list[list[int]] = []
    for base in (chain, list(reversed(chain))):
        out.extend(_rotations(base))
    return out
def _best_ring(chain: list[int], substituents: list) -> list[int]:
    best = chain
    for cand in _ring_candidates(chain):
        best = _prefer_chain(best, cand, substituents)
    return best
def _orient_cycloalkane(chain: list[int], parent: dict, substituents: list) -> list[int]:
    return chain if not chain or not substituents else _best_ring(chain, substituents)
def _orient_ring_fixed(
    chain: list[int], parent: dict, substituents: list, key: str
) -> list[int]:
    c = parent.get(key)
    if c is None or c not in chain:
        return chain
    base = _rotate_to_front(chain, c)
    rev = _rotate_to_front(list(reversed(chain)), c)
    return _prefer_chain(base, rev, substituents)
def _orient_cycloalcohol(chain: list[int], parent: dict, substituents: list) -> list[int]:
    return _orient_ring_fixed(chain, parent, substituents, "oh_c_idx")
def _orient_boronic(chain: list[int], parent: dict, substituents: list) -> list[int]:
    return _orient_ring_fixed(chain, parent, substituents, "c_attach")
def _better_ring_pair(best, best_locs, cand, cs, subs):
    locs = _pair_locs_on(cand, cs)
    if locs is None:
        return best, best_locs
    if best_locs is None or locs < best_locs:
        return cand, locs
    if locs == best_locs:
        return _prefer_chain(best, cand, subs), best_locs
    return best, best_locs
def _orient_ring_pair(chain: list[int], parent: dict, key: str, subs: list) -> list[int]:
    cs = parent.get(key)
    if not cs or not chain:
        return chain
    best, best_locs = chain, _pair_locs_on(chain, cs)
    for cand in _ring_candidates(chain):
        best, best_locs = _better_ring_pair(best, best_locs, cand, cs, subs)
    return best
def _orient_benzenediol(chain: list[int], parent: dict, substituents: list) -> list[int]:
    return _orient_ring_pair(chain, parent, "oh_c_idxs", substituents)
def _orient_benzenediamine(chain: list[int], parent: dict, substituents: list) -> list[int]:
    return _orient_ring_pair(chain, parent, "amine_c_idxs", substituents)
def _orient_diazine(chain: list[int], parent: dict, substituents: list) -> list[int]:
    return _orient_ring_pair(chain, parent, "n_idxs", substituents)
def _orient_pyrimidinamine(chain: list[int], parent: dict, substituents: list) -> list[int]:
    attach, virtual = parent.get("amine_c_idx"), list(substituents)
    if attach is not None:
        virtual = virtual + [{"attach_idx": attach, "en": ""}]
    return _orient_diazine(chain, parent, virtual)
def _orient_cycloketone(chain: list[int], parent: dict, substituents: list) -> list[int]:
    return _orient_ring_fixed(chain, parent, substituents, "ketone_c_idx")
def _orient_cycloamine(chain: list[int], parent: dict, substituents: list) -> list[int]:
    return _orient_ring_fixed(chain, parent, substituents, "amine_c_idx")
def _orient_benzoic(chain: list[int], parent: dict, substituents: list) -> list[int]:
    return _orient_ring_fixed(chain, parent, substituents, "ring_attach_idx")
def _orient_pyridine(chain: list[int], parent: dict, substituents: list) -> list[int]:
    return _orient_ring_fixed(chain, parent, substituents, "n_idx")
def _orient_hetero5(chain: list[int], parent: dict, substituents: list) -> list[int]:
    return _orient_ring_fixed(chain, parent, substituents, "hetero_idx")
def _orient_sat_hetero(chain: list[int], parent: dict, substituents: list) -> list[int]:
    """Mono hetero = 1; di-hetero use pair locs (O before N for morpholine)."""
    hs = parent.get("hetero_idxs") or []
    if len(hs) >= 2:
        return _orient_ring_pair(chain, parent, "hetero_idxs", substituents)
    return _orient_ring_fixed(chain, parent, substituents, "hetero_idx")
def _n_loc_on(chain: list[int], n: int | None) -> int:
    return chain.index(n) + 1 if n is not None and n in chain else 99
def _orient_imidazole(chain: list[int], parent: dict, substituents: list) -> list[int]:
    """NH = 1; choose direction so the other N has the lowest locant (3 vs 4; 2 vs 5)."""
    nh = parent.get("nh_idx")
    if nh is None or nh not in chain:
        return chain
    base = _rotate_to_front(chain, nh)
    rev = _rotate_to_front(list(reversed(chain)), nh)
    n = parent.get("n_idx")
    return base if _n_loc_on(base, n) <= _n_loc_on(rev, n) else rev
def _virtual_cooh_subs(parent: dict, substituents: list) -> list:
    attach, virtual = parent.get("ring_attach_idx"), list(substituents)
    if attach is not None:
        virtual = virtual + [{"attach_idx": attach, "en": ""}]
    return virtual
def _orient_pyridinecarboxylic(chain: list[int], parent: dict, substituents: list) -> list[int]:
    return _orient_pyridine(chain, parent, _virtual_cooh_subs(parent, substituents))
def _orient_hetero5carboxylic(chain: list[int], parent: dict, substituents: list) -> list[int]:
    """Hetero fixed as 1; COOH attach as virtual sub for direction."""
    return _orient_hetero5(chain, parent, _virtual_cooh_subs(parent, substituents))
def _orient_fixed_hetero(chain, parent, subs, h):
    return _orient_ring_fixed(chain, {**parent, "hetero_idx": h}, subs, "hetero_idx")
def _orient_sat_hetero_carboxylic(chain: list[int], parent: dict, substituents: list) -> list[int]:
    """Hetero=1; asym dihetero fix first (O); sym try both for lowest COOH."""
    hs, virt = parent.get("hetero_idxs") or [], _virtual_cooh_subs(parent, substituents)
    if len(hs) < 2 or parent.get("hetero_asym"):
        return _orient_fixed_hetero(chain, parent, virt, hs[0] if hs else None)
    a = _orient_fixed_hetero(chain, parent, virt, hs[0])
    b = _orient_fixed_hetero(chain, parent, virt, hs[1])
    return _prefer_chain(a, b, virt)
def _orient_diazolecarboxylic(chain: list[int], parent: dict, substituents: list) -> list[int]:
    return _orient_imidazole(chain, parent, substituents)
def _orient_pyridin_fg(chain: list[int], parent: dict, substituents: list, key: str) -> list[int]:
    attach, virtual = parent.get(key), list(substituents)
    if attach is not None:
        virtual = virtual + [{"attach_idx": attach, "en": ""}]
    return _orient_pyridine(chain, parent, virtual)
def _orient_pyridinamine(chain: list[int], parent: dict, substituents: list) -> list[int]:
    return _orient_pyridin_fg(chain, parent, substituents, "amine_c_idx")
def _orient_pyridinol(chain: list[int], parent: dict, substituents: list) -> list[int]:
    return _orient_pyridin_fg(chain, parent, substituents, "oh_c_idx")
def _naph_loc_on(chain: list[int], attach: int) -> int:
    if attach not in chain: return 99
    loc = _NAPH_LOCANTS[chain.index(attach)]
    return 99 if loc is None else loc
def _naph_loc_key(chain: list[int], substituents: list) -> tuple:
    locs = sorted(_naph_loc_on(chain, s["attach_idx"]) for s in substituents)
    return tuple(locs) if locs else ()
def _pick_naph_chain(cands: list, key_fn) -> list[int]:
    best = cands[0]
    for cand in cands[1:]:
        if key_fn(cand) < key_fn(best): best = cand
    return best
def _orient_naphthalene(chain: list[int], parent: dict, substituents: list) -> list[int]:
    cands = parent.get("naph_chains") or [chain]
    if not cands: return chain
    if not substituents: return cands[0]
    return _pick_naph_chain(cands, lambda c: _naph_loc_key(c, substituents))
def _orient_naphthalenecarboxylic(chain: list[int], parent: dict, substituents: list) -> list[int]:
    """Lowest sub set, then lowest principal COOH locant (P-14.4)."""
    cands = parent.get("naph_chains") or [chain]
    if not cands: return chain
    virt = _virtual_cooh_subs(parent, substituents)
    attach = parent.get("ring_attach_idx")
    def key(c):
        cooh = 99 if attach is None else _naph_loc_on(c, attach)
        return (_naph_loc_key(c, virt), cooh)
    return _pick_naph_chain(cands, key)
def _orient_indole(chain: list[int], parent: dict, substituents: list) -> list[int]:
    return chain
_SAT_HETERO_PLAIN = (
    "aziridine", "oxirane", "oxolane", "oxane", "pyrrolidine", "piperidine",
    "morpholine", "piperazine", "dioxolane", "dioxane", "thiolane",
)
_SAT_HETERO_COOH = (
    "piperidinecarboxylic", "pyrrolidinecarboxylic", "piperazinecarboxylic",
    "morpholinecarboxylic", "oxolanecarboxylic", "oxanecarboxylic",
    "thiolanecarboxylic", "aziridinecarboxylic",
)
def _sat_hetero_orienters() -> dict:
    d = {k: _orient_sat_hetero for k in _SAT_HETERO_PLAIN}
    d.update({k: _orient_sat_hetero_carboxylic for k in _SAT_HETERO_COOH})
    return {**d, "sat_hetero_repl": _orient_sat_hetero_repl}
def _aza_orienters() -> dict:
    return {
        "pyridine": _orient_pyridine, "pyridinecarboxylic": _orient_pyridinecarboxylic,
        "pyridinecarbonitrile": _orient_pyridinecarboxylic,
        "pyridinamine": _orient_pyridinamine, "pyridinol": _orient_pyridinol,
        "imidazole": _orient_imidazole, "pyrazole": _orient_imidazole,
        "oxazole": _orient_imidazole, "thiazole": _orient_imidazole,
        "imidazolecarboxylic": _orient_diazolecarboxylic,
        "pyrazolecarboxylic": _orient_diazolecarboxylic,
        "pyrimidine": _orient_diazine, "pyrazine": _orient_diazine,
        "pyridazine": _orient_diazine, "pyrimidinamine": _orient_pyrimidinamine}
def _fused_orienters() -> dict:
    d = {k: _orient_indole for k in _INDOLE_ORIENT_KINDS}
    d["naphthalene"] = _orient_naphthalene
    d["naphthalenecarboxylic"] = _orient_naphthalenecarboxylic
    return d
def _arene_orienters() -> dict:
    return {
        "furan": _orient_hetero5, "thiophene": _orient_hetero5, "pyrrole": _orient_hetero5,
        "furancarboxylic": _orient_hetero5carboxylic,
        "thiophenecarboxylic": _orient_hetero5carboxylic,
        "pyrrolecarboxylic": _orient_hetero5carboxylic,
        **_fused_orienters(), **_aza_orienters(), **_sat_hetero_orienters(),
    }
def _hetero_orienters() -> dict:
    return {
        "alcohol": _orient_alcohol, "thiol": _orient_thiol,
        "diol": _orient_polyol, "triol": _orient_polyol, "diamine": _orient_diamine,
        "cycloalcohol": _orient_cycloalcohol, "phenol": _orient_cycloalcohol,
        "boronic": _orient_boronic, "benzenediol": _orient_benzenediol,
        "benzenediamine": _orient_benzenediamine,
        "cycloamine": _orient_cycloamine, "aniline": _orient_cycloamine,
        "amine": _orient_amine, "sec_amine": _orient_amine, "tert_amine": _orient_amine,
        **_arene_orienters(),
    }
def _orient_polyene(chain: list[int], parent: dict, substituents: list) -> list[int]:
    return orient_polyene(chain, parent, substituents, _prefer_chain)
def _orient_alkenedioic(chain: list[int], parent: dict, substituents: list) -> list[int]:
    return orient_alkenedioic(chain, parent, substituents, _prefer_chain, _orient_alkene)
def _orient_diacid(chain: list[int], parent: dict, substituents: list) -> list[int]:
    if parent.get("double_bond") or parent.get("double_bonds"):
        return _orient_alkenedioic(chain, parent, substituents)
    return _orient_alkane(chain, substituents)
def _unsat_orienters() -> dict:
    b = _orient_benzoic
    return {
        "cycloketone": _orient_cycloketone, "alkene": _orient_alkene,
        "diacid": _orient_diacid, "diester": _orient_diacid, "polyene": _orient_polyene,
        "cyclopolyene": orient_cyclopolyene, "cycloalkene": _orient_cycloalkene,
        "alkyne": _orient_alkyne, "cycloalkane": _orient_cycloalkane,
        "benzene": _orient_cycloalkane, "benzoic": b, "benzaldehyde": b,
        "acetophenone": b, "benzoate": b, "benzonitrile": b,
        "benzoyl_chloride": b, "cycloalkanecarboxylic": b,
    }
def _kind_orienters() -> dict:
    return {**_hetero_orienters(), **_unsat_orienters(), **_carbonyl_orienters()}
def _orient_by_kind(kind: str, chain: list[int], parent: dict, subs: list) -> list[int]:
    fn = _kind_orienters().get(kind)
    return _orient_alkane(chain, subs) if fn is None else fn(chain, parent, subs)
def _orient_chain(parent: dict, substituents: list) -> list[int]:
    chain = list(parent.get("chain") or [])
    if not chain:
        return chain
    return _orient_by_kind(parent.get("kind"), chain, parent, substituents)
def _fg_locant(oriented: dict, kinds: tuple, key: str) -> int | None:
    if oriented.get("kind") not in kinds:
        return None
    return _pos_on(oriented.get("chain") or [], oriented.get(key))
def _oh_locant(oriented: dict) -> int | None:
    kind, chain, a = oriented.get("kind"), oriented.get("chain") or [], oriented.get("oh_c_idx")
    if kind == "quinolinol" and a is not None:
        return _naph_sub_locant(chain, a)
    if kind == "benzothiophenol" and a is not None:
        return _indole_sub_locant(chain, a)
    return _fg_locant(
        oriented, ("alcohol", "cycloalcohol", "pyridinol"), "oh_c_idx",
    )
def _sh_locant(oriented: dict) -> int | None:
    return _fg_locant(oriented, ("thiol",), "sh_c_idx")
def _pair_locants(oriented: dict, kinds, key: str) -> list[int] | None:
    if oriented.get("kind") not in kinds:
        return None
    locs = _pair_locs_on(oriented.get("chain") or [], oriented.get(key))
    return list(locs) if locs else None
def _oh_locants(oriented: dict) -> list[int] | None:
    return _pair_locants(oriented, ("diol", "triol", "benzenediol"), "oh_c_idxs")
def _amine_pair_locants(oriented: dict) -> list[int] | None:
    return _pair_locants(oriented, ("diamine", "benzenediamine"), "amine_c_idxs")
_FUSED_AM = (
    "benzofuranamine", "benzothiazolamine", "benzoxazolamine", "benzimidazolamine",
)
_AMINE_KINDS = (
    "amine", "cycloamine", "sec_amine", "tert_amine", "pyridinamine",
    "pyrimidinamine",
) + _FUSED_AM
def _amine_locant(oriented: dict) -> int | None:
    if oriented.get("kind") in _FUSED_AM:
        chain, a = oriented.get("chain") or [], oriented.get("amine_c_idx")
        return _indole_sub_locant(chain, a) if a is not None else None
    return _fg_locant(oriented, _AMINE_KINDS, "amine_c_idx")
def _ketone_locant(oriented: dict) -> int | None:
    return _fg_locant(oriented, ("ketone", "cycloketone"), "ketone_c_idx")
def _ketone_pair_locants(oriented: dict) -> list[int] | None:
    return _pair_locants(oriented, ("dione",), "ketone_c_idxs")
def _has_parent_ene(oriented: dict) -> bool:
    return bool(oriented.get("double_bond") or oriented.get("double_bonds"))
def _ene_locant(oriented: dict) -> int | None:
    kind = oriented.get("kind")
    if kind in ("alkene", "cycloalkene") or _has_parent_ene(oriented):
        return _ene_locant_of(
            _ene_ends_on(oriented.get("chain") or [], oriented.get("double_bond"))
        )
    return None
def _yne_locant(oriented: dict) -> int | None:
    if oriented.get("kind") == "alkyne" or oriented.get("triple_bond"):
        return _ene_locant_of(_ene_ends_on(oriented.get("chain") or [], oriented.get("triple_bond")))
    return None
def _has_parent_yne(oriented: dict) -> bool:
    return bool(oriented.get("triple_bond"))
def _omit_oh(oh_pos: int | None, n_carbons: int, kind: str | None = None, parent: dict | None = None) -> bool:
    if kind == "cycloalcohol":
        return True
    if kind == "alcohol" and parent and (_has_parent_ene(parent) or _has_parent_yne(parent)):
        return False
    return oh_pos == 1 and n_carbons <= 2
def _omit_sh(sh_pos: int | None, n_carbons: int) -> bool:
    return sh_pos == 1 and n_carbons <= 2
def _omit_amine(am_pos: int | None, n_carbons: int, kind: str | None = None) -> bool:
    if kind == "cycloamine":
        return True
    return am_pos == 1 and n_carbons <= 2
def _omit_unsat(n_carbons: int, kind: str | None = None, parent: dict | None = None) -> bool:
    if kind == "cycloalkene":
        return True
    if kind == "alcohol" and parent and _has_parent_yne(parent):
        return False
    if parent and _has_parent_ene(parent) and kind != "alkene":
        return False
    return n_carbons <= 3
def _naph_sub_locant(chain: list[int], attach: int) -> int:
    if attach not in chain or len(chain) != 10:
        return chain.index(attach) + 1 if attach in chain else 0
    loc = _NAPH_LOCANTS[chain.index(attach)]
    return loc if loc is not None else chain.index(attach) + 1
def _indole_sub_locant(chain: list[int], attach: int) -> int:
    if attach not in chain or len(chain) != 9:
        return chain.index(attach) + 1 if attach in chain else 0
    loc = _INDOLE_LOCANTS[chain.index(attach)]
    return loc if loc is not None else chain.index(attach) + 1
def _sub_locant(chain: list[int], attach: int, kind: str | None) -> int:
    plan = plan_from_chain(chain, kind)
    loc = effective_sub_locant(plan, attach) if plan else None
    if loc is not None:
        return loc
    if kind in _NAPH_KINDS or kind in _Q_KINDS:
        return _naph_sub_locant(chain, attach)
    if kind in _FUSED56_KINDS:
        return _indole_sub_locant(chain, attach)
    return chain.index(attach) + 1
def _with_locants(chain: list[int], substituents: list, kind: str | None = None) -> list:
    return [{**s, "locant": _sub_locant(chain, s["attach_idx"], kind)} for s in substituents]
def _unsat_locants(oriented: dict, n: int) -> dict:
    kind = oriented.get("kind")
    return {
        "ene_locant": _ene_locant(oriented),
        "ene_locants": ene_locants(oriented),
        "omit_ene_locant": _omit_unsat(n, kind, oriented),
        "yne_locant": _yne_locant(oriented),
        "omit_yne_locant": _omit_unsat(n, kind, oriented),
    }
def _oh_am_locants(oriented: dict, n: int) -> dict:
    oh, am = _oh_locant(oriented), _amine_locant(oriented)
    kind = oriented.get("kind")
    return {
        "oh_locant": oh, "oh_locants": _oh_locants(oriented),
        "omit_oh_locant": _omit_oh(oh, n, kind, oriented),
        "amine_locant": am, "amine_locants": _amine_pair_locants(oriented),
        "omit_amine_locant": _omit_amine(am, n, kind),
    }
def _sh_locants(oriented: dict, n: int) -> dict:
    sh = _sh_locant(oriented)
    return {"sh_locant": sh, "omit_sh_locant": _omit_sh(sh, n)}
def _fg_locants(oriented: dict) -> dict:
    n = oriented.get("n_carbons", 0)
    base = {
        **_oh_am_locants(oriented, n), **_sh_locants(oriented, n),
        "ketone_locant": _ketone_locant(oriented),
        "ketone_locants": _ketone_pair_locants(oriented),
    }
    return {**base, **_unsat_locants(oriented, n)}
def _pack(oriented: dict, substituents: list) -> dict:
    base = {"parent": oriented, "substituents": substituents}
    return {**base, **_fg_locants(oriented)}
def number(parent: dict, substituents: list) -> dict:
    chain, kind = _orient_chain(parent, substituents), parent.get("kind")
    oriented = {**parent, "chain": chain}
    plan = plan_from_chain(chain, kind)
    if plan is not None:
        oriented["numbering"] = plan
    return _pack(oriented, _with_locants(chain, substituents, kind))
