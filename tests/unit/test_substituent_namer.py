"""Ordered SubstituentNamer backends: retained → rooted_tree → recursive."""
from __future__ import annotations

from namepredict.layer3.claimable_block import ClaimedBlock, SideSlot
from namepredict.layer3.substituent_namer import SubstituentName, SubstituentNamer


def _claim(atoms=(1, 2, 3), root=1, attach=0) -> ClaimedBlock:
    return ClaimedBlock(
        slot=SideSlot.OTHER,
        attach_parent=attach,
        root=root,
        atoms=frozenset(atoms),
    )


class _FakeBackend:
    def __init__(self, name: str, result: SubstituentName | None, log: list[str]):
        self.name = name
        self._result = result
        self._log = log

    def try_name(self, mol, claim: ClaimedBlock, *, depth: int) -> SubstituentName | None:
        self._log.append(self.name)
        return self._result


def _ok(claim: ClaimedBlock, backend: str) -> SubstituentName:
    return SubstituentName(
        claim=claim,
        en=f"{backend}-en",
        zh=f"{backend}-zh",
        requires_parentheses=False,
        backend=backend,
    )


def test_retained_success_skips_later_backends():
    claim = _claim()
    log: list[str] = []
    namer = SubstituentNamer(
        backends=[
            _FakeBackend("retained", _ok(claim, "retained"), log),
            _FakeBackend("rooted_tree", _ok(claim, "rooted_tree"), log),
            _FakeBackend("recursive", _ok(claim, "recursive"), log),
        ]
    )
    result = namer.name(None, claim, depth=0)
    assert result is not None
    assert result.backend == "retained"
    assert result.claim is claim
    assert log == ["retained"]


def test_rooted_tree_runs_after_retained_failure():
    claim = _claim()
    log: list[str] = []
    namer = SubstituentNamer(
        backends=[
            _FakeBackend("retained", None, log),
            _FakeBackend("rooted_tree", _ok(claim, "rooted_tree"), log),
            _FakeBackend("recursive", _ok(claim, "recursive"), log),
        ]
    )
    result = namer.name(None, claim, depth=0)
    assert result is not None
    assert result.backend == "rooted_tree"
    assert log == ["retained", "rooted_tree"]


def test_all_backends_fail_returns_none():
    claim = _claim()
    log: list[str] = []
    namer = SubstituentNamer(
        backends=[
            _FakeBackend("retained", None, log),
            _FakeBackend("rooted_tree", None, log),
            _FakeBackend("recursive", None, log),
        ]
    )
    assert namer.name(None, claim, depth=0) is None
    assert log == ["retained", "rooted_tree", "recursive"]


def test_backend_none_never_yields_empty_name_or_mutates_claim():
    claim = _claim(atoms=(10, 11))
    before = (claim.slot, claim.attach_parent, claim.root, claim.atoms)
    namer = SubstituentNamer(
        backends=[_FakeBackend("retained", None, []), _FakeBackend("rooted_tree", None, [])]
    )
    result = namer.name(None, claim, depth=0)
    assert result is None
    assert (claim.slot, claim.attach_parent, claim.root, claim.atoms) == before
