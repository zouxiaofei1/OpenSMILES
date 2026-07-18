# Universal Claimable-Block Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make every selected parent either account for every heavy atom through typed, ordered side-block naming or fail/retry another parent; deliver this as a class-level shared mechanism rather than molecule-specific naming paths.

**Architecture:** L2 claims only ownership, boundary, and attachment: it never selects a naming mode or emits a name. The selected parent owns an immutable terminal `owned_atoms` set; outside connected components become typed `ClaimedBlock`s. L3's ordered `SubstituentNamer` tries retained, rooted-tree, then bounded recursive backends and returns a typed `SubstituentName`; L4/L5 assemble only successful names. `CoverageLedger` is the final acceptance gate: a candidate is usable only when it has neither gaps nor overlaps, otherwise parent selection retries the next candidate.

**Tech Stack:** Python 3.12, RDKit, pytest, `SMILESNNamer` dual output, `tools/structure_lint.py`, `benchmarks.benchmark_parallel`.

## Global Constraints

- Code belongs only in `src/namepredict/layer0` through `layer5`, `src/namepredict/namer.py`, and `tests/unit`; `data/*` is read-only.
- Keep source files at or below 500 lines and executable function bodies at or below 10 non-comment lines; split by responsibility before exceeding either limit.
- L2 may determine parent ownership, boundary, attachment, and topology only. L3 names substituents. L4/L5 numbers and assembles; no molecule-specific output-string patch is permitted.
- Do not add deep learning, unbounded caches, or benchmark-gold edits.
- TDD is mandatory: red test, observed failure, minimum implementation, green test, then commit. One independently reviewable class slice per task.
- Every naming change runs its affected unit tests; the final gate runs `python -m benchmarks.benchmark_parallel --data data/merged_benchmark.json --time`, with dual fallback regression no greater than 0.5%.
- IUPAC references: `docs/iupac/` P-29 (substituents), P-14 (locants/order), P-63.2 (ethers), and P-66.1 (amide N substituents).
- A recognized external block must either be named and covered or make the current candidate fail. It must never be silently omitted.
- No parent-specific copy of the recursive substituent engine is allowed in `benzamide.py`, ring producers, or chain producers.
- `parent_atom_set(parent, mol)` remains only a compatibility adapter. The terminal truth is `parent["owned_atoms"]`, created once per selected candidate and never recomputed from a later partial representation.
- Do not introduce or retain `mode="systematic_alkyl"`. Alkyl naming is one `SubstituentNamer` backend, not a claim mode or L2 concern.
- Under Windows Git Bash run `source .venv/Scripts/activate` and `export PYTHONPATH=src` in a worktree.

---

## Confirmed target model

### Ownership and claims

```python
# src/namepredict/layer2/claimable_block.py
from dataclasses import dataclass
from enum import Enum

class SideSlot(str, Enum):
    CHAIN_C = "chain_c"
    RING_C = "ring_c"
    AMIDE_N = "amide_n"
    ETHER_O = "ether_o"
    OTHER = "other"

@dataclass(frozen=True)
class ClaimedBlock:
    slot: SideSlot
    attach_parent: int
    root: int
    atoms: frozenset[int]


def claim_block(mol, *, owned_atoms: frozenset[int], attach_parent: int,
                root: int, slot: SideSlot) -> ClaimedBlock | None:
    """Return the complete outside heavy component and its attachment, or None."""


def iter_claims(mol, owned_atoms: frozenset[int]) -> list[ClaimedBlock]:
    """Return canonical-order claims for all outside heavy components."""
```

`claim_block` verifies only that `attach_parent` is owned, `root` is its outside heavy neighbour, and the full component reached without entering `owned_atoms` is nonempty. It performs no leaf recognition, no alkyl classification, no cut-submol naming, and has no `mode` field. A component touching more than one owned atom is returned as one claim with the canonical attachment only when the topology policy says it is a valid substituent; otherwise it returns `None` and invalidates the candidate.

### Typed L3 naming contract

```python
# src/namepredict/layer3/substituent_namer.py
from dataclasses import dataclass
from typing import Protocol
from namepredict.layer2.claimable_block import ClaimedBlock

@dataclass(frozen=True)
class SubstituentName:
    claim: ClaimedBlock
    en: str
    zh: str
    requires_parentheses: bool
    backend: str  # "retained", "rooted_tree", or "recursive"

class SubstituentBackend(Protocol):
    name: str
    def try_name(self, mol, claim: ClaimedBlock, *, depth: int) -> SubstituentName | None: ...

class SubstituentNamer:
    def name(self, mol, claim: ClaimedBlock, *, depth: int = 0) -> SubstituentName | None: ...
```

`SubstituentNamer` tries backends in this fixed order: retained names, saturated-carbon rooted tree, then bounded recursive `cut_block → name_as_substituent → yl_form`. The rooted-tree backend accepts only its explicit pure saturated carbon domain (`max_atoms=12`, `max_depth=3`); it is not a globally visible mode. All successful output carries exactly the source `ClaimedBlock`, including source atoms and attachment. Recursive failure returns `None`, never an empty prefix.

### Final coverage gate

```python
# src/namepredict/layer3/coverage.py
from dataclasses import dataclass

@dataclass(frozen=True)
class CoverageLedger:
    owned_atoms: frozenset[int]
    named_claims: tuple[SubstituentName, ...]
    gap: frozenset[int]
    overlap: frozenset[int]

    @property
    def complete(self) -> bool:
        return not self.gap and not self.overlap


def build_coverage_ledger(mol, *, owned_atoms: frozenset[int],
                          names: list[SubstituentName]) -> CoverageLedger: ...
```

For heavy atoms, `gap` is `all_heavy - owned_atoms - union(claim.atoms)` and `overlap` contains an atom appearing in more than one of `owned_atoms` and named claims. The ledger is evaluated after all claims are named and before assembly. A non-complete ledger rejects the candidate, including a parent that otherwise generated a plausible partial string.

### Candidate retry boundary

`select_parent` exposes ranked candidates, not only an eagerly accepted dict. The naming coordinator attempts candidates in existing rank order: finalize immutable `owned_atoms`, claim every external component, name all claims, build the ledger, then orient/assemble only if complete. Any claim, backend, ledger, orientation, or assembly failure moves to the next candidate. Public `SMILESNNamer.name(smiles)` remains unchanged and returns failure only after candidates are exhausted.

---

## File structure and migration map

| Path | Responsibility |
|---|---|
| `src/namepredict/layer2/block_cut.py` | Existing cut/submol primitives; consumes terminal `owned_atoms`, never derives naming policy. |
| `src/namepredict/layer2/claimable_block.py` | `SideSlot`, `ClaimedBlock`, complete-component claims, deterministic ordering. |
| `src/namepredict/layer2/parent_core.py` or new focused `parent_ownership.py` | Materialize `owned_atoms` for every candidate kind; retain `parent_atom_set` adapter. |
| `src/namepredict/layer2/parent_selector.py` | Return ordered candidates and preserve rank/retry metadata. |
| `src/namepredict/layer3/substituent_namer.py` | Typed orchestrator and backend protocol/order. |
| `src/namepredict/layer3/substituent_retained.py` | Existing true simple leaves behind the backend protocol. |
| `src/namepredict/layer2/side_alkyl_sys.py` and `src/namepredict/layer3/alkyl_sys_names.py` | Rooted saturated-carbon topology and dual naming backend. |
| `src/namepredict/layer3/as_substituent.py`, `yl_form.py` | Recursive backend only; bounded recursion and attachment-aware `-yl`. |
| `src/namepredict/layer3/coverage.py` | Final heavy-atom gap/overlap ledger. |
| `src/namepredict/layer3/substituent_extractor.py` | Convert `SubstituentName` to L3 dictionaries; remove parallel extraction paths. |
| `src/namepredict/namer.py` | Candidate attempt/retry coordinator and final gate. |
| `src/namepredict/layer2/alkenamide.py`, `aryl_sub.py`, `benzamide.py` | Remove false shortcuts and consume shared claims only. |

Delete, rather than preserve as fallbacks, legacy side paths that independently accept unowned atoms or create `n_alkyl`, `n_phenyl`, `_SIDE_PROBES`, `claimable_side`, or parent-local recursive names outside the typed orchestrator. Keep simple-name presentation adapters only when they construct a `SubstituentName` with the original claim. Search and remove all references to `systematic_alkyl` mode before declaring migration complete.

---

### Task 1: Establish class-level red bars and migration inventory

**Files:**
- Create: `tests/unit/test_universal_claimable_block.py`
- Create: `tests/unit/test_coverage_ledger.py`
- Modify: `docs/superpowers/plans/2026-07-18-universal-claimable-block.md` only if observed names require correcting this plan's test fixtures before implementation

**Consumes:** public `SMILESNNamer.name`.

**Produces:** failing class fixtures and an exact list of legacy side-path symbols to remove.

- [ ] **Step 1: Write the red-bar tests.**

```python
import pytest
from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    ("O=C(NC1CCCCCC1)c1ccccc1", "N-cycloheptylbenzamide", "N-环庚基苯甲酰胺"),
    ("O=C(Nc1cc(C)cc(C)c1)c1ccccc1", "N-(3,5-dimethylphenyl)benzamide", "N-(3,5-二甲基苯基)苯甲酰胺"),
    ("COC(C)CC", "2-methoxybutane", "2-甲氧基丁烷"),
    ("CCCCC(CC(C)CC)C(=O)O", "2-(2-methylbutyl)hexanoic acid", "2-(2-甲基丁基)己酸"),
    ("c1ccc(cc1)CC(C)CC", "(2-methylbutyl)benzene", "(2-甲基丁基)苯"),
]

@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_universal_claimable_classes(smiles, en, zh):
    result = SMILESNNamer().name(smiles)
    assert result.success
    assert normalize_en(result.en) == normalize_en(en)
    assert normalize_zh(result.zh) == normalize_zh(zh)
```

- [ ] **Step 2: Run the red bar.** Run `PYTHONPATH=src pytest tests/unit/test_universal_claimable_block.py -v`. Expected: at least the documented false-shortcut and evaporation cases fail; existing simple benzamide regressions remain green.
- [ ] **Step 3: Inventory removals.** Run `rg -n "systematic_alkyl|claimable_side|_SIDE_PROBES|n_alkyl|n_phenyl" src/namepredict` and record each producer/consumer in the task branch notes. Classify every hit as retained implementation, adapter, or deletion; no unclassified hit may survive Task 8.
- [ ] **Step 4: Commit.** `git add tests/unit/test_universal_claimable_block.py tests/unit/test_coverage_ledger.py && git commit -m "test(namepredict): add universal claimable-block red bars"`.

### Task 2: Materialize terminal parent ownership and ordered candidates

**Files:**
- Modify: `src/namepredict/layer2/parent_selector.py`
- Modify: `src/namepredict/layer2/block_cut.py`
- Modify: `src/namepredict/layer2/parent_core.py` or create `src/namepredict/layer2/parent_ownership.py`
- Create: `tests/unit/test_parent_ownership.py`

**Consumes:** existing ranked parent recognition facts.

**Produces:** `iter_parent_candidates(info) -> list[dict]`, each with immutable `owned_atoms`; `parent_atom_set` compatibility adapter.

- [ ] **Step 1: Write failing ownership tests.** Assert a benzamide candidate owns aryl core plus amide C/N/O but not N-phenyl atoms; acid owns its carboxyl C and both oxygens; each `owned_atoms` is a `frozenset` and is unchanged after extraction.
- [ ] **Step 2: Run.** `PYTHONPATH=src pytest tests/unit/test_parent_ownership.py -v`; expected import/assertion failures.
- [ ] **Step 3: Implement ownership finalization.** Add `finalize_parent_ownership(parent, mol) -> dict` that copies the candidate once with `owned_atoms=frozenset(...)`. Make `parent_atom_set(parent, mol)` return `parent["owned_atoms"]` when present, otherwise finalize only for legacy direct callers; do not use it as an alternative truth inside the new pipeline.
- [ ] **Step 4: Expose ranked candidates.** Preserve current ordering and make the legacy `select_parent` return the first finalized candidate only for untouched callers. The coordinator added in Task 7 must use `iter_parent_candidates`.
- [ ] **Step 5: Run and commit.** Run `PYTHONPATH=src pytest tests/unit/test_parent_ownership.py tests/unit/test_benzamide.py -v`; commit `feat(namepredict): finalize parent ownership for candidates`.

### Task 3: Add ownership-only claims and typed coverage ledger

**Files:**
- Create: `src/namepredict/layer2/claimable_block.py`
- Create: `src/namepredict/layer3/coverage.py`
- Create: `tests/unit/test_claimable_block_api.py`
- Modify: `tests/unit/test_coverage_ledger.py`

**Consumes:** `parent["owned_atoms"]` from Task 2.

**Produces:** `ClaimedBlock`, `iter_claims`, and `build_coverage_ledger` exactly as defined above.

- [ ] **Step 1: Write failing topology tests.** For N-phenyl benzamide assert one `AMIDE_N` claim with six atoms; for 2-methylbutylbenzene assert one `RING_C` claim containing all five side atoms; assert a two-parent-attachment component is rejected.
- [ ] **Step 2: Write failing ledger tests.** Construct typed test names from disjoint claims and assert complete; omit one atom and assert it is in `gap`; duplicate an atom across two names and assert it is in `overlap`.
- [ ] **Step 3: Run.** `PYTHONPATH=src pytest tests/unit/test_claimable_block_api.py tests/unit/test_coverage_ledger.py -v`; expected import failures.
- [ ] **Step 4: Implement L2 claims.** Derive `SideSlot` exclusively from the owned attachment atom's role. Traverse all heavy neighbours outside ownership, de-duplicate components, and sort by `(attach_parent, root, slot.value)`. Do not import L3 and do not add a mode field.
- [ ] **Step 5: Implement coverage.** Count membership across ownership and named claim atom sets. Include every RDKit heavy atom; hydrogen is excluded. Return immutable sets.
- [ ] **Step 6: Run and commit.** Run both test files plus `python tools/structure_lint.py --root src/namepredict/layer2/claimable_block.py --root src/namepredict/layer3/coverage.py`; commit `feat(namepredict): add ownership-only claims and coverage ledger`.

### Task 4: Implement typed, ordered substituent naming backends

**Files:**
- Create: `src/namepredict/layer3/substituent_namer.py`
- Modify: `src/namepredict/layer3/as_substituent.py`
- Modify: `src/namepredict/layer3/yl_form.py`
- Create: `tests/unit/test_substituent_namer.py`

**Consumes:** `ClaimedBlock`, existing retained probes, `cut_block`, `build_cut_submol`, and `_name_mol`.

**Produces:** `SubstituentName`, `SubstituentBackend`, `SubstituentNamer`; ordered retained/rooted-tree/recursive behavior.

- [ ] **Step 1: Write failing order tests.** Inject three fake backends that append their names to a list. Assert retained success prevents later calls; retained failure then rooted-tree success prevents recursive; all failures return `None`.
- [ ] **Step 2: Write recursive failure test.** A backend returning `None` must produce no empty `en`/`zh` object and must not mutate the source claim.
- [ ] **Step 3: Run.** `PYTHONPATH=src pytest tests/unit/test_substituent_namer.py -v`; expected import failures.
- [ ] **Step 4: Implement the protocol and orchestrator.** Construct each successful `SubstituentName` with the exact input claim. Give the implementations stable backend labels `retained`, `rooted_tree`, and `recursive`; no enum/string named `systematic_alkyl` is permitted.
- [ ] **Step 5: Adapt recursive naming.** `name_as_substituent` is called only by the recursive backend, receives the claim's complete atom set and root attachment, increments depth, and returns `None` at `depth >= 4` or on an unattached/unnumberable result.
- [ ] **Step 6: Run and commit.** Run `PYTHONPATH=src pytest tests/unit/test_substituent_namer.py tests/unit/test_as_substituent.py tests/unit/test_yl_form.py -v`; commit `feat(namepredict): add typed ordered substituent namer`.

### Task 5: Move saturated carbon naming behind the rooted-tree backend

**Files:**
- Create or modify: `src/namepredict/layer2/side_alkyl_sys.py`
- Create or modify: `src/namepredict/layer3/alkyl_sys_names.py`
- Modify: `src/namepredict/layer3/substituent_namer.py`
- Create: `tests/unit/test_rooted_alkyl_backend.py`

**Consumes:** Task 4 backend protocol.

**Produces:** rooted-tree backend for pure saturated carbon claims only.

- [ ] **Step 1: Write red tests.** Test 2-methylbutyl attached to an acid and benzene, plus an unsupported alkene and an over-limit 13-atom tree. Assert supported cases return backend `rooted_tree`; unsupported cases return `None` from this backend so recursion may try next.
- [ ] **Step 2: Run.** `PYTHONPATH=src pytest tests/unit/test_rooted_alkyl_backend.py -v`; expected failures.
- [ ] **Step 3: Implement topology then names.** Build a rooted tree from `claim.root`; accept only single-bond, non-ring, non-aromatic carbon atoms, at most 12 atoms and depth 3. Reuse retained isopentyl/isobutyl presentation only in the retained backend; otherwise generate the systematic dual `-yl` form in L3.
- [ ] **Step 4: Remove direct carbon bypasses.** Replace `_SIDE_PROBES`/extractor direct systematic calls with registration of this backend. Delete the old mode field and all checks for `"systematic_alkyl"`.
- [ ] **Step 5: Run and commit.** Run `PYTHONPATH=src pytest tests/unit/test_rooted_alkyl_backend.py tests/unit/test_universal_claimable_block.py -v`; commit `feat(namepredict): route saturated sides through rooted-tree backend`.

### Task 6: Migrate amide N and open-chain alkoxy to shared names

**Files:**
- Modify: `src/namepredict/layer2/alkenamide.py`
- Modify: `src/namepredict/layer2/aryl_sub.py`
- Modify: `src/namepredict/layer2/benzamide.py`
- Modify: `src/namepredict/layer3/alkoxy_names.py`
- Modify: `src/namepredict/layer3/substituent_extractor.py`
- Create: `tests/unit/test_amide_n_claims.py`
- Create: `tests/unit/test_open_chain_alkoxy_claims.py`

**Consumes:** shared claims/namer from Tasks 3–5.

**Produces:** no false `N-methyl`/`N-phenyl`; N and ether sides use `SubstituentName`.

- [ ] **Step 1: Write red shortcut tests.** Assert cycloheptyl does not set simple `n_alkyl`; 3,5-dimethylphenyl does not set simple `n_phenyl`; their final names are the class fixtures from Task 1.
- [ ] **Step 2: Write ether red tests.** Assert `COC(C)CC` names as 2-methoxybutane and that every named ether side is represented by a typed claim.
- [ ] **Step 3: Run.** `PYTHONPATH=src pytest tests/unit/test_amide_n_claims.py tests/unit/test_open_chain_alkoxy_claims.py -v`; expected failures.
- [ ] **Step 4: Remove false metadata acceptance.** Simple N-alkyl requires an entirely acyclic, saturated, non-aromatic carbon arm; simple phenyl requires no external ring substituent. On failure do not emit an alternative simple meta flag; let the normal claim enumeration create the N claim.
- [ ] **Step 5: Convert extract/assembly input.** Make N-block and alkoxy extraction consume `SubstituentName` and its `claim.atoms`; do not independently recut or rename. Preserve `N-(...)` formatting in L5 based on the typed result.
- [ ] **Step 6: Run and commit.** Run `PYTHONPATH=src pytest tests/unit/test_amide_n_claims.py tests/unit/test_open_chain_alkoxy_claims.py tests/unit/test_benzamide.py -v`; commit `feat(namepredict): migrate N and alkoxy sides to shared claims`.

### Task 7: Gate assembly with coverage and retry candidates

**Files:**
- Modify: `src/namepredict/namer.py`
- Modify: `src/namepredict/layer3/substituent_extractor.py`
- Create: `tests/unit/test_candidate_retry.py`
- Modify: `tests/unit/test_coverage_ledger.py`

**Consumes:** Tasks 2–6.

**Produces:** complete-ledger-only assembly and retry semantics.

- [ ] **Step 1: Write a candidate retry red test.** Monkeypatch ranked candidates so the first owns too little or has an unnameable claim and the second is complete; assert final success uses the second. Add an all-fail case asserting `success is False`, not a partial name.
- [ ] **Step 2: Run.** `PYTHONPATH=src pytest tests/unit/test_candidate_retry.py -v`; expected failure because the first candidate is accepted.
- [ ] **Step 3: Implement `try_candidate`.** Finalize ownership, enumerate claims, name every claim, return failure immediately on a missing name, build a ledger, reject any non-complete ledger, then call numbering and assembly. Do not call L4/L5 before ledger acceptance.
- [ ] **Step 4: Implement retry.** Iterate `iter_parent_candidates` in rank order. Catch only known per-candidate naming failures; do not hide programming errors. Preserve diagnostics containing candidate kind and `gap`/`overlap` for failed attempts.
- [ ] **Step 5: Run and commit.** Run `PYTHONPATH=src pytest tests/unit/test_candidate_retry.py tests/unit/test_coverage_ledger.py tests/unit/test_universal_claimable_block.py -v`; commit `feat(namepredict): retry candidates behind coverage gate`.

### Task 8: Remove old bypasses and add parent-class integration coverage

**Files:**
- Modify: all files identified in Task 1 inventory
- Create: `tests/unit/test_claimable_parent_batches.py`
- Modify: `workstate.md`

**Consumes:** all main-line tasks.

**Produces:** one shared route for listed parent classes; no legacy parallel side route.

- [ ] **Step 1: Write parameterized batch tests.** Cover at least two examples each for acid, alkane, ketone, alcohol, benzene, phenol, aniline, pyridine, amide, and benzamide. For every successful result assert `build_coverage_ledger(...).complete`; include at least one rooted-carbon and one recursive/heteroaryl claim.
- [ ] **Step 2: Run.** `PYTHONPATH=src pytest tests/unit/test_claimable_parent_batches.py -v`; expected parent-specific bypass failures.
- [ ] **Step 3: Delete legacy paths.** Remove mode branching, parent-local recursion, and any assembler path accepting a raw side dict without `SubstituentName` provenance. Replace remaining callers with the typed shared orchestrator. Remove dead imports and tests that assert obsolete modes.
- [ ] **Step 4: Prove deletion.** Run `rg -n "systematic_alkyl|claimable_side|_SIDE_PROBES" src/namepredict`; expected no production hits. Run `rg -n "parent_atom_set\(" src/namepredict`; each new-pipeline hit must be the compatibility adapter or an explicitly documented legacy boundary, never final coverage truth.
- [ ] **Step 5: Run and commit.** Run `PYTHONPATH=src pytest tests/unit/test_claimable_parent_batches.py tests/unit/test_universal_claimable_block.py -v`; update `workstate.md` to state the shared typed path is mainline and north-star fused-scaffold work is optional; commit `refactor(namepredict): remove legacy claimable-side bypasses`.

### Task 9: Full validation and benchmark

**Files:**
- Modify: `workstate.md`

- [ ] **Step 1: Run focused suite.**

```bash
source .venv/Scripts/activate
export PYTHONPATH=src
pytest tests/unit/test_universal_claimable_block.py \
  tests/unit/test_parent_ownership.py \
  tests/unit/test_claimable_block_api.py \
  tests/unit/test_coverage_ledger.py \
  tests/unit/test_substituent_namer.py \
  tests/unit/test_rooted_alkyl_backend.py \
  tests/unit/test_amide_n_claims.py \
  tests/unit/test_open_chain_alkoxy_claims.py \
  tests/unit/test_candidate_retry.py \
  tests/unit/test_claimable_parent_batches.py -v
```

Expected: all pass.

- [ ] **Step 2: Run structural and benchmark gates.**

```bash
python tools/structure_lint.py --root src/namepredict
python -m benchmarks.benchmark_parallel --data data/merged_benchmark.json --time
```

Expected: lint passes and dual fallback does not regress more than 0.5% from the recorded pre-change baseline.

- [ ] **Step 3: Record result and commit.** Record exact benchmark before/after values and any intentionally unsupported class in `workstate.md`; commit `test(namepredict): verify universal claimable-block coverage`.

---

## Optional appendix: recursive capability expansion after mainline

This appendix is strictly downstream of Tasks 1–9. It is not a dependency of coverage, candidate retry, N-cycloalkyl, substituted phenyl, rooted alkyl, or open-chain alkoxy. Do not keep mainline tasks red while attempting it.

1. Add class-level recursive retained support for cycloalkyl `-yl` and substituted phenyl `-yl`, with at least three examples per class.
2. Add a fused heterocycle family only after its free parent, attachment locants, and `yl_form` have independent tests; do not add a single north-star scaffold exception.
3. Add acyl heterocycle leaves as a separately typed retained/recursive backend capability, not an amide-specific bypass.
4. The prior recursive-substituent north-star tasks may be executed here only after being rewritten to consume `ClaimedBlock`, `SubstituentNamer`, `owned_atoms`, and `CoverageLedger`. They may not restore removed routes or become a prerequisite for Tasks 1–9.

## Explicitly out of scope

- Perfect complex substituent naming for all approximately 141 parent kinds.
- A universal fused-ring nomenclature engine.
- Polyether and amino side-chain trees.
- Benchmark gold changes.
- Single-SMILES special-case scaffolds.

## Dependency order

```text
Task 1 red bars/inventory
  -> Task 2 terminal ownership/candidates
  -> Task 3 ownership-only claims + ledger
  -> Task 4 typed ordered namer
  -> Task 5 rooted-carbon backend
  -> Task 6 N/alkoxy migration
  -> Task 7 coverage-gated retry
  -> Task 8 bypass deletion + parent batches
  -> Task 9 validation
  -> Optional appendix
```

## Self-review

- Confirmed design is represented directly: claim has ownership/boundary/attachment only; no claim naming mode exists.
- `CoverageLedger` checks both gap and overlap before any assembled success, and failed candidates retry in rank order.
- `owned_atoms` is terminal source of truth; `parent_atom_set` is expressly a legacy adapter.
- The backend order and typed `SubstituentNamer`/`SubstituentName` contracts are fixed and used by migration tasks.
- The old recursive north-star work is optional and explicitly downstream, resolving the former mainline/appendix contradiction.
- Class-level goals, TDD, commits, lint, and benchmark gates are retained.

## Execution Handoff

Plan complete and saved to `docs/superpowers/plans/2026-07-18-universal-claimable-block.md`. Two execution options:

1. **Subagent-Driven (recommended)** — dispatch a fresh subagent per task, with review between tasks.
2. **Inline Execution** — execute tasks in this session using `executing-plans`, with checkpoints.

Which approach?
