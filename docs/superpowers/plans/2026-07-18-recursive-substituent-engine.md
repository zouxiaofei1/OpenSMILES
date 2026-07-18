# Recursive Substituent Engine Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 在**母体已选定**的前提下，把贴在母体边界外的复杂块（尤其 benzamide 的 N-substituent 与环上取代）统一为 **`claimable_block` → 简单叶子快路径 / 剪开后全管线 `name_as_substituent` → IUPAC 收成 -yl**，使 `N-(heteroaryl)` 与更深层嵌套可组合；北极星分子为 4-methoxy-N-(4-(pyrrolidine-1-carbonyl)-5,6-dihydro-4H-cyclopenta[d]thiazol-2-yl)benzamide。

**Architecture:** 不把「任意子图递归」塞进 `AlkylTree`。总策略固定为：`select_parent` 锁定母体原子集 → 枚举外沿块 → 简单块走现有 probe /（并行 plan 的）碳树 / alkoxy / leaves → 否则 `cut_block` 得子图 + 连接点 → `name_as_substituent(submol, attach)` 再跑 L0–L5（带 depth）→ `yl_form` 按 IUPAC 收取代基形 → L3/L5 组装。裁掉母体后重原子严格减少，递归终止有保证。本 plan **不**用 `max_depth` 否定架构；depth 只作实现保险丝（默认 ≤4）。

**Tech Stack:** Python 3.12、RDKit、pytest、`SMILESNNamer` dual、`tools/structure_lint.py`、触及命名时 `benchmarks.benchmark_parallel`。

## Global Constraints

- 代码只落在 `src/namepredict/layer0`–`layer5`、`namer.py` 与 `tests/unit/`；`data/*` 只读。
- 单文件 ≤500 行；函数体非空非注释行 ≤10；逼近 500 行必须拆文件、不跨层。
- 层边界：
  - **L2**：母体原子集、块拓扑 cut、claim 成败；**不**拼最终多语种串。
  - **L3**：`name_as_substituent` / `yl_form` / N-block extract → `{kind, attach_idx, atoms, en, zh, paren?}`。
  - **L4/L5**：位次与组装；禁止为个案改字符串补丁。
- 禁止深度学习；禁止无界 cache；禁止改 benchmark 金标数据。
- TDD：先红测 → 最小实现 → 绿 → commit；每任务只做一个可验收切片。
- 触及命名：相关 unit 全绿；必要时 `python -m benchmarks.benchmark_parallel --data data/merged_benchmark.json --time`，dual 回退不得 >0.5%。
- IUPAC：`docs/iupac/` **P-29**（取代基 / -yl）、**P-66.1**（酰胺 / benzamide）、**P-22**（杂环母体氢化物）、**P-25**（稠环，融环切片）、**P-65**（酰基前缀，如 pyrrolidine-1-carbonyl）。
- **禁止** 在 `benzamide.py` 内复制整棵侧链/递归命名器；只放宽 gate + 调共享 API。
- **禁止** 「不认识侧块 → 当空气」；能 claim 则命名，否则 gate/extract 失败（显式），不得静默蒸发碳/杂原子。
- **禁止** 把开链饱和烷基塞进本 plan 的递归主路径（那是 `2026-07-18-side-chain-recursive-engine` 的叶子）；本 plan 消费「简单叶子已成功」或调用共享 `claimable_side`（若已存在）。
- Windows Git Bash：`source .venv/Scripts/activate` 后再 `pytest` / `python`。
- 与并行 side-chain plan 关系：本 plan **不**依赖其完成即可启动 P0–P2；若 `claimable_side` / systematic alkyl 已合并，P1 快路径优先调用，否则 retained probe 即可。

---

## 北极星与分阶段能力

**目标 dual（用户给定）：**

| 语 | 名 |
|----|-----|
| EN | `4-methoxy-N-(4-(pyrrolidine-1-carbonyl)-5,6-dihydro-4H-cyclopenta[d]thiazol-2-yl)benzamide` |
| ZH | `4-甲氧基-N-(4-(吡咯烷-1-羰基)-5,6-二氢-4H-环戊并[d]噻唑-2-基)苯甲酰胺` |

**参考 SMILES（实现时以 RDKit 规范化后断言；两式等价则均可）：**

```text
COc1ccc(C(=O)Nc2nc3c(s2)CCC3C(=O)N2CCCC2)cc1
```

**现状（实施前基线，2026-07-18 探测）：**

| SMILES / 结构 | 现状 |
|---------------|------|
| `c1ccccc1C(=O)N` | benzamide ✅ |
| `O=C(N)c1ccc(OC)cc1` | 4-methoxybenzamide ✅ |
| `c1ccccc1C(=O)Nc2ccccc2` | N-phenylbenzamide ✅ |
| `c1ccccc1C(=O)Nc2ncccc2` | phenylformamide ❌（应 N-(pyridin-2-yl)benzamide） |
| `c1ccccc1C(=O)Nc2nccs2` | phenylformamide ❌（应 N-(1,3-thiazol-2-yl)benzamide） |
| `CC(=O)N1CCCC1` | N,N-dimethylacetamide ❌（酰基/吡咯烷路径另修） |
| `c1nc2c(s1)CCC2` | methane ❌（融环未支持） |
| 北极星 SMILES | ethane ❌ |

**能力分层（本 plan 交付）：**

| ID | 能力 | 本 plan |
|----|------|---------|
| R0 | `parent_atom_set` / `cut_block` / 原子数递减不变量 | **必做** |
| R1 | `name_as_substituent` + depth 保险丝 + `yl_form`（表驱动 retained） | **必做** |
| R2 | benzamide / amide **gate**：N 上允许 claimable 复杂块（非仅 C1–4/Ph） | **必做** |
| R3 | L3 `n_block` extract：`N-(…yl)` 组装（复用 n_benzyl 括号风格） | **必做** |
| R4 | 第一批叶子：未取代 **pyridin-2-yl / 1,3-thiazol-2-yl** 挂 N | **必做** |
| R5 | 子图上简单环取代（如 4-methyl-1,3-thiazol-2-yl） | **必做** |
| R6 | **pyrrolidine-1-carbonyl** 酰基前缀（-C(=O)-N 环） | **必做** |
| R7 | **5,6-dihydro-4H-cyclopenta[d]thiazole** retained / scaffold 母体 | **必做** |
| R8 | 组合北极星：benzamide + 4-methoxy + N-(4-(pyrrolidine-1-carbonyl)-…-2-yl) | **必做** |
| R9 | 开链 amide 同步 N-block（非仅 benzamide） | **Phase 后段可选** |
| R10 | 通用任意融环 Hantzsch–Widman 自动导出 | **不在本 plan**（只做北极星所需一片） |
| R11 | 侧链 N 脂肪树 / polyether | **不在本 plan**（side-chain / 另开） |

---

## 目标态文件结构

```
src/namepredict/
  layer2/
    block_cut.py              # 新建：parent 原子集、外沿根、cut 连通块
    submol_build.py           # 新建：子图 Mol + atom_map + attach_new
    benzamide.py              # 改：gate 接受 n_block claimable
    alkenamide.py             # 改：_amide_n_meta 增加复杂 N 路径钩子
    # 融环切片：
    cyclopenta_thiazole.py    # 新建：5,6-dihydro-4H-cyclopenta[d]thiazole parent
  layer3/
    yl_form.py                # 新建：free PIN → …-yl @ attach
    as_substituent.py         # 新建：name_as_substituent 管道
    n_block_extract.py        # 新建：amide/benzamide N 复杂块 → n_block sub
    substituent_extractor.py  # 改：挂 n_block extract
  layer5/
    assembler_prefixes.py     # 改：n_block 与 n_benzyl 同样 omit 位次逻辑
  namer.py                    # 改：内部 _name_mol / depth 参数（保持 public API 稳定）
tests/unit/
  test_block_cut.py
  test_yl_form.py
  test_as_substituent.py
  test_n_heteroaryl_benzamide.py
  test_pyrrolidine_1_carbonyl.py
  test_cyclopenta_thiazole.py
  test_recursive_benzamide_north_star.py
```

---

## 核心接口（全任务共用，禁止改名）

```python
# layer2/block_cut.py
def parent_atom_set(parent: dict, mol) -> frozenset[int]:
    """Atoms owned by selected parent (ring/chain + FG atoms declared by kind)."""

def side_roots(mol, parent_atoms: frozenset[int]) -> list[int]:
    """Neighbor atoms of parent_atoms that are outside the set (block roots)."""

def cut_block(mol, root: int, parent_atoms: frozenset[int]) -> frozenset[int] | None:
    """Connected component from root not entering parent_atoms; None if empty/invalid."""

# layer2/submol_build.py
@dataclass(frozen=True)
class CutSubmol:
    mol: object                 # RDKit Mol, attach capped with H for free naming
    atom_map: dict[int, int]    # new_idx -> old_idx
    inv_map: dict[int, int]     # old_idx -> new_idx
    attach_new: int             # attach atom index in submol
    attach_old: int
    atoms_old: frozenset[int]

def build_cut_submol(mol, atoms: frozenset[int], attach_old: int) -> CutSubmol | None:
    """Induced submol on atoms; free valence at attach filled with H."""

# layer3/yl_form.py
def yl_form(en: str, zh: str, attach_locant: int | None, *, paren: bool = True) -> tuple[str, str, bool]:
    """Free parent name → substituent …-yl / …-基. IUPAC P-29."""

# layer3/as_substituent.py
def name_as_substituent(
    mol, attach_old: int, atoms: frozenset[int], *, depth: int = 0, max_depth: int = 4,
) -> tuple[str, str, bool] | None:
    """Cut → pipeline on submol → yl_form. None if unsupported. depth increases each recurse."""

# layer3/n_block_extract.py
def extract_n_blocks(info: dict, parent: dict) -> list[dict]:
    """For amide/benzamide: one complex N-substituent → kind=n_block, en like N-(…)."""
```

**递归不变量（写入测试）：**

1. `len(atoms_old) < mol.GetNumAtoms()`（相对全分子）或相对父调用 `atoms` 严格更小。  
2. `depth < max_depth`，否则返回 `None`（失败，不截断装成功）。  
3. 子调用 **禁用** `CommonNameCache` 写回全分子 cache 键污染（子 SMILES 可只读 cache，成功结果不 `put` 到外层 smiles；或使用 `cache=None`）。

**组装约定（N-block）：**

- 与 `n_benzyl` 一致：复杂则 `N-(en)` / `N-(zh)`；简单无位次可 `N-en`（本 plan 复杂块一律括号）。  
- `attach_idx`：benzamide 用 `ring_attach_idx`（与现有 n_alkyl 一致），便于 L4 省略 N 位次。

---

## Task 1: 红测 — N-heteroaryl benzamide + 北极星锁名

**Files:**

- Create: `tests/unit/test_n_heteroaryl_benzamide.py`
- Create: `tests/unit/test_recursive_benzamide_north_star.py`

**Interfaces:**

- Consumes: `SMILESNNamer`, `normalize_en`, `normalize_zh`
- Produces: 锁定 dual 名，供后续 task 变绿

- [ ] **Step 1: 写失败测试（MVP 路径）**

```python
# IUPAC: P-66.1.1 / P-29 / P-22.2.1
# Layer: L2,L3,L5
"""N-heteroaryl benzamide via recursive substituent path."""
from __future__ import annotations

import pytest
from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

MVP = [
    (
        "c1ccccc1C(=O)Nc2ncccc2",
        "N-(pyridin-2-yl)benzamide",
        "N-(吡啶-2-基)苯甲酰胺",
    ),
    (
        "c1ccccc1C(=O)Nc2nccs2",
        "N-(1,3-thiazol-2-yl)benzamide",
        "N-(1,3-噻唑-2-基)苯甲酰胺",
    ),
    (
        "O=C(Nc1nccs1)c1ccc(OC)cc1",
        "4-methoxy-N-(1,3-thiazol-2-yl)benzamide",
        "4-甲氧基-N-(1,3-噻唑-2-基)苯甲酰胺",
    ),
]

REGRESS = [
    ("c1ccccc1C(=O)N", "benzamide", "苯甲酰胺"),
    ("O=C(N)c1ccc(OC)cc1", "4-methoxybenzamide", "4-甲氧基苯甲酰胺"),
    ("c1ccccc1C(=O)Nc2ccccc2", "N-phenylbenzamide", "N-苯基苯甲酰胺"),
    ("c1ccccc1C(=O)NC", "N-methylbenzamide", "N-甲基苯甲酰胺"),
    ("CC(=O)N", "acetamide", "乙酰胺"),
]


@pytest.mark.parametrize("smiles,en,zh", MVP)
def test_n_heteroaryl_benzamide_mvp(smiles, en, zh):
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)


@pytest.mark.parametrize("smiles,en,zh", REGRESS)
def test_n_heteroaryl_benzamide_regress(smiles, en, zh):
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)
```

- [ ] **Step 2: 北极星红测（允许长期红到 Task 8）**

```python
# tests/unit/test_recursive_benzamide_north_star.py
# IUPAC: P-66.1 / P-29 / P-25 / P-65
# Layer: L2–L5
from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

NORTH = (
    "COc1ccc(C(=O)Nc2nc3c(s2)CCC3C(=O)N2CCCC2)cc1",
    "4-methoxy-N-(4-(pyrrolidine-1-carbonyl)-5,6-dihydro-4H-cyclopenta[d]thiazol-2-yl)benzamide",
    "4-甲氧基-N-(4-(吡咯烷-1-羰基)-5,6-二氢-4H-环戊并[d]噻唑-2-基)苯甲酰胺",
)


def test_north_star_benzamide():
    smiles, en, zh = NORTH
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)
```

- [ ] **Step 3: 跑测确认红**

```bash
source .venv/Scripts/activate
pytest tests/unit/test_n_heteroaryl_benzamide.py tests/unit/test_recursive_benzamide_north_star.py -v
```

Expected: MVP 与 NORTH FAIL；REGRESS PASS。

- [ ] **Step 4: Commit 红测**

```bash
git add tests/unit/test_n_heteroaryl_benzamide.py tests/unit/test_recursive_benzamide_north_star.py
git commit -m "$(cat <<'EOF'
test(namepredict): red-bar recursive N-heteroaryl benzamide

Lock MVP N-(pyridin/thiazol-2-yl) and north-star dual names.

Co-Authored-By: Claude <noreply@anthropic.com>
EOF
)"
```

---

## Task 2: L2 — `parent_atom_set` / `cut_block` / `build_cut_submol`

**Files:**

- Create: `src/namepredict/layer2/block_cut.py`
- Create: `src/namepredict/layer2/submol_build.py`
- Create: `tests/unit/test_block_cut.py`

**Interfaces:**

- Produces: `parent_atom_set`, `side_roots`, `cut_block`, `CutSubmol`, `build_cut_submol`

- [ ] **Step 1: 拓扑单测**

```python
# tests/unit/test_block_cut.py
from namepredict.layer0.preprocessor import preprocess
from namepredict.layer1.analyzer import analyze
from namepredict.layer2.parent_selector import select_parent
from namepredict.layer2.block_cut import parent_atom_set, cut_block, side_roots
from namepredict.layer2.submol_build import build_cut_submol


def test_cut_n_phenyl_from_benzamide():
    mol = preprocess("c1ccccc1C(=O)Nc2ccccc2")
    parent = select_parent(analyze(mol))
    assert parent["kind"] == "benzamide"
    patoms = parent_atom_set(parent, mol)
    # amide N is in parent; N-phenyl ring is outside
    am = analyze(mol)["amides"][0]
    n_idx = am["n_idx"]
    roots = [r for r in side_roots(mol, patoms) if r in set(am.get("n_c_idxs") or [])]
    assert len(roots) == 1
    block = cut_block(mol, roots[0], patoms)
    assert block is not None and len(block) == 6
    sub = build_cut_submol(mol, block, roots[0])
    assert sub is not None
    assert sub.mol.GetNumAtoms() >= 6
    assert sub.attach_new in sub.atom_map
```

- [ ] **Step 2: 跑测红**

```bash
pytest tests/unit/test_block_cut.py -v
```

Expected: import/FAIL。

- [ ] **Step 3: 实现 `block_cut.py`（每函数 ≤10 行）**

要点：

```python
def parent_atom_set(parent: dict, mol) -> frozenset[int]:
    """Union chain + kind-specific FG atoms (amide C/N/O for benzamide/amide)."""
    # benzamide: chain(ring) + amide_c_idx + n_idx + dbl O
    # do NOT include N-substituent carbons

def side_roots(mol, parent_atoms: frozenset[int]) -> list[int]:
    # heavy neighbors of any parent atom that are not in parent_atoms

def cut_block(mol, root: int, parent_atoms: frozenset[int]) -> frozenset[int] | None:
    # BFS from root; never enter parent_atoms; include all heavies in component
```

`submol_build.py`：

```python
def build_cut_submol(mol, atoms: frozenset[int], attach_old: int) -> CutSubmol | None:
    # EditableMol: copy atoms in `atoms`, bonds within set;
    # at attach_old, add H to satisfy valence for free naming;
    # build atom_map / inv_map; Sanitize
```

- [ ] **Step 4: 绿测 + lint**

```bash
pytest tests/unit/test_block_cut.py -v
python tools/structure_lint.py --root src/namepredict/layer2/block_cut.py
python tools/structure_lint.py --root src/namepredict/layer2/submol_build.py
```

- [ ] **Step 5: Commit**

```bash
git add src/namepredict/layer2/block_cut.py src/namepredict/layer2/submol_build.py tests/unit/test_block_cut.py
git commit -m "$(cat <<'EOF'
feat(namepredict): cut parent-boundary blocks into submols

parent_atom_set + cut_block + build_cut_submol for recursive sides.

Co-Authored-By: Claude <noreply@anthropic.com>
EOF
)"
```

---

## Task 3: L3 — `yl_form` + `name_as_substituent`（无宿主接线）

**Files:**

- Create: `src/namepredict/layer3/yl_form.py`
- Create: `src/namepredict/layer3/as_substituent.py`
- Modify: `src/namepredict/namer.py`（抽出 `_name_mol(mol, *, depth, use_cache)`，public `name(smiles)` 不变）
- Create: `tests/unit/test_yl_form.py`
- Create: `tests/unit/test_as_substituent.py`

**Interfaces:**

- Consumes: `build_cut_submol`, `CutSubmol`, pipeline layers
- Produces: `yl_form`, `name_as_substituent`

- [ ] **Step 1: `yl_form` 单测**

```python
# tests/unit/test_yl_form.py
from namepredict.layer3.yl_form import yl_form

def test_pyridine_to_yl():
    en, zh, paren = yl_form("pyridine", "吡啶", 2)
    assert en == "pyridin-2-yl"
    assert zh == "吡啶-2-基"
    assert paren is True

def test_thiazole_to_yl():
    en, zh, paren = yl_form("1,3-thiazole", "1,3-噻唑", 2)
    assert en == "1,3-thiazol-2-yl"
    assert "噻唑-2-基" in zh

def test_benzene_to_yl():
    en, zh, paren = yl_form("benzene", "苯", 1)
    # unsubstituted phenyl convention
    assert en in ("phenyl", "benzen-1-yl")  # pick ONE in implement; lock phenyl
```

**实现规则（本 plan 固定，写入 `yl_form` docstring + 测试）：**

| free EN | attach locant | yl EN |
|---------|---------------|-------|
| `benzene` | 1 | `phenyl`（retained 捷径） |
| `pyridine` | k | `pyridin-k-yl` |
| `1,3-thiazole` | k | `1,3-thiazol-k-yl` |
| `…ole` / `…ine` 等 | k | 去末端 e（若有）+ `-k-yl` |
| 已有前缀的 free 名 | k | 前缀保留 + 母体改 yl：`4-methyl-1,3-thiazol-2-yl` |

中文：`吡啶` → `吡啶-2-基`；`1,3-噻唑` → `1,3-噻唑-2-基`；`苯` → `苯基`。

- [ ] **Step 2: 实现 `yl_form`（表 + 小规则，禁止巨型 if 链散落）**

- [ ] **Step 3: `namer` 抽出 mol 入口**

```python
# namer.py 概念
def _name_mol(mol, *, depth: int = 0, cache: CommonNameCache | None = None) -> NameResult:
    # same layers as _pipeline but starts from mol; if depth>0: cache write disabled

def name_as_substituent(...):
    sub = build_cut_submol(...)
    if sub is None or depth >= max_depth:
        return None
    # optional: force numbering so attach atom gets known locant — MVP:
    # name free sub.mol; compute attach_locant from parent numbering in sub result meta
    # if meta lacks locant map, use ring-orient helpers for known kinds (pyridine N=1 etc.)
```

**MVP 位次策略（Task 3–4）：**

对 **未取代** pyridine / 1,3-thiazole：

- 子图 `select_parent` 已有固定杂原子位次（吡啶 N=1；噻唑 S=1,N=3）。  
- `attach_locant` = 子图 parent 编号下 attach 原子的位次（从 `number()` 结果或 ring orient 读）。  
- 若读不到位次 → `None`（失败），不要猜。

- [ ] **Step 4: `as_substituent` 单测**

```python
def test_name_thiazol_2_yl_from_benzamide_context():
    mol = preprocess("c1ccccc1C(=O)Nc2nccs2")
    # locate thiazole atoms + attach carbon bonded to amide N
    ...
    got = name_as_substituent(mol, attach_old, atoms, depth=0)
    assert got is not None
    en, zh, paren = got
    assert en == "1,3-thiazol-2-yl"
```

- [ ] **Step 5: 绿测 + commit**

```bash
pytest tests/unit/test_yl_form.py tests/unit/test_as_substituent.py -v
git add src/namepredict/layer3/yl_form.py src/namepredict/layer3/as_substituent.py \
  src/namepredict/namer.py tests/unit/test_yl_form.py tests/unit/test_as_substituent.py
git commit -m "$(cat <<'EOF'
feat(namepredict): name_as_substituent + yl_form core

Cut submol, rerun layers, emit P-29 -yl dual names.

Co-Authored-By: Claude <noreply@anthropic.com>
EOF
)"
```

---

## Task 4: benzamide gate 放宽 + `_amide_n_meta` 钩子

**Files:**

- Modify: `src/namepredict/layer2/benzamide.py`
- Modify: `src/namepredict/layer2/alkenamide.py`（或新建 `layer2/amide_n_block.py` 若行数爆）
- Modify: `src/namepredict/layer2/arene_carbonyl.py` 仅当 `_arene_subs_ok` 排除集需纳入 N-block 原子
- Test: `tests/unit/test_benzamide.py` 回归 + 新单测 gate

**Interfaces:**

- Produces: parent 上 `n_block: True`（或 `n_block_root: int`）当 N 上存在复杂可 claim 块
- 保持 simple N 路径不变（methyl / phenyl 仍走旧 meta）

- [ ] **Step 1: 行为规格**

```text
_amide_n_meta 顺序：
  1. 现有 n_benzyl / n_phenyl / n_alkyl（简单）
  2. else if 单一 n_c 且 cut_block 可 claim → {"n_block": True, "n_block_root": c}
  3. else {}  → benzamide gate 仍失败（与现网一致）

_is_simple_benzamide:
  - 简单 N：旧逻辑
  - n_block：要求 name_as_substituent 在 L2 可预检成功？ 
    **本 plan 选择：** L2 只做拓扑 claim（cut_block 非空 + 非简单已处理）；
    真正命名失败留给 L3（则 assemble 失败或回退 — MVP 要求 L3 成功，红测锁死）。
  - 环上取代：仍走 _arene_subs_ok；exclude 必须并入 N-block 全部 atoms，避免「环外未知」拒收
```

- [ ] **Step 2: 实现 exclude 打包**

```python
def _amide_pack_benz(mol, am):
    excl, allowed = ...  # existing
    # if complex N later flagged, caller unions cut_block atoms into excl
```

在 `_is_simple_benzamide` / `_benzamide_parent`：

```python
meta = _amide_n_meta(info)
# if meta.get("n_block"):
#   atoms = cut_block(mol, meta["n_block_root"], parent_core_atoms)
#   excl |= atoms
```

- [ ] **Step 3: 单测 parent kind**

```python
def test_benzamide_parent_accepts_n_thiazole():
    mol = preprocess("c1ccccc1C(=O)Nc2nccs2")
    parent = select_parent(analyze(mol))
    assert parent["kind"] == "benzamide"
    assert parent.get("n_block") is True
```

- [ ] **Step 4: 跑 `test_benzamide.py` 全绿 + commit**

```bash
pytest tests/unit/test_benzamide.py tests/unit/test_n_alkyl_amide.py tests/unit/test_n_phenyl_amide.py -v
git commit -m "$(cat <<'EOF'
feat(namepredict): benzamide gate allows claimable N-block

Topology claim for complex N; exclude block atoms from ring side check.

Co-Authored-By: Claude <noreply@anthropic.com>
EOF
)"
```

---

## Task 5: L3 `extract_n_blocks` + L5 组装

**Files:**

- Create: `src/namepredict/layer3/n_block_extract.py`
- Modify: `src/namepredict/layer3/substituent_extractor.py`（`extract_substituents` 串联）
- Modify: `src/namepredict/layer5/assembler_prefixes.py`（`n_block` 视同 N-前缀 omit 位次）
- Test: MVP 用例应变绿

**Interfaces:**

- Consumes: `name_as_substituent`, parent `n_block` / `n_block_root`
- Produces: `{"kind": "n_block", "attach_idx", "atoms", "en": "N-(…)", "zh": "N-(…)", "paren": True}`

- [ ] **Step 1: extract 实现**

```python
def extract_n_blocks(info, parent):
    if parent.get("kind") not in {"amide", "benzamide"} or not parent.get("n_block"):
        return []
    mol, root = info["mol"], parent["n_block_root"]
    patoms = parent_atom_set(parent, mol)  # without N-block
    atoms = cut_block(mol, root, patoms)
    if atoms is None:
        return []
    named = name_as_substituent(mol, root, atoms, depth=0)
    if named is None:
        return []
    en, zh, paren = named
    attach = parent.get("ring_attach_idx") if parent["kind"] == "benzamide" else parent.get("amide_c_idx")
    return [{
        "kind": "n_block", "n_carbons": 0, "attach_idx": attach,
        "atoms": list(atoms), "en": f"N-({en})", "zh": f"N-({zh})", "paren": True,
    }]
```

注意：`parent_atom_set` 在 extract 时**不得**已含 N-block 原子（与 gate 打包一致）。

- [ ] **Step 2: `substituent_extractor` 接入顺序**

在现有 `_extract_n_alkyl + _extract_n_phenyl + _extract_n_benzyl` 旁增加 `extract_n_blocks`；**互斥**：有 n_block 则不应再出 n_alkyl。

- [ ] **Step 3: L5**

```python
# assembler_prefixes.py
if kind in ("sec_amine", "tert_amine", "amide", "benzamide"):
    return {s.get("kind") for s in substituents} <= {
        "n_alkyl", "n_phenyl", "n_benzyl", "n_block",
    }
# _parts_for_stem: n_block → omit = True（同 n_benzyl）
```

确认 `benzene_names` / special stem 对 benzamide 前缀拼接不丢 `N-(…)`。

- [ ] **Step 4: 跑 MVP 绿**

```bash
pytest tests/unit/test_n_heteroaryl_benzamide.py tests/unit/test_benzamide.py -v
```

Expected: MVP 全绿；REGRESS 全绿；NORTH 仍红。

- [ ] **Step 5: Commit**

```bash
git commit -m "$(cat <<'EOF'
feat(namepredict): N-block extract and assemble for benzamide

Wire name_as_substituent into N-(heteroaryl)benzamide dual names.

Co-Authored-By: Claude <noreply@anthropic.com>
EOF
)"
```

---

## Task 6: 子图简单环取代（4-methyl-1,3-thiazol-2-yl）

**Files:**

- Modify: `as_substituent.py` / 子 pipeline 已有 mono 取代能力则只加测
- Modify: `azole13` / `pyridine` 仅当子图 gate 过严导致带甲基噻唑失败时放宽 **子图场景**（禁止为宿主特例改字符串）
- Test: 扩 `test_n_heteroaryl_benzamide.py`

- [ ] **Step 1: 红测**

```python
(
    "c1ccccc1C(=O)Nc2nc(C)cs2",
    "N-(4-methyl-1,3-thiazol-2-yl)benzamide",
    "N-(4-甲基-1,3-噻唑-2-基)苯甲酰胺",
),
```

（位次以 IUPAC 噻唑编号为准；实现前用 `SMILESNNamer` 对游离 `Cc1ncsc1` / `Cc1cscn1` 核对 free 名再锁。）

- [ ] **Step 2: 保证子调用 `extract_substituents` 在 thiazole parent 上能提出 methyl，再 `yl_form` 保留前缀**

`yl_form` 必须支持 **带前缀的 free 名**：

```python
yl_form("4-methyl-1,3-thiazole", "4-甲基-1,3-噻唑", 2)
→ ("4-methyl-1,3-thiazol-2-yl", "4-甲基-1,3-噻唑-2-基", True)
```

- [ ] **Step 3: 绿测 + commit**

```bash
git commit -m "$(cat <<'EOF'
feat(namepredict): substituted heteroaryl-yl via sub-pipeline

Keep ring prefixes when forming N-(…thiazol-2-yl).

Co-Authored-By: Claude <noreply@anthropic.com>
EOF
)"
```

---

## Task 7: `pyrrolidine-1-carbonyl` 酰基前缀（R6）

**Files:**

- Create: `src/namepredict/layer3/acyl_radical_names.py`（或 `layer2/side_acyl.py` 拓扑 + L3 名）
- Modify: 子图 `extract_substituents` / ring side 识别 `-C(=O)-N` 且 N 在饱和氮杂环
- Test: `tests/unit/test_pyrrolidine_1_carbonyl.py`

**化学定义（写死）：**

```text
连接点 = 羰基碳（挂在宿主环上）
结构 = 环碳/芳碳 — C(=O) — N（吡咯烷 N，另两键为环内 CH2）
EN: pyrrolidine-1-carbonyl
ZH: 吡咯烷-1-羰基
```

**不做：** 任意内酰胺、开链 -C(=O)NMe2 全称系统化（可另 task）；本 task 只认 **pyrrolidine**（五元饱和 NH 环，N 连羰基）。

- [ ] **Step 1: 游离/侧基红测**

```python
CASES = [
    # 最小：噻唑 4-位挂 pyrrolidine-1-carbonyl，2-位 H
    # 用可解析 SMILES；期望 free 或取代基路径稳定
    (
        "O=C(c1cscn1)N1CCCC1",
        "4-(pyrrolidine-1-carbonyl)-1,3-thiazole",  # 位次以 orient 为准，实现时锁定实测 free 名
        "4-(吡咯烷-1-羰基)-1,3-噻唑",
    ),
]
```

若 free 母体选择更倾向 amide 而非 thiazole，**取代基模式**下应优先 **保留连接点所在杂环为 parent**（Task 3 的 attach 约束增强）：

```python
# as_substituent: when orienting, prefer parents whose chain/ring contains attach_new
```

本 task 最小实现可二选一（文档锁定一种）：

**A（推荐）：** L2/L3 专用 `match_pyrrolidine_1_carbonyl(mol, c_carbonyl)` → 直接返回 en/zh，作环侧链 fact，**不**走错误的 acetamide 母体。  
**B：** 修 `CC(=O)N1CCCC1` 全局母体（更大，易回归）。

**本 plan 指定 A**，避免被错误 amide 母体绑架；另开 bugfix 可修 B。

- [ ] **Step 2: 实现 match + 接入杂环 / 子 pipeline 侧链 extract**

```python
def match_pyrrolidine_1_carbonyl(mol, c_idx: int) -> frozenset[int] | None:
    # c is C=O carbon; double O; single N; N in 5-sat ring all C else; return atoms set
```

- [ ] **Step 3: 与 `name_as_substituent` 组合测**

`N-(4-(pyrrolidine-1-carbonyl)-1,3-thiazol-2-yl)benzamide` 中间里程碑（无稠环）。

- [ ] **Step 4: Commit**

```bash
git commit -m "$(cat <<'EOF'
feat(namepredict): pyrrolidine-1-carbonyl ring prefix leaf

Direct match avoids false acetamide parent on N-acylpyrrolidine.

Co-Authored-By: Claude <noreply@anthropic.com>
EOF
)"
```

---

## Task 8: `5,6-dihydro-4H-cyclopenta[d]thiazole` 母体（R7）

**Files:**

- Create: `src/namepredict/layer2/cyclopenta_thiazole.py`
- Modify: `layer2/fg_producers.py` / `ring_producers.py` / `kind_registry.py` / L4 orient / L5 stems（按现网 fused56 注册模式）
- Test: `tests/unit/test_cyclopenta_thiazole.py`

**范围（故意窄）：**

```text
仅：噻唑与饱和环戊烷以 [d] 边稠合的 5,6-dihydro-4H-cyclopenta[d]thiazole
允许：2-位连接（宿主 N）、4-位 pyrrolidine-1-carbonyl、无其它杂取代
不做：任意 cyclopenta[d] 杂环自动命名引擎
```

- [ ] **Step 1: 游离母体红测**

```python
CASES = [
    ("c1nc2c(s1)CCC2", "5,6-dihydro-4H-cyclopenta[d]thiazole", "5,6-二氢-4H-环戊并[d]噻唑"),
    # 2-chloro or 2-methyl if needed to lock numbering
]
```

- [ ] **Step 2: 实现 parent 识别**

参考 `fused56` / `benzothiazole`：数环、芳香噻唑、饱和五元碳环共享边、固定位次（S/N 与 4/5/6 氢化标号按 IUPAC 稠环惯例写进测试与 orient）。

- [ ] **Step 3: 2-yl 形态**

```python
yl_form("5,6-dihydro-4H-cyclopenta[d]thiazole", "…", 2)
→ "5,6-dihydro-4H-cyclopenta[d]thiazol-2-yl"
```

- [ ] **Step 4: 带 4-(pyrrolidine-1-carbonyl) 的子图 free/substituent 名绿**

- [ ] **Step 5: Commit**

```bash
git commit -m "$(cat <<'EOF'
feat(namepredict): 5,6-dihydro-4H-cyclopenta[d]thiazole parent

Narrow fused scaffold for north-star N-block recursion.

Co-Authored-By: Claude <noreply@anthropic.com>
EOF
)"
```

---

## Task 9: 组合北极星 + 碳/杂原子覆盖断言

**Files:**

- Modify: 仅接线缺口（gate exclude、深度 2 递归：benzamide → N-block 子图 → 环上 acyl 叶）
- Create: `tests/unit/test_n_block_coverage.py`
- Test: `test_recursive_benzamide_north_star.py` 变绿

- [ ] **Step 1: 端到端**

```bash
pytest tests/unit/test_recursive_benzamide_north_star.py tests/unit/test_n_heteroaryl_benzamide.py -v
```

Expected: NORTH 绿。

- [ ] **Step 2: 覆盖断言**

```python
def heavy_coverage_gap(mol, parent, subst) -> set[int]:
    heavies = {a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1}
    covered = set(parent_atom_set(parent, mol))
    for s in subst:
        covered |= set(s.get("atoms") or [])
    return heavies - covered

def test_north_star_no_gap():
    ...
    assert heavy_coverage_gap(mol, parent, subst) == set()
```

- [ ] **Step 3: 递归深度测**

```python
def test_depth_guard_returns_none(monkeypatch):
    # force max_depth=0 on name_as_substituent → None, host must not emit empty N-
```

- [ ] **Step 4: Commit**

```bash
git commit -m "$(cat <<'EOF'
feat(namepredict): north-star recursive N-block benzamide

4-methoxy-N-(4-(pyrrolidine-1-carbonyl)-5,6-dihydro-4H-cyclopenta[d]thiazol-2-yl)benzamide.

Co-Authored-By: Claude <noreply@anthropic.com>
EOF
)"
```

---

## Task 10: 文档边界 + 回归锁 + workstate

**Files:**

- Modify: `workstate.md`（按仓库日志格式追加）
- Optional: `docs/superpowers/plans/2026-07-18-side-chain-recursive-engine.md` 顶部加 **See also** 链到本 plan（互指：叶子 vs 递归块）
- Run: unit 子集 + structure_lint + benchmark

- [ ] **Step 1: 跑验收**

```bash
source .venv/Scripts/activate
pytest tests/unit/test_block_cut.py tests/unit/test_yl_form.py \
  tests/unit/test_as_substituent.py tests/unit/test_n_heteroaryl_benzamide.py \
  tests/unit/test_pyrrolidine_1_carbonyl.py tests/unit/test_cyclopenta_thiazole.py \
  tests/unit/test_recursive_benzamide_north_star.py tests/unit/test_benzamide.py \
  tests/unit/test_n_phenyl_amide.py tests/unit/test_n_alkyl_amide.py -q
python tools/structure_lint.py --root src/namepredict
python -m benchmarks.benchmark_parallel --data data/merged_benchmark.json --time
```

- [ ] **Step 2:** dual 回退 ≤0.5%；记录 fails 变化  
- [ ] **Step 3:** workstate 一条：架构「母体已定 → cut → name_as_substituent → yl」，IUPAC P-29/P-66.1，与 side-chain plan 分工  
- [ ] **Step 4: Commit 文档**

```bash
git commit -m "$(cat <<'EOF'
docs(workstate): recursive substituent engine north-star

Co-Authored-By: Claude <noreply@anthropic.com>
EOF
)"
```

---

## 各「不会自动拥有」项（汇总）

| 缺口 | 做法 | Task |
|------|------|------|
| 切边 / 子图 | `parent_atom_set` + `cut_block` + `build_cut_submol` | 2 |
| free → -yl | `yl_form` 读 P-29 + 表驱动 | 3 |
| 整管线再跑 | `name_as_substituent` + `_name_mol(depth)` | 3 |
| benzamide 复杂 N | gate + exclude block atoms | 4 |
| 组装 N-(…) | `n_block` extract + L5 omit | 5 |
| 取代杂环基 | 子 pipeline 前缀 + yl_form | 6 |
| 吡咯烷-1-羰基 | **专用 match 叶**（不做错误 acetamide） | 7 |
| 稠环母体 | **窄** cyclopenta[d]thiazole | 8 |
| 北极星组合 | 深度 2 递归 + 覆盖断言 | 9 |
| 通用融环引擎 / 脂肪 N 树 / polyether | **不做** | — |

---

## 风险与缓解

| 风险 | 缓解 |
|------|------|
| 子调用选错母体（acetamide 抢走酰吡咯烷） | 酰基用专用 leaf；attach 约束优先含连接点的环系 |
| cache 污染 | `depth>0` 禁用 cache put |
| 递归装成功但半截名 | depth 耗尽 / yl 失败 → `None` → 上层失败，不蒸发 |
| benzamide 与 formamide 回归 | REGRESS 锁；gate 仍要求 Ph–C(=O)–N 拓扑 |
| 文件超 500 行 | 新逻辑只进 cut/submol/yl/as_substituent/n_block/cyclopenta_* |
| dual 回归 | 每 task 跑 benzamide 相关；Task 10 全 benchmark |
| 与 side-chain plan 抢 `_side_atoms` | 本 plan 复杂块走 cut 递归；简单烷基仍走叶子/并行 plan |
| 位次 2 vs 4 在稠环 | 单测锁 orient；先绿 free 母体再绿 -yl |

---

## 建议实施顺序（依赖图）

```
Task1 红测
  → Task2 cut/submol
  → Task3 yl_form + name_as_substituent
  → Task4 benzamide gate
  → Task5 n_block extract/assemble     ===== MVP: N-(pyridin/thiazol-2-yl)benzamide
  → Task6 取代杂环基
  → Task7 pyrrolidine-1-carbonyl leaf
  → Task8 cyclopenta[d]thiazole parent
  → Task9 北极星组合
  → Task10 验收
```

**最小可交付（MVP）：** Task 1–5 + 10 的子集（N-heteroaryl benzamide + 4-methoxy 组合）。  
**完整 plan 声明范围：** MVP + Task 6–9（北极星 dual）。

---

## Self-Review

1. **Spec coverage:** 讨论结论（母体已定、cut 变小、整管线再跑、yl 读 IUPAC、depth 非架构否定、专用酰基叶、窄融环）均有 Task；与 side-chain 叶子分工写明。  
2. **Placeholder scan:** 无 TBD；可选 R9/R10/R11 标不做。  
3. **类型一致:** `CutSubmol` / `name_as_substituent` / `yl_form` / `extract_n_blocks` / `n_block` 贯穿。  
4. **层纯度:** cut L2、命名 L3、组装 L5；benzamide 不内嵌递归器。  
5. **YAGNI:** 不做通用融环引擎；酰吡咯烷专用 match；不修全部错误 amide 母体 unless 挡路。

---

## Execution Handoff

Plan complete and saved to `docs/superpowers/plans/2026-07-18-recursive-substituent-engine.md`.

**两种执行方式：**

1. **Subagent-Driven（推荐）** — 每 Task 新开 subagent，Task 间 review  
2. **Inline Execution** — 本会话按 `executing-plans` 连续做，设 checkpoint  

你更想用哪一种？若只要 MVP，也可以说「先做 Task 1–5」。
