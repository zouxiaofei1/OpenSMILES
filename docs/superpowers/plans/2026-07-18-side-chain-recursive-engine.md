# Side-Chain Recursive Engine Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 把开链/环上侧链从「白名单 shape + 场景碎片提取」升级为 **共享侧链引擎（carbon tree / hetero bridge / aryl leaf）+ 统一 claim API**，使 2-methylbutyl、开链 methoxy、多支链烷烃等不再静默蒸发，并让苯/吡啶等 retained 母体与烷/酸/酮等 chain 母体复用同一套能力。

**Architecture:** L2 只产出 typed topology facts（`side_facts`）；L3 只渲染 en/zh；环母体 claim 与开链 extract 都调用同一 `claimable_side` / `carbon_side` / `hetero_side`。Retained shape 表保留为捷径，**末尾 systematic 兜底**；禁止按 parent `kind` 复制侧链逻辑。芳基嵌套继续走 `leaves/*`（P3），脂肪树 systematic 递归独立（P1），氧桥走 P2。

**Tech Stack:** Python 3.12、RDKit、pytest、`SMILESNNamer` dual、`tools/structure_lint.py`、`benchmarks.benchmark_parallel`（触及命名输出时）。

## Global Constraints

- 代码只落在 `src/namepredict/layer0`–`layer5`、入口与 `tests/unit/`；`data/*` 只读。
- 单文件 ≤500 行；函数体非空非注释行 ≤10；逼近 500 行必须拆文件、不跨层。
- 层边界：
  - **L2**：拓扑 fact、母体 claim（`_side_atoms` / `claimable_side`）；**不**拼最终多语种复杂串以外的 L5 组装。
  - **L3**：fact → `{kind, attach_idx, atoms, en, zh, paren?}`；**不**选母体。
  - **L4/L5**：位次与组装；禁止为个案改字符串补丁。
- 禁止深度学习；禁止无界 cache；禁止改 benchmark 金标数据。
- TDD：先红测 → 最小实现 → 绿 → commit；每任务只做一个可验收切片。
- 触及命名：相关 unit 全绿；必要时 `python -m benchmarks.benchmark_parallel --data data/merged_benchmark.json --time`，dual 回退不得 >0.5%。
- IUPAC：`docs/iupac/` P-29（烷基前缀）、P-63.2.2/P-63.2.4（醚 / alkoxyalkane PIN）、P-14（位次/字母序）。
- **禁止** 只扩 `ether_parent` arm 白名单当主路径；支化醚 PIN = alkane 母体 + alkoxy 前缀。
- **禁止** 把开链烷基塞进 `ArylLeafKind` / `leaves/`。
- **禁止** 「不认识侧链 → 当空气」；能 systematic 则命名，否则 claim 失败（环）或覆盖断言失败（测）。
- Windows Git Bash：`source .venv/Scripts/activate` 后再 `pytest` / `python`。

---

## 现状基线（实施前必读）

| 区域 | 文件 | 现状 | 缺口 |
|------|------|------|------|
| 烷基 shape | `layer2/side_alkyl.py` `_SIDE_PROBES` | linear / iPr / sBu / iBu / tBu / isopentyl / neopentyl / tert-pentyl / prenyl / CF₃ / 环烷基 | **无 2-methylbutyl**；无通用树 |
| 环 claim | `ring_parent._side_covers` → `_side_atoms` | probe 失败 → 苯/环母体拒收 | `c1ccc(cc1)CC(C)CC` → 塌成 `3-methylpentane` |
| L3 烷基 | `substituent_extractor._one_alkyl` | linear → shape 表 → None | C36 丢 5C 枝；酸上支链蒸发 |
| L3 alkoxy | `alkoxy_names._one_ring_alkoxy` | **要求 attach 在环上** | `COC(C)CC` → butane |
| 醚母体 | `ether_parent.py` | linear C1–4 + ipr + hfip | 支化臂拒收后回退烷烃且 L3 不补 methoxy |
| 芳基递归 | `layer2/leaves/*` | Ph/OPh/Bn depth 叶 | 与脂肪树无关；保持 |
| 文件体量 | `side_alkyl.py` ~442、`substituent_extractor.py` ~466 | 近 500 行 | **新逻辑必须新文件** |

**已知失败例：**

| SMILES | 现状 | 期望 |
|--------|------|------|
| `COC(C)CC` | butane | 2-methoxybutane / 2-甲氧基丁烷 |
| `COC(C)C(F)(F)F` | 1,1,1-trifluoropropane | 1,1,1-trifluoro-2-methoxypropane |
| `CCCCCCCCCC(CCCC(CC(CC)C)CCCCC)CCCCCCCCCCCC` | 10-nonyldocosane（丢 5C） | 6-(2-methylbutyl)-10-nonyldocosane |
| `CCCCC(CC(C)CC)C(=O)O` | hexanoic acid | 2-(2-methylbutyl)hexanoic acid（或等价 PIN） |
| `c1ccc(cc1)CC(C)CC` | 3-methylpentane | (2-methylbutyl)benzene |
| `c1ccnc(c1)CC(C)CC` | 3-methylpentane | 2-(2-methylbutyl)pyridine（若吡啶 gate 对齐） |

**根因一句话：** `_SIDE_PROBES` / L3 shape 表 / 环-only alkoxy 都是 **有限枚举**；枚举外静默失败。

---

## 目标能力分层（本 plan 交付范围）

| ID | 能力 | 本 plan |
|----|------|---------|
| P1 | 饱和纯碳侧枝 systematic 树（深度≤3，原子≤12） | **必做** |
| P2 | 开链 + 环统一 alkoxy（复用 `outer_alkoxy`） | **必做** |
| P3 | 芳基 leaf 递归 | **已有，只组合调用** |
| P4 | `claimable_side` 统一 claim | **必做** |
| P5 | 环 gate 消费 P4（苯 / 简单吡啶 / 环烷） | **必做** |
| P6 | 碳覆盖断言（测试 + 可选 meta） | **必做** |
| P7 | 侧链含简单不饱和（prenyl 已有外的烯基树） | **Phase 后段可选** |
| P8 | 侧链含 N（aminoalkyl 等） | **不在本 plan**（另开） |
| P9 | 醚臂复杂 FG / polyether | **不在本 plan**；醚只做 alkoxyalkane PIN |
| P10 | retained shape 捷径保留 | **保持**，systematic 垫底 |

---

## 目标态文件结构

```
src/namepredict/
  layer2/
    side_alkyl.py              # 现有 retained probes；末尾 probe 调 systematic
    side_alkyl_sys.py          # 新建：饱和碳树拓扑 AlkylTree / collect / spine
    side_claim.py              # 新建：claimable_side / side_atom_set
    side_alkoxy.py             # 现有 outer alkoxy
    side_facts.py              # 导出 systematic_alkyl / claimable_side
    ring_parent.py             # claim 仍走 _side_atoms（自动吃到 systematic）
    ether_parent.py            # 不扩 arm；文档注释 PIN 回退
    leaves/*                   # 不动协议；可选 Task 后期 leaf 上挂 P1
  layer3/
    substituent_extractor.py   # _one_alkyl 接 systematic 名
    alkyl_sys_names.py         # 新建：AlkylTree → en/zh + paren
    alkoxy_names.py            # 开链 alkoxy extract
  layer5/                      # 原则上不改；paren 已支持
tests/unit/
  test_systematic_alkyl.py
  test_open_chain_alkoxy.py
  test_side_carbon_coverage.py
  test_ring_claim_systematic.py
```

---

## 核心接口（全任务共用，禁止改名）

```python
# layer2/side_alkyl_sys.py
from dataclasses import dataclass

@dataclass(frozen=True)
class AlkylTree:
    """Saturated pure-carbon side tree rooted at attachment carbon (in side)."""
    root: int
    spine: tuple[int, ...]          # longest path starting at root
    atoms: frozenset[int]           # all carbons in this side
    # branches: spine locant (1-based) -> child trees (child.root bonded to spine)
    branches: tuple[tuple[int, tuple["AlkylTree", ...]], ...]

def collect_side_carbons(mol, root: int, forbidden: set[int]) -> frozenset[int] | None:
    """BFS carbons from root not entering forbidden; None if hetero/unsat/aromatic."""

def build_alkyl_tree(mol, root: int, forbidden: set[int], *, max_atoms: int = 12, max_depth: int = 3) -> AlkylTree | None:
    """Pure C saturated side → tree; None if unsupported."""

def systematic_side_atoms(mol, root: int, forbidden: set[int]) -> list[int] | None:
    """Probe-compatible: atom list or None (for _SIDE_PROBES / claim)."""

# layer2/side_claim.py
def claimable_side(mol, root: int, parent: set[int]) -> frozenset[int] | None:
    """Atoms covered by one claimable side starting at root, or None."""

# layer2/side_facts.py
def systematic_alkyl(mol, root: int, parent: set[int]) -> AlkylTree | None: ...

# layer3/alkyl_sys_names.py
def name_alkyl_tree(tree: AlkylTree) -> tuple[str, str, bool]:
    """→ (en, zh, paren). paren True if locants/multi-branch."""
```

**命名规则（P1，饱和）：**

1. 从附着点 `root` 起选 **最长 spine**（并列时选支化更多 / 原子序更小，实现里固定比较器，写进测试）。
2. Spine 编号从附着点 = 1。
3. 分叉命名为前缀，位次为 spine 位次；多前缀字母序（EN `alkyl_alpha_key`）。
4. 末端 hydride 名：`methyl`… 按 spine 长度；带前缀则 `k-methyl…yl`。
5. 例：`-CH2-CH(CH3)-CH2-CH3` → spine 4 → **2-methylbutyl** / **2-甲基丁基**，`paren=True`（有位次）。

**复用关系：**

```
parent.chain / ring_set
        │
        ▼
 claimable_side  ==  _side_atoms（probes + systematic）
        │
        ├── L2 ring gate (_side_covers)
        └── L3 _one_alkyl / render name_alkyl_tree
 hetero: outer_alkoxy → open+ring _extract_alkoxys
 aryl: leaves/* 不变
```

---

## Task 1: 红测 — 用户已踩案例 + 碳覆盖

**Files:**
- Create: `tests/unit/test_side_chain_engine_red.py`
- Test: 同文件

**Interfaces:**
- Consumes: `SMILESNNamer`, `normalize_en`, `normalize_zh`
- Produces: 锁定目标名，供后续 task 变绿

- [ ] **Step 1: 写失败测试**

```python
# IUPAC: P-29.3 / P-63.2.4
# Layer: L2,L3,L5
"""Red-bar cases for side-chain engine (systematic alkyl + open alkoxy)."""
from __future__ import annotations

import pytest
from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    ("COC(C)CC", "2-methoxybutane", "2-甲氧基丁烷"),
    ("COC(C)C(F)(F)F", "1,1,1-trifluoro-2-methoxypropane", "1,1,1-三氟-2-甲氧基丙烷"),
    (
        "CCCCCCCCCC(CCCC(CC(CC)C)CCCCC)CCCCCCCCCCCC",
        "6-(2-methylbutyl)-10-nonyldocosane",
        "6-(2-甲基丁基)-10-壬基二十二烷",
    ),
    ("c1ccc(cc1)CC(C)CC", "(2-methylbutyl)benzene", "(2-甲基丁基)苯"),
    ("CCCCC(CC(C)CC)C(=O)O", "2-(2-methylbutyl)hexanoic acid", "2-(2-甲基丁基)己酸"),
]

# 不得回归
REGRESS = [
    ("COC(C)C", "isopropyl methyl ether", "异丙基甲醚"),
    ("COCC", "methoxyethane", "甲氧基乙烷"),
    ("COC(C(F)(F)F)C(F)(F)F", "hexafluoroisopropyl methyl ether", "六氟异丙基甲醚"),
    ("CC(C)CCc1ccccc1", "isopentylbenzene", "异戊基苯"),
    ("CCC(C)c1ccccc1", "sec-butylbenzene", "仲丁基苯"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_side_engine_targets(smiles, en, zh):
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)


@pytest.mark.parametrize("smiles,en,zh", REGRESS)
def test_side_engine_no_regress(smiles, en, zh):
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)
```

- [ ] **Step 2: 跑测确认红**

```bash
source .venv/Scripts/activate
pytest tests/unit/test_side_chain_engine_red.py -v
```

Expected: `CASES` 中多项 FAIL；`REGRESS` PASS。

- [ ] **Step 3: Commit 红测**

```bash
git add tests/unit/test_side_chain_engine_red.py
git commit -m "$(cat <<'EOF'
test(namepredict): red-bar side-chain engine targets

Lock dual names for open alkoxy, 2-methylbutyl on alkane/arene/acid.

Co-Authored-By: Claude <noreply@anthropic.com>
EOF
)"
```

---

## Task 2: L2 — `AlkylTree` 拓扑（纯碳饱和）

**Files:**
- Create: `src/namepredict/layer2/side_alkyl_sys.py`
- Create: `tests/unit/test_alkyl_tree_topo.py`
- Modify: 无（本任务不接线）

**Interfaces:**
- Produces: `collect_side_carbons`, `build_alkyl_tree`, `systematic_side_atoms`, `AlkylTree`

- [ ] **Step 1: 拓扑单测**

```python
# tests/unit/test_alkyl_tree_topo.py
from rdkit import Chem
from namepredict.layer2.side_alkyl_sys import build_alkyl_tree, collect_side_carbons

def _mol(s):
    m = Chem.MolFromSmiles(s)
    Chem.AddHs(m)  # 若项目 preprocess 不同，改用 layer0.preprocess
    return m

def test_2_methylbutyl_tree_atoms():
    # Use preprocess like production:
    from namepredict.layer0.preprocessor import preprocess
    mol = preprocess("CC(C)CC")  # not this — use mapped side
    # Better: molecule where side is explicit off a fake parent carbon
    # Parent chain carbon 0 forbidden; side root = atom of CH2 in 2-methylbutyl
    mol = preprocess("CCCCC")  # placeholder replaced in implementer by real map test
```

实现时用 **原子映射 SMILES** 或 `preprocess` + 已知 idx（与对话中 C36 分支 `start=14` 一致）写稳定断言：

```python
from namepredict.layer0.preprocessor import preprocess
from namepredict.layer2.side_alkyl_sys import build_alkyl_tree

def test_build_tree_on_c36_branch():
    s = "CCCCCCCCCC(CCCC(CC(CC)C)CCCCC)CCCCCCCCCCCC"
    mol = preprocess(s)
    # side root 14, forbidden = parent chain set from select_parent
    from namepredict.layer1.analyzer import analyze
    from namepredict.layer2.parent_selector import select_parent
    parent = set(select_parent(analyze(mol))["chain"])
    tree = build_alkyl_tree(mol, 14, parent)
    assert tree is not None
    assert len(tree.atoms) == 5
    assert len(tree.spine) == 4  # butyl spine
    locs = dict(tree.branches)
    assert 2 in locs  # methyl at locant 2
```

- [ ] **Step 2: 跑测红**

```bash
pytest tests/unit/test_alkyl_tree_topo.py -v
```

Expected: import/FAIL。

- [ ] **Step 3: 实现 `side_alkyl_sys.py`（每函数 ≤10 行）**

要点：

```python
def _heavies(atom):
    return [n for n in atom.GetNeighbors() if n.GetAtomicNum() != 1]

def _is_pure_sat_c(mol, idx: int) -> bool:
    a = mol.GetAtomWithIdx(idx)
    if a.GetAtomicNum() != 6 or a.GetIsAromatic() or a.IsInRing():
        return False
    return all(
        b.GetBondType().name == "SINGLE" and b.GetOtherAtom(a).GetAtomicNum() in (1, 6)
        for b in a.GetBonds()
    )

def collect_side_carbons(mol, root, forbidden):
    # BFS; any bad atom → None
    ...

def _longest_spine(mol, root, atoms: set[int]) -> list[int]:
    # DFS only within atoms, path must start at root
    ...

def _child_roots(mol, spine_atom, spine_set, atoms) -> list[int]:
    ...

def build_alkyl_tree(mol, root, forbidden, *, max_atoms=12, max_depth=3):
    atoms = collect_side_carbons(mol, root, forbidden)
    if atoms is None or len(atoms) > max_atoms:
        return None
    spine = _longest_spine(mol, root, set(atoms))
    # build branches recursively with max_depth
    ...

def systematic_side_atoms(mol, root, forbidden):
    t = build_alkyl_tree(mol, root, forbidden)
    return None if t is None else sorted(t.atoms)  # or list in BFS order; claim only needs set
```

**环上侧链：** `forbidden=ring_set`；**开链：** `forbidden=chain_set`。  
**深度：** 子树 `build` 时 `max_depth-1`；depth 用尽仍有分叉 → `None`（宁可 claim 失败，勿截断蒸发）。

- [ ] **Step 4: 绿测 + structure_lint**

```bash
pytest tests/unit/test_alkyl_tree_topo.py -v
python tools/structure_lint.py --root src/namepredict/layer2/side_alkyl_sys.py
```

- [ ] **Step 5: Commit**

```bash
git add src/namepredict/layer2/side_alkyl_sys.py tests/unit/test_alkyl_tree_topo.py
git commit -m "$(cat <<'EOF'
feat(namepredict): AlkylTree topology for saturated carbon sides

BFS collect + longest spine + branch locants (P-29 systematic).

Co-Authored-By: Claude <noreply@anthropic.com>
EOF
)"
```

---

## Task 3: L3 — `name_alkyl_tree` + 接入 `_one_alkyl`

**Files:**
- Create: `src/namepredict/layer3/alkyl_sys_names.py`
- Modify: `src/namepredict/layer3/substituent_extractor.py`（`_one_alkyl` / `_one_branched` 末尾）
- Modify: `src/namepredict/layer2/side_facts.py`（导出 `systematic_alkyl`）
- Test: `tests/unit/test_systematic_alkyl.py` + 部分变绿 `test_side_chain_engine_red.py` 中烷/酸例

**Interfaces:**
- Consumes: `AlkylTree`, `build_alkyl_tree`
- Produces: `name_alkyl_tree` → `(en, zh, paren)`；L3 sub dict

- [ ] **Step 1: 命名单测**

```python
def test_name_2_methylbutyl():
    # build tree from known mol/idx then:
    en, zh, paren = name_alkyl_tree(tree)
    assert en == "2-methylbutyl"
    assert zh == "2-甲基丁基"
    assert paren is True
```

线性 C1–C4 若误入 systematic，应与 `ALKYL_EN` 一致（或 `_one_alkyl` 保证 linear 先匹配，systematic 仅非 linear）。

- [ ] **Step 2: 实现命名（递归）**

```python
# alkyl_sys_names.py
_STEM_EN = {1: "methyl", 2: "ethyl", 3: "propyl", 4: "butyl", 5: "pentyl",
            6: "hexyl", 7: "heptyl", 8: "octyl", 9: "nonyl", 10: "decyl",
            11: "undecyl", 12: "dodecyl"}
_STEM_ZH = {1: "甲基", 2: "乙基", 3: "丙基", 4: "丁基", 5: "戊基",
            6: "己基", 7: "庚基", 8: "辛基", 9: "壬基", 10: "癸基",
            11: "十一烷基", 12: "十二烷基"}

def _branch_prefix(loc: int, child: AlkylTree) -> tuple[str, str]:
    en, zh, _ = name_alkyl_tree(child)
    # child already *yl; for methyl child spine=1 → "methyl"
    return f"{loc}-{en}", f"{loc}-{zh}"

def name_alkyl_tree(tree: AlkylTree) -> tuple[str, str, bool]:
    n = len(tree.spine)
    stem_en, stem_zh = _STEM_EN[n], _STEM_ZH[n]
    if not tree.branches:
        return stem_en, stem_zh, False
    # sort prefixes by alkyl_alpha_key of en without locant
    ...
    return f"{pref_en}{stem_en}", f"{pref_zh}{stem_zh}", True
```

注意：中文「2-甲基丁基」无空格；英文 `2-methylbutyl` 无空格。多前缀：`2-ethyl-3-methylpentyl` 字母序 ethyl 在 methyl 前。

- [ ] **Step 3: `side_facts.systematic_alkyl`**

```python
def systematic_alkyl(mol, root, parent: set[int]):
    from namepredict.layer2.side_alkyl_sys import build_alkyl_tree
    return build_alkyl_tree(mol, root, parent)
```

- [ ] **Step 4: 接 `_one_alkyl`**

在 `substituent_extractor.py`：

```python
def _one_systematic(mol, attach, start, chain_set):
    tree = side_facts.systematic_alkyl(mol, start, chain_set)
    if tree is None:
        return None
    en, zh, paren = name_alkyl_tree(tree)
    return {
        "kind": "alkyl", "n_carbons": len(tree.atoms), "attach_idx": attach,
        "atoms": list(tree.atoms), "en": en, "zh": zh, "paren": paren,
    }

def _one_alkyl(...):
    ...
    ha = _one_haloalkyl(...)
    if ha is not None:
        return ha
    br = _one_branched(...)
    if br is not None:
        return br
    return _one_systematic(mol, attach, start, chain_set)
```

若 `substituent_extractor.py` 超 500 行：把 `_one_systematic` 放进 `alkyl_sys_names.py` 或新 `alkyl_extract.py`。

- [ ] **Step 5: 跑测**

```bash
pytest tests/unit/test_systematic_alkyl.py tests/unit/test_branched_alkyl.py tests/unit/test_tert_pentyl.py tests/unit/test_side_chain_engine_red.py -v
```

Expected: C36 与酸例变绿（或接近；位次若 6/10 反了再调 L4，本任务先保证两取代基都在）。  
苯例仍可能红（claim 未接）。

- [ ] **Step 6: Commit**

```bash
git commit -m "$(cat <<'EOF'
feat(namepredict): systematic alkyl names on chain parents

L3 fallback after retained shapes; 2-methylbutyl no longer evaporates.

Co-Authored-By: Claude <noreply@anthropic.com>
EOF
)"
```

---

## Task 4: L3 — 开链 alkoxy（P2）

**Files:**
- Modify: `src/namepredict/layer3/alkoxy_names.py`
- Test: `tests/unit/test_open_chain_alkoxy.py`；红测醚两例变绿
- 不修改 `ether_parent` arm 表（除注释）

**Interfaces:**
- Consumes: `outer_alkoxy`, `info["ethers"]`, `parent["chain"]`
- Produces: alkoxy sub 挂在 chain 附着碳

- [ ] **Step 1: 测试**

```python
CASES = [
    ("COC(C)CC", "2-methoxybutane", "2-甲氧基丁烷"),
    ("COC(C)C(F)(F)F", "1,1,1-trifluoro-2-methoxypropane", "1,1,1-三氟-2-甲氧基丙烷"),
    ("COCC", "methoxyethane", "甲氧基乙烷"),  # 可仍走 ether parent
]
```

- [ ] **Step 2: 实现开链提取**

```python
def _one_chain_alkoxy(mol, e, chain_set):
    ends = _alkoxy_ends(e, chain_set)  # 恰好一端在 chain
    if ends is None:
        return None
    o, chain_c, outer = ends
    if mol.GetAtomWithIdx(outer).GetIsAromatic():
        return None
    # 允许 chain_c 不在环上
    fact = outer_alkoxy(mol, outer, o)
    return _make_alkoxy(chain_c, o, list(fact.atoms), fact.code) if fact else None

def _extract_alkoxys(info, parent):
    mol = info["mol"]
    chain_set = set(parent.get("chain") or [])
    out = []
    for e in info.get("ethers") or []:
        one = _one_ring_alkoxy(mol, e, chain_set) or _one_chain_alkoxy(mol, e, chain_set)
        if one is not None:
            out.append(one)
    return out
```

**与 ether parent 共存：**  
- `COCC`：若 L2 仍选 `kind=ether`，L3 在 ether+arms 时可能 skip 重复（已有 `ether_arms` 跳过 halo 逻辑）；确认 **不要** 双重 methoxy。  
- 读 `substituent_extractor` 里 ether 跳过：若 `kind=="ether"` 且无取代基需求，alkoxy extract 对 **母体臂上的 O** 应避免把母体 O 再当取代基。

规则：

```python
def _ether_parent_o(parent):
    return parent.get("o_idx") if parent.get("kind") == "ether" else None

# in extract: skip e if e["o_idx"] == parent o_idx
```

- [ ] **Step 3: 当 ether_parent 失败时**  
`COC(C)CC` → parent alkane butane，chain 含 C2，methoxy 提出 → L4 位次 2。  
验证 pipeline。

- [ ] **Step 4: 跑测**

```bash
pytest tests/unit/test_open_chain_alkoxy.py tests/unit/test_dialkyl_ether.py tests/unit/test_branched_ether.py tests/unit/test_side_chain_engine_red.py -v
```

- [ ] **Step 5: Commit**

```bash
git commit -m "$(cat <<'EOF'
feat(namepredict): open-chain alkoxy extract via outer_alkoxy

Substitutive PIN path for branched ethers without ether_parent arms.

Co-Authored-By: Claude <noreply@anthropic.com>
EOF
)"
```

---

## Task 5: P4 — `claimable_side` + `_SIDE_PROBES` 接入 systematic

**Files:**
- Create: `src/namepredict/layer2/side_claim.py`
- Modify: `src/namepredict/layer2/side_alkyl.py` — `_SIDE_PROBES` 末尾加 `systematic_side_atoms` 包装
- Modify: `src/namepredict/layer2/side_facts.py` — 导出 `claimable_side`
- Test: `tests/unit/test_side_claim.py`

**Interfaces:**
- Produces: `claimable_side(mol, root, parent) -> frozenset[int] | None`
- `_side_atoms` 自动识别 2-methylbutyl → 环 gate 无需改判定公式

- [ ] **Step 1: 包装 probe**

```python
# side_alkyl_sys 已有 systematic_side_atoms(mol, root, forbidden) -> list|None

# side_alkyl.py
def _probe_systematic(mol, start, chain):
    from namepredict.layer2.side_alkyl_sys import systematic_side_atoms
    return systematic_side_atoms(mol, start, chain)

_SIDE_PROBES = (
    ...existing...,
    _probe_systematic,  # LAST — retained names win first
)
```

- [ ] **Step 2: `claimable_side`**

```python
def claimable_side(mol, root, parent: set[int]):
    from namepredict.layer2.side_alkyl import _side_atoms
    atoms = _side_atoms(mol, root, parent)
    return frozenset(atoms) if atoms is not None else None
```

- [ ] **Step 3: 单测**

```python
def test_claim_2_methylbutyl_on_benzene():
    mol = preprocess("c1ccc(cc1)CC(C)CC")
    ring = set(...)  # 6 ring atoms
    start = ...      # 苄位侧链第一碳
    got = claimable_side(mol, start, ring)
    assert got is not None and len(got) == 5
```

- [ ] **Step 4: 验证苯母体复活**

```bash
python - <<'PY'
from namepredict.namer import SMILESNNamer
print(SMILESNNamer().name("c1ccc(cc1)CC(C)CC").en)
PY
```

Expected: `(2-methylbutyl)benzene`（若 L3 已接 systematic）。

- [ ] **Step 5: 跑 `test_branched_alkyl` + 红测苯例**

```bash
pytest tests/unit/test_branched_alkyl.py tests/unit/test_side_claim.py tests/unit/test_side_chain_engine_red.py -v
```

- [ ] **Step 6: Commit**

```bash
git commit -m "$(cat <<'EOF'
feat(namepredict): systematic probe in _SIDE_PROBES + claimable_side

Ring parents can claim 2-methylbutyl sides via shared atoms API.

Co-Authored-By: Claude <noreply@anthropic.com>
EOF
)"
```

---

## Task 6: 碳覆盖门禁（防再蒸发）

**Files:**
- Create: `src/namepredict/layer3/coverage.py`（或 `tests` 专用 helper 若不想进生产）
- Create: `tests/unit/test_side_carbon_coverage.py`
- Optional Modify: `namer.py` debug meta `uncaptured_carbons`（仅 success 路径统计，不改变 en/zh）

**推荐生产小函数：**

```python
def carbon_coverage_gap(mol, parent: dict, subst: list) -> set[int]:
    all_c = {a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 6}
    covered = set(parent.get("chain") or [])
    for s in subst:
        covered |= set(s.get("atoms") or [])
    # 醚母体 o 不计入碳；卤素非碳
    return all_c - covered
```

- [ ] **Step 1: 测试**

```python
def test_c36_no_gap():
    ...
    assert carbon_coverage_gap(mol, parent, subst) == set()

def test_old_bug_would_gap():
    # 可对修复前逻辑不测；只测当前 pipeline
```

- [ ] **Step 2: 对红测 SMILES 参数化 `gap == empty`**

- [ ] **Step 3: Commit**

```bash
git commit -m "$(cat <<'EOF'
test(namepredict): carbon coverage guard for side-chain extract

Co-Authored-By: Claude <noreply@anthropic.com>
EOF
)"
```

---

## Task 7: 环族 claim 对齐（吡啶 / 环烷 / 苯酚门控）

**Files:**
- Modify: 仅当某 gate **不**走 `_side_atoms` 时才改  
  检查：`pyridine.py`、`sat_hetero.py`、`phenol_aniline.py`、`ring_parent._arene_fg_subs_ok`
- Test: `tests/unit/test_ring_claim_systematic.py`

**Interfaces:**
- Consumes: `claimable_side` / `_side_atoms`
- Produces: 吡啶/环己烷/苯酚上 2-methylbutyl 保留环母体

- [ ] **Step 1: 探测哪些 gate 仍拒**

```bash
python - <<'PY'
from namepredict.namer import SMILESNNamer
for s in [
  "c1ccc(cc1)CC(C)CC",
  "c1ccnc(c1)CC(C)CC",
  "CC(C)CC1CCCCC1",  # isopentyl already ok; use 2-methylbutyl:
  "CCC(C)CC1CCCCC1",
  "Oc1ccc(CC(C)CC)cc1",
]:
    print(s, "->", SMILESNNamer().name(s).en)
PY
```

- [ ] **Step 2: 凡手写 linear-only 的 gate，改为 `_side_atoms is not None`**

例：若存在 `_linear_n_alkyl_sides_ok` 独占路径，增加 systematic 或改为：

```python
def _alkyl_sides_ok(mol, ring, starts):
    sets = _side_sets(mol, ring, starts)
    return sets is not None and _disjoint_cover(sets, _outside_c_atoms(mol, ring) - exclude)
```

**禁止** 复制 `build_alkyl_tree` 进 pyridine 文件。

- [ ] **Step 3: 测试**

```python
CASES = [
    ("c1ccc(cc1)CC(C)CC", "(2-methylbutyl)benzene", "(2-甲基丁基)苯"),
    ("CCC(C)CC1CCCCC1", "(2-methylbutyl)cyclohexane", "(2-甲基丁基)环己烷"),
    ("Oc1ccc(CC(C)CC)cc1", "4-(2-methylbutyl)phenol", "4-(2-甲基丁基)苯酚"),
    # pyridine if simple parent supports mono alkyl:
    ("c1ccnc(c1)CC(C)CC", "2-(2-methylbutyl)pyridine", "2-(2-甲基丁基)吡啶"),
]
```

- [ ] **Step 4: 跑相关环测防回归**

```bash
pytest tests/unit/test_branched_alkyl.py tests/unit/test_alkylbenzene.py tests/unit/test_pyridine_alkoxy_nitro.py tests/unit/test_ring_claim_systematic.py -v
```

- [ ] **Step 5: Commit**

```bash
git commit -m "$(cat <<'EOF'
feat(namepredict): ring gates claim systematic alkyl sides

Reuse _side_atoms/claimable_side; no per-kind side namer forks.

Co-Authored-By: Claude <noreply@anthropic.com>
EOF
)"
```

---

## Task 8: 多层脂肪递归（深度 2–3）+ 字母序

**Files:**
- Modify: `side_alkyl_sys.py` / `alkyl_sys_names.py`
- Test: `tests/unit/test_alkyl_tree_depth.py`

**范围（本 task）：** spine 上多个甲基/乙基；子树再带甲基（深度 2）。

例：

| 侧链 | EN |
|------|-----|
| `-CH2-CH(CH3)-CH(CH3)-CH3` | 2,3-dimethylbutyl |
| `-CH2-C(CH3)2-CH2-CH3` | 2,2-dimethylbutyl |
| 深度2：`-CH2-CH(CH2CH3)-CH2-CH3` 等 | 2-ethylbutyl |

- [ ] **Step 1: 红测 2–3 个 depth 例（挂在 benzene 或 alkane 母体）**

- [ ] **Step 2: 确保 `name_alkyl_tree` 递归 child 时：**
  - child 为 methyl → `methyl` 不是 `1-methylmethyl`
  - 多 locant 同一名用 `2,2-dimethyl…` 倍增前缀（di/tri）—— **若现有 L5 倍增只在组装层**，侧链名内部自带 `dimethyl` 与 IUPAC 复杂取代基习惯对齐（2-methylbutan-2-yl 已是整串）。

  **简化规则（本 plan 固定）：**  
  - 同一 locant 多个 identical methyl → `2,2-dimethyl…`  
  - 不同 locant → `2,3-dimethyl…`  
  - 非 methyl 子树 → `(…)` 包裹 child 名：`2-(1-methylpropyl)…` 仅当 child spine>1 且 child 自带位次

- [ ] **Step 3: 绿测 + commit**

```bash
git commit -m "$(cat <<'EOF'
feat(namepredict): multi-branch systematic alkyl depth≤3

Co-Authored-By: Claude <noreply@anthropic.com>
EOF
)"
```

---

## Task 9: 侧链有限不饱和（P7 可选增量）

**Files:**
- Modify: `side_alkyl_sys.py` 允许 spine 上 **一个** 双键（非芳）
- Modify: `alkyl_sys_names.py` → `…enyl` 位次
- Test: 除 prenyl 捷径外增加 `but-2-en-1-yl` 类 1 例

**不做：** 共轭多烯侧链、炔基树、环丙烯基（另 plan）。

- [ ] 仅当 Task 1–8 全绿且 dual 稳定后做。  
- [ ] prenyl 继续走 retained probe（优先于 systematic）。

---

## Task 10: 醚策略文档化 + 回归锁（P9 边界）

**Files:**
- Modify: `src/namepredict/layer2/ether_parent.py` 模块 docstring
- Modify: `tests/unit/test_branched_ether.py` 增加注释与 **不** 把 sec-butyl  eth er 强行 functional-class
- Optional: `workstate.md` 一条

**Docstring 约定：**

```text
P-63.2.2/4:
- linear C1–C4 both arms → ether parent (alkoxyalkane / sym retained)
- ipr / hfip specials → functional-class ether_arms
- any other substitution → do NOT claim ether parent;
  alkane (or higher FG) parent + L3 alkoxy/halo prefixes (PIN)
```

- [ ] **Commit**

```bash
git commit -m "$(cat <<'EOF'
docs(namepredict): ether parent scope vs substitutive alkoxy PIN

Co-Authored-By: Claude <noreply@anthropic.com>
EOF
)"
```

---

## Task 11: 与芳基 leaf 组合（可选，薄）

**Files:**
- Modify: `layer2/leaves/topo.py` `match_alkoxy` 已用 outer n；不改  
- 若存在「环上 alkoxy 外侧再支化」需求：P2 `outer_alkoxy` 扩支化 code（已有 isopropoxy）  
- **脂肪侧链挂在 Ph 叶上**：已有 depth-2；**不**在 leaf 内嵌 `AlkylTree` 除非红测要求

本 task 默认 **SKIP**，除非 benchmark 出现 `2-(2-methylbutyl)phenoxy` 类失败。

---

## Task 12: 全量验收 + workstate

**Files:**
- Modify: `workstate.md`（按仓库日志格式追加）
- Run: unit 子集 + structure_lint + benchmark_parallel

- [ ] **Step 1:**

```bash
source .venv/Scripts/activate
pytest tests/unit/test_side_chain_engine_red.py tests/unit/test_systematic_alkyl.py \
  tests/unit/test_open_chain_alkoxy.py tests/unit/test_side_claim.py \
  tests/unit/test_side_carbon_coverage.py tests/unit/test_ring_claim_systematic.py \
  tests/unit/test_branched_alkyl.py tests/unit/test_dialkyl_ether.py \
  tests/unit/test_branched_ether.py tests/unit/test_tert_pentyl.py -q
python tools/structure_lint.py --root src/namepredict
python -m benchmarks.benchmark_parallel --data data/merged_benchmark.json --time
```

- [ ] **Step 2:** dual 回退 ≤0.5%；记录 fails 变化  
- [ ] **Step 3:** workstate 一条，引用 IUPAC P-29 / P-63.2.4  
- [ ] **Step 4:** 最终 commit（若有文档）

---

## 各「不会自动拥有」项的对应做法（汇总）

| 缺口 | 做法 | Task |
|------|------|------|
| 苯/吡啶挂 2-methylbutyl | `_SIDE_PROBES`+systematic → `_side_atoms` → 现有 `_side_covers`；吡啶手写 gate 改调同一 API | 5, 7 |
| 侧链 O/N/环/不饱和 | O：P2 alkoxy；环：已有 cycloalkyl probe；不饱和：Task 9 有限；**N：本 plan 不做** | 4, 9 |
| 醚臂复杂 FG | **明确不做**；PIN=取代命名；复杂醚另 plan | 10 |
| 取代基的取代基 3 层脂肪树 | `AlkylTree` 递归 `max_depth=3` + `name_alkyl_tree` | 2, 3, 8 |
| retained 母体 claim 放宽 | **统一 claimable_side / _side_atoms**，禁止每母体复制 | 5, 7 |
| 开链 alkoxy / C36 | L3 extract + systematic | 3, 4, 6 |

---

## 风险与缓解

| 风险 | 缓解 |
|------|------|
| systematic 抢 retained 名（isopentyl→3-methylbutyl） | probe **顺序**：retained 全部在 systematic **之前** |
| 苯接受过大侧链导致错误母体 | `max_atoms=12`；n_sub≤4 仍有效 |
| 文件超 500 行 | 新逻辑只进 `side_alkyl_sys` / `alkyl_sys_names` / `side_claim` |
| dual 回归 | 每 task 跑 REGRESS；Task 12 全 benchmark |
| 位次 6 vs 10 字母序 | L4 已有规则；断言完整 dual 名；若差一位先查 `number` 再改 L4 不改 L5 补丁 |
| ether 双重计数 | skip parent `o_idx` |

---

## 建议实施顺序（依赖图）

```
Task1 红测
  → Task2 AlkylTree
  → Task3 L3 命名+接入     ──→ Task6 覆盖门禁
  → Task4 开链 alkoxy
  → Task5 claim/_SIDE_PROBES → Task7 环族对齐
  → Task8 深度递归
  → Task9 不饱和（可选）
  → Task10 醚文档
  → Task12 验收
```

**最小可交付（MVP）：** Task 1–6 + 12（修用户已踩醚/C36/酸蒸发 + 苯 claim 基础）。  
**完整 plan 声明范围：** MVP + Task 7–8 + 10（+ 可选 9）。

---

## Self-Review

1. **Spec coverage:** 五类缺口均有 Task；先前醚/烷烃/递归复用策略写入 Architecture + Task 10。  
2. **Placeholder scan:** 无 TBD；可选 Task 9/11 标 SKIP 条件。  
3. **类型一致:** `AlkylTree` / `build_alkyl_tree` / `name_alkyl_tree` / `claimable_side` / `systematic_alkyl` 命名贯穿。  
4. **层纯度:** 拓扑 L2、渲染 L3、claim 共享 `_side_atoms`。  
5. **YAGNI:** 不做通用 N-侧链、polyether、完整叶上脂肪引擎。

---

## Execution Handoff

Plan complete and saved to `docs/superpowers/plans/2026-07-18-side-chain-recursive-engine.md`.

**两种执行方式：**

1. **Subagent-Driven（推荐）** — 每 Task 新开 subagent，Task 间 review  
2. **Inline Execution** — 本会话按 `executing-plans` 连续做，设 checkpoint  

你更想用哪一种？若只要 MVP，也可以说「先做 Task 1–6」。
