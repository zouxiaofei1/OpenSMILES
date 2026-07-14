# Complete CIP Stereo Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 分阶段把 NamePredict 的立体化学从「RDKit 读标签 + 白名单前缀」推进到可审计、可扩展的完整 CIP 管线（四面体 R/S、双键 E/Z、假不对称、轴向/面性、相对构型），每阶段 dual 可验收且不破层纯度。

**Architecture:** 立体分两层职责——**(A) 赋值/排序**（谁决定 R 还是 S）与 **(B) 名称前缀**（如何写成 `(2S,4R)-` / `(E,8R)-`）。近期 (A) 继续委托 RDKit，但经统一 **StereoProvider** 接口；远期用自研 digraph 替换 provider 实现而不改 L5。E/Z 与 R/S 统一为 `list[StereoPart]`，在 L5 末尾合并一次。编号层 (L4) 仅在「立体影响最低位次集」时消费同一接口。

**Tech Stack:** Python 3.12、RDKit（`AssignStereochemistry` / `BondStereo` / 远期 digraph）、pytest、现有 `SMILESNNamer` 端到端 dual 测试、`benchmarks.benchmark_parallel`。

## Global Constraints

- 代码只能落在 `src/namepredict/layer0`–`layer5` 或入口点；其他目录只读（`prompt.txt`）。
- 单文件 ≤500 行；每个函数 ≤10 行；函数式、尽量单输入/单输出。
- 不得混层：赋值可放 L1（结构事实），前缀组装只在 L5；L2 不写立体字母。
- 禁止深度学习；禁止大缓存；禁止改 benchmark 数据。
- 每任务结束跑相关 unit tests；触及命名输出时跑 `python -m benchmarks.benchmark_parallel --data data/merged_benchmark.json --time`；dual 回退不得 >0.5%。
- 每任务一个改进点；review PASS 后再 commit；`workstate.md` 按既有日志格式记一条。
- IUPAC 依据：`docs/iupac/` 与 `docs/iupac/cn/第7章_立体化学.md`（P-91～P-93 族）。
- 不自研完整 digraph 直到 Phase D；Phase A–C 必须 dual 正收益。

---

## 现状基线（实施前必读）

| 文件 | 职责 | 缺口 |
|---|---|---|
| `src/namepredict/layer5/stereo_rs.py` | 母体链 R/S 前缀；与已有 E/Z 字符串合并 | `_RS_KINDS` 仅 8 种；只扫 `parent.chain` |
| `src/namepredict/layer5/unsat_acid.py` | 烯酸/烯醇 E/Z（单/多键） | 未挂 alkene/alkenal/alkenoate/alkenenitrile |
| `src/namepredict/layer5/assembler.py` | `assemble` 末尾 `apply_rs_prefix` | E/Z 分散在各 `_names_*` |
| `tests/unit/test_rs_stereo.py` | 酸/醇/胺 + 一条 E+R/S | 无酮/酯/腈/简单烯烃 |

当前 dual 约 **12.8%**（日志 `#556ed16` R/S 后）；实施前先跑一轮 benchmark 记入 workstate。

---

## 文件结构（目标态）

```
src/namepredict/
  layer1/
    stereo_types.py      # StereoPart / StereoKind 数据类（无 RDKit 赋值逻辑可放此；或仅类型）
    stereo_provider.py   # Protocol + RdkitStereoProvider（≤500 行，函数≤10）
  layer5/
    stereo_format.py     # 解析/合并/格式化前缀 (E,8R)-  （从 stereo_rs 抽出）
    stereo_rs.py         # apply_rs_prefix → 调 provider + format（变薄）
    stereo_ez.py         # 统一 E/Z 收集（从 unsat_acid 抽出可复用部分）
    unsat_acid.py        # 只保留酸名组装，E/Z 调 stereo_ez
    assembler.py         # assemble: names → ez/rs 统一 apply
tests/unit/
  test_rs_stereo.py      # 扩母体白名单
  test_ez_stereo.py      # 新建：alkene/alkenal/… E/Z
  test_stereo_format.py  # 纯字符串合并
  test_stereo_provider.py
  test_cip_digraph.py    # Phase D+
  test_axial_planar.py   # Phase E
  test_relative_stereo.py# Phase F
```

**层边界：**

- L1 `StereoProvider`：输入 `Mol` + 关注原子/键索引 → 输出 `R|S|E|Z|…`
- L5：只做 locant 映射与字符串；不调用 `AssignStereochemistry` 以外的化学规则（经 provider）
- L4（后期）：编号候选评分可读「立体是否存在」，不发明描述符

---

## 阶段总览

| Phase | 名称 | dual 意图 | 是否自研 CIP |
|---|---|---|---|
| **A** | 覆盖扩展（白名单 + E/Z 挂载） | 主涨分 | 否，RDKit |
| **B** | 统一 Stereo 模型与 Provider | 持平/微涨，可扩展 | 否 |
| **C** | 取代基手性 + 环上 R/S | 中涨 | 否 |
| **D** | 可审计 CIP digraph（四面体） | 正确性；dual 防回归 | **是** |
| **E** | 轴向 / 面性 (M/P) | 小涨（长尾） | 是 |
| **F** | 相对/外消旋/假不对称 | 小涨 | 是 |
| **G** | L4 立体参与编号 | 边界修正 | 消费 D |

**执行原则：** 一次只做一个 Task；A 做完再开 B。Phase D 前若 dual 仍主要死在母体选择，优先插队非立体选题，不阻塞本计划文档。

---

### Task 1: 基线 benchmark + 立体 fail 桶

**Files:**
- Create: `tools/_stereo_fail_bucket.py`
- Modify: none（只读 benchmark 输出）
- Test: 脚本可运行即可

**Interfaces:**
- Consumes: `python -m benchmarks.benchmark_parallel` 的 JSON/文本输出习惯
- Produces: `tools/_stereo_fail_buckets.json` — `{bucket: [smiles, ...]}`

- [ ] **Step 1: 跑 dual 基线**

```bash
PYTHONPATH=src python -m benchmarks.benchmark_parallel --data data/merged_benchmark.json --time
```

Expected: 打印 dual 百分比与 ok 条数；记下 `dual_pct` / `ok_dual` 到工作笔记（暂不改 workstate 除非本轮有代码）。

- [ ] **Step 2: 写 fail 分桶脚本**

```python
# tools/_stereo_fail_bucket.py
"""Cluster dual fails that look stereo-related (heuristic)."""
from __future__ import annotations
import json, re, sys
from pathlib import Path

# 启发式：金标含 (R|S|E|Z|r|s) 而预测缺；或预测有但字母不同
_ST = re.compile(r"\((?:\d*[EZezRSrs],?)+\)-")

def load_cases(path: Path) -> list[dict]:
    data = json.loads(path.read_text(encoding="utf-8"))
    return data if isinstance(data, list) else data.get("cases") or data.get("data") or []

def bucket(gold: str, pred: str) -> str | None:
    g, p = _ST.search(gold or ""), _ST.search(pred or "")
    if g and not p:
        return "missing_stereo_prefix"
    if g and p and g.group(0) != p.group(0):
        return "wrong_stereo_letter"
    if (gold or "").count("@") or "/" in (gold or ""):
        return None
    return None

def main() -> None:
    # 用法：先用现有 bench 失败列表；若无，则对 merged 抽样 namer
    print("implement against local fail dump; see plan Task1 notes")

if __name__ == "__main__":
    main()
```

实现时对接项目里已有的 fail 导出（若 `benchmark_parallel` 已写 fails 文件则读该文件；否则对 `data/merged_benchmark.json` 用 `SMILESNNamer` 扫金标含立体前缀的子集）。

- [ ] **Step 3: 输出桶大小，决定 Task2 白名单顺序**

Expected: `missing_stereo_prefix` 按推测母体 kind 再细分；优先桶 ≥8 的 kind。

- [ ] **Step 4: Commit 工具脚本（若有价值）**

```bash
git add tools/_stereo_fail_bucket.py
git commit -m "chore(tools): stereo fail bucketing helper for CIP plan"
```

---

### Task 2: 扩展 `_RS_KINDS`（酮/酯/酰胺/腈/硫醇/醛）

**Files:**
- Modify: `src/namepredict/layer5/stereo_rs.py`（仅 `_RS_KINDS`）
- Modify: `tests/unit/test_rs_stereo.py`
- Test: `tests/unit/test_rs_stereo.py`

**Interfaces:**
- Consumes: 现有 `apply_rs_prefix(numbered, en, zh) -> tuple[str, str]`
- Produces: 更多 `parent.kind` 带链上 R/S

- [ ] **Step 1: 写失败测试（新 kind）**

在 `tests/unit/test_rs_stereo.py` 的 `CASES` 追加：

```python
    # ketone / ester / amide / nitrile / thiol / aldehyde
    (
        "CC[C@H](C)C(=O)C",
        "(3S)-3-methylpentan-2-one",
        "(3S)-3-甲基戊-2-酮",
    ),
    (
        "CC[C@H](C)C(=O)OC",
        "methyl (2S)-2-methylbutanoate",
        "(2S)-2-甲基丁酸甲酯",
    ),
    (
        "CC[C@H](C)C(=O)N",
        "(2S)-2-methylbutanamide",
        "(2S)-2-甲基丁酰胺",
    ),
    (
        "CC[C@H](C)C#N",
        "(2S)-2-methylbutanenitrile",
        "(2S)-2-甲基丁腈",
    ),
    (
        "CC[C@H](C)S",
        "(2S)-butane-2-thiol",
        "(2S)-丁-2-硫醇",
    ),
    (
        "CC[C@H](C)C=O",
        "(2S)-2-methylbutanal",
        "(2S)-2-甲基丁醛",
    ),
```

**注意：** 金标 R/S 必须以 RDKit `_CIPCode` + 本仓库 `normalize_en` 为准；写测试前先用：

```bash
PYTHONPATH=src python -c "from rdkit import Chem; m=Chem.MolFromSmiles('CC[C@H](C)C(=O)C'); Chem.AssignStereochemistry(m,force=True,cleanIt=True); print([(a.GetIdx(),a.GetProp('_CIPCode') if a.HasProp('_CIPCode') else None) for a in m.GetAtoms()])"
```

若字母与上表不符，**改测试期望为 RDKit 字母**，不要硬拧 CIP。

- [ ] **Step 2: 跑测试确认失败**

```bash
PYTHONPATH=src pytest tests/unit/test_rs_stereo.py -v
```

Expected: 新 case 失败（名称无 `(nX)-` 前缀）。

- [ ] **Step 3: 扩展白名单**

```python
_RS_KINDS = frozenset({
    "acid", "alkenoic_acid", "alcohol", "alkenol",
    "diol", "triol", "amine", "diamine",
    "ketone", "dione", "ester", "alkenoate",
    "amide", "nitrile", "alkenenitrile", "alkenal",
    "aldehyde", "thiol",
})
```

- [ ] **Step 4: 跑测试通过**

```bash
PYTHONPATH=src pytest tests/unit/test_rs_stereo.py -v
```

Expected: PASS。若酯名是 `methyl …` 而 R/S 被贴到整串最前，确认 `apply_rs_prefix` 仍正确（当前实现是整名 strip 前缀；若金标要 `methyl (2S)-…`，则本 Task 只保证酸根/后缀母体，酯专项放到 Task 3）。

- [ ] **Step 5: 单元通过后跑 benchmark**

```bash
PYTHONPATH=src python -m benchmarks.benchmark_parallel --data data/merged_benchmark.json --time
```

Expected: dual 不回退；`missing_stereo` 桶下降。

- [ ] **Step 6: Commit + workstate**

```bash
git add src/namepredict/layer5/stereo_rs.py tests/unit/test_rs_stereo.py
git commit -m "feat(namepredict): extend CIP R/S kinds to ketone/ester/amide/nitrile/thiol/aldehyde"
```

`workstate.md` 追加一行（填真实 hash 与 dual）：

```
[#HASH][IUPAC P-92 / P-93] 扩展母体链 R/S 白名单(酮/酯/酰胺/腈/硫醇/醛) [+N tests, dual X%→Y%]
```

---

### Task 3: 酯/酰胺名中 R/S 插入位置（`methyl (2S)-…`）

**Files:**
- Modify: `src/namepredict/layer5/stereo_rs.py`（`_with_rs` 或新增 `_with_rs_ester`）
- Modify: `src/namepredict/layer5/assembler.py`（若需 kind 提示）
- Test: `tests/unit/test_rs_stereo.py`

**Interfaces:**
- Consumes: `parent.kind == "ester" | "alkenoate"`，英文名形如 `methyl 2-methylbutanoate`
- Produces: `methyl (2S)-2-methylbutanoate`（R/S 在烷基名之后、酰基之前）

- [ ] **Step 1: 失败测试**

```python
    (
        "CC[C@H](C)C(=O)OC",
        "methyl (2S)-2-methylbutanoate",
        "(2S)-2-甲基丁酸甲酯",  # 中文以现网金标为准，可只 assert en
    ),
```

- [ ] **Step 2: 实现最小插入逻辑**

约束：函数 ≤10 行；可拆：

```python
def _split_ester_en(en: str) -> tuple[str, str] | None:
    """'methyl 2-methylbutanoate' → ('methyl ', '2-methylbutanoate')."""
    if " " not in en or en.startswith("("):
        return None
    i = en.find(" ")
    return en[: i + 1], en[i + 1 :]

def _with_rs_en(en: str, rs: list[tuple[int, str]], kind: str | None) -> str:
    if kind in ("ester", "alkenoate"):
        sp = _split_ester_en(en)
        if sp:
            head, tail = sp
            return head + _with_rs(tail, rs)
    return _with_rs(en, rs)
```

`apply_rs_prefix` 传入 `kind`：

```python
def apply_rs_prefix(numbered: dict, en: str, zh: str) -> tuple[str, str]:
    rs = _rs_parts(numbered)
    if not rs:
        return en, zh
    kind = (numbered.get("parent") or {}).get("kind")
    return _with_rs_en(en, rs, kind), _with_rs(zh, rs)
```

- [ ] **Step 3: 测试 + benchmark + commit**

```bash
PYTHONPATH=src pytest tests/unit/test_rs_stereo.py -v
PYTHONPATH=src python -m benchmarks.benchmark_parallel --data data/merged_benchmark.json --time
git add src/namepredict/layer5/stereo_rs.py tests/unit/test_rs_stereo.py
git commit -m "feat(namepredict): place R/S after ester alkyl word (methyl (2S)-…)"
```

---

### Task 4: 统一 E/Z 收集模块并挂到 alkene

**Files:**
- Create: `src/namepredict/layer5/stereo_ez.py`
- Modify: `src/namepredict/layer5/unsat_acid.py`（改为 re-export / 调用 stereo_ez）
- Modify: `src/namepredict/layer5/assembler.py`（`_alkene_names` 加 ez）
- Create: `tests/unit/test_ez_stereo.py`

**Interfaces:**
- Consumes: `parent.mol`, `parent.double_bond` 或 `parent.double_bonds`, `parent.chain`
- Produces: `ez_prefix(numbered) -> str` 如 `""` / `"(E)-"` / `"(2E,6Z)-"`

- [ ] **Step 1: 从 unsat_acid 抽出纯函数到 stereo_ez（行为不变）**

```python
# src/namepredict/layer5/stereo_ez.py
"""Double-bond E/Z prefixes (IUPAC P-91.2 / P-93.4)."""
from __future__ import annotations
from rdkit.Chem import BondStereo, Mol

def _stereo_tag(st) -> str:
    if st == BondStereo.STEREOE:
        return "(E)-"
    if st == BondStereo.STEREOZ:
        return "(Z)-"
    return ""

def _bond_stereo(mol: Mol | None, double_bond) -> str:
    if mol is None or not double_bond:
        return ""
    c1, c2 = double_bond
    bond = mol.GetBondBetweenAtoms(int(c1), int(c2))
    return _stereo_tag(bond.GetStereo()) if bond is not None else ""

# ... 迁移 _bond_min_loc, _ez_letter, _ez_bond_part, _ez_parts,
# _ez_multi_prefix, _ez_prefix, _ez_for_alkenol 保持签名
```

`unsat_acid.py`：

```python
from namepredict.layer5.stereo_ez import (
    _ez_for_alkenol, _ez_prefix,  # 或公开名 ez_prefix
)
```

- [ ] **Step 2: 失败测试 — 简单烯烃**

```python
# tests/unit/test_ez_stereo.py
import pytest
from namepredict.constants import normalize_en
from namepredict.namer import SMILESNNamer

CASES = [
    (r"C/C=C/C", "(E)-but-2-ene"),
    (r"C/C=C\C", "(Z)-but-2-ene"),
    (r"CC/C=C/CC", "(E)-hex-3-ene"),
]

@pytest.mark.parametrize("smi,en", CASES)
def test_alkene_ez(smi, en):
    r = SMILESNNamer().name(smi)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
```

- [ ] **Step 3: 改 `_alkene_names` 路径**

在 `assembler.py` 的 alkene 组装处（现有 `_alkene_or_poly` / `_alkene_names`）前缀拼接 `ez`：

```python
def _alkene_names_ez(n: int, numbered: dict) -> tuple[str, str] | None:
    from namepredict.layer5.stereo_ez import ez_for_parent
    base = _alkene_names(
        n, numbered.get("ene_locant"), numbered.get("omit_ene_locant", False),
    )
    if not base:
        return None
    ez = ez_for_parent(numbered)  # 公开包装 _ez_for_alkenol / _ez_prefix
    en, zh = base
    return f"{ez}{en}", f"{ez}{zh}"
```

确保 `parent` 对 `kind=="alkene"` 已带 `double_bond`；若 L2 未填，本 Task 只读已有字段，缺则返回无 ez（并在 workstate 记 follow-up）。若缺字段，最小补丁在 **L2 alkene parent 产出** 填 `double_bond`（仍属 Task 4，单独小函数 ≤10 行）。

- [ ] **Step 4: 测试 + 旧 unsat 测试回归**

```bash
PYTHONPATH=src pytest tests/unit/test_ez_stereo.py tests/unit/test_alkenoic_ez.py tests/unit/test_alkenol.py tests/unit/test_rs_stereo.py -v
```

Expected: 全 PASS。

- [ ] **Step 5: benchmark + commit**

```bash
PYTHONPATH=src python -m benchmarks.benchmark_parallel --data data/merged_benchmark.json --time
git add src/namepredict/layer5/stereo_ez.py src/namepredict/layer5/unsat_acid.py src/namepredict/layer5/assembler.py tests/unit/test_ez_stereo.py
git commit -m "feat(namepredict): unified E/Z helper and alkene stereo prefixes"
```

---

### Task 5: E/Z 挂到 alkenal / alkenoate / alkenenitrile

**Files:**
- Modify: `src/namepredict/layer5/assembler.py`（各 unsat 名函数）
- Modify: `tests/unit/test_ez_stereo.py`
- 可能 Modify: 对应 L2 parent 的 `double_bond` 字段

**Interfaces:**
- Consumes: `ez_for_parent(numbered) -> str`
- Produces: `(E)-…enal` / `methyl (E)-…enoate` 等

- [ ] **Step 1: 失败测试**

```python
CASES += [
    (r"C/C=C/C=O", "(E)-but-2-enal"),
    (r"C/C=C/C#N", "(E)-but-2-enenitrile"),
    (r"C/C=C/C(=O)OC", "methyl (E)-but-2-enoate"),
]
```

- [ ] **Step 2: 在 assembler 对应分支加 ez（与 Task4 同模式）**

- [ ] **Step 3: pytest + benchmark + commit**

```bash
PYTHONPATH=src pytest tests/unit/test_ez_stereo.py tests/unit/test_alkenal.py tests/unit/test_alkenoate.py tests/unit/test_alkenenitrile.py -v
git commit -m "feat(namepredict): E/Z prefixes on alkenal/alkenoate/alkenenitrile"
```

---

### Task 6: 抽出 `stereo_format` + 单测纯字符串合并

**Files:**
- Create: `src/namepredict/layer5/stereo_format.py`
- Modify: `src/namepredict/layer5/stereo_rs.py`（委托 format）
- Create: `tests/unit/test_stereo_format.py`

**Interfaces:**
- Produces:

```python
def parse_stereo_tag(tag: str) -> list[tuple[int | None, str]]: ...
def format_stereo_parts(parts: list[tuple[int | None, str]]) -> str: ...
def merge_stereo_parts(
    old: list[tuple[int | None, str]],
    new: list[tuple[int | None, str]],
) -> list[tuple[int | None, str]]: ...
def apply_stereo_tag(name: str, parts: list[tuple[int | None, str]]) -> str: ...
```

- [ ] **Step 1: 纯函数测试（无 RDKit）**

```python
from namepredict.layer5.stereo_format import (
    parse_stereo_tag, format_stereo_parts, merge_stereo_parts, apply_stereo_tag,
)

def test_merge_ez_then_rs():
    old = parse_stereo_tag("(E)-")
    new = [(8, "R")]
    parts = merge_stereo_parts(old, [(loc, let) for loc, let in new])
    assert format_stereo_parts(parts) == "(E,8R)-"

def test_apply_on_bare():
    assert apply_stereo_tag("butan-2-ol", [(2, "S")]) == "(2S)-butan-2-ol"
```

- [ ] **Step 2: 迁移 stereo_rs 中 strip/parse/format/merge**

- [ ] **Step 3: pytest 全立体相关 + commit**

```bash
PYTHONPATH=src pytest tests/unit/test_stereo_format.py tests/unit/test_rs_stereo.py -v
git commit -m "refactor(namepredict): extract stereo prefix format helpers"
```

---

### Task 7: L1 StereoProvider（RDKit 实现）

**Files:**
- Create: `src/namepredict/layer1/stereo_types.py`
- Create: `src/namepredict/layer1/stereo_provider.py`
- Modify: `src/namepredict/layer5/stereo_rs.py`（经 provider 读码）
- Create: `tests/unit/test_stereo_provider.py`

**Interfaces:**
- Produces:

```python
# stereo_types.py
from dataclasses import dataclass
from typing import Literal

StereoLetter = Literal["R", "S", "E", "Z", "r", "s", "M", "P"]

@dataclass(frozen=True)
class AtomStereo:
    atom_idx: int
    letter: StereoLetter

@dataclass(frozen=True)
class BondStereoPart:
    a_idx: int
    b_idx: int
    letter: StereoLetter

# stereo_provider.py
class StereoProvider(Protocol):
    def atom_codes(self, mol: Mol) -> dict[int, str]: ...
    def bond_codes(self, mol: Mol) -> dict[tuple[int, int], str]: ...

class RdkitStereoProvider:
    def atom_codes(self, mol: Mol) -> dict[int, str]: ...
    def bond_codes(self, mol: Mol) -> dict[tuple[int, int], str]: ...

DEFAULT_PROVIDER: StereoProvider = RdkitStereoProvider()
```

每个方法 ≤10 行：

```python
def _assign(mol: Mol) -> None:
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)

def _atom_letter(atom) -> str | None:
    if not atom.HasProp("_CIPCode"):
        return None
    c = atom.GetProp("_CIPCode")
    return c if c in ("R", "S", "r", "s") else None
```

- [ ] **Step 1: provider 单测（甘油酸、2-丁醇）**

```python
from rdkit import Chem
from namepredict.layer1.stereo_provider import RdkitStereoProvider

def test_glyceric_S():
    mol = Chem.MolFromSmiles("O=C(O)[C@@H](O)CO")
    codes = RdkitStereoProvider().atom_codes(mol)
    assert "S" in codes.values() or "R" in codes.values()
```

- [ ] **Step 2: stereo_rs 改为 `DEFAULT_PROVIDER.atom_codes(mol)`**

- [ ] **Step 3: 全回归 + benchmark + commit**

```bash
PYTHONPATH=src pytest tests/unit/test_stereo_provider.py tests/unit/test_rs_stereo.py tests/unit/test_ez_stereo.py -v
git commit -m "feat(namepredict): L1 RdkitStereoProvider behind stereo_rs"
```

---

### Task 8: 环母体链/环原子 R/S（cycloalcohol / cycloamine / cycloketone）

**Files:**
- Modify: `src/namepredict/layer5/stereo_rs.py`（`_cip_on_chain` 泛化为 `_cip_on_atoms`）
- Modify: L2 环 parent 确保 `chain` 或 `ring_atoms` 有序
- Test: `tests/unit/test_rs_stereo.py`

**Interfaces:**
- Consumes: `parent.chain` 或 `parent.ring_atoms: list[int]`（编号序）
- Produces: 环上位次 R/S

- [ ] **Step 1: 选 2–3 个 dual fail 中「环醇/环胺骨架已对只缺 R/S」的 SMILES 做金标测试**

（具体 SMILES 来自 Task 1 桶；示例）

```python
    (
        "C[C@H]1CCC[C@@H](C1)O",  # 示例：以 namer 当前骨架名为准
        "(1R,3S)-3-methylcyclohexan-1-ol",  # 期望以实际编号为准校准
        None,
    ),
```

- [ ] **Step 2: 实现 `_atoms_for_rs(parent) -> list[int]`**

```python
def _atoms_for_rs(parent: dict) -> list[int]:
    chain = parent.get("chain") or []
    if chain:
        return [int(x) for x in chain]
    ring = parent.get("ring_atoms") or parent.get("atoms") or []
    return [int(x) for x in ring]
```

- [ ] **Step 3: 将 cycloalcohol / cycloamine / cycloketone 加入 `_RS_KINDS`**

- [ ] **Step 4: pytest + benchmark + commit**

```bash
git commit -m "feat(namepredict): R/S on cycloalcohol/amine/ketone ring atoms"
```

---

### Task 9: 取代基手性描述符（侧链立体）

**Files:**
- Modify: `src/namepredict/layer3/` 取代基结构（若需携带 atom map）
- Modify: `src/namepredict/layer5/benzene_names.py` 或 `_build_prefix` 路径
- Create: `tests/unit/test_substituent_stereo.py`

**Interfaces:**
- Consumes: 取代基片段 mol + 连接点；`StereoProvider.atom_codes`
- Produces: 如 `(1-phenylethyl)` 带 `(R)` / `[(1R)-1-phenylethyl]` 形式（严格按 IUPAC P-93 括号层级；与现有 L5 括号风格对齐）

- [ ] **Step 1: 只做「开链烷基侧链一个手性中心」MVP**

金标例（以 benchmark 真实条目校准）：

```python
# 例如 sec-butyl 手性：骨架名已能出 1-methylpropyl 时再贴 (1R)
```

- [ ] **Step 2: 若 L3 无原子映射则最小扩展 `Substituent` 字典字段 `atom_ids: list[int]`**

- [ ] **Step 3: L5 前缀组装读 provider，仅当 codes 非空**

- [ ] **Step 4: 测试 + benchmark（防 regress 字母序）+ commit**

```bash
git commit -m "feat(namepredict): CIP descriptors on simple chiral alkyl substituents"
```

**范围控制：** 不做复杂嵌套 yl 的立体；嵌套留给 Phase D 后。

---

### Task 10: Phase D — CIP digraph 骨架（四面体，可切换）

**Files:**
- Create: `src/namepredict/layer1/cip_digraph.py`
- Create: `src/namepredict/layer1/cip_rank.py`
- Create: `src/namepredict/layer1/cip_provider.py`  # `DigraphStereoProvider`
- Modify: `src/namepredict/layer1/stereo_provider.py`（可选切换）
- Create: `tests/unit/test_cip_digraph.py`

**Interfaces:**
- Produces:

```python
def build_digraph(mol: Mol, center: int) -> object: ...
def rank_neighbors(mol: Mol, center: int) -> list[int]: ...  # 4 neighbor atom idx, high→low
def tetrahedral_letter(mol: Mol, center: int) -> str | None:  # "R"|"S"|None
```

- [ ] **Step 1: 规范用例（IUPAC 教材级，固定 SMILES）**

```python
# tests/unit/test_cip_digraph.py
from rdkit import Chem
from namepredict.layer1.cip_rank import tetrahedral_letter

def _code(smi: str, idx: int) -> str | None:
    mol = Chem.MolFromSmiles(smi)
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    # digraph letter must match RDKit on these easy cases
    from namepredict.layer1.cip_rank import tetrahedral_letter
    return tetrahedral_letter(mol, idx)

def test_easy_match_rdkit():
    smi = "C[C@H](O)Cl"  # 校准 idx
    mol = Chem.MolFromSmiles(smi)
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    for a in mol.GetAtoms():
        if a.HasProp("_CIPCode"):
            assert tetrahedral_letter(mol, a.GetIdx()) == a.GetProp("_CIPCode")
```

- [ ] **Step 2: 实现 digraph（复制原子、原子序数、双键 duplicate）**

拆文件保持 ≤500 行 / 函数 ≤10 行：
- `cip_digraph.py`: 节点扩展
- `cip_rank.py`: BFS 层级比较
- `cip_assign.py`（可选）: 手性体积 → R/S

- [ ] **Step 3: `DigraphStereoProvider.atom_codes` 与 RDKit 对齐测试集**

```python
class DigraphStereoProvider:
    def atom_codes(self, mol: Mol) -> dict[int, str]:
        out = {}
        for a in mol.GetAtoms():
            let = tetrahedral_letter(mol, a.GetIdx())
            if let:
                out[a.GetIdx()] = let
        return out
```

- [ ] **Step 4: 环境开关**

```python
# stereo_provider.py
import os
def get_provider() -> StereoProvider:
    if os.environ.get("NAMEPREDICT_CIP") == "digraph":
        from namepredict.layer1.cip_provider import DigraphStereoProvider
        return DigraphStereoProvider()
    return RdkitStereoProvider()
```

默认仍 RDKit；CI 可加 digraph job。

- [ ] **Step 5: 对齐集 100% + dual 默认路径不回退 + commit**

```bash
NAMEPREDICT_CIP=digraph PYTHONPATH=src pytest tests/unit/test_cip_digraph.py -v
PYTHONPATH=src python -m benchmarks.benchmark_parallel --data data/merged_benchmark.json --time
git commit -m "feat(namepredict): optional CIP digraph provider for tetrahedral R/S"
```

---

### Task 11: digraph 难例 — 同位素 / 多重键 / 环展开

**Files:**
- Modify: `src/namepredict/layer1/cip_digraph.py`, `cip_rank.py`
- Modify: `tests/unit/test_cip_digraph.py`

**Interfaces:** 同 Task 10

- [ ] **Step 1: 用例表**（从 `docs/iupac/cn/第7章_立体化学.md` 摘 5–10 例）

覆盖：C  vs 同位素 H/D；C=C 复制；小环。

- [ ] **Step 2: 逐例 TDD 修 digraph**

- [ ] **Step 3: 与 RDKit 分歧时记录 `tools/cip_rdkit_divergence.md`（只列 SMILES，不改默认 provider）**

- [ ] **Step 4: commit**

```bash
git commit -m "feat(namepredict): CIP digraph hard cases (isotope/unsaturation/rings)"
```

---

### Task 12: Phase E — 轴向手性 M/P（联烯 MVP）

**Files:**
- Create: `src/namepredict/layer1/cip_axial.py`
- Modify: `src/namepredict/layer1/cip_provider.py` 或 `stereo_ez` 旁路
- Create: `tests/unit/test_axial_planar.py`
- Modify: L5 前缀支持 `M`/`P` 字母（`stereo_format` 白名单）

**Interfaces:**
- Produces: `axial_letter(mol, cumulated_atoms) -> "M"|"P"|None`

- [ ] **Step 1: 失败测试**

```python
# 联烯例：以 IUPAC 例题 SMILES 为准
# 期望名称含 (M)- 或 (P)-；仅当 L2 已能给出 penta-2,3-diene 类母体
```

- [ ] **Step 2: 实现联烯两端取代基排序 + 螺旋性**

- [ ] **Step 3: L5 合并字母 `M/P` 到 stereo tag**

- [ ] **Step 4: pytest + benchmark（预期 dual 微涨或持平）+ commit**

```bash
git commit -m "feat(namepredict): allene axial M/P stereodescriptors"
```

**延后：** 联苯阻转、螺环轴性 — 单独立项，不在本 Task 塞入。

---

### Task 13: Phase F — 假不对称 r/s 与相对构型

**Files:**
- Modify: `src/namepredict/layer1/cip_rank.py`（r/s）
- Create: `src/namepredict/layer5/stereo_relative.py`
- Create: `tests/unit/test_relative_stereo.py`

**Interfaces:**
- Produces: 描述符 `r`/`s`；名称前缀 `rel-` / `rac-`（仅当输入无绝对立体或显式标记）

- [ ] **Step 1: r/s 用例（内消旋糖醇类简化结构）**

- [ ] **Step 2: 仅当 digraph 判定假不对称中心时输出小写**

- [ ] **Step 3: `rel-` 仅在 mol 无绝对立体且编号需相对关系时**（输入启发式：无 `@` 但有 cis/trans 环关系 — 可先 skip 自动 rel，只做 r/s）

- [ ] **Step 4: commit**

```bash
git commit -m "feat(namepredict): pseudoasymmetric r/s descriptors"
```

---

### Task 14: Phase G — L4 编号消费立体（最低位次集）

**Files:**
- Modify: `src/namepredict/layer4/` 编号评分
- Test: 专有 case「两端等价时立体描述符集合更低」

**Interfaces:**
- Consumes: 候选编号下的 `list[StereoPart]`
- Produces: 选中使立体位次集最低的编号

- [ ] **Step 1: 找 benchmark 中「编号镜像导致 R/S 位次不同」的失败例**

- [ ] **Step 2: 评分元组增加 stereo locant set（IUPAC 次序：先主 FG，再不饱和，再取代，再立体）**

- [ ] **Step 3: 回归 + benchmark + commit**

```bash
git commit -m "feat(namepredict): L4 numbering tie-break with stereo locant sets"
```

---

### Task 15: 默认切换 digraph + 文档收尾

**Files:**
- Modify: `src/namepredict/layer1/stereo_provider.py` 默认 provider
- Modify: `docs/RULES.md` 或 `workstate.md`「其他」节 ≤20 行说明
- Test: 全 unit + full benchmark

- [ ] **Step 1: 对齐集与 dual：digraph ≥ RDKit 默认 dual**

- [ ] **Step 2: 切换默认；保留 `NAMEPREDICT_CIP=rdkit` 回滚开关**

- [ ] **Step 3: 全量**

```bash
PYTHONPATH=src pytest tests/unit/test_rs_stereo.py tests/unit/test_ez_stereo.py tests/unit/test_stereo_format.py tests/unit/test_stereo_provider.py tests/unit/test_cip_digraph.py -v
PYTHONPATH=src python -m benchmarks.benchmark_parallel --data data/merged_benchmark.json --time
```

- [ ] **Step 4: commit**

```bash
git commit -m "feat(namepredict): default to digraph CIP provider with rdkit fallback"
```

---

## 依赖关系

```
Task1 基线桶
  → Task2 RS 白名单扩展
  → Task3 酯插入位置
  → Task4 E/Z 模块 + alkene
  → Task5 其它 unsat E/Z
  → Task6 format 抽取
  → Task7 Provider 接口
  → Task8 环 R/S
  → Task9 取代基立体
  → Task10 digraph MVP
  → Task11 digraph 难例
  → Task12 轴向 M/P
  → Task13 r/s rel
  → Task14 L4 编号
  → Task15 默认切换
```

可并行（在接口稳定后）：Task6 ∥ Task5；Task12 ∥ Task13（均依赖 Task10）。

---

## 每阶段 dual 验收门槛

| 阶段结束于 | 最低要求 |
|---|---|
| Task 2–5 (Phase A) | dual 相对阶段起点 **≥ +0.0**；争取 `missing_stereo_prefix` 桶明显下降 |
| Task 7 (Phase B) | dual 回退 **≤ 0**（持平可接受） |
| Task 8–9 (Phase C) | dual **≥ 0**；有明确 +N 用例 |
| Task 10–11 (Phase D) | digraph 对齐集通过；默认 RDKit dual 不降 |
| Task 12–13 | 新单测绿；dual 不降 |
| Task 14–15 | 切换后 dual **≥** 切换前 |

---

## 明确不做（YAGNI 边界）

- 不在 L2 写死 SMILES→R/S 特判表。
- 不做完整无机立体 / 配位几何。
- 不做糖完整绝对构型系统名（属天然产物 P-10x，另计划）。
- 不把 CIP digraph 塞进 L5。
- 不在 dual 仍被母体错误主导时强行推进 Task12–15；可暂停本计划插队母体选题。

---

## Self-Review

1. **Spec coverage:** 四面体 R/S、E/Z、假不对称、轴向、相对构型、L4 编号、可替换 provider 均有 Task；完整面性联苯仅声明延后。
2. **Placeholder scan:** 金标 SMILES 中环系/联烯需用 Task1 桶校准——已写明「以 benchmark 校准」，实施时不得留 TBD 提交。
3. **Type consistency:** `StereoLetter`、`StereoProvider.atom_codes`、`format_stereo_parts` 命名在 Task6–7 固定，后续 Task 沿用。
4. **Constraints:** 文件拆分满足 ≤500 行与 ≤10 行函数；层边界 L1 赋值 / L5 格式。

---

## 执行交接

Plan complete and saved to `docs/superpowers/plans/2026-07-14-complete-cip-stereo.md`.

**Two execution options:**

1. **Subagent-Driven (recommended)** — 每 Task 新 subagent，Task 间 review，快迭代  
2. **Inline Execution** — 本会话按 executing-plans 批量推进并设检查点  

**Which approach?**
