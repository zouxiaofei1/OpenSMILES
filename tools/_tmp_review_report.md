## 代码审查报告
- **范围**: `src/namepredict/layer1/analyzer.py`, `src/namepredict/layer2/parent_selector.py`, `src/namepredict/layer5/assembler.py`, `src/namepredict/layer5/stems.py`, `tests/unit/test_dialkyl_sulfide.py`（新）, `tests/unit/test_alkanethiol.py`（另：`prompt.txt` 编排改动，非本功能核心）
- **主规则/意图**: IUPAC **P-63.2.1** 开链简单二烷基硫醚（对称 di… sulfide / 二…硫醚；不对称 alphabetical alkyl alkyl sulfide / …基…基硫醚），C1–C4 直链，镜像既有 P-63.2.2 醚路径
- **选题同意**: 是
- **选题意见**: 与 `workstate.md`「下一步候选: 硫醚」一致；粒度是自然规则边界（dialkyl sulfide 一条主规则，非只做 methyl）；增益虽小（+1 dual）但有结构推进（L1 检测 + L2 母体 + L5 组装，并修正原 CCSC→ethane 假绿）。非特判、非只 L5 糊名。
- **共识状态**: 一致（往返轮次: 0）
- **Layer**: L1 检测 / L2 母体 / L5 组装；**无越界**。L4 无位次需求合理。优先级 alcohol > thiol > amine > ether > sulfide 与 BAD 互斥合理。
- **structure_lint**: pass（`python tools/structure_lint.py --root src/namepredict` → ok）
- **pytest**: pass（`tests/unit/test_dialkyl_sulfide.py` + `test_alkanethiol.py` → **17 passed**）
- **benchmark**: dual **4.4%(179/4062) → 4.4%(180/4062)**；fails **3883 → 3882**；**中英文双重准确: 4.4%**（en=4.4% 180/4062，zh=18.4% 139/756）。回退 0，无 ROLLBACK 阈值问题。
- **体量**: 函数体超限 **0** 处；文件超限 **0** 处（L1 448 / L2 476 / assembler 473 / stems 215，均 ≤500；L2/assembler 逼近上限，后续再动宜拆文件）
- **架构**: 大体无问题；醚/硫醚平行实现清晰。提示：`_ether_cs`/`_sulfide_cs` 同构可接受。未发现高层逻辑下沉或「只 L5 补同分」。
- **化学逻辑**: **有问题（须修）**
  1. **`_is_sulfide_sulfur` 过宽**：仅 `Z=16`、无 H、恰好 2 个 C 邻，未要求总配位数/邻居数=2。`CS(C)=O`（DMSO，金标 dimethyl sulfoxide）被标 `has_sulfide` 并命名为 **dimethyl sulfide**（实测）。醚侧有 `not any(_has_double_bonded_o(c)…)`，硫醚未对齐。
  2. **硫酯误检为 sulfide**：`CC(=O)SC` 同时 `has_sulfide=True` 且醛误检（醛抢母体→acetaldehyde，硫酯本就不在本轮范围）；应对齐醚：排除连在羰基 C 上的 S，避免污染 FG 图。
  3. 支化/芳基/多硫醚/二硫醚正确未吃进（`_arm_ok`、单 S、2C 约束）— 与声明的简单开链范围一致，可接受。
  4. 中文对称「二甲硫醚」与不对称「乙基甲基硫醚」（按 EN 字母序）与测例/金标 CCSC 一致，未见位次格式问题。
- **路径合规**: 生产改动在 `src/namepredict/` + `tests/unit/` 可写范围；未碰 `data/*` / `benchmarks/benchmark.py`。`prompt.txt` 非本规则必需，但不构成硬红旗。工作区另有无关未跟踪噪声（如 `layer4/polyene.py`、tools 临时文件）— 勿一并 commit。
- **代码质量**: **不合格（可局部改正）** — 选题与主路径正确，但 L1 检测缺 degree/羰基约束导致亚砜假阳性，不得按「已知缺口」直接 PASS。
- **裁决**: **FIX**
- **必须修改**（FIX 时逐条可执行）:
  1. 收紧 `src/namepredict/layer1/analyzer.py` 中 `_is_sulfide_sulfur`：在现有条件上增加 **仅二配位硫**（如 `atom.GetTotalDegree() == 2` 或显式邻居全为 C 且 `len(neighbors)==2`），排除亚砜/砜等。
  2. 对齐醚逻辑：对两臂碳增加 **`not any(_has_double_bonded_o(c) for c in cs)`**（或等价），避免硫酯 `CC(=O)SC` 进入 `sulfides`。
  3. 在 `tests/unit/test_dialkyl_sulfide.py` 增加负例至少：`CS(C)=O` **不得**为 dimethyl sulfide（可断言 en 不等于该串，或成功时不含 sulfide 保留名）；建议再加 `CC(=O)SC` 不得为 * methyl sulfide 类名。
  4. 改完后重跑：`pytest tests/unit/test_dialkyl_sulfide.py tests/unit/test_alkanethiol.py -q` 与 `python -m benchmarks.benchmark --data data/merged_benchmark.json`；确认 dual 不回退 >0.5%（预期 DMSO 仍 dual 失败但化学输出不再谎称 sulfide；CCSC 等正例保持）。
  5. commit 时 **不要** 混入 `prompt.txt` 以外的无关未跟踪文件；本轮功能 diff 仅限上述 namepredict + 单测（`prompt.txt` 若属编排另议，勿与硫醚声称绑死）。
- **回滚原因**（ROLLBACK 时）: （无）
- **挑战要点**（CHALLENGE 时）: （无）
- **可提交摘要**（仅 PASS 时，供 workstate）: （FIX 通过前不可提交）