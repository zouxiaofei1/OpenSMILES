## 代码审查报告
- **范围**: `src/namepredict/layer2/ring_parent.py`, `tests/unit/test_aminophenol.py`
- **主规则/意图**: IUPAC P-63.1.4 / P-62.5：氨基酚以 phenol 为母体、amino 为前缀（OH 优先于 NH2）；L2 放行环上伯胺并计入 arene 取代预算
- **选题同意**: 是
- **选题意见**: 同意。与 `workstate.md`「下一步候选: 氨基酚」一致；在既有 phenol 母体 + L3 amino 前缀 + L4 OH=1 编号上只放开 L2 门闩，是自然规则边界而非 SMILES 特判。dual +2（5.5%→5.6%）虽小，但与近期芳环切片量级一致，且同一实现覆盖 o/m/p，并自然容纳 diaminophenol / amino+halo/methyl/nitro（仍受 cap）。非「只做某一异构体」的过窄刷 commit。
- **共识状态**: 一致（往返轮次: 0）
- **Layer**: L2（`ring_parent` 判定 phenol 可带环上伯胺）；L3/L4/L5 沿用既有 amino/phenol 路径，无越界、无 L5 糊名
- **structure_lint**: pass（`python tools/structure_lint.py --root src/namepredict` → ok）
- **pytest**: pass（`test_aminophenol` + `test_phenol_aniline` + `test_nitrobenzene`：29 passed）
- **benchmark**: dual 5.5%(224) → 5.6%(226)；fails 3838→3836；en=5.6%(226) zh=23.0%(174)；**中英文双重准确: 5.6%**（本轮复跑确认，与实现方一致）
- **体量**: 函数体超限 0 处；文件超限 0 处（新增 `_ring_primary_amines` / `_phenol_amines_ok` / `_phenol_allowed` / `_phenol_subs_ok` 与改写后 `_is_simple_phenol` 均 ≤10）
- **架构**: 无问题。仅 L2 门控；`select_parent` 仍 alcohol→amine（酚优先于苯胺）；`_is_simple_aniline` 仍拒 `hydroxyls`；L3 `_extract_aminos` 在 parent=phenol 时已提取 amino
- **化学逻辑**: 无问题（要点已核）:
  - 仅环上伯胺：`degree==1` 且 `c_idx in ring_set`；存在非环/非伯胺则 `_phenol_amines_ok`→None
  - amino 计入预算：`_arene_fg_subs_ok(..., n_amino)`，`halo+alkyl+nitro+amino ≤ 2`
  - N 原子进 `allowed`（与 nitro 模式一致），避免 outside hetero 误杀
  - 叔胺酚 `Oc1ccccc1N(C)C`→失败（不误收）
  - phenol / aniline / 4-nitrophenol 不回归
- **路径合规**: 可写范围内（仅 `src/namepredict/layer2/` + `tests/unit/`）
- **代码质量**: 合格（与选题解耦：代码合格）
- **裁决**: PASS
- **必须修改**（FIX 时逐条可执行）:
  （无）
- **回滚原因**（ROLLBACK 时）: （无）
- **挑战要点**（CHALLENGE 时）: （无）
- **可提交摘要**（仅 PASS 时，供 workstate）:
  [#pending][IUPAC P-63.1.4 / P-62.5] 氨基酚：phenol 母体 + 环上伯胺 amino 前缀（OH 优先） [+aminophenol tests, dual 5.5%→5.6% (224→226)]

**说明（非阻断）**: 测例 docstring 写「exactly one … NH2」，实现实际允许 0–N 个环上伯胺（受 cap）；o/m 未强制 zh、也未断言 diaminophenol 等已能出的名——属覆盖偏窄，不构成 FIX/红旗。