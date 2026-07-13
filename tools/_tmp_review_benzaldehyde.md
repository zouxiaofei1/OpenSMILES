## 代码审查报告
- **范围**: `src/namepredict/layer2/ring_parent.py`, `src/namepredict/layer2/parent_selector.py`, `src/namepredict/layer4/numbering.py`, `src/namepredict/layer5/assembler.py`, `src/namepredict/layer5/benzene_names.py`, `tests/unit/test_benzaldehyde.py`
- **主规则/意图**: IUPAC P-66.6.1 保留母体 benzaldehyde/苯甲醛（CHO 挂苯环；环上 ≤2 卤/甲基/酚羟基；连接碳=1）
- **选题同意**: 是
- **选题意见**: workstate 下一步候选已写「苯甲醛」；与刚落地的 benzoic 同构切片（≤2 取代），属最小可证伪步且有结构推进；dual +3 例、fails −3，非特判气质。不必再拆「只做未取代」或扩到 acetophenone。
- **共识状态**: 一致（往返轮次: 0）
- **Layer**: L2 选母体 + L4 复用 `ring_attach` 定向 + L5 保留名/位次保留；无越界、无「只 L5 糊名」
- **structure_lint**: pass（`python tools/structure_lint.py --root src/namepredict` → ok）
- **pytest**: pass（`test_benzaldehyde`+`test_benzoic_acid`+`test_phenol_aniline`+`test_benzenediol` 41 passed；全量 unit 忽略无关 API 收集错误后 608 passed）
- **benchmark**: dual 5.3% → 5.4%；**中英文双重准确: 5.4%**（en=5.4% 218/4062，zh=22.4% 169/756，fails 3847→3844）
- **体量**: 函数体超限 0 处；文件超限 0 处（`parent_selector.py` 493 行，`ring_parent.py` 492 行，均 ≤500）
- **架构**: 无问题。`_aldehyde_parent` 先 `_try_benzaldehyde_parent` 再链醛；L4 `"benzaldehyde": _orient_benzoic` 走 `ring_attach_idx`；L5 `arene_fg_parent_names` + `_omit_sub_locants` 不省略位次。
- **化学逻辑**: 无问题（审查重点通过）。
  - **环外醛碳 attach**：`_aldehyde_ring_c`/`_fg_ring_c` 要求唯一醛且恰 1 个环碳邻接 → `ring_attach_idx`，L4 定 1。
  - **outside 排除**：`_benzoic_alkyl_ok`/`_benzoic_extra_n` 以 `exclude_c=aldehyde_c_idx` 从 side starts 与 outside carbons 剔除 CHO 碳，避免当甲基侧链。
  - **杂原子**：`_cooh_oxygen_idxs(fg_c)` 对醛碳返回羰基 O，允许环外 CHO 氧；胺等仍被 `_arene_carbonyl_conflict` 挡住。
  - **回归抽检**：`O=Cc1ccccc1`→benzaldehyde；羟基/卤/甲基位次正确；`OC(=O)c1ccccc1` 仍 benzoic；`CC=O`/`CCC=O` 仍链醛；`O=CCc1ccccc1`（苯乙醛）不进 benzaldehyde；`O=Cc1ccc(N)cc1` 不进本母体。
  - 命名债（非阻断）：`_benzoic_*`/`_cooh_oxygen_idxs` 已泛化参数但函数名仍 benzoic 味，可后续改名，不挡本轮。
- **路径合规**: 可写范围内；未触 `data/*` / `benchmarks/benchmark.py`；`test_benzaldehyde.py` 为新增单测（未跟踪属预期）
- **代码质量**: 合格（与选题解耦：代码合格且选题同意）
- **裁决**: PASS
- **必须修改**（FIX 时逐条可执行）:
  （无）
- **回滚原因**（ROLLBACK 时）: （无）
- **挑战要点**（CHALLENGE 时）: （无）
- **可提交摘要**（仅 PASS 时，供 workstate）:
  [#pending][IUPAC P-66.6.1] 保留母体 benzaldehyde/苯甲醛：环上≤2卤/甲基/羟基；CHO 连接碳=1 [+tests, dual 5.3%→5.4% (215→218), fails 3847→3844]