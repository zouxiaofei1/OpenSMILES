## 代码审查报告
- **范围**: `src/namepredict/layer2/parent_selector.py`, `src/namepredict/layer5/assembler.py`, `src/namepredict/layer5/unsat_acid.py`, `tests/unit/test_alkenoic_acid.py`, `tests/unit/test_alkenoic_ez.py`（新）
- **主规则/意图**: IUPAC P-93.4 / P-31.1 — 开链一元烯酸在母体 C=C 有 RDKit 立体时加 `(E)-`/`(Z)-` 前缀
- **选题同意**: 是
- **选题意见**: 同意。`workstate` 下一步候选已写「单烯酸 E/Z」；与既有烯二酸 E/Z 同构扩展，走 RDKit `BondStereo` 系统规则而非分子特判；增益虽小（+3 dual）但补齐主路径缺口，粒度合适
- **共识状态**: 一致（往返轮次: 0）
- **Layer**: L2 给 parent 挂 `mol`+`double_bond`；L5 组装前缀与词干。无越界（L5 不重选母体）
- **structure_lint**: pass
- **pytest**: pass（`test_alkenoic_ez` + `test_alkenoic_acid` + `test_alkenedioic`：24 passed）
- **benchmark**: dual 5.7%(233) → 5.8%(236)；**中英文双重准确: 5.8%**；en=5.8%(237) zh=23.5%(178) fails 3829→3826
- **体量**: 函数体超限 0 处；文件超限 0 处（`parent_selector.py`=500、`assembler.py`=496，贴顶未破）
- **架构**: 无问题。`_ez_prefix` / `_bond_stereo` 复用烯二酸路径；`alkenoic_acid_names` 从 assembler 下沉到 `unsat_acid.py` 与 `alkenedioic_names` 对齐；`_try_unsat_fg` 统一挂 `mol`（顺带影响 alkenal 等 unsat FG，仅元数据，合理）
- **化学逻辑**: 无问题。E/Z 来自 `BondStereo.STEREOE/Z`；无立体不加前缀；烯二酸回归保留；负例用无立体 SMILES 避免与 E/Z 测冲突
- **路径合规**: 可写范围内
- **代码质量**: 合格
- **裁决**: PASS
- **必须修改**（FIX 时逐条可执行）:
  （无）
- **回滚原因**（ROLLBACK 时）: （无）
- **挑战要点**（CHALLENGE 时）: （无）
- **可提交摘要**（仅 PASS 时，供 workstate）:
  [#pending][IUPAC P-93.4 / P-31.1] 一元烯酸 alkenoic acid (E)/(Z) 前缀 [+test_alkenoic_ez, dual 5.7%→5.8%]