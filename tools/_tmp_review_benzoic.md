## 代码审查报告
- **范围**: `src/namepredict/layer2/ring_parent.py`, `src/namepredict/layer2/parent_selector.py`, `src/namepredict/layer4/numbering.py`, `src/namepredict/layer5/assembler.py`, `src/namepredict/layer5/benzene_names.py`, `tests/unit/test_benzoic_acid.py`
- **主规则/意图**: IUPAC P-65.1.1.1 保留母体 benzoic acid/苯甲酸；环上≤2 卤/甲基/酚羟基；COOH 连接碳=1
- **选题同意**: 是
- **选题意见**: 同意。自然规则边界清晰（单羧基直连苯环 + 允许 ≤2 简单环取代），非「只做一分子」特判；dual +7（5.1%→5.3%）有结构推进；hydroxy 走酸优先+前缀而非 phenol 母体，层优先级正确。workstate「下一步候选」未列本条，但本轮用户/主 Agent 已指定 P-65.1.1.1，可接受。非过窄切片（未拆 methyl/halo/hydroxy 分轮）。
- **共识状态**: 一致（往返轮次: 0）
- **Layer**: L2 选母体 + L4 定向 + L5 词干/位次保留；无越界。L3 未改，因既有 hydroxy/halo/alkyl 提取对 `kind=benzoic` 已适用
- **structure_lint**: pass
- **pytest**: pass（`test_benzoic_acid` + `phenol_aniline` + `benzenediol` + benzene/酸回归相关 78 passed）
- **benchmark**: dual 5.1% → 5.3%；**中英文双重准确: 5.3%**（en=5.3% 215/4062，zh=22.2% 168/756，fails 3854→3847）
- **体量**: 函数体超限 0 处；文件超限 0 处（`parent_selector.py` 恰好 500 行，贴顶可接受）
- **架构**: 无问题。`_try_benzoic_parent` 挂在 `_acid_parent`（酸优先于醇）；L4 抽出 `_orient_ring_fixed` 复用 cycloalcohol/ketone/amine/benzoic；L5 仅词干 + 保留取代基位次，无 L5 糊同分
- **化学逻辑**: 无阻断问题
  - **COOH 不在环上 attach 定向**：`chain`=苯环；`cooh_c_idx`=羧基碳；`ring_attach_idx`=`_carboxyl_ring_c`（羧基碳唯一环碳邻点）；L4 固定 `ring_attach_idx` 为 1 — 正确
  - **hydroxy 前缀**：COOH 优先级使母体走 benzoic；`benzoic ∉ _PARENT_OH_KINDS` → L3 抽 hydroxy；烟测 2/3/4-hydroxybenzoic acid 正确
  - **phenol/链酸不回归**：phenol、acetic、propanedioic、benzenediol 仍对
  - **非特判**：结构谓词（单羧基、直连环、甲基-only 侧链、杂原子允许集、extra_n≤2），无 SMILES 映射
  - 范围外预存：`OC(=O)Cc1ccccc1`（苯乙酸）仍误为链酸 — 正确未进 benzoic，不属本轮
- **路径合规**: 可写范围内（`src/namepredict/layer2|4|5` + `tests/unit`）；未碰 `data/*` / 计分
- **代码质量**: 合格（与选题解耦亦合格）
- **裁决**: PASS
- **必须修改**:
  （无）
- **回滚原因**: （无）
- **挑战要点**: （无）
- **可提交摘要**（仅 PASS 时，供 workstate）:
  [#pending][IUPAC P-65.1.1.1] 保留名 benzoic acid/苯甲酸（环上≤2 卤/甲基/羟基；COOH 连接碳=1） [+tests, dual 5.1%→5.3% (208→215), fails 3854→3847]