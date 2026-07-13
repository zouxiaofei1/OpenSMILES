## 代码审查报告
- **范围**: `src/namepredict/layer2/arene_carbonyl.py`（新）, `src/namepredict/layer2/ring_parent.py`, `src/namepredict/layer2/parent_selector.py`, `src/namepredict/layer4/numbering.py`, `src/namepredict/layer5/assembler.py`, `src/namepredict/layer5/benzene_names.py`, `tests/unit/test_acetophenone.py`（新）, `tests/unit/test_benzaldehyde.py`
- **主规则/意图**: IUPAC P-64.1.1 保留母体 acetophenone/苯乙酮（Ph–CO–CH3；环上≤2 卤/甲基/羟基；连接碳=1）；并将 benzoic/benzaldehyde 逻辑抽至同文件
- **选题同意**: 是
- **选题意见**: workstate 下一步候选明确为「苯乙酮/acetophenone」；与 benzoic/benzaldehyde 同族保留芳羰基母体，粒度是自然规则边界而非过窄 methyl 特判；propiophenone 以 acetyl-methyl 结构规则拒绝；dual 仅 +1 但属结构性母体推进，非金标糊名
- **共识状态**: 一致（往返轮次: 0）
- **Layer**: L2 母体选择（`arene_carbonyl` 合法）+ L4 定向 + L5 保留名组装；无越界
- **structure_lint**: pass（`structure_lint: ok`）
- **pytest**: pass（`test_acetophenone`+`test_benzaldehyde`+`test_benzoic_acid` 31 passed；全量 unit 618 passed，API 收集错误与本轮无关）
- **benchmark**: dual 5.4%(218) → 5.4%(219)；**中英文双重准确: 5.4%**；en=5.4%(219) zh=22.5%(170) fails=3843（未回退）
- **体量**: 函数体超限 0 处；文件超限 0 处（`parent_selector.py` 498、`assembler.py` 492、`numbering.py` 489、`ring_parent.py` 375、`arene_carbonyl.py` 192，均 ≤500）
- **架构**: 无问题。从 `ring_parent` 抽出 arene 羰基保留母体到 L2 新文件合理；L4 复用 `_orient_benzoic`（`ring_attach_idx=1`）；L5 仅词表扩展；L3 靠 `_is_pure_alkyl_c` 自然不把酮碳当烷基侧链，无需 SMILES 特判
- **化学逻辑**: 无问题。`_acetyl_methyl_c` + `_is_methyl_carbon` 保证仅 acetyl（拒 propiophenone/苯丙酮/二苯酮）；`exclude={ket_c, me}` 不把乙酰甲基计入 ≤2 取代；`_arene_fg_conflict` 对 benzoic/benzaldehyde 显式传 `has_ketone`，对 acetophenone 传 `has_acid/has_aldehyde`，冲突面正确；benzoic/benzaldehyde 探针无回归
- **路径合规**: 可写范围内（仅 `src/namepredict/layer2|4|5` 与 `tests/unit`）
- **代码质量**: 合格（与选题解耦）
- **裁决**: PASS
- **必须修改**（FIX 时逐条可执行）:
- **回滚原因**（ROLLBACK 时）:
- **挑战要点**（CHALLENGE 时）:
- **可提交摘要**（仅 PASS 时，供 workstate）:
  [#pending][IUPAC P-64.1.1] 保留母体 acetophenone/苯乙酮（环上≤2卤/甲基/羟基；乙酰连接=1；抽 arene_carbonyl） [+acetophenone tests, dual 5.4%(218)→5.4%(219), fails 3844→3843]

**备注（非阻断）**: `parent["acetyl_methyl_idx"]` 目前仅 L2 写入、下游未读，可留作文档/后续；propiophenone 负例仍落到错误链酮名 `nonan-3-one`，与本轮「非 acetophenone」断言一致，不扩 scope。