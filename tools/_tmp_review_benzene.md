## 代码审查报告
- **范围**: `src/namepredict/layer2/ring_parent.py`, `src/namepredict/layer2/parent_selector.py`, `src/namepredict/layer4/numbering.py`, `src/namepredict/layer5/assembler.py`, `tests/unit/test_simple_benzene.py`
- **主规则/意图**: IUPAC P-22.1.3 简单苯（未取代 + 单卤 + 单 C1–C2 正烷基；甲基保留 toluene/甲苯）
- **选题同意**: 是
- **选题意见**: 同意。`workstate.md` 下一步候选含「芳环单取代」；本轮打开芳环母体路径，粒度（未取代/单卤/C1–C2 + toluene 保留名）是合理最小可证伪步，非过窄刷 commit；层优先级正确（L2 母体 → L4 环定向 → L5 组装/保留名）；dual +0.2%（+7）有结构收益。C3+ 烷基苯可作 follow-up，不构成反对。
- **共识状态**: 一致（往返轮次: 0）
- **Layer**: L2 / L4 / L5；无越界（L2 选 benzene 母体，L4 复用环旋转，L5 词干与省略位次/甲苯）
- **structure_lint**: pass（`python tools/structure_lint.py --root src/namepredict` → ok）
- **pytest**: pass（`test_simple_benzene.py` 12 passed；`tests/unit` 忽略无关 API 收集错误后 551 passed）
- **benchmark**: dual 4.4%(180/4062) → 4.6%(187/4062)；**中英文双重准确: 4.6%**；en=4.6% zh=19.3% fails=3875（本审查复跑一致）
- **体量**: 函数体超限 0 处；文件超限 0 处（`ring_parent.py` 267、`parent_selector.py` 486、`numbering.py` 468、`assembler.py` 488；新增/改动函数体均 ≤10）
- **架构**: 无问题。`_ring_parent` 收束环烃母体；`benzene` 编号复用 `_orient_cycloalkane`；甲苯前缀抑制在 L5，未在 L5 重选母体。
- **化学逻辑**: 有问题  
  1. **不饱和侧链被当成乙基**：`C=Cc1ccccc1`（styrene）、`C#Cc1ccccc1` 成功输出 `ethylbenzene`/`乙基苯`。根因是 `_outside_ok`→`_pure_alkyl_outside` 与 L3 `_is_pure_alkyl_c` 只看 C/H 原子序数、不校验键级；L2 把其收进 simple benzene，L3 再抽成 `n_carbons=2` 烷基。这与「mono C1–C2 **n-烷基**苯」声明不符，属本轮新暴露的错误成功路径（环己烷上乙烯基→`ethylcyclohexane` 为同源既有洞，本轮芳环路径扩大了危害面）。  
  2. 测试负例不足：仅有饱和环己烷/乙醇，未钉死乙烯基/乙炔基苯不得为 ethylbenzene。  
  3. 范围外说明（不阻断本轮主路径）：`CCCc1ccccc1` 等未进 simple benzene 而落到错误链烷名，属既有「芳环未支持→乱选链」问题，非本轮宣称成功集。
- **路径合规**: 可写范围内（无 `data/*`、无 benchmark 计分改动）
- **代码质量**: 不合格（选题同意；结构/体量/门禁绿，但化学假阳性须修）
- **裁决**: FIX
- **必须修改**（FIX 时逐条可执行）:
  1. 在 simple-benzene 准入（优先修共享「环外侧纯烷基」判定，或等价的 benzene 局部检查）中 **拒绝含 C=C/C≡C 的侧链**：`C=Cc1ccccc1`、`C#Cc1ccccc1` 不得再命名为 `ethylbenzene`/`乙基苯`（本轮可 `success=False` 或走未支持，不得假绿）。
  2. 单元测试增加负例并双语严格断言：至少 `C=Cc1ccccc1`、`C#Cc1ccccc1` 的 `normalize_en` **不等于** `ethylbenzene`（若仍 success，也不得是该错误名）；保持既有正例与环己烷回归。
  3. 改完后重跑：`structure_lint`、`pytest tests/unit/test_simple_benzene.py`、以及 dual bench；确认 dual 相对改前 4.4% 回退 ≤0.5%（预期仍应 ≥ 本轮 4.6% 或仅小幅波动）。
- **回滚原因**（ROLLBACK 时）: （无）
- **挑战要点**（CHALLENGE 时）: （无）
- **可提交摘要**（仅 PASS 时，供 workstate）: （无；FIX 完成并复检通过前禁止 commit）