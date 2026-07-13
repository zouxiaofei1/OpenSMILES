## 代码审查报告
- **范围**: `src/namepredict/layer2/ring_parent.py`, `src/namepredict/layer5/assembler.py`, `src/namepredict/layer5/benzene_names.py`（新）, `tests/unit/test_multi_benzene.py`（新）
- **主规则/意图**: IUPAC P-14.3.4 / P-22.1.3 简单多取代苯（环上 2–3 个卤素与/或甲基；英文 xylene 保留名 + 位次；中文系统「二甲基苯」）
- **选题同意**: 是
- **选题意见**: 与 workstate「下一步候选: 多取代苯(二卤/二烷基位次)」及本轮明确范围（卤/甲基）一致；mono 已支持 C1–C2，multi 先收紧为甲基是合理最小可证伪步（xylene 保留名、侧链总碳=起点数的结构门控）；非「只改 3 个 SMILES」特判选题。follow-up 可做 multi-ethyl/C3+。
- **共识状态**: 一致（往返轮次: 0）
- **Layer**: L2（母体门控）+ L5（保留名/前缀组装）；L4 复用既有环旋转最低位次；无越界
- **structure_lint**: pass（`structure_lint: ok`）
- **pytest**: pass（`test_multi_benzene`+`test_simple_benzene` 28；全 unit 忽略无关 API 收集错误后 567 passed）
- **benchmark**: dual 4.6%(187/4062) → 4.8%(196/4062)；fails 3875→3866；**中英文双重准确: 4.8%**（en=4.8% zh=20.5%）
- **体量**: 函数体超限 0 处；文件超限 0 处（`ring_parent.py` 291、`assembler.py` 482、`benzene_names.py` 36，均 ≤500；新/改函数体均 ≤10）
- **架构**: 无问题。甲苯/二甲苯逻辑从 assembler 下沉到同层 `benzene_names.py`；L2 只做「是否简单苯母体」门控，不拼名；无 L5 糊同分、无跨层 import
- **化学逻辑**: 无阻断问题
  - **xylene 特判**: 按 `substituents` 的 `kind==alkyl` 且 `n_carbons==1` 且恰好 2 个，**非 SMILES 表**；中文走系统「{locs}-二甲基」+ 母体「苯」
  - **多取代假阳性**: multi 要求 `len(outside)==len(starts)` 且全为甲基；`CCc1ccc(CC)cc1` 等 **不会** 出 diethylbenzene（实测回落到既有链烷错误路径，非本轮引入的苯假阳性）
  - **单取代回归**: toluene / ethylbenzene / monohalo / benzene 仍绿；mono 仍允许 outside∈{1,2}（乙基）
  - 混合前缀字母序（如 1-bromo-4-chlorobenzene）依赖既有 L4/L5 分组，探针正确
- **路径合规**: 可写范围内（仅 layer2/layer5 + tests/unit）；未碰 `data/*` / benchmark 计分
- **代码质量**: 合格（与选题解耦）
- **裁决**: PASS
- **必须修改**（FIX 时逐条可执行）:
  - （无）
- **回滚原因**（ROLLBACK 时）: （无；dual 未回退）
- **挑战要点**（CHALLENGE 时）: （无）
- **可提交摘要**（仅 PASS 时，供 workstate）:
  [#pending][IUPAC P-14.3.4 / P-22.1.3] 简单多取代苯(2–3卤/甲基；en xylene保留；zh系统二甲基苯) [+tests multi_benzene, dual 4.6%→4.8%(187→196), fails 3875→3866]

### 审查备注（不阻断）
1. 建议后续补 **负例**：multi-ethyl / ethyl+methyl 不得命名为 `*benzene`（当前仅靠 L2 门控，无直接断言）。
2. 单测仅覆盖 `1,4-xylene`；`1,2-`/`1,3-` 探针正确，可补测加固。
3. `_benzene_alkyl_ns` 用 outside 计数推断全甲基（返回 `[1]*n`），对当前范围正确；扩展 multi-ethyl 时应改为真实链长而非硬编码 1。
4. 工作区另有无关改动/未跟踪（`prompt.txt`、`layer4/polyene.py` 等）——**本轮 commit 应只纳入上述 4 个路径**，勿捎带。