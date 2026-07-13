补充确认：无 SMILES 特判；`structure_lint` 绿；相关单测 64 全绿、全量 namepredict 单测 579 绿；bench dual 4.8%→5.0%（+9，fails 3866→3857）。

## 代码审查报告
- **范围**: `src/namepredict/layer2/ring_parent.py`, `src/namepredict/layer2/parent_selector.py`, `src/namepredict/layer3/substituent_extractor.py`, `src/namepredict/layer4/numbering.py`, `src/namepredict/layer5/assembler.py`, `src/namepredict/layer5/benzene_names.py`, `tests/unit/test_phenol_aniline.py`
- **主规则/意图**: IUPAC P-63.1.4 / P-62.2.1.1.1 保留母体 phenol / aniline；OH/NH2 定位为 1；允许 0–2 个环上卤/甲基
- **选题同意**: 是
- **选题意见**: 同意。`workstate` 下一步候选含「酚/苯胺」；酚+胺同一套 arene FG 骨架（`_arene_fg_subs_ok`）属自然单规则切片，非过窄 methyl 刷 commit；dual +0.2%（196→205）、fails −9，有实质推进
- **共识状态**: 一致（往返轮次: 0）
- **Layer**: L2 选 phenol/aniline 母体；L3 父 FG 过滤；L4 复用环醇/环胺定向；L5 保留名组装。无越界、非纯 L5 糊名
- **structure_lint**: pass（ok）
- **pytest**: pass（`test_phenol_aniline`+simple/multi benzene+mono cycloalcohol/amine：64；忽略无关 API 收集错误后 unit 579 passed）
- **benchmark**: dual 4.8%(196/4062) → 5.0%(205/4062)；**中英文双重准确: 5.0%**；en=5.0% zh=21.4% fails=3857（与改前 fails=3866 一致改进；回退未发生）
- **体量**: 函数体超限 0 处；文件超限 0 处。软提醒：`parent_selector.py` 490 行、`assembler.py` 488 行逼近 500，后续再塞规则前宜先拆
- **架构**: 无问题。母体 kind 贯通 L2→L5；编号复用 `_orient_cycloalcohol` / `_orient_cycloamine`；名称落在 `benzene_names.arene_fg_parent_names`，职责清晰
- **化学逻辑**: 无阻断问题
  - **仲胺 c_idx KeyError**：已修。`_mono_amine_on_ring` 用 `am.get("c_idx")`；L1 对 deg≥2 只给 `c_idxs`，与 aniline 路径兼容，避免 N-烷基苯胺/环上仲胺探测时崩溃
  - **误伤 benzene/cyclohexanol**：负例与探针均正确（benzene/chlorobenzene/toluene/cyclohexanol/cyclohexanamine/ethanol 不回归）
  - **氨基酚边界**：`_is_simple_phenol` 拒 `amines`、`_is_simple_aniline` 拒 `hydroxyls` → 双侧 retained 不命中，结果 `success=False`（安全失败，非错名）。IUPAC 更优是 phenol 母体 + amino 前缀，属明确 follow-up，本轮可不扩
  - **N-甲基苯胺**：`degree!=1` 正确排除 aniline；落到错误 `hexane` 为既有开链/芳香仲胺缺口，非本轮引入
  - 乙基酚、邻苯二酚等超范围拒绝；二氯酚在 ≤2 额外取代内正确
  - 无 SMILES `==` / 全名映射特判
- **路径合规**: 可写范围内（`src/namepredict/layer2–5` + `tests/unit`）；未碰 `data/*` / bench 计分
- **代码质量**: 合格（与选题解耦）
- **裁决**: **PASS**
- **必须修改**（FIX 时逐条可执行）:
  - （无）
- **回滚原因**（ROLLBACK 时）:
  - （无）
- **挑战要点**（CHALLENGE 时）:
  - （无）
- **可提交摘要**（仅 PASS 时，供 workstate）:
  [#pending][IUPAC P-63.1.4 / P-62.2.1.1.1] 保留母体 phenol/aniline（环上≤2 卤/甲基，FG=1；仲胺 c_idx 防 KeyError） [+phenol_aniline tests, dual 4.8%→5.0% (196→205), fails 3866→3857]

**提交注意（非否决）**: 工作区另有无关改动/未跟踪（如 `prompt.txt`、`layer4/polyene.py`、其它 unit 等）。本轮 commit 应只纳入上列 7 个审查文件；氨基酚系统命名与 N-烷基苯胺勿顺手塞入。