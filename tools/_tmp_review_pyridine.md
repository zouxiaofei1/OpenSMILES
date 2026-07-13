## 代码审查报告
- **范围**: `src/namepredict/layer2/pyridine.py`（新）, `src/namepredict/layer2/parent_selector.py`, `src/namepredict/layer4/numbering.py`, `src/namepredict/layer5/assembler.py`, `src/namepredict/layer5/benzene_names.py`, `tests/unit/test_pyridine.py`（新）
- **主规则/意图**: IUPAC P-22.2.1 保留母体 pyridine：未取代 / 单甲基·单卤（N=1 编号）+ pyridinecarboxylic（吡啶-n-羧酸）
- **选题同意**: 是
- **选题意见**: 同意。苯系保留名（benzene/phenol/benzoic…）之后切入首个杂芳保留母体，主规则边界清晰（N 固定为 1、COOH 位次、与 benzoic 平行的 L2/L4/L5 通路）。虽 workstate「下一步候选」偏 arene+nitro / sec-butyl，但本轮是杂芳主路径探针而非错层补丁；实现经 `_benzene_subs_ok` 实际可覆盖多取代，并非「只写 methyl 硬表」。dual +3 / fails −3 属合理首刀增益。
- **共识状态**: 一致（往返轮次: 0）
- **Layer**: L2 母体判定 + L4 环定向 + L5 词干组装；未越界（L5 未重选母体；L2 未拼最终名）
- **structure_lint**: pass（`python tools/structure_lint.py --root src/namepredict` → ok）
- **pytest**: pass（`test_pyridine + test_simple_benzene + test_benzoic_acid` → 36 passed；全 namepredict unit 687 passed，排除无关 API 收集错误）
- **benchmark**: dual 5.8%(237/4062) → 5.9%(240/4062)；**中英文双重准确: 5.9%**；en=6.0%(243) zh=24.1%(182) fails=3822（与实现方报告一致；回退 0）
- **体量**: 函数体超限 0 处；文件超限 0 处（`parent_selector.py` 491、`assembler.py` 491、`numbering.py` 486、`pyridine.py` 100，均 ≤500）
- **架构**: 无问题。`_ring_parent` 插 pyridine、`_acid_parent` 优先 pyridinecarboxylic 再 benzoic，与 FG 优先序一致；编号复用 `_orient_ring_fixed(..., "n_idx")`，羧酸用虚拟 `ring_attach` 破平局。`pyridine` 词条挂在 `benzene_names._ARENE_FG` / 羧酸名同文件，属命名归类债，非层混用。
- **化学逻辑**: 无阻断问题。抽查：2/3/4-methylpyridine、2-chloro、pyridine-2/3/4-carboxylic acid、2-ethyl、多取代 4-chloro-2-methyl 均正确；2- 而非 6-；benzene/toluene/benzoic 不回归。无 `smiles ==` 特判。已知范围外：氨基吡啶被 `_pyridine_fg_block` 挡掉后掉进开链 amine（本轮未声称覆盖），可接受。
- **路径合规**: 可写范围内（`src/namepredict/**`、`tests/unit/**`）；未改 `data/*` / benchmark 计分
- **代码质量**: 合格（与选题解耦）
- **裁决**: PASS
- **必须修改**（FIX 时逐条可执行）:
  - （无）
- **回滚原因**（ROLLBACK 时）:
  - （无）
- **挑战要点**（CHALLENGE 时）:
  - （无）
- **可提交摘要**（仅 PASS 时，供 workstate）:
  [#pending][IUPAC P-22.2.1] 保留 pyridine + 单甲基/单卤（N=1）+ pyridinecarboxylic [+test_pyridine, dual 5.8%→5.9% (237→240), fails 3825→3822]

**非阻断备注**（不挡 commit）：测例可再补 pyridine-4-carboxylic / 多取代回归；`benzene_names.py` 后续可改名或拆 arene 词表；`parent_selector` 有压缩式无关重排，无害。