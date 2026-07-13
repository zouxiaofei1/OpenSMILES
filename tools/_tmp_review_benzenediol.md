## 代码审查报告
- **范围**: `src/namepredict/layer2/ring_parent.py`, `src/namepredict/layer2/parent_selector.py`, `src/namepredict/layer3/substituent_extractor.py`, `src/namepredict/layer4/numbering.py`, `src/namepredict/layer5/assembler.py`, `src/namepredict/layer5/benzene_names.py`, `tests/unit/test_benzenediol.py`
- **主规则/意图**: 未取代苯二酚系统名 benzene-a,b-diol（1,2/1,3/1,4）；条款宜标 **P-63.1.2**（ol 后缀/乘积 di），非 P-63.1.4（更高特征基时 hydroxy 前缀）
- **选题同意**: 是
- **选题意见**: workstate 已列「苯二酚」；一轮覆盖三异构体为自然边界；L2 母体 + L4 最低位次集 + L5 组装，非 L5 糊名、非 SMILES 特判。窄在「未取代」可接受（烷基/卤代 follow-up）。条款编号宜在 commit/workstate 中改为 P-63.1.2（或 P-63.1.2 + P-14.3.3）
- **共识状态**: 一致（往返轮次: 0）
- **Layer**: L2→L5；无越界。L2 判 `benzenediol` 并写 `oh_c_idxs`；L3 将 kind 纳入 `_PARENT_OH_KINDS` 滤父 FG；L4 环向 + 位次；L5 仅拼 en/zh
- **structure_lint**: pass
- **pytest**: pass（`test_benzenediol` + `test_phenol_aniline` + `test_simple_benzene`：34 passed）；全量 unit 本轮未重跑，相关簇已绿
- **benchmark**: dual 5.0%(205/4062) → 5.1%(206/4062)；**中英文双重准确: 5.1%**；fails 3857→3856；en=5.1%(208)，zh=21.4%(162)（复核命令与报告一致）
- **体量**: 函数体超限 0 处；文件超限 0 处（`parent_selector.py` **499/500**，下一刀必炸；本轮 `_ring_alcohol_parent` 拆法正确）
- **架构**: 无问题。未取代约束用 `_outside_carbons` + `_hetero_allowed`（严于 phenol 的 `_arene_fg_subs_ok`），与「本轮无烷基/卤」一致
- **化学逻辑**: 有问题（见下，可 FIX）
  1. **中文后缀与金标不一致**：实现/单测 `苯-{a},{b}-二醇`；`data/merged_benchmark.json` 金标为 `苯-1,2-二酚` / `苯-1,3-二酚`（tiers-55/98）
  2. **1,2/1,3/1,4 位次**：`_orient_benzenediol` + `_pair_locs_on` 最低对位次集正确；探针 1,2/1,3/1,4 与 1,3 另一 SMILES 均 OK
  3. **phenol / 链 diol 不回归**：`phenol`、`ethane-1,2-diol`、`propane-1,2,3-triol`、卤/甲基酚英文路径正常
  4. **无 SMILES 特判**
- **路径合规**: 可写范围内；未碰 `data/*` / 计分
- **代码质量**: 不合格（中文 stem 与项目金标错位；架构/体量/en 路径合格）
- **裁决**: FIX
- **必须修改**（FIX 时逐条可执行）:
  1. `benzenediol_names` 中文由 `苯-{loc}-二醇` 改为 `苯-{loc}-二酚`（与金标及「酚」族一致；链状 diol 仍用「二醇」）
  2. 同步 `tests/unit/test_benzenediol.py` 期望：`苯-1,2-二酚`、`苯-1,3-二酚`（1,4 无 zh 金标可继续 `None` 或同规则写二酚）
  3. 改后重跑：`pytest tests/unit/test_benzenediol.py tests/unit/test_phenol_aniline.py tests/unit/test_alkanediol.py -q` 与 full bench；预期 dual 再 +2（→208 量级），en 不变，zh +2
  4. workstate/提交摘要条款改为 **P-63.1.2**（勿写 P-63.1.4）
- **回滚原因**（ROLLBACK 时）: （无）
- **挑战要点**（CHALLENGE 时）: （无）
- **可提交摘要**（仅 PASS 时，供 workstate）: （FIX 关闭前勿 commit）

### 为何 dual 仅 +1 而 en +3（zh 匹配？）
| 行 | SMILES | en | zh 金标 | pred zh | dual |
|----|--------|----|---------|---------|------|
| tiers-55 | Oc1ccccc1O | ✓ benzene-1,2-diol | 苯-1,2-**二酚** | 苯-1,2-**二醇** | ✗ |
| tiers-98 | Oc1cc(O)ccc1 | ✓ benzene-1,3-diol | 苯-1,3-**二酚** | 苯-1,3-**二醇** | ✗ |
| chebi-1357 | Oc1ccc(O)cc1 | ✓ benzene-1,4-diol | （eval_zh=false） | 任意 | ✓（仅 en） |

dual 定义：`all(适用字段)`。三例 en 全对 → en +3；仅 1,4 无 zh 考核 → dual +1；1,2/1,3 因 **酚/醇一字之差** zh 失败。zh 总数仍 162，无新增 zh 命中。修 stem 后 dual 应追上 en 增益。