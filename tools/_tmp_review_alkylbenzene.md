I have enough evidence for the review. Dual is unchanged at 5.8% with +1 correct; the main defect is L2 mono acceptance by raw outside-carbon count, which lets isobutyl/sec-butyl/tert-butylbenzene become success with bare `benzene`.

## 代码审查报告
- **范围**: `src/namepredict/layer2/ring_parent.py`；`src/namepredict/layer3/substituent_extractor.py`；`tests/unit/test_alkylbenzene.py`
- **主规则/意图**: IUPAC P-29.3.1 / P-22.1.3 单取代烷基苯：L2 放宽 mono outside 至 C1–C4；L3 识别 isopropyl 支链前缀（保留名 isopropyl/异丙基）
- **选题同意**: 是
- **选题意见**: 同意。workstate「下一步候选」已写明 C3+烷基苯/异丙苯；粒度是自然规则边界（C1–C4 直链 + 唯一常见支链 isopropyl），非「只做甲基」过窄切片；落在 L2 母体门槛 + L3 取代基识别，层优先级正确。dual +1 / fails −1 虽小，但是结构性推进。
- **共识状态**: 一致（往返轮次: 0）
- **Layer**: L2（`_mono_benzene_ok`）+ L3（`_is_isopropyl` / `_one_alkyl`）；无越界；L5 未糊名
- **structure_lint**: pass（`structure_lint: ok`）
- **pytest**: pass（`test_alkylbenzene` + simple/multi benzene 及相关 monoalkyl 共 80 passed；全量 `tests/unit -k` 被无关 API 收集错误打断，与本 diff 无关）
- **benchmark**: dual **5.8%(236)** → **5.8%(237)**；fails 3826→3825；en=5.9%(238) zh=23.7%(179)；**中英文双重准确: 5.8%**；无 >0.5% 回退
- **体量**: 函数体超限 0 处；文件超限 0 处  
  新函数约：`_heavies_c` 2、`_is_terminal_methyl` 4、`_is_isopropyl` 8、`_make_isopropyl` 9、`_one_alkyl` 5、`_mono_benzene_ok` 3（均 ≤10）
- **架构**: 有问题  
  1. L2 mono 门槛只数 `_outside_carbons ∈ {1,2,3,4}`，**不校验侧链拓扑**；L3 仅能抽 **直链 C1–C4** 或 **isopropyl**。二者不对齐 → 母体选中 benzene 但取代基为空，L5 输出残缺名。  
  2. `_heavies_c` 与既有 `_c_neighbors` 同义重复（次要）。
- **化学逻辑**: 有问题  
  1. **硬缺陷（本轮引入）**：`CC(C)Cc1ccccc1` / `CCC(C)c1ccccc1` / `CC(C)(C)c1ccccc1`（isobutyl / sec-butyl / tert-butylbenzene）在改后 `success=True` 且 en=`benzene`（缺取代基）。改前 outside 仅 (1,2)，这些分子不会进 simple mono benzene。  
  2. isopropyl 判定本身合理：附着碳为纯烷基 C，相对母体链恰 2 个 free 碳且均为端甲基；与 n-propyl 不混淆；`CC(C)C`→`2-methylpropane`、`CC(C)CC`→`2-methylbutane` 未坏。  
  3. multi-methyl 仍走 `_multi_benzene_ok`（每侧链 outside 计数强制 n=1），xylene/trimethyl 回归 OK。  
  4. 开链副作用：`CCCC(C(C)C)CCC`→`4-isopropylheptane`（保留名，一般命名可接受；PIN 应为 propan-2-yl）。非阻断，但属 L3 全局支链能力外溢，宜在报告知悉，不必本轮 PIN 化。  
  5. 无 SMILES `==` 特判；`branched` 字段写入后暂无读取方，无害。
- **路径合规**: 可写范围内（layer2/layer3 + tests/unit）；未碰 `data/*` / benchmark 计分
- **代码质量**: 不合格（与选题解耦：选题可做，L2 门槛过宽引入假绿残缺名）
- **裁决**: **FIX**
- **必须修改**（FIX 时逐条可执行）:
  1. **收紧 L2 mono 门槛，与 L3 可提取拓扑对齐**：`_mono_benzene_ok` 不得仅用 outside 碳数 ∈{1,2,3,4}。在 `n_sub<=1` 且存在 side start 时，只接受「单一起点 + 可命名侧链」——即 outside 碳构成 **一条直链 C1–C4**，或 **isopropyl（中心 C + 两枚端甲基）**；isobutyl / sec-butyl / tert-butyl 等应 **拒绝** simple benzene（`_is_simple_benzene` 为 False），而不是命名成裸 `benzene`。实现上优先复用/下沉与 L3 一致的 walk/isopropyl 结构判定（避免两套化学定义漂移），且保持函数体 ≤10。
  2. **补负例单测**（`tests/unit/test_alkylbenzene.py` 或并列）：至少覆盖  
     - `CC(C)Cc1ccccc1`、`CCC(C)c1ccccc1`、`CC(C)(C)c1ccccc1`：**不得** `success and en=="benzene"`（期望失败或非残缺苯名）；  
     保留现有正例 toluene/ethyl/propyl/butyl/isopropyl 与 `CC(C)C`→2-methylpropane。
  3. **小清理（建议同轮）**：删除 `_heavies_c`，改用 `_c_neighbors`；未使用的 `branched` 键可删或等 L5 真正消费再留。
- **回滚原因**（ROLLBACK 时）: （无）
- **挑战要点**（CHALLENGE 时）: （无）
- **可提交摘要**（仅 PASS 时，供 workstate）: （FIX 未通过，禁止 commit）

**复核重点对照**  
| 关注点 | 结论 |
|--------|------|
| isopropyl 结构判定 | 拓扑正确；端甲基双重检查 OK |
| mono 放宽不破坏 multi-methyl | 是（multi 仍限甲基） |
| toluene/ethyl 回归 | 是 |
| 开链 2-methylpropane | 是 |
| 函数 ≤10 | 是 |
| 无 SMILES 特判 | 是 |
| L2/L3 对齐 | **否 → FIX** |