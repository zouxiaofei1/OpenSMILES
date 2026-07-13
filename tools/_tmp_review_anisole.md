## 代码审查报告
- **范围**: `src/namepredict/layer2/ring_parent.py`, `src/namepredict/layer3/substituent_extractor.py`, `src/namepredict/layer5/benzene_names.py`, `tests/unit/test_anisole.py`
- **主规则/意图**: IUPAC P-63.2.2 芳环烷氧基前缀（methoxy/ethoxy C1–C2）+ 未取代甲氧基苯保留名 anisole/苯甲醚
- **选题同意**: 是
- **选题意见**: 同意。workstate 已列「烷氧基苯/anisole」为下一步候选；本轮以保留 anisole 为探针、同步 C1–C2 直链烷氧基，属最小可证伪步且未把 C3+ 硬塞进同一轮；dual +4 / fails −4 有结构推进，非 SMILES 特判选题
- **共识状态**: 一致（往返轮次: 0）
- **Layer**: L2（环母体放行 alkoxy 并 exclude 外侧 C/O）+ L3（提取 alkoxy 前缀）+ L5（anisole 保留名与前缀省略）；无越界；未只在 L5 糊名
- **structure_lint**: pass（`python tools/structure_lint.py --root src/namepredict`）
- **pytest**: pass（相关 49；全量 namepredict 单测 646 passed，忽略无关 fastapi 收集错误 `test_api_name`/`test_settings_api`）
- **benchmark**: dual 5.6% (226/4062) → 5.7% (230/4062)；fails 3836→3832；en=5.7% (231) zh=23.3% (176)；**中英文双重准确: 5.7%**
- **体量**: 函数体超限 0 处；文件超限 0 处（`ring_parent.py` 485 行，逼近 500 但未超）
- **架构**: 无问题。L2 仅做母体可识别性与 `exclude` 防外侧碳/杂原子双计；L3 产出 `kind=alkoxy`+en/zh；L5 仅保留名分支。开链 `ether` 母体仍由 L2 ether 路径独占，芳环侧不吸收为 dialkyl ether。L2/L3 烷氧行走逻辑有镜像重复（类似既有 nitro），属既有模式，不构成层混用
- **化学逻辑**: 无问题（在宣称范围内）。验证：`COc1ccccc1`→anisole；`CCOc1ccccc1`→ethoxybenzene；`COc1ccc(O)cc1`→4-methoxyphenol；`COc1ccccc1C` 甲基与甲氧基分计无双计；`CCOCC`/`CCOC` 开链醚不回归；C3+ propoxy 不进 simple arene（落回失败路径，属本轮 C1–C2 边界）。未见 SMILES `==` 特判/居所捷径
- **路径合规**: 可写范围内（仅 layer2/3/5 + tests/unit；未碰 `data/*` / benchmark 计分）
- **代码质量**: 合格（与选题解耦）。次要：`ring_parent._ring_alkoxy_n` 已定义未引用（死代码，不阻断）；`ring_parent` 行数偏高，后续 arene 扩展宜再拆文件
- **裁决**: PASS
- **必须修改**（FIX 时逐条可执行）:
  1. （无）
- **回滚原因**（ROLLBACK 时）: （无）
- **挑战要点**（CHALLENGE 时）: （无）
- **可提交摘要**（仅 PASS 时，供 workstate）:
  [#pending][IUPAC P-63.2.2] 芳环烷氧基 methoxy/ethoxy + anisole/苯甲醚 [+test_anisole, dual 5.6%→5.7% (226→230), fails 3836→3832]