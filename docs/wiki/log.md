# Repo Wiki 操作日志

> Append-only — 记录所有 wiki 生成/更新/修改操作。

| 时间 | 操作 | 涉及页面 | 触发者 |
|------|------|---------|--------|
| 2026-07-20 | init | `index.md`, `wiki_plan.yaml`, 目录结构 | chem-dev |
| 2026-07-20 | generate | 6 architecture + 7 supplementary pages | chem-dev |
| 2026-08-06 | update | `architecture/layer2-parent-selector.md`, `architecture/overview.md` — P-44 规则驱动主链管线、纯烃表达、苯环保留名(benzoate)、fg 注册层删除 | chem-dev |
| 2026-08-06 | update | layer2 大重构后 wiki 同步：死代码删除 + 8 组合并（`principal_selection+principal_registry→principal.py`、`aliph_fg+parent_selector_common→fg_helpers.py`、`cyclo_ene_fg+cyclo_poly_fg→cyclo_fg.py`、`spiro_parent+bridged_parent→polycyclic_parent.py`、`benzothiazole+benzoxazole→benzazole.py`、`benzofuran+benzothiophene→fused56_mono.py`、`anthracene→anthraquinone`、`ortho_benzoquinone→benzoquinone`）；注册表（fg/ring/unsat producers）整体删除，骨架识别改由 `scaffold/ring_core`（@_register）+ `scaffold/specs`；涉及 `architecture/layer2-parent-selector.md`, `architecture/overview.md`, `concepts/functional-group-priority.md`, `concepts/atom-ownership.md`, `guides/adding-new-ring-system.md`, `guides/adding-new-functional-group.md`, `reference/core-data-contracts.md`, `index.md` | chem-dev |
| 2026-08-07 | update | Layer3 取代基命名引擎重构后 wiki 同步：侧链拓扑事实层（`side_facts`/`side_alkyl`/`side_alkoxy`/`aryl_sub`/`aryl_depth2`/`leaves/`）自 layer2/tools 迁入 layer3，L2 仅反向借用（carboxyalkyl_arms/_walk_linear/aryl_sub）；`tools/anchored_table.py` 锚定 canonical-SMILES 查表取代 `_BRANCH_CHECKS` 模板 + `_try_*` 一串尝试器；删除 `ring_namer.py`(recursive_ph_name)、`heteroaryl_sub.py`、`cycloalkyl_names.py`、`n_side_extract.py`/`n_block_extract.py`；`substituent_extractor` 445→211、`substituent_namer` 334→128、`make_match` 改显式参数满足 layer3 契约；yl 转换迁 `layer5/free_to_yl.py`，新增 `submol_build.py`(build_anchor_submol)；重写 `architecture/layer3-substituents.md`，修正 `architecture/layer2-parent-selector.md`/`architecture/overview.md`/`guides/*` 中已删模块引用，更新 `index.md` 统计 | chem-dev |
