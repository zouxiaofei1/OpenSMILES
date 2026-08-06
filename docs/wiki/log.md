# Repo Wiki 操作日志

> Append-only — 记录所有 wiki 生成/更新/修改操作。

| 时间 | 操作 | 涉及页面 | 触发者 |
|------|------|---------|--------|
| 2026-07-20 | init | `index.md`, `wiki_plan.yaml`, 目录结构 | chem-dev |
| 2026-07-20 | generate | 6 architecture + 7 supplementary pages | chem-dev |
| 2026-08-06 | update | `architecture/layer2-parent-selector.md`, `architecture/overview.md` — P-44 规则驱动主链管线、纯烃表达、苯环保留名(benzoate)、fg 注册层删除 | chem-dev |
| 2026-08-06 | update | layer2 大重构后 wiki 同步：死代码删除 + 8 组合并（`principal_selection+principal_registry→principal.py`、`aliph_fg+parent_selector_common→fg_helpers.py`、`cyclo_ene_fg+cyclo_poly_fg→cyclo_fg.py`、`spiro_parent+bridged_parent→polycyclic_parent.py`、`benzothiazole+benzoxazole→benzazole.py`、`benzofuran+benzothiophene→fused56_mono.py`、`anthracene→anthraquinone`、`ortho_benzoquinone→benzoquinone`）；注册表（fg/ring/unsat producers）整体删除，骨架识别改由 `scaffold/ring_core`（@_register）+ `scaffold/specs`；涉及 `architecture/layer2-parent-selector.md`, `architecture/overview.md`, `concepts/functional-group-priority.md`, `concepts/atom-ownership.md`, `guides/adding-new-ring-system.md`, `guides/adding-new-functional-group.md`, `reference/core-data-contracts.md`, `index.md` | chem-dev |
