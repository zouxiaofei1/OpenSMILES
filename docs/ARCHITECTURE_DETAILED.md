# NamePredict V2 架构概要

## 处理流水线

```
SMILES → molecule_parser.py → skeleton_identifier_v3.py → fg_detector.py 
       → parent_selector_v3.py → numbering_engine.py → substituent_analyzer.py 
       → name_assembler.py → IUPACName{en, zh}
```

## 核心模块

**1. 入口与常量**: `namer.py`, `constants.py`, `salt_handler.py`

**2. 骨架识别**: `skeleton_identifier_v3.py`, `skeleton_identifier_v2.py`, `skeleton_identifier_v2_chain.py`, `skeleton_identifier_v2_special.py`, `skeleton_identifier_v2_n_oxide.py`, `skeleton_models.py`, `skeleton_creator.py`, `ring_skeleton_creator.py`, `special_skeleton_handler.py`, `special_types.py`, `chain_finder.py`, `chain_utils.py`, `chain_selection_helper.py`, `ring_analyzer.py`, `ring_utils.py`, `spiro_analyzer.py`, `bridged_ring_analyzer.py`

**3. 官能团检测**: `fg_detector.py`, `functional_groups.py`, `functional_group_utils.py`, `fg_position.py`, `_detect_mixin.py`, `_detect_oxygen_mixin.py`, `_detect_nitrogen_mixin.py`, `_detect_sulfur_mixin.py`, `_detect_phosphorus_mixin.py`, `_detect_quinone_mixin.py`, `_detect_cyanate_mixin.py`, `_detect_salt_mixin.py`, `_detect_inorganic_mixin.py`

**4. 母体选择**: `parent_selector_v3.py`, `parent_selector_v2_helper.py`, `parent_selector_with_fg.py`, `parent_naming.py`

**5. 编号引擎**: `numbering_engine.py`, `numbering_cyclic.py`, `numbering_utils.py`, `numbering_types.py`, `fused_ring_numbering.py`, `fused_ring_numbering_base.py`, `fused_ring_numbering_benzo.py`, `fused_ring_numbering_indole.py`, `fused_ring_numbering_quinoline.py`, `fused_ring_numbering_quinazoline.py`, `fused_ring_numbering_phthalide.py`, `fused_ring_types.py`

**6. 取代基分析**: `substituent_analyzer.py`, `substituent_formatter.py`, `substituent_constants.py`, `substituent_branched.py`, `substituent_element_base.py`, `substituent_element_analyzers.py`, `substituent_element_halogen.py`, `substituent_element_nitrogen.py`, `substituent_element_sulfur.py`, `substituent_element_phenyl.py`, `substituent_element_fused_heterocyclic.py`, `substituent_element_oxygen_helpers.py`, `substituent_element_oxygen_benzyloxy.py`, `substituent_fg_amide.py`, `substituent_fg_amino.py`, `substituent_fg_oxo.py`, `substituent_fg_isocyanato.py`, `substituent_fg_simple.py`, `substituent_heterocyclic.py`, `substituent_n_alkyl.py`, `substituent_sulfonyl_handler.py`

**7. 名称组装**: `name_assembler.py`, `assembler_core.py`, `base_name_builder.py`, `iupac_name.py`, `trivial_names.py`, `unsaturation_utils.py`, `alkyl_utils.py`

**8. 杂环系统**: `heterocyclic_detector.py`, `heterocyclic_naming.py`, `heterocyclic_assembler.py`, `heterocyclic_numbering.py`, `heterocyclic_types.py`, `_heterocycle_identifier.py`, `_heterocyclic_utils.py`, `_saturated_heterocycle.py`, `_oxadiazole_numbering.py`

**9. 特殊处理器**: `ester_assembler.py`, `carbonate_handler.py`, `orthoester_handler.py`, `sulfonic_assembler.py`, `sulfone_handler.py`, `sulfide_sulfoxide_handler.py`, `phosphine_handler.py`, `phosphate_handler.py`, `carbamate_urea_assembler.py`, `boronic_acid_assembler.py`, `quaternary_ammonium_handler.py`, `amine_n_oxide_assembler.py`, `ether_handler.py`, `halide_detector.py`, `halide_detection.py`, `salicylic_pattern.py`, `salt_detection.py`, `arsine_handler.py`, `bismuthine_handler.py`, `silane_handler.py`, `stibine_handler.py`, `pnictogen_base.py`

**10. 缓存**: `cache/common_names.py`, `cache/exceptions.py`

**11. 规则**: `rules/functional_groups.py`, `rules/iupac_priority.py`, `rules/naming_utils.py`

**12. 测试**: `tests/`, `tests-inorganic/`, `tests-pending/`, `data/smiles_tiers.json`, `data/smiles_tiers_inorganic.json`

## 优先级规则

羧酸 > 酐 > 酯 > 酰卤 > 酰胺 > 腈 > 醛 > 酮 > 醇 > 胺 > 炔 > 烯 > 醚 > 卤素
