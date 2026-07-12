##日志
[#243ece8][IUPAC P-15.1.7 / 中文位次格式] 一元醇中文位次补连字符(丙-1-醇) [+3 zh dual, benchmark dual 0.3%→0.4%]
[#143333f][IUPAC P-44 母体链选择] 醇母体最长链穿过OH碳(propan-2-ol/butan-2-ol) [+5 dual, benchmark dual 0.4%→0.5%]
[#3352d9a][IUPAC P-29.3.1] 通用直链单烷基C1–C4侧链前缀贯通L3/L4/L5 [+10 tests, dual 0.5%→0.9%]
[#ff557f5][IUPAC P-61.3.1] 通用卤素前缀fluoro/chloro/bromo/iodo贯通L2/L3/L5 [+16 tests, dual 0.9%→1.2%]
[#4359897][IUPAC P-14.2.1 / P-16.3] 倍数词头扩至penta–deca(五–十) [+9 tests, dual 1.2%(48)→1.2%(49), fails 4014→4013]
[#9b6c9b5][IUPAC P-65.1.1] 一元羧酸：L1羧基检测+排除羧基OH；L2/L4/L5 acid母体与词干 [+8 tests, dual 1.2%(49)→1.5%(60), fails 4013→4002]

## 其他
- 基线(本轮前): en=1.2%(49/4062) zh=5.4%(41/756) dual=1.2%(49/4062)
- 本次: en=1.5%(60/4062) zh=6.5%(49/756) dual=1.5%(60/4062) fails=4002
- 改动: L1 carboxyls/has_acid；羟基排除羧基OH；L2 kind=acid链穿羧基碳；L4 carboxyl=1；L5 ACID stems(含formic/acetic保留)
- 验证: CC(=O)O→acetic acid/乙酸；CCC(=O)O→propanoic acid/丙酸；CC(C)C(=O)O→2-methylpropanoic acid；CCO乙醇不回退
- 已知缺口: 酯/酸酐/烯酸/二元酸；支化烷基isopropyl；酮/醛/胺
- 下一步候选: 简单酮(alkanone)；或一元醛；或酯(需酸依赖)；或支化烷基
