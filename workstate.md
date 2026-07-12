##日志
[#243ece8][IUPAC P-15.1.7 / 中文位次格式] 一元醇中文位次补连字符(丙-1-醇) [+3 zh dual, benchmark dual 0.3%→0.4%]
[#143333f][IUPAC P-44 母体链选择] 醇母体最长链穿过OH碳(propan-2-ol/butan-2-ol) [+5 dual, benchmark dual 0.4%→0.5%]
[#3352d9a][IUPAC P-29.3.1] 通用直链单烷基C1–C4侧链前缀贯通L3/L4/L5 [+10 tests, dual 0.5%→0.9%]
[#ff557f5][IUPAC P-61.3.1] 通用卤素前缀fluoro/chloro/bromo/iodo贯通L2/L3/L5 [+16 tests, dual 0.9%→1.2%]
[#4359897][IUPAC P-14.2.1 / P-16.3] 倍数词头扩至penta–deca(五–十) [+9 tests, dual 1.2%(48)→1.2%(49), fails 4014→4013]
[#9b6c9b5][IUPAC P-65.1.1] 一元羧酸：L1羧基检测+排除羧基OH；L2/L4/L5 acid母体与词干 [+8 tests, dual 1.2%(49)→1.5%(60), fails 4013→4002]
[#1bbe283][IUPAC P-64.2.1] 一元酮alkanone：L1酮羰基+L2/L4/L5 …-n-one/…-n-酮 [+11 tests, dual 1.5%(60)→1.7%(68), fails 4002→3994]
[#pending][IUPAC P-66.6.1] 一元醛alkanal：L1醛羰基+L2/L4/L5 保留formaldehyde/acetaldehyde与…anal/…醛 [+10 tests, dual 1.7%(68)→1.8%(74), fails 3994→3988]

## 其他
- 基线(本轮前): en=1.7%(68/4062) zh=7.3%(55/756) dual=1.7%(68/4062)
- 本次: en=1.8%(74/4062) zh=8.1%(61/756) dual=1.8%(74/4062) fails=3988
- 改动: L1 aldehydes/has_aldehyde（=O且恰1碳邻、非羧基）；L2 acid>aldehyde>ketone>alcohol>alkane；L4醛碳=1；L5 ALDEHYDE stems
- 验证: C=O→formaldehyde/甲醛；CC=O→acetaldehyde/乙醛；CCC=O→propanal/丙醛；酮/酸/醇不回退
- 已知缺口: 羟基醛前缀；烯醛/二醛；酯/酸酐；支化烷基isopropyl；胺
- 下一步候选: 酯；或支化烷基；或胺；或羟基醛前缀
