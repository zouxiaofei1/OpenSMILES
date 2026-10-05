# OpenSMILES
基于规则的 SMILES → IUPAC化学命名工具。手动输入，绘制分子或粘贴SMILES，获得IUPAC系统命名

## 安装
python版
使用/src目录

或询问Agent

Webjs版
打开index.min.html可用，基于python代码翻译


## 开发服务器

```bash
uvicorn server.backend.app:app --host 127.0.0.1 --port 8666 --reload
```
或双击restart-console.bat

## Benchmark 与 性能参数

benchmark命令
python -m benchmarks.benchmark_parallel --data benchmarks/merged_benchmark.json --time  --timeout 1 

数据集 | 匹配率 | 平均每分子耗时(测试于9950x@5.2Ghz)
merged_benchmark.json | 82.7% | 20.6ms
pubchem 10000抽样 | 97.3% | 19.3ms

## License

MIT 许可证，详见 [LICENSE](LICENSE)。

## 命名示例

SMILES | EN | ZH | explanation |
BrC1=C(C=C(C=C1)C(C)(C)C)[N+](=O)[O-] | 1-bromo-4-tert-butyl-2-nitrobenzene | 1-溴-4-叔丁基-2-硝基苯 | 简单化合物 |
C1SC21CC2 | 1-thiaspiro[2.2]pentane | 1-硫杂螺[2.2]戊烷 | |
CN1C=NC2=C1C(=O)N(C(=O)N2C)C | 1,3,7-trimethylpurine-2,6-dione | 1,3,7-三甲基嘌呤-2,6-二酮 | 咖啡因
CC1=C2[C@@]([C@]([C@H]([C@@H]3[C@]4([C@H](OC4)C[C@@H]([C@]3(C(=O)[C@@H]2OC(=O)C)C)O)OC(=O)C)OC(=O)c5ccccc5)(C[C@@H]1OC(=O)[C@H](O)[C@@H](NC(=O)c6ccccc6)c7ccccc7)O)(C)C | [(1S,2S,3R,4S,7R,9S,10S,12R,15S)-4,12-diacetyloxy-15-[(2R,3S)-3-benzamido-2-hydroxy-3-phenylpropanoyl]oxy-1,9-dihydroxy-10,14,17,17-tetramethyl-11-oxo-6-oxatetracyclo[11.3.1.03,10.04,7]heptadec-13-en-2-yl] benzoate | 苯甲酸(1S,2S,3R,4S,7R,9S,10S,12R,15S)-[4,12-二乙酰氧基-15-[(2R,3S)-3-苯甲酰胺基-2-羟基-3-苯基丙酰氧基]-1,9-二羟基-10,14,17,17-四甲基-11-氧代-6-氧杂四环[11.3.1.03,10.04,7]十七-13-烯-2-基]酯 | 紫杉醇(和Pubchem结果相同) |