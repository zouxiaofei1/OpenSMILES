# OpenSMILES

基于规则的 SMILES → 双语 IUPAC 命名引擎。

## 安装

Ask the Agent.

## 开发服务器

```bash
uvicorn server.backend.app:app --host 127.0.0.1 --port 8765 --reload
```
或双击restart-console.bat

## Benchmark 与 Pytest

pytest
python -m benchmarks.benchmark_parallel --data benchmarks/merged_benchmark.json --time  --timeout 1 


## License

本项目采用 MIT 许可证，详见 [LICENSE](LICENSE)。