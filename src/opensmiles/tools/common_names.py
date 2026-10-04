"""常用名结果缓存：以 SMILES 为键的内存字典，容量满时抛错。"""

from __future__ import annotations
from opensmiles.types import NameResult


class CommonNameCache:
    """SMILES → NameResult 的内存缓存，达上限后拒绝写入。"""

    def __init__(self, max_entries: int = 100) -> None:
        """初始化缓存容量上限与底层字典。"""
        self.max_entries = max_entries
        self._data: dict[str, NameResult] = {}

    def get(self, smiles: str) -> NameResult | None:
        """按 SMILES 键取缓存结果，未命中返回 None。"""
        return self._data.get(smiles)

    def put(self, smiles: str, result: NameResult) -> None:
        """写入缓存；已存在则覆盖，容量满时抛 ValueError。"""
        if smiles in self._data:
            self._data[smiles] = result
            return
        if len(self._data) >= self.max_entries:
            raise ValueError(f"cache exceeds max_entries={self.max_entries}")
        self._data[smiles] = result

    def clear(self) -> None:
        """清空全部缓存条目（并行评测每行独立命名时用于隔离跨行状态）。"""
        self._data.clear()