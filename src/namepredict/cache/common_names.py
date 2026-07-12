from __future__ import annotations
from namepredict.types import NameResult


class CommonNameCache:
    def __init__(self, max_entries: int = 100) -> None:
        self.max_entries = max_entries
        self._data: dict[str, NameResult] = {}

    def get(self, smiles: str) -> NameResult | None:
        return self._data.get(smiles)

    def put(self, smiles: str, result: NameResult) -> None:
        if smiles in self._data:
            self._data[smiles] = result
            return
        if len(self._data) >= self.max_entries:
            raise ValueError(f"cache exceeds max_entries={self.max_entries}")
        self._data[smiles] = result

    def __len__(self) -> int:
        return len(self._data)
