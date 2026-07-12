import pytest
from namepredict.cache.common_names import CommonNameCache
from namepredict.types import NameResult

def test_cache_rejects_over_100():
    c = CommonNameCache(max_entries=100)
    for i in range(100):
        c.put(f"C{i}", NameResult(en=f"n{i}", zh=f"名{i}", success=True, source="cache"))
    with pytest.raises(ValueError):
        c.put("overflow", NameResult(en="x", zh="x", success=True, source="cache"))
