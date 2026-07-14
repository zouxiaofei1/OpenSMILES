#!/usr/bin/env python
import os
import sys
os.environ.setdefault("PYTHONPATH", "src")
sys.path.insert(0, "src")
import pytest
sys.exit(pytest.main(["tests/unit/test_arene_side_extend.py", "-q", "--tb=short"]))
