"""跨层通用的化学拓扑工具，供 L2（母体选择）与 L3（取代基命名）共享。
烷基侧链与芳基叶注册表已移至 layer3；本模块不依赖流水线层。"""
from __future__ import annotations
