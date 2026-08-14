"""跨层通用的化学拓扑工具，供 L2（母体选择）与 L3（取代基命名）共享。
烷基侧链拓扑（side_alkyl）与芳基叶注册表现在位于 layer3；本模块不依赖任何流水线层。"""
from __future__ import annotations
