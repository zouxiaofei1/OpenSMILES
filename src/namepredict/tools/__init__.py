"""Layer-agnostic chemical topology tools shared by L2 (parent selection)
and L3 (substituent naming).

Some side-topology helpers now live in layer3 (side_alkoxy / aryl_depth2 /
aryl_sub / side_alkoxy) and are re-imported here: side_alkyl imports
layer3.side_alkoxy, and tools.leaves imports layer3.aryl_depth2.  There is
no dependency on layer2/layer4/layer5.
"""
from __future__ import annotations
