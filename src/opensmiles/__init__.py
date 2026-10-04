"""基于规则的 SMILES → 双语 IUPAC 命名。"""

from opensmiles.tools.rdkit_fast import install as _install_rdkit_fast

_install_rdkit_fast()

__version__ = "0.1.0"
