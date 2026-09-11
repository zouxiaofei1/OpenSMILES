from __future__ import annotations

import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))  # E:\chem\src\namepredict
_ROOT = os.path.dirname(os.path.dirname(_HERE))     # E:\chem
sys.path = [p for p in sys.path if p != _HERE]
for _p in (_ROOT, os.path.join(_ROOT, "src")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from namepredict.namer import SMILESNNamer  # noqa: E402


def show(namer: SMILESNNamer, s: str) -> None:
    r = namer.name(s)
    print(f"SMILES : {s}")
    print(f"  EN   : {r.en or '<fail>'}")
    print(f"  ZH   : {r.zh or '<fail>'}")
    if not r.success:
        print(f"  meta : {r.meta}")
    print()


def main() -> None:
    namer = SMILESNNamer()
    args = [a for a in sys.argv[1:] if a.strip()]
    if args:
        for s in args:
            show(namer, s)
    else:
        print("输入 SMILES：")
        for line in sys.stdin:
            s = line.strip()
            if s:
                show(namer, s)


if __name__ == "__main__":
    main()
