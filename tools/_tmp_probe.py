from namepredict.namer import SMILESNNamer
n = SMILESNNamer()
for s in [
    "BrCCCC1=CC=CC=C1",
    "BrC1=CC(=C(C=C1)OCCOC)F",
    "BrC1=CC=C2C=CC(=NC2=C1)C(F)(F)F",
    "ClC1=NC2=CC=CC=C2C(=C1)C(F)(F)F",
    "CCc1ccccc1",
    "FC(F)(F)c1ccccc1",
    "NCCN",
    "COc1ccccc1",
    "CCOc1ccccc1",
]:
    r = n.name(s)
    print(repr(s), "->", getattr(r, "en", None), "/", getattr(r, "zh", None), "ok", r.success)
