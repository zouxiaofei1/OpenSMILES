"""神经网络命名引擎封装: knowledgator/SMILES2IUPAC-canonical-base(本地 tools/ml_models)。

mT5 版 SMILES→IUPAC, 仅产英文。该模型对输入 SMILES 拼写敏感(不吃小写芳香记号、
羰基尾缀形式不同会误读), 统一喂 RDKit canonical + Kekulé。进程内单例 lazy 加载,
generate 用锁串行化(CPU 无并发增益)。结果规整成 NameResult 同构 dict, zh 恒空
(模型不产中文), 供 routes_name 统一追加 engine/gold。
"""

from __future__ import annotations

import os
import threading
import time
from pathlib import Path
from typing import Any

_ML_ROOT = Path(__file__).resolve().parents[2] / "tools" / "ml_models"
_ALOAD_LOCK = threading.Lock()
_INFER_LOCK = threading.Lock()
_engine: Any = None


def _ensure() -> Any:
    """进程内单例 lazy 加载本地模型与 tokenizer(首次约数秒)。"""
    global _engine
    if _engine is not None:
        return _engine
    with _ALOAD_LOCK:
        if _engine is None:
            os.environ.setdefault("HF_HUB_OFFLINE", "1")  # 全部本地, 禁止访问 HF
            from chemicalconverters.model_utils import MT5ForConditionalGeneration
            from rdkit import Chem
            from transformers import AutoTokenizer

            class _E:
                pass

            e = _E()
            e.model = MT5ForConditionalGeneration.from_pretrained(
                str(_ML_ROOT / "SMILES2IUPAC-canonical-base"))
            e.stok = AutoTokenizer.from_pretrained(str(_ML_ROOT / "SMILES-FAST-TOKENIZER"))
            e.itok = AutoTokenizer.from_pretrained(str(_ML_ROOT / "IUPAC-FAST-TOKENIZER"))
            e.chem = Chem
            _engine = e
    return _engine


def _preprocess(smiles: str, chem) -> str | None:
    """canonical kekule SMILES(去立体/同位素、无芳香小写); 解析失败返回 None。"""
    mol = chem.MolFromSmiles((smiles or "").strip())
    if mol is None:
        return None
    mol = chem.RWMol(mol)
    try:
        chem.Kekulize(mol, clearAromaticFlags=True)
    except Exception:
        pass
    return chem.MolToSmiles(mol)


def name_result(smiles: str) -> dict[str, Any]:
    """生成该结构的英文 IUPAC 名(仅 en; zh 恒空), 返回 NameResult 同构 dict。"""
    t0 = time.perf_counter()
    try:
        e = _ensure()
        p = _preprocess(smiles, e.chem)
        if p is None:
            return {"en": "", "zh": "", "success": False, "source": "error",
                    "time_ms": (time.perf_counter() - t0) * 1000.0,
                    "meta": {"error": "Invalid SMILES"}}
        with _INFER_LOCK:
            enc = e.stok("<BASE>" + p, return_tensors="pt")
            out = e.model.generate(input_ids=enc["input_ids"], attention_mask=enc["attention_mask"],
                                   max_new_tokens=156, num_beams=1)
            en = e.itok.decode(out[0], skip_special_tokens=True).strip()
    except Exception as exc:
        return {"en": "", "zh": "", "success": False, "source": "error",
                "time_ms": (time.perf_counter() - t0) * 1000.0,
                "meta": {"error": f"{type(exc).__name__}: {exc}"}}
    return {
        "en": en,
        "zh": "",
        "success": bool(en),
        "source": "ml" if en else "error",
        "time_ms": (time.perf_counter() - t0) * 1000.0,
        "meta": {},
    }
