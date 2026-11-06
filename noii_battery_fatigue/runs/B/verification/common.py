"""Compilation and paths for verification; never changes the supplied host."""
from pathlib import Path
import hashlib
import importlib
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "host"))
import fehost


def identical_bytes(a, b):
    """Compare scalar bit patterns, including signed zero, in a fixed order."""
    left, right = np.asarray(a), np.asarray(b)
    return (left.dtype == right.dtype and left.shape == right.shape
            and left.tobytes(order='C') == right.tobytes(order='C'))


def source_identity():
    """Hash the exact compiled-source inputs and the runtime contract."""
    contract = fehost.load_contract(str(ROOT))
    sources = contract['source']
    if isinstance(sources, str):
        sources = [sources]
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
            for name in sources+['element/uel_contract.json']}


def compiled():
    contract = fehost.load_contract(str(ROOT))
    sources = contract["source"]
    if isinstance(sources, str):
        sources = [sources]
    digest = hashlib.sha256()
    for source in sources:
        digest.update((ROOT / source).read_bytes())
    digest.update(str(contract.get("compile_flags", [])).encode())
    name = "verified_uel_" + digest.hexdigest()[:12]
    work = ROOT / "tmp" / name
    work.mkdir(parents=True, exist_ok=True)
    if list(work.glob(name + "*.so")):
        sys.path.insert(0, str(work))
        module = importlib.import_module(name)
    else:
        module = fehost.build_module(str(ROOT), contract, str(work), name)
    return module, contract


def direct(module, contract, X, U, previous=None, state=None, dt=.2, time=1., flags=None, props=None):
    if previous is None:
        previous = U
    if state is None:
        state = np.zeros(contract["nsvars"])
    return module.drive_uel(
        np.asarray(state, float).copy(), np.asfortranarray(np.asarray(X, float).T),
        np.asarray(U, float), np.asarray(U-previous, float),
        np.asarray(contract["props"] if props is None else props, float),
        np.asarray(contract["jprops"], np.int32), np.array([time, time]), dt, 1.,
        np.asarray(contract["lflags"] if flags is None else flags, np.int32),
        np.zeros(3), 1, 1, 1, 1, 1.)
