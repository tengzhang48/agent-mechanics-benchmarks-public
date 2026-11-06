"""Compile once, entirely inside workspace tmp/, using the unmodified host."""
from pathlib import Path
import sys
import hashlib
import importlib
import json
import fcntl

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "host"))
import fehost


def compiled(contract=None, checked=False):
    c = contract or fehost.load_contract(str(ROOT))
    if checked:
        c = {**c, "compile_flags": c["compile_flags"]+["-fcheck=all", "-fbacktrace"]}
    sources = c["source"] if isinstance(c["source"], list) else [c["source"]]
    digest = hashlib.sha256()
    for file in sources:
        digest.update((ROOT / file).read_bytes())
    digest.update(json.dumps(c.get("compile_flags", [])).encode())
    name = "challenge_"+digest.hexdigest()[:12]
    workdir = ROOT / "tmp" / name
    sys.path.insert(0, str(workdir))
    lockfile=ROOT/"tmp"/(name+".lock")
    with lockfile.open("w") as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        if workdir.exists() and list(workdir.glob(name+"*.so")):
            return importlib.import_module(name)
        return fehost.build_module(str(ROOT), c, str(workdir), name)


def properties(c):
    return dict(zip(c["prop_names"], c["props"]))
