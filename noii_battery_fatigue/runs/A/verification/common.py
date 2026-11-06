"""Workspace-local utilities; use the unchanged host exclusively as a host."""
import hashlib
import importlib
import json
import os
from pathlib import Path
import sys
import fcntl

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
os.environ['TMPDIR'] = str(ROOT / 'tmp')
sys.path.insert(0, str(ROOT / 'host'))
import fehost
import meshes


def load_element(n_elem=1, dt=0.2, parameters=None, variant=None):
    contract = fehost.load_contract(str(ROOT))
    if variant:
        contract.update(variant)
    if parameters:
        contract['props'] = [parameters.get(k, v) for k, v in
                             zip(contract['prop_names'], contract['props'])]
    sources = contract['source']
    if isinstance(sources, str):
        sources = [sources]
    digest = hashlib.sha256(b''.join((ROOT / s).read_bytes() for s in sources)
                            + json.dumps(contract['compile_flags']).encode()).hexdigest()[:12]
    name = 'uel_local_' + digest
    build = ROOT / 'tmp' / 'build' / name
    sys.path.insert(0, str(build))
    build.parent.mkdir(parents=True,exist_ok=True)
    with (build.parent/(name+'.lock')).open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        if list(build.glob(name + '*.so')):
            mod = importlib.import_module(name)
        else:
            mod = fehost.build_module(str(ROOT), contract, str(build), name)
    return fehost.ContractUEL(mod, contract, n_elem, dt), contract


def raw_call(el, x, u, old, sv=None, dt=None):
    if sv is None:
        sv = el.svars[0].copy()
    saved_dt = el.dt
    if dt is not None:
        el.dt = dt
    try:
        rhs, k, sn, pn = el._call(0, x, u, u - old, sv.copy())
        return -rhs, k, sn, pn
    finally:
        el.dt = saved_dt


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
