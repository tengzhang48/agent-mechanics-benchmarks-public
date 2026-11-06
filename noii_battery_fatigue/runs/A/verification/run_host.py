"""Run the original host CLI, keep build/temp files under tmp/, copy results.

Example: python verification/run_host.py D0_s8 --n 8 --cycles 20
Remaining arguments go unchanged to host/run_problem.py.
"""
from pathlib import Path
import shutil
import subprocess
import sys
import os
import json
import hashlib
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[1]
name = sys.argv[1]
out = ROOT / 'tmp' / 'runs' / name
dest = ROOT / 'verification' / 'runs' / name
out.mkdir(parents=True, exist_ok=True)
env = dict(os.environ, TMPDIR=str(ROOT/'tmp'), PYTHONDONTWRITEBYTECODE='1', OPENBLAS_NUM_THREADS='1')
contract=json.loads((ROOT/'element'/'uel_contract.json').read_text())
metadata={'started_utc':datetime.now(timezone.utc).isoformat(),'contract':contract,
          'source_sha256':hashlib.sha256((ROOT/'element'/'uel.f90').read_bytes()).hexdigest(),
          'arguments':sys.argv[2:]}
(out/'input_manifest.json').write_text(json.dumps(metadata,indent=2)+'\n')
cmd = [sys.executable, str(ROOT/'host'/'run_problem.py'), str(ROOT),
       '--out', str(out), *sys.argv[2:]]
with (out/'console.log').open('w') as log:
    completed = subprocess.run(cmd, cwd=ROOT, env=env, stdout=log, stderr=subprocess.STDOUT)
dest.mkdir(parents=True, exist_ok=True)
for p in out.iterdir():
    if p.is_file():
        shutil.copy2(p, dest/p.name)
print((out/'console.log').read_text(), flush=True)
if completed.returncode:
    raise SystemExit(completed.returncode)
