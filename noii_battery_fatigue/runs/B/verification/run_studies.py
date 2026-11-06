"""Run public D0 studies via the supplied CLI without editing the host.

All compiler and intermediate output is under tmp/. Only durable scientific
outputs and source/contract provenance are copied into verification/runs/final.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import errno
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
VERSION="selected"
FROZEN_ROOT=ROOT
RUN_MANIFEST=None
CASE_TIMEOUT=None


def publish_record(staged, destination):
    """Expose the final record after scientific files and provenance exist.

    The managed workspace can mount tmp/ on a different filesystem. Keep
    staging there as required, and use a final copy if rename cannot cross it.
    """
    try:
        staged.replace(destination)
    except OSError as failure:
        if failure.errno != errno.EXDEV:
            raise
        shutil.copyfile(staged,destination)


def source_manifest():
    contract=json.loads((ROOT/"element/uel_contract.json").read_text())
    names=contract["source"]+["element/uel_contract.json","element/model.py","element/generate.py"]
    if (ROOT/"element/lump_reaction.py").exists():
        names.append("element/lump_reaction.py")
    return {name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in names}


def execute(case):
    mesh,n,cycles,spc,seed=case
    label=f"D0_{mesh}_{n}_s{seed}_c{cycles}_spc{spc}"
    scratch=ROOT/"tmp/public_runs"/VERSION/label
    durable=ROOT/"verification/runs"/VERSION/label
    durable.mkdir(parents=True,exist_ok=True)
    manifest=RUN_MANIFEST or source_manifest()
    cached=durable/"run.json"
    if cached.exists():
        old=json.loads(cached.read_text())
        if old.get("element_source_sha256")==manifest and old["status"]=="ok":
            return label,"cached",old
    staged=scratch/'durable_run.json'
    if staged.exists():
        prior=json.loads(staged.read_text())
        if prior.get('element_source_sha256')==manifest and prior['status']=='ok':
            # This record is written only after a successful CLI process and
            # complete scientific-output copy. Recover publication, not a solve.
            assert (durable/'observables.csv').exists()
            assert len(list(durable.glob('fields_cycle_*.npz')))==cycles
            publish_record(staged,durable/'run.json')
            print(f'RECOVER {label}: completed simulation record published',flush=True)
            return label,'ok',prior
    scratch.mkdir(parents=True,exist_ok=True)
    cmd=[sys.executable,str(ROOT/"host/run_problem.py"),str(FROZEN_ROOT),"--mesh",mesh,
         "--n",str(n),"--seed",str(seed),"--cycles",str(cycles),
         "--steps-per-cycle",str(spc),"--out",str(scratch)]
    env={**os.environ,"TMPDIR":str(ROOT/"tmp"),"PYTHONDONTWRITEBYTECODE":"1",
         "XDG_CACHE_HOME":str(ROOT/"tmp/cache"),"OPENBLAS_NUM_THREADS":"1"}
    print(f"START {label}",flush=True)
    start=time.monotonic()
    with (durable/"stdout.txt").open("w") as log:
        try:
            result=subprocess.run(cmd,cwd=ROOT,env=env,stdout=log,stderr=subprocess.STDOUT,
                                  timeout=CASE_TIMEOUT)
        except subprocess.TimeoutExpired:
            record=dict(status='process_timeout',timeout_seconds=CASE_TIMEOUT,
                        element_source_sha256=manifest,command=cmd)
            (durable/'controller_timeout.json').write_text(json.dumps(record,indent=2)+'\n')
            return label,'process_timeout',record
    if result.returncode or not (scratch/"run.json").exists():
        return label,"process_failed",dict(returncode=result.returncode)
    for file in scratch.iterdir():
        if file.is_file() and file.name not in ('run.json','durable_run.json') and file.suffix in (".csv",".npz",".json"):
            shutil.copyfile(file,durable/file.name)
    summary=json.loads((scratch/"run.json").read_text())
    summary["element_source_sha256"]=manifest
    summary["command"]=cmd
    summary["elapsed_including_compile_s"]=time.monotonic()-start
    # Publish only after fields/CSV and source identity are ready. Even the
    # staging file stays under tmp/; cross-mount publication uses a final copy.
    staged=scratch/'durable_run.json'
    staged.write_text(json.dumps(summary,indent=1)+"\n")
    publish_record(staged,durable/'run.json')
    print(f'END {label}: {summary["status"]} ({summary["wall_s"]} solve seconds)',flush=True)
    return label,summary["status"],summary


def main():
    global VERSION, FROZEN_ROOT, RUN_MANIFEST, CASE_TIMEOUT
    ap=argparse.ArgumentParser()
    ap.add_argument("--group",default="all",choices=["mesh","time","fine_time","long","long_time","all"])
    ap.add_argument("--workers",type=int,default=2)
    ap.add_argument("--version",default="selected")
    ap.add_argument("--n64-cycles",type=int,default=5)
    ap.add_argument("--long-time-spcs",type=int,nargs='+',default=[80,160,320,640,1280,2560])
    ap.add_argument("--case-timeout-seconds",type=float)
    a=ap.parse_args()
    VERSION=a.version
    CASE_TIMEOUT=a.case_timeout_seconds
    # Freeze this run's element sources before launching subprocesses. The
    # host itself is used directly from its original location, unchanged.
    RUN_MANIFEST=source_manifest()
    digest=hashlib.sha256(json.dumps(RUN_MANIFEST,sort_keys=True).encode()).hexdigest()[:16]
    FROZEN_ROOT=ROOT/"tmp/study_source"/digest
    (FROZEN_ROOT/"element").mkdir(parents=True,exist_ok=True)
    for name,expected in RUN_MANIFEST.items():
        target=FROZEN_ROOT/name
        if not target.exists():
            shutil.copyfile(ROOT/name,target)
        assert hashlib.sha256(target.read_bytes()).hexdigest()==expected
    groups={
        "mesh":[(m,n,5 if n<64 else a.n64_cycles,80,1) for n in (8,16,32,64) for m in ("structured","distorted")],
        "time":[("structured",16,5,s,1) for s in (20,40,160,320,10,640,1280)]
               +[("structured",32,5,s,1) for s in (40,160)],
        "long":[(m,16,20,40,1) for m in ("structured","distorted")]
               +[("structured",32,20,80,1)],
    }
    groups['fine_time']=[('structured',16,5,s,1) for s in (640,1280)]
    groups['long_time']=[('structured',16,20,s,1) for s in a.long_time_spcs]
    cases=sum((groups[k] for k in ('mesh','time','long')),[]) if a.group=="all" else groups[a.group]
    results={}
    out=ROOT/f'verification/studies_{VERSION}_{a.group}.json'
    with ThreadPoolExecutor(max_workers=a.workers) as executor:
        futures={executor.submit(execute,case):case for case in cases}
        for future in as_completed(futures):
            label,status,summary=future.result()
            results[label]=dict(status=status,summary=summary)
            out.write_text(json.dumps(results,indent=1)+"\n")
    if any(v["status"] not in ("ok","cached") for v in results.values()):
        raise SystemExit("one or more studies failed; inspect the recorded evidence")


if __name__=="__main__":
    main()
