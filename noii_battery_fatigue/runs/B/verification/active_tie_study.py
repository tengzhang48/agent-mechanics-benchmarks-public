"""Prepare and test the fully implicit active-H-tie candidate under tmp only."""
import argparse
import hashlib
import importlib
import json
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np

from common import ROOT, direct, fehost
from run_case import run

CANDIDATE=ROOT/"tmp/active_tie_candidate"
OUTPUT=ROOT/"verification/active_tie_study"


def prepare():
    (CANDIDATE/"element").mkdir(parents=True,exist_ok=True)
    (CANDIDATE/"verification").mkdir(exist_ok=True)
    text=(ROOT/"element/model.py").read_text()
    assert text.count("if candidate.real > hhist_old.real:")==1
    assert text.count("afat_new = afat_old + w_new - wprev_old")==1
    text=text.replace("if candidate.real > hhist_old.real:","if candidate.real >= hhist_old.real:")
    text=text.replace("afat_new = afat_old + w_new - wprev_old",
                      "afat_new = afat_old + (w_new - wprev_old)")
    (CANDIDATE/"element/model.py").write_text(text)
    for name in ("uel.f90","uel_contract.json"):
        shutil.copyfile(ROOT/"element"/name,CANDIDATE/"element"/name)
    generator=(ROOT/"element/generate.py").read_text()
    generator=generator.replace('str(ROOT / "tools/abaqus_ufl")',f'str(Path({str(ROOT / "tools/abaqus_ufl")!r}))')
    (CANDIDATE/"element/generate.py").write_text(generator)
    subprocess.run([sys.executable,str(CANDIDATE/"element/generate.py")],cwd=ROOT,check=True)
    manifest=dict(base_generated_sha256=hashlib.sha256((ROOT/"element/generated_model.for").read_bytes()).hexdigest(),
                  candidate_generated_sha256=hashlib.sha256((CANDIDATE/"element/generated_model.for").read_bytes()).hexdigest(),
                  changes=["At exact equality select current H candidate derivative (>=)",
                           "Add the nonnegative W increment to alpha with explicit parentheses"],
                  physics="unchanged fully implicit trial-H phase equation",tool_modifications=[])
    OUTPUT.mkdir(parents=True,exist_ok=True)
    (OUTPUT/"candidate_manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")


def module_contract():
    contract=fehost.load_contract(str(CANDIDATE))
    digest=hashlib.sha256()
    for source in contract["source"]:
        digest.update((CANDIDATE/source).read_bytes())
    name="active_tie_"+digest.hexdigest()[:12]
    work=CANDIDATE/"build"
    if list(work.glob(name+"*.so")):
        sys.path.insert(0,str(work));module=importlib.import_module(name)
    else:
        module=fehost.build_module(str(CANDIDATE),contract,str(work),name)
    return module,contract


def uniform_history_checks(module,contract):
    X=np.array([[0.,0.],[1.1,-.05],[.95,.85],[-.12,1.03]])
    ee=np.array([[.03,.012],[.012,-.01]])
    p=dict(zip(contract["prop_names"],contract["props"]))
    props=np.asarray(contract["props"]).copy()
    props[8],props[9]=.002,.0005
    base=np.array([1.4492080900029878,.0020754914424242423,.0020754914424242423])
    unload=base.copy();unload[2]*=.25
    reload=np.array([2.1337751328249275,.003632110024242424,base[2]])
    floor=np.array([6.918181818181818e14,.10107549144242424,base[2]])
    results={}
    for name,scale,old,expected in [
        ("loading",1.,[.3,.001,.001],base),
        ("unloading",.5,base,unload),
        ("reloading",1.,unload,reload),
        ("floor",1.,[.3,.1,.001],floor)]:
        U=np.zeros((4,6));U[:,:2]=X@(scale*ee+p["beta"]*(.55-p["c0"])*np.eye(2)).T
        U[:,2]=.55;U[:,3]=.4;U[:,4]=.3600064*.11837121212121211*scale;U[:,5]=.6
        rhs,K,sv,pn=direct(module,contract,X,U.ravel(),state=np.tile(old,4),props=props)
        error=np.max(abs(sv.reshape(4,3)-expected)/np.maximum(1.,abs(expected)))
        assert pn==1. and error<1e-12,(name,pn,error)
        results[name]=dict(max_state_scaled_error=float(error),expected=expected.tolist())
    # Explicitly exercise the chosen derivative AT an exact H equality.
    # A pilot call is clearly active (H_old=0, W_old=1 suppresses fatigue).
    # Seed its exact per-point output back as old history. The new values
    # and RHS stay the same, and >= must retain the pilot's active tangent.
    U=np.zeros((4,6));U[:,:2]=X@(ee+p["beta"]*(.55-p["c0"])*np.eye(2)).T
    U[:,2]=.55;U[:,3]=.4;U[:,4]=.3600064*.11837121212121211;U[:,5]=.6
    pilot=direct(module,contract,X,U.ravel(),state=np.tile([0.,0.,1.],4),props=props)
    tied=direct(module,contract,X,U.ravel(),state=pilot[2],props=props)
    damage_rows=np.arange(3,24,6)
    displacement_columns=np.array([6*k+j for k in range(4) for j in (0,1)])
    active_block=pilot[1][np.ix_(damage_rows,displacement_columns)]
    tied_block=tied[1][np.ix_(damage_rows,displacement_columns)]
    block_error=float(np.max(abs(active_block-tied_block)))
    assert pilot[3]==1. and tied[3]==1.
    assert np.array_equal(pilot[0],tied[0])
    assert np.max(abs(active_block))>.01 and block_error<1e-10,block_error
    results["exact_H_equality_active_derivative"]=dict(
        rhs_bitwise_unchanged=True,active_damage_u_block_inf=float(np.max(abs(active_block))),
        equality_vs_active_block_error=block_error,
        deliberately_frozen_zero_block_error=float(np.max(abs(active_block))))
    (OUTPUT/"uniform_history_checks.json").write_text(json.dumps(results,indent=2)+"\n")


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("case",choices=["prepare","checks","notch_structured_s1",
                                       "notch_distorted_s1","notch_distorted_s3","D0","floor"])
    args=parser.parse_args()
    if args.case=="prepare":
        prepare(); module_contract();return
    module,contract=module_contract()
    if args.case=="checks":
        import check_element
        target=OUTPUT/"assembled_checks"
        (target/"verification").mkdir(parents=True,exist_ok=True)
        check_element.compiled=lambda:(module,contract)
        check_element.ROOT=target
        check_element.main()
        uniform_history_checks(module,contract)
        return
    if args.case.startswith("notch_"):
        mesh,seed=("structured",1) if args.case=="notch_structured_s1" else ("distorted",int(args.case[-1]))
        run(ROOT/"verification/problems/notch_sides.json",mesh=mesh,n=16,seed=seed,
            cycles=2,spc=80,out=OUTPUT/args.case,module_contract=(module,contract))
    elif args.case=="D0":
        run(ROOT/"host/dev_problem/d0.json",n=8,cycles=1,spc=80,
            out=OUTPUT/args.case,module_contract=(module,contract))
    else:
        run(ROOT/"verification/problems/uniform_fatigue_floor.json",n=8,cycles=2,spc=80,
            uniform_strain=.03,out=OUTPUT/args.case,module_contract=(module,contract))


if __name__=="__main__":
    main()
