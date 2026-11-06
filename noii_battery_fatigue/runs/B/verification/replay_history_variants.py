"""Isolated exact-equality history-branch replays; production files unchanged."""
import argparse
import hashlib
import importlib
from pathlib import Path
import sys

from common import ROOT, compiled, fehost
from run_case import run


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("variant",choices=["geH","geW","geBoth"])
    args=parser.parse_args()
    _,contract=compiled()
    # Freeze the original fully implicit source used by the captured failure.
    frozen=ROOT/"tmp/history_tangent_geW_cs10/generated_model.for"
    source=frozen.read_text().replace(
        "IF (DBLE(w_new) .GE. DBLE(wprev_old)) THEN",
        "IF (DBLE(w_new) .GT. DBLE(wprev_old)) THEN")
    assert hashlib.sha256(source.encode()).hexdigest()=="744b5303512e26f261cd5d3158f85e2d8492fb20712bc2a901cb66062812e1ca"
    if args.variant in ("geH","geBoth"):
        source=source.replace("IF (DBLE(candidate) .GT. DBLE(hhist_old)) THEN",
                              "IF (DBLE(candidate) .GE. DBLE(hhist_old)) THEN")
    if args.variant in ("geW","geBoth"):
        source=source.replace("IF (DBLE(w_new) .GT. DBLE(wprev_old)) THEN",
                              "IF (DBLE(w_new) .GE. DBLE(wprev_old)) THEN")
    folder=ROOT/"tmp"/f"notch_replay_{args.variant}"
    folder.mkdir(parents=True,exist_ok=True)
    (folder/"generated_model.for").write_text(source)
    name=f"notch_replay_{args.variant}"
    work=folder/"build"
    custom={**contract,"source":[str((folder/"generated_model.for").relative_to(ROOT)),"element/uel.f90"]}
    if list(work.glob(name+"*.so")):
        sys.path.insert(0,str(work));module=importlib.import_module(name)
    else:
        module=fehost.build_module(str(ROOT),custom,str(work),name)
    run(ROOT/"verification/problems/notch_sides.json",mesh="distorted",n=16,seed=3,
        out=ROOT/"verification/notch_tie_replays"/args.variant,
        module_contract=(module,custom))


if __name__=="__main__":
    main()
