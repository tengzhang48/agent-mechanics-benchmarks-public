"""Converged-state DTIME-limit test of the optional H accuracy controller.

Restarts the controller's last accepted state using the uncontrolled lumped
lag element. This distinguishes an unavoidable lag relaxation from an
exaggerated unconverged Newton trial. No production files are modified.
"""
import csv
import importlib.util
import json

import numpy as np

from common import ROOT, fehost


def main():
    directory=ROOT/"tmp/lumped_lag_h_control_candidate/notch_distorted32_s3_c1_spc80"
    z=np.load(directory/"last_state.npz")
    rows=list(csv.DictReader((directory/"increments.csv").open()))
    time0=float(rows[-1]["t"])
    nodes,elems,initial,sv=z["nodes"],z["elems"],z["U"],z["svars"]
    spec=importlib.util.spec_from_file_location("review_lag_common",ROOT/"tmp/lumped_lagged_candidate/common.py")
    common=importlib.util.module_from_spec(spec);spec.loader.exec_module(common)
    module,contract=common.compiled()
    parameters=json.loads((ROOT/"verification/problems/notch_sides.json").read_text())["parameters"]
    contract={**contract,"props":[parameters.get(k,v)for k,v in zip(contract["prop_names"],contract["props"])]}
    x,y=nodes.T
    Hslots=np.asarray(contract["history_slots"]["H"])-1
    Aslots=np.asarray(contract["history_slots"]["alpha"])-1
    a=1-parameters["kappa"]
    results=[]
    for dt in (.00625,.0001,.00001,.000001):
        el=fehost.ContractUEL(module,contract,len(elems),dt);el.svars=sv.copy()
        def bc(fraction):
            t=time0+fraction*dt
            concentration=parameters["c0"]+.6*(1-np.cos(2*np.pi*t/8))/2
            b={6*k+2:concentration for k in np.where(np.isclose(x,0)|np.isclose(x,1))[0]}
            b.update({6*k+3:1. for k in np.where(np.isclose(y,.5)&(x<=.125+1e-12))[0]})
            for xx in (0.,1.):
                k=np.where(np.isclose(x,xx)&np.isclose(y,0))[0][0]
                b[6*k]=b[6*k+1]=0.
            return b
        U,report=fehost.solve(nodes,elems,el,bc,1,dt,U0=initial,tol=1e-9,
                              residual_scale=fehost.residual_weights(el,contract),max_cutbacks=0)
        hold,hnew=sv[:,Hslots],el.svars[:,Hslots]
        proxy=1/(1+a*hold)-1/(1+a*hnew)
        ei,gp=np.unravel_index(np.argmax(proxy),proxy.shape)
        row=dict(dt=dt,physical_time_start=time0,history_proxy_max=float(proxy[ei,gp]),
                 element=int(ei),gp=int(gp),H_old=float(hold[ei,gp]),H_new=float(hnew[ei,gp]),
                 alpha_old_over_threshold=float(sv[ei,Aslots[gp]]/parameters["alpha_T"]),
                 alpha_new_over_threshold=float(el.svars[ei,Aslots[gp]]/parameters["alpha_T"]),
                 nodal_damage_change_max=float(np.max(abs(U[3::6]-initial[3::6]))),
                 solver=report.summary())
        results.append(row);print(json.dumps(row),flush=True)
    (ROOT/"verification/lag_controller_deadlock_review.json").write_text(json.dumps(results,indent=2)+"\n")


if __name__=="__main__":
    main()
