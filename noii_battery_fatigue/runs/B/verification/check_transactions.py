"""Force the same host cutback before UEL evaluation or after trial commit.

The element is compiled unchanged. The fault injection belongs to the test
driver; the unmodified host performs its own state/solution rollback.
"""
import json

import numpy as np

from common import ROOT, compiled, fehost, source_identity, identical_bytes
import meshes
from hostfe.bounds import IncrementTooLarge


def scenario(after):
    module, contract=compiled()
    prob=json.loads((ROOT/"host/dev_problem/d0.json").read_text())
    params=prob["parameters"]
    contract={**contract,"props":[params.get(k,v) for k,v in zip(contract["prop_names"],contract["props"])]}
    nodes,elems=meshes.distorted(8,seed=5)
    dt=.4;spc=20;nsteps=2*spc
    el=fehost.ContractUEL(module,contract,len(elems),dt)
    dpn=el.dof_per_node; comp={n:el.comp(n) for n in fehost.REQUIRED_DOFS}
    top=np.where(np.isclose(nodes[:,1],1.))[0]
    a=meshes.node_at(nodes,0,0);b=meshes.node_at(nodes,1,0)
    U0=np.zeros(len(nodes)*dpn);U0[comp["c"]::dpn]=params["c0"]
    forced=False;rows=[];accepted=[]
    slots={k:np.asarray(v)-1 for k,v in contract["history_slots"].items()}
    old_sv=el.svars.copy()
    minimum={k:0. for k in slots}

    def is_target():
        return (not forced and np.isclose(el._time[1],2.0,atol=1e-13) and
                np.isclose(el.dt,.4,atol=1e-13))

    def bc(frac):
        t=frac*nsteps*dt
        cb=params["c0"]+.45*(1-np.cos(2*np.pi*t/8))/2
        out={int(dpn*k+comp["c"]):cb for k in top}
        out[dpn*a+comp["ux"]]=0.;out[dpn*a+comp["uy"]]=0.;out[dpn*b+comp["uy"]]=0.
        return out

    def batch(coords,U,DU):
        nonlocal forced
        if not after and is_target():
            forced=True
            raise IncrementTooLarge("verification injection before UEL trial evaluation")
        return el.element_rk_batch(coords,U,DU)

    def commit(nodes_,elems_,dpn_,U,Uprev):
        nonlocal forced,old_sv
        el.commit_step(nodes_,elems_,dpn_,U,Uprev)
        if after and is_target():
            forced=True
            raise IncrementTooLarge("verification injection after converged UEL trial state write")
        for k,s in slots.items():
            minimum[k]=min(minimum[k],float(np.min(el.svars[:,s]-old_sv[:,s])))
        old_sv=el.svars.copy()
        accepted.append((el._time[1],el.dt))

    def observe(step,frac,U,metrics=None):
        rows.append(U.copy())

    U,report=fehost.solve(nodes,elems,el,bc,nsteps,dt,U0=U0,
                           residual_scale=fehost.residual_weights(el,contract),
                           batch_fn=batch,commit_fn=commit,step_callback=observe)
    assert forced and report.cutbacks>=1
    assert all(v>=0 for v in minimum.values())
    return U,el.svars.copy(),np.asarray(rows),np.asarray(accepted),report.summary(),minimum


def main():
    before=scenario(False);after=scenario(True)
    assert all(identical_bytes(a,b) for a,b in zip(before[:4],after[:4]))
    summary=dict(bitwise_identical_final_U=True,bitwise_identical_final_svars=True,
                 bitwise_identical_all_output_fields=True,bitwise_identical_accepted_times=True,
                 before_summary=before[4],after_summary=after[4],history_min_change=after[5],
                 forced_reject_at_time=2.,forced_reject_dt=.4,
                 mesh="distorted",n=8,seed=5,cycles=2,steps_per_cycle=20,
                 element_source_sha256=source_identity(),
                 comparison_method='Identical dtype, shape, and C-ordered scalar bytes, including signed zero')
    (ROOT/"verification/transactions.json").write_text(json.dumps(summary,indent=2)+"\n")
    print(json.dumps(summary,indent=2))


if __name__=="__main__":
    main()
