"""Same forced rejection before/after actual UEL trials on the unchanged host."""
from pathlib import Path
import sys
import json
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from verification.common import compiled
import fehost
import meshes
from hostfe.bounds import IncrementTooLarge


def experiment(after_trials):
    contract=fehost.load_contract(str(ROOT))
    module=compiled(contract)
    nodes,elems=meshes.distorted(8,seed=3)
    el=fehost.ContractUEL(module,contract,len(elems),.4)
    top=np.where(np.isclose(nodes[:,1],1))[0]
    n0=meshes.node_at(nodes,0,0);n1=meshes.node_at(nodes,1,0)
    def bc(frac):
        cb=.35+.45*(1-np.cos(2*np.pi*frac*2))/2
        vals={6*k+2:cb for k in top}
        vals.update({6*n0:0.,6*n0+1:0.,6*n1+1:0.})
        return vals
    control=dict(injected=False,calls=0,seen_trial_states=0)
    def batch(coords,U,DU):
        target=not control["injected"] and abs(el._time[1]-.8)<1e-10 and el.dt>.39
        if target:
            if after_trials:
                # The second call includes a Newton trial; returned SVARS are
                # discarded in assembly, just as in a normal host evaluation.
                result=el.element_rk_batch(coords,U,DU)
                control["seen_trial_states"]+=1
                control["calls"]+=1
                if control["calls"]<2:
                    return result
            control["injected"]=True
            raise IncrementTooLarge("verification forced rejection")
        return el.element_rk_batch(coords,U,DU)
    U0=np.zeros(len(nodes)*6);U0[2::6]=.35
    rows=[];history=[]
    def output(step,frac,U,metrics=None):
        rows.append(U.copy());history.append(el.svars.copy())
    U,report=fehost.solve(nodes,elems,el,bc,40,.4,U0=U0,batch_fn=batch,
        residual_scale=fehost.residual_weights(el,contract),step_callback=output)
    assert control["injected"]
    return U,el.svars.copy(),np.array(rows),np.array(history),report,control


def main():
    before=experiment(False)
    after=experiment(True)
    comparisons={name:bool(np.array_equal(a,b)) for name,a,b in zip(
        ["final_nodal","final_svars","all_output_nodal","all_output_svars"],before[:4],after[:4])}
    assert all(comparisons.values()),comparisons
    assert after[5]["seen_trial_states"]==2
    report=dict(status="passed",exact_equality=comparisons,
        before=before[4].summary(),after=after[4].summary(),
        after_trial_calls=after[5]["seen_trial_states"],
        scope="distorted n=8 seed=3, 2 cycles x 20 steps, forced rejection at t=.8, dt=.4")
    (ROOT/"verification/rollback_checks.json").write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(report,indent=2))


if __name__=="__main__":
    main()
