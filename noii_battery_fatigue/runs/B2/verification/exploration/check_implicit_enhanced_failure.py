"""Check the current-H alternative's saved failure against the separate oracle."""
from pathlib import Path
import sys
import json
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from verification.exploration.implicit_enhanced import oracle
from verification.common import compiled,properties
import fehost


if __name__=='__main__':
    out=ROOT/'verification/exploration/implicit_enhanced_severe32'
    z=np.load(out/'failed_trial.npz',allow_pickle=False)
    c=json.loads((out/'contract.json').read_text());c['props']=z['props'].tolist()
    el=fehost.ContractUEL(compiled(c,checked=True),c,len(z['U']),float(z['dt']))
    el._time=z['time'];p=properties(c)
    rhs,K,sn,pn=el._call_batch(z['coords'],z['U'],z['DU'],z['svars'].copy())
    R=-rhs.T;states=sn.T
    chosen=np.argsort(np.max(abs(R[:,3::6]),axis=1))[-12:]
    errors=[];state_errors=[]
    for i in chosen:
        rr,ss=oracle(z['coords'][i],z['U'][i],z['U'][i]-z['DU'][i],z['svars'][i],float(z['dt']),p)
        errors.append(float(np.max(abs(rr-R[i]))))
        state_errors.append(float(np.max(abs(ss-states[i]))))
        assert np.allclose(rr,R[i],rtol=1e-9,atol=2e-11),(i,errors[-1])
        assert np.allclose(ss,states[i],rtol=1e-9,atol=2e-9),(i,state_errors[-1])
    result=dict(status='passed',checked_elements=chosen.tolist(),pnewdt=float(pn),
        residual_max_abs_error=max(errors),state_max_abs_error=max(state_errors),
        nodal_intactness_min=float(np.min(z['U'][:,5::6])),
        note='Independent residual/state at twelve large local phase-residual elements in the actual failed trial. This is not a claim that a classical tangent exists at a history switch.')
    (out/'independent_failure_check.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
