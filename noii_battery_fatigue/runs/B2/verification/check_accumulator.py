"""Compiled regression for accumulator roundoff at tiny positive W increments."""
from pathlib import Path
import sys
import json
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from verification.common import compiled
import fehost


def main():
    archive=ROOT/"verification/exploration/alpha_roundoff_before"
    old=json.loads((archive/"uel_contract.json").read_text())
    old["source"]=[str((archive/Path(s).name).relative_to(ROOT)) for s in old["source"]]
    (archive/"uel_contract.json").write_text(json.dumps(old,indent=2)+"\n")
    new=fehost.load_contract(str(ROOT))
    mods=[compiled(old),compiled(new)]
    X=np.array([[0.,0.],[1.,0.],[1.,1.],[0.,1.]])
    U=np.zeros((4,6));U[:,0]=.01*X[:,1];U[:,2]=.35;U[:,5]=1.
    def call(m,sv):
        return m.drive_uel(sv.copy(),X.T,U.ravel(),np.zeros(24),np.array(new["props"]),
            np.array([],np.int32),np.array([1.,1.]),.1,1.,np.array(new["lflags"],np.int32),
            np.zeros(3),1,1,1,1,1.)
    W=call(mods[0],np.zeros(12))[2][2::3]
    rng=np.random.default_rng(12091)
    state=np.zeros((4,3));state[:,0]=1e6
    broken=[]
    for k,w in enumerate(W):
        wp=np.nextafter(w,-np.inf)
        found=False
        for a in rng.uniform(w,100*w,100000):
            b=(a+w)-wp
            if b<a:
                state[k,1]=a;state[k,2]=wp
                broken.append(float(b-a));found=True
                break
        assert found,"No roundoff witness found"
    prior=state.ravel()
    ro,ko,so,po=call(mods[0],prior)
    rn,kn,sn,pn=call(mods[1],prior)
    assert po<1 and pn==1,(po,pn)
    assert np.all(sn[0::3]>=prior[0::3])
    assert np.all(sn[1::3]>=prior[1::3])
    expected=prior[1::3]+np.maximum(sn[2::3]-prior[2::3],0)
    assert np.array_equal(sn[1::3],expected)
    result=dict(status="passed",old_pnewdt=float(po),corrected_pnewdt=float(pn),
        old_alpha_change_before_guard=broken,old_alpha=prior[1::3].tolist(),
        W_old=prior[2::3].tolist(),W_new=sn[2::3].tolist(),
        corrected_alpha_change=(sn[1::3]-prior[1::3]).tolist(),
        correction="alpha_old + (W_new-W_old) on the increasing branch; preserve operation order")
    (ROOT/"verification/accumulator_checks.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))


if __name__=="__main__":
    main()
