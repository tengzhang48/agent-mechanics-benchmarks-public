"""Inspect a captured failed UEL call without altering the solver or element."""
from pathlib import Path
import sys
import json
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from verification.common import compiled
from verification.reference import GP,shape
import fehost


def main(file):
    data=np.load(file,allow_pickle=False)
    c=json.loads((ROOT/"verification/exploration/gauss_intact/uel_contract.json").read_text())
    c={**c,"props":data["props"].tolist()}
    el=fehost.ContractUEL(compiled(c),c,len(data["U"]),float(data["dt"]))
    el._time=data["time"];el.svars=data["svars"].copy()
    rhs,K,snew,pn=el._call_batch(data["coords"],data["U"],data["DU"],el.svars)
    H=snew.T[:,0::3]
    sn=data["U"][:,5::6]
    r64=np.zeros((len(sn),4));r80=np.zeros((len(sn),4),np.longdouble)
    largest={};maxcond=0.
    for k,(xi,eta) in enumerate(GP):
        N,dn=shape(xi,eta)
        sq=sn@N
        s80=sn.astype(np.longdouble)@N.astype(np.longdouble)
        J=np.einsum("mai,ar->mir",data["coords"],dn)
        det=np.linalg.det(J)
        Aq=1+(1-c["props"][4])*H[:,k]
        r64+=(Aq*sq-1)[:,None]*N[None,:]*det[:,None]
        r80+=(Aq.astype(np.longdouble)*s80-1)[:,None]*N.astype(np.longdouble)[None,:]*det.astype(np.longdouble)[:,None]
        cond=(abs(sn)@abs(N))/np.maximum(abs(sq),1e-300)
        m=int(np.argmax(cond))
        if cond[m]>maxcond:
            maxcond=float(cond[m]);largest=dict(element=m,gp=k,H=float(H[m,k]),
                nodal_intactness=sn[m].tolist(),interpolated_s_64=float(sq[m]),
                interpolated_s_80=str(s80[m]),cancellation_condition=maxcond)
    err=np.array(r64.astype(np.longdouble)-r80,float)
    result=dict(pnewdt=float(pn),dt=float(data["dt"]),time=data["time"].tolist(),
        H_min=float(H.min()),H_max=float(H.max()),
        old_H_max=float(el.svars[:,0::3].max()),
        nodal_s_min=float(sn.min()),nodal_s_max=float(sn.max()),
        max_64_vs_80_reaction_error=float(abs(err).max()),
        norm_64_vs_80_reaction_error=float(np.linalg.norm(err)),
        worst_interpolation=largest,
        element_reaction_norm=float(np.linalg.norm(rhs[3::6])),
        note="80-bit interpolation is a diagnostic, not a revised element or an accuracy oracle for H")
    target=Path(file).with_name("failure_diagnosis.json")
    target.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))


if __name__=="__main__":
    main(sys.argv[1])
