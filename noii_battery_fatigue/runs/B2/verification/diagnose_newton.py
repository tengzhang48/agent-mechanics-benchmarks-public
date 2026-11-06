"""Global replay of a captured trial: Newton direction and branch transitions."""
from pathlib import Path
import sys
import json
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from verification.common import compiled
import fehost
import meshes
from hostfe.driver import assemble_batched,_apply_dirichlet_rows,_scipy_solve


def main(file,contractfile=None):
    z=np.load(file,allow_pickle=False);ne=len(z["U"]);n=round(np.sqrt(ne))
    nodes,elems=meshes.structured(n)
    c=json.loads(Path(contractfile).read_text()) if contractfile else json.loads((ROOT/"verification/exploration/gauss_intact/uel_contract.json").read_text())
    c={**c,"props":z["props"].tolist()}
    el=fehost.ContractUEL(compiled(c),c,len(elems),float(z["dt"]))
    el.svars=z["svars"].copy();el._time=z["time"]
    gm=(elems[:,:,None]*6+np.arange(6)).reshape(ne,24)
    U=np.zeros(len(nodes)*6);DU=np.zeros_like(U)
    U[gm]=z["U"];DU[gm]=z["DU"];Uold=U-DU
    con=np.array(sorted([6*k+2 for k in np.where(np.isclose(nodes[:,1],1))[0]]+[0,1,6*n+1]))
    free=np.ones(len(U),bool);free[con]=False
    weights=np.tile([100.,100.,1.,1.,100.,1.],len(nodes))
    def assemble(u):
        R,K=assemble_batched(nodes,elems,6,el.element_rk_batch,u,u-Uold)
        R[con]=0
        return R,K
    R,K=assemble(U);Kc,Rc=_apply_dirichlet_rows(K,R,con);du=_scipy_solve(Kc,-Rc);du[con]=0
    rows=[]
    for j in range(18):
        alpha=2.**(-j)
        rt,_=assemble(U+alpha*du)
        rows.append(dict(alpha=alpha,weighted_norm=float(np.linalg.norm(weights*rt)),
            damage_norm=float(np.linalg.norm(rt[3::6]))))
    fd=[];direction=du/np.max(abs(du))
    target=K@direction;target[con]=0
    for step in (1e-5,1e-7,1e-9,1e-11,1e-13):
        rp,_=assemble(U+step*direction)
        rm,_=assemble(U-step*direction)
        error=(rp-rm)/(2*step)-target
        fd.append(dict(step=step,relative_error=float(np.linalg.norm(error)/np.linalg.norm(target))))
    result=dict(residual_norm=float(np.linalg.norm(weights*R)),
        component_raw_norms=[float(np.linalg.norm(R[j::6])) for j in range(6)],
        direction_max_by_component=[float(abs(du[j::6]).max()) for j in range(6)],
        line_search=rows,directional_FD=fd)
    target="newton_diagnosis_"+Path(contractfile).parent.name+".json" if contractfile else "newton_diagnosis.json"
    Path(file).with_name(target).write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))


if __name__=="__main__":
    main(*sys.argv[1:])
