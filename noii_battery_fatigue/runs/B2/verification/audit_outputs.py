"""Audit persisted field coverage and independently recompute host observables."""
from pathlib import Path
import sys
import json
import hashlib
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from verification.reference import GP,shape


def audit(folder):
    summary=json.loads((folder/"run.json").read_text())
    assert summary["status"]=="ok", summary["status"]
    for name,expected in summary.get("source_sha256",{}).items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==expected,name
    csv=np.genfromtxt(folder/"observables.csv",delimiter=",",names=True)
    csv=np.atleast_1d(csv)
    assert len(csv)==summary["cycles"]*summary["steps_per_cycle"]
    assert np.all(np.isfinite(np.column_stack([csv[k] for k in csv.dtype.names])))
    assert np.all(np.diff(csv["t"])>0)
    params=summary.get("parameters")
    if params is None:
        params=json.loads((ROOT/"host/dev_problem/d0.json").read_text())["parameters"]
    specifications=[json.loads(p.read_text()) for p in
        [ROOT/"host/dev_problem/d0.json"]+sorted((ROOT/"verification").glob("*problem.json"))]
    problem=next(p for p in specifications if p["id"]==summary["problem"])
    E,nu,beta,c0,kappa=(params[k] for k in ("E","nu","beta","c0","kappa"))
    mu=E/(2*(1+nu));lam=E*nu/((1+nu)*(1-2*nu))
    keys=["D_int","C_int","E_el","A_damage","d_max"]
    maximum=np.zeros(5);complement=0.;nfields=0;previous=None
    files=sorted(folder.glob("fields_cycle_*.npz"))
    assert len(files)==summary["cycles"]
    for file in files:
        z=np.load(file,allow_pickle=False)
        nodes,elems,U,sv=z["nodes"],z["elems"],z["U"],z["svars"]
        names=z["dof_names"].tolist();dpn=len(names)
        assert len(set(names))==dpn and all(k in names for k in ("ux","uy","c","d"))
        assert U.shape==(len(nodes)*dpn,)
        assert len(np.unique(nodes,axis=0))==len(nodes)
        assert len(np.unique(elems,axis=0))==len(elems)
        assert np.array_equal(np.unique(elems),np.arange(len(nodes)))
        assert elems.shape==(summary["n"]**2,4)
        assert sv.shape==(len(elems),12)
        assert all(np.all(np.isfinite(a)) for a in (nodes,U,sv))
        assert np.all(sv>=0)
        if previous is not None:
            assert np.all(sv[:,0::3]>=previous[:,0::3])
            assert np.all(sv[:,1::3]>=previous[:,1::3])
        previous=sv
        nodal=U.reshape(-1,dpn)
        if "s" in names:
            complement=max(complement,float(np.max(abs(nodal[:,names.index("d")]+nodal[:,names.index("s")]-1))))
        idx=[names.index(k) for k in ("ux","uy","c","d")]
        vals=nodal[:,idx][elems];X=nodes[elems]
        totals=np.zeros(3)
        for xi,eta in GP:
            N,dn=shape(xi,eta)
            J=np.einsum("mai,ar->mir",X,dn)
            det=np.linalg.det(J);assert np.all(det>0)
            B=np.einsum("ar,mri->mai",dn,np.linalg.inv(J))
            gu=np.einsum("mai,maj->mij",vals[:,:,:2],B)
            cg,dg=(np.einsum("ma,a->m",vals[:,:,j],N) for j in (2,3))
            ee=(gu+np.swapaxes(gu,1,2))/2-beta*(cg-c0)[:,None,None]*np.eye(2)
            trace=np.trace(ee,axis1=1,axis2=2)
            psi=mu*np.sum(ee*ee,axis=(1,2))+lam*trace**2/2
            energy=((1-kappa)*(1-dg)**2+kappa)*psi
            totals+=np.array([np.sum(dg*det),np.sum(cg*det),np.sum(energy*det)])
        area=0.
        sub=(np.arange(4)+.5)/2-1
        for xi in sub:
            for eta in sub:
                N,dn=shape(xi,eta)
                det=np.linalg.det(np.einsum("mai,ar->mir",X,dn))
                area+=float(np.sum(((vals[:,:,3]@N)>=.5)*det/4))
        dmax=nodal[:,names.index("d")].copy()
        if problem.get("notch"):
            prescribed=np.isclose(nodes[:,1],.5)&(nodes[:,0]<=problem["notch"]["x_end"]+1e-12)
            dmax[prescribed]=-np.inf
        totals=np.r_[totals,area,float(dmax.max())]
        cycle=int(file.stem.rsplit("_",1)[1]);target=cycle*problem["loading"]["T"]
        line=csv[np.argmin(abs(csv["t"]-target))]
        stored=np.array([line[k] for k in keys])
        maximum=np.maximum(maximum,abs(totals-stored))
        assert np.allclose(totals,stored,rtol=3e-10,atol=2e-12),(folder,cycle,totals,stored)
        nfields+=len(nodes)*dpn+sv.size
    assert complement<2e-8,complement
    return dict(path=str(folder.relative_to(ROOT)),status="passed",cycles=len(files),
        checked_scalar_values=nfields,observable_max_abs_errors=dict(zip(keys,maximum.tolist())),
        complement_max_error=complement)


if __name__=="__main__":
    contract=json.loads((ROOT/"element/uel_contract.json").read_text())
    expected=sorted(hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in contract["source"])
    rows=[];excluded=[];failed=[]
    for path in sorted((ROOT/"verification/runs").glob("*/run.json")):
        summary=json.loads(path.read_text())
        record=dict(path=str(path.parent.relative_to(ROOT)),status=summary["status"])
        if sorted(summary.get("source_sha256",{}).values())!=expected:
            record["reason"]="source differs from final submission; exploration evidence only"
            excluded.append(record)
        elif summary["status"]!="ok":
            failed.append(record)
        else:
            rows.append(audit(path.parent))
    result=dict(status="passed" if not failed else "failed",runs=rows,
        failed_final_source_runs=failed,excluded_predecessor_runs=excluded,
        scope="completed runs with the exact final source hashes; reference Qp uses separate observables")
    (ROOT/"verification/output_audit.json").write_text(json.dumps(result,indent=2)+"\n")
    host={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
          for p in sorted((ROOT/"host").rglob("*")) if p.is_file()}
    manifest=ROOT/"verification/host_manifest.json"
    if manifest.exists():
        assert json.loads(manifest.read_text())==host,"Host files differ from the earlier recorded manifest"
    else:
        manifest.write_text(json.dumps(host,indent=2)+"\n")
    print(json.dumps(result,indent=2))
