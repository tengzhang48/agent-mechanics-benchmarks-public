"""Host-level rollback, large-step, saturation and Newton-work tests.

The same one-block Newton/line search and transaction as fehost.solve are used;
the linear solve is counted and failures can be injected before/after UEL calls.
"""
import argparse
import hashlib
import json
import time

import numpy as np
from scipy.sparse.linalg import spsolve

from common import ROOT, fehost, meshes, load_element, write_json
from hostfe import CompiledElementState, solve_steps_staggered
from hostfe.bounds import IncrementTooLarge
from observables import observables
from run_problem import edge_nodes


def run(name,n=8,spc=40,cycles=2,distorted=False,problem=None,force_reject=None,lag=False,lump=1):
    prob=json.loads((ROOT/'host'/'dev_problem'/'d0.json').read_text()) if problem is None else problem
    p,load=prob['parameters'],prob['loading']
    dt=load['T']/spc
    nodes,elems=meshes.distorted(n,seed=3) if distorted else meshes.structured(n)
    variant={'source':'verification/alternatives/uel_lagged.f90'} if lag else None
    el,contract=load_element(len(elems),dt,p,variant=variant)
    el._jprops=np.array([lump],dtype=np.int32)
    contract['jprops']=[lump]
    dpn=el.dof_per_node
    cnodes=edge_nodes(nodes,prob['c_edges'])
    notch=meshes.notch_nodes(nodes,prob['notch']['x_end']) if prob.get('notch') else np.array([],int)
    constraints=[(meshes.node_at(nodes,x,y),comps) for x,y,comps in prob['constraints']]
    U0=np.zeros(len(nodes)*dpn); U0[2::dpn]=p['c0']; U0[dpn*notch+3]=1.
    def bc(frac):
        t=frac*cycles*load['T']
        cb=p['c0']+load['A']*(1-np.cos(2*np.pi*t/load['T']))/2
        result={dpn*a+2:cb for a in cnodes}
        result.update({dpn*a+3:1. for a in notch})
        for a,comps in constraints:
            for comp in comps:
                result[dpn*a+el.comp(comp)]=0.
        return result
    rows=[]; commits=[]; attempts=[]
    latest={'U':U0.copy()}
    free_damage=np.ones(len(nodes),bool); free_damage[notch]=False
    work={'solves':0,'calls':0,'injected':False}
    previous={'H':el.svars[:,0::3].copy(),'alpha':el.svars[:,1::3].copy(),'d':U0[3::dpn].copy()}
    minima={'H':0.,'alpha':0.,'d':0.}
    def batch(coords,u,du):
        work['calls']+=1
        if force_reject and not work['injected'] and abs(el._time[1]-.8)<1e-10:
            work['injected']=True
            if force_reject=='after':
                # Force the UEL to see nonconverged, high-history trial states.
                trial=u.copy()
                trial[:,0::dpn]+=.02*coords[:,:,0]
                trial[:,2::dpn]+=.1
                el.element_rk_batch(coords,trial,du+(trial-u))
                el.element_rk_batch(coords,u,du)
            raise IncrementTooLarge('verification: forced rejection '+force_reject)
        return el.element_rk_batch(coords,u,du)
    def linear(k,b):
        work['solves']+=1
        return spsolve(k,b)
    def attempt(rec):
        attempts.append({**vars(rec),'linear_solves':work['solves'],'batch_calls':work['calls']})
        work['solves']=0; work['calls']=0
    def commit(nodes,elems,dpn,U,old):
        el.commit_step(nodes,elems,dpn,U,old)
        latest['U']=U.copy()
        current={'H':el.svars[:,0::3],'alpha':el.svars[:,1::3],'d':U[3::dpn]}
        for k,now in current.items():
            minima[k]=min(minima[k],float((now-previous[k]).min()))
            previous[k]=now.copy()
        commits.append({'time':float(el._time[1]),'dt':float(el.dt),
                        'd_min':float(U[3::dpn].min()),'d_max':float(U[3::dpn].max()),
                        'd_max_free':float(U[3::dpn][free_damage].max()),
                        'H_max':float(current['H'].max()),'alpha_max':float(current['alpha'].max())})
    def output(step,frac,U,metrics=None):
        rows.append({'step':step,'t':frac*cycles*load['T'],
                     'd_max':float(U[3::dpn].max()),**observables(nodes,elems,U,dpn,
                     {k:el.comp(k) for k in fehost.REQUIRED_DOFS},p)})
    start=time.time(); status='ok'; report=None; U=U0
    try:
        U,report=solve_steps_staggered(nodes,elems,dpn,el.element_rk,bc,
            [list(range(dpn))],cycles*spc,dt=dt,state=CompiledElementState([el]),
            batch_fn=batch,commit_fn=commit,U0=U0,tol=1e-8,max_field_newton=25,
            max_cutbacks=10,cutback_on_element_rejection=True,
            residual_scale=fehost.residual_weights(el,contract),linear_solve=linear,
            attempt_callback=attempt,step_callback=output,return_report=True)
    except Exception as exc:
        status=f'{type(exc).__name__}: {exc}'
        U=latest['U']
    result={'name':name,'status':status,'n':n,'spc':spc,'cycles':cycles,
            'distorted':distorted,'force_reject':force_reject,'injected':work['injected'],
            'wall_s':round(time.time()-start,3),'history_min_change':minima,
            'max_linear_solves_per_attempt':max((r['linear_solves'] for r in attempts),default=0),
            'mean_linear_solves_per_attempt':float(np.mean([r['linear_solves'] for r in attempts])) if attempts else 0.,
            'rejected_attempts':sum(not r['converged'] for r in attempts),
            'd_max':max((r['d_max'] for r in commits),default=0.),
            'd_max_free':max((r['d_max_free'] for r in commits),default=0.),
            'd_min':min((r['d_min'] for r in commits),default=0.),
            'H_max':max((r['H_max'] for r in commits),default=0.),
            'alpha_max':max((r['alpha_max'] for r in commits),default=0.),
            'accepted_subincrements':len(commits),
            'contract':contract,'problem':prob,
            'source_sha256':hashlib.sha256((ROOT/contract['source']).read_bytes()).hexdigest()}
    if report is not None:
        result.update(report.summary())
    dest=ROOT/'verification'/'global_checks'/name
    write_json(dest/'summary.json',result)
    # Failed attempts have infinity norm by design; serialize those as null.
    clean=[{k:(None if isinstance(v,float) and not np.isfinite(v) else v) for k,v in a.items()} for a in attempts]
    write_json(dest/'attempts.json',clean)
    write_json(dest/'commits.json',commits)
    write_json(dest/'observables.json',rows)
    np.savez(dest/'final.npz',nodes=nodes,elems=elems,U=U,svars=el.svars)
    print(json.dumps({k:v for k,v in result.items() if k not in ['contract','problem']}),flush=True)
    return U,el.svars.copy(),rows,result


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('mode',choices=['rollback','large','saturation','extreme','statistics'])
    a=ap.parse_args()
    if a.mode=='rollback':
        before=run('reject_before',force_reject='before')
        after=run('reject_after',force_reject='after')
        equal=all(np.array_equal(a,b) for a,b in zip(before[:2],after[:2]))
        assert equal and before[2]==after[2]
        assert before[3]['injected'] and after[3]['injected']
        write_json(ROOT/'verification'/'rollback.json',{'final_U_bitwise_equal':True,
                   'final_SVARS_bitwise_equal':True,'observables_identical':True})
    elif a.mode=='large':
        run('large_step_10',n=16,spc=10,cycles=20,distorted=True)
    elif a.mode=='statistics':
        run('D0_statistics',n=16,spc=40,cycles=20)
    else:
        prob=json.loads((ROOT/'host'/'dev_problem'/'d0.json').read_text())
        extreme=a.mode=='extreme'
        prob['id']='saturated_fatigue' if extreme else 'constrained_notch'
        prob['parameters'].update(G_c=1e-8 if extreme else 1e-9,alpha_T=1e-15 if extreme else 1e-4)
        prob['loading']['A']=.65
        prob['c_edges']=['left','top']
        prob['notch']={'x_end':.125}
        prob['constraints']=[[0.,0.,['ux','uy']],[1.,0.,['ux','uy']],
                             [1.,1.,['ux','uy']],[0.,1.,['ux','uy']]]
        write_json(ROOT/'verification'/('extreme_problem.json' if extreme else 'saturation_problem.json'),prob)
        run(prob['id'],n=8,spc=40,cycles=5,distorted=True,problem=prob)


if __name__=='__main__':
    main()
