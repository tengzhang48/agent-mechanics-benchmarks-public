"""Build inspectable result tables from persisted runs with final-source hashes."""
from pathlib import Path
import hashlib
import json
import math
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
KEYS=('d_max','D_int','C_int','E_el')


def csv(path):
    return np.atleast_1d(np.genfromtxt(path/'observables.csv',delimiter=',',names=True))


def table(headers,rows):
    return ['| '+' | '.join(headers)+' |','| '+' | '.join(['---']*len(headers))+' |']+[
        '| '+' | '.join(str(v) for v in row)+' |' for row in rows]


def fmt(v):
    return f'{v:.7g}' if isinstance(v,(float,np.floating)) else str(v)


def discrepancy(coarse,fine):
    return {k:float(np.max(abs(coarse[k]-np.interp(coarse['t'],fine['t'],fine[k])))/
              max(np.max(abs(fine[k])),1e-30)) for k in KEYS}


def audit_reference(path,r):
    data=csv(path)
    assert len(data)==r['cycles']*r['spc'],path
    assert all(np.all(np.isfinite(data[k])) for k in data.dtype.names),path
    assert np.all(np.diff(data['t'])>0),path
    files=sorted(path.glob('fields_cycle_*.npz'));assert len(files)==r['cycles'],path
    previous=None;count=0
    ng=r.get('quadrature_order',r['order']+1)**2
    for j,file in enumerate(files):
        z=np.load(file,allow_pickle=False)
        assert z['elems'].shape==(r['n']**2,(r['order']+1)**2)
        assert z['U'].size==len(z['nodes'])*4
        state=z['svars'];assert state.shape==(r['n']**2,3*ng)
        assert np.all(state>=0) and all(np.all(np.isfinite(z[k])) for k in ('nodes','U','svars'))
        if previous is not None:
            assert np.all(state[:,0::3]>=previous[:,0::3])
            assert np.all(state[:,1::3]>=previous[:,1::3])
        previous=state
        assert abs(z['U'][3::4].max()-data['d_max'][(j+1)*r['spc']-1])<2e-12
        count+=z['U'].size+state.size
    return dict(status='passed',saved_cycles=len(files),finite_scalar_values=count,
                scope='Coverage, finiteness, history monotonicity and nodal-maximum output bridge; Qp integral formula is separately checked analytically')


def main():
    c=json.loads((ROOT/'element/uel_contract.json').read_text())
    expected=sorted(hashlib.sha256((ROOT/s).read_bytes()).hexdigest() for s in c['source'])
    runs=[];pending=[];other=[];failed=[]
    for path in sorted((ROOT/'verification/runs').iterdir()):
        if not path.is_dir():continue
        if not (path/'run.json').exists():
            if c['formulation_id'] in path.name:
                pending.append(dict(path=str(path.relative_to(ROOT)),
                    saved_cycle_files=len(list(path.glob('fields_cycle_*.npz')))))
            continue
        r=json.loads((path/'run.json').read_text())
        if sorted(r.get('source_sha256',{}).values())!=expected:
            other.append(str(path.relative_to(ROOT)));continue
        record=dict(path=str(path.relative_to(ROOT)),**r)
        if r['status']!='ok':failed.append(record);continue
        rows=csv(path);assert len(rows)==r['cycles']*r['steps_per_cycle']
        record['final_observables']={k:float(rows[k][-1]) for k in KEYS}
        record['cycle_observables']=[dict(cycle=j+1,**{k:float(rows[k][(j+1)*r['steps_per_cycle']-1]) for k in KEYS})
            for j in range(r['cycles'])]
        runs.append(record)
    manufactured=[]
    for mesh in ('structured','distorted'):
        prior=None
        for n in (8,16,32,64):
            r=json.loads((ROOT/f'verification/enhanced_manufactured/{mesh}_{n}_s1.json').read_text())
            rec=dict(r)
            if prior:
                rec['rates']={k:math.log(prior['L2'][k]/r['L2'][k],2) for k in ('ux','uy','c','d','p')}
                rec['rates']['flux']=math.log(prior['flux_L2']/r['flux_L2'],2)
            manufactured.append(rec);prior=r
    references=[]
    pending_references=[dict(path=str(p.relative_to(ROOT)),saved_cycle_files=len(list(p.glob('fields_cycle_*.npz'))))
        for p in sorted((ROOT/'verification/reference_runs').iterdir()) if p.is_dir() and not (p/'run.json').exists()]
    for path in sorted((ROOT/'verification/reference_runs').glob('*/run.json')):
        r=json.loads(path.read_text())
        rec=dict(path=str(path.parent.relative_to(ROOT)),**r)
        if r['status']=='ok':
            rec['output_coverage']=audit_reference(path.parent,r)
            rec['final_observables']={k:float(csv(path.parent)[k][-1]) for k in KEYS}
        references.append(rec)
    temporal=[]
    for window,spcs in ((5,(40,80,160,320)),(20,(40,80,160,320))):
        selected=[]
        for spc in spcs:
            matches=[r for r in runs if r['problem']=='D0' and r['mesh']=='structured' and r['n']==32
                     and r['steps_per_cycle']==spc and r['cycles']>=window]
            if not matches:continue
            r=min(matches,key=lambda r:r['cycles'])
            rows=csv(ROOT/r['path']);rows=rows[rows['t']<=8*window+1e-8]
            selected.append((r,rows))
        if len(selected)>1:
            fine,fr=selected[-1]
            for r,rows in selected[:-1]:
                temporal.append(dict(cycles=window,coarse_spc=r['steps_per_cycle'],fine_spc=fine['steps_per_cycle'],
                    normalized_max_history_difference=discrepancy(rows,fr)))
    comparisons=[];seen_comparisons=set()
    for p in sorted((ROOT/'verification/field_comparisons').glob('*.json')):
        rec=json.loads(p.read_text());details=rec['comparison']
        key=(details['production'],details['reference'],details['cycle'])
        if key in seen_comparisons:continue
        seen_comparisons.add(key);comparisons.append(rec)
    control_path=ROOT/'verification/enhancement_comparison.json'
    control=json.loads(control_path.read_text()) if control_path.exists() else None
    diffusion_path=ROOT/'verification/damaged_diffusion_resolution.json'
    damaged_diffusion=json.loads(diffusion_path.read_text()) if diffusion_path.exists() else None
    tests={}
    for filename in ('element_checks','rollback_checks','limit_checks','history_time_checks','energy_bridge_checks','accumulator_checks','mutation_checks',
                     'reference_checks','reference_checks_q6','regeneration','output_audit'):
        p=ROOT/'verification'/f'{filename}.json'
        if p.exists():tests[filename]=json.loads(p.read_text()).get('status')
    result=dict(status='failed' if failed or any(v!='passed' for v in tests.values()) else 'recorded',
        source_sha256=dict(zip(c['source'],[hashlib.sha256((ROOT/s).read_bytes()).hexdigest() for s in c['source']])),
        tests=tests,completed_final_source_runs=runs,pending_final_source_runs=pending,
        failed_final_source_runs=failed,excluded_predecessor_runs=other,
        manufactured=manufactured,independent_references=references,pending_independent_references=pending_references,
        time_refinement=temporal,field_comparisons=comparisons,matched_enhancement_comparison=control,
        damaged_diffusion_resolution=damaged_diffusion)
    (ROOT/'verification/results.json').write_text(json.dumps(result,indent=2)+'\n')
    lines=['# Measured results','',
        'Generated by `source env.sh; python verification/summarize.py`. Only runs whose four source hashes match the submission enter the primary tables. A completed solve is evidence of robustness; reference comparisons and refinement assess accuracy.','',
        f'Completed final-source runs: {len(runs)}. Failed final-source runs: {len(failed)}. Pending runs: {len(pending)}.','',
        '## Compiled and transactional checks','']
    lines+=table(['Check','Result'],tests.items())
    lines+=['','## Global runs','',
        'Values below are at the final cycle; `d_max` excludes prescribed notch nodes. Integrals and energy are computed by the unmodified host from the nodal Q1 fields.','']
    rows=[]
    for r in runs:
        v=r['final_observables'];stats=r.get('newton_linear_solves',{})
        rows.append([Path(r['path']).name,r['accepted_output_steps'],r.get('rejected_attempts',0),
                     f"{stats.get('min','—')}–{stats.get('max','—')}"]+[fmt(v[k]) for k in KEYS])
    lines+=table(['Run','Outputs','Rejected','Newton updates','d_max','D_int','C_int','E_el'],rows)
    if control:
        lines+=['','## Matched enhancement comparison','',
            'Structured n=16, twenty cycles, 40 steps/cycle. Properties, history split, reaction, guard and residual scales are identical. The generated core is byte-identical after removing only the enhancement hooks. This isolates the strain approximation; both solves finish.','']
        records=[('Ordinary Q1 strain',control['basic']),('Enhanced strain',control['enhanced'])]
        if 'independent_Q4_acceptedH_reference' in control:
            records.append(('Independent Q4, same history split',control['independent_Q4_acceptedH_reference']))
        lines+=table(['Variant']+list(KEYS),[[name]+[fmt(rec['final'][k]) for k in KEYS] for name,rec in records])
    lines+=['','## Manufactured full-PDE fields','',
        'History H=.7 is an inactive accepted prehistory; analytical continuum forcing is differentiated independently with SymPy. L2 norms use 5×5 error quadrature. No forcing is included in the submitted UEL.','']
    rows=[]
    for r in manufactured:
        rows.append([r['mesh'],r['n']]+[fmt(r['L2'][k]) for k in ('ux','c','d')]+[fmt(r['grad_p_L2']),fmt(r['flux_L2']),fmt(r.get('rates',{}).get('c','—'))])
    lines+=table(['Mesh','n','u_x L2','c L2','d L2','grad p L2','J L2','c rate'],rows)
    lines+=['','The distorted raw pressure-gradient error does not decrease strongly. Concentration is about first order with the selected enhanced strain; displacement and damage are about second order. The full flux decreases.','',
        '## Time refinement at n=32','',
        'Each number is max absolute observable difference over common physical output times, divided by the maximum absolute value in the finer run. This measures time sensitivity of the spatial discretization; it is not an error bound for the PDE.','']
    lines+=table(['Cycles','Coarse spc','Finer spc']+list(KEYS),[
        [r['cycles'],r['coarse_spc'],r['fine_spc']]+[fmt(r['normalized_max_history_difference'][k]) for k in KEYS] for r in temporal])
    if damaged_diffusion:
        lines+=['','## Low-mobility concentration resolution','',
            'A fully damaged beta=0 step-charging limit has diffusivity 1.4e-6 and diffusion length .003347 at t=8. Compiled host n=8/16 solutions agree with the independently assembled 1D reduction within 6e-15. The finer refinement below uses that 1D reduction, not full 2D host runs. A gamma-time erfc formula gives the continuous-space backward-Euler reference; its resolvent recurrence and integration order are verified independently.','',
            'Consistent-mass concentration undershoots the physical initial value c0=.35 on unresolved meshes. No positivity guarantee is claimed. This demonstrates a resolution mechanism; it does not attribute all coupled severe-case errors to it.','']
        lines+=table(['1D n','Final c_min','C_int','Exact BE C_int','Spatial L2 error'],[
            [r['n']]+[fmt(r[k]) for k in ('c_min_final','C_int','continuous_space_BE_C_int','L2_spatial_error')]
            for r in damaged_diffusion['spatial_study']])
    lines+=['','## Independent Qp runs','',
        'The independent NumPy implementation uses physical displacement Hessians and direct grad(g p_eff), no generator or projected pressure. Standard runs use current endpoint H; `_acceptedH` runs use the same first-order history split as the UEL but a consistent high-order reaction. Finite reference discretization and time errors remain.','']
    lines+=table(['Run','Status']+list(KEYS),[
        [Path(r['path']).name,r['status']]+[fmt(r.get('final_observables',{}).get(k,'—')) for k in KEYS] for r in references])
    if pending_references:
        lines+=['','Independent reference runs still pending:','']+[
            f"- `{r['path']}`: {r['saved_cycle_files']} saved cycles, no completed summary yet." for r in pending_references]
    if comparisons:
        matched=[]
        for rec in comparisons:
            r=rec.get('comparison',{})
            if Path(r.get('reference','')).name!='D0_Q4_structured_8_spc80_cy20_acceptedH':continue
            prod=json.loads((ROOT/r['production']/'run.json').read_text())
            matched.append((prod,r))
        if matched:
            matched.sort(key=lambda pair:(pair[0]['mesh'],pair[0]['n']))
            lines+=['','## Spatial comparisons at the same history split','',
                'Both UEL and independent Q4 reference use accepted H and dt=.1. This isolates spatial effects from the finite-step difference to endpoint-history references. The reference still has finite spatial/integration error. Relative observable differences are signed fractions.','']
            lines+=table(['Mesh n','Cycle','u_x L2','c L2','d L2','d_max relative','D_int relative','C_int relative','E_el relative'],[
                [f"{p['mesh']} {p['n']}",r['cycle']]+[fmt(r['primary_L2_errors'][k]) for k in ('ux','c','d')]+
                [fmt(r['observable_relative_errors'][k]) for k in KEYS] for p,r in matched])
        lines+=['','## Physical-space field comparisons','',
            'Qp interpolation is checked against physical quadratic polynomials on structured and distorted meshes. L2 errors integrate production Q1 fields against the independent interpolated reference. The enhanced-energy recovery is additional diagnostic evidence; the scored observable remains the host Q1 energy. The dt column gives UEL/reference step lengths. Differing steps or history schemes include temporal effects and do not isolate spatial error.','']
        rows=[]
        for rec in comparisons:
            if 'comparison' not in rec:continue
            r=rec['comparison']
            rows.append([Path(r['production']).name,Path(r['reference']).name,r['cycle'],
                f"{fmt(r['production_dt'])}/{fmt(r['reference_dt'])}"]+
                [fmt(r['primary_L2_errors'][k]) for k in ('ux','c','d')]+
                [fmt(r[k]) for k in ('host_Q1_energy','independently_recovered_enhanced_energy','reference_energy')]+
                [fmt(r.get('reference_sampled_on_production_nodes_Q1_energy','—'))])
        lines+=table(['UEL run','Reference','Cycle','dt UEL/ref','u_x L2','c L2','d L2','Host energy','Enhanced energy','Reference energy','Reference sampled Q1 energy'],rows)
    if pending:
        lines+=['','## Pending at table-generation time','']+[
            f"- `{r['path']}`: {r['saved_cycle_files']} saved cycles, no completed summary yet." for r in pending]
    if failed:
        lines+=['','## Failed final-source runs','']+[f"- `{r['path']}`: {r['status']}" for r in failed]
    lines+=['','Predecessor and failed experimental formulations remain under `exploration/` or are explicitly excluded by source hash. `results.json` retains the full run summaries, history minima, source hashes, time differences and exclusions. Hidden tests and Abaqus are unverified.','']
    (ROOT/'verification/results.md').write_text('\n'.join(lines))
    print(json.dumps(dict(completed=len(runs),failed=len(failed),pending=len(pending),tests=tests),indent=2))


if __name__=='__main__':main()
