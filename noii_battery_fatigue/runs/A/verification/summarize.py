"""Regenerate inspectable convergence tables from saved numerical artifacts."""
import json
from pathlib import Path
import numpy as np
from common import ROOT,write_json

OUT=ROOT/'verification'
metrics=['d_max','D_int','C_int','E_el']
data={}; text=['# Verification results','',
    'Generated from saved artifacts by `python verification/summarize.py`. '
    'Host observables use the host\'s own Q4 quadrature. '
    'Fine-grid comparisons are convergence studies, not an independent 2D answer key.','']

def table(headers,rows):
    text.append('| '+' | '.join(headers)+' |')
    text.append('| '+' | '.join(['---']*len(headers))+' |')
    for row in rows:
        text.append('| '+' | '.join(f'{v:.6g}' if isinstance(v,(float,np.floating)) else str(v) for v in row)+' |')
    text.append('')

text+=['## Compiled element and transaction checks','']
for filename in ['element_tests.json','reference_tests.json','rollback.json']:
    path=OUT/filename
    if path.exists():
        data[filename]=json.loads(path.read_text())
        text.append(f'- [{filename}]({filename})')
text+=['','The element tests include analytic forces on a nonaffine element, '
       'positive-trace and deviatoric energy, rigid motion, local mass balance, '
       'finite zero-time rejection, invalid orientation, and H=1e22. '
       'The global rejected-step comparison is bitwise identical.','']

text+=['## Manufactured 2D continuum solution','',
       'Analytical body sources are differentiated from the continuum laws '
       'using Sympy. H is a frozen spatial field. All physical fields have '
       'analytical Dirichlet values; q,z have no prescribed values.','']
mm=[]
for mesh in ['structured','distorted']:
    for n in [8,16,32,64]:
        path=OUT/'manufactured'/f'{mesh}_n{n}_lump1_p0.json'
        if path.exists():
            a=json.loads(path.read_text()); mm.append(a)
table(['mesh','n','L2 ux','L2 c','L2 d','L2 q','L2 grad q','L2 flux'],
      [[a['mesh'],a['n'],a['L2']['ux'],a['L2']['c'],a['L2']['d'],a['L2']['q'],
        a['grad_q_L2'],a['flux_L2']] for a in mm])
data['manufactured']=mm
text+=['The structured field errors approach order two. On randomly distorted '
       'meshes concentration/pressure convergence is slower. Isolated pressure '
       'gradient error does not decrease through n=64; total flux does. '
       'This limitation also appears with the tested mixed-volumetric-stress alternative.','']

text+=['## Independent coupled transient reference','',
       'The cell-centred FV solver analytically eliminates displacement under '
       'ux=0 and sigma_yy=0. It uses direct stress differences in the flux, '
       'and the original damage ODE, with no q or z. The UEL solves a narrow '
       '2D strip on the unchanged host. All integrals are divided by strip width. '
       'Both use five cycles; G_c=1e-4 and alpha_T=5e-4.','']
refdir=OUT/'reference_1d'
ref=np.load(refdir/'lag_fv_n512_spc80'/'result.npz')
reference=ref['history']; refscale=np.max(np.abs(reference[:,1:]),axis=0)
comparison=[]
for n in [16,32,64,128]:
    path=refdir/f'uel_n{n}_spc80'/'result.npz'
    if path.exists():
        a=np.load(path); hist=a['history']
        diff=np.max(np.abs(hist[:,1:]-reference[:,1:]),axis=0)
        field=[float(np.sqrt(np.mean((np.interp(ref['y'],a['y'],a[k])-ref[k])**2))) for k in ['c','d']]
        comparison.append({'n':n,'max_abs_observable_error':dict(zip(metrics,diff.tolist())),
                           'max_normalized_observable_error':dict(zip(metrics,(diff/refscale).tolist())),
                           'final_c_L2':field[0],'final_d_L2':field[1]})
table(['UEL rows','max abs dmax error','max abs D_int error','max abs C_int error',
       'max abs E_el error','final c L2','final d L2'],
      [[a['n'],*a['max_abs_observable_error'].values(),a['final_c_L2'],a['final_d_L2']] for a in comparison])
data['1d_spatial_against_FV512']=comparison
alpha=float(ref['svars'][:,1].max())
text+=[f'The FV reference reaches alpha={alpha:.8g}, above alpha_T. '
       'The two-dimensional tensile split is separately tested by the element patch.','']
fvc=[]
for method,ns in [('lag_fv',[128,256,512]),('uel',[16,32,64,128])]:
    for n in ns:
        path=refdir/f'{method}_n{n}_spc80'/'summary.json'
        if path.exists():
            a=json.loads(path.read_text()); fvc.append(a)
table(['method','n','final dmax','final D_int','final C_int','final E_el'],
      [[a['method'],a['n'],*a['final'][1:]] for a in fvc])
data['1d_spatial_final']=fvc
temporal=[]
for method,n in [('fv',256),('lag_fv',256),('uel',128)]:
    for spc in [40,80,160,320]:
        path=refdir/f'{method}_n{n}_spc{spc}'/'summary.json'
        if path.exists():
            a=json.loads(path.read_text()); temporal.append(a)
table(['method','steps/cycle','final dmax','final D_int','final C_int','final E_el'],
      [[a['method'],a['spc'],*a['final'][1:]] for a in temporal])
data['1d_temporal_final']=temporal
text+=['The spatial error table uses the same lagged-history split in the '
       'independent FV solver. The temporal table additionally includes the '
       'implicit-history FV solution: both schemes approach the same limit. '
       'Successive 80/160/320-step changes in the integrated observables '
       'roughly halve. A remaining finite spatial error is visible.','']
uniformpath=OUT/'uniform_history'/'summary.json'
if uniformpath.exists():
    uniform=json.loads(uniformpath.read_text()); data['uniform_history']=uniform
    text+=['## Uniform fatigue time benchmark','',
           'Uniform concentration is prescribed everywhere, ux is constrained '
           'and sigma_yy=0. This isolates constitutive/history/damage integration; '
           'it is not a diffusion solve. An independent scalar bracketed root '
           'solves the implicit damage/history law at 5120 steps/cycle, with '
           'a 2560-step reference self-check. There is no spatial discretization '
           'error in the homogeneous Q1 element.','']
    table(['steps/cycle','max abs d error','max abs H error','max abs alpha error','final d','final alpha'],
          [[a['spc'],a['max_abs_error']['d'],a['max_abs_error']['H'],a['max_abs_error']['alpha'],
            a['final'][1],a['final'][3]] for a in uniform['runs']])
    text+=['The full error histories and reference resolution check are in '
           '`uniform_history/summary.json`. Early damage peaks can require many '
           'steps to resolve even though Newton itself converges.','']

text+=['## D0 on the original host CLI','',
       'Final formulation: lumped damage source, direct stress, 3x3 history '
       'points. All listed runs use the original host/run_problem.py. '
       'Each full run is 20 cycles unless explicitly marked five.','']
names=[('structured',8,40,'selected_D0_s8_spc40'),
       ('structured',16,40,'selected_D0_s16_spc40'),
       ('structured',32,40,'selected_D0_s32_spc40'),
       ('structured',64,40,'selected_D0_s64_spc40_5cyc'),
       ('distorted',8,40,'selected_D0_d8_spc40'),
       ('distorted',16,40,'selected_D0_d16_spc40'),
       ('distorted',32,40,'selected_D0_d32_spc40'),
       ('structured',16,80,'selected_D0_s16_spc80'),
       ('structured',16,160,'selected_D0_s16_spc160')]
names.extend([('distorted',16,80,'selected_D0_d16_spc80'),
              ('distorted',16,160,'selected_D0_d16_spc160')])
runs=[]; pending=[]
for mesh,n,spc,name in names:
    path=OUT/'runs'/name
    if not (path/'run.json').exists():
        pending.append(name); continue
    a=json.loads((path/'run.json').read_text())
    hist=np.genfromtxt(path/'observables.csv',delimiter=',',names=True)
    row={'name':name,**a,'cycle5':None,'cycle20':None}
    for cycle in [5,20]:
        if len(hist)>=cycle*spc:
            row[f'cycle{cycle}']={k:float(hist[k][cycle*spc-1]) for k in metrics}
    runs.append(row)
table(['mesh','n','steps/cycle','cycles reached','status','cutbacks'],
      [[r['mesh'],r['n'],r['steps_per_cycle'],r['accepted_output_steps']/r['steps_per_cycle'],
        r['status'],r.get('cutbacks','?')] for r in runs])
text+=['At the end of cycle five:','']
table(['mesh','n','steps/cycle',*metrics],
      [[r['mesh'],r['n'],r['steps_per_cycle'],*r['cycle5'].values()] for r in runs if r['cycle5']])
text+=['At the end of cycle twenty:','']
table(['mesh','n','steps/cycle',*metrics],
      [[r['mesh'],r['n'],r['steps_per_cycle'],*r['cycle20'].values()] for r in runs if r['cycle20']])
data['D0_runs']=runs; data['pending_D0_runs']=pending
if pending:
    text+=['Pending when this table was generated: '+', '.join(pending)+'.','']
text+=['The smallest damage-area indicator A_damage can be zero throughout a '
       'fine run while coarse meshes exceed d=0.5. That threshold observable '
       'is particularly mesh sensitive; the full host CSVs retain it. '
       'These tables do not establish an independent two-dimensional reference '
       'error or spatial asymptotic rate for the long-cycle problem.','']
text+=['The earlier `final_D0_*` and `D0_lumped_s8_spc40` directories belong '
       'to the fully implicit candidate before the final pivot. Their names '
       'are retained for audit continuity; they are not the selected time '
       'integration. Their archived source hash is fcd4c37a1a12.','']
fieldpath=OUT/'field_convergence.json'
if fieldpath.exists():
    fields=json.loads(fieldpath.read_text()); data['D0_field_convergence']=fields
    text+=['Q1 field differences at common physical quadrature points, relative '
           'to the selected structured n=64 solution at cycle five:','']
    table(['mesh','n','L2 u','L2 c','L2 d','L2 q'],
          [[a['mesh'],a['n'],*a['L2_difference_from_n64'].values()] for a in fields['spatial_cycle_5']])
    text+=['Field changes at cycle twenty for time refinement on n=16:','']
    table(['mesh','coarse/fine steps','L2 u','L2 c','L2 d','L2 q','max delta H','max delta alpha'],
          [[a['mesh'],str(a['coarse_spc'])+'/'+str(a['fine_spc']),*a['L2_change'].values(),
            a['H_max_change'],a['alpha_max_change']] for a in fields['temporal_cycle_20']])

text+=['## Solver robustness and known failures','',
       'Newton work is counted as sparse linear solves in the unchanged host '
       'one-block Newton path. Cutbacks use its ordinary transaction and limits.','']
checks=[]
for path in sorted((OUT/'global_checks').glob('*/summary.json')):
    a=json.loads(path.read_text()); checks.append(a)
selected_hash=__import__('hashlib').sha256((ROOT/'element'/'uel.f90').read_bytes()).hexdigest()
table(['case','history','reaction','status','max free d','max H','cutbacks/rejections','mean Newton updates','max updates'],
      [[a['name'],'lagged' if a['source_sha256']==selected_hash else 'implicit',
        'lumped' if a['contract']['jprops'] and a['contract']['jprops'][0]==1 else 'full',
        'ok' if a['status']=='ok' else 'FAILED',a.get('d_max_free',a['d_max']),a['H_max'],
        a.get('cutbacks',a['rejected_attempts']),a['mean_linear_solves_per_attempt'],
        a['max_linear_solves_per_attempt']] for a in checks])
data['global_checks']=checks
text+=['Failed cases are retained with their accepted partial trajectories '
       'and per-attempt records. The two old notch failures predate the addition '
       'of a free-damage diagnostic, so their tabulated maximum may include '
       'the prescribed notch d=1. They do not demonstrate successful propagation. '
       'The free-damage toughness sweep has no prescribed damage.','',
       'The implicit extreme fatigue-floor case fails in its first nominal '
       'increment. The selected split completes that same physical problem '
       'with H near 1e20 and free d=1. A separate fixed-H=1e22 element call '
       'also passes.','',
       'All accepted H and alpha values in the reported runs are nondecreasing. '
       'Small nodal d decreases can remain at global-solve accuracy (the old '
       'implicit D0 statistics run records about 2.8e-7). No separate nodal irreversibility '
       'constraint or post-convergence clipping was imposed.','',
       'Actual Abaqus execution, hidden problems, arbitrary quadrilateral '
       'families, and universal Newton convergence are unverified.','']

q2runs=[]
for path in sorted((OUT/'q2_reference').glob('n*_spc*/summary.json')):
    q2runs.append(json.loads(path.read_text()))
q2runs.sort(key=lambda a:(a['n'],a['spc']))
if q2runs:
    data['Q2_reference_runs']=q2runs
    text+=['## Separate higher-order 2D comparison','',
           'This self-authored reference uses affine biquadratic displacement, '
           'concentration and damage, 5x5 quadrature, direct displacement '
           'Hessians in grad(g*p), and the original consistent damage equation. '
           'It has no projected-pressure or damage-source auxiliary fields. '
           'It uses the same lagged history in time to isolate spatial differences. '
           'Its nine-node macroelements are reference-only; the delivered element '
           'retains the required four-node topology. The generic host is unchanged.','']
    table(['Q2 macro n','steps/cycle','status','final dmax','final D_int','final C_int','final E_el'],
          [[a['n'],a['spc'],a['status'],*a['final'][1:]] for a in q2runs if a['final']])
    q2test=OUT/'q2_reference'/'element_tests.json'
    if q2test.exists():
        qa=json.loads(q2test.read_text());data['Q2_reference_element_tests']=qa
        text+=[f"The separate Q2 tangent has relative FD error {qa['Jacobian_relative_error']:.3g}; "
               f"the quadratic-displacement stress-flux residual error is {qa['quadratic_displacement_flux_residual_max_error']:.3g}.",'']
q2path=OUT/'q2_comparison.json'
if q2path.exists():
    q2=json.loads(q2path.read_text());data['Q2_comparison']=q2
    text+=['Cycle-five primary fields compared at common physical quadrature '
           'points against Q2 n=32:','']
    table(['Q1 mesh','n','L2 u','L2 c','L2 d','L2 flux','L2 grad q'],
          [[a['mesh'],a['n'],*a['L2_field_difference'].values(),a['flux_L2_difference'],
            a['grad_q_L2_difference']] for a in q2['Q1_comparison']])
    text+=['Maximum observable differences over the first five cycles:','']
    table(['Q1 mesh','n','dmax','D_int','C_int','E_el'],
          [[a['mesh'],a['n'],*a['max_abs_observable_difference'].values()] for a in q2['Q1_comparison']])
    text+=['Q2 reference family self-convergence against its n=32 mesh:','']
    table(['Q2 macro n','L2 u','L2 c','L2 d','max dmax diff','max D_int diff','max C_int diff','max E_el diff'],
          [[a['n'],*a['L2_field_difference'].values(),*a['max_abs_observable_difference'].values()] for a in q2['Q2_self_convergence']])
    text+=['This is an independently discretized numerical comparison from '
           'the stated equations, with its own finite reference error. It is '
           'not the organizer\'s hidden reference or an external Abaqus result.','']

write_json(OUT/'results.json',data)
(OUT/'RESULTS.md').write_text('\n'.join(text)+'\n')
print('Wrote verification/results.json and verification/RESULTS.md')
print('Pending D0 runs:',pending)
