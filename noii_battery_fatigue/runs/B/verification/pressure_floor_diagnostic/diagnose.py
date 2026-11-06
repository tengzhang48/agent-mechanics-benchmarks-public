"""Frozen-state pressure recovery diagnostic. No constitutive/path solve or edits.

Independent Q1 quadrature assembles two global consistent L2 projections on
exactly the stored u,c,d,v fields. Neither recovery is an exact gradient oracle.
"""
from pathlib import Path
import csv,datetime,hashlib,json,platform,shutil
import numpy as np
import scipy
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import splu

ROOT=Path(__file__).resolve().parents[2]
CASE=ROOT/'verification/custom_runs/selected/spatial_floor_distorted16_s3_c2_spc80'
OUT=ROOT/'verification/pressure_floor_diagnostic'
CORNERS=np.array([[-1.,-1.],[1.,-1.],[1.,1.],[-1.,1.]])

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def quad(X):
 for xi,eta in CORNERS/np.sqrt(3):
  s,t=CORNERS.T;N=(1+s*xi)*(1+t*eta)/4
  DN=np.array([s*(1+t*eta),t*(1+s*xi)]).T/4
  J=X.T@DN;B=DN@np.linalg.inv(J)
  yield N,B,float(np.linalg.det(J)),N@X

def globalize(elems,Ue,n):
 values=Ue.reshape(len(elems),4,6);U=np.empty((n,6));seen=np.zeros(n,bool);spread=0.
 for en,ue in zip(elems,values):
  for i,u in zip(en,ue):
   if seen[i]:spread=max(spread,float(np.max(abs(U[i]-u))))
   else:U[i]=u;seen[i]=True
 assert seen.all()
 return U,spread

def weighted(values,w,mask=None):
 x=np.asarray(values);m=np.ones(len(w),bool) if mask is None else np.asarray(mask,bool)
 if not m.any():return dict(count=0,area=0.)
 ww=w[m];xx=x[m];norm=ww.sum();order=np.argsort(xx);cw=np.cumsum(ww[order])/norm
 return dict(count=int(m.sum()),area=float(norm),min=float(xx.min()),max=float(xx.max()),mean=float(ww@xx/norm),rms=float(np.sqrt(ww@(xx*xx)/norm)),q50=float(np.interp(.5,cw,xx[order])),q90=float(np.interp(.9,cw,xx[order])),q99=float(np.interp(.99,cw,xx[order])))

def main():
 OUT.mkdir(parents=True,exist_ok=True)
 data=np.load(CASE/'failed_trial.npz');old=np.load(CASE/'last_state.npz');run=json.loads((CASE/'run.json').read_text());prov=json.loads((CASE/'provenance.json').read_text());pr=run['problem']['parameters']
 nodes,elems=data['nodes'],data['elems'];U,spread=globalize(elems,data['Ue'],len(nodes));Up,oldspread=globalize(elems,data['Ue']-data['DUe'],len(nodes));sv=data['svars'].reshape(-1,4,4);dt=float(data['dt'])
 E,nu,beta,c0,kap,D,chi=[pr[x] for x in ['E','nu','beta','c0','kappa','D','chi']];a=1-kap;mu=E/(2*(1+nu));lam=E*nu/((1+nu)*(1-2*nu))
 ii=[];jj=[];mm=[];bp=np.zeros(len(nodes));bg=np.zeros(len(nodes));points=[];rq_local=np.zeros((len(elems),4));rc_local=np.zeros_like(rq_local)
 for ei,en in enumerate(elems):
  uu=U[en];prev=Up[en]
  for gp,(N,B,w,xy) in enumerate(quad(nodes[en])):
   assert w>0
   grad_u=uu[:,:2].T@B;c,d,q,v=N@uu[:,2:];strain=(grad_u+grad_u.T)/2-beta*(c-c0)*np.eye(2);stress_eff=2*mu*strain+lam*np.trace(strain)*np.eye(2);p=float(np.trace(stress_eff)/3);g=float(kap+a*v*v);gc=uu[:,2]@B;gv=uu[:,5]@B;gg=2*a*v*gv;gq=uu[:,4]@B
   J=g*D*gc-chi*c*gq
   rq_local[ei]+=N*(q-g*p)*w;rc_local[ei]+=(N*(c-N@prev[:,2])/dt+B@J)*w
   M=np.outer(N,N)*w;ii.extend(np.repeat(en,4));jj.extend(np.tile(en,4));mm.extend(M.ravel());np.add.at(bp,en,N*p*w);np.add.at(bg,en,N*g*p*w)
   points.append(dict(ei=ei,gp=gp,N=N,B=B,w=w,xy=xy,c=float(c),d=float(d),v=float(v),g=g,p=p,q=float(q),grad_c=gc,grad_v=gv,grad_g=gg,grad_q_trial=gq,history_old=sv[ei,gp].copy()))
 M=coo_matrix((mm,(ii,jj)),shape=(len(nodes),len(nodes))).tocsr();factor=splu(M.tocsc());qe,qd=factor.solve(np.column_stack([bp,bg])).T
 qtrial=U[:,4];rq=M@qtrial-bg;rd=M@qd-bg;re=M@qe-bp
 # All following comparisons hold u,c,d,v and old state fixed. Only pressure recovery varies.
 rc={'trial':np.zeros(len(nodes)),'projected_degraded':np.zeros(len(nodes)),'effective_then_product':np.zeros(len(nodes))};gp_values=[]
 for x in points:
  ei,gp,N,B,w=x['ei'],x['gp'],x['N'],x['B'],x['w'];en=elems[ei];qed=float(N@qe[en]);qdd=float(N@qd[en]);ge=qe[en]@B;gd=qd[en]@B;ga=x['g']*ge+qed*x['grad_g'];diff=x['g']*D*x['grad_c'];drifts={key:-chi*x['c']*gg for key,gg in [('trial',x['grad_q_trial']),('projected_degraded',gd),('effective_then_product',ga)]};storage=N*(x['c']-N@Up[en,2])/dt
  for key in rc:np.add.at(rc[key],en,(storage+B@(diff+drifts[key]))*w)
  gp_values.append(dict(element=ei,gp=gp,xy=x['xy'].tolist(),weight=w,c=x['c'],d=x['d'],vfrac=x['v'],g=x['g'],p_eff=x['p'],q_trial=x['q'],q_degraded=qdd,q_eff=qed,product_g_q_eff=x['g']*qed,q_trial_over_g=x['q']/x['g'],q_degraded_over_g=qdd/x['g'],grad_c=x['grad_c'].tolist(),grad_g=x['grad_g'].tolist(),grad_q_trial=x['grad_q_trial'].tolist(),grad_q_degraded=gd.tolist(),grad_q_eff=ge.tolist(),grad_g_q_eff=ga.tolist(),diffusive_flux=diff.tolist(),drift_trial=drifts['trial'].tolist(),drift_projected_degraded=drifts['projected_degraded'].tolist(),drift_effective_then_product=drifts['effective_then_product'].tolist(),flux_trial=(diff+drifts['trial']).tolist(),flux_projected_degraded=(diff+drifts['projected_degraded']).tolist(),flux_effective_then_product=(diff+drifts['effective_then_product']).tolist(),history_old_H_alpha_W_psi=x['history_old'].tolist()))
 get=lambda k:np.asarray([x[k] for x in gp_values]);w=get('weight');g=get('g');p=get('p_eff');qt=get('q_trial');qD=get('q_degraded');qE=get('q_eff');Gq=get('grad_q_trial');GD=get('grad_q_degraded');GA=get('grad_g_q_eff');JT=get('drift_trial');JD=get('drift_projected_degraded');JA=get('drift_effective_then_product');diff=get('diffusive_flux');area=w.sum();norm=lambda x:float(np.sqrt(np.sum(w*np.sum(x*x,axis=1)))) if np.ndim(x)==2 else float(np.sqrt(w@(x*x)));mag=lambda x:np.sqrt(np.sum(x*x,axis=1));p_rms=float(np.sqrt(w@(p*p)/area));p_mask=abs(p)>=.01*p_rms
 peq=qD/g-p;deltaG=GD-GA;deltaJ=JD-JA;deltaTrialG=Gq-GD
 metrics={
  'pressure_q_over_g_error_squared':peq*peq,
  'potential_gradient_difference_squared':np.sum(deltaG*deltaG,axis=1),
  'species_drift_difference_squared':np.sum(deltaJ*deltaJ,axis=1),
  'current_projected_drift_squared':np.sum(JD*JD,axis=1),
  'alternative_product_drift_squared':np.sum(JA*JA,axis=1)}
 totals={k:float(w@x) for k,x in metrics.items()}
 def group(mask,label):
  row=dict(label=label,g=weighted(g,w,mask),p_eff=weighted(p,w,mask),q_over_g=weighted(qD/g,w,mask),q_over_g_minus_p=weighted(peq,w,mask),q_eff_minus_p=weighted(qE-p,w,mask),grad_q_degraded_norm=weighted(mag(GD),w,mask),grad_product_norm=weighted(mag(GA),w,mask),gradient_difference_norm=weighted(mag(deltaG),w,mask),drift_degraded_norm=weighted(mag(JD),w,mask),drift_product_norm=weighted(mag(JA),w,mask),qD_minus_raw_degraded_target=weighted(qD-g*p,w,mask),product_minus_raw_degraded_target=weighted(g*(qE-p),w,mask),g_grad_qE_term_norm=weighted(mag(g[:,None]*get('grad_q_eff')),w,mask),qE_grad_g_term_norm=weighted(mag(qE[:,None]*get('grad_g')),w,mask),diffusive_flux_norm=weighted(mag(diff),w,mask))
  area_fraction=float(w[mask].sum()/area);row['area_fraction']=area_fraction
  diff_rms=row['diffusive_flux_norm'].get('rms',0.);row['current_drift_over_diffusive_rms']=row['drift_degraded_norm'].get('rms',0.)/diff_rms if diff_rms else None;row['alternative_drift_over_diffusive_rms']=row['drift_product_norm'].get('rms',0.)/diff_rms if diff_rms else None
  row['squared_norm_fractions']={k:float(w[mask]@v[mask]/totals[k]) if totals[k]>0 else 0. for k,v in metrics.items()};row['error_enrichment_per_area']={k:v/area_fraction if area_fraction>0 else None for k,v in row['squared_norm_fractions'].items()}
  meaningful=mask&p_mask;ratios=np.abs(qD/g)/np.maximum(abs(p),np.finfo(float).tiny);row['amplification_abs_qD_over_g_p_where_p_meaningful']=weighted(ratios,w,meaningful);row['sign_disagreement_area_fraction_where_p_meaningful']=float(w[meaningful&(np.sign(qD/g)!=np.sign(p))].sum()/w[meaningful].sum()) if meaningful.any() else None
  return row
 thresholds=[kap,2*kap,10*kap,100*kap,.01,.1,np.inf];bins=[];lower=-np.inf
 for upper in thresholds:
  mask=(g>lower)&(g<=upper);bins.append(group(mask,f'{lower:g}<g<={upper:g}'));lower=upper
 cumul=[group(g<=cut,f'g<={cut:g}') for cut in [kap*(1+1e-8),2*kap,10*kap,100*kap,.01,.1]]
 # Match the stored element residuals, before any Dirichlet elimination.
 rq_ind=np.zeros(len(nodes));np.add.at(rq_ind,elems.ravel(),rq_local.ravel());RcStored=np.zeros(len(nodes));np.add.at(RcStored,elems.ravel(),data['Re'][:,2::6].ravel());free_c=~(np.isclose(nodes[:,0],0)|np.isclose(nodes[:,0],1))
 verify=dict(shared_nodal_trial_max_spread=spread,shared_nodal_previous_max_spread=oldspread,previous_U_vs_last_accepted_max_abs=float(np.max(abs(Up.ravel()-old['U']))),old_svars_equal_last_accepted=bool(np.array_equal(data['svars'],old['svars'])),stored_q_element_R_max_abs_difference=float(np.max(abs(rq_local-data['Re'][:,4::6]))),stored_c_element_R_max_abs_difference=float(np.max(abs(rc_local-data['Re'][:,2::6]))),assembled_q_R_vs_Mq_minus_b_max_abs=float(np.max(abs(rq_ind-rq))),pressure_projection_q_trial_residual_l2=float(np.linalg.norm(rq)),pressure_projection_q_trial_residual_inf=float(np.max(abs(rq))),exact_qD_projection_residual_inf=float(np.max(abs(rd))),exact_qE_projection_residual_inf=float(np.max(abs(re))),mass_matrix_symmetry_error_inf=float(np.max(abs((M-M.T).data))) if (M-M.T).nnz else 0.,minimum_mass_diagonal=float(M.diagonal().min()))
 projection=dict(raw_degraded_pressure_target_l2=norm(g*p),qD_minus_raw_degraded_target_l2=norm(qD-g*p),product_minus_raw_degraded_target_l2=norm(g*(qE-p)),trial_minus_qD_over_g_l2=norm((qt-qD)/g),g_grad_qE_term_l2=norm(g[:,None]*get('grad_q_eff')),qE_grad_g_term_l2=norm(qE[:,None]*get('grad_g')),q_trial_minus_qD_nodal_max_abs=float(np.max(abs(qtrial-qd))),q_trial_minus_qD_mass_norm=float(np.sqrt((qtrial-qd)@(M@(qtrial-qd)))),qD_mass_norm=float(np.sqrt(qd@(M@qd))),trial_projection_gradient_error_l2=norm(deltaTrialG),qD_gradient_l2=norm(GD),trial_vs_qD_gradient_relative_l2=norm(deltaTrialG)/norm(GD),qD_over_g_minus_p_l2=norm(peq),p_eff_l2=norm(p),q_eff_minus_p_l2=norm(qE-p),qD_minus_g_q_eff_l2=norm(qD-g*qE),gradient_difference_l2=norm(deltaG),gradient_product_l2=norm(GA),gradient_difference_relative_to_qD_l2=norm(deltaG)/norm(GD),gradient_difference_relative_to_product_l2=norm(deltaG)/norm(GA),species_drift_trial_l2=norm(JT),species_drift_qD_l2=norm(JD),species_drift_product_l2=norm(JA),species_diffusive_flux_l2=norm(diff),species_drift_difference_l2=norm(deltaJ),species_drift_difference_relative_to_qD_l2=norm(deltaJ)/norm(JD))
 species={k:dict(free_residual_l2=float(np.linalg.norm(v[free_c])),free_residual_inf=float(np.max(abs(v[free_c]))),all_node_residual_sum=float(v.sum())) for k,v in rc.items()};species['current_trial_vs_stored_free_R_l2_difference']=float(np.linalg.norm((rc['trial']-RcStored)[free_c]))
 source={str(p.relative_to(ROOT)):sha(p) for p in [ROOT/'element/generated_model.for',ROOT/'element/uel.f90',ROOT/'element/model.py',ROOT/'element/uel_contract.json']};inputs={str((CASE/x).relative_to(ROOT)):sha(CASE/x) for x in ['failed_trial.npz','last_state.npz','run.json','provenance.json']}
 assert inputs[str((CASE/'failed_trial.npz').relative_to(ROOT))]==prov['artifact_sha256']['failed_trial.npz']
 exact_source=all(source[str(ROOTp)]==digest for ROOTp,digest in [(Path('element/generated_model.for'),prov['source_sha256']['element/generated_model.for']),(Path('element/uel.f90'),prov['source_sha256']['element/uel.f90']),(Path('element/uel_contract.json'),prov['source_sha256']['element/uel_contract.json'])])
 for i,x in enumerate(gp_values):
  x.update(qD_over_g_minus_p_eff=float(peq[i]),gradient_choice_difference=deltaG[i].tolist(),gradient_choice_difference_norm=float(mag(deltaG)[i]),meaningful_p_eff=bool(p_mask[i]),amplification_abs_qD_over_g_p=(float(abs(qD[i]/g[i]/p[i])) if p_mask[i] else None))
 top=lambda vals:[gp_values[int(i)] for i in np.argsort(vals)[-20:][::-1]]
 summary=dict(recorded_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),scope='Read-only frozen-state recovery diagnostic; only two independent linear global L2 projection solves, no nonlinear or time-path solve.',inputs_sha256=inputs,current_final_source_sha256=source,compiled_core_wrapper_contract_match_failed_source=exact_source,failed_source_sha256=prov['source_sha256'],diagnostic_script_sha256=sha(Path(__file__)),environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),time=data['time'].tolist(),dt=dt,n_nodes=len(nodes),n_elements=len(elems),n_gp=len(points),total_area=float(area),parameters=pr,definitions=dict(qD='consistent Q1 L2(g(v_h)*p_eff(u_h,c_h)), same2x2 quadrature as the UEL',qE='consistent Q1 L2(p_eff(u_h,c_h)), same frozen state and quadrature',alternative_gradient='g*grad(qE)+qE*2*(1-kappa)*v_h*grad(v_h)',species_flux='g*D*grad(c)-chi*c*grad(potential)',meaningful_pressure_threshold=.01*p_rms,p_eff_rms=p_rms),verification=verify,projection_comparison=projection,species_frozen_residual=species,bins=bins,cumulative_low_g_groups=cumul,overall=group(np.ones(len(g),bool),'all'),largest_q_over_g_error=top(abs(peq)),largest_gradient_difference=top(mag(deltaG)),limitations=['Unconverged failed Newton state only; neither recovered gradient is an exact physical/reference gradient.','No path solve, mesh/time convergence study, replacement formulation, or evidence that pressure recovery caused the nonlinear failure.','q/g is a diagnostic quotient rather than a required constitutive unknown; normalized ratios exclude |p_eff|<.01 global p_rms.','Two projections use the existing2x2 quadrature and fixed Q1 mesh; finite-h recovery, interpolation, and quadrature effects are not separated.'])
 (OUT/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n');(OUT/'gp_values.json').write_text(json.dumps(gp_values,indent=1,allow_nan=False)+'\n')
 np.savez(OUT/'projection_fields.npz',nodes=nodes,elems=elems,U=U,previous_U=Up,old_svars=data['svars'],q_trial=qtrial,q_degraded=qd,q_eff=qe,M_data=M.data,M_indices=M.indices,M_indptr=M.indptr,b_degraded=bg,b_effective=bp,g=get('g'),p_eff=get('p_eff'),weights=w,grad_q_trial=Gq,grad_q_degraded=GD,grad_g_q_eff=GA)
 if Path(__file__).resolve()!=(OUT/'diagnose.py').resolve():shutil.copy2(Path(__file__),OUT/'diagnose.py')
 names=['diagnose.py','summary.json','gp_values.json','projection_fields.npz']
 (OUT/'artifact_manifest.json').write_text(json.dumps(dict(recorded_utc=summary['recorded_utc'],artifacts_sha256={str((OUT/x).relative_to(ROOT)):sha(OUT/x) for x in names}),indent=2)+'\n')
 audit_source={str(x.relative_to(ROOT)):sha(x) for x in [ROOT/'element/generated_model.for',ROOT/'element/uel.f90',ROOT/'element/uel_contract.json']}
 assert all(audit_source[k]==source[k] for k in audit_source)
 low=cumul[2]
 results=dict(recorded_utc=summary['recorded_utc'],scope=summary['scope'],source_sha256=audit_source,input_sha256=inputs,artifact_sha256={str((OUT/x).relative_to(ROOT)):sha(OUT/x) for x in names+['artifact_manifest.json']},compiled_source_matches_failed_trial_provenance=exact_source,assessment=dict(local_low_g_equivalent_pressure_amplification_supported=True,gradient_discrepancy_concentrated_at_low_g=False,causal_role_in_failed_global_solve_established=False,alternative_recovery_validated_as_solver_fix=False),key_numbers=dict(trial_q_mass_equation_residual_l2=verify['pressure_projection_q_trial_residual_l2'],trial_vs_exact_qD_gradient_relative_l2=projection['trial_vs_qD_gradient_relative_l2'],low_g_threshold=10*kap,low_g_area_fraction=low['area_fraction'],low_g_pressure_quotient_error_squared_fraction=low['squared_norm_fractions']['pressure_q_over_g_error_squared'],low_g_gradient_difference_squared_fraction=low['squared_norm_fractions']['potential_gradient_difference_squared'],low_g_species_drift_difference_squared_fraction=low['squared_norm_fractions']['species_drift_difference_squared'],maximum_meaningful_qD_over_g_p_amplification=summary['overall']['amplification_abs_qD_over_g_p_where_p_meaningful']['max'],current_drift_l2=projection['species_drift_qD_l2'],alternative_drift_l2=projection['species_drift_product_l2'],alternative_over_current_drift_l2=projection['species_drift_product_l2']/projection['species_drift_qD_l2']),conclusion='At this frozen unconverged state, projection creates large q/g amplification near the degradation floor. However, the largest absolute gradient/flux differences between the two recovery choices occur outside the lowest-g region, and effective-pressure-first recovery gives a much larger global drift norm. Neither recovery is an exact gradient reference; this does not establish the cause of nonlinear failure or validate a replacement.',limitations=summary['limitations'])
 (OUT/'results.json').write_text(json.dumps(results,indent=2,allow_nan=False)+'\n')
 print(json.dumps(dict(verification=verify,projection_comparison=projection,species_frozen_residual=species,cumulative_low_g_groups=[dict(label=x['label'],area_fraction=x['area_fraction'],q_over_g_error_rms=x['q_over_g_minus_p'].get('rms'),gradient_diff_rms=x['gradient_difference_norm'].get('rms'),fractions=x['squared_norm_fractions'],amplification=x['amplification_abs_qD_over_g_p_where_p_meaningful']) for x in cumul]),indent=2,allow_nan=False))
if __name__=='__main__':main()
