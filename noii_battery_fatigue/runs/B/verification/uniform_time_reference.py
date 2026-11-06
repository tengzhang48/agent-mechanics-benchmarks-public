"""Independent quasistatic uniform-fatigue time reference.

All displacement/concentration DOFs are prescribed as isotropic strain
e(t)=.03*sin(pi*t/8)**2 and c=c0. The exact spatial solution is uniform,
with psi+=2*K*e**2, v=1/(1+(1-kappa)*H), and d=1-v.  At each time the
implicit alpha equation is a scalar monotone root because W(alpha) is
nonincreasing in alpha. Bracketing avoids the monolithic UEL algorithm.
No element or code-generator module is imported.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path
import math

import numpy as np
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[1]
PARAMETERS = dict(E=10., nu=.28, kappa=1e-5, ell=.12, G_c=3.5e-6,
                  alpha_T=1e-10, k_fg=.5, strain_amplitude=.03, T=8.)


def fatigue(alpha, p):
    if alpha <= p['alpha_T']:
        return 1.
    return max(1e-15, max(0., 1-p['k_fg']*math.log10(alpha/p['alpha_T']))**2)


def update(psi, state, p, implicit=True, psi_old=0., predictor_floor=0.,
           intact_old=1., relaxation=0., relaxation_cap=1000., relaxation_form='history'):
    H0, alpha0, W0 = state
    a = 1-p['kappa']
    def gamma_at(phaseH):
        if relaxation_form == 'accepted_stiffness':
            return relaxation*(1+min(a*H0, relaxation_cap))
        if relaxation_form == 'phase_stiffness':
            return relaxation*(1+min(a*phaseH, relaxation_cap))
        return relaxation*min(H0, relaxation_cap)

    def history_and_work(alpha):
        H = max(H0, 2*p['ell']*psi/(p['G_c']*fatigue(alpha, p)))
        if implicit == 'psilag':
            driver = min(psi, max(psi_old, predictor_floor*psi))
            phaseH = max(H0, 2*p['ell']*driver/(p['G_c']*fatigue(alpha, p)))
        else:
            phaseH = H if implicit else H0
        gamma = gamma_at(phaseH)
        v = (1+gamma*intact_old)/(1+gamma+a*phaseH)
        return H, (p['kappa']+a*v*v)*psi

    _, Wa = history_and_work(alpha0)
    if Wa <= W0:
        alpha = alpha0
    elif not implicit or fatigue(alpha0, p) == 1e-15:
        alpha = alpha0+(Wa-W0)
    else:
        # The damage at frozen H0 bounds W(alpha) from above.
        gamma = gamma_at(H0)
        vmax = (1+gamma*intact_old)/(1+gamma+a*H0)
        upper_increment = max(0., (p['kappa']+a*vmax*vmax)*psi-W0)
        scale = p['alpha_T']

        def equation(x):
            aa = alpha0+x*scale
            _, W = history_and_work(aa)
            return x-max(W-W0, 0.)/scale

        upper = upper_increment/scale
        if equation(upper) < 0.:
            upper = np.nextafter(upper, np.inf)
        alpha = alpha0+scale*brentq(equation, 0., upper, xtol=1e-12, rtol=1e-14)
    H, W = history_and_work(alpha)
    assert H >= H0 and alpha >= alpha0
    return H, alpha, W


def solve(spc, p=PARAMETERS, cycles=2, implicit=True, times=None, predictor_floor=0.,
          relaxation=0., relaxation_cap=1000., relaxation_form='history'):
    mu = p['E']/(2*(1+p['nu']))
    lam = p['E']*p['nu']/((1+p['nu'])*(1-2*p['nu']))
    K = lam+2*mu/3
    if times is None:
        times = np.arange(1, cycles*spc+1)*p['T']/spc
    state = (0., 0., 0.)
    rows = []
    onset = None
    psi_old = 0.
    intact_old = 1.
    time_old = 0.
    for t in times:
        dt = float(t-time_old)
        if dt <= 0.:
            raise ValueError('Time samples must strictly increase from zero')
        step_relaxation = relaxation/math.sqrt(dt) if relaxation_form == 'sqrt_dt_phase_stiffness' else relaxation
        step_form = 'phase_stiffness' if relaxation_form == 'sqrt_dt_phase_stiffness' else relaxation_form
        e = p['strain_amplitude']*math.sin(math.pi*t/p['T'])**2
        psi = 2*K*e*e
        oldH = state[0]
        state = update(psi, state, p, implicit, psi_old, predictor_floor,
                       intact_old, step_relaxation, relaxation_cap, step_form)
        H, alpha, W = state
        if implicit == 'psilag':
            driver = min(psi, max(psi_old,predictor_floor*psi))
            phaseH = max(oldH, 2*p['ell']*driver/(p['G_c']*fatigue(alpha,p)))
        else:
            phaseH = H if implicit else oldH
        if step_form == 'accepted_stiffness':
            gamma = step_relaxation*(1+min((1-p['kappa'])*oldH, relaxation_cap))
        elif step_form == 'phase_stiffness':
            gamma = step_relaxation*(1+min((1-p['kappa'])*phaseH, relaxation_cap))
        else:
            gamma = step_relaxation*min(oldH, relaxation_cap)
        v = (1+gamma*intact_old)/(1+gamma+(1-p['kappa'])*phaseH)
        if onset is None and fatigue(alpha, p) == 1e-15:
            onset = float(t)
        rows.append((t, H, alpha, W, 1-v, fatigue(alpha, p)))
        psi_old = psi
        intact_old = v
        time_old = float(t)
    method = ('implicit fatigue, lagged psi' if implicit == 'psilag' else
              'implicit scalar bracketing' if implicit else 'accepted-H lag')
    return dict(spc=spc, method=method, predictor_floor=predictor_floor,
                relaxation=relaxation, relaxation_cap=relaxation_cap,
                relaxation_form=relaxation_form,
                cycles=cycles, final_H=state[0], final_alpha=state[1], final_W=state[2],
                floor_onset=onset, rows=np.asarray(rows))


def main():
    out = ROOT/'verification/uniform_time_reference'
    out.mkdir(exist_ok=True)
    records = []
    step_counts = (20, 40, 80, 160, 320, 640, 1280, 2560, 10240, 40960, 163840, 655360)
    for method in (True, False):
        for spc in step_counts:
            result = solve(spc, implicit=method)
            rows = result.pop('rows')
            records.append(result)
            if spc in (80, 40960):
                np.savetxt(out/f'{"implicit" if method else "lag"}_spc{spc}.csv', rows,
                           delimiter=',', header='t,H,alpha,W,d,f', comments='')
    reference = next(r for r in records if r['spc'] == max(step_counts)
                     and r['method'] == 'implicit scalar bracketing')
    for row in records:
        row['final_alpha_relative_difference_to_finest_implicit'] = row['final_alpha']/reference['final_alpha']-1
    data = dict(parameters=PARAMETERS,
                scope='Independent uniform quasistatic solution; finest time is a reference, not exact continuum',
                records=records)
    (out/'results.json').write_text(json.dumps(data, indent=2)+'\n')
    with (out/'summary.csv').open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]))
        writer.writeheader(); writer.writerows(records)
    print(json.dumps(data, indent=2))


if __name__ == '__main__':
    main()
