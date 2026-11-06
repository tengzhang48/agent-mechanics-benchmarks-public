"""Independent scalar Newton experiments for high-H damage representations."""
import json
import math
from pathlib import Path

import numpy as np


A_FACTOR = 1-1e-5


def solve(H, old_H=0.0, method="unscaled", initial_d=None, max_iter=25, initial_aux=None):
    B = A_FACTOR*H
    A = 1+B
    Aold = 1+A_FACTOR*old_H
    dold = (Aold-1)/Aold if initial_d is None else initial_d
    oldscale = max(1e-12, min(1.0, abs(1-dold)))
    if method == "gap_aux":
        value = np.array([dold, 1/Aold])
    elif method == "logH_aux":
        value = np.array([dold, math.log(Aold)])
    else:
        value = np.array([dold])
    if initial_aux is not None:
        value[1] = initial_aux

    def residual_jacobian(x):
        d = x[0]
        raw = d-B*(1-d)
        if method == "unscaled":
            return np.array([raw]), np.array([[A]])
        if method == "old_gap":
            return np.array([oldscale*raw]), np.array([[oldscale*A]])
        if method == "current_gap":
            scale = math.sqrt((1-d)**2+1e-24)
            ds = -(1-d)/scale
            return np.array([scale*raw]), np.array([[scale*A+ds*raw]])
        if method == "gap_aux":
            z = x[1]
            return np.array([d-B*z, (d-1)+z]), np.array([[1., -B], [1., 1.]])
        if method == "logH_aux":
            s = x[1]
            scale = math.exp(-s)
            r = scale*raw
            return np.array([r, s-math.log(A)]), np.array([[scale*A, -r], [0., 1.]])
        raise ValueError(method)

    trace = []
    converged = False
    for iteration in range(max_iter+1):
        r, jac = residual_jacobian(value)
        norm = float(np.linalg.norm(r))
        trace.append(dict(iteration=iteration, value=value.tolist(), residual=norm))
        if norm < 1e-8:
            converged = True
            break
        if iteration == max_iter:
            break
        step = np.linalg.solve(jac, -r)
        accepted = False
        for k in range(30):
            trial = value+step*2.**(-k)
            trial_r, _ = residual_jacobian(trial)
            if np.linalg.norm(trial_r) < norm:
                value = trial
                accepted = True
                break
        if not accepted:
            break
    return dict(method=method,H=H,old_H=old_H,initial_d=dold,
                converged=converged,iterations=trace[-1]["iteration"],
                solution=value.tolist(),residual=trace[-1]["residual"],
                d_error=abs(value[0]-B/A),trace=trace)


if __name__ == "__main__":
    results=[]
    for H, old_H in [(1.,0.), (1e8,0.), (1e12,0.), (1e15,0.), (1e20,0.),
                      (1e15,1e6), (1e15,1e8), (1e20,1e15)]:
        for method in ("unscaled","old_gap","current_gap","logH_aux","gap_aux"):
            results.append(solve(H,old_H,method))
    results.append(solve(1.,0.,"current_gap",initial_d=.999))
    for H in (1., 1e12, 1e15, 1e20):
        row=solve(H,0.,"gap_aux",initial_aux=0.0)
        row["initial_aux"]=0.0
        results.append(row)
    path=Path(__file__).resolve().parent/"high_history_representation.json"
    path.write_text(json.dumps(results,indent=2)+"\n")
    for row in results:
        print(json.dumps({key:value for key,value in row.items() if key!="trace"}))
