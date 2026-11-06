"""Exact steady zero-flux transport driven by prescribed u_x=A*x^2.

All displacement and damage DOFs are prescribed; c is prescribed only on
left/right. The pressure recovery therefore sees I_h u and the transport
equation sees its recovered gradient. No material or generator is imported.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.special import lambertw

from common import ROOT, fehost
import meshes
from psilag_phase_oracle import compiled


SIGNS = np.array([[-1., -1.], [1., -1.], [1., 1.], [-1., 1.]])


def shape(xi, eta):
    sx, sy = SIGNS.T
    return ((1+sx*xi)*(1+sy*eta)/4,
            np.column_stack((sx*(1+sy*eta), sy*(1+sx*xi)))/4)


def exact(x, p, amplitude):
    mu = p['E']/(2*(1+p['nu']))
    lam = p['E']*p['nu']/((1+p['nu'])*(1-2*p['nu']))
    P = 4*amplitude*(mu+lam)/3
    B = 4*p['beta']*(mu+lam)/3
    aa, bb = p['chi']*B/p['D'], p['chi']*P/p['D']
    cL = p['c0']
    # The principal real Lambert W branch inverts the chemical potential.
    argument = aa*cL*np.exp(aa*cL+bb*np.asarray(x))
    c = np.real(lambertw(argument))/aa
    dc = p['chi']*P*c/(p['D']+p['chi']*B*c)
    q = P*np.asarray(x)-B*(c-cL)
    dq = P-B*dc
    assert np.max(abs(p['D']*np.log(c/cL)+p['chi']*B*(c-cL)-p['chi']*P*np.asarray(x))) < 1e-13
    return c, dc, q, dq


def integrals(nodes, elems, values, index, p, amplitude):
    xy = nodes[elems]
    cv, qv = values[elems, index['c']], values[elems, index['q']]
    points, weights = np.polynomial.legendre.leggauss(5)
    total = dict(c_L2=0., grad_c_L2=0., q_L2=0., grad_q_L2=0., flux_L2=0.)
    Cint = flux_average = 0.
    for xi, wx in zip(points, weights):
        for eta, wy in zip(points, weights):
            N, D = shape(xi, eta)
            J = np.einsum('eai,aj->eij', xy, D)
            det = np.linalg.det(J)
            assert np.all(det > 0.)
            B = np.einsum('ai,eij->eaj', D, np.linalg.inv(J))
            x = np.einsum('a,ea->e', N, xy[:, :, 0])
            ce, gce, qe, gqe = exact(x, p, amplitude)
            ch, qh = cv@N, qv@N
            gc = np.einsum('ea,eai->ei', cv, B)
            gq = np.einsum('ea,eai->ei', qv, B)
            flux = p['D']*gc-p['chi']*ch[:, None]*gq
            w = wx*wy*det
            total['c_L2'] += np.sum(w*(ch-ce)**2)
            total['q_L2'] += np.sum(w*(qh-qe)**2)
            total['grad_c_L2'] += np.sum(w*((gc[:, 0]-gce)**2+gc[:, 1]**2))
            total['grad_q_L2'] += np.sum(w*((gq[:, 0]-gqe)**2+gq[:, 1]**2))
            total['flux_L2'] += np.sum(w*np.sum(flux*flux, axis=1))
            Cint += np.sum(w*ch)
            flux_average += np.sum(w*flux[:, 0])
    result = {name: float(math.sqrt(value)) for name, value in total.items()}
    mu = p['E']/(2*(1+p['nu']))
    lam = p['E']*p['nu']/((1+p['nu'])*(1-2*p['nu']))
    P, B = 4*amplitude*(mu+lam)/3, 4*p['beta']*(mu+lam)/3
    cL, cR = exact(np.array([0., 1.]), p, amplitude)[0]
    Cexact = (p['D']*(cR-cL)+p['chi']*B*(cR*cR-cL*cL)/2)/(p['chi']*P)
    result.update(C_int=float(Cint), C_int_exact=float(Cexact),
                  C_int_error=float(abs(Cint-Cexact)), flux_x_average=float(flux_average))
    return result


def run(module, contract, n, mesh_name, seed, amplitude, dt):
    started = time.monotonic()
    nodes, elems = (meshes.structured(n) if mesh_name == 'structured'
                    else meshes.distorted(n, seed=seed))
    p = dict(zip(contract['prop_names'], contract['props']))
    index = {name: i for i, name in enumerate(contract['dof_names'])}
    dpn = len(index)
    element = fehost.ContractUEL(module, contract, len(elems), dt)
    c, _, q, _ = exact(nodes[:, 0], p, amplitude)
    cL, cR = exact(np.array([0., 1.]), p, amplitude)[0]
    initial = np.zeros((len(nodes), dpn))
    initial[:, index['ux']] = amplitude*nodes[:, 0]**2
    initial[:, index['c']] = cL+(cR-cL)*nodes[:, 0]
    initial[:, index['q']] = q
    initial[:, index['vfrac']] = 1.
    bc = {}
    for i, (x, y) in enumerate(nodes):
        bc[i*dpn+index['ux']] = amplitude*x*x
        bc[i*dpn+index['uy']] = 0.
        bc[i*dpn+index['d']] = 0.
        if np.isclose(x, 0.):
            bc[i*dpn+index['c']] = cL
        elif np.isclose(x, 1.):
            bc[i*dpn+index['c']] = cR
    solved, report = fehost.solve(nodes, elems, element, lambda fraction: bc, 1, dt,
                                  U0=initial.ravel(), tol=1e-11,
                                  residual_scale=fehost.residual_weights(element, contract))
    values = solved.reshape(len(nodes), dpn)
    assert np.isfinite(values).all()
    constraint = float(np.max(abs((values[:, index['d']]-1)+values[:, index['vfrac']])))
    assert constraint < 1e-11
    nodal = float(np.max(abs(values[:, index['c']]-c)))
    return dict(n=n, mesh=mesh_name, seed=seed, nodal_c_max_error=nodal,
                intact_constraint_max=constraint, wall_s=time.monotonic()-started,
                **integrals(nodes, elems, values, index, p, amplitude), solver=report.summary())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--element-root', type=Path, default=ROOT)
    parser.add_argument('--out', type=Path, default=ROOT/'verification/quadratic_drift_results.json')
    parser.add_argument('--n', type=int, nargs='+', default=[8, 16, 32, 64])
    parser.add_argument('--seed', type=int, default=3)
    args = parser.parse_args()
    root = args.element_root.resolve()
    assert root == ROOT or ROOT in root.parents
    module, contract = compiled(root)
    records = []
    result = dict(scope='independent analytical zero-flux transport with prescribed quadratic displacement',
                  amplitude=.01, dt=1e12, seed=args.seed,
                  source_sha256={name: hashlib.sha256((root/name).read_bytes()).hexdigest()
                                 for name in contract['source']+['element/uel_contract.json']}, results=records)
    for mesh_name in ('structured', 'distorted'):
        last = None
        for n in args.n:
            row = run(module, contract, n, mesh_name, args.seed, .01, 1e12)
            if last:
                for name in ('c_L2', 'q_L2', 'grad_c_L2', 'grad_q_L2', 'flux_L2', 'C_int_error'):
                    row[name+'_order'] = math.log(last[name]/row[name])/math.log(n/last['n'])
            records.append(row)
            print(json.dumps(row), flush=True)
            args.out.write_text(json.dumps(result, indent=2)+'\n')
            last = row


if __name__ == '__main__':
    main()
