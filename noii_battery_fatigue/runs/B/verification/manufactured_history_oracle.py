"""Independent elliptic PDE oracle for spatially varying, frozen history.

d*=.4+.1*cos(pi*x)*cos(pi*y) has zero normal derivative on the unit square.
H*=(d*-ell**2*laplacian(d*))/((1-kappa)*(1-d*)) makes it an exact solution
of the specified damage PDE without an added source.  Prescribed u=0,c=c0
give psi+=0, so this positive initial H remains frozen.  This checks the
spatial reaction discretization against an analytical PDE solution; it does
not import the material model, code generator, or generated assembly.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib
import json
import math
from pathlib import Path
import sys
import time

import numpy as np

from common import ROOT, fehost
import meshes


def shape(xi, eta):
    signs = np.array([[-1., -1.], [1., -1.], [1., 1.], [-1., 1.]])
    s, t = signs.T
    N = (1+s*xi)*(1+t*eta)/4
    D = np.column_stack((s*(1+t*eta), t*(1+s*xi)))/4
    return N, D


def exact(xy, ell, kappa):
    x, y = np.asarray(xy).T
    mode = np.cos(np.pi*x)*np.cos(np.pi*y)
    d = .4+.1*mode
    gradient = -.1*np.pi*np.column_stack((np.sin(np.pi*x)*np.cos(np.pi*y),
                                        np.cos(np.pi*x)*np.sin(np.pi*y)))
    laplacian = -.2*np.pi**2*mode
    H = (d-ell**2*laplacian)/((1-kappa)*(1-d))
    return d, gradient, H


def error_integrals(nodes, elems, d, ell, kappa):
    # Five-point integration is independent of the element's 2x2 quadrature.
    points, weights = np.polynomial.legendre.leggauss(5)
    l2 = h1 = integral = 0.
    xy, dv = nodes[elems], d[elems]
    for xi, wx in zip(points, weights):
        for eta, wy in zip(points, weights):
            N, D = shape(xi, eta)
            J = np.einsum('eai,aj->eij', xy, D)
            det = np.linalg.det(J)
            assert np.all(det > 0)
            B = np.einsum('ai,eij->eaj', D, np.linalg.inv(J))
            gp = np.einsum('a,eai->ei', N, xy)
            de, ge, _ = exact(gp, ell, kappa)
            dh = dv @ N
            gh = np.einsum('ea,eai->ei', dv, B)
            w = wx*wy*det
            l2 += np.sum(w*(dh-de)**2)
            h1 += np.sum(w*np.sum((gh-ge)**2, axis=1))
            integral += np.sum(w*dh)
    return dict(d_L2=math.sqrt(l2), grad_d_L2=math.sqrt(h1), D_int=float(integral),
                D_int_error=float(abs(integral-.4)))


def run(module, contract, mesh_name, n, seed, history_scale=1.):
    started = time.monotonic()
    nodes, elems = (meshes.structured(n) if mesh_name == 'structured'
                    else meshes.distorted(n, seed=seed))
    params = dict(zip(contract['prop_names'], contract['props']))
    # Use c=c0=0 for exactly zero strain energy, including in floating point.
    # The frozen-history damage PDE and its analytical solution are unchanged.
    props = list(contract['props'])
    props[contract['prop_names'].index('c0')] = 0.
    params['c0'] = 0.
    el = fehost.ContractUEL(module, {**contract, 'props': props}, len(elems), .1)
    ell, kappa = params['ell'], params['kappa']
    slots = np.asarray(contract['history_slots']['H'])-1
    assert len(slots) == 4, 'oracle initializes cyclic 2x2 history points'
    for gp, (xi, eta) in enumerate(np.array([[-1., -1.], [1., -1.],
                                            [1., 1.], [-1., 1.]])/np.sqrt(3)):
        N, _ = shape(xi, eta)
        xy = np.einsum('a,eai->ei', N, nodes[elems])
        el.svars[:, slots[gp]] = history_scale*exact(xy, ell, kappa)[2]
    sv0 = el.svars.copy()
    dpn, ci, di, vi = el.dof_per_node, el.comp('c'), el.comp('d'), el.comp('vfrac')
    initial = np.zeros((len(nodes), dpn))
    initial[:, ci] = params['c0']
    # The exact solution is stationary. Initialize old physical damage from
    # its nodal samples so any declared vanishing time relaxation has zero
    # continuum contribution; interpolation errors are still measured.
    initial[:, di] = exact(nodes, ell, kappa)[0]
    bc = {dpn*i+el.comp(name): (params['c0'] if name == 'c' else 0.)
          for i in range(len(nodes)) for name in ('ux', 'uy', 'c')}
    U, report = fehost.solve(nodes, elems, el, lambda fraction: bc, 1, .1,
                             U0=initial.ravel(), tol=1e-11,
                             residual_scale=fehost.residual_weights(el, contract))
    values = U.reshape(len(nodes), dpn)
    assert np.isfinite(U).all() and np.array_equal(sv0, el.svars)
    d, _, _ = exact(nodes, ell, kappa)
    nodal = float(np.max(abs(values[:, di]-d)))
    constraint = float(np.max(abs((values[:, di]-1)+values[:, vi])))
    assert constraint < 1e-11
    return dict(mesh=mesh_name, n=n, seed=seed, history_scale=history_scale,
                nodal_d_max_error=nodal, intact_constraint_max=constraint,
                history_bitwise_unchanged=True, wall_s=time.monotonic()-started,
                **error_integrals(nodes, elems, values[:, di], ell, kappa),
                solver=report.summary())


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--element-root', type=Path, default=ROOT)
    ap.add_argument('--out', type=Path, default=ROOT/'verification/manufactured_history_results')
    ap.add_argument('--n', type=int, nargs='+', default=[8, 16, 32, 64])
    ap.add_argument('--mesh', nargs='+', default=['structured', 'distorted'])
    ap.add_argument('--seed', type=int, default=7)
    args = ap.parse_args()
    root = args.element_root.resolve()
    assert root == ROOT or ROOT in root.parents
    contract = fehost.load_contract(str(root))
    sources = contract['source'] if isinstance(contract['source'], list) else [contract['source']]
    manifest = {name: hashlib.sha256((root/name).read_bytes()).hexdigest() for name in sources}
    digest = hashlib.sha256(json.dumps(manifest, sort_keys=True).encode()).hexdigest()[:12]
    name = 'history_oracle_'+digest
    build = ROOT/'tmp'/name
    if list(build.glob(name+'*.so')):
        sys.path.insert(0, str(build))
        module = importlib.import_module(name)
    else:
        module = fehost.build_module(str(root), contract, str(build), name)
    args.out.mkdir(parents=True, exist_ok=True)
    rows = []
    result = dict(scope='independent analytical elliptic damage PDE with frozen spatial H',
                  exact='d=.4+.1*cos(pi*x)*cos(pi*y)',
                  source_sha256={**manifest,'element/uel_contract.json':hashlib.sha256((root/'element/uel_contract.json').read_bytes()).hexdigest()},
                  results=rows)
    for mesh in args.mesh:
        for n in args.n:
            row = run(module, contract, mesh, n, args.seed)
            rows.append(row)
            print(json.dumps(row), flush=True)
            (args.out/'results.json').write_text(json.dumps(result, indent=2)+'\n')
    rates = {}
    for mesh in args.mesh:
        sequence = sorted((r for r in rows if r['mesh'] == mesh), key=lambda r: r['n'])
        rates[mesh] = []
        for before, after in zip(sequence, sequence[1:]):
            log_ratio = math.log(after['n']/before['n'])
            p2 = math.log(before['d_L2']/after['d_L2'])/log_ratio
            p1 = math.log(before['grad_d_L2']/after['grad_d_L2'])/log_ratio
            assert p2 > 1.6 and p1 > .8, (mesh, before, after)
            rates[mesh].append(dict(coarse=before['n'], fine=after['n'], d_L2_order=p2,
                                    grad_d_L2_order=p1))
    broken = run(module, contract, 'structured', max(args.n), args.seed, history_scale=1.1)
    assert broken['d_L2'] > .002 and broken['d_L2'] > 50*max(r['d_L2'] for r in rows if r['n']==max(args.n))
    result.update(rates=rates, broken_control=dict(history_scaled_1p1='detected', result=broken))
    (args.out/'results.json').write_text(json.dumps(result, indent=2)+'\n')
    columns = [key for key in rows[0] if key != 'solver']
    with (args.out/'results.csv').open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        writer.writerows({key: row[key] for key in columns} for row in rows)


if __name__ == '__main__':
    main()
