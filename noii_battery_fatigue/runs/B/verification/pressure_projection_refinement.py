"""Independent recovery check for q = L2 projection of d(I_h x^2)/dx.

The exact pressure potential is q=2x, with gradient (2,0). This isolates the
chosen projection of a displacement gradient; no element/model or generator
code is imported. It is a gradient-recovery check, not a coupled PDE solution.
"""
import json

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve

from common import ROOT
import meshes


SIGNS = np.array([[-1., -1.], [1., -1.], [1., 1.], [-1., 1.]])


def shape(xi, eta):
    sx, sy = SIGNS.T
    return ((1+sx*xi)*(1+sy*eta)/4,
            np.column_stack((sx*(1+sy*eta), sy*(1+sx*xi)))/4)


def run(n, seed):
    nodes, elems = meshes.structured(n) if seed is None else meshes.distorted(n, seed=seed)
    xy = nodes[elems]
    prescribed_u = xy[:, :, 0]**2
    cell_mass = np.zeros((len(elems), 4, 4))
    cell_load = np.zeros((len(elems), 4))
    for xi, eta in SIGNS/np.sqrt(3):
        N, D = shape(xi, eta)
        J = np.einsum('eai,aj->eij', xy, D)
        det = np.linalg.det(J)
        assert np.all(det > 0.)
        B = np.einsum('ai,eij->eaj', D, np.linalg.inv(J))
        broken_p = np.sum(prescribed_u*B[:, :, 0], axis=1)
        cell_mass += np.einsum('a,b,e->eab', N, N, det)
        cell_load += np.einsum('a,e->ea', N, det*broken_p)
    ri = np.broadcast_to(elems[:, :, None], cell_mass.shape)
    ci = np.broadcast_to(elems[:, None, :], cell_mass.shape)
    M = coo_matrix((cell_mass.ravel(), (ri.ravel(), ci.ravel())),
                   shape=(len(nodes), len(nodes))).tocsr()
    load = np.zeros(len(nodes))
    np.add.at(load, elems.ravel(), cell_load.ravel())
    q = spsolve(M, load)
    assert np.max(abs(M@q-load)) < 1e-12
    points, weights = np.polynomial.legendre.leggauss(5)
    qerr = gerr = 0.
    for xi, wx in zip(points, weights):
        for eta, wy in zip(points, weights):
            N, D = shape(xi, eta)
            J = np.einsum('eai,aj->eij', xy, D)
            det = np.linalg.det(J)
            B = np.einsum('ai,eij->eaj', D, np.linalg.inv(J))
            qp = q[elems]@N
            gp = np.einsum('ea,eai->ei', q[elems], B)
            xe = np.einsum('a,ea->e', N, xy[:, :, 0])
            qerr += np.sum(wx*wy*det*(qp-2*xe)**2)
            gerr += np.sum(wx*wy*det*np.sum((gp-np.array([2., 0.]))**2, axis=1))
    return dict(mesh='structured' if seed is None else 'distorted', n=n, seed=seed,
                q_L2=float(np.sqrt(qerr)), grad_q_L2=float(np.sqrt(gerr)))


def main():
    rows = []
    for seed in (None, 1, 3):
        last = None
        for n in (8, 16, 32, 64, 128):
            row = run(n, seed)
            if last:
                row['q_L2_order'] = float(np.log2(last['q_L2']/row['q_L2']))
                row['grad_q_L2_order'] = float(np.log2(last['grad_q_L2']/row['grad_q_L2']))
            rows.append(row)
            print(json.dumps(row), flush=True)
            last = row
    result = dict(scope='independent displacement-gradient L2 recovery, not a coupled PDE solve',
                  exact='u_x=x^2, p=du_x/dx=2x, grad p=(2,0)', results=rows)
    (ROOT/'verification/pressure_projection_refinement.json').write_text(json.dumps(result, indent=2)+'\n')


if __name__ == '__main__':
    main()
