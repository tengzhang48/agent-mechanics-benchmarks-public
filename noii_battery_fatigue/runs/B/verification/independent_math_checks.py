"""Independent analytical checks of the proposed discretization, not its code.

No abaqus_ufl or challenge element is imported.  The pressure-projection study
is the tensor-product-reducible case u_x=x**2, u_y=0 on a Cartesian grid,
with all material prefactors divided out: p=u_x,x=2*x.  It shows why checking
only a homogeneous patch does not establish accuracy of grad(project(p)).
"""

from __future__ import annotations

import json
import math

import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve


def projection_study():
    result = []
    points, weights = np.polynomial.legendre.leggauss(4)
    for n in (8, 16, 32, 64, 128):
        h = 1.0 / n
        x = np.linspace(0.0, 1.0, n + 1)
        u = x**2
        p_broken = np.diff(u) / h
        diagonal = np.full(n + 1, 2.0 * h / 3.0)
        diagonal[[0, -1]] = h / 3.0
        mass = diags((np.full(n, h / 6), diagonal, np.full(n, h / 6)),
                     (-1, 0, 1), format="csc")
        load = np.zeros(n + 1)
        load[:-1] += p_broken * h / 2
        load[1:] += p_broken * h / 2
        q = spsolve(mass, load)
        q_error_sq = grad_error_sq = 0.0
        for element in range(n):
            xx = x[element] + h * (points + 1.0) / 2.0
            qq = ((1.0 - points) * q[element]
                  + (1.0 + points) * q[element + 1]) / 2.0
            dq = (q[element + 1] - q[element]) / h
            q_error_sq += h / 2 * np.sum(weights * (qq - 2 * xx)**2)
            grad_error_sq += h * (dq - 2.0)**2
        result.append(dict(n=n, q_L2=math.sqrt(q_error_sq),
                           grad_q_L2=math.sqrt(grad_error_sq),
                           q_boundary=q[0]))
    return result


def scalar_oracles():
    E, nu, beta, kappa = 10.0, 0.28, 0.035, 1e-5
    mu = E / (2 * (1 + nu))
    lam = E * nu / ((1 + nu) * (1 - 2 * nu))
    bulk = lam + 2 * mu / 3
    exx, eyy, exy, d = 0.03, -0.01, 0.012, 0.4
    tr = exx + eyy
    psi = bulk / 2 * max(tr, 0)**2 + mu * ((exx - eyy)**2 / 2 + 2 * exy**2)
    g = (1 - kappa) * (1 - d)**2 + kappa
    p = (2 * mu + 2 * lam) * tr / 3
    alpha0, W0, alpha_T, G_c, ell = 0.001, 0.001, 0.0005, 0.002, 0.12
    W = g * psi
    alpha = alpha0 + max(W - W0, 0)
    f = 1.0
    if alpha > alpha_T:
        f = max(1e-15, max(0.0, 1 - 0.5 * math.log10(alpha / alpha_T))**2)
    H = max(0.3, 2 * ell * psi / (G_c * f))
    roundoff = []
    for history in (0.0, 1.0, 1e4, 1e8, 1e12, 1e15, 1e20):
        B = (1 - kappa) * history
        damage = B / (1 + B)
        reaction = damage - B * (1 - damage)
        roundoff.append(dict(H=history, d=damage, direct_reaction=reaction))
    return dict(mu=mu, lambda_=lam, K=bulk, psi_plus=psi,
                g=g, p_eff=p, W=W, alpha=alpha, f=f, H=H,
                compression_drift_B=4 * beta * (mu + lam) / 3,
                uniform_damage_roundoff=roundoff)


if __name__ == "__main__":
    print(json.dumps(dict(projection=projection_study(), scalar=scalar_oracles()),
                     indent=2))
