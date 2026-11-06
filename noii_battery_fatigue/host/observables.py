"""Host-side observables, computed from nodal values with Q4 interpolation (independent of the element).

D_int = int d, C_int = int c and E_el = int g(d) psi(eps_e) use 2x2 Gauss points; A_damage is the area
where d >= 0.5, from a 4x4 sub-sampling of each element. d_max (the nodal maximum) is taken by the runner.
"""
import numpy as np

GP = np.array([-1.0, 1.0]) / np.sqrt(3.0)
SUB = (np.arange(4) + 0.5) / 4 * 2 - 1  # 4x4 sub-sampling points for the damage-area indicator


def shape(xi, eta):
    s = np.array([-1, 1, 1, -1.0])
    t = np.array([-1, -1, 1, 1.0])
    N = 0.25 * (1 + s * xi) * (1 + t * eta)
    dN = np.column_stack([0.25 * s * (1 + t * eta), 0.25 * t * (1 + s * xi)])
    return N, dN


_GEOM = {}


def _geometry(nodes, elems):
    """Per-element Q4 data at the 2x2 Gauss points and the 4x4 sub-points, computed once per mesh."""
    key = (id(nodes), id(elems), len(nodes), len(elems))
    if key in _GEOM:
        return _GEOM[key]
    X = nodes[elems]                                         # (M, 4, 2)
    gp = []
    for xi in GP:
        for eta in GP:
            N, dN = shape(xi, eta)
            J = np.einsum("ai,maj->mij", dN, X)              # (M, 2, 2), J[i, j] = dX_j/dxi_i
            detJ = np.linalg.det(J)
            B = np.einsum("ai,mji->maj", dN, np.linalg.inv(J))   # (M, 4, 2), dN/dX
            gp.append((N, B, detJ))
    sub = []
    for xi in SUB:
        for eta in SUB:
            N, dN = shape(xi, eta)
            detJ = np.linalg.det(np.einsum("ai,maj->mij", dN, X))
            sub.append((N, detJ * (2.0 / 4) ** 2))
    _GEOM[key] = (gp, sub)
    return _GEOM[key]


def observables(nodes, elems, U, dpn, comp, prm):
    """D_int, C_int, E_el and A_damage of the nodal state U (vectorized over elements)."""
    E, nu, beta, c0, kappa = prm["E"], prm["nu"], prm["beta"], prm["c0"], prm["kappa"]
    mu = E / (2 * (1 + nu))
    lam = E * nu / ((1 + nu) * (1 - 2 * nu))
    ux, uy = U[comp["ux"]::dpn][elems], U[comp["uy"]::dpn][elems]   # (M, 4)
    c, d = U[comp["c"]::dpn][elems], U[comp["d"]::dpn][elems]
    gp, sub = _geometry(nodes, elems)
    D_int = C_int = E_el = A_dam = 0.0
    for N, B, detJ in gp:
        gx = np.einsum("ma,mak->mk", ux, B)                  # (M, 2): dux/dX, dux/dY
        gy = np.einsum("ma,mak->mk", uy, B)
        cg, dg = c @ N, d @ N
        sw = beta * (cg - c0)
        exx, eyy, exy = gx[:, 0] - sw, gy[:, 1] - sw, 0.5 * (gx[:, 1] + gy[:, 0])
        psi = mu * (exx ** 2 + eyy ** 2 + 2 * exy ** 2) + 0.5 * lam * (exx + eyy) ** 2
        g = (1 - kappa) * (1 - dg) ** 2 + kappa
        D_int += float(np.sum(dg * detJ))
        C_int += float(np.sum(cg * detJ))
        E_el += float(np.sum(g * psi * detJ))
    for N, w in sub:
        A_dam += float(np.sum((d @ N >= 0.5) * w))
    return {"D_int": D_int, "C_int": C_int, "E_el": E_el, "A_damage": A_dam}
