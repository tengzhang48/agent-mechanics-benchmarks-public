"""Independent 2D matrix implementation of the public equations and Q4 residual.

This module does not import the model declaration or generator. Its purpose is
to catch transcription, trace, sign, mapping and state-order errors.
"""
import numpy as np

GP = np.array([[-1, -1], [1, -1], [1, 1], [-1, 1]]) / np.sqrt(3.0)


def shape(xi, eta):
    signs = np.array([[-1., -1.], [1., -1.], [1., 1.], [-1., 1.]])
    a, b = signs.T
    N = (1+a*xi)*(1+b*eta)/4
    dN = np.column_stack((a*(1+b*eta), b*(1+a*xi)))/4
    return N, dN


def constitutive(grad_u, c, d, old, prm):
    E, nu = prm["E"], prm["nu"]
    mu = E / (2*(1+nu))
    lam = E*nu/((1+nu)*(1-2*nu))
    ee = (grad_u+grad_u.T)/2 - prm["beta"]*(c-prm["c0"])*np.eye(2)
    trace = np.trace(ee)
    dev = ee - trace/2*np.eye(2)
    psip = (lam+2*mu/3)/2*max(trace, 0)**2 + mu*np.sum(dev*dev)
    g = (1-prm["kappa"])*(1-d)**2+prm["kappa"]
    seff = 2*mu*ee+lam*trace*np.eye(2)
    stress = g*seff
    W = g*psip
    H0, alpha0, W0 = old
    alpha = alpha0+max(W-W0, 0)
    if alpha <= prm["alpha_T"]:
        fatigue = 1.
    else:
        fatigue = max(1e-15, max(0., 1-prm["k_fg"]*np.log10(alpha/prm["alpha_T"]))**2)
    H = max(H0, 2*prm["ell"]*psip/(prm["G_c"]*fatigue))
    return stress, np.trace(stress)/3, np.array([H, alpha, W]), g


def enhanced_strain(coords,nodal,prm):
    """Independent local mode fit by weighted least squares, not a 4x4 inverse."""
    mu=prm["E"]/(2*(1+prm["nu"]))
    lam=prm["E"]*prm["nu"]/((1+prm["nu"])*(1-2*prm["nu"]))
    C=np.array([[lam+2*mu,lam,0.],[lam,lam+2*mu,0.],[0.,0.,mu]])
    root=np.linalg.cholesky(C).T
    raw=[];weights=[];eps=[];factors=[]
    for xi,eta in GP:
        N,dn=shape(xi,eta);J=coords.T@dn;inv=np.linalg.inv(J);B=dn@inv
        raw.append(np.array([[-2*xi,0.],[0.,-2*eta]])@inv)
        weights.append(np.linalg.det(J))
        gu=nodal[:,:2].T@B;sw=prm["beta"]*(N@nodal[:,2]-prm["c0"])
        eps.append(np.array([gu[0,0]-sw,gu[1,1]-sw,gu[0,1]+gu[1,0]]))
        factors.append((1-prm["kappa"])*(N@nodal[:,5])**2+prm["kappa"])
    raw=np.array(raw);weights=np.array(weights)
    raw-=np.einsum("q,qij->ij",weights,raw)[None]/weights.sum()
    modes=[];lhs=[];rhs=[]
    for k in range(4):
        Q=np.zeros((3,4))
        for j in range(2):
            gx,gy=raw[k,j]
            Q[:,j]=[gx,0.,gy];Q[:,j+2]=[0.,gy,gx]
        modes.append(Q)
        lhs.append(np.sqrt(weights[k]*factors[k])*root@Q)
        rhs.append(-np.sqrt(weights[k]*factors[k])*root@eps[k])
    q=np.linalg.lstsq(np.vstack(lhs),np.concatenate(rhs),rcond=None)[0]
    corrections=[]
    for Q in modes:
        a,b,c=Q@q
        corrections.append(np.array([[a,c/2],[c/2,b]]))
    return corrections


def residual(coords, U, Uold, svars, dt, prm, rule="gauss"):
    nodal, previous = np.reshape(U, (4, 6)), np.reshape(Uold, (4, 6))
    R = np.zeros((4, 6))
    state = np.zeros((4, 3))
    corrections=enhanced_strain(coords,nodal,prm) if rule=="enhanced_lagged" else [0.]*4
    for k, (xi, eta) in enumerate(GP):
        N, dn = shape(xi, eta)
        jac = coords.T@dn
        B = dn@np.linalg.inv(jac)
        det = np.linalg.det(jac)
        grad_u = nodal[:, :2].T@B+corrections[k]
        cg, dg, pg, sg = N@nodal[:, 2:]
        co = N@previous[:, 2]
        stress, pressure, state[k], g = constitutive(grad_u, cg, 1-sg, svars.reshape(4, 3)[k], prm)
        gc, gs, gp = (nodal[:, j]@B for j in [2, 5, 4])
        Jc = g*prm["D"]*gc-prm["chi"]*cg*gp
        R[:, :2] += B@stress.T*det
        R[:, 2] += (N*(cg-co)/dt+B@Jc)*det
        R[:, 3] += prm["ell"]**2*(B@gs)*det
        if rule in ("gauss","scaled_gauss"):
            R[:,3]+=N*((1+(1-prm["kappa"])*state[k,0])*sg-1)*det
        elif rule=="gauss_lumped":
            R[:,3]+=N*((1+(1-prm["kappa"])*state[k,0])*nodal[:,5]-1)*det
        elif rule in ("lagged_gauss_lumped","enhanced_lagged"):
            oldH=svars.reshape(4,3)[k,0]
            R[:,3]+=N*((1+(1-prm["kappa"])*oldH)*nodal[:,5]-1)*det
        elif rule=="lagged_gauss":
            oldH=svars.reshape(4,3)[k,0]
            R[:,3]+=N*((1+(1-prm["kappa"])*oldH)*sg-1)*det
        R[:, 4] += N*(pg-pressure)*det
        R[:, 5] += N*(dg+sg-1)*det
    if rule=="lumped":
        for k,(xi,eta) in enumerate([(-1,-1),(1,-1),(1,1),(-1,1)]):
            N,dn=shape(xi,eta);J=coords.T@dn;B=dn@np.linalg.inv(J)
            cg=N@nodal[:,2];sg=N@nodal[:,5]
            stress,pressure,state[k],g=constitutive(nodal[:,:2].T@B,cg,1-sg,svars.reshape(4,3)[k],prm)
            R[:,3]+=N*((1+(1-prm["kappa"])*state[k,0])*sg-1)*np.linalg.det(J)
    elif rule not in ("gauss","gauss_lumped","scaled_gauss","lagged_gauss_lumped","lagged_gauss","enhanced_lagged"):
        raise ValueError(rule)
    if rule=="scaled_gauss":
        R[:,3]*=np.sqrt((1-prm["kappa"])*nodal[:,5]**2+prm["kappa"])
    return R.ravel(), state.ravel()
