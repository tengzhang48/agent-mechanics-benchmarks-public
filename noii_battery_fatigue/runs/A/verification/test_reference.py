"""Check the independent FV Jacobian and demonstrate broken-gradient failure."""
import json
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve
from common import ROOT,write_json
from reference_1d import fv_evaluate

n=12; h=1/n; rng=np.random.default_rng(92)
c=.45+.1*np.sin(np.pi*(np.arange(n)+.5)/n)
d=.12+.02*rng.random(n); v=np.r_[c,d]
old=v.copy(); old[:n]-=.015
diag=np.full(n,-2.); diag[[0,-1]]=-1
L=sp.diags([np.ones(n-1),diag,np.ones(n-1)],[-1,0,1],format='csr')/h**2
sv=np.tile([.05,.0008,1e-5],(n,1))
r,k,sn=fv_evaluate(v,old,sv,.1,.65,L)
fd=np.zeros((2*n,2*n)); eps=1e-7
for j in range(2*n):
    inc=np.zeros_like(v); inc[j]=eps
    fd[:,j]=(fv_evaluate(v+inc,old,sv,.1,.65,L)[0]-
             fv_evaluate(v-inc,old,sv,.1,.65,L)[0])/(2*eps)
relative=float(np.linalg.norm(fd-k.toarray())/np.linalg.norm(k.toarray()))
assert relative<2e-8,relative
r,k,_=fv_evaluate(v,old,sv,.1,.65,L,lag=True)
for j in range(2*n):
    inc=np.zeros_like(v); inc[j]=eps
    fd[:,j]=(fv_evaluate(v+inc,old,sv,.1,.65,L,True)[0]-
             fv_evaluate(v-inc,old,sv,.1,.65,L,True)[0])/(2*eps)
lag_relative=float(np.linalg.norm(fd-k.toarray())/np.linalg.norm(k.toarray()))
assert lag_relative<2e-8,lag_relative

# u_x=x^2/2, c=c0, d=0 => p_eff=a*x. Q1 displacement on rectangles
# has cellwise constant dux/dx, so its broken pressure gradient is zero.
projection=[]
for n in [8,16,32,64,128]:
    h=1/n; x=np.linspace(0,1,n+1); mid=.5*(x[:-1]+x[1:])
    diag=np.full(n+1,2*h/3); diag[[0,-1]]=h/3
    M=sp.diags([np.full(n,h/6),diag,np.full(n,h/6)],[-1,0,1],format='csr')
    b=np.zeros(n+1); b[:-1]+=.5*h*mid; b[1:]+=.5*h*mid
    q=spsolve(M,b)
    graderror=float(np.sqrt(h*np.sum((np.diff(q)/h-1)**2)))
    qerror=0.
    for xi,w in zip(*np.polynomial.legendre.leggauss(3)):
        alpha=.5*(1+xi)
        interp=(1-alpha)*q[:-1]+alpha*q[1:]
        exact=(1-alpha)*x[:-1]+alpha*x[1:]
        qerror+=float(.5*h*w*np.sum((interp-exact)**2))
    projection.append({'n':n,'projected_pressure_L2':qerror**.5,
                       'projected_gradient_L2':graderror,'broken_gradient_L2':1.})
assert projection[-1]['projected_gradient_L2']<projection[0]['projected_gradient_L2']
result={'independent_FV_jacobian_relative_error':relative,
        'lagged_FV_jacobian_relative_error':lag_relative,
        'quadratic_displacement_pressure_counterexample':projection}
write_json(ROOT/'verification'/'reference_tests.json',result)
print(json.dumps(result,indent=2))
