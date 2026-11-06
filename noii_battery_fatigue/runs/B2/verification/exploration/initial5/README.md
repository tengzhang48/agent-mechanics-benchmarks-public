# Rejected five-field prototype

Archived solely to reproduce the known failure in `check_limits.py`.
This is not the submitted element; the scoring contract is in `element/`.

This prototype solves (ux,uy,c,d,p), evaluating the damage reaction as
`d-(1-kappa)*(1-d)*H`. At alpha_old=.02 and shear .01, fatigue reaches its
floor and H is about 1.34e16. Double precision cannot resolve 1-d accurately
enough to reduce the residual. The eight Newton residual norms in
`verification/limit_checks.json` reproduce the failure in compiled Fortran.

`model.py` is the original generator declaration. The kernel was generated
with the unmodified pinned abaqus_ufl standard Quad4 generator, with the sole
text change UEL -> CHEMO_KERNEL to allow the handwritten entry wrapper.
The included contract is historical and is overridden with archive-relative
source paths by `check_limits.py`. Do not use it to score the final element.
