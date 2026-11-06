# Localized notch failure and equality-policy experiments

This diagnoses the original fully implicit six-field source with generated
SHA256 `744b5303512e26f261cd5d3158f85e2d8492fb20712bc2a901cb66062812e1ca`.
It does not change the host, generator, or production element. Experimental
source copies and builds are confined to workspace `tmp/`.

## Captured failure

The n=16 distorted seed-3 notch/side-charging variant stalled at time
2.3149414062500013 with DTIME=9.765624999991118e-5. The saved accepted histories
had maximum H around 2350 and alpha below its threshold, so fatigue reduction
was inactive. Damage alone accounts for the residual 1.02442777015e-8; other
field norms are at most 2e-14. The intact-fraction representation did not have
its high-H subtraction problem here.

The global Newton solve itself is accurate: relative linear residual
3.10e-14 and max damage update 9.14e-8. The smallest computed weighted
Jacobian eigenvalue is 2.09e-4. Nonetheless the actual residual increases
along that Newton direction, from 1.0244e-8 to 1.2871e-8 for a full step;
every tested positive line-search length fails to reduce it. The forward
linearization discrepancy is 1.28706726737e-8.

An independent Gauss-point calculation localizes essentially the entire
discrepancy to element 127, point 1 (both zero-based). Its old H is
201.3098415122905 and the exact production candidate is
201.3098415122875, lower by 2.98427949019e-12, or 105 ulps. The Newton direction
would raise that candidate by .00141159556745, crossing immediately from
the frozen side of the max-history kink to the active side. The missing
active-side contribution predicts the forward discrepancy to 2.20e-14.

This is not evidence of an incorrect derivative within a smooth branch.
The finite difference crosses a nearby nonsmooth history surface. Reducing
the generated complex-step size from 1e-10 to 1e-30 preserves the residual
bitwise and changes element tangent entries by at most 4.44e-15. Thus the
complex-step magnitude is not the cause of this captured failure.

Selecting the active derivative for this nearly tied point also fails as a
local remedy. Its Newton direction lowers the H candidate by
.00550599247600, crossing back to the frozen side, and the full-step residual
becomes 5.02e-8. The two one-sided tangent matrices differ by a rank-one
update, with singular values [3.6969207682,1.41e-15,6.72e-16,2.05e-16] on
its nonzero block. Its Sherman-Morrison coefficient is -1.25637440908;
an interpolated generalized Jacobian becomes singular at active fraction
.795941076777. Neither one-sided linearized root is consistent with its
own branch. This supports a nonsmooth loss of local regularity of this
discrete path. It does not prove that no more distant equilibrium exists.

Reproduction and numerical evidence: `debug_notch_failure.py`,
`notch_failure_review.json`, and the original capture under
`custom_runs/notch_distorted_16_debug/`.

## Exact-equality controls and the earlier Newton path

Changing H's strict `candidate > H_old` branch to `candidate >= H_old`
does not change any real max-history value. It selects the active candidate
derivative at exact equality. At the already-stalled capture there is no
exact equality, so this change leaves its tangent and failed Newton step
unchanged. It is not a universal repair for arbitrary near ties.

Nevertheless, a full replay with this change alone completed the n=16
distorted seed-3 notch problem: all 240 requested outputs over three cycles,
one cutback, minimum DTIME .05. Choosing the active derivative at earlier
exact ties steered the Newton trajectory away from the captured stalled
path. This is direct solver evidence for that particular loading and mesh,
not a proof that all softening paths will converge. Results are in
`notch_tie_replays/geH/`.

Changing W's strict comparison to `>=` alone, or together with H's `>=`,
instead failed at step 3 through the wrapper's monotone-history guard.
A diagnostic-only wrapper copy temporarily exposed the rejected trial's
state: 42 integration points decreased alpha by one ulp while W_new=W_old.
For example alpha_old=3.021006044529564e-9 became
3.0210060445295635e-9, a change of -4.1359030627651384e-25. The expression
`alpha_old + W_new - W_old` is left-associated, so the mathematically zero
increment can round downward. H and W stayed unchanged. The final wrapper's
strict monotonicity guard remains enabled.

Evaluating `alpha_old + (W_new-W_old)` first forms the nonnegative increment;
adding a nonnegative number cannot reduce a finite nonnegative alpha under
ordinary rounding. This preserves the exact model and is the preferable
arithmetic order. W's strict comparison can remain unchanged.

## Fully implicit candidate under test

`active_tie_study.py` prepares `tmp/active_tie_candidate` from the unmodified
generator, with only the H equality choice and parenthesized alpha increment.
Its generated SHA256 is
`3a2932ebc621ffaca412212be1985fc3cb433a88d5a2c1e42502536d7e70f3ae`.
The phase equation continues to use current trial H. No history lag, damage
clamp, row multiplier, or additional physics is introduced.

All five material verification states, independent assembled residual and
finite-difference checks, invariants, invalid-state checks, discarded-trial
repeatability, and large-H limits passed. The four independent frozen
hand-history cases match to at most 1.28e-14 scaled state error. Host results
for this candidate are recorded under `active_tie_study/`.

| case | requested cycles | outputs | cutbacks | minimum DTIME | maximum d |
|---|---:|---:|---:|---:|---:|
| D0, structured n=8 | 1 | 80 | 0 | .1 | .562214 |
| uniform fatigue floor, n=8 | 2 | 160 | 0 | .1 | 1+2.22e-16 |
| notch, structured n=16 | 2 | 160 | 2 | .05 | 1.048406 |
| notch, distorted n=16, seed 1 | 2 | 160 | 1 | .05 | 1.057148 |
| notch, distorted n=16, seed 3 | 2 | 160 | 1 | .05 | 1.055400 |

Every run completed with nondecreasing stored H and alpha. The floor run
reached H=9.350649350649437e17. The maximum intact constraint error across
these runs was 4.88e-16. Thus the notch overshoots are actual Q1 discrete
damage overshoots, not disagreement between d and vfrac. Nodal damage also
decreases at some points (largest per-increment decrease .01214 in the
seed-1 notch run). The equality choice fixes neither the discrete maximum
principle nor nodal irreversibility; mesh-refinement evidence remains needed
for their effect on continuum observables.

## History-lag alternative: consistency envelope

A frozen accepted H_n in the damage equation is a first-order partitioned
time discretization on smooth regular paths, while end-state H, alpha, and W
can still use the exact accepted-state formulas. Its phase residual evaluated
with the new history differs by

\[
 - (1-\kappa)\int N_a\,v_{n+1}(H_{n+1}-H_n).
\]

This defect is O(DTIME) only when the relevant fields and histories are
sufficiently smooth. At quasistatic snap-through a lag can select a different
jump branch and irreversibly accumulate a pre-damage W peak into alpha.
Uniform first-order accuracy through such jumps is not established. A bounded
linear predictor can reduce a smooth-path H defect to O(DTIME^2), but retains
first-order switch behavior and needs its own tests. These alternatives were
assessed, not substituted into the fully implicit candidate reported above.
