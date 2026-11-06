# Results

## Official results (hidden problem, 21 October 2026)

Each element was scored by `../scoring/score.py` on the host the agents had, on a hidden problem of the same kind as the
practice problem D0: other boundary conditions, loading and parameters, structured and distorted meshes, the first 10
cycles, and a stress variant that drives the damage to 1. The hidden problem and its reference stay private (see
`../README.md`); the numbers below are the scorer's output.

| Check | A | B | B2 |
|---|---|---|---|
| R0 runs complete | pass | pass | pass |
| C1 key observables: d_max, ∫d, ∫c (limits 4 %, 4 %, 0.5 %) | **fail**: 6.0 %, 8.5 %, 0.020 % | **fail**: 8.3 %, 9.3 %, 0.029 % | pass: 0.66 %, 0.55 %, 0.0041 % |
| C1 key fields: u_x, u_y, c, d (limit 1 %) | **fail**: 0.056 %, 0.033 %, 0.043 %, 3.3 % | **fail**: 0.071 %, 0.032 %, 0.048 %, 3.7 % | pass: 0.0091 %, 0.013 %, 0.0060 %, 0.20 % |
| C2 converges to the key | pass | pass | pass |
| C3 distortion (largest of d_max, ∫d, ∫c over two meshes) | pass: 1.5 % | pass: 1.5 % | pass: 0.89 % |
| R1 history never decreases | pass | pass | pass |
| R2 hidden state (U and SVARS differences) | pass: 0, 0 | pass: 0, 0 | pass: 0, 0 |
| R3 DTIME = 0 | pass: PNEWDT 0.5 | pass: PNEWDT 0.5 | pass: PNEWDT 0.5 |
| R4 large step | pass, 0 cutbacks | pass, 0 cutbacks | pass, 0 cutbacks |
| R5 damage near 1 (largest d) | pass: 0.985 | pass: 0.998 | pass: 0.999 |
| **Checks passed** | **8 of 10** | **8 of 10** | **10 of 10** |

Errors in C1 are those of the element's time-step limit against the key's, on the 32 × 32 mesh (max over the cycle ends,
divided by the key's maximum). All three elements were robust; A and B predicted too much damage. The cause is their
discretization: with equal-order bilinear fields on these meshes, the swelling strain turns into false stress. Run A's own
quadratic solution agrees with the key to 0.1 % on D0 over 5 cycles (`../reference/README.md`), and the damage errors of A
and B fall 2.9 to 3.7 times from the 16 × 16 to the 32 × 32 mesh (C2), toward the key. B2 removed the false stress with four condensed enhanced-strain modes per element.

### Points (rubric v1.1)

| Part | Points | A | B | B2 |
|---|---:|---:|---:|---:|
| Convergence (C1 per quantity 15, C2 10, C3 5) | 30 | 23.1 | 23.1 | 30 |
| History (R1 5, R2 5, time-step behaviour 5) | 15 | 15 | 15 | 15 |
| Robustness (R0, R3, R4, R5 4 each; cutbacks 4) | 20 | 20 | 20 | 20 |
| Design, judged from `design.md` against the run's evidence | 20 | 17 | 17 | 20 |
| Honesty, judged from `STATUS.md` against what was run | 10 | 9 | 10 | 10 |
| Effort: 2.5 × (best wall time / own) + 2.5 × (best tokens / own) | 5 | 5 | 0.67 | 0.90 |
| **Total** | **100** | **89** | **86** | **96** |

How the parts were fixed: the weights were drafted before the runs. The split of the convergence points and the rubric for
design and honesty were fixed after run A's first score and before run B's (rubric v1, 20 October 2026, 21:40 UTC). Rubric
v1.1 (22:30 UTC) changed one design item: the stress-gradient term's formulation (6 points) became interpolation order (3:
full credit for mixed order, meaning u one order higher than the fields fed by its gradient, for example through condensed
modes; 0 for equal-order bilinear) and the projection or consistency of the gradient term with a convergence test (3); the
number of integration points is not scored. The design and honesty points are the organizer's judgment; the reasons for each
deduction:

- **A, design 17:** equal-order bilinear fields (−3).
- **A, honesty 9:** it reports its Q2 damage numbers but does not state the coarse-mesh damage bias as a limitation (−1).
- **B, design 17:** equal-order bilinear fields (−3).
- **B, honesty 10:** its `STATUS.md` states the coarse-mesh damage gap itself ("structured n32 versus n64 still differs
  9.35%") and both of its procedural deviations.


## Practice results (D0, rerunnable)

`practice_D0/A`, `practice_D0/B` and `practice_D0/B2` hold `score_test1.json` and `score_test1.md` from

```bash
python score.py run   ../runs/<run> out/<run> --problem ../host/dev_problem/d0.json --key reference_D0 --jobs 5
python score.py score ../runs/<run> out/<run> --problem ../host/dev_problem/d0.json --key reference_D0
```

on the released host (public abaqus_ufl `691c4eb`), 6 November 2026. R5 is not run on D0 (its stress variant belongs to the
hidden problem).

| Check | A | B | B2 |
|---|---|---|---|
| R0 runs complete | pass | pass | pass |
| C1 key observables: d_max, ∫d, ∫c (limits 4 %, 4 %, 0.5 %) | **fail**: 11.5 %, 11.4 %, 0.25 % | **fail**: 18.2 %, 13.7 %, 0.33 % | pass: 0.90 %, 0.70 %, 0.0088 % |
| C1 key fields: u_x, u_y, c, d (limit 1 %) | **fail**: 0.23 %, 0.49 %, 0.23 %, 5.0 % | **fail**: 0.30 %, 0.65 %, 0.31 %, 6.2 % | pass: 0.021 %, 0.015 %, 0.013 %, 0.28 % |
| C2 converges to the key | pass | pass | pass |
| C3 distortion (largest of d_max, ∫d, ∫c over two meshes) | pass: 2.4 % | pass: 2.5 % | pass: 1.5 % |
| R1 history never decreases | pass | pass | pass |
| R2 hidden state (U and SVARS differences) | pass: 0, 0 | pass: 0, 0 | pass: 0, 0 |
| R3 DTIME = 0 | pass: PNEWDT 0.5 | pass: PNEWDT 0.5 | pass: PNEWDT 0.5 |
| R4 large step | pass, 0 cutbacks | pass, 0 cutbacks | pass, 0 cutbacks |
| R5 damage near 1 | not run | not run | not run |
| **Checks passed** | **7 of 9** | **7 of 9** | **9 of 9** |

The practice problem shows the same pattern as the hidden one: A and B are robust but predict too much damage on the
32 × 32 mesh; B2 passes every check that runs.

## The released host gives the official numbers

To check the port, run B2's element was scored again on the released host with the same scorer, on the hidden problem
(privately). The scorer's output (`score_test1.json`) is identical to the official one in every check and number; only the
workspace field differs, and the new `not_run` list is empty. With the three comparison runs in `../host/README.md` (88
exact checks), this is the evidence that the released host and scorer reproduce the official evaluation. A and B were not
rescored on the hidden problem; their official numbers are from the host the agents had.
