# Verification results

Generated from saved artifacts by `python verification/summarize.py`. Host observables use the host's own Q4 quadrature. Fine-grid comparisons are convergence studies, not an independent 2D answer key.

## Compiled element and transaction checks

- [element_tests.json](element_tests.json)
- [reference_tests.json](reference_tests.json)
- [rollback.json](rollback.json)

The element tests include analytic forces on a nonaffine element, positive-trace and deviatoric energy, rigid motion, local mass balance, finite zero-time rejection, invalid orientation, and H=1e22. The global rejected-step comparison is bitwise identical.

## Manufactured 2D continuum solution

Analytical body sources are differentiated from the continuum laws using Sympy. H is a frozen spatial field. All physical fields have analytical Dirichlet values; q,z have no prescribed values.

| mesh | n | L2 ux | L2 c | L2 d | L2 q | L2 grad q | L2 flux |
| --- | --- | --- | --- | --- | --- | --- | --- |
| structured | 8 | 4.78876e-05 | 0.000664561 | 0.00243215 | 0.000493149 | 0.0137994 | 0.00197757 |
| structured | 16 | 1.21661e-05 | 0.000167449 | 0.000620251 | 0.000120748 | 0.0069019 | 0.000987044 |
| structured | 32 | 3.05778e-06 | 4.19399e-05 | 0.000155847 | 2.99742e-05 | 0.00345143 | 0.000493344 |
| structured | 64 | 7.65915e-07 | 1.04893e-05 | 3.90111e-05 | 7.47406e-06 | 0.00172579 | 0.000246652 |
| distorted | 8 | 5.28411e-05 | 0.00147072 | 0.00266122 | 0.000882345 | 0.023943 | 0.00231521 |
| distorted | 16 | 1.57233e-05 | 0.000774628 | 0.000692103 | 0.000452235 | 0.0254238 | 0.00133923 |
| distorted | 32 | 4.25591e-06 | 0.000353228 | 0.000176036 | 0.000208375 | 0.0239027 | 0.000877506 |
| distorted | 64 | 9.1913e-07 | 0.000192734 | 4.42919e-05 | 0.000112349 | 0.0255691 | 0.000532212 |

The structured field errors approach order two. On randomly distorted meshes concentration/pressure convergence is slower. Isolated pressure gradient error does not decrease through n=64; total flux does. This limitation also appears with the tested mixed-volumetric-stress alternative.

## Independent coupled transient reference

The cell-centred FV solver analytically eliminates displacement under ux=0 and sigma_yy=0. It uses direct stress differences in the flux, and the original damage ODE, with no q or z. The UEL solves a narrow 2D strip on the unchanged host. All integrals are divided by strip width. Both use five cycles; G_c=1e-4 and alpha_T=5e-4.

| UEL rows | max abs dmax error | max abs D_int error | max abs C_int error | max abs E_el error | final c L2 | final d L2 |
| --- | --- | --- | --- | --- | --- | --- |
| 16 | 0.00251494 | 0.0027545 | 0.00151682 | 1.55103e-06 | 0.000846137 | 0.00277319 |
| 32 | 0.000603513 | 0.00113312 | 0.000616403 | 5.38547e-07 | 0.000332804 | 0.000952031 |
| 64 | 0.000699952 | 0.000506218 | 0.000268622 | 2.05973e-07 | 0.000140796 | 0.000392155 |
| 128 | 0.000459612 | 0.000238029 | 0.000123988 | 8.71405e-08 | 6.52859e-05 | 0.000181627 |

The FV reference reaches alpha=0.00067813256, above alpha_T. The two-dimensional tensile split is separately tested by the element patch.

| method | n | final dmax | final D_int | final C_int | final E_el |
| --- | --- | --- | --- | --- | --- |
| lag_independent_FV | 128 | 0.452335 | 0.29755 | 0.467933 | 5.69646e-05 |
| lag_independent_FV | 256 | 0.45239 | 0.297545 | 0.467928 | 5.69624e-05 |
| lag_independent_FV | 512 | 0.452404 | 0.297543 | 0.467927 | 5.69618e-05 |
| UEL_strip | 16 | 0.45093 | 0.299873 | 0.467669 | 5.67642e-05 |
| UEL_strip | 32 | 0.452621 | 0.298402 | 0.467781 | 5.68541e-05 |
| UEL_strip | 64 | 0.452881 | 0.297896 | 0.467851 | 5.69061e-05 |
| UEL_strip | 128 | 0.452728 | 0.297699 | 0.467887 | 5.69331e-05 |

| method | steps/cycle | final dmax | final D_int | final C_int | final E_el |
| --- | --- | --- | --- | --- | --- |
| independent_FV | 40 | 0.45092 | 0.296077 | 0.467111 | 5.62871e-05 |
| independent_FV | 80 | 0.451736 | 0.297422 | 0.467884 | 5.69389e-05 |
| independent_FV | 160 | 0.452097 | 0.29811 | 0.468274 | 5.72675e-05 |
| independent_FV | 320 | 0.452269 | 0.298457 | 0.46847 | 5.74328e-05 |
| lag_independent_FV | 40 | 0.452389 | 0.296342 | 0.467203 | 5.6336e-05 |
| lag_independent_FV | 80 | 0.45239 | 0.297545 | 0.467928 | 5.69624e-05 |
| lag_independent_FV | 160 | 0.452406 | 0.298168 | 0.468295 | 5.72788e-05 |
| lag_independent_FV | 320 | 0.452419 | 0.298485 | 0.46848 | 5.74384e-05 |
| UEL_strip | 40 | 0.452738 | 0.296499 | 0.467163 | 5.63072e-05 |
| UEL_strip | 80 | 0.452728 | 0.297699 | 0.467887 | 5.69331e-05 |
| UEL_strip | 160 | 0.452742 | 0.298322 | 0.468254 | 5.72493e-05 |
| UEL_strip | 320 | 0.452751 | 0.298639 | 0.46844 | 5.74088e-05 |

The spatial error table uses the same lagged-history split in the independent FV solver. The temporal table additionally includes the implicit-history FV solution: both schemes approach the same limit. Successive 80/160/320-step changes in the integrated observables roughly halve. A remaining finite spatial error is visible.

## Uniform fatigue time benchmark

Uniform concentration is prescribed everywhere, ux is constrained and sigma_yy=0. This isolates constitutive/history/damage integration; it is not a diffusion solve. An independent scalar bracketed root solves the implicit damage/history law at 5120 steps/cycle, with a 2560-step reference self-check. There is no spatial discretization error in the homogeneous Q1 element.

| steps/cycle | max abs d error | max abs H error | max abs alpha error | final d | final alpha |
| --- | --- | --- | --- | --- | --- |
| 40 | 0.238954 | 54.5094 | 2.64298e-06 | 0.993816 | 5.45773e-06 |
| 80 | 0.117316 | 18.9883 | 9.00056e-07 | 0.992344 | 3.92333e-06 |
| 160 | 0.0571817 | 8.60063 | 4.05272e-07 | 0.991824 | 3.51804e-06 |
| 320 | 0.0280133 | 4.13223 | 1.93888e-07 | 0.99159 | 3.35346e-06 |
| 640 | 0.0139192 | 2.02799 | 9.46757e-08 | 0.991478 | 3.27814e-06 |

The full error histories and reference resolution check are in `uniform_history/summary.json`. Early damage peaks can require many steps to resolve even though Newton itself converges.

## D0 on the original host CLI

Final formulation: lumped damage source, direct stress, 3x3 history points. All listed runs use the original host/run_problem.py. Each full run is 20 cycles unless explicitly marked five.

| mesh | n | steps/cycle | cycles reached | status | cutbacks |
| --- | --- | --- | --- | --- | --- |
| structured | 8 | 40 | 20 | ok | 0 |
| structured | 16 | 40 | 20 | ok | 0 |
| structured | 32 | 40 | 20 | ok | 0 |
| structured | 64 | 40 | 5 | ok | 0 |
| distorted | 8 | 40 | 20 | ok | 0 |
| distorted | 16 | 40 | 20 | ok | 0 |
| distorted | 32 | 40 | 20 | ok | 0 |
| structured | 16 | 80 | 20 | ok | 0 |
| structured | 16 | 160 | 20 | ok | 0 |
| distorted | 16 | 80 | 20 | ok | 0 |
| distorted | 16 | 160 | 20 | ok | 0 |

At the end of cycle five:

| mesh | n | steps/cycle | d_max | D_int | C_int | E_el |
| --- | --- | --- | --- | --- | --- | --- |
| structured | 8 | 40 | 0.392487 | 0.163629 | 0.52604 | 1.36632e-06 |
| structured | 16 | 40 | 0.250576 | 0.0922006 | 0.513934 | 7.22538e-07 |
| structured | 32 | 40 | 0.201954 | 0.0717942 | 0.510253 | 5.22032e-07 |
| structured | 64 | 40 | 0.186747 | 0.0664073 | 0.509317 | 4.69154e-07 |
| distorted | 8 | 40 | 0.434605 | 0.17674 | 0.527967 | 1.4507e-06 |
| distorted | 16 | 40 | 0.258727 | 0.0973481 | 0.514767 | 7.65767e-07 |
| distorted | 32 | 40 | 0.205228 | 0.0731475 | 0.510488 | 5.3459e-07 |
| structured | 16 | 80 | 0.256282 | 0.0947201 | 0.515726 | 7.25597e-07 |
| structured | 16 | 160 | 0.259108 | 0.0959803 | 0.516637 | 7.27092e-07 |
| distorted | 16 | 80 | 0.26509 | 0.100087 | 0.516591 | 7.69778e-07 |
| distorted | 16 | 160 | 0.268128 | 0.101453 | 0.517514 | 7.71703e-07 |

At the end of cycle twenty:

| mesh | n | steps/cycle | d_max | D_int | C_int | E_el |
| --- | --- | --- | --- | --- | --- | --- |
| structured | 8 | 40 | 0.514202 | 0.175698 | 0.530063 | 1.37747e-06 |
| structured | 16 | 40 | 0.349495 | 0.0974245 | 0.515779 | 7.44548e-07 |
| structured | 32 | 40 | 0.283021 | 0.0755862 | 0.511634 | 5.38879e-07 |
| distorted | 8 | 40 | 0.598531 | 0.191101 | 0.532915 | 1.43884e-06 |
| distorted | 16 | 40 | 0.365222 | 0.10329 | 0.516876 | 7.89698e-07 |
| distorted | 32 | 40 | 0.288944 | 0.077119 | 0.511939 | 5.52056e-07 |
| structured | 16 | 80 | 0.368056 | 0.100552 | 0.517825 | 7.50914e-07 |
| structured | 16 | 160 | 0.378025 | 0.102144 | 0.518876 | 7.54214e-07 |
| distorted | 16 | 80 | 0.387064 | 0.106778 | 0.519013 | 7.97434e-07 |
| distorted | 16 | 160 | 0.399178 | 0.108571 | 0.520117 | 8.01499e-07 |

The smallest damage-area indicator A_damage can be zero throughout a fine run while coarse meshes exceed d=0.5. That threshold observable is particularly mesh sensitive; the full host CSVs retain it. These tables do not establish an independent two-dimensional reference error or spatial asymptotic rate for the long-cycle problem.

The earlier `final_D0_*` and `D0_lumped_s8_spc40` directories belong to the fully implicit candidate before the final pivot. Their names are retained for audit continuity; they are not the selected time integration. Their archived source hash is fcd4c37a1a12.

Q1 field differences at common physical quadrature points, relative to the selected structured n=64 solution at cycle five:

| mesh | n | L2 u | L2 c | L2 d | L2 q |
| --- | --- | --- | --- | --- | --- |
| structured | 8 | 0.000462417 | 0.0180357 | 0.120402 | 0.000440062 |
| structured | 16 | 0.000130282 | 0.00492374 | 0.0322698 | 0.000175761 |
| structured | 32 | 2.67964e-05 | 0.000999538 | 0.00676574 | 5.14149e-05 |
| distorted | 8 | 0.000516853 | 0.0201582 | 0.135846 | 0.000528214 |
| distorted | 16 | 0.000154382 | 0.00582496 | 0.0385505 | 0.000290794 |
| distorted | 32 | 3.35855e-05 | 0.00125292 | 0.00838699 | 0.000121198 |

Field changes at cycle twenty for time refinement on n=16:

| mesh | coarse/fine steps | L2 u | L2 c | L2 d | L2 q | max delta H | max delta alpha |
| --- | --- | --- | --- | --- | --- | --- | --- |
| structured | 40/80 | 8.14996e-05 | 0.00233431 | 0.00382269 | 1.10659e-05 | 0.250286 | 5.1028e-06 |
| structured | 80/160 | 4.17983e-05 | 0.00119788 | 0.00196732 | 5.67386e-06 | 0.14433 | 2.27832e-06 |
| distorted | 40/80 | 8.2594e-05 | 0.00239692 | 0.0043797 | 1.1183e-05 | 0.319857 | 5.51343e-06 |
| distorted | 80/160 | 4.2424e-05 | 0.0012343 | 0.00229292 | 5.75171e-06 | 0.19338 | 2.48066e-06 |

## Solver robustness and known failures

Newton work is counted as sparse linear solves in the unchanged host one-block Newton path. Cutbacks use its ordinary transaction and limits.

| case | history | reaction | status | max free d | max H | cutbacks/rejections | mean Newton updates | max updates |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| D0_Gc_3.5e-07 | implicit | lumped | ok | 0.984759 | 155.025 | 0 | 2.52 | 5 |
| D0_Gc_3.5e-08 | implicit | lumped | FAILED | 0.998701 | 1816.1 | 11 | 9.25 | 25 |
| D0_Gc_3.5e-09 | implicit | lumped | ok | 0.999922 | 41657.1 | 0 | 2.515 | 7 |
| D0_statistics | lagged | lumped | ok | 0.349495 | 2.34483 | 0 | 2.02875 | 3 |
| constrained_notch | lagged | lumped | ok | 0.999993 | 981952 | 0 | 2.035 | 3 |
| D0_statistics | implicit | lumped | ok | 0.369808 | 2.63611 | 0 | 2.81 | 5 |
| large_step_10 | implicit | lumped | ok | 0.415706 | 3.64772 | 25 | 6.13091 | 25 |
| lag_D0_Gc_3.5e-08 | lagged | lumped | ok | 0.998895 | 1801.61 | 0 | 2.035 | 3 |
| lag_D0_s8_spc40 | lagged | lumped | ok | 0.514202 | 4.89021 | 0 | 2.03125 | 3 |
| lag_D0_s8_spc80 | lagged | lumped | ok | 0.528594 | 5.23515 | 0 | 2.00875 | 3 |
| lag_saturation | lagged | lumped | ok | 0.999993 | 981952 | 0 | 2.035 | 3 |
| lag_saturation_local | lagged | lumped | ok | 1 | 9.28106e+19 | 0 | 1.62 | 3 |
| large_step_10 | lagged | lumped | ok | 0.293615 | 1.82972 | 0 | 2.67 | 3 |
| reject_after | lagged | lumped | ok | 0.356921 | 1.54634 | 1 | 2.22892 | 3 |
| reject_before | lagged | lumped | ok | 0.356921 | 1.54634 | 1 | 2.22892 | 3 |
| saturated_fatigue | lagged | lumped | ok | 1 | 9.28106e+19 | 0 | 1.62 | 3 |
| saturation | implicit | lumped | FAILED | 1 | 11173 | 11 | 12.4545 | 23 |
| saturation | implicit | lumped | FAILED | 1 | 4.96478e+09 | 11 | 21.0667 | 25 |
| selected_full_integration | lagged | full | ok | 1.15981 | 72.4364 | 0 | 2.0575 | 3 |

Failed cases are retained with their accepted partial trajectories and per-attempt records. The two old notch failures predate the addition of a free-damage diagnostic, so their tabulated maximum may include the prescribed notch d=1. They do not demonstrate successful propagation. The free-damage toughness sweep has no prescribed damage.

The implicit extreme fatigue-floor case fails in its first nominal increment. The selected split completes that same physical problem with H near 1e20 and free d=1. A separate fixed-H=1e22 element call also passes.

All accepted H and alpha values in the reported runs are nondecreasing. Small nodal d decreases can remain at global-solve accuracy (the old implicit D0 statistics run records about 2.8e-7). No separate nodal irreversibility constraint or post-convergence clipping was imposed.

Actual Abaqus execution, hidden problems, arbitrary quadrilateral families, and universal Newton convergence are unverified.

## Separate higher-order 2D comparison

This self-authored reference uses affine biquadratic displacement, concentration and damage, 5x5 quadrature, direct displacement Hessians in grad(g*p), and the original consistent damage equation. It has no projected-pressure or damage-source auxiliary fields. It uses the same lagged history in time to isolate spatial differences. Its nine-node macroelements are reference-only; the delivered element retains the required four-node topology. The generic host is unchanged.

| Q2 macro n | steps/cycle | status | final dmax | final D_int | final C_int | final E_el |
| --- | --- | --- | --- | --- | --- | --- |
| 4 | 40 | ok | 0.196139 | 0.0686475 | 0.509535 | 4.81773e-07 |
| 8 | 40 | ok | 0.184937 | 0.0653972 | 0.509216 | 4.56914e-07 |
| 16 | 40 | ok | 0.178269 | 0.0644102 | 0.509017 | 4.51073e-07 |
| 16 | 80 | ok | 0.185217 | 0.066327 | 0.510779 | 4.50786e-07 |
| 16 | 160 | ok | 0.189475 | 0.0673365 | 0.511695 | 4.50787e-07 |
| 32 | 40 | ok | 0.177744 | 0.0642636 | 0.508995 | 4.50239e-07 |

The separate Q2 tangent has relative FD error 6.51e-11; the quadratic-displacement stress-flux residual error is 1.14e-17.

Cycle-five primary fields compared at common physical quadrature points against Q2 n=32:

| Q1 mesh | n | L2 u | L2 c | L2 d | L2 flux | L2 grad q |
| --- | --- | --- | --- | --- | --- | --- |
| structured | 8 | 0.000471883 | 0.0183845 | 0.123104 | 0.00879946 | 0.00839719 |
| structured | 16 | 0.000139853 | 0.00527597 | 0.0350127 | 0.00248212 | 0.00671669 |
| structured | 32 | 3.63835e-05 | 0.00135363 | 0.00954275 | 0.000790481 | 0.004875 |
| structured | 64 | 9.60026e-06 | 0.000357609 | 0.0028055 | 0.000353792 | 0.003439 |
| distorted | 8 | 0.00052628 | 0.0205066 | 0.138551 | 0.0099568 | 0.0117865 |
| distorted | 16 | 0.000163937 | 0.00617701 | 0.0412812 | 0.00289395 | 0.0147639 |
| distorted | 32 | 4.31696e-05 | 0.00160599 | 0.011156 | 0.000901838 | 0.0138368 |

Maximum observable differences over the first five cycles:

| Q1 mesh | n | dmax | D_int | C_int | E_el |
| --- | --- | --- | --- | --- | --- |
| structured | 8 | 0.214742 | 0.0993651 | 0.0225419 | 1.69533e-06 |
| structured | 16 | 0.0728319 | 0.027937 | 0.00590405 | 4.7942e-07 |
| structured | 32 | 0.0242101 | 0.00753057 | 0.00150805 | 1.24177e-07 |
| structured | 64 | 0.009003 | 0.00214365 | 0.000397268 | 3.24399e-08 |
| distorted | 8 | 0.256861 | 0.112476 | 0.0265986 | 1.88767e-06 |
| distorted | 16 | 0.0809825 | 0.0330845 | 0.00694513 | 5.59187e-07 |
| distorted | 32 | 0.0274835 | 0.00888386 | 0.00178904 | 1.46831e-07 |

Q2 reference family self-convergence against its n=32 mesh:

| Q2 macro n | L2 u | L2 c | L2 d | max dmax diff | max D_int diff | max C_int diff | max E_el diff |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 4 | 2.43068e-05 | 0.000721797 | 0.00546194 | 0.018395 | 0.00438383 | 0.000780815 | 6.70228e-08 |
| 8 | 7.55109e-06 | 0.000247606 | 0.0016477 | 0.00719292 | 0.00113359 | 0.000258843 | 1.39028e-08 |
| 16 | 9.74092e-07 | 3.04678e-05 | 0.000205885 | 0.00105119 | 0.000180011 | 3.91876e-05 | 2.31018e-09 |

This is an independently discretized numerical comparison from the stated equations, with its own finite reference error. It is not the organizer's hidden reference or an external Abaqus result.

