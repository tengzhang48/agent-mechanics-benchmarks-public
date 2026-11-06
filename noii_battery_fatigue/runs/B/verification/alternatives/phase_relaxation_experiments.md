# Archived phase-relaxation experiments

2026-10-20T21:07:23.423383+00:00

The clean bounded-psi candidate was retained. All tested relaxation variants failed the strongly heterogeneous spatial fatigue-floor case: notch_sides with alpha_T=1e-10, distorted n16 seed3, two cycles and80 steps per cycle. The physical floor stayed1e-15 and k_fg=.5. No main, host, or package files were changed by these experiments.

Each variant passed independent assembled residual/state checks and all24 finite-difference Jacobian columns. The independent compiled oracles also checked the exact added gamma terms. Every failed global run preserved its accepted state and its actual failed trial inputs. A passing short restart did not predict whole-path robustness.

| Relaxation gamma | Coefficient | Failed nominal step | Failed trial time | Wall seconds |
|---|---:|---:|---:|---:|
| eta*Hold | 0.1 | 8 | 0.720410156 | 82.23 |
| eta*Hold | 0.5 | 14 | 1.39423828 | 200.17 |
| eta*Hold | 1.0 | 8 | 0.70859375 | 140.85 |
| eta*(1+min(aHold,1000)) | 0.5 | 7 | 0.670507812 | 102.75 |
| eta*(1+min(aHold,1000)) | 1.0 | 7 | 0.64609375 | 97.75 |
| eta*(1+min(aHold,1000)) | 2.0 | 13 | 1.28642578 | 385.89 |
| eta*(1+min(aHold,1000)) | 5.0 | 3 | 0.225097656 | 38.55 |
| eta*(1+min(aHphase,1000)) | 0.5 | 14 | 1.34404297 | 88.35 |
| eta*(1+min(aHphase,1000)) | 1.0 | 14 | 1.34589844 | 110.56 |
| eta*(1+min(aHphase,1000)) | 2.0 | 3 | 0.234472656 | 65.53 |
| C/sqrt(dt)*(1+min(aHphase,1000)) | 0.1 | 14 | 1.37851562 | 347.23 |
| C/sqrt(dt)*(1+min(aHphase,1000)) | 0.3 | 3 | 0.235742187 | 73.19 |
| C/sqrt(dt)*(1+min(aHphase,1000)) | 1.0 | 2 | 0.157714844 | 43.49 |

The frozen-H eta=.5 variant nevertheless passed uniform fatigue-floor2cycles at20/80/320 steps per cycle against the independent scalar solution. All four heterogeneous fixed-H cases also passed, with residual<=2.89e-14 and intact-fraction oracle relative error<=1.01e-15. The accepted-full-stiffness eta=.5 variant passed D0 n16 fivecycles and notch n32 distortedseed3 onecycle with zero cutbacks. These successes were insufficient to justify adding relaxation.

Exact source hashes, modification UTC, generation counts, and reference changes are in each candidate verification/frozen_source.json, verification/generation.json, and check_element.py. Full global diagnostics are in phase_relaxation_global_runs.json and each individual run folder. The clean candidate remains tmp/lumped_psilag_candidate, generated SHA90b72e74a82b697696bbd5ab1ee5f1c6fea2f87d7789bb554dbcb7a2066f88b0.
