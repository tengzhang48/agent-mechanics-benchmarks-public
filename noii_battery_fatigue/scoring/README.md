# Scoring

`score.py` and `battery.py` are the scripts that produced the official scores, with the changes listed at the top of each
file and no change to any check or tolerance. The official scores came from a hidden problem of the same kind as the
practice problem D0, which stays private (see `../README.md`). Here the same checks run on D0, against its reference.

## Run it

With the host set up (`../host/README.md`), from this folder:

```bash
python score.py run   ../runs/B2 out/B2 --problem ../host/dev_problem/d0.json --key reference_D0 --jobs 8
python score.py score ../runs/B2 out/B2 --problem ../host/dev_problem/d0.json --key reference_D0
```

`run` makes 7 host runs and 2 robustness tests, each in its own process (about 15 minutes for B2 on 8 cores; B's element
is about 3 times slower per call than A's); `score` writes `out/B2/score_test1.json` and `score_test1.md`. To score your own element,
point both commands at a folder with your `element/` and `element/uel_contract.json` (format in
`../task/host_README_as_given.md` §2).

## The checks

All on the first 10 cycles of the problem; "the key" is the reference solution (`reference_D0/`). The element's history
update is its own choice, so both are compared in the time-step limit, v = 2 v(T/80) − v(T/40), at the cycle ends.

| Check | Passes when |
|---|---|
| R0 runs complete | all 7 host runs finish every step |
| C1 key observables | on the 32 × 32 mesh, d_max and ∫d within 4 %, ∫c within 0.5 % of the key (max over cycle ends / max of the key) |
| C1 key fields | on the 32 × 32 mesh, the RMS difference of the nodal u_x, u_y, c and d within 1 % of the key's maximum |
| C2 converges to key | each of these errors drops at least 1.5 times from n = 16 to n = 32, or is already below 0.5 % |
| C3 distortion | d_max, ∫d and ∫c on two distorted 32 × 32 meshes agree with the structured mesh within the C1 tolerances |
| R1 history | H and α are declared in the contract and never decrease between accepted increments |
| R2 hidden state | two runs that reject the same increment, one after the element saw the trial states and one before, agree in U and SVARS to 1e-12 |
| R3 DTIME = 0 | one direct call with DTIME = 0 returns finite arrays or PNEWDT < 1 |
| R4 large step | a run with a quarter of the steps per cycle completes |
| R5 damage near 1 | a stress variant of the problem completes and d reaches 0.9 within [−0.01, 1.01] |

The elastic energy E_el is reported, not scored: a consistent 4-node element can be 8 % off in it on these meshes.

R5 needs the stress variant of the hidden problem, which is not released, so on D0 it is reported as not run. The
point scores of the official results also include the organizer's judgment of each run's design and honesty and its
effort; see `../results/README.md`.

## The reference subset

`reference_D0/F1_u64_div40` and `F1_u64_div80` hold what `score.py` reads from the reference: the 20 cycle ends of
`observables.csv` (columns t, D_int, C_int, E_el, d_max_host16, d_max_host32) and of the nodal values on the 16 × 16 and
32 × 32 host meshes (`host_nodes.npz`), with the values unchanged. The reference code is in `../reference/`; the full output
(944 MB) is not in this repository.
