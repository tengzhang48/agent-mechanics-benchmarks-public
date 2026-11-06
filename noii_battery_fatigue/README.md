# Battery fatigue: a chemo-mechanical user element (Noii et al. 2024)

An AI coding agent writes an Abaqus/Standard user element (UEL) for a published model of lithium diffusion with swelling and
fatigue phase-field fracture: a 4-node quadrilateral with displacement, lithium concentration and damage. The model is
fixed; the formulation is the agent's choice. The element is run on a host that plays the role of Abaqus/Standard and is
compared with an independent reference solution.

The same agent (OpenAI Codex, gpt-6.1-sol) built the element three times, on 20 to 21 October 2026: without our
code-generation framework abaqus_ufl (A), with it and a team of sub-agents (B), and with it as a single agent (B2). All
three elements were robust. A and B predicted 8.5 % and 9.3 % too much total damage on the 32 × 32 mesh (the limit is 4 %);
B2 was within 0.55 %.

| Run | Arm | Official checks passed | Points (of 100) | Element final after | Tokens |
|---|---|---|---:|---|---|
| A | without abaqus_ufl | 8 of 10 (damage accuracy failed) | 89 | 35 min | 10.2 M |
| B | abaqus_ufl and sub-agents | 8 of 10 (damage accuracy failed) | 86 | 87 min | 158.2 M |
| B2 | abaqus_ufl, single agent | 10 of 10 | 96 | 74 min | 81.8 M |

These official scores were measured on a hidden problem, which stays private (below). On the public practice problem D0,
rescored on the released host on 6 November 2026, the pattern is the same: A and B fail damage accuracy, and B2 passes
every check that runs. Details, the rubric and both sets of results: `results/README.md`. The case study written for the source paper's authors,
without the hidden problem: `case_study.html`.

## Rerun it

You need Python 3.10 or newer and gfortran. Install the public abaqus_ufl with its finite-element extras (NumPy, SciPy,
Meson, Ninja) at the commit this folder was checked with:

```bash
pip install "abaqus_ufl[fe] @ git+https://github.com/tengzhang48/abaqus_ufl@691c4eb"
```

Then score a saved element on the practice problem D0 (about 15 minutes for B2 on 8 cores):

```bash
cd scoring
python score.py run   ../runs/B2 out/B2 --problem ../host/dev_problem/d0.json --key reference_D0 --jobs 8
python score.py score ../runs/B2 out/B2 --problem ../host/dev_problem/d0.json --key reference_D0
```

`out/B2/score_test1.md` should match `results/practice_D0/B2/score_test1.md`. To test your own element, put it in a folder
with `element/` and `element/uel_contract.json` (format: `task/host_README_as_given.md` §2), and run the same two commands
on that folder.

## In this folder

| Folder | What it holds |
|---|---|
| `task/` | The task, the run conditions and the host description, as the agents received them; how each run was started |
| `host/` | The host on the public abaqus_ufl, the meshes, the observables and the practice problem D0 |
| `scoring/` | The scorer, the robustness tests, and the part of D0's reference the scorer reads |
| `reference/` | The reference solution's code (dolfinx) for D0, and its verification records |
| `runs/A`, `runs/B`, `runs/B2` | Each run's element, design note, status, research log and verification |
| `results/` | The official results, the points, and the practice results on D0 |
| `case_study.html` | The case study shared with the source paper's authors on 24 October 2026, without the hidden problem |
| `RELEASE_RECORD.md` | Every text replacement made for this release (machine paths and one private repository name) |

## What is not here, and why

- **The hidden problem** (with its stress variant and their reference results), on which the official scores were
  measured. It is kept private so that future agents can be scored on the same problem as these three runs. The practice
  problem, its reference and the scorer are here, so the method can be checked and rerun. To have an element scored on the
  hidden problem, open an issue.
- **The agents' session logs.** They are kept and will be released only after a review for credentials, personal
  information and private correspondence.
- **Large outputs:** NumPy field arrays and other files over 1 MB from the runs' verification folders (about 340 MB). Each
  run's `DATA_MANIFEST.tsv` lists them with their SHA-256.
- **The abaqus_ufl copy given to B and B2:** it is the public abaqus_ufl at commit `a02ed7c`.
- **The source papers and the authors' code.** They are cited below; none of their code, meshes, figures or PDFs is here.

## Source papers and credit

- The model: N. Noii, D. Milijasevic, H. Waisman, A. Khodadadian (2024), Fatigue failure theory for lithium diffusion
  induced fracture in lithium-ion battery electrode particles, *Comput. Methods Appl. Mech. Eng.* 428, 117068,
  doi:10.1016/j.cma.2024.117068.
- Its FEniCS implementation: N. Noii, D. Milijasevic, A. Khodadadian, T. Wick (2024), An efficient FEniCS implementation for
  coupling lithium-ion battery charge/discharge processes with fatigue phase-field fracture, *Eng. Fract. Mech.* 306,
  110251, doi:10.1016/j.engfracmech.2024.110251, and the authors' public code, github.com/noiiG/Battery-Li-Ion-Simulation-Fatigue
  (revision 8bcec27).
- **Taken from them:** the model's equations and its 2D conventions, which follow the authors' script (stated in
  `task/PUBLIC_TASK.md` §1). **Ours:** the parameters and problems, the task, the host and scorer, the reference solution,
  and the analysis. Where an element disagrees with the reference, that is the element's error, not an error in the papers.
- Prepared by Claude (organizer: task, host, scorer, reference solution and analysis) with Teng Zhang, who set the
  direction and made the decisions recorded in `results/README.md` and the case study. The files in `runs/` were written
  by the agents. The results were shared with the source paper's authors on 24 October 2026.

## Licence

Code under the MIT License, text, figures and data under CC BY 4.0 (see `../LICENSE` and `../LICENSE-CONTENT`). The files in
`runs/` were written by OpenAI Codex during runs set up by Teng Zhang, and are released under the same terms.
