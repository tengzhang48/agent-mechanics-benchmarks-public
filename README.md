# MechanicsHarness: agent mechanics benchmarks

MechanicsHarness evaluates AI agents on mechanics problems with reproducible tasks, independent reference solutions and
documented human verification. Website: https://mechanicsharness.org

This repository holds the released benchmarks: for each, the task as the agents received it, the host that runs a submitted
element, the scorer, the reference solution, the agents' submissions and the results, with a practice problem that anyone can
rerun from a clean checkout.

| Benchmark | What the agent writes | Status |
|---|---|---|
| [`noii_battery_fatigue/`](noii_battery_fatigue/) | An Abaqus/Standard user element (4-node quadrilateral; displacement, lithium concentration and damage) for the battery fatigue model of Noii et al. (2024), *Comput. Methods Appl. Mech. Eng.* 428, 117068, and its FEniCS implementation, *Eng. Fract. Mech.* 306, 110251 | Released: task, host, scorer, reference code, three runs (20 to 21 October 2026), results, case study |
| Rock fracture | An Abaqus/Standard user element (8-node brick; displacement and damage) for the gradient-damage rock model of L. Anand (2026), *Int. J. Solids Struct.* 325, 113739, tested on the one-element compression test of its Fig. 17 | In preparation, under review by the source-paper author |

## Start here

```bash
git clone https://github.com/tengzhang48/agent-mechanics-benchmarks-public
cd agent-mechanics-benchmarks-public/noii_battery_fatigue
```

Then follow `noii_battery_fatigue/README.md`: install the public abaqus_ufl and gfortran, and score one of the saved
elements on the practice problem.

## In this repository

- `EVALUATION_POLICY.md`: hidden tests, traceable verdicts, agent and evaluator errors, human verification, source papers.
- `CONTRIBUTING.md`: how to report a problem or propose a case.
- `CLAIMS_AND_EVIDENCE.md`: each statement in this README and its evidence.
- `RELEASE_NOTES.md`, `CITATION.cff`, `SHA256SUMS` (hashes of every file here).
- `LICENSE` (MIT, for code) and `LICENSE-CONTENT` (CC BY 4.0, for text, figures and data).

## Not in this repository

- The source papers and any files their authors shared with us. Papers are cited and linked.
- Tests reserved for future blind rounds and their reference results.
- Agent logins, home-folder state and the agents' session logs (kept, and released only after review).
- Internal correspondence and working notes.

## Licence

Code (Python, Fortran, shell scripts) is under the MIT License (`LICENSE`); text, figures and data are under CC BY 4.0
(`LICENSE-CONTENT`). Copyright 2026 Teng Zhang. Files written by AI agents during the benchmark runs are in each benchmark's
`runs/` folder and are labelled there. The source papers are cited, not included; no code, mesh, figure or PDF from their
authors is in this repository.
