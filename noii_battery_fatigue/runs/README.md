# The three runs (20 to 21 October 2026)

The same agent built the element three times: OpenAI Codex CLI 0.162.1, model gpt-6.1-sol, maximum reasoning effort, the
same task, conditions and budget (about 6 hours) for each. Everything in these folders was written by the agent during its
run, except the placeholders listed in `../RELEASE_RECORD.md`.

| Run | Arm | Element final | Whole run | Tokens (whole run) |
|---|---|---|---|---|
| `A` | without abaqus_ufl | 35 min (20:14 UTC) | 61 min | 10.2 M |
| `B` | with abaqus_ufl (public `a02ed7c`) and its AI skill, which advises a team of sub-agents | 87 min (21:06 UTC) | 301 min, then stopped while still running extra studies | 158.2 M |
| `B2` | B rerun as a single agent (sub-agents switched off) | 74 min (00:23 UTC) | 262 min | 81.8 M |

Times are from the launch; tokens include sub-agents and cached input, from the Codex session files.

## What each run built

- **A:** a hand-written 4-node element with six bilinear nodal fields (u_x, u_y, c, d, q = g p_eff, and z = (1 − κ)H(1 − d)
  for the damage equation near d = 1), 3 × 3 Gauss points and a lumped damage reaction, and its own quadratic (Q2)
  reference to check convergence.
- **B:** an element generated with abaqus_ufl (`a02ed7c`): six bilinear nodal fields (u_x, u_y, c, d, a pressure
  projection q and an intact fraction), a lumped damage reaction and a bounded lag of the tensile-energy numerator in the
  damage equation, chosen after studies of many alternatives (`verification/alternatives/`). In its `STATUS.md`, B
  declared its own run invalid: one broad process listing it ran (`ps`) showed unrelated jobs on the machine, against
  the workspace-only condition. The organizer checked the visible listing: it named no organizer file, hidden problem,
  answer key or other run, and A and B2 never listed processes. The organizer recommended scoring the run with this
  deviation recorded, and Teng Zhang ruled so on 6 November 2026; its results are reported on that basis.
- **B2:** an element generated with abaqus_ufl plus 223 hand-written lines (3,666 of its 3,889 Fortran lines are generated):
  six bilinear nodal fields (u_x, u_y, c, d, a pressure p and an intactness s) and four statically condensed
  enhanced-strain modes inside each element, which remove the false swelling stress of the bilinear element.
  It found the problem with its own higher-order reference, which showed three times too much damage before the fix.

## In each folder

| Path | What it is |
|---|---|
| `element/` | The submitted UEL sources and `uel_contract.json`, byte-identical to the element that was scored |
| `design.md` | The formulation, the alternatives considered, and the reasons |
| `STATUS.md` | What works, what does not, and what is unverified, in the agent's words |
| `research_log.md` | Times, failed attempts and pivots |
| `verification/` | The agent's own tests and studies, without NumPy field arrays and files over 1 MB |
| `DATA_MANIFEST.tsv` | The files left out of `verification/` (path, size, SHA-256), for a later data deposit |
| `env.sh` (B2 only) | B2's environment script |

Not here: the agents' session logs (to be reviewed before any release), their `tmp/` folders, build products, and the
abaqus_ufl copy given to B and B2 (public commit `a02ed7c`).
