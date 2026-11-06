# The task, as the agents received it

Each run started in its own workspace with these files and an empty `element/`, `verification/` and `tmp/`
(written by the organizer's build scripts on 20 October 2026):

| File here | In the workspace | What it is |
|---|---|---|
| `PUBLIC_TASK.md` | `PUBLIC_TASK.md` | The challenge: the model, what is the agent's to decide, the host, how the element is judged, the deliverables. Identical in all three runs. |
| `ORGANIZER_NOTE_A.md`, `ORGANIZER_NOTE_B.md`, `ORGANIZER_NOTE_B2.md` | `ORGANIZER_NOTE.md` | The run conditions (isolation, environment, time budget, deliverables) and the arm line. |
| `host_README_as_given.md` | `host/README.md` | How the host calls the element, the contract file, the practice problem D0, and what the scorer adds. |

The workspaces also held the host code (`host/`, the solver part of abaqus_ufl.fe at a development revision of 20 October
2026), and, for arms B and B2, a copy of the public abaqus_ufl at commit `a02ed7c` (18 October 2026) under
`tools/abaqus_ufl/`. Neither copy is repeated here: the released host in `../host/` runs the same calculation on the public
abaqus_ufl (see `../host/README.md`), and `a02ed7c` is in the public abaqus_ufl repository.

## How each run was started

Each run was one non-interactive Codex session (OpenAI Codex CLI 0.162.1, model gpt-6.1-sol, reasoning effort max),
started in its workspace with a 7-hour hard limit (the stated budget was about 6 hours), with TMPDIR set to the workspace's
`tmp/`:

```
codex exec --skip-git-repo-check -s workspace-write -C <workspace> -o <final message file> "<prompt>"
```

Run B2 was started the same way with sub-agents switched off (`codex exec --disable multi_agent ...`). The prompt, the
same for all three runs:

> Read ORGANIZER_NOTE.md in this folder first, then PUBLIC_TASK.md and host/README.md, and carry out the challenge under
> the organizer's conditions. Work only inside this folder. End with a short summary of what you reached and what remains
> unverified.

Launch times: A and B 2026-10-20 19:38:37 UTC; B2 2026-10-20 23:09:55 UTC.

## What was changed for this release

Only placeholders replace machine paths (`<env>` for the Python environment, `<workspace>` for the run's own folder) and
one private repository name in `host_README_as_given.md`; every replacement is counted in `../RELEASE_RECORD.md`. The
wording is otherwise unchanged, including the sentences about the hidden problems, which the scorer used and which stay
private (see `../README.md`).
