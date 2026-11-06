# Organizer note (read this first, then `PUBLIC_TASK.md` and `host/README.md`)

You are a test agent in a controlled trial. The challenge is in `PUBLIC_TASK.md`. These conditions
are pinned by the organizer and are identical for every run of this challenge, except the arm line at
the end.

## Conditions

- **Visibility:** this workspace only. No internet, no package installs. Do not read files outside
  this workspace, except the Python environment and the system compilers named below. Do not search
  the disk for other copies of this challenge, answer keys or other agents' work. Any such access
  invalidates the run. Your commands are logged and audited.
- **Temporary files:** only under `tmp/` in this workspace (TMPDIR points there). Do not write to
  /tmp or your home directory.
- **Environment:** use `<env>`. Put its `bin/` first on PATH:
  `export PATH=<env>/bin:$PATH`.
  It has Python 3.13, numpy, scipy, sympy, gfortran, meson and ninja (for `numpy.f2py`).
- **The host** in `host/` is part of the challenge. Use it as it is; do not modify it. Your element is
  scored on an unmodified copy.
- **Abaqus is not available.** Do not claim Abaqus execution.
- **No human help.** Nobody will answer questions. Make reasonable decisions and record them.
- **Time budget:** about 6 hours of wall-clock time. Stop when done or when the budget is used, and
  report what you actually reached.

## Deliverables

As listed in `PUBLIC_TASK.md` §5: `element/` (source and `uel_contract.json`), `design.md`,
`verification/`, `research_log.md` and `STATUS.md`. Do not hard-code answers: the organizer runs your
compiled element on hidden problems on the same host and compares with an independent reference.

## Your arm: B (with abaqus_ufl as a tool)

You may use `tools/abaqus_ufl/`, the `abaqus_ufl` code generator at the pinned public commit a02ed7c
(18 October 2026), with its documentation and its AI skill (`ai_skills/abaqus-ufl/`). Import it by
putting `tools/abaqus_ufl` on `PYTHONPATH`. Use it to create the UEL where it fits. Where the task
needs something it does not support, extend it or write that part yourself, and document it. Record
what you used, every modification, and how much of the Fortran was generated versus hand-written. Do
not imply capabilities it does not have.
