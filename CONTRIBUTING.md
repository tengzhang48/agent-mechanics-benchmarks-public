# Contributing

MechanicsHarness is maintained by Teng Zhang (Syracuse University). The website is https://mechanicsharness.org.

## Report a problem

Open an issue in this repository. Say which benchmark and which file, what you ran (command, operating system, Python,
NumPy and gfortran versions), what you expected, and what happened. A wrong statement in a case study is a problem too;
please point to the sentence and to the evidence.

## Have an element scored on a hidden problem

For the battery benchmark (`noii_battery_fatigue/`), you can ask for your own element to be scored on its hidden problem.
First score it on the practice problem (`noii_battery_fatigue/README.md`, "Rerun it"), then open an issue that links your
element folder (`element/` with `uel_contract.json`) and that score. We run the same scorer privately and report the result of
each check, but not the hidden problem's setup or its reference values.

## Propose a case

We plan to accept proposed cases once the evaluation workflow has been exercised on a contributed submission. Until then,
open an issue with:

- the published model and the paper (citation and DOI) and any public code from its authors;
- the test you have in mind and where its conditions are stated;
- whether the source material may be redistributed (we cite and link papers; we include third-party files only where their
  terms permit).

## Rules that protect the evaluation

- Do not post the content of hidden problems, or reference results for them, in issues or pull requests.
- Released historical tests are no longer blind; future rounds use new private cases (see `EVALUATION_POLICY.md`).
- Credit the reference work exactly: what was taken from a paper or its code, and what is new.
