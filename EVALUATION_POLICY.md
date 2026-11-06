# Evaluation policy

This policy follows the MechanicsHarness release plan (Teng Zhang, November 2026). The same text is on
https://mechanicsharness.org/evaluation.html.

## Hidden tests

Tests reserved for future blind evaluation, and their reference results, stay private. We publish the properties they test
and this policy. Completed historical tests may be released as evaluation cases; once released, they are no longer hidden, and
later rounds use new private cases.

## Traceable verdicts

Every published verdict will identify the task, the agent and model version, the tools, the resource budget, the submission,
the host and scorer versions, and the reference data. Where a fix to the host or the scorer affects a result, the saved
submissions are rescored on a common corrected evaluator, the original results are kept, and evaluation changes made after the
runs are labelled. A featured result may appear on the website's home page before its case page; it then states the
versions of the key and checker it rests on, and links the case page once that is published.

## Agent errors and evaluator errors

Case studies distinguish defects of the agent, defects of our own host or scorer, and outcomes that remain unresolved. The runs
form a small collection and are described as case studies; we avoid broad rankings across different task conditions.

## Human verification

Each case study states who specified the physics, who interpreted the source paper, who wrote or checked the reference, who
reviewed discrepancies, and who made each final judgment.

## Source papers and their authors

Source papers are cited and linked. For the research-model benchmarks (level 2), results are shared with the authors of the
source paper before publication. At every level, factual corrections from the source authors are recorded. Third-party
material is included only where its redistribution terms permit, and files that authors share with us privately stay private.

## Agent records

Session records are kept for every run. Transcripts or excerpts are released only after review for credentials, personal
information and private correspondence.
