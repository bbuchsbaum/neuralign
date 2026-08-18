# Goal: Stage 1 truth-preserving objects

This file is the completion contract for the current harness.
It covers **Stage 1 only** (`nl-s1`). Stages 2–4 stay closed until this goal is met.

## What finished means

`AlignmentResult` and `AlignmentModel` tell the same truth the analysis-facing layer already states:

- one model means one target coordinate system
- aligned outputs were produced by the attached model, or the object is a resample/crossfit artifact
- fold-specific spaces never look like a global aligned matrix
- operator properties belong to the fitted map

No new aligners. No Stage 2 lifecycle work. No compact-operator rewrite.

## Completion contract

The Stage 1 goal is complete only when **all** of the following are true:

1. Every required checkbox in `PLAN.md` is checked.
2. `scripts/verify-goal` exits 0.
3. `tests/testthat/test-correctness-contracts.R` passes without skip/fail of Stage 1 cases.
4. Related suites still pass: `aligned-study`, `alignment_result`, `apply_alignment`, `fit_alignment`, `quality`, `compose`.
5. No Stage 1 ticket remains ready/open in mote (`nl-s1-*`), or PLAN.md records that mote is unavailable and the PLAN checkboxes plus tests are the source of truth.
6. The agent has not weakened tests or acceptance criteria to obtain a green run.

## Non-goals

- `AlignmentBundle`, lazy storage, `harmonize_shared_spaces`
- CONE re-enable or new aligners
- `OrthogonalUpdate` / compact Procrustes (Stage 3)
- Provider schema-only prepare, `calibrate_subject`, identity split (Stage 2)
- Freezing 1.0 or “fixing” hosted CI (Stage 4)

## Review gate

Human review is required for these contract diffs before Stage 2 or Cloud Stage 3 starts:

- `get_aligned()` / `as_aligned_matrix()` / `alignment_quality()` hard-error on fold-specific results
- authoritative model/result validator
- observation-CV evaluation artifact vs deployment refit split
