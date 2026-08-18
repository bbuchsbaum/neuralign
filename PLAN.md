# PLAN: Stage 1 truth-preserving objects

Authoritative runbook for the current Agent. Mote ids in parentheses are the durable tracker (`nl-s1-*`).
Update checkboxes only after the listed tests pass. Do not start Stage 2.

See `GOAL.md` for the completion contract.

## Completion contract

The task is complete only when:

- [ ] Steps 1–7 below are complete
- [ ] `tests/testthat/test-correctness-contracts.R` Stage 1 cases pass
- [ ] `devtools::test(filter="aligned-study|alignment_result|apply_alignment|fit_alignment|quality|compose")` passes
- [ ] `scripts/verify-goal` exits 0
- [ ] no TODO/FIXME introduced by this work in touched files
- [ ] no new aligners or Stage 2/3 APIs were added
- [ ] implementation has been reviewed against every acceptance criterion in `GOAL.md`

The agent must not mark an item complete until it has verified it.

## Steps

### 1. Regression contracts first (`nl-s1-regression-tests`)

- [ ] Add `tests/testthat/test-correctness-contracts.R`
- [ ] Low-rank dims: `U` 7×2, `V` 5×2 → target 7, source 5
- [ ] Scaled inverse: several `s ≠ 1`, `A %*% inverse(A)` is identity; auto inverse never raw transpose
- [ ] Fold-space safety: data-driven CV cannot produce or stack a global aligned matrix
- [ ] Result/model identity: every retained fold model reproduces its stored assessment output
- [ ] Model validity: duplicate names, `NA`, `Inf`, character operators, heterogeneous codomains rejected
- [ ] Composition: composed apply equals sequential apply; unknown or unequal intermediate spaces error

Write failing tests first if the code is not yet fixed. Do not skip Stage 1 cases to go green.

### 2. Low-rank target dimension (`nl-s1-lowrank-dim`)

- [ ] `.transform_target_dim()` returns `nrow(U)`, not `ncol(U)`
- [ ] Audit other rank-vs-dimension helpers
- [ ] Regression test green

Files: `R/shared_feature_space.R`

### 3. Scaled Procrustes inverse (`nl-s1-scaled-procrustes`)

- [ ] Scale is retained on the fitted map (not stripped)
- [ ] Auto inverse uses `(1/s) Q^T` when `|s| ≠ 1`
- [ ] `restrict_to_identified` rejects scaled maps unless it handles scale
- [ ] Regression test green

Files: `R/aligner_procrustes.R`, `R/apply_alignment.R`, `R/fit_alignment.R`

### 4. Authoritative validator (`nl-s1-validator`)

- [ ] One validator used by constructors, S4 validity, provider returns, `add_transform()`, subset, compose, load, apply
- [ ] Every map shares one identified target coordinate system
- [ ] Reject duplicate/empty names, non-finite ops, unsupported classes, mixed target dims
- [ ] `AlignmentResult` checks subject/codomain agreement, or is a resample artifact
- [ ] Regression test green

Files: `R/alignment_model.R`, `R/alignment_result.R`, `R/validate_alignment_objects.R`, `R/serialize.R`

### 5. Hard-disable fold-specific global aligned data (`nl-s1-fold-hard-disable`)

- [ ] When `anchor_common` is not strictly `TRUE`, no global `aligned` matrix
- [ ] No cross-fold pairwise/reconstruction quality
- [ ] Primary return is `AlignedResampleSet`
- [ ] `get_aligned()`, `as_aligned_matrix()`, `alignment_quality()`, `compose_alignment()` hard-error
- [ ] Regression test green

Files: `R/fit_alignment_internal.R`, `R/alignment_result.R`, `R/quality.R`, `R/compose.R`

### 6. Evaluation vs deployment (`nl-s1-eval-vs-deploy`)

- [ ] CV evaluation outputs are not attached to a full-data refit model
- [ ] Optional deployment refit is a separate artifact
- [ ] Each retained fold model reproduces its stored assessment output
- [ ] Regression test green

Files: `R/fit_alignment_internal.R`, `R/obs_crossfit.R`, `R/aligned_resample_set.R`

### 7. Fail-closed extract / quality / compose / stack (`nl-s1-compose-guards`)

- [ ] Unknown or unequal intermediate space errors (expert override only, stamped)
- [ ] Partial subject drop requires explicit `allow_partial` / `subjects = intersect(...)`
- [ ] Metamorphic: `apply(T2 ∘ T1, X) = apply(T2, apply(T1, X))`
- [ ] Declared exact inverse satisfies `T^{-1}(TX) = X`
- [ ] Regression test green

Files: `R/compose.R`, `R/alignment_result.R`, `R/quality.R`

## Order

`1 → {2, 3, 4} → 5 → {6, 7}`

4 should land before 5 so mixed-codomain models are rejected rather than merely warned.

## Out of scope

Stage 2 (`nl-s2-*`), Stage 3 (`nl-s3-*`), Stage 4 (`nl-s4-*`).
