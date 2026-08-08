# neuralign Architecture

`neuralign` is an orchestration layer for cross-subject alignment. It centralizes
shared alignment plumbing (data containers, reference selection, CV scaffolding,
method registry, apply/compose/serialize) while keeping domain semantics in
downstream packages.

The goal is a small, stable core where packages can “bring their own features”
and reuse the same alignment machinery.

## Core invariants

### Data orientation + operators

- **Data matrices are `features × observations`** throughout the core.
- **Transforms are left-multiply operators** with shape `target_features × source_features`.
  - Applying a transform to subject `s` is `Q_s %*% X_s`.

The helper `procrustes_rotation(..., convention = "right")` exists for packages
that naturally work in an `observations × features` / right-multiply convention.

### Domain-agnostic

Core functions do not assume fMRI-specific concepts (“TR”, “run”, “HRF”, “MNI”).
Those belong in client packages (e.g., feature extraction, run splitting policy,
and object-specific transform application).

## Main objects

### `AlignmentData`

S4 container for multi-subject alignment inputs.

- `data`: named list of per-subject matrices
- `subjects`: subject ids
- `obs_labels`: optional column labels for sparse/partial overlap workflows
- `space`: optional space identifier (character or `gds_space`)
- optional guidance (see `R/guidance.R`)

Coercion is via `as_alignment_data()`.

### `AlignmentModel`

S4 object representing a fitted model.

- `transforms`: named list `subject → operator` (`target × source`)
- `reference` / `reference_data`: reference spec and realized reference
- `method`: registry key (e.g. `"procrustes"`, `"gw"`)
- `space_from` / `space_to`: space bookkeeping
- `provenance`: params/versions/hashes
- `method_state`: method-specific state for apply
- `train_subjects`: subjects used for training (subject-axis CV)

### `AlignmentResult`

S4 object representing aligned outputs.

- `model`: `AlignmentModel`
- `aligned`: named list of aligned matrices
- `quality`: summary diagnostics (optional)
- `cv_info`: cross-validation metadata (optional)

## Aligner registry

Aligners are registered at package load and accessed by name.

- `register_aligner(name, fit_fn, apply_fn, capabilities, ...)`
- `available_aligners()`, `aligner_capabilities()`

Aligner contract (`fit_fn`):

```
fit_fn <- function(data, reference, train_idx = NULL, ...) {
  list(
    transforms = list(subject_id = operator, ...),
    reference_data = <matrix/template>,
    space_from = <ANY>,
    space_to = <ANY>,
    method_state = list(...)
  )
}
```

The core orchestrator (`fit_alignment()`) is responsible for:

- validating inputs
- resolving references (fixed subject, template, data-driven)
- coordinating CV
- assembling `AlignmentModel` / `AlignmentResult`

## Core workflows

### Fit/apply

- `fit_alignment(data, method = ..., reference = ..., cv = ...)`
- `apply_alignment(model_or_result, new_data, ...)`

### Composition

- `compose_alignment(model1, model2)` (and `%*%` for models/results)

Composition is purely linear-algebraic: it does not introduce domain semantics.

### Reference selection

`select_reference()` provides standard reference selection policies (e.g.,
`"medoid"`, `"centroid"`), using distance metrics that can be label-aware via
`obs_labels`.

### Cross-validation

Two axes are supported:

- **Subject-axis CV** (LOSO / k-fold): `create_cv_folds()` + `fit_alignment(cv=...)`
- **Observation-axis folds** (run / blocked-time patterns): `create_obs_folds()`

Observation-axis CV orchestration lives in `fit_alignment()`’s observation CV
path and related utilities in `R/obs_crossfit.R` / `R/cv.R`.

## Feature blocks (portable correspondence signals)

Block utilities are a key abstraction for making downstream packages thin:

- `alignment_feature_block()`: (matrix + name + weight + feature_names + meta)
- `harmonize_feature_blocks()`: intersection/reorder across subjects
- `stack_feature_blocks()`: sqrt-weight scaling + row-bind
- `build_alignment_features()`: harmonize + stack + coverage/identifiability guardrails
- `block_alignment_report()`: unified coverage/rank/quality report

The block system is domain-agnostic. Client packages decide what blocks mean
(task betas, contrasts, ROI fingerprints, priors, etc.) and only supply matrices
and labels.

Optional metadata (not interpreted by core methods):

- `meta$source_type`: `"supervised"`, `"unsupervised"`, `"prior"`, …
- `meta$requires_independence`: `TRUE/FALSE` (used for advisory warnings in
  `build_alignment_features()` when `obs_crossfit=FALSE`)

## Guidance channels (optional)

`neuralign` includes optional guidance utilities that remain generic:

- ROI/parcel projectors (`R/roi_anchors.R`)
- intrinsic geometry coordinates (`R/intrinsic_geometry.R`)

These are inputs a method may use as priors; they are not “alignment by anatomy”
in themselves.

## Serialization + interop

- `save_alignment()` / `load_alignment()`
- `export_alignment()` / `import_alignment()`
- `as_map_family()` / `from_map_family()` for `fmrigds` integration

## Aligned representation layer (analysis-facing)

`AlignmentResult` is an **algorithm execution artifact**. Long-lived scientific
workflows should consume an analysis-facing object:

- `SharedFeatureSpace`: typed shared-coordinate identity with stable digest id
- `AlignedBlock`: one subject/session block (`observations × shared features`)
- `AlignedStudy`: one common shared space + blocks + model + lineage + safety
- `AlignedResampleSet`: retained fold models plus fold-specific `AlignedStudy`s
  (no raw stacking)

Orientation boundary (explicit, tested):

- Algorithm layer (`AlignmentData` / transforms): **features × observations**,
  left-multiply `Q %*% X`
- Analysis layer (`AlignedStudy`): **observations × shared features**,
  with one explicit transpose at the representation boundary

Core entrypoints:

- `as_aligned_study(result, ...)` — common-space results only
- `as_aligned_resample_set(result, ...)` — fold-specific spaces
- `align_study(model, data, mode = ...)` — frozen application / other modes
- `as_hyperdesign(aligned_study)` — adapter toward multidesign-style workflows

Invariant: an aligned numerical matrix is never separated from its observation
metadata, shared-space identity, model reference, and lineage.

`AlignedStudy` is a **generic scientific-data contract**, not an fMRI object.
It knows about blocks/subjects, observations, shared features, space identity,
provenance, CV role, and eager storage — not domain-specific acquisition or
file-format semantics.

## Non-goals

- No fMRI-specific feature extraction (GLMs, run/TR semantics, atlas loading).
- No pipeline DSL or multi-stage execution framework in core.
- No package-specific object application (e.g., “rotate NEF gamma”) in core.

Downstream packages should:

1) construct alignment-ready matrices (possibly via feature blocks), then
2) call `fit_alignment()` / `apply_alignment()`, then
3) either wrap outputs with `as_aligned_study()` / `align_study()` for analysis,
   or apply resulting operators to their own domain objects.
