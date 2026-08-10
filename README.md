# neuralign

`neuralign` is an R package for fitting, applying, composing, and evaluating
cross-subject alignment transforms through one consistent interface.

It acts as the orchestration layer around alignment methods: the package
standardizes data containers, model/result objects, reference selection,
cross-validation, leakage checks, feature harmonization, and serialization,
while allowing the actual alignment algorithm to come from built-in methods or
companion packages such as `manifoldalign`.

## Why this package exists

Alignment workflows tend to repeat the same infrastructure work:

- choosing a reference subject or template
- fitting on training data and applying to held-out data
- keeping transform conventions straight
- computing quality summaries
- storing provenance and serialized artifacts

`neuralign` centralizes that shared machinery so downstream packages can focus
on domain-specific features and alignment backends instead of rewriting the same
plumbing.

## What `neuralign` standardizes

- `AlignmentData`: a common input container for multi-subject matrices and
  optional metadata such as observation labels, geometry, or design structure
- `AlignmentModel`: fitted transforms, reference information, method metadata,
  and provenance
- `AlignmentResult`: aligned outputs, quality summaries, and CV metadata
- `AlignedStudy`: analysis-facing observations in one identified shared space
- `AlignedResampleSet`: retained fold artifacts when each split has a different
  shared coordinate system
- `fit_alignment()` / `apply_alignment()`: one fit/apply workflow across methods
- `compose_alignment()`: linear composition of alignment models
- `register_aligner()`: a registry for extending the ecosystem with new methods
- subject-axis CV (`loso`, `kfold`) and observation-axis crossfit utilities
- leakage warnings and audit helpers
- feature-block utilities for portable correspondence signals

## Core conventions

Two conventions matter throughout the package:

- subject data are stored as `features x observations`
- transforms are left-multiply operators with shape `target x source`, so
  aligned data are computed as `Q %*% X`

These rules are enforced consistently across fitting, applying, composing, and
quality assessment.

## Installation

`neuralign` is not currently on CRAN.

Install from GitHub with either `pak` or `remotes`:

```r
pak::pak("bbuchsbaum/neuralign")
```

```r
remotes::install_github("bbuchsbaum/neuralign")
```

Some alignment methods depend on optional suggested packages. The package works
with its built-in methods out of the box, and registers additional backends
when supporting packages are installed.

## Quick start

```r
library(neuralign)

set.seed(42)
n_features <- 20
n_obs <- 30

# Shared signal seen by every subject
shared <- matrix(rnorm(n_features * n_obs), n_features, n_obs)

data_list <- lapply(1:4, function(i) {
  Q <- qr.Q(qr(matrix(rnorm(n_features^2), n_features)))
  Q %*% shared + 0.3 * matrix(rnorm(n_features * n_obs), n_features, n_obs)
})
names(data_list) <- paste0("sub-", sprintf("%02d", 1:4))

# Wrap input matrices
adat <- AlignmentData(data_list)

# Fit a built-in aligner
fit <- fit_alignment(
  adat,
  method = "procrustes",
  reference = "consensus"
)

fit
get_quality(fit)

# Apply the fitted model to a new subject
new_subject <- list(
  "sub-05" = {
    Q <- qr.Q(qr(matrix(rnorm(n_features^2), n_features)))
    Q %*% shared + 0.3 * matrix(rnorm(n_features * n_obs), n_features, n_obs)
  }
)

applied <- apply_alignment(fit, AlignmentData(new_subject))
names(get_aligned(applied))
```

## Choosing a method

Methods are discovered through the aligner registry. In a plain installation,
the built-in methods are:

- `procrustes`
- `procrustes_graph`
- `kprocrustes`

When `manifoldalign` is installed, `neuralign` also attempts to register
additional methods including:

- `gw`
- `fpgw`
- `kema`
- `coupled_diag`
- `gpca`
- `grasp`
- `lowrank`

`cone` is intentionally disabled: the upstream estimator does not currently
meet neuralign's deterministic graph-correspondence accuracy contract. Use
`grasp` for the supported graph-alignment path.

Use `available_aligners()` or `available_aligners(details = TRUE)` to inspect
what is registered in your current session.

Capabilities are method-specific. In particular, `kema` returns nonlinear
training embeddings and deliberately does not support new-subject or new-data
application. Registration should never be read as a promise that every method
supports the same fit/apply lifecycle.

## From a fit to an analysis artifact

`AlignmentResult` records one algorithm run. Convert a non-fold-specific result
to an `AlignedStudy` when downstream code needs analysis-facing
`observations x shared features` matrices plus coordinate identity, lineage,
and safety metadata:

```r
study <- as_aligned_study(fit, source_data = adat)
aligned_matrix(study, "sub-01")
shared_space(study)$id
```

For cross-validation with fold-specific coordinates, refit with
`return_resample_artifacts = TRUE` and use `as_aligned_resample_set()`. Raw
coordinates from different splits cannot be stacked; aggregate metrics or
predictions instead.

## Saving artifacts

`save_alignment()` and `load_alignment()` support models, results,
`AlignedStudy` objects, and `AlignedResampleSet` objects. The RDS envelope is
versioned and integrity checked. Saving an `AlignmentResult` stores its model by
default; set `include_data = TRUE` to retain aligned matrices.

```r
path <- tempfile(fileext = ".rds")
save_alignment(study, path)
restored <- load_alignment(path)
unlink(path)
```

## Cross-validation and leakage control

`neuralign` treats evaluation as a first-class concern.

- subject-axis CV is available through `fit_alignment(cv = "loso")` and
  `fit_alignment(cv = "kfold")`
- observation-axis crossfit is available through `create_obs_folds()` and
  `run_obs_crossfit_from_data()`
- reference selection is performed within training folds to reduce leakage
- `assess_leakage_risk()` provides a more explicit audit when needed

This makes it easier to compare methods without silently evaluating on training
information.

## Feature blocks and guidance channels

Many workflows do not start from one clean shared matrix per subject. To support
that, `neuralign` includes utilities for building and harmonizing heterogeneous
correspondence signals:

- `alignment_feature_block()`
- `harmonize_feature_blocks()`
- `stack_feature_blocks()`
- `build_alignment_features()`
- `block_alignment_report()`

Optional guidance channels such as ROI projectors and intrinsic geometry are
also supported for methods that can use them as priors.

## Extending the ecosystem

Packages can register new alignment methods and immediately inherit the core
workflow around fitting, applying, CV, composition, quality summaries, and
serialization.

The aligner contract is intentionally small: a method provides a `fit_fn`,
optionally an `apply_fn`, and a capability description via
`register_aligner()`. Aligner API v2 is the only supported provider contract.
Every fit callback receives the exact fold or full-data context; an optional
preflight callback can validate the complete resampling plan once before
fitting begins. See the extension vignette for the expected interfaces.

## Documentation map

The package ships several vignettes covering common workflows:

- `vignette("neuralign")`: end-to-end getting started guide
- `vignette("introduction")`: package overview and method-selection framing
- `vignette("cv_leakage")`: subject-axis CV, observation-axis crossfit, and
  leakage prevention
- `vignette("advanced_blocks")`: feature blocks and subspace workflows
- `vignette("anchors")`: anchors, guidance, and correspondence regimes
- `vignette("extending")`: writing and registering custom aligners
- `vignette("fmri_profile")`: fMRI-oriented usage patterns

Additional project context:

- [ARCHITECTURE.md](ARCHITECTURE.md)
- [PRD.md](PRD.md)
- [API_STABILITY.md](API_STABILITY.md)
- [RELEASE.md](RELEASE.md)

## Scope

`neuralign` is intentionally not an fMRI preprocessing package and not a
feature-extraction framework. It expects downstream code to construct
alignment-ready matrices, then uses the package to fit and apply the resulting
transform family in a consistent, auditable way.
