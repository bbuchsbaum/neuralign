# neuralign product requirements

## Product statement

`neuralign` is the contract and orchestration layer for fitting, applying,
evaluating, and retaining cross-subject alignment artifacts. It lets analysis
code use one set of data, model, result, cross-validation, safety, and
serialization conventions while algorithm packages continue to own their
estimators.

Version 0.1 targets research developers who need explicit scientific semantics
more than a single all-purpose alignment algorithm. It is not a preprocessing
pipeline, a feature extractor, or a claim that every registered backend has the
same lifecycle.

## The three contracts users should remember

1. Input matrices use `features x observations` orientation.
2. Operator-returning methods use target-by-source transforms and apply them as
   `aligned = transform %*% input`.
3. A shared coordinate system is an identified artifact. Data may be combined
   only when its `SharedFeatureSpace` identities agree; fold-specific spaces are
   retained in an `AlignedResampleSet` and are not stacked as raw coordinates.

## Public artifact model

### Algorithm artifacts

- `AlignmentData` owns per-subject inputs and optional observation labels,
  design, geometry, and guidance.
- `AlignmentModel` owns fitted transforms, reference state, method state,
  provenance, and training membership.
- `AlignmentResult` owns a model plus aligned outputs, diagnostics, and
  cross-validation evidence.

An aligner declares one of two truthful return contracts:

- `returns = "operator"`: the model stores reusable target-by-source operators.
- `returns = "embedding"`: the model stores fitted embedding artifacts. Such a
  method must set `supports_new_data = FALSE` unless it implements a genuine
  out-of-sample extension.

### Analysis artifacts

- `AlignedStudy` stores analysis-facing `observations x shared features`
  blocks, one shared-space identity, lineage, and safety declarations.
- `AlignedResampleSet` stores the exact model and coordinate system for each
  retained split. Only metrics and predictions may be aggregated across
  incompatible split coordinates.

The conversion boundary between algorithm and analysis orientations is one
explicit transpose.

## Backend ownership and lifecycle

`neuralign` owns adapter semantics, capability declarations, shape validation,
fit/apply routing, and mandatory conformance fixtures. Backend packages own the
numerical estimator.

| Backend path | Lifecycle in 0.1 | Contract |
|---|---|---|
| `procrustes`, `procrustes_graph`, `kprocrustes` | built in | reusable operators |
| `gw`, `fpgw` | optional `manifoldalign` | transport operators |
| `gpca`, `coupled_diag`, `lowrank` | optional `manifoldalign` | latent operators; optional reference-space lift |
| `grasp` | optional `manifoldalign` | graph assignment operator; mandatory permutation oracle |
| `kema` | optional `manifoldalign` | nonlinear training embedding; no new-subject or new-data application |
| `cone`, `cone_align` | disabled | fails the frozen graph-correspondence oracle |
| external providers such as `dkge` | plugin | provider registers against the versioned aligner API |

Registration is not evidence of correctness. Every supported optional backend
must have a deterministic contract fixture that runs in the backend CI lane.
Methods that fail their contract are fixed or disabled; their tests are not
made opt-in.

## Safety requirements

- Reference selection occurs within the analysis fold for cross-validation.
- In-sample, frozen-application, and cross-fitted outputs carry distinct safety
  records.
- Confirmatory cross-subject analysis fails closed unless the artifact contains
  verified cross-fit evidence or the caller explicitly accepts a declared mode.
- Fold-specific coordinates are not represented as one global reference.
- Unknown shared-space identities fail compatibility checks.

## Persistence requirements

`save_alignment()` writes a versioned envelope for `AlignmentModel`,
`AlignmentResult`, `AlignedStudy`, and `AlignedResampleSet`.

- The envelope records its format version, object kind, package/R versions, and
  a SHA-256 object hash.
- Integrity mismatches and newer unsupported formats fail closed.
- The version-1 model/result envelope and raw legacy RDS objects remain readable
  through explicit compatibility paths.
- A future object-layout change must add a tested migration before incrementing
  the serialization format.

## Release requirements

A release tag may be created only when all locally executable gates pass and
the external evidence status is written down rather than inferred.

Required gates:

1. full tests and `R CMD check` on the exact source tree;
2. a clean source build/install smoke test;
3. a minimal-Suggests check;
4. Linux, macOS, and Windows CI checks;
5. mandatory manifoldalign semantic and numerical conformance;
6. scheduled downstream provider compatibility;
7. generated documentation matches roxygen source; and
8. `NEWS.md`, release evidence, tracker state, Git commit, tag, and remote refs
   identify the same release.

Local passage is not cross-platform evidence. Workflow configuration is not a
served CI receipt. See `RELEASE.md` for the live gate ledger.

## API policy

The core workflow is stabilized for the 0.1 series. Advanced representation,
guidance, subspace, and integration surfaces remain provisional and must be
identified as such. The aligner plugin contract is versioned independently by
`NEURALIGN_ALIGNER_API_VERSION`.

No exported function is removed silently. Renames require an alias and warning
for at least one 0.1 minor release unless the old behavior is scientifically
incorrect or unsafe. The current tiers and migration rules are in
`API_STABILITY.md`.

## Current roadmap

The remaining planned representation extensions are deliberately separate:

- a sparse/lazy aligned-value provider contract; and
- aligner-owned decoder capabilities.

Neither is implied by eager `AlignedStudy` storage or by a reference-space
low-rank operator. Each needs its own ownership, serialization, and validation
contract before implementation.

## Non-goals for 0.1

- hiding input orientation or reference semantics;
- fabricating linear operators for nonlinear embeddings;
- decoding without an aligner-owned decoder contract;
- treating raw coordinates from different CV folds as interchangeable;
- importing every optional estimator into the core package; or
- claiming CRAN, platform, performance, or downstream compatibility without a
  corresponding receipt.
