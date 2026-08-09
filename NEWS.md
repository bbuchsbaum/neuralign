# neuralign 0.1.0

## Core contracts

- Added aligner API v2 with a one-shot provider preflight hook, detached exact
  subject/observation fit contexts, and explicit API-v1 callback compatibility.
- Made observation-axis train/test overlap an engine-level error before any
  provider callback, and materialized generated subject folds once for both
  preflight and execution.
- Added a unified fit/apply workflow with target-by-source operator semantics,
  fold-aware reference handling, capability validation, and leakage checks.
- Added analysis-facing `AlignedStudy`, `SharedFeatureSpace`, and
  `AlignedResampleSet` artifacts with explicit coordinate identity and safety
  metadata.
- Added version-2 serialization for models, results, studies, and retained
  resample sets, with SHA-256 verification and legacy version-1 loading.

## Backend correctness

- Corrected KEMA to return manifoldalign's nonlinear training embeddings rather
  than presenting primal-vector diagnostics as reusable linear operators. KEMA
  now rejects new-subject, new-data, reference-space, and CV application until
  a truthful out-of-sample kernel extension exists.
- Disabled CONE and its alias after the upstream estimator failed the frozen
  graph-correspondence oracle. GRASP remains the supported graph-assignment
  path.
- Made manifoldalign semantic, dimensional, and numerical conformance tests a
  mandatory CI lane instead of an environment-variable opt-in.

## Release engineering

- Added Linux, macOS, and Windows package checks, clean-install and
  minimal-Suggests lanes, mandatory backend conformance, and a scheduled `dkge`
  provider smoke test.
- Added an API stability policy, release evidence ledger, deterministic backend
  reproductions, and current product requirements.
