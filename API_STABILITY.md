# API stability policy

This policy applies to the 0.2 release series. `neuralign` is still young, but
its ordinary workflow should not require users to chase unannounced renames.

## Stable core

These entry points define the ordinary workflow and are source-stable within
0.2.x:

- construction: `AlignmentData()`, `AlignmentModel()`, `AlignmentResult()`;
- fitting and application: `fit_alignment()`, `apply_alignment()`;
- inspection: `get_model()`, `get_aligned()`, `get_transform()`,
  `get_transforms()`, `get_quality()`, `get_cv_info()`;
- method discovery: `available_aligners()`, `aligner_capabilities()`;
- persistence: `save_alignment()`, `load_alignment()`; and
- analysis representation: `as_aligned_study()`, `AlignedStudy()`,
  `AlignedResampleSet()`, and shared-space compatibility checks.

Stable means argument and return semantics will not change incompatibly in a
patch release. It does not mean every optional backend supports new subjects,
new observations, inverse transforms, or decoding. Callers must inspect method
capabilities.

## Versioned extension contract

External method providers register through `register_aligner()` and declare an
`api_version`. Compatibility is governed by
`NEURALIGN_ALIGNER_API_VERSION`, independently of the package version.
neuralign supports one provider API at a time: a provider whose declaration
does not equal the current version fails during registration. The 0.2.0
release uses API 2; API-1 compatibility is not retained.

## Provisional surfaces

The following expert surfaces are public so current research pipelines can use
them, but their ergonomics may change in a minor 0.2 release with a documented
migration:

- guidance and intrinsic-geometry channels;
- feature-block harmonization and reporting;
- subspace restriction/lifting helpers;
- observation-axis cross-fit orchestration;
- `fmrigds` and hyperdesign conversions; and
- portable CSV/JSON/MATLAB import/export.

The serialized RDS envelope is not provisional: its format version and
migration behavior are part of the persistence contract.

## Deprecation and correction

- A convenience rename keeps a forwarding alias and emits a deprecation warning
  for at least one minor release.
- A behavior that returns the wrong mathematical object may be corrected or
  disabled immediately. The release notes must name the semantic break.
- Disabled methods remain discoverable through an actionable error rather than
  silently registering a known-bad estimator.
- Export count is reviewed before each minor release. New aliases are not added
  merely to mirror an upstream package's vocabulary.

For 0.2.0, KEMA is corrected from a linear-proxy claim to a nonlinear training
embedding contract, and CONE is disabled after failing the mandatory graph
oracle.
