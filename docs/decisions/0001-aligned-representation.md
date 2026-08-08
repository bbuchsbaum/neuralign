# Aligned representation contract

Status: accepted for the `0.1.x` implementation

## Context

`neuralign` provides common support for alignment algorithms. It owns the
matrix conventions, fit and apply orchestration, model and result objects,
cross-validation metadata, aligner registry, and serialization used by those
algorithms.

An `AlignmentResult` records one execution of an algorithm. Downstream code
also needs a stable description of the aligned values: which coordinate system
contains the values, which observations produced each row, and whether the
values came from training or held-out data. That description must work for GPA,
SRM, XPA-R, and future aligners without importing their method-specific
diagnostics.

## Decision

`neuralign` owns a method-neutral aligned-output contract:

- `SharedFeatureSpace` identifies one fitted or externally defined coordinate
  system.
- `AlignedBlock` stores one named block of observations in that coordinate
  system.
- `AlignedStudy` collects blocks that are proven to share one coordinate
  system.
- `AlignedResampleSet` retains fold-specific models and aligned assessment
  values when resampling produces more than one coordinate system.

The analysis-facing matrix convention is observations by shared features. The
algorithm-facing convention remains features by observations. Conversion
between them performs one explicit transpose.

## Required invariants

### Coordinate identity

`SharedFeatureSpace$id` identifies coordinates, not a method configuration.
Two fitted models with different coordinate-defining artifacts must not receive
the same ID merely because they used the same method, subjects, and dimension.

An aligner may provide an explicit coordinate identifier. Otherwise
`neuralign` derives a conservative fingerprint from the fitted model, including
its reference artifact, target-space description, method state, and fitted
operators. Conservative false negatives are acceptable: treating two equal
spaces as different prevents aggregation. False positives are not acceptable:
treating different spaces as equal corrupts an analysis.

Compatibility checks fail when either identity is missing. They never treat an
unknown space as compatible with a known space.

### Resampling

One `AlignedStudy` contains one coordinate system. Reassembled observations
from fold-specific spaces cannot be relabeled as a new common space.

Each `AlignedResampleSet` split retains the actual fold model, its
`SharedFeatureSpace`, and its aligned assessment values. It may also retain the
aligned analysis values. A result that discarded these fold artifacts cannot
be converted after the fact; the conversion reports what must be retained at
fit time.

Raw coordinates are not stacked across splits. Downstream code aggregates
predictions or metrics unless every split was explicitly mapped to one external
coordinate system.

### Analysis safety

Safety is an evidence state, not a caller-selected label. The supported states
are:

- `verified_safe`: recorded folds or observation identifiers establish the
  required separation;
- `declared`: the caller describes a safe use, but the object cannot prove it;
- `unknown`: the object lacks a sufficient declaration and evidence;
- `unsafe`: the recorded construction violates the requested use.

In-sample aligned values are `unsafe` for confirmatory cross-subject
assessment. `assert_analysis_safe()` rejects `declared`, `unknown`, and
`unsafe` values by default. A caller may explicitly allow a declared result,
but cannot promote it to verified status.

### Data and metadata

Every aligned block has a non-empty block name, a subject identifier, a numeric
finite matrix, one observation-data row per matrix row, and the study's shared
space ID. Subject/session metadata use block names as their primary key.

The initial implementation accepts in-memory base matrices and `Matrix`
objects. A storage list describes this eager representation; it is not a lazy
backend API.

## Package boundaries

`multidesign` owns `multidesign` and `hyperdesign`. The optional adapter in
`neuralign` calls its public constructors and returns genuine objects.

`multivarious` owns generic projector, score, reconstruction, and component
semantics. The aligned-output contract does not add a competing effect or
projector type.

`fmrilatent` may later provide an optional values backend. Core aligned objects
do not depend on it.

Aligner packages own estimators, diagnostics, simulations, and any method state
needed to fit new data. They may contribute a coordinate descriptor through
the aligner contract.

## Deferred work

The first release does not expose generic effect decoding, lazy storage,
cross-version schema migration, or an XPA-specific class. These require
separate provider or aligner contracts and conformance tests.
