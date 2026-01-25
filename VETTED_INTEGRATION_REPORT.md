# Vetting: External Aligner Integration Report (dkge / fmrireg.*)

Date: 2026-01-25

This document vets the provided “integration potential” report against the current `neuralign` repository state (R package `neuralign` v0.1.0) and highlights what is **verified in-code**, what is **unverified**, and what appears **misaligned with the current core contracts**.

## Executive summary

- **Verified in `neuralign` today:** a working aligner registry (`register_aligner()`, lazy auto-load mapping), reference selection (incl. “medoid”), subject-level cross-validation scaffolding (LOSO / k-fold), and leakage warnings.
- **Not present in `neuralign` today:** built-in adapters/wrappers for `dkge`, `fmrireg.gnef`, or `fmrireg.garrr`. If those methods are to be usable, the *external packages* must register themselves via `.onLoad()` (or `neuralign` must add wrappers).
- **Key contract mismatch to watch:** `apply_alignment()` assumes per-subject **matrix operators** and applies them via left-multiply (`transform %*% data`). “Embedding-only” methods are **not currently plug-compatible** unless they can provide operators or a custom `apply_fn`.
- **Local environment note:** in this session, `topofmri` is installed, but `dkge`, `fmrireg.gnef`, `fmrireg.garrr`, `fmrireg.lowrank`, and `manifoldalign` are **not installed**, so package-level API claims could not be validated here.

## What `neuralign` actually provides (relevant to the report)

### Registry + lazy loading

`neuralign` implements an internal registry plus an autoload hook that maps method names to packages:

- Registry implementation: `R/registry.R`
- Lazy-load mapping includes:
  - `"dkge" = "dkge"`
  - `"nef" = "fmrireg.gnef"`

Important implication: **having a mapping is not the same as being integrated**. Integration happens only if loading the provider package results in `register_aligner("<method>", ...)` being called (typically via `.onLoad()` in that provider package).

### Reference selection (incl. medoid)

`neuralign` has framework-level reference selection, including medoid selection:

- `select_reference()` supports `"medoid"`, `"centroid"`, `"first"`, `"random"`: `R/reference.R`

This medoid is computed in the **raw data space** (based on per-feature correlations or Euclidean/Frobenius distances), which may or may not be meaningful for methods that operate in a derived space (e.g., “effect-space”).

### CV scaffolding + leakage warnings

`fit_alignment()` supports subject-level cross-validation modes:

- `cv = "none" | "loso" | "kfold"`: `R/fit_alignment.R`
- Fold generation: `R/cv.R`
- Leakage checks and workflow assessment: `R/leakage.R`

This is general plumbing. Whether a specific method is *actually safe* under CV depends on its `capabilities` and implementation.

### Apply semantics (critical for “embedding” claims)

Core apply is operator-based:

- `apply_alignment()` applies existing transforms as `transform %*% subj_data`: `R/apply_alignment.R`
- Default “apply to new subjects” path re-calls `fit_fn()` and expects `fit_result$transforms[[subject]]` to be a matrix-like operator: `R/apply_alignment.R`

While the registry capability schema includes `returns = "operator" | "embedding"` (`R/registry.R`), **the current apply path does not branch on `returns`**, and therefore does not implement an “embedding-only” workflow.

## Vetting the provided package assessments

The table below focuses on what can be confirmed from the `neuralign` repo itself. Claims about *external packages’ internal algorithms* are marked unverified unless their code/docs are present in this repo or installed locally.

### 1) `dkge` — “HIGH Integration Potential”

| Claim from report | Vetting verdict | Notes / evidence |
|---|---|---|
| “Plugin registry ready — neuralign already has ‘dkge’ in lazy-load map” | **Verified** | `R/registry.R` maps `"dkge" = "dkge"`; `DESCRIPTION` lists `dkge` under `Suggests`. |
| “Multi-subject data handling / effect alignment with partial overlap support” | **Unverified** | `dkge` package not present in this repo and not installed in this environment. |
| “CV workflows — LOSO/K-fold contrast engine” | **Unverified** | `neuralign` has LOSO/k-fold *alignment* CV plumbing (`R/fit_alignment.R`, `R/cv.R`), but no “contrast engine” appears in core; `dkge` internals not verified. |
| “Reference selection — Needs adaptation (no built-in medoid)” | **Partially incorrect / needs nuance** | `neuralign` provides framework-level medoid selection (`R/reference.R`). However, if `dkge` needs medoid selection in *effect/kernel space*, then adaptation may indeed be required. |
| “Mark as `needs_design = TRUE`, `returns = \"embedding\"`” | **Mixed** | `needs_design` is supported as a capability and enforced (`R/registry.R`). However, `returns = "embedding"` is **not currently supported** by `apply_alignment()`; embedding-only methods would require either operators, a custom `apply_fn`, or core changes. |

**Bottom line:** `dkge` can be “high potential” **only if** it can supply a per-subject operator consistent with `neuralign`’s apply semantics (or it ships a compatible `apply_fn`). If it is fundamentally “embedding-returning”, integration requires extending core apply/results contracts.

### 2) `fmrireg.garrr` — “MODERATE (Complementary Role)”

| Claim from report | Vetting verdict | Notes / evidence |
|---|---|---|
| “Produces transforms — projectors (voxel→anchor), not subject-to-subject alignment” | **Unverified** | Package not present/installed here. |
| “neurotransform integration — excellent MNI/template warping” | **Unverified** | Not verifiable from this repo. |
| “Recommended role: downstream analysis engine, not a registered aligner” | **Plausible** | This aligns with `neuralign`’s current scope (alignment/orchestration) vs. downstream stats, but it’s a product decision, not something `neuralign` enforces today. |

**Bottom line:** nothing in `neuralign` currently auto-loads or registers a `garrr` method. Treating it as downstream tooling is reasonable, but it’s outside the current aligner contract either way.

### 3) `fmrireg.gnef` — “HIGH Integration Potential”

| Claim from report | Vetting verdict | Notes / evidence |
|---|---|---|
| “Procrustes alignment / medoid selection / cross-fit alignment exist” | **Unverified** | Package not present/installed here. |
| “Transform output — K×K orthogonal Q matrices per subject” | **Unverified** | If true, this would be a strong fit for `neuralign` because operators can be applied by left-multiply. |
| “Just needs adapter wrapper” | **Partially true** | `neuralign` has a lazy-load mapping for `"nef" = "fmrireg.gnef"` (`R/registry.R`), but `neuralign` itself does not include a wrapper. Integration still requires `fmrireg.gnef` to call `register_aligner("nef", ...)` on load (or `neuralign` to add a wrapper). |

**Additional discrepancy:** `PRD.md` discusses `fmrireg.gnef` as a suggested provider package, but the current `DESCRIPTION` does **not** list it under `Suggests`. If `nef` integration is intended, that packaging metadata likely needs updating.

### 4) `fmrireg.lowrank` — “NOT VIABLE”

| Claim from report | Vetting verdict | Notes / evidence |
|---|---|---|
| “Single-subject only; not an aligner” | **Unverified** | Package not present/installed here. |
| “Recommended role: upstream first-level engine feeding coefficients into neuralign” | **Plausible** | Consistent with `neuralign` being alignment/orchestration, but not validated against package docs here. |

## Revised integration priority (based on `neuralign`’s current contracts)

1. **`fmrireg.gnef` (as `nef`)** — Highest *if* it truly provides per-subject orthogonal operators; minimal mismatch with current apply semantics.
2. **`dkge`** — High only if it can be expressed as operators (or ships a custom `apply_fn`). If it is embedding-only, core changes are required and effort rises substantially.
3. **`fmrireg.garrr`** — Treat as downstream analysis; not an aligner unless a clear operator/apply story exists.
4. **`fmrireg.lowrank`** — Likely not an aligner; upstream preprocessing/inference tooling.

## Action items / open questions for a follow-up iteration

- Decide whether `neuralign` should officially support `capabilities$returns = "embedding"` (this likely requires extending `AlignmentModel`, `AlignmentResult`, and `apply_alignment()` semantics).
- Confirm that each provider package (`dkge`, `fmrireg.gnef`, etc.) actually registers its method in `.onLoad()` via `neuralign::register_aligner()`.
- If `nef` is intended as a supported provider, add `fmrireg.gnef` to `DESCRIPTION` `Suggests` (and potentially to any CI install matrix).
- Define and document a **standard design format** for `AlignmentData@design` if `dkge`-like supervised methods are first-class.

