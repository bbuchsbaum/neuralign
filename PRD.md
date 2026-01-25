# Neuralign Product Requirements Document

## North Star

> neuralign is the place where an alignment method becomes a reusable, auditable, composable, cross-validated transform family—independent of which package invented the workflow.

---

## The Problem: Reinventing the Wheel

The ecosystem has multiple packages implementing alignment (manifoldalign, topofmri, dkge, fmrireg.gnef, fmrireg.garrr). The duplication isn't in the core algorithms (Procrustes, OT, FUGW)—it's in **all the glue** that keeps getting rewritten:

- Reference/medoid selection
- "Fit on training, apply to held-out" workflows
- Operator application conventions (some use `X %*% Q`, others `Q %*% X`)
- Serialization of alignment artifacts
- Quality assessment
- fmrigds integration

Each package invents bespoke "fit objects" and bespoke "apply()" code.

---

## The Solution: One Contract, Many Methods

neuralign provides:

### 1. One Canonical Representation of "A Learned Alignment"

Every method returns the same shape:
- Per-subject operators in a declared `from_space → to_space`
- Traits: `orthogonal`, `mass_preserving`, `invertible`
- Provenance: method, params, package versions, hashes

### 2. One Operator Convention

- Operators are always `(target × source)` dimensions
- Application is always `Y = A %*% X` (left-multiply)
- OT couplings, Procrustes rotations, and linear maps all follow the same rule
- Eliminates "transpose/flip dims" bugs

### 3. Centralized Utilities Everyone Uses

| Utility | What It Does |
|---------|--------------|
| `select_reference()` | Medoid/centroid/template with consistent distance bookkeeping |
| `fit_alignment(cv=...)` | Reference selection happens inside fold (no leakage) |
| Leakage warnings | Automatic when fit/apply use overlapping subjects |
| `save_alignment()` / `load_alignment()` | Versioned metadata, reproducible |
| `alignment_quality()` | Standard summary dispatching to method-specific metrics |

### 4. Common Methods Become Shared Wrappers

- **Procrustes**: neuralign wraps manifoldalign, provides canonical orthogonal operator family
- **OT**: neuralign wraps manifoldalign GW/FPGW, standardizes coupling + variance propagation
- Other packages call these instead of rolling their own

### 5. No Dependency Cycles

```
neuralign (core + registry)
    ├── Imports: manifoldalign, Matrix, digest
    └── Suggests: topofmri, dkge, fmrireg.gnef, fmrigds

topofmri, dkge, fmrireg.gnef (method providers)
    └── Import neuralign
    └── Call register_aligner() in .onLoad()
```

### 6. Standard Downstream Integration

Output is convertible to `fmrigds::MapFamily`, so every pipeline consumes alignment the same way.

---

## What Neuralign Owns vs. Delegates

### Owns (Contract + Orchestration)
- Single alignment interface: fit, apply, compose, serialize, QA, CV hooks
- Standardized objects: `AlignmentData`, `AlignmentModel`, `AlignmentResult`
- Plugin/registry system: `register_aligner()` for ecosystem extension
- Reference selection utilities
- Operator composition (structural × functional)
- Cross-validation plumbing + leakage warnings

### Delegates (No Reimplementation)
- **Algorithms** → manifoldalign (Procrustes, KEMA, GW, FPGW)
- **Geometric transforms** → neurotransform/neurofunctor
- **Group statistics** → fmrigds/fmrireg

---

## Object Model

### AlignmentData
Input container for data to be aligned.
```r
AlignmentData(
  data,           # per-subject: list of matrices/NeuroVec
  subjects,       # character: subject IDs
  space,          # gds_space or compatible
  design = NULL,  # optional: task structure for dkge-like methods
  geometry = NULL # optional: adjacency for graph methods
)
```

### AlignmentModel
What you fit (separating fit from apply enables CV).
```r
AlignmentModel(
  transforms,     # per-subject operators (target × source)
  reference,      # medoid subject, template object, or "consensus"
  method,         # method name
  space_from,     # source space
  space_to,       # target space
  provenance      # params, package versions, hashes
)
```

### AlignmentResult
What you get after applying.
```r
AlignmentResult(
  model,          # AlignmentModel
  aligned,        # aligned data (lazy or realized)
  quality,        # diagnostics
  cv_info = NULL  # fold info if cross-validated
)
```

---

## Core API

### Naming (Avoids Conflict with fmrigds::align)
```r
# Primary verbs
fit_alignment(data, method = "procrustes", reference = "medoid", cv = "none", ...)
apply_alignment(model, new_data, ...)
compose_alignment(model1, model2)

# Registry
register_aligner(name, fit_fn, capabilities)
available_aligners()
aligner_capabilities(name)

# QA
alignment_quality(result, metrics = c("reconstruction", "correlation"))

# I/O
save_alignment(model, path)
load_alignment(path)

# fmrigds integration
as_map_family(model)
```

---

## Aligner Contract

### The `fit_fn` Signature

Every aligner must implement:
```r
fit_fn <- function(data, reference, train_idx, ...) {
  # MUST return:
  list(
    transforms = list(         # Named: subject_id → operator (target × source)
      "sub-01" = matrix(...),
      "sub-02" = matrix(...)
    ),
    reference_data = ...,      # The actual reference
    space_from = ...,
    space_to = ...,
    method_state = list()      # Method-specific state for apply
  )
}
```

### Capability Metadata

```r
register_aligner(
  name = "procrustes",
  fit_fn = procrustes_fit,
  capabilities = list(
    supports_cv = TRUE,
    needs_geometry = FALSE,
    needs_design = FALSE,
    returns_invertible = TRUE,
    transform_type = "orthogonal",  # orthogonal | linear | ot | permutation
    mass_preserving = TRUE,
    reference_types = c("subject", "consensus", "template")
  ),
  package = "manifoldalign",
  description = "Orthogonal Procrustes via manifoldalign"
)
```

**Note on embeddings:** neuralign is currently operator-based. Methods that
naturally produce embeddings should still expose an `(target × source)` operator
(often a projection) and store any extra embedding artifacts in `method_state`.

### dkge integration contract (operator-based)

`dkge` operates in *effect space* with a design-kernel metric `K` (a `q×q` PSD kernel over effects).
Its K-Procrustes alignment solves for an orthogonal `R` that maximizes `tr((U_ref^T K U) R)`, where
`U` is a `q×r` K-orthonormal basis (`U^T K U = I`).

To integrate with neuralign’s operator convention (`Y = A %*% X`), the dkge adapter should:

- **AlignmentData$data**: per-subject matrices `X_s` with shape `r×q` (components × effects), typically `X_s = t(U_s)`.
  Column order must match the effect ordering of `K`.
- **AlignmentData@design**: a list containing at least:
  - `K`: `q×q` design kernel in effect space.
  - `effects`: character vector length `q` matching row/col names of `K` and the columns of `X_s`.
- **Reference semantics**:
  - `reference="consensus"`: compute a consensus basis via `dkge_consensus_basis_K(U_list, K)` on training subjects.
  - `reference=<subject_id>`: use that subject’s basis as `U_ref`.
  - `reference=<matrix>`: accept a template basis (either `r×q` like `X_s`, or `q×r` like `U`).
  - Caveat: `reference="medoid"` in neuralign is Euclidean; if dkge needs a kernel-weighted medoid, pass an explicit subject id (computed in dkge) or use `consensus`.
- **Transforms returned to neuralign**:
  - For each subject, compute `R_s` using `dkge_procrustes_K(U_ref, U_s, K)`.
  - Return operator `A_s = t(R_s)` (an `r×r` orthogonal matrix) so that:
    - dkge’s right-multiply alignment: `U_s_aligned = U_s %*% R_s`
    - neuralign’s left-multiply convention on `X_s = t(U_s)`: `t(U_s_aligned) = t(R_s) %*% t(U_s) = A_s %*% X_s`
- **spaces**:
  - Recommend `space_from = space_to = "dkge_effect_basis"` (or a more specific string incorporating a kernel/effect hash).
- **Partial overlap (different subjects, different effect sets)**:
  - v1 recommendation: preprocess to a shared effect index/order in dkge before calling neuralign (e.g., via `dkge_align_effects()`), so `X_s` all share the same `q` and `effects`.
  - If the adapter supports partial overlap, it should own the `dkge_align_effects()` call and define how effect unions are built fold-safely (train-only union, optional priors), then fit bases on the completed kernels.

Recommended dkge `register_aligner()` capabilities:
`needs_design=TRUE`, `supports_cv=TRUE`, `transform_type="orthogonal"`, `returns_invertible=TRUE`,
`supports_new_subject=TRUE`, `reference_types=c("subject","consensus","template")`, `returns="operator"`.

---

## MVP Implementation Order

### Phase 1: Minimal Core
1. `AlignmentData` class + coercions
2. `AlignmentModel` class (thin wrapper around MapFamily + provenance)
3. `register_aligner()`, `available_aligners()`
4. `fit_alignment()` with method dispatch
5. `apply_alignment()` for held-out/new data
6. `save_alignment()` / `load_alignment()`

### Phase 2: First Method + Reference Selection
7. **Procrustes wrapper via manifoldalign** (one real method end-to-end)
8. `select_reference(data, method = "medoid")` - fold-safe
9. `as_map_family()` - first-class fmrigds pathway

### Phase 3: CV Infrastructure
10. CV hooks in `fit_alignment(cv = "loso" | "kfold")`
11. Leakage warnings

### Phase 4: Polish + Second Method
12. `alignment_quality()` metrics
13. OT wrapper via manifoldalign
14. Documentation + vignettes

### Phase 5: Ecosystem Migration
15. topofmri registers "fugw"
16. dkge registers "dkge"
17. fmrireg.gnef registers "nef"

---

## Package Structure

```
neuralign/
├── R/
│   ├── alignment_data.R      # AlignmentData class + coercions
│   ├── alignment_model.R     # AlignmentModel class
│   ├── alignment_result.R    # AlignmentResult class
│   ├── fit_alignment.R       # Main fit generic
│   ├── apply_alignment.R     # Apply to new data
│   ├── compose.R             # Operator composition
│   ├── registry.R            # register_aligner(), available_aligners()
│   ├── reference.R           # select_reference()
│   ├── quality.R             # alignment_quality()
│   ├── cv.R                  # Cross-validation utilities
│   ├── leakage.R             # Leakage detection/warnings
│   ├── serialize.R           # save/load
│   ├── compat_fmrigds.R      # as_map_family()
│   ├── aligner_procrustes.R  # Wrapper for manifoldalign
│   ├── aligner_ot.R          # Wrapper for manifoldalign OT
│   └── zzz.R                 # .onLoad registration
├── DESCRIPTION
├── NAMESPACE
└── tests/testthat/
```

---

## Practical Gotchas

**Operator semantics**: fmrigds uses `map %*% X` where map is `(n_target × n_source)`. Neuralign matches this exactly.

**topofmri exports**: Functions like `graph_from_neurovec()` and `align_to_template()` are NOT exported. Use high-level entry points (`align_epifmri()`, `fit_topo_template()`) or coordinate exports upstream.

**Naming**: fmrigds already exports `align()`. Use `fit_alignment()` / `apply_alignment()` to avoid masking.

---

## Verification Plan

1. **Unit tests**: Each class, each registry operation, each method wrapper
2. **Integration test**: `fit_alignment()` → `apply_alignment()` → `as_map_family()` → fmrigds pipeline
3. **CV test**: Verify train/test isolation, check leakage warning fires correctly
4. **Round-trip test**: `save_alignment()` → `load_alignment()` produces identical results

---

## Net Effect

You still have multiple algorithms and domain-specific workflows, but they all stop duplicating the same scaffolding. Procrustes/OT become "just another plugin that returns an operator family," and all the boring, failure-prone glue lives once in neuralign.
