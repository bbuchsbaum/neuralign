# Neuralign Synthesis Report: Cross-Package Analysis and Architecture Design

## Executive Summary

This report synthesizes findings from 9 R packages in the fMRI analysis ecosystem to inform the design of `neuralign`, a unified framework for functional neural alignment. The analysis reveals significant duplication of alignment approaches across packages, inconsistent abstractions, and clear opportunities for a unifying layer.

**Key Finding**: The ecosystem has mature implementations of alignment algorithms (manifoldalign), geometric transforms (neurotransform/neurofunctor), and group analysis (fmrigds/fmrireg), but lacks a coherent abstraction layer that bridges these components for functional alignment workflows.

---

## 1. Package Ecosystem Overview

### 1.1 Package Roles and Dependencies

```
                    ┌─────────────────────────────────────────────┐
                    │              USER APPLICATIONS              │
                    └─────────────────────────────────────────────┘
                                         │
          ┌──────────────────────────────┼──────────────────────────────┐
          │                              │                              │
          ▼                              ▼                              ▼
   ┌─────────────┐              ┌─────────────────┐            ┌─────────────┐
   │   fmrireg   │              │    fmrigds      │            │    dkge     │
   │  First-Level│              │  Group Data     │            │Design-Kernel│
   │     GLM     │              │   Structures    │            │  Embedding  │
   └──────┬──────┘              └────────┬────────┘            └──────┬──────┘
          │                              │                            │
          ├──────────────────────────────┼────────────────────────────┤
          │                              │                            │
   ┌──────┴──────┐              ┌────────┴────────┐            ┌──────┴──────┐
   │fmrireg.garrr│              │  fmrireg.gnef   │            │  topofmri   │
   │ Anchor-Space│              │Neural Factorize │            │ Topological │
   │   Low-Rank  │              │  Latent Align   │            │  FUGW Align │
   └─────────────┘              └─────────────────┘            └─────────────┘
          │                              │                            │
          └──────────────────────────────┼────────────────────────────┘
                                         │
                              ┌──────────┴──────────┐
                              │    manifoldalign    │
                              │  Alignment Library  │
                              │ KEMA, GW, Procrustes│
                              └──────────┬──────────┘
                                         │
          ┌──────────────────────────────┼──────────────────────────────┐
          │                              │                              │
          ▼                              ▼                              ▼
   ┌─────────────┐              ┌─────────────────┐            ┌─────────────┐
   │neurofunctor │              │  neurotransform │            │   neuroim2  │
   │   Domain    │              │    Geometric    │            │   Volume    │
   │   Graphs    │              │   Transforms    │            │   I/O       │
   └─────────────┘              └─────────────────┘            └─────────────┘
```

### 1.2 Package Summary Table

| Package | Lines | Primary Role | Alignment Approach | Group Support |
|---------|-------|--------------|-------------------|---------------|
| **fmrireg** | ~15K | First-level GLM | Assumes pre-aligned (MNI) | group_data + meta-analysis |
| **fmrireg.garrr** | ~12K | Low-rank engines | Template projectors | Anchor-space group GLM |
| **fmrireg.gnef** | ~6.5K | Neural factorization | Procrustes in latent space | Medoid-based group stats |
| **fmrigds** | ~10K | Group data structures | MapFamily abstraction | Meta-analysis reducers |
| **topofmri** | ~8K | Topological alignment | FUGW optimal transport | Template training |
| **dkge** | ~9K | Design-kernel embedding | OT + anchor alignment | LOSO cross-validation |
| **manifoldalign** | ~12K | Alignment algorithms | KEMA, GW, Procrustes, etc. | Multi-domain sync |
| **neurofunctor** | ~7K | Spatial composition | Projector paths | Per-subject graphs |
| **neurotransform** | ~4K | Transform primitives | Morphism composition | Domain hashing |

---

## 2. Alignment Paradigms Identified

### 2.1 Geometric/Anatomical Alignment
**Packages**: neurotransform, neurofunctor

**Approach**: Transform coordinates between anatomical spaces (native → MNI → template)

**Key Abstractions**:
- `Morphism`: Typed transforms (Affine3D, Warp3D, VolToSurf, SurfToSurf)
- `MorphismPath`: Composable sequences with pullback semantics
- `Projector`: Sparse matrix operators compiled from paths
- `BrainGraph`: Domain nodes + morphism edges

**Strengths**:
- Single interpolation (eliminates cumulative blur)
- Exact/approximate/adjoint invertibility taxonomy
- Provenance tracking
- Tool interoperability (ANTs, FSL, AFNI, FreeSurfer)

**Limitations**:
- Assumes anatomical correspondence is sufficient
- No functional adaptation

### 2.2 Latent-Space Alignment
**Packages**: fmrireg.gnef, fmrireg.garrr, dkge

**Approach**: Project to low-dimensional space, align there, project back

**Key Abstractions**:
- Anchor-space projectors P (V × M) mapping voxels ↔ anchors
- Neural drivers S (T × K) in component space
- Orthogonal Procrustes Q for latent alignment

**Strengths**:
- Dramatic dimensionality reduction (V >> K)
- Denoising through projection
- Fast group statistics in low-D space

**Limitations**:
- Loss of fine spatial detail
- Alignment quality depends on K selection
- Different packages use incompatible projectors

### 2.3 Optimal Transport Alignment
**Packages**: manifoldalign, topofmri, dkge

**Approach**: Find optimal correspondence via transport plans

**Key Methods**:
- **Sinkhorn**: Entropy-regularized transport (fast)
- **Gromov-Wasserstein**: Structure-only transport (features incomparable)
- **Fused GW (FUGW)**: Combines features + structure with α parameter
- **Partial transport**: Handles mass mismatch

**Strengths**:
- Principled handling of non-bijective correspondences
- Flexible feature/structure trade-off
- Mass-preserving properties for valid statistics

**Limitations**:
- Computational cost O(n²) or worse
- Requires careful regularization selection

### 2.4 Manifold/Kernel Alignment
**Packages**: manifoldalign

**Approach**: Preserve manifold geometry while aligning class structure

**Key Methods**:
- **KEMA**: Kernel manifold alignment with semi-supervised support
- **Generalized Procrustes**: Orthogonal alignment for partial observations
- **GRASP/CONE-Align**: Graph spectral alignment
- **Coupled Diagonalization**: Multi-modal spectral alignment

**Strengths**:
- Preserves intrinsic data geometry
- Handles partial labels (semi-supervised)
- Multiple paradigms for different data types

**Limitations**:
- Requires kernel/similarity specification
- u parameter (geometry/class trade-off) needs tuning

### 2.5 Topological Alignment
**Packages**: topofmri

**Approach**: Align via functional connectivity graph topology

**Key Concepts**:
- Laplacian eigenbases as functional coordinate system
- Weak orientation anchors (LR/AP/SI) for disambiguation
- FUGW barycenter for template construction

**Strengths**:
- No structural image needed at inference
- Directly targets functional organization
- Handles native-space fMRI

**Limitations**:
- Novel approach, less validated
- Requires MNI data for template training

### 2.6 Design-Aware Alignment
**Packages**: dkge

**Approach**: Incorporate experimental design structure into alignment

**Key Concepts**:
- Design kernels encoding factorial structure
- Effect-space alignment respecting task similarity
- LOSO cross-validation for unbiased inference

**Strengths**:
- Task-aware grouping of conditions
- Principled handling of partial task overlap
- Integrated classification framework

**Limitations**:
- Requires design specification
- Less applicable to resting-state

---

## 3. Common Data Structures

### 3.1 Projector/Operator Abstractions

| Package | Class | Structure | Operations |
|---------|-------|-----------|------------|
| neurofunctor | `Projector` | dgCMatrix (rows=target, cols=source) | apply, transpose, compose |
| fmrireg.garrr | Numeric matrix P | V × M (voxels × anchors) | forward, pseudoinverse |
| fmrigds | `map_linear` | operator matrix or function | apply, variance propagation |
| topofmri | OT coupling P | N_s × N_t | soft alignment |

**Opportunity**: Unified `Projector` interface with:
- `forward(x)`: Apply transform
- `backward(x)`: Apply inverse/adjoint
- `variance_propagate(var)`: Uncertainty handling
- `compose(other)`: Chain projectors

### 3.2 Space/Domain Abstractions

| Package | Class | Contents |
|---------|-------|----------|
| neurofunctor | `Domain` = Space + Geometry | Hash, coordinates, masks |
| neurotransform | Hash strings (opaque) | Source/target identifiers |
| fmrigds | `gds_space` hierarchy | space_voxel, space_parcels, space_surface, space_basis |
| topofmri | `TopoGraph` | Adjacency, features, indices |

**Opportunity**: Unified `Space` with:
- Type discriminant (voxel/surface/parcel/basis/graph)
- Coordinate accessor
- Mask support
- Hash for identity

### 3.3 Transform/Map Abstractions

| Package | Class | Subject-Specific? | Invertibility |
|---------|-------|------------------|---------------|
| neurotransform | `Morphism` hierarchy | Via domain hash | exact/approximate/adjoint/none |
| neurofunctor | `MorphismPath` | Per-subject BrainGraph | Quality flags |
| fmrigds | `MapFamily` | `by_subject` list | UncertaintyRule |
| manifoldalign | `align_transform` | O/GL/perm types | compose/invert methods |

**Opportunity**: Unified `AlignmentMap` with:
- Subject-specific operators
- Invertibility metadata
- Variance propagation rules
- Provenance

### 3.4 Group Data Containers

| Package | Class | Dimensions | Formats |
|---------|-------|------------|---------|
| fmrireg | `group_data` | Varies by format | CSV, NIFTI, H5, GDS |
| fmrigds | `gds` | [sample × subject × contrast] | Unified lazy |
| dkge | `dkge_data` | betas + designs | Lists |
| manifoldalign | `hyperdesign` | x + design per domain | Lists |

**Opportunity**: `fmrigds` already provides the most sophisticated abstraction; neuralign should integrate with it via `MapFamily`.

---

## 4. Duplicated Functionality

### 4.1 Procrustes Alignment

**Implementations**:
- `manifoldalign::generalized_procrustes()` - General, partial observations
- `fmrireg.gnef::nef_procrustes_align()` - Latent timecourses
- `topofmri` (basis alignment) - Laplacian eigenbases
- `dkge` (internal) - Transport alignment

**Recommendation**: Use manifoldalign as canonical implementation; provide adapters for package-specific data structures.

### 4.2 Optimal Transport

**Implementations**:
- `manifoldalign::gromov_wasserstein()`, `fpgw()`, `parrot()` - General algorithms
- `topofmri::rcpp_fugw_*()` - Fused GW for graphs
- `dkge::dkge_transport_*()` - Sinkhorn for cluster alignment

**Recommendation**: Consolidate into manifoldalign; topofmri/dkge use as dependency or re-export.

### 4.3 Group Statistics

**Implementations**:
- `fmrireg::fmri_meta()`, `fmri_ttest()` - Traditional meta-analysis
- `fmrigds` reducers - Lazy pipeline meta-analysis
- `fmrireg.gnef::nef_group_stats()` - Latent-space group contrasts
- `dkge::dkge_group_*()` - Design-kernel group inference

**Recommendation**: fmrigds provides the most flexible infrastructure; other packages should register custom reducers.

### 4.4 Template/Anchor Construction

**Implementations**:
- `fmrireg.garrr::build_anchors_fps()` - Farthest-point sampling
- `topofmri::fit_topo_template()` - FUGW barycenter
- `dkge::dkge_make_anchors()` - k-means anchors

**Recommendation**: Consolidate anchor strategies into neuralign with method dispatch.

---

## 5. Identified Gaps

### 5.1 No Unified Alignment Interface

Each package defines its own:
- Input format (matrices, NeuroVec, hyperdesign, group_data)
- Output format (aligned matrices, transport plans, projectors)
- Method invocation (different function signatures)

**Need**: Common `align(data, method, ...)` interface returning standardized results.

### 5.2 No Hybrid Alignment Pipelines

Current packages operate in isolation:
- neurofunctor handles geometric alignment
- manifoldalign handles functional alignment
- No easy way to compose: native → MNI → hyperalign → template

**Need**: Pipeline composition that chains geometric + functional alignment.

### 5.3 No Cross-Validation Aware Alignment

Only fmrireg.gnef (Mode B) and dkge (LOSO) address leakage:
- Most alignment methods don't prevent double-dipping
- Group statistics may be inflated if alignment uses test data

**Need**: Built-in cross-validation support with leakage warnings.

### 5.4 No Alignment Quality Metrics

Scattered diagnostic functions:
- `manifoldalign::alignment_quality()` - Label-free metrics
- `topofmri::qa_*()` - Template-specific QA
- No standard quality reporting

**Need**: Unified alignment quality assessment framework.

### 5.5 No Alignment Persistence

Each package has ad-hoc serialization:
- neurofunctor: projector cache
- topofmri: save_topo_template()
- dkge: transport specs

**Need**: Standard alignment artifact format for sharing/reuse.

---

## 6. Proposed Neuralign Architecture

### 6.1 Design Principles

1. **Orchestration, not reimplementation**: Wrap existing packages
2. **Unified interfaces**: Common API across alignment methods
3. **Composable pipelines**: Chain geometric + functional alignment
4. **Statistical safety**: Built-in cross-validation, leakage prevention
5. **Provenance tracking**: Full audit trail
6. **Format agnostic**: Work with any input via adapters

### 6.2 Core Abstractions

```r
# 1. Alignment Method Registry
register_aligner("procrustes", manifoldalign::generalized_procrustes, ...)
register_aligner("kema", manifoldalign::kema, ...)
register_aligner("fugw", topofmri_fugw_wrapper, ...)
register_aligner("hyperalign", hyperalign_wrapper, ...)

# 2. Unified Alignment Interface
align(
  data,                    # AlignmentData or coercible

  method = "procrustes",   # Registered method name
  reference = "medoid",    # Reference selection strategy
  cv = "loso",             # Cross-validation scheme (none/loso/kfold)
  ...                      # Method-specific parameters
) → AlignmentResult

# 3. AlignmentData Container
AlignmentData(
  data,                    # List of per-subject data (matrices/NeuroVec)
  space,                   # gds_space or neurofunctor Domain
  subjects,                # Subject identifiers
  design = NULL,           # Optional: task structure (fmridesign)
  covariates = NULL        # Optional: subject-level covariates
)

# 4. AlignmentResult Container
AlignmentResult(
  transforms,              # Per-subject alignment transforms

  aligned_data,            # Aligned data (lazy or realized)
  reference,               # Reference space/subject

  quality,                 # Alignment quality metrics
  provenance,              # Full audit trail
  cv_folds = NULL          # Cross-validation info if applicable
)

# 5. Alignment Pipeline
pipeline <- alignment_pipeline(
  step("geometric", neurofunctor_align, target = "MNI"),
  step("functional", "procrustes", reference = "medoid"),
  step("smooth", spatial_smooth, fwhm = 6)
)
result <- execute(pipeline, data)
```

### 6.3 Module Structure

```
neuralign/
├── R/
│   ├── # Core abstractions
│   ├── alignment_data.R       # AlignmentData class
│   ├── alignment_result.R     # AlignmentResult class
│   ├── alignment_transform.R  # Unified transform representation
│   ├── alignment_space.R      # Space abstraction (wraps gds_space)
│   │
│   ├── # Method registry
│   ├── aligner_registry.R     # register_aligner(), get_aligner()
│   ├── align.R                # Main align() generic
│   │
│   ├── # Method wrappers (adapters to existing packages)
│   ├── aligner_procrustes.R   # Wraps manifoldalign
│   ├── aligner_kema.R         # Wraps manifoldalign
│   ├── aligner_gw.R           # Wraps manifoldalign
│   ├── aligner_fugw.R         # Wraps topofmri
│   ├── aligner_nef.R          # Wraps fmrireg.gnef
│   ├── aligner_dkge.R         # Wraps dkge
│   │
│   ├── # Pipeline composition
│   ├── pipeline.R             # alignment_pipeline(), step(), execute()
│   ├── compose.R              # Compose alignment stages
│   │
│   ├── # Cross-validation
│   ├── cv.R                   # Cross-validation schemes
│   ├── leakage.R              # Leakage detection/prevention
│   │
│   ├── # Quality assessment
│   ├── quality.R              # Unified quality metrics
│   ├── diagnostics.R          # Diagnostic plots
│   │
│   ├── # Integration
│   ├── compat_fmrigds.R       # MapFamily integration
│   ├── compat_fmrireg.R       # group_data integration
│   ├── compat_neurofunctor.R  # Projector integration
│   │
│   ├── # I/O
│   ├── serialize.R            # Save/load alignment artifacts
│   └── provenance.R           # Audit trail management
│
├── src/
│   └── # C++ acceleration if needed
│
├── tests/
│   └── testthat/
│
├── vignettes/
│   ├── introduction.Rmd
│   ├── alignment_methods.Rmd
│   ├── pipelines.Rmd
│   └── group_analysis.Rmd
│
├── DESCRIPTION
├── NAMESPACE
└── README.md
```

### 6.4 Integration Points

#### With fmrigds (Primary)
```r
# Convert AlignmentResult to MapFamily
map_family <- as_map_family(alignment_result)

# Use in fmrigds pipeline
gds("data.h5") %>%
  align(map_family) %>%
  reduce("meta:fe") %>%
  compute()
```

#### With fmrireg
```r
# Create aligned group_data
gd <- as_group_data(alignment_result, format = "gds")

# Use standard fmrireg analysis
fmri_meta(gd, formula = ~ 1 + group)
```

#### With neurofunctor
```r
# Compose geometric + functional alignment
geometric_proj <- compile_projector(graph, native, mni)
functional_align <- align(data_mni, "procrustes")

# Combined projector
combined <- compose_projector(geometric_proj, functional_align$projector)
```

### 6.5 Method Coverage

| Method | Underlying Package | CV Support | Use Case |
|--------|-------------------|------------|----------|
| `procrustes` | manifoldalign | Yes | Orthogonal rotation, partial obs |
| `kema` | manifoldalign | Yes | Semi-supervised manifold |
| `gw` | manifoldalign | Yes | Structure-only transport |
| `fpgw` | manifoldalign | Yes | Fused feature+structure |
| `fugw` | topofmri | Yes | Graph topology alignment |
| `nef` | fmrireg.gnef | Yes (Mode B) | Neural driver alignment |
| `dkge` | dkge | Yes (LOSO) | Design-kernel alignment |
| `cone` | manifoldalign | Yes | Graph alignment |
| `grasp` | manifoldalign | Yes | Spectral graph matching |

---

## 7. Implementation Roadmap

### Phase 1: Foundation (Core Abstractions)
1. Define AlignmentData, AlignmentResult, AlignmentTransform classes
2. Implement aligner registry with method dispatch
3. Create Procrustes wrapper (simplest case)
4. Basic I/O and provenance

### Phase 2: Method Integration
1. Wrap manifoldalign methods (KEMA, GW, FPGW)
2. Wrap topofmri FUGW
3. Wrap fmrireg.gnef NEF alignment
4. Wrap dkge alignment

### Phase 3: Pipeline Composition
1. Implement pipeline DSL
2. Integrate with neurofunctor for geometric alignment
3. Compose geometric + functional stages

### Phase 4: Statistical Safety
1. Cross-validation infrastructure
2. Leakage detection
3. Quality metrics framework

### Phase 5: Integration & Polish
1. fmrigds MapFamily integration
2. fmrireg group_data compatibility
3. Documentation and vignettes

---

## 8. Recommendations

### 8.1 For neuralign Development
1. **Start with fmrigds integration**: MapFamily provides the best existing abstraction
2. **Use manifoldalign as algorithm library**: Don't reimplement; wrap
3. **Adopt neurofunctor patterns**: Provenance, lazy evaluation, projector algebra
4. **Enforce cross-validation by default**: Statistical safety first

### 8.2 For Ecosystem Coordination
1. **Consolidate OT implementations**: manifoldalign should own Sinkhorn/GW
2. **Standardize space representations**: Converge on fmrigds gds_space
3. **Share alignment quality metrics**: Common diagnostic vocabulary
4. **Document alignment assumptions**: Each package should state MNI/native expectations

### 8.3 For Users
1. **Use neuralign for new workflows**: Unified interface, statistical safety
2. **Access underlying packages via wrappers**: Method-specific tuning available
3. **Save alignment artifacts**: Enable reproducibility and sharing

---

## Appendix: Package-Specific Key Functions

### manifoldalign
- `kema()`, `generalized_procrustes()`, `grasp()`, `cone_align()`
- `gromov_wasserstein()`, `fpgw()`, `parrot()`
- `align_many()`, `rotation_sync()`, `permutation_sync()`

### topofmri
- `fit_topo_template()`, `align_to_template()`
- `rcpp_fugw_*()`, `compute_basis()`
- `qa_template_*()`, `qa_alignment_*()

### fmrireg.gnef
- `nef_fit()`, `nef_fit_group()`, `nef_crossfit()`
- `nef_procrustes_align()`, `nef_select_medoid()`
- `nef_group_stats()`, `nef_group_stats_modeb()`

### dkge
- `dkge()`, `dkge_fit()`, `dkge_contrast()`
- `dkge_transport_*()`, `dkge_anchor_*()`, `dkge_align_*()

### fmrigds
- `gds()`, `align()`, `map_to()`, `reduce()`, `compute()`
- `MapFamily()`, `OrthogonalFamily()`, `OTFamily()`
- `UncertaintyRule()`

### neurofunctor
- `Domain()`, `BrainGraph()`, `compile_projector()`
- `Field()`, `to()`, `get_data()`, `backproject()`

### neurotransform
- `Affine3DMorphism()`, `Warp3DMorphism()`, `VolToSurfMorphism()`
- `transform()`, `compose()`, `invert()`

---

*Report generated by Claude Code analysis agents*
*Date: 2025-01-25*
