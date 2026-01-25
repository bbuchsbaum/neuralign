# Neuralign Architecture Design

## Vision

**neuralign** is the orchestration layer for functional neural alignment in R. It provides a unified interface to multiple alignment algorithms while ensuring statistical safety, proper variance propagation, and seamless integration with the fMRI analysis ecosystem.

**Core Principle**: Orchestrate, don't reimplement.

---

## 1. Design Goals

### 1.1 Unified Interface
- Single `align()` generic dispatches to any registered method
- Consistent input/output contracts across all alignment algorithms
- Method-agnostic downstream analysis

### 1.2 Composable Pipelines
- Chain geometric alignment (neurofunctor) with functional alignment
- DSL for multi-stage alignment workflows
- Lazy evaluation where possible

### 1.3 Statistical Safety
- Built-in cross-validation to prevent leakage
- Provenance tracking for audit trails
- Quality metrics for alignment validation

### 1.4 Ecosystem Integration
- First-class fmrigds MapFamily support
- fmrireg group_data compatibility
- neurofunctor Projector composition

---

## 2. Core Abstractions

### 2.1 AlignmentData

Container for data to be aligned. Wraps heterogeneous inputs into standard form.
```r
setClass("AlignmentData",
  slots = c(
    data = "list",           # Per-subject data (matrices or NeuroVec)
    space = "ANY",           # gds_space, Domain, or character
    subjects = "character",  # Subject identifiers
    design = "ANY",          # Optional fmridesign object
    covariates = "ANY",      # Optional subject-level data.frame
    metadata = "list"        # Creation timestamp, source info
  )
)

# Coercion methods
as_alignment_data.list <- function(x, subjects = NULL, space = NULL, ...)
as_alignment_data.hyperdesign <- function(x, ...)  # From manifoldalign
as_alignment_data.group_data <- function(x, ...)   # From fmrireg
as_alignment_data.gds <- function(x, ...)          # From fmrigds
```

### 2.2 AlignmentTransform

Unified representation of per-subject alignment transforms.

```r
setClass("AlignmentTransform",
  slots = c(
    type = "character",        # "orthogonal", "linear", "ot", "warp"
    by_subject = "list",       # Named list: subject_id → operator
    from_space = "ANY",        # Source space
    to_space = "ANY",          # Target space (may be same)
    invertible = "logical",    # Can be inverted?
    mass_preserving = "logical", # For OT-based transforms
    uncertainty_rule = "ANY",  # Variance propagation strategy
    metadata = "list"
  )
)

# Core operations
apply_transform(transform, data, subject_id, direction = "forward")
invert_transform(transform)
compose_transforms(t1, t2)  # t2 after t1
```

### 2.3 AlignmentResult

Output from any alignment method.

```r
setClass("AlignmentResult",
  slots = c(
    transforms = "AlignmentTransform",
    aligned_data = "ANY",      # Lazy or realized aligned data
    reference = "ANY",         # Reference subject/space/template
    method = "character",      # Method name used
    quality = "list",          # Alignment quality metrics
    cv_info = "ANY",           # Cross-validation details if applicable
    provenance = "list",       # Full audit trail
    call = "call"              # Original call for reproducibility
  )
)

# Accessors
transforms(result)
aligned(result, realize = TRUE)  # Get aligned data, optionally force realization
reference(result)
quality(result)
provenance(result)

# Coercion
as_map_family(result)     # For fmrigds integration
as_projector_list(result) # For neurofunctor integration
```

### 2.4 Aligner (Method Interface)

Protocol for alignment methods.

```r
# Registration
register_aligner(

  name,                  # Method identifier

  fit_fn,                # Function: (data, reference, ...) → AlignmentResult
  capabilities = list(  # Method properties
    supports_cv = TRUE,
    supports_partial = FALSE,
    supports_weights = FALSE,
    transform_type = "orthogonal"
  ),
  package = NULL,        # Underlying package
  description = NULL
)

# Retrieval
get_aligner(name)
list_aligners()
aligner_capabilities(name)
```

---

## 3. Main API

### 3.1 The align() Generic

```r
align(
  data,                      # AlignmentData or coercible
  method = "procrustes",     # Registered aligner name
  reference = "medoid",      # Reference selection: "medoid", "first", subject_id, or template
  cv = "none",               # Cross-validation: "none", "loso", "kfold", or custom
  cv_folds = NULL,           # Fold specification if cv = "kfold"
  weights = NULL,            # Optional subject weights
  return_aligned = TRUE,     # Include aligned data in result?
  ...                        # Method-specific parameters
) → AlignmentResult
```

**Reference Selection Strategies**:
- `"medoid"`: Subject with minimum mean distance to others (default)
- `"first"`: First subject in list
- `"centroid"`: Synthetic average reference
- `subject_id`: Specific subject as reference
- `template`: External template object

**Cross-Validation Schemes**:
- `"none"`: Full-data alignment (fast, but potential leakage)
- `"loso"`: Leave-one-subject-out (recommended for group inference)
- `"kfold"`: K-fold with user-specified folds
- Custom function: `cv_fn(data, fold) → list(train, test)`

### 3.2 Convenience Functions

```r
# Quick alignment with defaults
quick_align(data, method = "procrustes")

# Select reference subject
select_reference(data, method = "medoid", space = "scores")

# Alignment with automatic method selection
auto_align(data, goal = "group_stats")  # Chooses method based on data characteristics

# Apply pre-computed alignment to new data
apply_alignment(result, new_data, subject_id)
```

### 3.3 Pipeline DSL

```r
# Define multi-stage pipeline
pipeline <- alignment_pipeline(
  # Stage 1: Geometric alignment via neurofunctor
  stage("geometric",
    method = "neurofunctor",
    graph = brain_graph,
    target = "MNI152"
  ),

  # Stage 2: Functional alignment via Procrustes
  stage("functional",
    method = "procrustes",
    reference = "medoid",
    cv = "loso"
  ),

  # Stage 3: Optional smoothing
  stage("smooth",
    method = "spatial_smooth",
    fwhm = 6,
    enabled = TRUE
  )
)

# Execute pipeline
result <- execute(pipeline, data)

# Access intermediate results
intermediate <- result$stages$functional
```

---

## 4. Method Wrappers

### 4.1 Procrustes (manifoldalign)

```r
# R/aligner_procrustes.R
.aligner_procrustes <- function(data, reference, partial = FALSE, ...) {
  # Convert to hyperdesign if needed
  hd <- as_hyperdesign(data)

  # Delegate to manifoldalign
  if (partial) {
    fit <- manifoldalign::generalized_procrustes(hd, ...)
  } else {
    fit <- manifoldalign::orthogonal_align(hd, reference, ...)
  }

  # Convert result
  transforms <- AlignmentTransform(
    type = "orthogonal",
    by_subject = fit$O_matrices,
    from_space = space(data),
    to_space = space(data),  # Same space, different orientation
    invertible = TRUE,
    mass_preserving = TRUE
  )

  AlignmentResult(
    transforms = transforms,
    aligned_data = fit$aligned,
    reference = fit$reference,
    method = "procrustes",
    quality = list(
      mean_distance = fit$mean_distance,
      per_subject_distance = fit$distances
    ),
    provenance = list(
      package = "manifoldalign",
      function_called = "generalized_procrustes",
      timestamp = Sys.time()
    )
  )
}

register_aligner("procrustes", .aligner_procrustes,
  capabilities = list(
    supports_cv = TRUE,
    supports_partial = TRUE,
    supports_weights = FALSE,
    transform_type = "orthogonal"
  ),
  package = "manifoldalign",
  description = "Orthogonal Procrustes alignment with partial observation support"
)
```

### 4.2 KEMA (manifoldalign)

```r
# R/aligner_kema.R
.aligner_kema <- function(data, reference, y = NULL, u = 0.5, ncomp = 10, ...) {
  hd <- as_hyperdesign(data)

  fit <- manifoldalign::kema(hd, y = y, u = u, ncomp = ncomp, ...)

  # KEMA returns scores and loadings, not direct transforms
  transforms <- AlignmentTransform(
    type = "linear",
    by_subject = extract_kema_transforms(fit),
    from_space = space(data),
    to_space = "kema_latent",
    invertible = FALSE,  # Projection, not bijection
    mass_preserving = FALSE
  )

  AlignmentResult(
    transforms = transforms,
    aligned_data = fit$s,  # Aligned scores
    reference = "consensus",
    method = "kema",
    quality = list(
      variance_explained = fit$sdev^2 / sum(fit$sdev^2)
    ),
    provenance = list(
      package = "manifoldalign",
      params = list(u = u, ncomp = ncomp)
    )
  )
}
```

### 4.3 FUGW (topofmri)

```r
# R/aligner_fugw.R
.aligner_fugw <- function(data, reference = NULL, alpha = 0.5, epsilon = 0.01, ...) {
  # Convert to TopoGraph objects
  graphs <- lapply(data$data, topofmri::graph_from_neurovec, ...)

  if (is.null(reference)) {
    # Train template from data
    template <- topofmri::fit_topo_template(graphs, alpha = alpha, ...)
    reference <- template
  }

  # Align each subject to template
  embeddings <- lapply(graphs, function(g) {
    topofmri::align_to_template(g, reference, alpha = alpha, epsilon = epsilon, ...)
  })

  transforms <- AlignmentTransform(
    type = "ot",
    by_subject = lapply(embeddings, function(e) e$coupling),
    from_space = "native_graph",
    to_space = "template_graph",
    invertible = FALSE,
    mass_preserving = TRUE,
    uncertainty_rule = "ot_variance"
  )

  AlignmentResult(
    transforms = transforms,
    aligned_data = lapply(embeddings, function(e) e$Z_nodes),
    reference = reference,
    method = "fugw",
    quality = list(
      coupling_entropy = sapply(embeddings, coupling_entropy),
      coord_rmse = sapply(embeddings, function(e) e$metadata$coord_rmse)
    ),
    provenance = list(
      package = "topofmri",
      params = list(alpha = alpha, epsilon = epsilon)
    )
  )
}
```

### 4.4 NEF Alignment (fmrireg.gnef)

```r
# R/aligner_nef.R
.aligner_nef <- function(data, reference = "medoid", K = 100, mode = "A", ...) {
  # Fit NEF per subject
  fits <- lapply(data$data, fmrireg.gnef::nef_fit, K = K, ...)
  names(fits) <- data$subjects

  # Select reference
  if (reference == "medoid") {
    medoid_info <- fmrireg.gnef::nef_select_medoid(fits)
    ref_id <- medoid_info$medoid_id
  } else {
    ref_id <- reference
  }

  # Procrustes alignment in latent space
  aligned <- fmrireg.gnef::nef_procrustes_align(fits, ref_id)

  transforms <- AlignmentTransform(
    type = "orthogonal",
    by_subject = aligned$Q,
    from_space = "nef_latent",
    to_space = "nef_latent",
    invertible = TRUE,
    mass_preserving = TRUE,
    metadata = list(K = K, mode = mode)
  )

  AlignmentResult(
    transforms = transforms,
    aligned_data = aligned$S_aligned,
    reference = ref_id,
    method = "nef",
    quality = list(
      procrustes_distance = aligned$residuals
    ),
    provenance = list(
      package = "fmrireg.gnef",
      params = list(K = K, mode = mode),
      nef_fits = fits  # Store for downstream GLM
    )
  )
}
```

---

## 5. Cross-Validation Infrastructure

### 5.1 CV Schemes

```r
# R/cv.R
cv_loso <- function(n_subjects) {
  # Leave-one-subject-out
  lapply(seq_len(n_subjects), function(i) {
    list(train = setdiff(seq_len(n_subjects), i), test = i)
  })
}

cv_kfold <- function(n_subjects, k = 5, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  folds <- sample(rep(1:k, length.out = n_subjects))
  lapply(1:k, function(i) {
    list(train = which(folds != i), test = which(folds == i))
  })
}

cv_custom <- function(fold_list) {
  # Validate and return
  validate_cv_folds(fold_list)
  fold_list
}
```

### 5.2 Cross-Validated Alignment

```r
# R/align_cv.R
align_cv <- function(data, method, cv, ...) {
  aligner <- get_aligner(method)
  folds <- get_cv_folds(cv, n_subjects(data))

  results_by_fold <- lapply(folds, function(fold) {
    # Fit on training subjects
    train_data <- subset_subjects(data, fold$train)
    fit <- aligner(train_data, ...)

    # Apply to test subjects
    test_data <- subset_subjects(data, fold$test)
    test_aligned <- apply_alignment(fit, test_data)

    list(
      fold = fold,
      train_fit = fit,
      test_aligned = test_aligned
    )
  })

  # Aggregate results
  AlignmentResult(
    transforms = aggregate_cv_transforms(results_by_fold),
    aligned_data = aggregate_cv_aligned(results_by_fold),
    reference = results_by_fold[[1]]$train_fit$reference,
    method = method,
    quality = aggregate_cv_quality(results_by_fold),
    cv_info = list(
      scheme = cv,
      folds = folds,
      per_fold_results = results_by_fold
    ),
    provenance = list(
      cross_validated = TRUE,
      n_folds = length(folds)
    )
  )
}
```

### 5.3 Leakage Detection

```r
# R/leakage.R
check_leakage <- function(result) {
  warnings <- character()

  # Check if alignment used same data as downstream analysis
  if (is.null(result$cv_info) && !is.null(result$metadata$used_for_inference)) {
    warnings <- c(warnings,
      "Alignment computed on full data may cause leakage in group inference")
  }

  # Check if reference selection used test data
  if (!is.null(result$provenance$reference_selection_used_all_subjects)) {
    warnings <- c(warnings,
      "Reference selection used all subjects; consider cross-validated selection")
  }

  if (length(warnings) > 0) {
    warning("Potential leakage detected:\n", paste("- ", warnings, collapse = "\n"))
  }

  invisible(warnings)
}
```

---

## 6. Quality Assessment

### 6.1 Alignment Quality Metrics

```r
# R/quality.R
alignment_quality <- function(result, metrics = "default") {
  if (metrics == "default") {
    metrics <- c("reconstruction", "correspondence", "variance_preserved")
  }

  quality <- list()

  if ("reconstruction" %in% metrics) {
    quality$reconstruction_error <- compute_reconstruction_error(result)
  }

  if ("correspondence" %in% metrics) {
    quality$correspondence_quality <- compute_correspondence_quality(result)
  }

  if ("variance_preserved" %in% metrics) {
    quality$variance_preserved <- compute_variance_preserved(result)
  }

  if ("cross_subject_correlation" %in% metrics) {
    quality$cross_subject_correlation <- compute_cross_subject_correlation(result)
  }

  class(quality) <- "alignment_quality"
  quality
}

# Print method
print.alignment_quality <- function(x, ...) {
  cat("Alignment Quality Assessment\n")
  cat("============================\n")
  for (name in names(x)) {
    cat(sprintf("%-25s: %.3f\n", name, mean(x[[name]])))
  }
}
```

### 6.2 Diagnostic Plots

```r
# R/diagnostics.R
plot_alignment_quality <- function(result, type = "overview") {
  switch(type,
    "overview" = plot_quality_overview(result),
    "per_subject" = plot_per_subject_quality(result),
    "convergence" = plot_convergence(result),
    "correspondence" = plot_correspondence_matrix(result)
  )
}

plot_quality_overview <- function(result) {
  quality <- alignment_quality(result)
  # Create multi-panel ggplot visualization
  # ...
}
```

---

## 7. Ecosystem Integration

### 7.1 fmrigds Integration

```r
# R/compat_fmrigds.R

# Convert AlignmentResult to MapFamily
as_map_family.AlignmentResult <- function(x, name = NULL, ...) {
  if (is.null(name)) name <- paste0("neuralign_", x$method)

  fmrigds::MapFamily(
    name = name,
    from_space = x$transforms@from_space,
    to_space = x$transforms@to_space,
    type = x$transforms@type,
    by_subject = x$transforms@by_subject,
    traits = list(
      orthogonal = x$transforms@type == "orthogonal",
      mass_preserving = x$transforms@mass_preserving
    ),
    uncertainty = x$transforms@uncertainty_rule
  )
}

# Use alignment in fmrigds pipeline
gds_with_alignment <- function(gds_plan, alignment_result) {
  map_family <- as_map_family(alignment_result)
  fmrigds::align(gds_plan, map_family)
}
```

### 7.2 fmrireg Integration

```r
# R/compat_fmrireg.R

# Create aligned group_data
as_group_data.AlignmentResult <- function(x, contrast = NULL, ...) {
  # Extract aligned data as matrices
  aligned_matrices <- lapply(x$aligned_data, as.matrix)

  # Create group_data structure
  fmrireg::group_data_from_matrices(
    betas = aligned_matrices,
    subjects = names(x$transforms@by_subject),
    contrasts = contrast,
    ...
  )
}

# Extend fmri_meta for aligned data
fmri_meta_aligned <- function(alignment_result, formula = ~ 1, ...) {
  gd <- as_group_data(alignment_result)
  fmrireg::fmri_meta(gd, formula = formula, ...)
}
```

### 7.3 neurofunctor Integration

```r
# R/compat_neurofunctor.R

# Compose geometric projector with functional alignment
compose_with_projector <- function(projector, alignment_result) {
  # Create alignment projector
  align_proj <- as_projector(alignment_result)

  # Compose: geometric first, then functional
  neurofunctor::compose_projectors(projector, align_proj)
}

# Use neurofunctor graph for geometric stage
geometric_align <- function(data, graph, target_domain) {
  projectors <- lapply(names(data$data), function(subj) {
    source <- get_subject_domain(graph, subj)
    neurofunctor::compile_projector(graph, source, target_domain)
  })

  AlignmentResult(
    transforms = projectors_to_transforms(projectors),
    aligned_data = apply_projectors(projectors, data$data),
    reference = target_domain,
    method = "neurofunctor_geometric"
  )
}
```

---

## 8. Serialization

### 8.1 Save/Load

```r
# R/serialize.R

save_alignment <- function(result, path, format = "rds") {
  # Add metadata
  result$provenance$saved_at <- Sys.time()
  result$provenance$neuralign_version <- packageVersion("neuralign")

  switch(format,
    "rds" = saveRDS(result, path),
    "h5" = save_alignment_h5(result, path),
    stop("Unknown format: ", format)
  )

  invisible(path)
}

load_alignment <- function(path) {
  ext <- tools::file_ext(path)
  result <- switch(ext,
    "rds" = readRDS(path),
    "h5" = load_alignment_h5(path),
    stop("Unknown format: ", ext)
  )

  validate_alignment_result(result)
  result
}
```

### 8.2 Artifact Format (HDF5)

```
/neuralign/
├── version                    # Package version
├── method                     # Alignment method name
├── reference/                 # Reference information
│   ├── type                   # "subject", "template", "consensus"
│   └── data                   # Reference data if applicable
├── transforms/
│   ├── type                   # Transform type
│   ├── from_space             # Source space descriptor
│   ├── to_space               # Target space descriptor
│   └── by_subject/            # Per-subject transforms
│       ├── sub-01/
│       │   └── operator       # Transform matrix/function
│       └── sub-02/
│           └── operator
├── quality/
│   └── [metrics]              # Quality assessment results
├── cv_info/                   # Cross-validation details (if applicable)
│   ├── scheme
│   ├── n_folds
│   └── folds/
└── provenance/
    ├── timestamp
    ├── call
    └── parameters
```

---

## 9. Dependencies

### 9.1 Required
- `methods`: S4 classes
- `Matrix`: Sparse matrices
- `digest`: Hashing for provenance

### 9.2 Suggested (for specific methods)
- `manifoldalign`: Procrustes, KEMA, GW, FPGW, etc.
- `topofmri`: FUGW, topological alignment
- `fmrireg.gnef`: NEF alignment
- `dkge`: Design-kernel alignment
- `fmrigds`: Group data integration
- `neurofunctor`: Geometric alignment composition

### 9.3 Optional (for enhanced functionality)
- `ggplot2`: Diagnostic plots
- `hdf5r`: HDF5 serialization
- `future`: Parallel execution

---

## 10. Testing Strategy

### 10.1 Unit Tests
- AlignmentData coercion methods
- AlignmentTransform operations (apply, invert, compose)
- Each aligner wrapper in isolation

### 10.2 Integration Tests
- Full pipeline execution
- Cross-validation workflows
- Package interoperability (fmrigds, fmrireg, neurofunctor)

### 10.3 Statistical Tests
- Verify variance propagation correctness
- Verify cross-validation prevents leakage
- Compare to reference implementations

### 10.4 Performance Tests
- Scaling with number of subjects
- Scaling with data dimensionality
- Memory usage profiling

---

*Architecture designed based on synthesis of 9 package analyses*
*neuralign v0.1.0 target*
