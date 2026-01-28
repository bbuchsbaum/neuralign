test_that("fit_alignment works with procrustes", {
  # Register procrustes for testing
  ensure_test_aligner("procrustes")

  # Create test data
  set.seed(123)
  n_features <- 20
  n_obs <- 30

  # Create data with shared structure
  reference_pattern <- matrix(rnorm(n_features * n_obs), n_features, n_obs)

  data_list <- lapply(1:4, function(i) {
    # Add subject-specific noise (procrustes will align in feature space)
    reference_pattern + matrix(rnorm(n_features * n_obs, sd = 0.5), n_features, n_obs)
  })
  names(data_list) <- paste0("sub-0", 1:4)

  adat <- AlignmentData(data_list)

  # Fit alignment
  result <- fit_alignment(adat, method = "procrustes", reference = "consensus")

  expect_s4_class(result, "AlignmentResult")
  expect_equal(length(result@aligned), 4)
  expect_equal(names(result@aligned), names(data_list))

  # Check model
  model <- get_model(result)
  expect_s4_class(model, "AlignmentModel")
  expect_equal(model@method, "procrustes")
})

test_that("fit_alignment threads obs_labels into AlignmentData coercion", {
  ensure_test_aligner("procrustes")

  set.seed(321)
  data_list <- make_test_data_list(n_subjects = 2, n_features = 5, n_obs = 10)
  labs <- paste0("reg", seq_len(10))

  res <- fit_alignment(data_list, method = "procrustes", obs_labels = labs, reference = "sub-01")
  expect_s4_class(res, "AlignmentResult")
  expect_equal(get_obs_labels(as_alignment_data(data_list, obs_labels = labs)), labs)

  # Mismatched labels should error (validation auto-enables when obs_labels present)
  expect_error(
    fit_alignment(data_list, method = "procrustes", obs_labels = c("a", "b"), reference = "sub-01"),
    "obs_labels length mismatch"
  )
})

test_that("fit_alignment with medoid reference works", {
  ensure_test_aligner("procrustes")

  set.seed(456)
  data_list <- make_test_data_list(n_subjects = 3, n_features = 10, n_obs = 10)
  adat <- AlignmentData(data_list)

  result <- fit_alignment(adat, method = "procrustes", reference = "medoid")

  expect_s4_class(result, "AlignmentResult")
  expect_true(result@model@reference %in% adat@subjects)
})

test_that("fit_alignment with specific reference subject works", {
  ensure_test_aligner("procrustes")

  set.seed(789)
  data_list <- make_test_data_list(n_subjects = 3, n_features = 10, n_obs = 10)
  adat <- AlignmentData(data_list)

  result <- fit_alignment(adat, method = "procrustes", reference = "sub-02")

  expect_equal(result@model@reference, "sub-02")
})

test_that("fit_alignment errors on unknown method", {
  data_list <- list(
    "sub-01" = matrix(1, 10, 10),
    "sub-02" = matrix(1, 10, 10)
  )
  adat <- AlignmentData(data_list)

  expect_error(
    fit_alignment(adat, method = "nonexistent_method"),
    "Unknown method"
  )
})

test_that("fit_alignment produces transforms with correct dimensions", {
  ensure_test_aligner("procrustes")

  n_features <- 15
  data_list <- make_test_data_list(n_subjects = 3, n_features = n_features, n_obs = 20)
  adat <- AlignmentData(data_list)

  result <- fit_alignment(adat, method = "procrustes")
  model <- get_model(result)

  for (subj in names(model@transforms)) {
    transform <- model@transforms[[subj]]
    # Transforms should be (target x source) = (n_features x n_features)
    expect_equal(dim(transform), c(n_features, n_features))
  }
})

test_that("fit_alignment with train_idx works", {
  ensure_test_aligner("procrustes")

  set.seed(111)
  data_list <- make_test_data_list(n_subjects = 4, n_features = 10, n_obs = 10)
  adat <- AlignmentData(data_list)

  # Fit on first 3 subjects only
  result <- fit_alignment(
    adat,
    method = "procrustes",
    train_idx = 1:3
  )

  # Should still produce transforms for all subjects
  expect_equal(length(result@model@transforms), 4)

  # Train subjects should be recorded
  expect_equal(result@model@train_subjects, c("sub-01", "sub-02", "sub-03"))
})

test_that("fit_alignment computes quality metrics", {
  ensure_test_aligner("procrustes")

  set.seed(222)
  data_list <- make_test_data_list(n_subjects = 3, n_features = 10, n_obs = 10)
  adat <- AlignmentData(data_list)

  result <- fit_alignment(adat, method = "procrustes", compute_quality = TRUE)

  expect_true(length(result@quality) > 0)
  expect_true("mean_pairwise_correlation" %in% names(result@quality))
})


# ---------- New tests appended below ----------

test_that("fit_alignment errors when AlignmentData has different obs_labels", {
  ensure_test_aligner("procrustes")

  set.seed(42)
  n_feat <- 10
  n_obs <- 5
  data_list <- make_test_data_list(n_subjects = 2, n_features = n_feat, n_obs = n_obs)

  # Create AlignmentData with obs_labels already set
  adat <- AlignmentData(data_list, obs_labels = make_test_obs_labels(n_obs, prefix = "obs"))

  # Passing conflicting obs_labels should error

  expect_error(
    fit_alignment(adat, method = "procrustes", obs_labels = paste0("different", 1:n_obs)),
    "already has different obs_labels"
  )
})

test_that("fit_alignment sets obs_labels when AlignmentData has none", {
  ensure_test_aligner("procrustes")

  set.seed(42)
  n_feat <- 10
  n_obs <- 5
  data_list <- make_test_data_list(n_subjects = 2, n_features = n_feat, n_obs = n_obs)

  # Create AlignmentData without obs_labels
  adat <- AlignmentData(data_list)
  expect_null(adat@obs_labels)

  # Pass obs_labels via fit_alignment; should succeed and set them
  labs <- make_test_obs_labels(n_obs, prefix = "obs")
  result <- fit_alignment(adat, method = "procrustes", obs_labels = labs, reference = "sub-01")
  expect_s4_class(result, "AlignmentResult")
})

test_that("fit_alignment warns when method does not support CV", {
  saved <- as.list(neuralign:::.aligner_registry)
  neuralign:::.clear_registry()
  on.exit({
    neuralign:::.clear_registry()
    for (nm in names(saved)) {
      assign(nm, saved[[nm]], envir = neuralign:::.aligner_registry)
    }
  }, add = TRUE)

  # Register a test aligner with supports_cv = FALSE
  test_fit <- function(data, reference, train_idx = NULL, ...) {
    n_feat <- nrow(data@data[[1]])
    n_obs <- ncol(data@data[[1]])
    if (is.null(train_idx)) train_idx <- seq_along(data@subjects)
    transforms <- lapply(data@subjects, function(s) diag(n_feat))
    names(transforms) <- data@subjects
    list(
      transforms = transforms,
      reference_data = matrix(0, n_feat, n_obs),
      space_from = NULL,
      space_to = NULL
    )
  }

  register_aligner(
    "no_cv_method",
    test_fit,
    capabilities = list(supports_cv = FALSE)
  )
  on.exit(unregister_aligner("no_cv_method"))

  set.seed(42)
  n_feat <- 10
  n_obs <- 5
  data_list <- lapply(1:4, function(i) matrix(rnorm(n_feat * n_obs), n_feat, n_obs))
  names(data_list) <- paste0("sub-0", 1:4)
  adat <- AlignmentData(data_list)

  # Use "consensus" reference to avoid the validate_cv_setup leakage error
  expect_warning(
    fit_alignment(adat, method = "no_cv_method", cv = "loso", reference = "consensus"),
    "may not fully support CV"
  )
})

test_that("fit_alignment restrict_to_identified preserves aligned data for orthogonal methods", {
  ensure_test_aligner("procrustes")

  set.seed(999)
  n_feat <- 20
  n_obs <- 5 # rank-deficient relative to n_feat
  data_list <- make_test_data_list(n_subjects = 3, n_features = n_feat, n_obs = n_obs)
  adat <- AlignmentData(data_list)

  res0 <- fit_alignment(
    adat,
    method = "procrustes",
    reference = "sub-01",
    compute_quality = FALSE
  )
  res1 <- fit_alignment(
    adat,
    method = "procrustes",
    reference = "sub-01",
    compute_quality = FALSE,
    restrict_to_identified = TRUE,
    restrict_tol = 1e-12
  )

  for (subj in names(res0@aligned)) {
    expect_equal(res1@aligned[[subj]], res0@aligned[[subj]], tolerance = 1e-10)
  }

  model1 <- get_model(res1)
  expect_true(isTRUE(model1@method_state$restrict_to_identified))
})

test_that("fit_alignment restrict_to_identified errors for non-orthogonal operator methods", {
  saved <- as.list(neuralign:::.aligner_registry)
  neuralign:::.clear_registry()
  on.exit({
    neuralign:::.clear_registry()
    for (nm in names(saved)) {
      assign(nm, saved[[nm]], envir = neuralign:::.aligner_registry)
    }
  }, add = TRUE)

  test_fit <- function(data, reference, train_idx = NULL, ...) {
    n_feat <- nrow(data@data[[1]])
    transforms <- lapply(data@subjects, function(s) diag(n_feat))
    names(transforms) <- data@subjects
    list(transforms = transforms, reference_data = data@data[[1]])
  }

  register_aligner(
    "linear_method",
    test_fit,
    capabilities = list(returns = "operator", transform_type = "linear", supports_cv = TRUE)
  )

  data_list <- make_test_data_list(n_subjects = 2, n_features = 6, n_obs = 3)
  adat <- AlignmentData(data_list)

  expect_error(
    fit_alignment(adat, method = "linear_method", restrict_to_identified = TRUE),
    "only supported for orthogonal operators"
  )
})

test_that(".validate_cv_folds_spec errors on invalid fold specs", {
  n_subjects <- 4

  # Not a list / missing $folds
  expect_error(
    neuralign:::.validate_cv_folds_spec("not_a_list", n_subjects),
    "must be a fold spec list"
  )
  expect_error(
    neuralign:::.validate_cv_folds_spec(list(something = "else"), n_subjects),
    "must be a fold spec list"
  )

  # Fewer than 2 folds
  expect_error(
    neuralign:::.validate_cv_folds_spec(
      list(folds = list(fold1 = list(train = 1:3, test = 4))),
      n_subjects
    ),
    "must contain >= 2 folds"
  )

  # Fold missing $train or $test
  expect_error(
    neuralign:::.validate_cv_folds_spec(
      list(folds = list(
        fold1 = list(train = 1:3),
        fold2 = list(train = c(1, 4), test = c(2, 3))
      )),
      n_subjects
    ),
    "must contain \\$train and \\$test"
  )

  # NA indices
  expect_error(
    neuralign:::.validate_cv_folds_spec(
      list(folds = list(
        fold1 = list(train = c(1, NA), test = 3),
        fold2 = list(train = c(1, 3), test = 2)
      )),
      n_subjects
    ),
    "has NA indices"
  )

  # Out-of-range indices
  expect_error(
    neuralign:::.validate_cv_folds_spec(
      list(folds = list(
        fold1 = list(train = c(1, 2), test = 99),
        fold2 = list(train = c(1, 3), test = 2)
      )),
      n_subjects
    ),
    "has out-of-range indices"
  )

  # Overlapping train/test
  expect_error(
    neuralign:::.validate_cv_folds_spec(
      list(folds = list(
        fold1 = list(train = c(1, 2, 3), test = c(3, 4)),
        fold2 = list(train = c(1, 3), test = 2)
      )),
      n_subjects
    ),
    "has overlapping train/test"
  )
})

test_that(".reference_kind classifies references correctly", {
  # Matrix -> "template"
  expect_equal(neuralign:::.reference_kind(matrix(1:4, 2, 2)), "template")

  # Subject name (not one of the special values) -> "fixed_subject"
  expect_equal(neuralign:::.reference_kind("sub-01"), "fixed_subject")

  # "medoid", "centroid", "consensus" -> "data_driven"
  expect_equal(neuralign:::.reference_kind("medoid"), "data_driven")
  expect_equal(neuralign:::.reference_kind("centroid"), "data_driven")
  expect_equal(neuralign:::.reference_kind("consensus"), "data_driven")

  # Numeric -> "unknown"
  expect_equal(neuralign:::.reference_kind(42), "unknown")
})

test_that(".fit_new_subject uses apply_fn when available", {
  saved <- as.list(neuralign:::.aligner_registry)
  neuralign:::.clear_registry()
  on.exit({
    neuralign:::.clear_registry()
    for (nm in names(saved)) {
      assign(nm, saved[[nm]], envir = neuralign:::.aligner_registry)
    }
  }, add = TRUE)

  set.seed(42)
  n_feat <- 10
  n_obs <- 5
  data_list <- lapply(1:3, function(i) matrix(rnorm(n_feat * n_obs), n_feat, n_obs))
  names(data_list) <- paste0("sub-0", 1:3)
  adat <- AlignmentData(data_list)

  # Create a mock aligner with apply_fn
  test_fit <- function(data, reference, train_idx = NULL, ...) {
    transforms <- lapply(data@subjects, function(s) diag(n_feat))
    names(transforms) <- data@subjects
    list(
      transforms = transforms,
      reference_data = matrix(0, n_feat, n_obs),
      space_from = NULL,
      space_to = NULL
    )
  }

  apply_fn_called <- FALSE
  test_apply <- function(fit_result, new_data, ...) {
    # Return identity transform scaled by 2 to distinguish from fit path
    subj <- new_data@subjects[1]
    list(transforms = setNames(list(diag(n_feat) * 2), subj))
  }

  # Build a mock aligner list (as returned by get_aligner)
  aligner <- list(
    name = "apply_test",
    fit_fn = test_fit,
    apply_fn = test_apply,
    capabilities = list(supports_cv = TRUE)
  )

  fit_result <- test_fit(adat, reference = "consensus")
  reference <- "consensus"

  result <- neuralign:::.fit_new_subject(aligner, fit_result, adat, 3, reference)

  # Should be the 2x identity (from apply_fn), not the 1x identity (from fit_fn)
  expect_equal(result, diag(n_feat) * 2)
})

test_that(".fit_new_subject falls back to fit_fn without apply_fn and uses reference_data", {
  saved <- as.list(neuralign:::.aligner_registry)
  neuralign:::.clear_registry()
  on.exit({
    neuralign:::.clear_registry()
    for (nm in names(saved)) {
      assign(nm, saved[[nm]], envir = neuralign:::.aligner_registry)
    }
  }, add = TRUE)

  set.seed(42)
  n_feat <- 10
  n_obs <- 5
  data_list <- lapply(1:3, function(i) matrix(rnorm(n_feat * n_obs), n_feat, n_obs))
  names(data_list) <- paste0("sub-0", 1:3)
  adat <- AlignmentData(data_list)

  ref_data <- matrix(rnorm(n_feat * n_obs), n_feat, n_obs)

  # fit_fn that records which reference was passed
  captured_ref <- NULL
  test_fit <- function(data, reference, train_idx = NULL, ...) {
    captured_ref <<- reference
    transforms <- lapply(data@subjects, function(s) diag(n_feat))
    names(transforms) <- data@subjects
    list(
      transforms = transforms,
      reference_data = ref_data,
      space_from = NULL,
      space_to = NULL
    )
  }

  aligner <- list(
    name = "no_apply_test",
    fit_fn = test_fit,
    apply_fn = NULL,
    capabilities = list(supports_cv = TRUE)
  )

  fit_result <- test_fit(adat, reference = "consensus")
  captured_ref <- NULL  # reset

  result <- neuralign:::.fit_new_subject(aligner, fit_result, adat, 3, "consensus")

  # The fallback path should use reference_data from fit_result, not "consensus"
  expect_true(is.matrix(captured_ref))
  expect_equal(captured_ref, ref_data)
  expect_equal(result, diag(n_feat))
})

test_that(".fit_new_subject falls back to original reference when reference_data is NULL", {
  saved <- as.list(neuralign:::.aligner_registry)
  neuralign:::.clear_registry()
  on.exit({
    neuralign:::.clear_registry()
    for (nm in names(saved)) {
      assign(nm, saved[[nm]], envir = neuralign:::.aligner_registry)
    }
  }, add = TRUE)

  set.seed(42)
  n_feat <- 10
  n_obs <- 5
  data_list <- lapply(1:3, function(i) matrix(rnorm(n_feat * n_obs), n_feat, n_obs))
  names(data_list) <- paste0("sub-0", 1:3)
  adat <- AlignmentData(data_list)

  captured_ref <- NULL
  test_fit <- function(data, reference, train_idx = NULL, ...) {
    captured_ref <<- reference
    transforms <- lapply(data@subjects, function(s) diag(n_feat))
    names(transforms) <- data@subjects
    list(
      transforms = transforms,
      reference_data = NULL,  # no reference_data
      space_from = NULL,
      space_to = NULL
    )
  }

  aligner <- list(
    name = "null_ref_test",
    fit_fn = test_fit,
    apply_fn = NULL,
    capabilities = list(supports_cv = TRUE)
  )

  fit_result <- test_fit(adat, reference = "consensus")
  captured_ref <- NULL  # reset

  result <- neuralign:::.fit_new_subject(aligner, fit_result, adat, 3, "consensus")

  # Should fall back to the original reference since reference_data is NULL
  expect_equal(captured_ref, "consensus")
  expect_equal(result, diag(n_feat))
})


# ---------- return_aligned = FALSE (model-only workflow) ----------

test_that("fit_alignment with return_aligned=FALSE returns empty aligned list", {
  ensure_test_aligner("procrustes")

  set.seed(42)
  n_feat <- 10
  n_obs <- 15
  data_list <- list(
    "sub-01" = matrix(rnorm(n_feat * n_obs), n_feat, n_obs),
    "sub-02" = matrix(rnorm(n_feat * n_obs), n_feat, n_obs),
    "sub-03" = matrix(rnorm(n_feat * n_obs), n_feat, n_obs)
  )
  adat <- AlignmentData(data_list)

  result <- fit_alignment(adat, method = "procrustes", reference = "consensus",
                          return_aligned = FALSE)

  expect_s4_class(result, "AlignmentResult")

  # Model should be fully populated
  model <- get_model(result)
  expect_equal(length(model@transforms), 3)

  # Aligned list should be empty
  expect_equal(length(result@aligned), 0)
})

test_that("fit_alignment with return_aligned=TRUE returns aligned data", {
  ensure_test_aligner("procrustes")

  set.seed(42)
  n_feat <- 10
  n_obs <- 15
  data_list <- list(
    "sub-01" = matrix(rnorm(n_feat * n_obs), n_feat, n_obs),
    "sub-02" = matrix(rnorm(n_feat * n_obs), n_feat, n_obs)
  )
  adat <- AlignmentData(data_list)

  result <- fit_alignment(adat, method = "procrustes", reference = "consensus",
                          return_aligned = TRUE)

  expect_equal(length(result@aligned), 2)
  expect_equal(names(result@aligned), c("sub-01", "sub-02"))
})

test_that("fit_alignment cv='loso' with return_aligned=FALSE returns empty aligned", {
  ensure_test_aligner("procrustes")

  set.seed(42)
  n_feat <- 10
  n_obs <- 15
  data_list <- list(
    "sub-01" = matrix(rnorm(n_feat * n_obs), n_feat, n_obs),
    "sub-02" = matrix(rnorm(n_feat * n_obs), n_feat, n_obs),
    "sub-03" = matrix(rnorm(n_feat * n_obs), n_feat, n_obs)
  )
  adat <- AlignmentData(data_list)

  result <- fit_alignment(adat, method = "procrustes", reference = "consensus",
                          cv = "loso", return_aligned = FALSE)

  expect_s4_class(result, "AlignmentResult")
  expect_equal(length(result@aligned), 0)
  # Model should still exist
  expect_equal(length(get_model(result)@transforms), 3)
})


# ---------- Observation-axis CV ----------

test_that(".subset_obs works with shared observation indices", {
  set.seed(42)
  data_list <- list(
    "sub-01" = matrix(1:30, 5, 6),
    "sub-02" = matrix(31:60, 5, 6)
  )
  adat <- AlignmentData(data_list)

  result <- neuralign:::.subset_obs(adat, c(1, 3, 5))

  result_list <- get_data_list(result)
  expect_equal(ncol(result_list[["sub-01"]]), 3)
  expect_equal(ncol(result_list[["sub-02"]]), 3)
  expect_equal(result_list[["sub-01"]], data_list[["sub-01"]][, c(1, 3, 5)])
})

test_that(".subset_obs works with per-subject observation indices", {
  set.seed(42)
  data_list <- list(
    "sub-01" = matrix(1:30, 5, 6),
    "sub-02" = matrix(31:60, 5, 6)
  )
  adat <- AlignmentData(data_list)

  obs_idx <- list(
    "sub-01" = c(1, 2),
    "sub-02" = c(3, 4, 5)
  )
  result <- neuralign:::.subset_obs(adat, obs_idx)

  result_list <- get_data_list(result)
  expect_equal(ncol(result_list[["sub-01"]]), 2)
  expect_equal(ncol(result_list[["sub-02"]]), 3)
})

test_that(".subset_obs errors on empty observation index", {
  data_list <- list(
    "sub-01" = matrix(1:30, 5, 6),
    "sub-02" = matrix(31:60, 5, 6)
  )
  adat <- AlignmentData(data_list)

  obs_idx <- list(
    "sub-01" = integer(0),
    "sub-02" = c(1, 2)
  )
  expect_error(
    neuralign:::.subset_obs(adat, obs_idx),
    "Empty observation index"
  )
})

test_that("fit_alignment routes observation-axis folds to .fit_cv_obs_folds", {
  ensure_test_aligner("procrustes")

  set.seed(42)
  n_feat <- 10
  n_obs <- 20
  data_list <- list(
    "sub-01" = matrix(rnorm(n_feat * n_obs), n_feat, n_obs),
    "sub-02" = matrix(rnorm(n_feat * n_obs), n_feat, n_obs),
    "sub-03" = matrix(rnorm(n_feat * n_obs), n_feat, n_obs)
  )
  adat <- AlignmentData(data_list)

  # Build observation-axis fold spec
  obs_folds <- list(
    axis = "observation",
    method = "custom",
    folds = list(
      fold1 = list(train_idx = 1:10, test_idx = 11:20),
      fold2 = list(train_idx = 11:20, test_idx = 1:10)
    )
  )

  result <- fit_alignment(adat, method = "procrustes", reference = "sub-01",
                          cv_folds = obs_folds)

  expect_s4_class(result, "AlignmentResult")
  expect_equal(result@cv_info$axis, "observation")
  expect_equal(result@cv_info$n_folds, 2)
  expect_equal(length(get_model(result)@transforms), 3)
})

test_that("observation-axis CV warns when method doesn't declare observation in cv_axes", {
  saved <- as.list(neuralign:::.aligner_registry)
  neuralign:::.clear_registry()
  on.exit({
    neuralign:::.clear_registry()
    for (nm in names(saved)) {
      assign(nm, saved[[nm]], envir = neuralign:::.aligner_registry)
    }
  }, add = TRUE)

  test_fit <- function(data, reference, train_idx = NULL, ...) {
    n_feat <- nrow(data@data[[1]])
    transforms <- lapply(data@subjects, function(s) diag(n_feat))
    names(transforms) <- data@subjects
    list(
      transforms = transforms,
      reference_data = matrix(0, n_feat, 1),
      space_from = NULL,
      space_to = NULL
    )
  }

  register_aligner(
    "subj_only_method",
    test_fit,
    capabilities = list(cv_axes = c("subject"))
  )
  on.exit(unregister_aligner("subj_only_method"))

  set.seed(42)
  data_list <- list(
    "sub-01" = matrix(rnorm(40), 4, 10),
    "sub-02" = matrix(rnorm(40), 4, 10)
  )
  adat <- AlignmentData(data_list)

  obs_folds <- list(
    axis = "observation",
    method = "custom",
    folds = list(
      fold1 = list(train_idx = 1:5, test_idx = 6:10),
      fold2 = list(train_idx = 6:10, test_idx = 1:5)
    )
  )

  expect_warning(
    fit_alignment(adat, method = "subj_only_method", reference = "sub-01",
                  cv_folds = obs_folds),
    "observation-axis CV may not be meaningful"
  )
})

test_that("fit_alignment validates reference against method reference_types", {
  ensure_test_aligner("procrustes_graph")

  set.seed(123)
  x <- matrix(rnorm(4 * 6), 4, 6)
  adat <- AlignmentData(
    list(s1 = x, s2 = x),
    obs_labels = paste0("o", seq_len(ncol(x)))
  )

  expect_error(
    fit_alignment(adat, method = "procrustes_graph", reference = "consensus"),
    "does not support reference type 'consensus'"
  )

  expect_error(
    fit_alignment(adat, method = "procrustes_graph", reference = matrix(0, 4, 6)),
    "does not support reference type 'template'"
  )
})


# ---------- .reference_type_for_validation edge cases ----------

test_that(".reference_type_for_validation classifies matrix as template", {
  adat <- AlignmentData(list(s1 = matrix(1, 3, 3)))
  ref_type <- neuralign:::.reference_type_for_validation(matrix(1, 3, 3), adat)
  expect_equal(ref_type, "template")
})

test_that(".reference_type_for_validation classifies medoid/centroid as subject", {
  adat <- AlignmentData(list(s1 = matrix(1, 3, 3)))
  expect_equal(neuralign:::.reference_type_for_validation("medoid", adat), "subject")
  expect_equal(neuralign:::.reference_type_for_validation("centroid", adat), "subject")
})

test_that(".reference_type_for_validation classifies consensus", {
  adat <- AlignmentData(list(s1 = matrix(1, 3, 3)))
  expect_equal(neuralign:::.reference_type_for_validation("consensus", adat), "consensus")
})

test_that(".reference_type_for_validation returns unknown for barycenter", {
  adat <- AlignmentData(list(s1 = matrix(1, 3, 3)))
  expect_equal(neuralign:::.reference_type_for_validation("barycenter", adat), "unknown")
})

test_that(".reference_type_for_validation classifies subject ID", {
  adat <- AlignmentData(list(s1 = matrix(1, 3, 3), s2 = matrix(1, 3, 3)))
  expect_equal(neuralign:::.reference_type_for_validation("s1", adat), "subject")
  expect_equal(neuralign:::.reference_type_for_validation("s2", adat), "subject")
})

test_that(".reference_type_for_validation returns unknown for unrecognized string", {
  adat <- AlignmentData(list(s1 = matrix(1, 3, 3)))
  expect_equal(neuralign:::.reference_type_for_validation("bogus_ref", adat), "unknown")
})

test_that(".reference_type_for_validation returns unknown for non-character non-matrix", {
  adat <- AlignmentData(list(s1 = matrix(1, 3, 3)))
  expect_equal(neuralign:::.reference_type_for_validation(42, adat), "unknown")
})


# ---------- .validate_reference_for_aligner ----------

test_that(".validate_reference_for_aligner passes when no reference_types constraint", {
  caps <- list(reference_types = NULL)
  adat <- AlignmentData(list(s1 = matrix(1, 3, 3)))
  expect_invisible(
    neuralign:::.validate_reference_for_aligner("s1", adat, "test", caps)
  )
})

test_that(".validate_reference_for_aligner errors on unknown string reference", {
  caps <- list(reference_types = c("subject"))
  adat <- AlignmentData(list(s1 = matrix(1, 3, 3)))
  expect_error(
    neuralign:::.validate_reference_for_aligner("not_a_subject", adat, "test", caps),
    "Unknown reference"
  )
})

test_that(".validate_reference_for_aligner errors on non-character non-matrix reference", {
  caps <- list(reference_types = c("subject"))
  adat <- AlignmentData(list(s1 = matrix(1, 3, 3)))
  expect_error(
    neuralign:::.validate_reference_for_aligner(42, adat, "test", caps),
    "Invalid reference specification"
  )
})

test_that(".validate_reference_for_aligner errors when ref type not in supported set", {
  caps <- list(reference_types = c("subject"))
  adat <- AlignmentData(list(s1 = matrix(1, 3, 3)))
  expect_error(
    neuralign:::.validate_reference_for_aligner("consensus", adat, "test_method", caps),
    "does not support reference type.*consensus"
  )
})


# ---------- Embedding + low-rank transforms ----------

test_that("fit_alignment supports embedding-returning aligners (cv='none' only)", {
  with_temp_registry(code = {
    embed_fit <- function(data, reference, train_idx = NULL, ...) {
      k <- 3
      aligned <- lapply(data@subjects, function(s) {
        X <- get_subject_data(data, s)
        X[seq_len(k), , drop = FALSE]
      })
      names(aligned) <- data@subjects
      list(
        aligned = aligned,
        reference_data = NULL,
        space_from = data@space,
        space_to = "latent"
      )
    }

    register_aligner(
      "embedder",
      embed_fit,
      capabilities = list(
        returns = "embedding",
        supports_new_data = FALSE
      )
    )

    set.seed(1)
    adat <- make_test_alignment_data(n_subjects = 2, n_features = 6, n_obs = 5)
    res <- fit_alignment(adat, method = "embedder", reference = "sub-01", compute_quality = FALSE)

    expect_s4_class(res, "AlignmentResult")
    expect_equal(names(res@aligned), adat@subjects)
    expect_equal(dim(res@aligned[[1L]]), c(3, 5))

    model <- get_model(res)
    expect_true(all(vapply(model@transforms, inherits, logical(1), "neuralign_embedding_transform")))

    expect_error(
      fit_alignment(adat, method = "embedder", cv = "loso"),
      "cross-validation"
    )
  })
})

test_that("fit_alignment can apply low-rank operator transforms", {
  with_temp_registry(code = {
    low_rank_fit <- function(data, reference, train_idx = NULL, ...) {
      k <- 2
      transforms <- lapply(data@subjects, function(s) {
        p <- nrow(get_subject_data(data, s))
        U <- matrix(runif(p * k), p, k)
        V <- matrix(runif(p * k), p, k)
        neuralign:::.new_low_rank_transform(U, V)
      })
      names(transforms) <- data@subjects
      list(
        transforms = transforms,
        reference_data = NULL,
        space_from = data@space,
        space_to = data@space
      )
    }

    register_aligner(
      "low_rank_method",
      low_rank_fit,
      capabilities = list(returns = "operator")
    )

    set.seed(2)
    adat <- make_test_alignment_data(n_subjects = 2, n_features = 5, n_obs = 4)
    res <- fit_alignment(adat, method = "low_rank_method", compute_quality = FALSE)

    expect_s4_class(res, "AlignmentResult")
    expect_true(inherits(get_model(res)@transforms[[1L]], "neuralign_low_rank_transform"))
    expect_equal(dim(res@aligned[[1L]]), c(5, 4))
  })
})
