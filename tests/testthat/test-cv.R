test_that("create_cv_folds with loso works", {
  data_list <- list(
    "sub-01" = matrix(1, 5, 5),
    "sub-02" = matrix(1, 5, 5),
    "sub-03" = matrix(1, 5, 5)
  )
  adat <- AlignmentData(data_list)

  folds <- create_cv_folds(adat, method = "loso")

  expect_equal(folds$method, "loso")
  expect_equal(folds$n_folds, 3)
  expect_equal(length(folds$folds), 3)

  # Each fold should have 1 test subject
  for (fold in folds$folds) {
    expect_equal(length(fold$test), 1)
    expect_equal(length(fold$train), 2)
  }
})

test_that("create_cv_folds with kfold works", {
  data_list <- list(
    "sub-01" = matrix(1, 5, 5),
    "sub-02" = matrix(1, 5, 5),
    "sub-03" = matrix(1, 5, 5),
    "sub-04" = matrix(1, 5, 5)
  )
  adat <- AlignmentData(data_list)

  folds <- create_cv_folds(adat, method = "kfold", k = 2, seed = 123)

  expect_equal(folds$method, "kfold")
  expect_equal(folds$n_folds, 2)

  # All subjects should be assigned
  expect_equal(length(folds$assignments), 4)
  expect_true(all(folds$assignments %in% 1:2))

  # Each fold should have some subjects
  for (fold in folds$folds) {
    expect_true(length(fold$test) > 0)
    expect_true(length(fold$train) > 0)
  }
})

test_that("create_cv_folds falls back to loso when k > n", {
  data_list <- list(
    "sub-01" = matrix(1, 5, 5),
    "sub-02" = matrix(1, 5, 5)
  )
  adat <- AlignmentData(data_list)

  expect_warning(
    folds <- create_cv_folds(adat, method = "kfold", k = 5),
    "LOSO"
  )

  expect_equal(folds$method, "loso")
})

test_that("create_cv_folds is reproducible with seed", {
  data_list <- list(
    "sub-01" = matrix(1, 5, 5),
    "sub-02" = matrix(1, 5, 5),
    "sub-03" = matrix(1, 5, 5),
    "sub-04" = matrix(1, 5, 5)
  )
  adat <- AlignmentData(data_list)

  folds1 <- create_cv_folds(adat, method = "kfold", k = 2, seed = 42)
  folds2 <- create_cv_folds(adat, method = "kfold", k = 2, seed = 42)

  expect_equal(folds1$assignments, folds2$assignments)
})

test_that("create_cv_folds with stratified works", {
  data_list <- lapply(1:6, function(i) matrix(1, 5, 5))
  names(data_list) <- paste0("sub-0", 1:6)
  adat <- AlignmentData(data_list)

  # Group assignment
  groups <- c("A", "A", "A", "B", "B", "B")

  folds <- create_cv_folds(adat, method = "stratified", k = 2,
    groups = groups, seed = 123
  )

  expect_equal(folds$method, "stratified")
})

test_that("fit_alignment with cv = loso works", {
  ensure_test_aligner("procrustes")

  set.seed(999)
  data_list <- make_test_data_list(n_subjects = 3, n_features = 10, n_obs = 10)
  adat <- AlignmentData(data_list)

  result <- fit_alignment(adat, method = "procrustes", cv = "loso")

  expect_equal(result@cv_info$method, "loso")
  expect_equal(result@cv_info$n_folds, 3)
  expect_false(isTRUE(result@cv_info$anchor_common))
  expect_equal(length(result@aligned), 0)
})

test_that("fit_alignment with cv = kfold works", {
  ensure_test_aligner("procrustes")

  set.seed(998)
  data_list <- make_test_data_list(n_subjects = 4, n_features = 10, n_obs = 10)
  adat <- AlignmentData(data_list)

  result <- fit_alignment(adat, method = "procrustes", cv = "kfold", cv_folds = 2)

  expect_equal(result@cv_info$method, "kfold")
  expect_equal(result@cv_info$n_folds, 2)
})

test_that("is_cv_result correctly identifies CV results", {
  ensure_test_aligner("procrustes")

  set.seed(997)
  data_list <- make_test_data_list(n_subjects = 2, n_features = 10, n_obs = 10)
  adat <- AlignmentData(data_list)

  result_no_cv <- fit_alignment(adat, method = "procrustes", cv = "none")
  result_cv <- fit_alignment(adat, method = "procrustes", cv = "loso")

  expect_false(is_cv_result(result_no_cv))
  expect_true(is_cv_result(result_cv))
})

test_that("get_fold_assignments returns fold info", {
  ensure_test_aligner("procrustes")

  set.seed(996)
  data_list <- make_test_data_list(n_subjects = 3, n_features = 10, n_obs = 10)
  adat <- AlignmentData(data_list)

  result <- fit_alignment(adat, method = "procrustes", cv = "loso")

  assignments <- get_fold_assignments(result)

  expect_equal(length(assignments), 3)
  expect_true(all(names(assignments) %in% adat@subjects))
})

test_that("create_cv_folds accepts numeric input", {
  folds <- create_cv_folds(5, method = "loso")

  expect_equal(folds$n_folds, 5)
  expect_equal(length(folds$assignments), 5)
  expect_true(all(grepl("^sub-", names(folds$assignments))))
})

test_that("create_cv_folds accepts character vector", {
  subjects <- c("alice", "bob", "charlie")
  folds <- create_cv_folds(subjects, method = "loso")

  expect_equal(folds$n_folds, 3)
  expect_equal(names(folds$assignments), subjects)
})

test_that("create_cv_folds rejects invalid input", {
  expect_error(
    create_cv_folds(list(a = 1)),
    "must be AlignmentData"
  )
})

test_that("stratified CV requires groups", {
  expect_error(
    create_cv_folds(5, method = "stratified"),
    "groups.*required"
  )
})

test_that("stratified CV validates groups length", {
  expect_error(
    create_cv_folds(5, method = "stratified", groups = c("A", "B")),
    "same length"
  )
})

test_that("run_cv_alignment works with loso string", {
  ensure_test_aligner("procrustes")

  set.seed(42)
  data_list <- make_test_data_list(n_subjects = 3, n_features = 10, n_obs = 5)
  adat <- AlignmentData(data_list)

  cv_result <- run_cv_alignment(adat, method = "procrustes", cv_folds = "loso")

  expect_true(is.list(cv_result))
  expect_s4_class(cv_result$result, "AlignmentResult")
  expect_null(cv_result$fold_results)
  expect_equal(cv_result$cv_info$method, "loso")
})

test_that("run_cv_alignment works with precomputed folds", {
  ensure_test_aligner("procrustes")

  set.seed(42)
  data_list <- make_test_data_list(n_subjects = 4, n_features = 10, n_obs = 5)
  adat <- AlignmentData(data_list)

  # Precompute folds
  folds <- create_cv_folds(adat, method = "kfold", k = 2, seed = 123)

  cv_result <- run_cv_alignment(adat, method = "procrustes", cv_folds = folds)

  expect_equal(cv_result$cv_info$method, "kfold")
  expect_equal(cv_result$cv_info$n_folds, 2)
})

test_that("run_cv_alignment accepts list input", {
  ensure_test_aligner("procrustes")

  set.seed(42)
  data_list <- make_test_data_list(n_subjects = 2, n_features = 10, n_obs = 5)

  # Pass list directly, not AlignmentData
  cv_result <- run_cv_alignment(data_list, method = "procrustes", cv_folds = "loso")

  expect_s4_class(cv_result$result, "AlignmentResult")
})

test_that("get_fold_assignments returns NULL for non-CV result", {
  ensure_test_aligner("procrustes")

  set.seed(42)
  data_list <- make_test_data_list(n_subjects = 2, n_features = 5, n_obs = 5)
  adat <- AlignmentData(data_list)

  result <- fit_alignment(adat, method = "procrustes", cv = "none")

  assignments <- get_fold_assignments(result)
  expect_null(assignments)
})


# --- Additional tests for uncovered lines ---

test_that("run_cv_alignment errors when obs_labels conflict with AlignmentData", {
  ensure_test_aligner("procrustes")

  set.seed(42)
  n_feat <- 10
  n_obs <- 5
  data_list <- make_test_data_list(n_subjects = 3, n_features = n_feat, n_obs = n_obs)

  # Create AlignmentData with obs_labels already set
  adat <- AlignmentData(data_list, obs_labels = make_test_obs_labels(n_obs))

  # Pass conflicting obs_labels to run_cv_alignment
  expect_error(
    run_cv_alignment(adat, method = "procrustes", cv_folds = "loso",
                     obs_labels = paste0("different_", 1:n_obs)),
    "obs_labels supplied but AlignmentData already has different obs_labels"
  )
})

test_that("run_cv_alignment sets obs_labels when AlignmentData has none", {
  ensure_test_aligner("procrustes")

  set.seed(42)
  n_feat <- 10
  n_obs <- 5
  data_list <- make_test_data_list(n_subjects = 3, n_features = n_feat, n_obs = n_obs)

  # Create AlignmentData without obs_labels
  adat <- AlignmentData(data_list)
  expect_null(adat@obs_labels)

  new_labels <- make_test_obs_labels(n_obs)

  # Should not error; obs_labels get set on the internal copy
  cv_result <- run_cv_alignment(adat, method = "procrustes", cv_folds = "loso",
                                obs_labels = new_labels)

  expect_s4_class(cv_result$result, "AlignmentResult")
  expect_null(cv_result$fold_results)
})

test_that("run_cv_alignment with data-driven reference yields fold-specific anchor", {
  ensure_test_aligner("procrustes")

  set.seed(42)
  n_feat <- 10
  n_obs <- 5
  data_list <- make_test_data_list(n_subjects = 4, n_features = n_feat, n_obs = n_obs)
  adat <- AlignmentData(data_list)

  # Default reference = "medoid" is data-driven
  cv_result <- run_cv_alignment(adat, method = "procrustes", cv_folds = "loso",
                                reference = "medoid")

  expect_s4_class(cv_result$result, "AlignmentResult")
  expect_equal(cv_result$result@cv_info$reference_kind, "data_driven")
  expect_false(cv_result$result@cv_info$anchor_common)

  # Model reference should be "fold_specific" for data-driven references
  expect_equal(cv_result$result@model@reference, "fold_specific")
})

test_that("run_cv_alignment with template matrix reference yields common anchor", {
  ensure_test_aligner("procrustes")

  set.seed(42)
  n_feat <- 10
  n_obs <- 5
  data_list <- make_test_data_list(n_subjects = 4, n_features = n_feat, n_obs = n_obs)
  adat <- AlignmentData(data_list)

  # Use a template matrix as reference (external anchor)
  template <- make_test_matrix(n_features = n_feat, n_obs = n_obs)
  cv_result <- run_cv_alignment(adat, method = "procrustes", cv_folds = "loso",
                                reference = template)

  expect_s4_class(cv_result$result, "AlignmentResult")
  expect_equal(cv_result$result@cv_info$reference_kind, "template")
  expect_true(cv_result$result@cv_info$anchor_common)
  expect_null(cv_result$result@cv_info$anchor_note)
})

test_that("has_common_anchor works for AlignmentModel", {
  # Model with fold_specific reference -> FALSE
  model_fold <- AlignmentModel(
    transforms = list("sub-01" = diag(3)),
    reference = "fold_specific",
    method = "test"
  )
  expect_false(has_common_anchor(model_fold))

  # Model with a specific subject reference -> TRUE
  model_subj <- AlignmentModel(
    transforms = list("sub-01" = diag(3)),
    reference = "sub-01",
    method = "test"
  )
  expect_true(has_common_anchor(model_subj))

  # Model with consensus reference -> TRUE
  model_consensus <- AlignmentModel(
    transforms = list("sub-01" = diag(3)),
    reference = "consensus",
    method = "test"
  )
  expect_true(has_common_anchor(model_consensus))
})

test_that("has_common_anchor errors on bad input", {
  expect_error(
    has_common_anchor("not_a_model"),
    "must be an AlignmentResult or AlignmentModel"
  )

  expect_error(
    has_common_anchor(42),
    "must be an AlignmentResult or AlignmentModel"
  )

  expect_error(
    has_common_anchor(list(a = 1)),
    "must be an AlignmentResult or AlignmentModel"
  )
})

test_that("validate_common_anchor with error action throws error", {
  model_fold <- AlignmentModel(
    transforms = list("sub-01" = diag(3)),
    reference = "fold_specific",
    method = "test"
  )

  expect_error(
    validate_common_anchor(model_fold, action = "error"),
    "No common anchor detected"
  )
})

test_that("validate_common_anchor with warn action gives warning", {
  model_fold <- AlignmentModel(
    transforms = list("sub-01" = diag(3)),
    reference = "fold_specific",
    method = "test"
  )

  expect_warning(
    result <- validate_common_anchor(model_fold, action = "warn"),
    "No common anchor detected"
  )
  expect_false(result)
})

test_that("validate_common_anchor with silent action returns FALSE invisibly", {
  model_fold <- AlignmentModel(
    transforms = list("sub-01" = diag(3)),
    reference = "fold_specific",
    method = "test"
  )

  # No error, no warning
  expect_no_warning(
    expect_no_error(
      result <- validate_common_anchor(model_fold, action = "silent")
    )
  )
  expect_false(result)
})

test_that("validate_common_anchor returns TRUE for common-anchor result", {
  model_common <- AlignmentModel(
    transforms = list("sub-01" = diag(3)),
    reference = "sub-01",
    method = "test"
  )

  result_common <- AlignmentResult(
    model = model_common,
    aligned = list("sub-01" = diag(3)),
    cv_info = list(method = "loso", anchor_common = TRUE)
  )

  result <- validate_common_anchor(result_common, action = "error")
  expect_true(result)
})

test_that("validate_common_anchor includes context in messages", {
  model_fold <- AlignmentModel(
    transforms = list("sub-01" = diag(3)),
    reference = "fold_specific",
    method = "test"
  )

  expect_error(
    validate_common_anchor(model_fold, action = "error", context = "ISC computation"),
    "ISC computation"
  )

  expect_warning(
    validate_common_anchor(model_fold, action = "warn", context = "decoding"),
    "decoding"
  )
})


# ---------- More CV coverage tests ----------

test_that("has_common_anchor for AlignmentResult without anchor_common field", {
  # Result with cv_info but no anchor_common field -> fallback to reference check
  model <- AlignmentModel(
    transforms = list("sub-01" = diag(3)),
    reference = "fold_specific",
    method = "test"
  )
  result <- AlignmentResult(
    model = model,
    aligned = list(),
    cv_info = list(method = "loso")  # no anchor_common
  )
  expect_false(has_common_anchor(result))
})

test_that("has_common_anchor for AlignmentResult with no CV info", {
  model <- AlignmentModel(
    transforms = list("sub-01" = diag(3)),
    reference = "sub-01",
    method = "test"
  )
  result <- AlignmentResult(
    model = model,
    aligned = list("sub-01" = diag(3)),
    cv_info = list()  # empty cv_info
  )
  expect_true(has_common_anchor(result))
})

test_that("run_cv_alignment with fixed subject reference yields common anchor", {
  ensure_test_aligner("procrustes")

  set.seed(42)
  data_list <- make_test_data_list(n_subjects = 3, n_features = 10, n_obs = 5)
  adat <- AlignmentData(data_list)

  expect_warning(
    cv_result <- run_cv_alignment(
      adat, method = "procrustes", cv_folds = "loso",
      reference = "sub-01"
    ),
    "Reference subject 'sub-01' appears in a test fold"
  )

  expect_equal(cv_result$result@cv_info$reference_kind, "fixed_subject")
  expect_true(cv_result$result@cv_info$anchor_common)

  # The fixed reference subject is not evaluated as held-out, but it is
  # represented in the aligned output with identity transform.
  expect_true("sub-01" %in% names(cv_result$result@model@transforms))
  expect_equal(cv_result$result@model@transforms[["sub-01"]], diag(10))
  expect_true("sub-01" %in% names(cv_result$result@aligned))

  expect_true("sub-01" %in% names(cv_result$result@cv_info$folds))
  expect_length(cv_result$result@cv_info$folds[["sub-01"]]$test, 0)
})

test_that("run_cv_alignment with kfold string works", {
  ensure_test_aligner("procrustes")

  set.seed(42)
  data_list <- make_test_data_list(n_subjects = 6, n_features = 10, n_obs = 5)
  adat <- AlignmentData(data_list)

  cv_result <- run_cv_alignment(
    adat, method = "procrustes", cv_folds = "kfold", k = 2
  )

  expect_s4_class(cv_result$result, "AlignmentResult")
  expect_equal(cv_result$cv_info$method, "kfold")
})

test_that("run_cv_alignment sets matching obs_labels on AlignmentData", {
  ensure_test_aligner("procrustes")

  set.seed(42)
  data_list <- make_test_data_list(n_subjects = 3, n_features = 10, n_obs = 3)
  labs <- make_test_obs_labels(3)
  adat <- AlignmentData(data_list, obs_labels = labs)

  # Same labels -> should not error
  cv_result <- run_cv_alignment(
    adat, method = "procrustes", cv_folds = "loso",
    obs_labels = labs
  )
  expect_s4_class(cv_result$result, "AlignmentResult")
})
