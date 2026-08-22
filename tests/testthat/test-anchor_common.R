test_that("has_common_anchor detects fold-specific anchors under CV", {
  ensure_test_aligner("procrustes")

  set.seed(123)
  data_list <- make_test_data_list(n_subjects = 3, n_features = 10, n_obs = 10)
  adat <- AlignmentData(data_list)

  res_cv <- fit_alignment(adat, method = "procrustes", cv = "loso")
  expect_false(has_common_anchor(res_cv))

  res_full <- fit_alignment(adat, method = "procrustes", cv = "none")
  expect_true(has_common_anchor(res_full))
})

test_that("template reference yields a common anchor under CV", {
  ensure_test_aligner("procrustes")

  set.seed(124)
  data_list <- make_test_data_list(n_subjects = 3, n_features = 10, n_obs = 10)
  adat <- AlignmentData(data_list)

  template <- make_test_matrix(n_features = 10, n_obs = 10)
  res_cv <- fit_alignment(adat, method = "procrustes", cv = "kfold", cv_folds = 3, reference = template)

  expect_true(has_common_anchor(res_cv))
  expect_true(isTRUE(res_cv@cv_info$anchor_common))
})

test_that("validate_common_anchor can error on fold-specific CV results", {
  ensure_test_aligner("procrustes")

  set.seed(125)
  data_list <- make_test_data_list(n_subjects = 3, n_features = 10, n_obs = 10)
  adat <- AlignmentData(data_list)

  res_cv <- fit_alignment(adat, method = "procrustes", cv = "loso")
  expect_error(
    validate_common_anchor(res_cv, action = "error"),
    "No common anchor"
  )
})

test_that("fit_alignment accepts explicit fold specs", {
  ensure_test_aligner("procrustes")

  set.seed(127)
  data_list <- make_test_data_list(n_subjects = 4, n_features = 10, n_obs = 10)
  adat <- AlignmentData(data_list)

  folds <- create_cv_folds(adat, method = "kfold", k = 2, seed = 1)
  res <- fit_alignment(adat, method = "procrustes", cv = "kfold", cv_folds = folds)

  expect_equal(res@cv_info$method, "kfold")
  expect_equal(res@cv_info$n_folds, 2)
  expect_equal(res@cv_info$fold_assignments, folds$assignments)
  expect_true(is.list(res@cv_info$folds))
})

test_that("apply_alignment refuses to fit new subjects for fold-specific anchor models", {
  ensure_test_aligner("procrustes")

  set.seed(126)
  data_list <- make_test_data_list(n_subjects = 3, n_features = 10, n_obs = 10)
  adat <- AlignmentData(data_list)

  res_cv <- fit_alignment(adat, method = "procrustes", cv = "loso")

  new_subject <- make_test_alignment_data(n_subjects = 1, n_features = 10, n_obs = 10, subject_ids = "sub-04")
  expect_error(
    apply_alignment(res_cv, new_subject),
    "not available for fold-specific models"
  )
})
