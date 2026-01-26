test_that("has_common_anchor detects fold-specific anchors under CV", {
  neuralign:::.register_procrustes()

  set.seed(123)
  data_list <- list(
    "sub-01" = matrix(rnorm(100), 10, 10),
    "sub-02" = matrix(rnorm(100), 10, 10),
    "sub-03" = matrix(rnorm(100), 10, 10)
  )
  adat <- AlignmentData(data_list)

  res_cv <- fit_alignment(adat, method = "procrustes", cv = "loso")
  expect_false(has_common_anchor(res_cv))

  res_full <- fit_alignment(adat, method = "procrustes", cv = "none")
  expect_true(has_common_anchor(res_full))
})

test_that("template reference yields a common anchor under CV", {
  neuralign:::.register_procrustes()

  set.seed(124)
  data_list <- list(
    "sub-01" = matrix(rnorm(100), 10, 10),
    "sub-02" = matrix(rnorm(100), 10, 10),
    "sub-03" = matrix(rnorm(100), 10, 10)
  )
  adat <- AlignmentData(data_list)

  template <- matrix(rnorm(100), 10, 10)
  res_cv <- fit_alignment(adat, method = "procrustes", cv = "kfold", cv_folds = 3, reference = template)

  expect_true(has_common_anchor(res_cv))
  expect_true(isTRUE(res_cv@cv_info$anchor_common))
})

test_that("validate_common_anchor can error on fold-specific CV results", {
  neuralign:::.register_procrustes()

  set.seed(125)
  data_list <- list(
    "sub-01" = matrix(rnorm(100), 10, 10),
    "sub-02" = matrix(rnorm(100), 10, 10),
    "sub-03" = matrix(rnorm(100), 10, 10)
  )
  adat <- AlignmentData(data_list)

  res_cv <- fit_alignment(adat, method = "procrustes", cv = "loso")
  expect_error(
    validate_common_anchor(res_cv, action = "error"),
    "No common anchor"
  )
})

test_that("fit_alignment accepts explicit fold specs", {
  neuralign:::.register_procrustes()

  set.seed(127)
  data_list <- list(
    "sub-01" = matrix(rnorm(100), 10, 10),
    "sub-02" = matrix(rnorm(100), 10, 10),
    "sub-03" = matrix(rnorm(100), 10, 10),
    "sub-04" = matrix(rnorm(100), 10, 10)
  )
  adat <- AlignmentData(data_list)

  folds <- create_cv_folds(adat, method = "kfold", k = 2, seed = 1)
  res <- fit_alignment(adat, method = "procrustes", cv = "kfold", cv_folds = folds)

  expect_equal(res@cv_info$method, "kfold")
  expect_equal(res@cv_info$n_folds, 2)
  expect_equal(res@cv_info$fold_assignments, folds$assignments)
  expect_true(is.list(res@cv_info$folds))
})

test_that("apply_alignment refuses to fit new subjects for fold-specific anchor models", {
  neuralign:::.register_procrustes()

  set.seed(126)
  data_list <- list(
    "sub-01" = matrix(rnorm(100), 10, 10),
    "sub-02" = matrix(rnorm(100), 10, 10),
    "sub-03" = matrix(rnorm(100), 10, 10)
  )
  adat <- AlignmentData(data_list)

  res_cv <- fit_alignment(adat, method = "procrustes", cv = "loso")

  new_subject <- AlignmentData(list("sub-04" = matrix(rnorm(100), 10, 10)))
  expect_error(
    apply_alignment(res_cv, new_subject),
    "fold-specific anchors"
  )
})
