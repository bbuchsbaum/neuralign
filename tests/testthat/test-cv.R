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
  neuralign:::.register_procrustes()

  set.seed(999)
  data_list <- list(
    "sub-01" = matrix(rnorm(100), 10, 10),
    "sub-02" = matrix(rnorm(100), 10, 10),
    "sub-03" = matrix(rnorm(100), 10, 10)
  )
  adat <- AlignmentData(data_list)

  result <- fit_alignment(adat, method = "procrustes", cv = "loso")

  expect_equal(result@cv_info$method, "loso")
  expect_equal(result@cv_info$n_folds, 3)
  expect_equal(length(result@aligned), 3)
})

test_that("fit_alignment with cv = kfold works", {
  neuralign:::.register_procrustes()

  set.seed(998)
  data_list <- list(
    "sub-01" = matrix(rnorm(100), 10, 10),
    "sub-02" = matrix(rnorm(100), 10, 10),
    "sub-03" = matrix(rnorm(100), 10, 10),
    "sub-04" = matrix(rnorm(100), 10, 10)
  )
  adat <- AlignmentData(data_list)

  result <- fit_alignment(adat, method = "procrustes", cv = "kfold", cv_folds = 2)

  expect_equal(result@cv_info$method, "kfold")
  expect_equal(result@cv_info$n_folds, 2)
})

test_that("is_cv_result correctly identifies CV results", {
  neuralign:::.register_procrustes()

  set.seed(997)
  data_list <- list(
    "sub-01" = matrix(rnorm(100), 10, 10),
    "sub-02" = matrix(rnorm(100), 10, 10)
  )
  adat <- AlignmentData(data_list)

  result_no_cv <- fit_alignment(adat, method = "procrustes", cv = "none")
  result_cv <- fit_alignment(adat, method = "procrustes", cv = "loso")

  expect_false(is_cv_result(result_no_cv))
  expect_true(is_cv_result(result_cv))
})

test_that("get_fold_assignments returns fold info", {
  neuralign:::.register_procrustes()

  set.seed(996)
  data_list <- list(
    "sub-01" = matrix(rnorm(100), 10, 10),
    "sub-02" = matrix(rnorm(100), 10, 10),
    "sub-03" = matrix(rnorm(100), 10, 10)
  )
  adat <- AlignmentData(data_list)

  result <- fit_alignment(adat, method = "procrustes", cv = "loso")

  assignments <- get_fold_assignments(result)

  expect_equal(length(assignments), 3)
  expect_true(all(names(assignments) %in% adat@subjects))
})
