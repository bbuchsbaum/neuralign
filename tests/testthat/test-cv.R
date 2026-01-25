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
  neuralign:::.register_procrustes()

  set.seed(42)
  data_list <- list(
    "sub-01" = matrix(rnorm(50), 10, 5),
    "sub-02" = matrix(rnorm(50), 10, 5),
    "sub-03" = matrix(rnorm(50), 10, 5)
  )
  adat <- AlignmentData(data_list)

  cv_result <- run_cv_alignment(adat, method = "procrustes", cv_folds = "loso")

  expect_true(is.list(cv_result))
  expect_s4_class(cv_result$result, "AlignmentResult")
  expect_equal(length(cv_result$fold_results), 3)
  expect_equal(cv_result$cv_info$method, "loso")
})

test_that("run_cv_alignment works with precomputed folds", {
  neuralign:::.register_procrustes()

  set.seed(42)
  data_list <- list(
    "sub-01" = matrix(rnorm(50), 10, 5),
    "sub-02" = matrix(rnorm(50), 10, 5),
    "sub-03" = matrix(rnorm(50), 10, 5),
    "sub-04" = matrix(rnorm(50), 10, 5)
  )
  adat <- AlignmentData(data_list)

  # Precompute folds
  folds <- create_cv_folds(adat, method = "kfold", k = 2, seed = 123)

  cv_result <- run_cv_alignment(adat, method = "procrustes", cv_folds = folds)

  expect_equal(cv_result$cv_info$method, "kfold")
  expect_equal(cv_result$cv_info$n_folds, 2)
})

test_that("run_cv_alignment accepts list input", {
  neuralign:::.register_procrustes()

  set.seed(42)
  data_list <- list(
    "sub-01" = matrix(rnorm(50), 10, 5),
    "sub-02" = matrix(rnorm(50), 10, 5)
  )

  # Pass list directly, not AlignmentData
  cv_result <- run_cv_alignment(data_list, method = "procrustes", cv_folds = "loso")

  expect_s4_class(cv_result$result, "AlignmentResult")
})

test_that("get_fold_assignments returns NULL for non-CV result", {
  neuralign:::.register_procrustes()

  set.seed(42)
  data_list <- list(
    "sub-01" = matrix(rnorm(25), 5, 5),
    "sub-02" = matrix(rnorm(25), 5, 5)
  )
  adat <- AlignmentData(data_list)

  result <- fit_alignment(adat, method = "procrustes", cv = "none")

  assignments <- get_fold_assignments(result)
  expect_null(assignments)
})
