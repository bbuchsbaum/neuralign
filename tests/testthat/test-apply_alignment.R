test_that("apply_alignment works with new subjects", {
  neuralign:::.register_procrustes()

  set.seed(333)
  # Training data
  train_data <- list(
    "sub-01" = matrix(rnorm(100), 10, 10),
    "sub-02" = matrix(rnorm(100), 10, 10),
    "sub-03" = matrix(rnorm(100), 10, 10)
  )
  train_adat <- AlignmentData(train_data)

  # Fit model
  result <- fit_alignment(train_adat, method = "procrustes")

  # New subject data
  new_data <- list(
    "sub-04" = matrix(rnorm(100), 10, 10)
  )
  new_adat <- AlignmentData(new_data)

  # Apply to new subject
  new_result <- apply_alignment(result, new_adat, warn_leakage = FALSE)

  expect_s4_class(new_result, "AlignmentResult")
  expect_true("sub-04" %in% names(new_result@aligned))
})

test_that("apply_alignment uses existing transforms for known subjects", {
  neuralign:::.register_procrustes()

  set.seed(444)
  data_list <- list(
    "sub-01" = matrix(rnorm(100), 10, 10),
    "sub-02" = matrix(rnorm(100), 10, 10)
  )
  adat <- AlignmentData(data_list)

  result <- fit_alignment(adat, method = "procrustes")
  model <- get_model(result)

  # Apply to same subjects (should use existing transforms)
  applied <- suppressWarnings(apply_alignment(model, adat, warn_leakage = FALSE))

  # Check that transforms are the same
  expect_equal(
    applied@model@transforms[["sub-01"]],
    model@transforms[["sub-01"]]
  )
})

test_that("apply_alignment warns about leakage", {
  neuralign:::.register_procrustes()

  set.seed(555)
  data_list <- list(
    "sub-01" = matrix(rnorm(100), 10, 10),
    "sub-02" = matrix(rnorm(100), 10, 10)
  )
  adat <- AlignmentData(data_list)

  result <- fit_alignment(adat, method = "procrustes")

  # Should warn when applying to training subjects
  expect_warning(
    apply_alignment(result, adat, warn_leakage = TRUE),
    "leakage"
  )
})

test_that("apply_alignment respects fit_new = FALSE", {
  neuralign:::.register_procrustes()

  set.seed(666)
  train_data <- list(
    "sub-01" = matrix(rnorm(100), 10, 10),
    "sub-02" = matrix(rnorm(100), 10, 10)
  )
  train_adat <- AlignmentData(train_data)

  result <- fit_alignment(train_adat, method = "procrustes")

  new_data <- list(
    "sub-03" = matrix(rnorm(100), 10, 10)
  )
  new_adat <- AlignmentData(new_data)

  # With fit_new = FALSE, should warn about missing transforms
  expect_warning(
    apply_alignment(result, new_adat, fit_new = FALSE, warn_leakage = FALSE),
    "No transforms"
  )
})

test_that("apply_transform validates dimensions", {
  transform <- matrix(1, 5, 10)  # 5x10
  data <- matrix(1, 10, 20)      # 10x20

  # Should work: transform cols (10) matches data rows (10)
  result <- apply_transform(transform, data)
  expect_equal(dim(result), c(5, 20))

  # Should fail: dimension mismatch
  bad_data <- matrix(1, 8, 20)  # 8x20, doesn't match 10
  expect_error(
    apply_transform(transform, bad_data),
    "mismatch"
  )
})

test_that("inverse_transform works for orthogonal transforms", {
  neuralign:::.register_procrustes()

  set.seed(777)
  data_list <- list(
    "sub-01" = matrix(rnorm(100), 10, 10),
    "sub-02" = matrix(rnorm(100), 10, 10)
  )
  adat <- AlignmentData(data_list)

  result <- fit_alignment(adat, method = "procrustes")
  model <- get_model(result)

  # Get forward and inverse transforms
  forward <- get_transform(model, "sub-01")
  inverse <- inverse_transform(model, "sub-01")

  # For orthogonal: forward %*% inverse should be identity
  identity_approx <- forward %*% inverse
  expect_equal(identity_approx, diag(nrow(forward)), tolerance = 1e-10)
})
