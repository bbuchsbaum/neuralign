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

test_that("apply_alignment validates model input", {
  expect_error(
    apply_alignment("not a model", list("sub-01" = diag(5))),
    "must be an AlignmentModel"
  )
})

test_that("apply_alignment coerces list input to AlignmentData", {
  neuralign:::.register_procrustes()

  set.seed(888)
  train_data <- list(
    "sub-01" = matrix(rnorm(50), 10, 5),
    "sub-02" = matrix(rnorm(50), 10, 5)
  )
  train_adat <- AlignmentData(train_data)

  result <- fit_alignment(train_adat, method = "procrustes")

  # Apply with raw list (should be coerced)
  new_list <- list("sub-03" = matrix(rnorm(50), 10, 5))
  new_result <- apply_alignment(result, new_list, warn_leakage = FALSE)

  expect_s4_class(new_result, "AlignmentResult")
  expect_true("sub-03" %in% names(get_aligned(new_result)))
})

test_that("apply_alignment errors for unregistered method", {
  # Create a model with unregistered method
  transforms <- list("sub-01" = diag(5))
  model <- AlignmentModel(
    transforms = transforms,
    reference = "consensus",
    method = "nonexistent_method_xyz",
    train_subjects = "sub-01"
  )

  new_data <- list("sub-02" = matrix(1, 5, 3))

  expect_error(
    apply_alignment(model, new_data, warn_leakage = FALSE),
    "not registered"
  )
})

test_that("apply_alignment handles method that doesn't support new subjects", {
  # Register a method that doesn't support new subjects
  register_aligner(
    name = "no_new_subj",
    fit_fn = function(data, reference, ...) {
      list(
        transforms = list("x" = diag(5)),
        reference_data = NULL
      )
    },
    capabilities = list(
      supports_new_subject = FALSE,
      returns = "operator"
    ),
    package = "neuralign"
  )

  transforms <- list("sub-01" = diag(5))
  model <- AlignmentModel(
    transforms = transforms,
    reference = "consensus",
    method = "no_new_subj",
    train_subjects = "sub-01"
  )

  new_data <- list("sub-02" = matrix(1, 5, 3))

  expect_error(
    apply_alignment(model, new_data, warn_leakage = FALSE),
    "does not support fitting transforms"
  )

  unregister_aligner("no_new_subj")
})

test_that("apply_transform handles non-matrix input", {
  transform <- diag(5)
  data_vec <- 1:5  # Vector, not matrix

  # Should be coerced to matrix and work
  result <- apply_transform(transform, data_vec)
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(5, 1))
})

test_that("apply_transform rejects invalid transform", {
  expect_error(
    apply_transform("not a matrix", diag(5)),
    "must be a matrix"
  )
})

test_that("inverse_transform with explicit method = transpose", {
  neuralign:::.register_procrustes()

  set.seed(999)
  data_list <- list(
    "sub-01" = matrix(rnorm(50), 10, 5),
    "sub-02" = matrix(rnorm(50), 10, 5)
  )
  adat <- AlignmentData(data_list)

  result <- fit_alignment(adat, method = "procrustes")
  model <- get_model(result)

  inverse <- inverse_transform(model, "sub-01", method = "transpose")
  forward <- get_transform(model, "sub-01")

  expect_equal(inverse, t(forward))
})

test_that("inverse_transform with method = solve", {
  transforms <- list("sub-01" = matrix(c(1,2,0,1), 2, 2))  # Invertible but not orthogonal
  model <- AlignmentModel(
    transforms = transforms,
    reference = "consensus",
    method = "test"
  )

  inverse <- inverse_transform(model, "sub-01", method = "solve")
  forward <- get_transform(model, "sub-01")

  # forward %*% inverse should be identity
  expect_equal(forward %*% inverse, diag(2), tolerance = 1e-10)
})

test_that("inverse_transform fails for singular matrix", {
  transforms <- list("sub-01" = matrix(0, 2, 2))  # Singular
  model <- AlignmentModel(
    transforms = transforms,
    reference = "consensus",
    method = "test"
  )

  expect_error(
    inverse_transform(model, "sub-01", method = "solve"),
    "not invertible"
  )
})

test_that("inverse_transform rejects unknown method", {
  transforms <- list("sub-01" = diag(5))
  model <- AlignmentModel(
    transforms = transforms,
    reference = "consensus",
    method = "test"
  )

  expect_error(
    inverse_transform(model, "sub-01", method = "unknown"),
    "Unknown inverse method"
  )
})

test_that("inverse_transform with method = pinv works for non-square operators", {
  transforms <- list("sub-01" = matrix(c(1, 0, 0, 1, 1, 1), nrow = 2)) # 2 x 3
  model <- AlignmentModel(
    transforms = transforms,
    reference = "consensus",
    method = "test"
  )

  forward <- get_transform(model, "sub-01")
  inverse <- inverse_transform(model, "sub-01", method = "pinv")

  expect_equal(dim(inverse), c(ncol(forward), nrow(forward)))
  expect_equal(forward %*% inverse %*% forward, forward, tolerance = 1e-10)
})

test_that("inverse_transform with method = ridge works for non-square operators", {
  transforms <- list("sub-01" = matrix(c(1, 0, 0, 1, 1, 1), nrow = 2)) # 2 x 3
  model <- AlignmentModel(
    transforms = transforms,
    reference = "consensus",
    method = "test"
  )

  forward <- get_transform(model, "sub-01")
  inverse <- inverse_transform(model, "sub-01", method = "ridge", lambda = 1e-3)

  expect_equal(dim(inverse), c(ncol(forward), nrow(forward)))
})

test_that("inverse_transform auto errors for non-square operator", {
  transforms <- list("sub-01" = matrix(c(1, 0, 0, 1, 1, 1), nrow = 2)) # 2 x 3
  model <- AlignmentModel(
    transforms = transforms,
    reference = "consensus",
    method = "test"
  )

  expect_error(
    inverse_transform(model, "sub-01", method = "auto"),
    "not square"
  )
})

test_that("inverse_transform rejects OT transforms", {
  register_aligner(
    name = "ot_dummy",
    fit_fn = function(data, reference, ...) {
      list(transforms = list("x" = diag(2)), reference_data = NULL)
    },
    capabilities = list(transform_type = "ot", returns_invertible = FALSE),
    package = "neuralign"
  )

  model <- AlignmentModel(
    transforms = list("sub-01" = diag(2)),
    reference = "consensus",
    method = "ot_dummy"
  )

  expect_error(
    inverse_transform(model, "sub-01", method = "pinv"),
    "OT-style"
  )

  unregister_aligner("ot_dummy")
})
