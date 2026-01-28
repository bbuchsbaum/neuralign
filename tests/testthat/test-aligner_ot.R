# Tests for OT aligner functions
# Note: Full GW/FPGW tests require manifoldalign package

test_that(".coupling_to_operator converts coupling correctly", {
  # Simple coupling matrix (source x target)
  P <- matrix(c(
    0.4, 0.1,
    0.1, 0.4
  ), 2, 2, byrow = TRUE)

  op <- neuralign:::.coupling_to_operator(P)

  # Should be transposed and row-normalized
  expect_equal(dim(op), c(2, 2))

  # Each row should sum to 1 (stochastic)
  expect_equal(rowSums(op), c(1, 1), tolerance = 1e-10)
})

test_that(".coupling_to_operator handles zero rows", {
  # Coupling with a zero column (becomes zero row after transpose)
  P <- matrix(c(
    0.5, 0,
    0.5, 0
  ), 2, 2, byrow = TRUE)

  op <- neuralign:::.coupling_to_operator(P)

  # Should not have NaN (division by zero handled)
  expect_false(any(is.nan(op)))
})

test_that(".gw_capabilities has correct structure", {
  caps <- neuralign:::.gw_capabilities

  expect_true(caps$supports_cv)
  expect_equal(caps$transform_type, "ot")
  expect_false(caps$mass_preserving)
  expect_false(caps$returns_invertible)
  expect_equal(caps$returns, "operator")
})

test_that(".fpgw_capabilities has correct structure", {
  caps <- neuralign:::.fpgw_capabilities

  expect_true(caps$supports_cv)
  expect_equal(caps$transform_type, "ot")
  expect_false(caps$mass_preserving)
})

test_that("GW fit function requires manifoldalign", {
  skip_if(requireNamespace("manifoldalign", quietly = TRUE),
          "manifoldalign is available")

  data_list <- make_test_data_list(n_subjects = 2, n_features = 10, n_obs = 2)
  adat <- AlignmentData(data_list)

  expect_error(
    neuralign:::.gw_fit(adat),
    "manifoldalign"
  )
})

test_that("FPGW fit function requires manifoldalign", {
  skip_if(requireNamespace("manifoldalign", quietly = TRUE),
          "manifoldalign is available")

  data_list <- make_test_data_list(n_subjects = 2, n_features = 10, n_obs = 2)
  adat <- AlignmentData(data_list)

  expect_error(
    neuralign:::.fpgw_fit(adat),
    "manifoldalign"
  )
})

test_that("GW consensus reference uses arithmetic mean", {
  skip_if_not_installed("manifoldalign")

  data_list <- list(
    s1 = matrix(1:4, 2, 2),
    s2 = matrix(5:8, 2, 2)
  )
  adat <- AlignmentData(data_list)

  result <- neuralign:::.gw_fit(
    adat,
    reference = "consensus",
    epsilon = 0.1,
    max_iter = 2
  )

  expect_equal(result$reference_data, (data_list$s1 + data_list$s2) / 2)
})

test_that("GW barycenter reference errors (not implemented)", {
  skip_if_not_installed("manifoldalign")

  data_list <- make_test_data_list(n_subjects = 2, n_features = 6, n_obs = 3)
  adat <- AlignmentData(data_list)

  expect_error(
    neuralign:::.gw_fit(adat, reference = "barycenter", epsilon = 0.1, max_iter = 2),
    "barycenter reference is not implemented"
  )
})

test_that("GW fit errors on non-finite inputs", {
  skip_if_not_installed("manifoldalign")

  X1 <- matrix(rnorm(12), 4, 3)
  X1[1, 1] <- NA_real_
  X2 <- matrix(rnorm(12), 4, 3)
  adat <- AlignmentData(list(s1 = X1, s2 = X2))

  expect_error(neuralign:::.gw_fit(adat, reference = "consensus", epsilon = 0.1, max_iter = 1))
})

test_that(".extract_pair_plan indexes packed pairs correctly", {
  # For n=4, packed order indices are:
  # 1:(1,2), 2:(1,3), 3:(1,4), 4:(2,3), 5:(2,4), 6:(3,4)
  plans <- as.list(seq_len(6))
  expect_equal(neuralign:::.extract_pair_plan(plans, 1, 2, 4), 1L)
  expect_equal(neuralign:::.extract_pair_plan(plans, 1, 4, 4), 3L)
  expect_equal(neuralign:::.extract_pair_plan(plans, 2, 4, 4), 5L)
  # Order-invariant
  expect_equal(neuralign:::.extract_pair_plan(plans, 4, 2, 4), 5L)
})

# Integration tests if manifoldalign is available
test_that("GW aligner works with manifoldalign", {
  skip_if_not_installed("manifoldalign")

  set.seed(42)
  data_list <- make_test_data_list(n_subjects = 3, n_features = 10, n_obs = 5)
  adat <- AlignmentData(data_list)

  # Test fit
  result <- neuralign:::.gw_fit(
    adat,
    reference = "medoid",
    epsilon = 0.1,
    max_iter = 10
  )

  expect_true(is.list(result))
  expect_true("transforms" %in% names(result))
  expect_equal(length(result$transforms), 3)

  # Each transform should be a valid operator
  for (t in result$transforms) {
    expect_true(is.matrix(t))
    expect_equal(nrow(t), 10)
  }
})

test_that("GW registration via fit_alignment works", {
  skip_if_not_installed("manifoldalign")

  set.seed(42)
  data_list <- make_test_data_list(n_subjects = 2, n_features = 10, n_obs = 3)
  adat <- AlignmentData(data_list)

  result <- fit_alignment(
    adat,
    method = "gw",
    reference = "medoid",
    epsilon = 0.1,
    max_iter = 10
  )

  expect_s4_class(result, "AlignmentResult")
  expect_equal(length(get_aligned(result)), 2)
})

test_that("apply_alignment fits GW transforms for new subjects using stored args", {
  skip_if_not_installed("manifoldalign")
  ensure_test_aligner("gw")

  set.seed(101)
  train_list <- make_test_data_list(
    n_subjects = 2,
    n_features = 8,
    n_obs = 4,
    subject_ids = c("s1", "s2")
  )
  train_data <- AlignmentData(train_list)

  res <- fit_alignment(
    train_data,
    method = "gw",
    reference = "s1",
    epsilon = 0.2,
    max_iter = 5,
    tol = 1e-6,
    compute_quality = FALSE
  )
  model <- get_model(res)

  new_data <- AlignmentData(list(s3 = make_test_matrix(8, 4)))
  applied <- apply_alignment(res, new_data)
  A_apply <- get_transform(get_model(applied), "s3")

  expected <- neuralign:::.gw_fit(
    new_data,
    reference = get_reference(model),
    epsilon = 0.2,
    max_iter = 5,
    tol = 1e-6
  )
  A_expected <- expected$transforms[["s3"]]

  expect_equal(A_apply, A_expected, tolerance = 1e-10)
})

test_that("apply_alignment fits FPGW transforms for new subjects using stored args", {
  skip_if_not_installed("manifoldalign")
  ensure_test_aligner("fpgw")

  set.seed(102)
  train_list <- make_test_data_list(
    n_subjects = 2,
    n_features = 8,
    n_obs = 4,
    subject_ids = c("s1", "s2")
  )
  train_data <- AlignmentData(train_list)

  res <- fit_alignment(
    train_data,
    method = "fpgw",
    reference = "s1",
    omega1 = 0.05,
    epsilon = 0.2,
    max_iter = 5,
    tol = 1e-6,
    compute_quality = FALSE
  )
  model <- get_model(res)

  new_data <- AlignmentData(list(s3 = make_test_matrix(8, 4)))
  applied <- apply_alignment(res, new_data)
  A_apply <- get_transform(get_model(applied), "s3")

  expected <- neuralign:::.fpgw_fit(
    new_data,
    reference = get_reference(model),
    omega1 = 0.05,
    epsilon = 0.2,
    max_iter = 5,
    tol = 1e-6
  )
  A_expected <- expected$transforms[["s3"]]

  expect_equal(A_apply, A_expected, tolerance = 1e-10)
})

test_that(".register_gw registers aligner correctly", {
  # Unregister first if exists
  tryCatch(unregister_aligner("gw"), error = function(e) NULL)

  # Register - expected to succeed
  ensure_test_aligner("gw")

  # Verify registration
  expect_true("gw" %in% available_aligners())

  # Check capabilities via registry
  caps <- aligner_capabilities("gw")
  expect_equal(caps$transform_type, "ot")
  expect_true(caps$supports_cv)
  expect_equal(caps$returns, "operator")

  # Cleanup
  unregister_aligner("gw")
})

test_that(".register_fpgw registers aligner correctly", {
  # Unregister first if exists
  tryCatch(unregister_aligner("fpgw"), error = function(e) NULL)

  # Register
  ensure_test_aligner("fpgw")

  # Verify registration
  expect_true("fpgw" %in% available_aligners())

  # Check capabilities
  caps <- aligner_capabilities("fpgw")
  expect_equal(caps$transform_type, "ot")
  expect_false(caps$mass_preserving)

  # Cleanup
  unregister_aligner("fpgw")
})

test_that(".coupling_to_operator handles larger matrices", {
  # 5x5 coupling matrix
  set.seed(123)
  P <- matrix(runif(25), 5, 5)
  P <- P / sum(P)  # Normalize to valid coupling

  op <- neuralign:::.coupling_to_operator(P)

  expect_equal(dim(op), c(5, 5))
  # Each row should sum to 1
  expect_equal(rowSums(op), rep(1, 5), tolerance = 1e-10)
})

test_that(".coupling_to_operator handles uniform coupling", {
  # Uniform coupling (all equal)
  P <- matrix(1/9, 3, 3)

  op <- neuralign:::.coupling_to_operator(P)

  expect_equal(dim(op), c(3, 3))
  # Each row should sum to 1
  expect_equal(rowSums(op), rep(1, 3), tolerance = 1e-10)
  # All elements should be equal (1/3)
  expect_equal(unique(as.vector(op)), 1/3, tolerance = 1e-10)
})

test_that(".coupling_to_operator handles sparse coupling", {
  # Sparse coupling (permutation-like)
  P <- matrix(0, 3, 3)
  P[1, 2] <- 1
  P[2, 3] <- 1
  P[3, 1] <- 1

  op <- neuralign:::.coupling_to_operator(P)

  expect_equal(dim(op), c(3, 3))
  # Result should be a permutation matrix (transposed)
  expect_equal(rowSums(op), rep(1, 3))
  expect_equal(colSums(op), rep(1, 3))
})

test_that(".gw_capabilities has all required fields", {
  caps <- neuralign:::.gw_capabilities

  # Check all expected fields exist
  expect_true("supports_cv" %in% names(caps))
  expect_true("cv_axes" %in% names(caps))
  expect_true("needs_geometry" %in% names(caps))
  expect_true("needs_design" %in% names(caps))
  expect_true("returns_invertible" %in% names(caps))
  expect_true("transform_type" %in% names(caps))
  expect_true("mass_preserving" %in% names(caps))
  expect_true("returns" %in% names(caps))
  expect_true("supports_new_subject" %in% names(caps))
  expect_true("supports_new_data" %in% names(caps))
  expect_true("reference_types" %in% names(caps))

  # Check specific values
  expect_equal(caps$cv_axes, c("subject"))
  expect_false(caps$needs_geometry)
  expect_false(caps$needs_design)
  expect_true(caps$supports_new_subject)
  expect_true(caps$supports_new_data)
  expect_true("consensus" %in% caps$reference_types)
})

test_that(".fpgw_capabilities has all required fields", {
  caps <- neuralign:::.fpgw_capabilities

  # Check all expected fields exist
  expect_true("supports_cv" %in% names(caps))
  expect_true("cv_axes" %in% names(caps))
  expect_true("needs_geometry" %in% names(caps))
  expect_true("needs_design" %in% names(caps))
  expect_true("returns_invertible" %in% names(caps))
  expect_true("transform_type" %in% names(caps))
  expect_true("mass_preserving" %in% names(caps))
  expect_true("returns" %in% names(caps))
  expect_true("supports_new_subject" %in% names(caps))
  expect_true("supports_new_data" %in% names(caps))
  expect_true("reference_types" %in% names(caps))

  # Check specific values
  expect_equal(caps$cv_axes, c("subject"))
  expect_false(caps$needs_geometry)
  expect_false(caps$needs_design)
  expect_false(caps$returns_invertible)
  expect_true(caps$supports_new_subject)
})
