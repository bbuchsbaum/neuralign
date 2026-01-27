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

# ---------- NEW TESTS ----------

test_that(".fit_transform_for_subject uses custom apply_fn when provided", {
  # Register an aligner that has a custom apply_fn
  register_aligner(
    name = "custom_apply_test",
    fit_fn = function(data, reference, train_idx = NULL, ...) {
      if (is.null(train_idx)) train_idx <- seq_along(data@subjects)
      train_data <- data[train_idx]
      subjects <- train_data@subjects
      transforms <- setNames(
        lapply(subjects, function(s) diag(5)),
        subjects
      )
      list(
        transforms = transforms,
        reference_data = matrix(0, 5, 3),
        space_from = NULL,
        space_to = NULL,
        method_state = list(custom_state = TRUE)
      )
    },
    apply_fn = function(fit_result, new_data, ...) {
      # Custom apply: return identity scaled by 2
      n <- nrow(get_subject_data(new_data, new_data@subjects[1]))
      list(transforms = list(diag(n) * 2))
    },
    capabilities = list(
      supports_new_subject = TRUE,
      returns = "operator"
    ),
    package = "neuralign"
  )

  set.seed(1234)
  train_data <- list(
    "sub-01" = matrix(rnorm(15), 5, 3),
    "sub-02" = matrix(rnorm(15), 5, 3)
  )
  train_adat <- AlignmentData(train_data)

  result <- fit_alignment(train_adat, method = "custom_apply_test")

  # New subject to trigger apply_fn path
  new_data <- list("sub-03" = matrix(rnorm(15), 5, 3))
  new_adat <- AlignmentData(new_data)

  applied <- apply_alignment(result, new_adat, warn_leakage = FALSE)

  # The custom apply_fn returns identity * 2, so the transform should be 2*I
  new_transform <- get_transform(get_model(applied), "sub-03")
  expect_equal(new_transform, diag(5) * 2)

  unregister_aligner("custom_apply_test")
})

test_that("inverse_transform method=auto errors when caps say returns_invertible=FALSE", {
  register_aligner(
    name = "noninvertible_test",
    fit_fn = function(data, reference, ...) {
      list(transforms = list("x" = diag(3)), reference_data = NULL)
    },
    capabilities = list(
      returns_invertible = FALSE,
      transform_type = "linear",
      returns = "operator"
    ),
    package = "neuralign"
  )

  model <- AlignmentModel(
    transforms = list("sub-01" = diag(3)),
    reference = "consensus",
    method = "noninvertible_test"
  )

  expect_error(
    inverse_transform(model, "sub-01", method = "auto"),
    "does not declare invertible"
  )

  unregister_aligner("noninvertible_test")
})

test_that("inverse_transform method=auto uses solve for linear invertible transforms", {
  register_aligner(
    name = "linear_invertible_test",
    fit_fn = function(data, reference, ...) {
      list(transforms = list("x" = diag(3)), reference_data = NULL)
    },
    capabilities = list(
      returns_invertible = TRUE,
      transform_type = "linear",
      returns = "operator"
    ),
    package = "neuralign"
  )

  # Non-orthogonal but invertible
  mat <- matrix(c(1, 2, 0, 0, 1, 0, 0, 0, 3), 3, 3)
  model <- AlignmentModel(
    transforms = list("sub-01" = mat),
    reference = "consensus",
    method = "linear_invertible_test"
  )

  inv <- inverse_transform(model, "sub-01", method = "auto")
  # "auto" should resolve to "solve" for linear invertible

  expect_equal(mat %*% inv, diag(3), tolerance = 1e-10)

  unregister_aligner("linear_invertible_test")
})

test_that("inverse_transform method=solve errors on non-square transform", {
  # 2x3 non-square matrix
  transforms <- list("sub-01" = matrix(1:6, 2, 3))
  model <- AlignmentModel(
    transforms = transforms,
    reference = "consensus",
    method = "test"
  )

  expect_error(
    inverse_transform(model, "sub-01", method = "solve"),
    "not square"
  )
})

test_that(".pseudoinverse on a 1x1 zero matrix returns zero matrix", {
  # The !length(d) guard in .pseudoinverse is unreachable for truly 0x0 matrices

  # because R's svd() errors on zero-dimension inputs. Test the near-zero path
  # where all singular values fall below the cutoff.
  zero_mat <- matrix(0, 1, 1)
  result <- neuralign:::.pseudoinverse(zero_mat)
  expect_equal(dim(result), c(1, 1))
  expect_equal(result[1, 1], 0)
})

test_that(".pseudoinverse on a rank-deficient matrix zeroes small singular values", {
  # 3x3 rank-1 matrix
  v <- c(1, 2, 3)
  rank1 <- outer(v, v)
  result <- neuralign:::.pseudoinverse(rank1)
  expect_equal(dim(result), c(3, 3))
  # rank1 %*% result %*% rank1 should equal rank1 (Moore-Penrose condition)
  expect_equal(rank1 %*% result %*% rank1, rank1, tolerance = 1e-10)
})

test_that(".ridge_inverse errors on invalid lambda values", {
  mat <- diag(3)

  # Negative lambda
  expect_error(
    neuralign:::.ridge_inverse(mat, lambda = -1),
    "lambda.*must be a single positive number"
  )

  # Non-numeric lambda
  expect_error(
    neuralign:::.ridge_inverse(mat, lambda = "abc"),
    "lambda.*must be a single positive number"
  )

  # NA lambda
  expect_error(
    neuralign:::.ridge_inverse(mat, lambda = NA_real_),
    "lambda.*must be a single positive number"
  )

  # Zero lambda
  expect_error(
    neuralign:::.ridge_inverse(mat, lambda = 0),
    "lambda.*must be a single positive number"
  )

  # Vector lambda
  expect_error(
    neuralign:::.ridge_inverse(mat, lambda = c(0.1, 0.2)),
    "lambda.*must be a single positive number"
  )

  # Inf lambda
  expect_error(
    neuralign:::.ridge_inverse(mat, lambda = Inf),
    "lambda.*must be a single positive number"
  )
})

test_that(".ridge_inverse wide matrix path (n_out < n_in)", {
  # Wide matrix: 2 rows x 5 cols => n_out=2, n_in=5, so n_out < n_in
  set.seed(42)
  wide_mat <- matrix(rnorm(10), 2, 5)

  result <- neuralign:::.ridge_inverse(wide_mat, lambda = 0.01)
  # Result should be 5 x 2 (n_in x n_out)
  expect_equal(dim(result), c(5, 2))

  # Verify the inverse is approximately correct: wide_mat %*% result ~ I_2
  product <- wide_mat %*% result
  expect_equal(dim(product), c(2, 2))
  # With ridge regularization, should be approximately identity
  expect_true(all(abs(diag(product) - 1) < 0.1))
})

test_that(".ridge_inverse tall matrix path (n_out >= n_in)", {
  # Tall matrix: 5 rows x 2 cols => n_out=5, n_in=2, so n_out >= n_in
  set.seed(43)
  tall_mat <- matrix(rnorm(10), 5, 2)

  result <- neuralign:::.ridge_inverse(tall_mat, lambda = 0.01)
  # Result should be 2 x 5 (n_in x n_out)
  expect_equal(dim(result), c(2, 5))

  # Verify: result %*% tall_mat ~ I_2
  product <- result %*% tall_mat
  expect_equal(dim(product), c(2, 2))
  expect_true(all(abs(diag(product) - 1) < 0.1))
})

test_that(".as_dense_matrix converts sparse Matrix to dense", {
  skip_if_not_installed("Matrix")
  sparse <- Matrix::sparseMatrix(i = c(1, 2, 3), j = c(1, 2, 3), x = c(1, 2, 3))
  result <- neuralign:::.as_dense_matrix(sparse)
  expect_true(is.matrix(result))
  expect_false(inherits(result, "Matrix"))
  expect_equal(result, matrix(c(1, 0, 0, 0, 2, 0, 0, 0, 3), 3, 3))
})

test_that(".as_dense_matrix converts data.frame to matrix", {
  df <- data.frame(a = 1:3, b = 4:6)
  result <- neuralign:::.as_dense_matrix(df)
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(3, 2))
  expect_equal(result[, 1], 1:3)
  expect_equal(result[, 2], 4:6)
})

test_that(".as_dense_matrix returns plain matrix unchanged", {
  mat <- matrix(1:6, 2, 3)
  result <- neuralign:::.as_dense_matrix(mat)
  expect_true(is.matrix(result))
  expect_identical(result, mat)
})


# ---------- Space mismatch warning in apply_alignment ----------

test_that("apply_alignment warns on space mismatch", {
  neuralign:::.register_procrustes()

  set.seed(42)
  data_list <- list(
    "sub-01" = matrix(rnorm(50), 10, 5),
    "sub-02" = matrix(rnorm(50), 10, 5)
  )
  adat <- AlignmentData(data_list)

  result <- fit_alignment(adat, method = "procrustes", reference = "sub-01")
  model <- get_model(result)

  # Force different space_from on the model
  model@space_from <- "MNI152"

  new_data <- list("sub-03" = matrix(rnorm(50), 10, 5))
  new_adat <- AlignmentData(new_data, space = "native")

  expect_warning(
    apply_alignment(model, new_adat, warn_leakage = FALSE),
    "Space mismatch"
  )
})


# ---------- fold_specific anchor and inverse_transform ----------

test_that("apply_alignment errors on fold_specific reference with new subjects", {
  neuralign:::.register_procrustes()

  model <- new("AlignmentModel",
    transforms = list("sub-01" = diag(5)),
    reference = "fold_specific",
    reference_data = NULL,
    method = "procrustes",
    space_from = NULL,
    space_to = NULL,
    provenance = list(),
    method_state = list(),
    train_subjects = "sub-01"
  )

  new_data <- AlignmentData(list("sub-99" = matrix(rnorm(25), 5, 5)))
  expect_error(apply_alignment(model, new_data), "fold-specific anchors")
})

test_that("apply_alignment with fit_new=FALSE warns for new subjects", {
  neuralign:::.register_procrustes()

  set.seed(55)
  d <- 5; n <- 4
  model <- new("AlignmentModel",
    transforms = list("sub-01" = diag(d)),
    reference = "sub-01",
    reference_data = matrix(rnorm(d * n), d, n),
    method = "procrustes",
    space_from = NULL,
    space_to = NULL,
    provenance = list(),
    method_state = list(),
    train_subjects = "sub-01"
  )

  new_data <- AlignmentData(list("sub-99" = matrix(rnorm(d * n), d, n)))
  expect_warning(
    apply_alignment(model, new_data, fit_new = FALSE, warn_leakage = FALSE),
    "fit_new=FALSE"
  )
})

test_that("inverse_transform errors on OT-style transform_type", {
  neuralign:::.clear_registry()
  dummy_fit <- function(data, reference, ...) {
    list(transforms = list(), reference_data = NULL)
  }
  register_aligner("ot_method", dummy_fit,
    capabilities = list(transform_type = "ot"))

  model <- new("AlignmentModel",
    transforms = list("s1" = diag(5)),
    reference = "s1", reference_data = NULL,
    method = "ot_method",
    space_from = NULL, space_to = NULL,
    provenance = list(), method_state = list(),
    train_subjects = "s1"
  )

  expect_error(inverse_transform(model, "s1"), "OT-style couplings")
  neuralign:::.clear_registry()
})

test_that("inverse_transform auto errors when not invertible", {
  neuralign:::.clear_registry()
  dummy_fit <- function(data, reference, ...) {
    list(transforms = list(), reference_data = NULL)
  }
  register_aligner("non_inv", dummy_fit,
    capabilities = list(returns_invertible = FALSE, transform_type = "linear"))

  model <- new("AlignmentModel",
    transforms = list("s1" = diag(5)),
    reference = "s1", reference_data = NULL,
    method = "non_inv",
    space_from = NULL, space_to = NULL,
    provenance = list(), method_state = list(),
    train_subjects = "s1"
  )

  expect_error(inverse_transform(model, "s1", method = "auto"), "does not declare invertible")
  neuralign:::.clear_registry()
})

test_that("inverse_transform pinv works for non-square matrix", {
  set.seed(56)
  W <- matrix(rnorm(50), 10, 5)
  model <- new("AlignmentModel",
    transforms = list("s1" = W),
    reference = "s1", reference_data = NULL,
    method = "test_method",
    space_from = NULL, space_to = NULL,
    provenance = list(), method_state = list(),
    train_subjects = "s1"
  )

  inv <- inverse_transform(model, "s1", method = "pinv")
  expect_equal(nrow(inv), 5L)
  expect_equal(ncol(inv), 10L)
})

test_that("inverse_transform ridge works", {
  set.seed(57)
  W <- matrix(rnorm(50), 10, 5)
  model <- new("AlignmentModel",
    transforms = list("s1" = W),
    reference = "s1", reference_data = NULL,
    method = "test_method",
    space_from = NULL, space_to = NULL,
    provenance = list(), method_state = list(),
    train_subjects = "s1"
  )

  inv <- inverse_transform(model, "s1", method = "ridge", lambda = 0.01)
  expect_equal(nrow(inv), 5L)
  expect_equal(ncol(inv), 10L)
})

test_that("apply_transform errors on dimension mismatch", {
  expect_error(apply_transform(diag(5), matrix(1, 3, 4)), "dimension mismatch")
})

test_that("apply_transform errors on non-matrix transform", {
  expect_error(apply_transform("not_a_matrix", matrix(1, 3, 3)), "must be a matrix")
})


# ---------- Additional edge-case coverage ----------

test_that("apply_alignment errors on non-operator returns method", {
  neuralign:::.clear_registry()
  dummy_fit <- function(data, reference, ...) {
    list(transforms = list(), reference_data = NULL)
  }
  register_aligner("embed_method", dummy_fit)

  # Manually override returns capability
  reg_env <- get(".aligner_registry", envir = asNamespace("neuralign"))
  entry <- reg_env[["embed_method"]]
  entry$capabilities$returns <- "embedding"
  assign("embed_method", entry, envir = reg_env)

  model <- new("AlignmentModel",
    transforms = list("s1" = diag(5)),
    reference = "s1", reference_data = NULL,
    method = "embed_method",
    space_from = NULL, space_to = NULL,
    provenance = list(), method_state = list(),
    train_subjects = "s1"
  )

  new_data <- AlignmentData(list("s1" = matrix(rnorm(25), 5, 5)))
  expect_error(
    apply_alignment(model, new_data, warn_leakage = FALSE),
    "does not return operator transforms"
  )
  neuralign:::.clear_registry()
})

test_that("apply_alignment rejects non-model input", {
  expect_error(
    apply_alignment("not_a_model", AlignmentData(list(s = matrix(1, 3, 3)))),
    "must be an AlignmentModel"
  )
})

test_that("inverse_transform auto fallback to solve for square unregistered method", {
  # Unregistered method: caps is NULL, square transform → solve
  model <- new("AlignmentModel",
    transforms = list("s1" = diag(5)),
    reference = "s1", reference_data = NULL,
    method = "totally_unknown_method_xyz",
    space_from = NULL, space_to = NULL,
    provenance = list(), method_state = list(),
    train_subjects = "s1"
  )

  inv <- inverse_transform(model, "s1", method = "auto")
  expect_equal(inv, diag(5))
})

test_that("inverse_transform auto errors for non-square with unregistered method", {
  # Non-square + no capabilities → error
  model <- new("AlignmentModel",
    transforms = list("s1" = matrix(rnorm(15), 3, 5)),
    reference = "s1", reference_data = NULL,
    method = "totally_unknown_method_xyz",
    space_from = NULL, space_to = NULL,
    provenance = list(), method_state = list(),
    train_subjects = "s1"
  )

  expect_error(
    inverse_transform(model, "s1", method = "auto"),
    "not square"
  )
})

test_that("inverse_transform solve errors on singular matrix", {
  # Singular matrix
  singular <- matrix(c(1, 2, 2, 4), 2, 2)
  model <- new("AlignmentModel",
    transforms = list("s1" = singular),
    reference = "s1", reference_data = NULL,
    method = "totally_unknown_method_xyz",
    space_from = NULL, space_to = NULL,
    provenance = list(), method_state = list(),
    train_subjects = "s1"
  )

  expect_error(
    inverse_transform(model, "s1", method = "solve"),
    "not invertible"
  )
})

test_that("inverse_transform solve errors on non-square matrix", {
  model <- new("AlignmentModel",
    transforms = list("s1" = matrix(1, 3, 5)),
    reference = "s1", reference_data = NULL,
    method = "totally_unknown_method_xyz",
    space_from = NULL, space_to = NULL,
    provenance = list(), method_state = list(),
    train_subjects = "s1"
  )

  expect_error(
    inverse_transform(model, "s1", method = "solve"),
    "not square"
  )
})

test_that("inverse_transform unknown method errors", {
  model <- new("AlignmentModel",
    transforms = list("s1" = diag(3)),
    reference = "s1", reference_data = NULL,
    method = "test_method",
    space_from = NULL, space_to = NULL,
    provenance = list(), method_state = list(),
    train_subjects = "s1"
  )

  expect_error(
    inverse_transform(model, "s1", method = "bogus"),
    "Unknown inverse method"
  )
})

test_that("inverse_transform transpose returns correct result", {
  set.seed(58)
  Q <- qr.Q(qr(matrix(rnorm(25), 5, 5)))  # orthogonal
  model <- new("AlignmentModel",
    transforms = list("s1" = Q),
    reference = "s1", reference_data = NULL,
    method = "test_method",
    space_from = NULL, space_to = NULL,
    provenance = list(), method_state = list(),
    train_subjects = "s1"
  )

  inv <- inverse_transform(model, "s1", method = "transpose")
  expect_equal(inv, t(Q))
  # For orthogonal: Q^T * Q = I
  expect_equal(inv %*% Q, diag(5), tolerance = 1e-10)
})

test_that(".ridge_inverse errors on invalid lambda", {
  expect_error(neuralign:::.ridge_inverse(diag(3), lambda = -1), "positive number")
  expect_error(neuralign:::.ridge_inverse(diag(3), lambda = 0), "positive number")
  expect_error(neuralign:::.ridge_inverse(diag(3), lambda = Inf), "positive number")
  expect_error(neuralign:::.ridge_inverse(diag(3), lambda = "abc"), "positive number")
  expect_error(neuralign:::.ridge_inverse(diag(3), lambda = c(0.1, 0.2)), "positive number")
})

test_that(".ridge_inverse handles wide matrix (n_out < n_in)", {
  set.seed(59)
  # Wide: 3 rows, 8 cols
  W <- matrix(rnorm(24), 3, 8)
  inv <- neuralign:::.ridge_inverse(W, lambda = 0.01)
  expect_equal(nrow(inv), 8L)
  expect_equal(ncol(inv), 3L)
  # inv %*% W should approximate identity(3)
  product <- W %*% inv
  expect_equal(dim(product), c(3L, 3L))
})

test_that(".ridge_inverse handles square matrix", {
  set.seed(60)
  W <- matrix(rnorm(25), 5, 5)
  inv <- neuralign:::.ridge_inverse(W, lambda = 0.001)
  expect_equal(dim(inv), c(5L, 5L))
  # Should approximate the true inverse for well-conditioned matrix
  product <- inv %*% W
  expect_equal(product, diag(5), tolerance = 0.05)
})

test_that(".fit_transform_for_subject uses custom apply_fn when available", {
  neuralign:::.clear_registry()

  # Register aligner with a custom apply_fn
  dummy_fit <- function(data, reference, ...) {
    list(transforms = list(), reference_data = NULL)
  }
  custom_apply <- function(fit_result, new_data, ...) {
    # Return a known transform for any new subject
    list(transforms = list(diag(3) * 7))
  }
  register_aligner("custom_apply_method", dummy_fit,
    apply_fn = custom_apply)

  model <- new("AlignmentModel",
    transforms = list("s1" = diag(3)),
    reference = "s1",
    reference_data = matrix(rnorm(9), 3, 3),
    method = "custom_apply_method",
    space_from = NULL, space_to = NULL,
    provenance = list(), method_state = list(),
    train_subjects = "s1"
  )

  new_data <- AlignmentData(list("s_new" = matrix(rnorm(9), 3, 3)))
  result <- apply_alignment(model, new_data, warn_leakage = FALSE)

  # Should have used custom_apply, which returns 7 * identity
  expect_equal(result@model@transforms[["s_new"]], diag(3) * 7)
  neuralign:::.clear_registry()
})
