# Algebra tests keep the old NULL-space fixtures; Stage 1 compose is
# fail-closed unless the caller opts into unverified spaces.
compose_alignment <- function(model1, model2, ...) {
  neuralign::compose_alignment(
    model1,
    model2,
    allow_unverified_spaces = TRUE,
    ...
  )
}
check_composition <- function(model1, model2, ...) {
  neuralign::check_composition(
    model1,
    model2,
    allow_unverified_spaces = TRUE,
    ...
  )
}

test_that("compose_alignment combines models", {
  # Create two simple models with compatible transforms
  transforms1 <- list(
    "sub-01" = matrix(c(1, 0, 0, 1), 2, 2),  # Identity
    "sub-02" = matrix(c(0, 1, 1, 0), 2, 2)   # Swap rows
  )

  transforms2 <- list(
    "sub-01" = matrix(c(2, 0, 0, 2), 2, 2),  # Scale by 2
    "sub-02" = matrix(c(1, 0, 0, 1), 2, 2)   # Identity
  )

  model1 <- AlignmentModel(
    transforms = transforms1,
    reference = "consensus",
    method = "method1"
  )

  model2 <- AlignmentModel(
    transforms = transforms2,
    reference = "consensus",
    method = "method2"
  )

  composed <- compose_alignment(model1, model2)

  expect_s4_class(composed, "AlignmentModel")
  expect_equal(composed@method, "method1+method2")

  # Check composed transform: T2 %*% T1
  # For sub-01: scale(2) %*% identity = scale(2)
  expect_equal(
    composed@transforms[["sub-01"]],
    matrix(c(2, 0, 0, 2), 2, 2)
  )

  # For sub-02: identity %*% swap = swap
  expect_equal(
    composed@transforms[["sub-02"]],
    matrix(c(0, 1, 1, 0), 2, 2)
  )
})

test_that("compose_alignment supports low-rank operator transforms", {
  # Build a rank-1 operator T = U V^T
  U <- matrix(c(1, 2, 3), 3, 1)
  V <- matrix(c(4, 5, 6), 3, 1)
  lr <- neuralign:::.new_low_rank_transform(U, V)

  mat <- diag(3) * 2

  X <- matrix(1:15, 3, 5)

  # matrix %*% low-rank (model1 first, model2 second)
  m1 <- AlignmentModel(list("sub-01" = lr), reference = NULL, method = "m1")
  m2 <- AlignmentModel(list("sub-01" = mat), reference = NULL, method = "m2")
  composed <- compose_alignment(m1, m2)
  expect_true(inherits(composed@transforms[["sub-01"]], "neuralign_low_rank_transform"))
  expect_equal(
    apply_transform(composed@transforms[["sub-01"]], X),
    apply_transform(mat, apply_transform(lr, X))
  )

  # low-rank %*% matrix (swap order)
  m3 <- AlignmentModel(list("sub-01" = mat), reference = NULL, method = "m3")
  m4 <- AlignmentModel(list("sub-01" = lr), reference = NULL, method = "m4")
  composed2 <- compose_alignment(m3, m4)
  expect_true(inherits(composed2@transforms[["sub-01"]], "neuralign_low_rank_transform"))
  expect_equal(
    apply_transform(composed2@transforms[["sub-01"]], X),
    apply_transform(lr, apply_transform(mat, X))
  )
})

test_that("compose_alignment supports low-rank × low-rank composition", {
  U1 <- matrix(c(1, 2, 3), 3, 1)
  V1 <- matrix(c(4, 5, 6), 3, 1)
  lr1 <- neuralign:::.new_low_rank_transform(U1, V1)

  U2 <- matrix(c(1, 0,
                 0, 1,
                 0, 0), 3, 2, byrow = TRUE)
  V2 <- matrix(c(0, 0,
                 1, 0,
                 0, 1), 3, 2, byrow = TRUE)
  lr2 <- neuralign:::.new_low_rank_transform(U2, V2)

  X <- matrix(1:12, 3, 4)

  m1 <- AlignmentModel(list("sub-01" = lr1), reference = NULL, method = "m1")
  m2 <- AlignmentModel(list("sub-01" = lr2), reference = NULL, method = "m2")
  composed <- compose_alignment(m1, m2)
  expect_true(inherits(composed@transforms[["sub-01"]], "neuralign_low_rank_transform"))
  expect_equal(
    apply_transform(composed@transforms[["sub-01"]], X),
    apply_transform(lr2, apply_transform(lr1, X))
  )
})

test_that("compose_alignment via %*% operator works", {
  transforms1 <- list(
    "sub-01" = diag(3)
  )
  transforms2 <- list(
    "sub-01" = diag(3) * 2
  )

  model1 <- AlignmentModel(
    transforms1, reference = NULL, method = "m1",
    space_from = "A", space_to = "B"
  )
  model2 <- AlignmentModel(
    transforms2, reference = NULL, method = "m2",
    space_from = "B", space_to = "C"
  )

  # model2 %*% model1 means model1 first, then model2
  composed <- model2 %*% model1

  expect_equal(composed@transforms[["sub-01"]], diag(3) * 2)
})

test_that("compose_alignment requires common subjects", {
  transforms1 <- list("sub-01" = diag(3))
  transforms2 <- list("sub-02" = diag(3))

  model1 <- AlignmentModel(transforms1, reference = NULL, method = "m1")
  model2 <- AlignmentModel(transforms2, reference = NULL, method = "m2")

  expect_error(
    compose_alignment(model1, model2),
    "no subjects in common"
  )
})

test_that("compose_alignment errors on partial overlap unless allow_partial", {
  transforms1 <- list(
    "sub-01" = diag(3),
    "sub-02" = diag(3)
  )
  transforms2 <- list(
    "sub-01" = diag(3),
    "sub-03" = diag(3)
  )

  model1 <- AlignmentModel(transforms1, reference = NULL, method = "m1")
  model2 <- AlignmentModel(transforms2, reference = NULL, method = "m2")

  expect_error(
    compose_alignment(model1, model2),
    "Partial subject drop|allow_partial"
  )

  composed <- compose_alignment(model1, model2, allow_partial = TRUE)
  expect_equal(names(composed@transforms), "sub-01")
})

test_that("compose_alignment checks dimension compatibility", {
  transforms1 <- list("sub-01" = matrix(1, 3, 5))  # 3x5
  transforms2 <- list("sub-01" = matrix(1, 4, 2))  # 4x2, incompatible

  model1 <- AlignmentModel(transforms1, reference = NULL, method = "m1")
  model2 <- AlignmentModel(transforms2, reference = NULL, method = "m2")

  expect_error(
    compose_alignment(model1, model2),
    "mismatch"
  )
})

test_that("check_composition returns compatibility info", {
  # Compatible models
  t1 <- list("sub-01" = matrix(1, 5, 10))  # 5x10
  t2 <- list("sub-01" = matrix(1, 3, 5))   # 3x5, output of t1 (5) matches input of t2 (5)

  m1 <- AlignmentModel(t1, reference = NULL, method = "m1")
  m2 <- AlignmentModel(t2, reference = NULL, method = "m2")

  result <- check_composition(m1, m2)
  expect_true(result$compatible)

  # Incompatible models
  t3 <- list("sub-01" = matrix(1, 3, 7))  # Input expects 7, not 5

  m3 <- AlignmentModel(t3, reference = NULL, method = "m3")
  result2 <- check_composition(m1, m3)
  expect_false(result2$compatible)
})

test_that("AlignmentModel %*% matrix applies transform", {
  transforms <- list(
    "sub-01" = matrix(c(2, 0, 0, 2), 2, 2)  # Scale by 2
  )

  model <- AlignmentModel(transforms, reference = NULL, method = "test")

  data <- matrix(c(1, 2, 3, 4), 2, 2)
  result <- model %*% data

  # Should apply scaling
  expected <- transforms[["sub-01"]] %*% data
  expect_equal(result, expected)
})

test_that("AlignmentModel %*% matrix errors when model has multiple transforms", {
  transforms <- list(
    "sub-01" = diag(2),
    "sub-02" = diag(2) * 2
  )
  model <- AlignmentModel(transforms, reference = NULL, method = "test")
  expect_error(
    model %*% matrix(1, 2, 2),
    "multiple subjects"
  )
})


# --- Additional tests for uncovered lines ---

test_that("compose_alignment accepts AlignmentResult arguments", {
  transforms1 <- list(
    "sub-01" = diag(3),
    "sub-02" = diag(3)
  )
  transforms2 <- list(
    "sub-01" = diag(3) * 2,
    "sub-02" = diag(3) * 3
  )

  model1 <- AlignmentModel(transforms1, reference = "consensus", method = "m1")
  model2 <- AlignmentModel(transforms2, reference = "consensus", method = "m2")

  # Wrap in AlignmentResult
  result1 <- AlignmentResult(
    model = model1,
    aligned = list("sub-01" = diag(3), "sub-02" = diag(3))
  )
  result2 <- AlignmentResult(
    model = model2,
    aligned = list("sub-01" = diag(3) * 2, "sub-02" = diag(3) * 3)
  )

  # Both arguments are AlignmentResult — should extract models internally

  composed <- compose_alignment(result1, result2)
  expect_s4_class(composed, "AlignmentModel")
  expect_equal(composed@method, "m1+m2")
  expect_equal(composed@transforms[["sub-01"]], diag(3) * 2)
  expect_equal(composed@transforms[["sub-02"]], diag(3) * 3)

  # One AlignmentResult and one AlignmentModel
  composed2 <- compose_alignment(result1, model2)
  expect_s4_class(composed2, "AlignmentModel")
  expect_equal(composed2@transforms[["sub-01"]], diag(3) * 2)
})

test_that("compose_alignment errors on non-model inputs", {
  model1 <- AlignmentModel(
    list("sub-01" = diag(3)),
    reference = NULL,
    method = "m1"
  )

  # Pass a string as model1

  expect_error(
    compose_alignment("not_a_model", model1),
    "must be an AlignmentModel"
  )

  # Pass a number as model2
  expect_error(
    compose_alignment(model1, 42),
    "must be an AlignmentModel"
  )

  # Pass a list (not an AlignmentModel or AlignmentResult)
  expect_error(
    compose_alignment(list(a = 1), model1),
    "must be an AlignmentModel"
  )
})

test_that("compose_alignment errors on non-operator transforms", {
  expect_error(
    AlignmentModel(list("sub-01" = function(x) x), reference = NULL, method = "m1"),
    "unsupported|operator|transform"
  )
  expect_error(
    AlignmentModel(list("sub-01" = "not_a_matrix"), reference = NULL, method = "m4"),
    "unsupported|operator|transform"
  )
})

test_that("check_composition handles AlignmentResult inputs", {
  transforms1 <- list("sub-01" = matrix(1, 5, 10))
  transforms2 <- list("sub-01" = matrix(1, 3, 5))

  model1 <- AlignmentModel(transforms1, reference = NULL, method = "m1")
  model2 <- AlignmentModel(transforms2, reference = NULL, method = "m2")

  result1 <- AlignmentResult(
    model = model1,
    aligned = list("sub-01" = matrix(1, 5, 10))
  )
  result2 <- AlignmentResult(
    model = model2,
    aligned = list("sub-01" = matrix(1, 3, 5))
  )

  # Both AlignmentResult
  check <- check_composition(result1, result2)
  expect_true(check$compatible)

  # Mixed: AlignmentResult and AlignmentModel
  check2 <- check_composition(result1, model2)
  expect_true(check2$compatible)
})

test_that("check_composition detects no subjects in common", {
  transforms1 <- list("sub-01" = diag(3))
  transforms2 <- list("sub-99" = diag(3))

  model1 <- AlignmentModel(transforms1, reference = NULL, method = "m1")
  model2 <- AlignmentModel(transforms2, reference = NULL, method = "m2")

  result <- check_composition(model1, model2)
  expect_false(result$compatible)
  expect_match(result$message, "No subjects in common")
})

test_that("check_composition detects non-operator transforms", {
  expect_error(
    AlignmentModel(list("sub-01" = function(x) x), reference = NULL, method = "m1"),
    "unsupported|operator|transform"
  )
})

test_that("check_composition detects dimension mismatch", {
  # model1 output rows = 5, model2 input cols = 7 -> mismatch
  transforms1 <- list("sub-01" = matrix(1, 5, 10))
  transforms2 <- list("sub-01" = matrix(1, 3, 7))

  model1 <- AlignmentModel(transforms1, reference = NULL, method = "m1")
  model2 <- AlignmentModel(transforms2, reference = NULL, method = "m2")

  result <- check_composition(model1, model2)
  expect_false(result$compatible)
  expect_match(result$message, "Dimension mismatch")
})

test_that("check_composition returns info message for compatible models", {
  transforms1 <- list(
    "sub-01" = matrix(1, 5, 10),
    "sub-02" = matrix(1, 5, 10)
  )
  transforms2 <- list(
    "sub-01" = matrix(1, 3, 5),
    "sub-02" = matrix(1, 3, 5)
  )

  model1 <- AlignmentModel(transforms1, reference = NULL, method = "m1")
  model2 <- AlignmentModel(transforms2, reference = NULL, method = "m2")

  result <- check_composition(model1, model2)
  expect_true(result$compatible)
  expect_match(result$message, "compatible")
  expect_match(result$message, "2 common subjects")
})


# ---------- Provenance preservation ----------

test_that("compose_alignment preserves provenance from both models", {
  transforms1 <- list("sub-01" = diag(3))
  transforms2 <- list("sub-01" = diag(3) * 2)

  model1 <- AlignmentModel(
    transforms1, reference = "consensus", method = "m1",
    params = list(scale = TRUE)
  )
  model2 <- AlignmentModel(
    transforms2, reference = "consensus", method = "m2",
    params = list(alpha = 0.5)
  )

  composed <- compose_alignment(model1, model2)

  # Provenance should have composed_from with both models
  prov <- composed@provenance
  expect_true(!is.null(prov$composed_from))
  expect_true(!is.null(prov$composed_from$model1))
  expect_true(!is.null(prov$composed_from$model2))
  expect_equal(prov$composed_from$model1$method, "m1")
  expect_equal(prov$composed_from$model2$method, "m2")
  expect_true(!is.null(prov$composed_at))
  expect_true(!is.null(prov$neuralign_version))
})

test_that("compose_alignment errors on space chain mismatch", {
  transforms1 <- list("sub-01" = diag(3))
  transforms2 <- list("sub-01" = diag(3))

  model1 <- AlignmentModel(
    transforms1, reference = "consensus", method = "m1",
    space_from = "native", space_to = "MNI"
  )
  model2 <- AlignmentModel(
    transforms2, reference = "consensus", method = "m2",
    space_from = "talairach", space_to = "functional"
  )

  expect_error(
    neuralign::compose_alignment(model1, model2),
    "Space chain mismatch"
  )
})

test_that("compose_alignment sets correct space_from and space_to", {
  transforms <- list("sub-01" = diag(3))

  model1 <- AlignmentModel(
    transforms, reference = "consensus", method = "m1",
    space_from = "native", space_to = "MNI"
  )
  model2 <- AlignmentModel(
    transforms, reference = "consensus", method = "m2",
    space_from = "MNI", space_to = "functional"
  )

  composed <- compose_alignment(model1, model2)
  expect_equal(composed@space_from, "native")
  expect_equal(composed@space_to, "functional")
})


# ---------- More compose coverage tests ----------

test_that("compose_alignment with unregistered method skips capability check", {
  transforms1 <- list("sub-01" = diag(3))
  transforms2 <- list("sub-01" = diag(3) * 2)

  # Use unregistered method names - should skip capability check
  model1 <- AlignmentModel(transforms1, reference = "consensus", method = "unregistered_m1")
  model2 <- AlignmentModel(transforms2, reference = "consensus", method = "unregistered_m2")

  composed <- compose_alignment(model1, model2)
  expect_s4_class(composed, "AlignmentModel")
  expect_equal(composed@method, "unregistered_m1+unregistered_m2")
})

test_that("compose_alignment method_state merges both models", {
  transforms <- list("sub-01" = diag(3))

  model1 <- AlignmentModel(
    transforms, reference = "consensus", method = "m1",
    method_state = list(scale = TRUE)
  )
  model2 <- AlignmentModel(
    transforms, reference = "consensus", method = "m2",
    method_state = list(reflection = FALSE)
  )

  composed <- compose_alignment(model1, model2)
  expect_equal(composed@method_state$model1_state$scale, TRUE)
  expect_equal(composed@method_state$model2_state$reflection, FALSE)
})

test_that("compose_alignment train_subjects is common subjects", {
  transforms1 <- list("sub-01" = diag(3), "sub-02" = diag(3))
  transforms2 <- list("sub-01" = diag(3), "sub-02" = diag(3))

  model1 <- AlignmentModel(transforms1, reference = "consensus", method = "m1")
  model2 <- AlignmentModel(transforms2, reference = "consensus", method = "m2")

  composed <- compose_alignment(model1, model2)
  expect_setequal(composed@train_subjects, c("sub-01", "sub-02"))
})


# ---------- compose_alignment non-operator capability error ----------

test_that("compose_alignment errors when model1 has non-operator returns", {
  neuralign:::.clear_registry()

  # Register a method that claims it returns something other than "operator"
  # We can't register "embedding" (rejected), so modify directly after registration
  dummy_fit <- function(data, reference, train_idx = NULL, fit_context = NULL, provider_plan = NULL, ...) {
    list(transforms = list(), reference_data = NULL)
  }
  register_aligner("op_method", dummy_fit)

  # Manually override the returns capability for testing
  reg_env <- get(".aligner_registry", envir = asNamespace("neuralign"))
  entry <- reg_env[["op_method"]]
  entry$capabilities$returns <- "some_other_type"
  assign("op_method", entry, envir = reg_env)

  model1 <- AlignmentModel(
    transforms = list("s1" = diag(3)),
    reference = "s1",
    method = "op_method"
  )
  model2 <- AlignmentModel(
    transforms = list("s1" = diag(3)),
    reference = "s1",
    method = "m2_unregistered"
  )

  expect_error(
    compose_alignment(model1, model2),
    "does not return operator transforms"
  )

  neuralign:::.clear_registry()
})

test_that("compose_alignment errors when model2 has non-operator returns", {
  neuralign:::.clear_registry()

  dummy_fit <- function(data, reference, train_idx = NULL, fit_context = NULL, provider_plan = NULL, ...) {
    list(transforms = list(), reference_data = NULL)
  }
  register_aligner("op_method2", dummy_fit)

  reg_env <- get(".aligner_registry", envir = asNamespace("neuralign"))
  entry <- reg_env[["op_method2"]]
  entry$capabilities$returns <- "some_other_type"
  assign("op_method2", entry, envir = reg_env)

  model1 <- AlignmentModel(
    transforms = list("s1" = diag(3)),
    reference = "s1",
    method = "m1_unregistered"
  )
  model2 <- AlignmentModel(
    transforms = list("s1" = diag(3)),
    reference = "s1",
    method = "op_method2"
  )

  expect_error(
    compose_alignment(model1, model2),
    "does not return operator transforms"
  )

  neuralign:::.clear_registry()
})

test_that("AlignmentModel %*% matrix errors on empty model", {
  model <- new("AlignmentModel",
    transforms = list(),
    reference = NULL,
    method = "test",
    space_from = NULL,
    space_to = NULL,
    provenance = list(),
    method_state = list(),
    train_subjects = character(0)
  )

  expect_error(
    model %*% matrix(1, 3, 3),
    "no transforms"
  )
})
