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

test_that("compose_alignment via %*% operator works", {
  transforms1 <- list(
    "sub-01" = diag(3)
  )
  transforms2 <- list(
    "sub-01" = diag(3) * 2
  )

  model1 <- AlignmentModel(transforms1, reference = NULL, method = "m1")
  model2 <- AlignmentModel(transforms2, reference = NULL, method = "m2")

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

test_that("compose_alignment warns on partial overlap", {
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

  expect_warning(
    compose_alignment(model1, model2),
    "in common"
  )
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
