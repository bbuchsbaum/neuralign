test_that("lift_aligned lifts into a target subject space via pinv", {
  k <- 3
  p <- 7
  n <- 11

  A_ref <- matrix(rnorm(k * p), k, p)
  A_s2 <- matrix(rnorm(k * p), k, p)

  Z_ref <- matrix(rnorm(k * n), k, n)
  Z_s2 <- matrix(rnorm(k * n), k, n)

  model <- AlignmentModel(
    transforms = list(s1 = A_ref, s2 = A_s2),
    reference = "s1",
    method = "dummy"
  )
  res <- AlignmentResult(model = model, aligned = list(s1 = Z_ref, s2 = Z_s2))

  inv_ref <- inverse_transform(model, "s1", method = "pinv")
  expected <- list(
    s1 = apply_transform(inv_ref, Z_ref),
    s2 = apply_transform(inv_ref, Z_s2)
  )

  out <- lift_aligned(res, to = "subject", subject = "s1", inverse = "pinv")
  expect_equal(out, expected)
})

test_that("lift_aligned can lift each aligned matrix into its own subject space", {
  k <- 2
  p <- 5
  n <- 6

  A_s1 <- matrix(rnorm(k * p), k, p)
  A_s2 <- matrix(rnorm(k * p), k, p)

  Z_s1 <- matrix(rnorm(k * n), k, n)
  Z_s2 <- matrix(rnorm(k * n), k, n)

  model <- AlignmentModel(
    transforms = list(s1 = A_s1, s2 = A_s2),
    reference = "s1",
    method = "dummy"
  )
  res <- AlignmentResult(model = model, aligned = list(s1 = Z_s1, s2 = Z_s2))

  out <- lift_aligned(res, to = "each", inverse = "pinv")
  expect_equal(out$s1, apply_transform(inverse_transform(model, "s1", method = "pinv"), Z_s1))
  expect_equal(out$s2, apply_transform(inverse_transform(model, "s2", method = "pinv"), Z_s2))
})

test_that("lift_aligned to='reference' uses model reference subject", {
  k <- 3
  p <- 4
  n <- 8

  A_ref <- matrix(rnorm(k * p), k, p)
  Z <- matrix(rnorm(k * n), k, n)

  model <- AlignmentModel(
    transforms = list(ref = A_ref),
    reference = "ref",
    method = "dummy"
  )
  res <- AlignmentResult(model = model, aligned = list(ref = Z))

  out <- lift_aligned(res, to = "reference", inverse = "pinv")
  expect_equal(out$ref, apply_transform(inverse_transform(model, "ref", method = "pinv"), Z))
})

