test_that("diag_cov_orthogonal identity and permutation behave as expected", {
  v <- c(a = 1, b = 4, c = 9)
  I <- diag(3)
  expect_equal(diag_cov_orthogonal(v, I, convention = "left"), v)
  expect_equal(diag_cov_orthogonal(v, I, convention = "right"), v)

  perm <- c(2, 3, 1)
  invperm <- match(seq_along(perm), perm)
  P <- diag(3)[perm, , drop = FALSE]

  out_left <- diag_cov_orthogonal(v, P, convention = "left")
  expect_equal(names(out_left), names(v))
  expect_equal(unname(out_left), unname(v[perm]))

  out_right <- diag_cov_orthogonal(v, P, convention = "right")
  expect_equal(names(out_right), names(v))
  expect_equal(unname(out_right), unname(v[invperm]))
})

test_that("diag_cov_orthogonal matches analytic 2D rotation case", {
  theta <- 0.3
  Q <- matrix(
    c(cos(theta), -sin(theta), sin(theta), cos(theta)),
    nrow = 2,
    byrow = TRUE
  )
  v <- c(1, 4)
  expected <- c(
    cos(theta)^2 * v[[1]] + sin(theta)^2 * v[[2]],
    sin(theta)^2 * v[[1]] + cos(theta)^2 * v[[2]]
  )

  out <- diag_cov_orthogonal(v, Q, convention = "left")
  expect_equal(out, expected, tolerance = 1e-12)
})

test_that("diag_cov_orthogonal supports multiple Q via list", {
  v <- c(1, 4, 9)
  Qs <- list(
    id = diag(3),
    perm = diag(3)[c(2, 3, 1), , drop = FALSE]
  )
  out <- diag_cov_orthogonal(v, Qs, convention = "left")
  expect_true(is.matrix(out))
  expect_equal(colnames(out), names(Qs))
  expect_equal(out[, "id"], v)
  expect_equal(unname(out[, "perm"]), unname(v[c(2, 3, 1)]))
})

test_that("diag_cov_orthogonal validates inputs", {
  expect_error(diag_cov_orthogonal(c(-1, 2), diag(2)), "non-negative")
  expect_error(diag_cov_orthogonal(c(1, 2, 3), diag(2)), "Dimension mismatch")
  expect_error(diag_cov_orthogonal(matrix(1, 2, 2), list(diag(2))), "must be a numeric vector")
})
