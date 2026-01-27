test_that("identified_subspace_basis returns expected rank and orthonormal basis", {
  x <- outer(c(1, 2), c(1, 1, 1)) # 2 x 3, rank 1
  out <- identified_subspace_basis(x, convention = "right", tol = 1e-12)
  expect_s3_class(out, "identified_subspace")
  expect_equal(out$transform_dim, 3L)
  expect_equal(out$numeric_rank, 1L)
  expect_equal(dim(out$basis), c(3, 1))
  expect_equal(crossprod(out$basis), matrix(1, 1, 1), tolerance = 1e-10)
})

test_that("project_to_subspace is idempotent", {
  set.seed(1)
  d <- 8
  r <- 3
  basis <- qr.Q(qr(matrix(rnorm(d * r), d, r)))
  x <- rnorm(d)

  p1 <- project_to_subspace(x, basis)
  p2 <- project_to_subspace(p1, basis)
  expect_equal(p2, p1, tolerance = 1e-12)
})

test_that("restrict_operator_to_subspace matches full-space application on the subspace", {
  set.seed(2)
  d <- 7
  r <- 3
  Q <- qr.Q(qr(matrix(rnorm(d * d), d, d)))
  basis <- qr.Q(qr(matrix(rnorm(d * r), d, r)))
  z <- rnorm(r)
  x <- basis %*% z

  Qr <- restrict_operator_to_subspace(Q, basis)
  lhs <- crossprod(basis, Q %*% x)
  rhs <- Qr %*% z
  expect_equal(as.numeric(lhs), as.numeric(rhs), tolerance = 1e-10)

  Q_lift <- lift_operator_from_subspace(Qr, basis, fill = "zero")
  expect_equal(as.numeric(Q_lift %*% x), as.numeric(basis %*% rhs), tolerance = 1e-10)
})

