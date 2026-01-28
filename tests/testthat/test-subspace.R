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

test_that("canonicalize_orthogonal_operator preserves action on identified subspace and removes complement ambiguity", {
  set.seed(3)
  d <- 8
  r <- 3

  # Rank-deficient correspondence matrix with column space span(B_source)
  B_source <- qr.Q(qr(matrix(rnorm(d * r), d, r)))
  X <- B_source %*% matrix(rnorm(r * 5), r, 5)

  # Build a full basis that starts with B_source
  Z <- matrix(rnorm(d * (d - r)), d, d - r)
  Z <- Z - B_source %*% crossprod(B_source, Z)
  B_perp <- qr.Q(qr(Z))
  B_full <- cbind(B_source, B_perp)

  Q_base <- qr.Q(qr(matrix(rnorm(d * d), d, d)))

  make_Q <- function() {
    R_perp <- qr.Q(qr(matrix(rnorm((d - r) * (d - r)), d - r, d - r)))
    Rcomp <- diag(d)
    Rcomp[(r + 1):d, (r + 1):d] <- R_perp
    R_full <- B_full %*% Rcomp %*% t(B_full)
    Q_base %*% R_full
  }

  Q1 <- make_Q()
  Q2 <- make_Q()

  # Q1 and Q2 agree on span(X) but differ on the orthogonal complement
  expect_equal(as.numeric(Q1 %*% X), as.numeric(Q2 %*% X), tolerance = 1e-10)

  Q1c <- canonicalize_orthogonal_operator(Q1, X, convention = "left", tol = 1e-12)
  Q2c <- canonicalize_orthogonal_operator(Q2, X, convention = "left", tol = 1e-12)

  # Canonicalization should remove the complement ambiguity
  expect_equal(Q1c, Q2c, tolerance = 1e-10)

  # And preserve action on the identified subspace
  expect_equal(as.numeric(Q1c %*% X), as.numeric(Q1 %*% X), tolerance = 1e-10)
  expect_equal(crossprod(Q1c), diag(d), tolerance = 1e-8)
})
