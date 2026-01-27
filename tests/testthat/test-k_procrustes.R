test_that("k_kernel_roots recovers K and produces inverse roots", {
  set.seed(1)
  q <- 6
  A <- matrix(rnorm(q * q), q, q)
  K <- t(A) %*% A + diag(0.1, q)

  roots <- k_kernel_roots(K, jitter = 1e-12)
  expect_true(all(roots$evals > 0))

  Iq <- diag(q)
  expect_equal(roots$Khalf %*% roots$Khalf, K, tolerance = 1e-8)
  expect_equal(roots$Khalf %*% roots$Kihalf, Iq, tolerance = 1e-8)

  roots_I <- k_kernel_roots(Iq)
  expect_equal(roots_I$Khalf, Iq)
  expect_equal(roots_I$Kihalf, Iq)
  expect_equal(roots_I$evals, rep(1, q))
})

test_that("k_orthonormalize produces K-orthonormal columns", {
  set.seed(2)
  q <- 8
  r <- 3
  A <- matrix(rnorm(q * q), q, q)
  K <- t(A) %*% A + diag(0.1, q)
  W <- matrix(rnorm(q * r), q, r)

  U <- k_orthonormalize(W, K)
  expect_equal(t(U) %*% K %*% U, diag(r), tolerance = 1e-8)
  expect_equal(ncol(U), r)

  Ui <- k_orthonormalize(W, diag(q))
  qrw <- qr(W)
  Q <- qr.Q(qrw)
  R <- qr.R(qrw)
  sgn <- sign(diag(R))
  sgn[!is.finite(sgn) | sgn == 0] <- 1
  Q <- sweep(Q, 2, sgn, `*`)
  expect_equal(Ui, Q, tolerance = 1e-10)
})

test_that("k_procrustes recovers known rotation and enforces reflection control", {
  set.seed(3)
  q <- 7
  r <- 4
  A <- matrix(rnorm(q * q), q, q)
  K <- t(A) %*% A + diag(0.1, q)

  Uref <- k_orthonormalize(matrix(rnorm(q * r), q, r), K)
  Rtrue <- qr.Q(qr(matrix(rnorm(r * r), r, r)))
  if (det(Rtrue) < 0) Rtrue[, r] <- -Rtrue[, r]

  U <- Uref %*% Rtrue
  pr <- k_procrustes(Uref, U, K, allow_reflection = TRUE)

  expect_equal(pr$U_aligned, Uref, tolerance = 1e-8)
  expect_equal(pr$R, t(Rtrue), tolerance = 1e-8)
  expect_equal(t(pr$R) %*% pr$R, diag(r), tolerance = 1e-8)
  expect_true(all(pr$cosines >= 0 & pr$cosines <= 1))

  # Reflection control
  Rneg <- Rtrue
  Rneg[, 1] <- -Rneg[, 1]
  expect_lt(det(Rneg), 0)
  Uneg <- Uref %*% Rneg

  pr_reflect <- k_procrustes(Uref, Uneg, K, allow_reflection = TRUE)
  expect_lt(pr_reflect$determinant, 0)
  expect_equal(pr_reflect$U_aligned, Uref, tolerance = 1e-8)

  pr_no_reflect <- k_procrustes(Uref, Uneg, K, allow_reflection = FALSE)
  expect_gte(pr_no_reflect$determinant, 0)
})

test_that("k_align_bases aligns identical bases to identity", {
  set.seed(4)
  q <- 6
  r <- 3
  K <- diag(q)
  U <- k_orthonormalize(matrix(rnorm(q * r), q, r), K)

  res <- k_align_bases(list(U, U), K, ref = 1L)
  expect_equal(res$R[[1]], diag(r), tolerance = 1e-10)
  expect_equal(res$R[[2]], diag(r), tolerance = 1e-10)
  expect_equal(res$U_aligned[[1]], U, tolerance = 1e-10)
  expect_equal(res$U_aligned[[2]], U, tolerance = 1e-10)
})

test_that("k_consensus_basis returns a K-orthonormal consensus", {
  set.seed(5)
  q <- 6
  r <- 3
  A <- matrix(rnorm(q * q), q, q)
  K <- t(A) %*% A + diag(0.1, q)
  U <- k_orthonormalize(matrix(rnorm(q * r), q, r), K)

  res <- k_consensus_basis(list(U, U, U), K, max_iter = 10, tol = 1e-10)
  expect_true(isTRUE(res$converged))
  expect_equal(t(res$U) %*% K %*% res$U, diag(r), tolerance = 1e-8)
  expect_equal(res$U, U, tolerance = 1e-8)
})
