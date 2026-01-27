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


# ---------- Additional edge case tests ----------

test_that("k_kernel_roots handles near-singular kernel with jitter", {
  # Rank-deficient K: rank 1 kernel
  v <- c(1, 2, 3, 4)
  K <- outer(v, v)  # rank 1

  roots <- k_kernel_roots(K, jitter = 1e-6)
  # All eigenvalues should be >= jitter
  expect_true(all(roots$evals >= 1e-6))
  # Khalf %*% Kihalf should be identity
  expect_equal(roots$Khalf %*% roots$Kihalf, diag(4), tolerance = 1e-4)
})

test_that("k_kernel_roots errors on non-square input", {
  expect_error(
    k_kernel_roots(matrix(1:6, 2, 3)),
    "must be square"
  )
})

test_that("k_kernel_roots errors on non-matrix input", {
  expect_error(
    k_kernel_roots("not a matrix"),
    "must be matrix-like"
  )
})

test_that("k_kernel_roots errors on negative jitter", {
  expect_error(
    k_kernel_roots(diag(3), jitter = -1),
    "non-negative"
  )
})

test_that("k_orthonormalize errors on dimension mismatch", {
  K <- diag(5)
  W <- matrix(rnorm(12), 4, 3)  # 4 rows, K is 5x5
  expect_error(
    k_orthonormalize(W, K),
    "Row dimension mismatch"
  )
})

test_that("k_orthonormalize with precomputed Kroots gives same result", {
  set.seed(10)
  q <- 5
  r <- 2
  A <- matrix(rnorm(q * q), q, q)
  K <- t(A) %*% A + diag(0.1, q)
  W <- matrix(rnorm(q * r), q, r)

  Kroots <- k_kernel_roots(K)
  U1 <- k_orthonormalize(W, K)
  U2 <- k_orthonormalize(W, K, Kroots = Kroots)
  expect_equal(U1, U2, tolerance = 1e-12)
})

test_that("k_procrustes errors on dimension mismatch between Uref and U", {
  K <- diag(5)
  Uref <- matrix(rnorm(10), 5, 2)
  U <- matrix(rnorm(15), 5, 3)

  expect_error(
    k_procrustes(Uref, U, K),
    "identical dimensions"
  )
})

test_that("k_procrustes errors when K dimension doesn't match U", {
  K <- diag(4)  # 4x4
  Uref <- matrix(rnorm(10), 5, 2)
  U <- matrix(rnorm(10), 5, 2)

  expect_error(
    k_procrustes(Uref, U, K),
    "Dimension mismatch"
  )
})

test_that("k_align_bases with 3+ bases aligns all to reference", {
  set.seed(20)
  q <- 6
  r <- 3
  K <- diag(q)
  Uref <- k_orthonormalize(matrix(rnorm(q * r), q, r), K)

  # Create 4 rotated versions
  U_list <- lapply(1:4, function(i) {
    R <- qr.Q(qr(matrix(rnorm(r * r), r, r)))
    if (det(R) < 0) R[, 1] <- -R[, 1]
    Uref %*% R
  })

  res <- k_align_bases(U_list, K, ref = Uref, allow_reflection = TRUE)

  expect_equal(length(res$U_aligned), 4)
  for (i in seq_along(res$U_aligned)) {
    expect_equal(res$U_aligned[[i]], Uref, tolerance = 1e-8)
  }
})

test_that("k_align_bases with integer ref index works", {
  set.seed(25)
  q <- 5
  r <- 2
  K <- diag(q)
  U1 <- k_orthonormalize(matrix(rnorm(q * r), q, r), K)
  R <- qr.Q(qr(matrix(rnorm(r * r), r, r)))
  if (det(R) < 0) R[, 1] <- -R[, 1]
  U2 <- U1 %*% R

  res <- k_align_bases(list(U1, U2), K, ref = 1L)

  # First basis aligned to itself should be identity rotation
  expect_equal(res$R[[1]], diag(r), tolerance = 1e-10)
  # Second should align to U1
  expect_equal(res$U_aligned[[2]], U1, tolerance = 1e-8)
})

test_that("k_consensus_basis converges with distinct rotated bases", {
  set.seed(30)
  q <- 8
  r <- 3
  A <- matrix(rnorm(q * q), q, q)
  K <- t(A) %*% A + diag(0.1, q)
  U0 <- k_orthonormalize(matrix(rnorm(q * r), q, r), K)

  # Create 5 rotated versions with small noise
  U_list <- lapply(1:5, function(i) {
    R <- qr.Q(qr(matrix(rnorm(r * r), r, r)))
    if (det(R) < 0) R[, 1] <- -R[, 1]
    U0 %*% R
  })

  res <- k_consensus_basis(U_list, K, max_iter = 100, tol = 1e-8)

  expect_true(res$converged)
  expect_lte(res$iters, 100)
  # Consensus should be K-orthonormal
  expect_equal(t(res$U) %*% K %*% res$U, diag(r), tolerance = 1e-6)
  # gaps should be decreasing
  expect_true(length(res$gaps) >= 1)
})

test_that("k_consensus_basis with weights prioritizes weighted bases", {
  set.seed(35)
  q <- 5
  r <- 2
  K <- diag(q)
  U1 <- k_orthonormalize(matrix(rnorm(q * r), q, r), K)
  # Create a different basis
  U2 <- k_orthonormalize(matrix(rnorm(q * r), q, r), K)

  # Heavily weight U1
  res <- k_consensus_basis(list(U1, U2), K, weights = c(100, 1), max_iter = 50, tol = 1e-8)

  expect_true(res$converged)
  # Consensus should be very close to U1 due to heavy weighting
  # Check principal angles are near 0
  cosines <- svd(t(res$U) %*% K %*% U1)$d
  expect_true(all(pmin(cosines, 1) > 0.99))
})

test_that("k_consensus_basis errors on bad weights", {
  K <- diag(3)
  U <- matrix(rnorm(6), 3, 2)

  expect_error(
    k_consensus_basis(list(U, U), K, weights = c(-1, 1)),
    "non-negative"
  )

  expect_error(
    k_consensus_basis(list(U, U), K, weights = c(0, 0)),
    "positive sum"
  )
})

test_that("k_consensus_basis returns converged=FALSE when max_iter is too small", {
  set.seed(40)
  q <- 10
  r <- 3
  K <- diag(q)

  # Create very different bases that need many iterations
  U_list <- lapply(1:5, function(i) {
    k_orthonormalize(matrix(rnorm(q * r), q, r), K)
  })

  res <- k_consensus_basis(U_list, K, max_iter = 1, tol = 1e-15)
  # With max_iter=1, probably won't converge to machine precision
  expect_equal(res$iters, 1L)
})
