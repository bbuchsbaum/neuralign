test_that("procrustes_rotation recovers known rotation (left)", {
  set.seed(1)
  d <- 10
  n <- 25
  Q_true <- qr.Q(qr(matrix(rnorm(d * d), d, d)))
  if (det(Q_true) < 0) Q_true[, d] <- -Q_true[, d]

  X <- matrix(rnorm(d * n), d, n)
  Y <- Q_true %*% X

  res <- procrustes_rotation(X, Y, convention = "left")
  expect_equal(res$convention, "left")
  expect_equal(res$residual, 0, tolerance = 1e-8)
  expect_equal(res$Q, Q_true, tolerance = 1e-8)
})

test_that("procrustes_rotation recovers known rotation (right)", {
  set.seed(2)
  d <- 8
  n <- 20
  Q_true <- qr.Q(qr(matrix(rnorm(d * d), d, d)))
  if (det(Q_true) < 0) Q_true[, d] <- -Q_true[, d]

  X <- matrix(rnorm(d * n), d, n)
  Xr <- t(X)
  Yr <- Xr %*% Q_true

  res <- procrustes_rotation(Xr, Yr, convention = "right")
  expect_equal(res$convention, "right")
  expect_equal(res$residual, 0, tolerance = 1e-8)
  expect_equal(res$Q, Q_true, tolerance = 1e-8)
})

test_that("left/right convention duality holds", {
  set.seed(3)
  d <- 6
  n <- 15
  Q_true <- qr.Q(qr(matrix(rnorm(d * d), d, d)))
  if (det(Q_true) < 0) Q_true[, d] <- -Q_true[, d]

  X <- matrix(rnorm(d * n), d, n)
  Y <- Q_true %*% X

  Q_left <- procrustes_rotation(X, Y, "left")$Q
  Q_right <- procrustes_rotation(t(X), t(Y), "right")$Q
  expect_equal(Q_left, t(Q_right), tolerance = 1e-10)
})

test_that("scale recovery returns expected factor", {
  set.seed(4)
  d <- 7
  n <- 21
  Q_true <- qr.Q(qr(matrix(rnorm(d * d), d, d)))
  if (det(Q_true) < 0) Q_true[, d] <- -Q_true[, d]

  X <- matrix(rnorm(d * n), d, n)
  Y <- 2 * (Q_true %*% X)

  res <- procrustes_rotation(X, Y, "left", scale = TRUE)
  expect_equal(res$scale_factor, 2, tolerance = 1e-8)
  expect_equal(res$Q / res$scale_factor, Q_true, tolerance = 1e-8)
  expect_equal(res$residual, 0, tolerance = 1e-8)
})

test_that("reflection control enforces det(Q) >= 0 when reflection=FALSE", {
  set.seed(5)
  d <- 9
  n <- 30
  Q_true <- qr.Q(qr(matrix(rnorm(d * d), d, d)))
  if (det(Q_true) > 0) Q_true[, 1] <- -Q_true[, 1] # force det < 0
  expect_lt(det(Q_true), 0)

  X <- matrix(rnorm(d * n), d, n)
  Y <- Q_true %*% X

  res_no_reflect <- procrustes_rotation(X, Y, "left", reflection = FALSE)
  expect_gte(det(res_no_reflect$Q), 0)

  res_reflect <- procrustes_rotation(X, Y, "left", reflection = TRUE)
  expect_lt(det(res_reflect$Q), 0)
  expect_equal(res_reflect$residual, 0, tolerance = 1e-8)
})

test_that("procrustes_distance basic properties hold", {
  set.seed(6)
  d <- 5
  n <- 12
  X <- matrix(rnorm(d * n), d, n)
  Q_true <- qr.Q(qr(matrix(rnorm(d * d), d, d)))
  if (det(Q_true) < 0) Q_true[, d] <- -Q_true[, d]
  Y <- Q_true %*% X

  expect_equal(procrustes_distance(X, X, "left"), 0, tolerance = 1e-10)
  expect_equal(procrustes_distance(X, Y, "left"), 0, tolerance = 1e-8)
  expect_equal(procrustes_distance(X, Y, "left"), procrustes_distance(Y, X, "left"), tolerance = 1e-10)
})

test_that("procrustes_rotation supports label intersection", {
  set.seed(7)
  d <- 6
  Q_true <- qr.Q(qr(matrix(rnorm(d * d), d, d)))
  if (det(Q_true) < 0) Q_true[, d] <- -Q_true[, d]

  labs_x <- paste0("img-", 1:10)
  labs_y <- paste0("img-", 6:15)

  X <- matrix(rnorm(d * length(labs_x)), d, length(labs_x))
  # Construct Y so that the overlapping part is exactly aligned
  X_overlap <- X[, 6:10, drop = FALSE]
  Y_overlap <- Q_true %*% X_overlap
  Y <- matrix(rnorm(d * length(labs_y)), d, length(labs_y))
  Y[, 1:5] <- Y_overlap

  res <- procrustes_rotation(
    X, Y,
    convention = "left",
    obs_labels_source = labs_x,
    obs_labels_target = labs_y
  )
  expect_equal(res$matched_labels, labs_x[6:10])
  expect_equal(res$residual, 0, tolerance = 1e-8)
})

test_that("procrustes_rotation errors on incompatible dimensions", {
  X <- matrix(rnorm(10), 2, 5)
  Y <- matrix(rnorm(12), 3, 4)
  expect_error(procrustes_rotation(X, Y, "left"))
})

