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


# ---------- Additional procrustes primitive tests ----------

test_that("procrustes_rotation errors on non-matrix input", {
  expect_error(procrustes_rotation(1:10, matrix(1, 2, 5), "left"), "matrix-like")
  expect_error(procrustes_rotation(matrix(1, 2, 5), "text", "left"), "matrix-like")
})

test_that("procrustes_rotation errors on invalid min_overlap", {
  X <- matrix(1, 3, 5)
  Y <- matrix(1, 3, 5)
  expect_error(procrustes_rotation(X, Y, "left", min_overlap = 0L), "positive integer")
  expect_error(procrustes_rotation(X, Y, "left", min_overlap = -1L), "positive integer")
})

test_that("procrustes_rotation right convention errors on dimension mismatch", {
  X <- matrix(rnorm(12), 4, 3)  # 4 obs x 3 features
  Y <- matrix(rnorm(15), 5, 3)  # 5 obs x 3 features - different observations
  expect_error(procrustes_rotation(X, Y, "right"), "matching observations")

  X2 <- matrix(rnorm(12), 3, 4)  # 3 obs x 4 features
  Y2 <- matrix(rnorm(9), 3, 3)   # 3 obs x 3 features - different features
  expect_error(procrustes_rotation(X2, Y2, "right"), "matching feature dimensions")
})

test_that("procrustes_rotation left convention errors on feature mismatch", {
  X <- matrix(rnorm(10), 2, 5)  # 2 features x 5 obs
  Y <- matrix(rnorm(15), 3, 5)  # 3 features x 5 obs
  expect_error(procrustes_rotation(X, Y, "left"), "matching feature dimensions")
})

test_that("procrustes_rotation right convention recovers rotation", {
  set.seed(10)
  d <- 5
  n <- 12
  Q_true <- qr.Q(qr(matrix(rnorm(d * d), d, d)))
  if (det(Q_true) < 0) Q_true[, d] <- -Q_true[, d]

  Xr <- matrix(rnorm(n * d), n, d)
  Yr <- Xr %*% Q_true

  res <- procrustes_rotation(Xr, Yr, "right", scale = TRUE)
  expect_equal(res$convention, "right")
  expect_equal(res$scale_factor, 1, tolerance = 1e-6)
  expect_equal(res$residual, 0, tolerance = 1e-8)
})

test_that("procrustes_rotation with scale=TRUE right convention", {
  set.seed(11)
  d <- 4
  n <- 10
  Q_true <- qr.Q(qr(matrix(rnorm(d * d), d, d)))
  if (det(Q_true) < 0) Q_true[, d] <- -Q_true[, d]

  Xr <- matrix(rnorm(n * d), n, d)
  Yr <- 3 * (Xr %*% Q_true)

  res <- procrustes_rotation(Xr, Yr, "right", scale = TRUE)
  expect_equal(res$scale_factor, 3, tolerance = 1e-6)
  expect_equal(res$residual, 0, tolerance = 1e-6)
})

test_that("procrustes_rotation with label matching in right convention", {
  set.seed(12)
  d <- 4
  Q_true <- qr.Q(qr(matrix(rnorm(d * d), d, d)))
  if (det(Q_true) < 0) Q_true[, d] <- -Q_true[, d]

  labs_x <- paste0("obs", 1:8)
  labs_y <- paste0("obs", 4:11)

  Xr <- matrix(rnorm(8 * d), 8, d)
  overlap_x <- Xr[4:8, , drop = FALSE]
  Yr <- matrix(rnorm(8 * d), 8, d)
  Yr[1:5, ] <- overlap_x %*% Q_true

  res <- procrustes_rotation(
    Xr, Yr, "right",
    obs_labels_source = labs_x,
    obs_labels_target = labs_y
  )
  expect_equal(res$matched_labels, paste0("obs", 4:8))
  expect_equal(res$residual, 0, tolerance = 1e-8)
})

test_that("procrustes_distance with right convention", {
  set.seed(13)
  d <- 4
  n <- 8
  Xr <- matrix(rnorm(n * d), n, d)
  Q <- qr.Q(qr(matrix(rnorm(d * d), d, d)))
  if (det(Q) < 0) Q[, d] <- -Q[, d]
  Yr <- Xr %*% Q

  expect_equal(procrustes_distance(Xr, Yr, "right"), 0, tolerance = 1e-8)
  expect_equal(procrustes_distance(Xr, Xr, "right"), 0, tolerance = 1e-10)
})

test_that("GPA builtin converges for multiple subjects", {
  neuralign:::.register_procrustes()

  set.seed(14)
  d <- 5
  n <- 10
  Z <- matrix(rnorm(d * n), d, n)

  data_list <- lapply(1:4, function(i) {
    Q <- qr.Q(qr(matrix(rnorm(d * d), d, d)))
    if (det(Q) < 0) Q[, 1] <- -Q[, 1]
    Q %*% Z
  })
  names(data_list) <- paste0("s", 1:4)

  adat <- AlignmentData(data_list)
  res <- fit_alignment(adat, method = "procrustes", reference = "consensus", compute_quality = FALSE)
  expect_s4_class(res, "AlignmentResult")

  model <- get_model(res)
  expect_equal(length(model@transforms), 4)
})

