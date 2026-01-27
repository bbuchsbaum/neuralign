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


# ---------- More procrustes coverage tests ----------

test_that("procrustes fit generates held-out subject transforms", {
  neuralign:::.register_procrustes()

  set.seed(15)
  d <- 4
  n <- 6
  Z <- matrix(rnorm(d * n), d, n)

  Q1 <- qr.Q(qr(matrix(rnorm(d * d), d, d)))
  if (det(Q1) < 0) Q1[, 1] <- -Q1[, 1]
  Q2 <- qr.Q(qr(matrix(rnorm(d * d), d, d)))
  if (det(Q2) < 0) Q2[, 1] <- -Q2[, 1]
  Q3 <- qr.Q(qr(matrix(rnorm(d * d), d, d)))
  if (det(Q3) < 0) Q3[, 1] <- -Q3[, 1]

  data_list <- list(
    s1 = Q1 %*% Z,
    s2 = Q2 %*% Z,
    s3 = Q3 %*% Z
  )
  adat <- AlignmentData(data_list)

  # Fit using train_idx for only 2 subjects, but expect all 3 to get transforms
  res <- fit_alignment(
    adat, method = "procrustes", reference = "s1",
    cv = "none", train_idx = 1:2, compute_quality = FALSE
  )

  model <- get_model(res)
  # All 3 subjects should have transforms (s3 is held-out, fitted via fallback)
  expect_equal(length(model@transforms), 3)
  expect_true("s3" %in% names(model@transforms))

  # Verify held-out transform is reasonable (orthogonal)
  T3 <- model@transforms[["s3"]]
  expect_equal(dim(T3), c(d, d))
  expect_equal(T3 %*% t(T3), diag(d), tolerance = 1e-6)
})

test_that("procrustes fit with template matrix reference", {
  neuralign:::.register_procrustes()

  set.seed(16)
  d <- 4
  n <- 6
  Z <- matrix(rnorm(d * n), d, n)

  Q <- qr.Q(qr(matrix(rnorm(d * d), d, d)))
  if (det(Q) < 0) Q[, 1] <- -Q[, 1]

  data_list <- list(s1 = Q %*% Z, s2 = Z)
  adat <- AlignmentData(data_list)

  res <- fit_alignment(
    adat, method = "procrustes", reference = Z,
    cv = "none", compute_quality = FALSE
  )

  model <- get_model(res)
  expect_equal(length(model@transforms), 2)
  # s2 data == Z == reference, so transform should be identity
  expect_equal(model@transforms[["s2"]], diag(d), tolerance = 1e-6)
})

test_that("procrustes fit with scale=TRUE recovers scale", {
  neuralign:::.register_procrustes()

  set.seed(17)
  d <- 4
  n <- 8
  Z <- matrix(rnorm(d * n), d, n)

  Q <- qr.Q(qr(matrix(rnorm(d * d), d, d)))
  if (det(Q) < 0) Q[, 1] <- -Q[, 1]

  data_list <- list(
    s1 = Z,
    s2 = 2 * (Q %*% Z)
  )
  adat <- AlignmentData(data_list)

  res <- fit_alignment(
    adat, method = "procrustes", reference = "s1",
    cv = "none", compute_quality = FALSE, scale = TRUE
  )

  model <- get_model(res)
  # Scale should be stored in method_state
  expect_true(model@method_state$scale)
})

test_that("GPA builtin with reflection parameter", {
  neuralign:::.register_procrustes()

  set.seed(18)
  d <- 4
  n <- 6
  Z <- matrix(rnorm(d * n), d, n)

  # Create subjects with a reflection
  Q_reflect <- qr.Q(qr(matrix(rnorm(d * d), d, d)))
  if (det(Q_reflect) > 0) Q_reflect[, 1] <- -Q_reflect[, 1]
  expect_lt(det(Q_reflect), 0)

  data_list <- list(
    s1 = Z,
    s2 = Q_reflect %*% Z
  )
  adat <- AlignmentData(data_list)

  # Without reflection
  res_no_ref <- fit_alignment(
    adat, method = "procrustes", reference = "s1",
    cv = "none", compute_quality = FALSE, reflection = FALSE
  )
  T_no_ref <- get_model(res_no_ref)@transforms[["s2"]]
  expect_gte(det(T_no_ref), 0)

  # With reflection
  res_ref <- fit_alignment(
    adat, method = "procrustes", reference = "s1",
    cv = "none", compute_quality = FALSE, reflection = TRUE
  )
  T_ref <- get_model(res_ref)@transforms[["s2"]]
  expect_lt(det(T_ref), 0)
})

test_that("procrustes_rotation scale with zero data returns scale=1", {
  # When source data is all zeros, denom = 0, scale defaults to 1
  X <- matrix(0, 3, 5)
  Y <- matrix(rnorm(15), 3, 5)

  res <- procrustes_rotation(X, Y, "left", scale = TRUE)
  expect_equal(res$scale_factor, 1)
})

