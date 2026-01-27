random_rotation <- function(d) {
  q <- qr.Q(qr(matrix(rnorm(d * d), d, d)))
  if (det(q) < 0) q[, 1] <- -q[, 1]
  q
}

.ensure_procrustes_graph_registered <- function() {
  if (!is_aligner_registered("procrustes_graph")) {
    if (!is_aligner_registered("procrustes")) {
      neuralign:::.register_procrustes()
    }
    neuralign:::.register_procrustes_graph()
  }
}

test_that("procrustes_graph aligns with empty global label intersection", {
  set.seed(1)
  .ensure_procrustes_graph_registered()

  d <- 5
  labels <- paste0("stim", 1:10)
  Z <- matrix(rnorm(d * length(labels)), d, length(labels))

  A1 <- random_rotation(d)
  A2 <- random_rotation(d)
  A3 <- random_rotation(d)

  X1 <- A1 %*% Z[, 1:5, drop = FALSE]
  X2 <- A2 %*% Z[, 1:10, drop = FALSE]
  X3 <- A3 %*% Z[, 6:10, drop = FALSE]

  adat <- AlignmentData(
    data = list(s1 = X1, s2 = X2, s3 = X3),
    obs_labels = list(
      s1 = labels[1:5],
      s2 = labels[1:10],
      s3 = labels[6:10]
    )
  )

  res <- fit_alignment(
    adat,
    method = "procrustes_graph",
    reference = "s1",
    min_overlap = d,
    compute_quality = FALSE
  )

  expect_equal(res@model@transforms[["s1"]], diag(d), tolerance = 1e-6)
  expect_equal(res@model@transforms[["s2"]], A1 %*% t(A2), tolerance = 1e-6)
  expect_equal(res@model@transforms[["s3"]], A1 %*% t(A3), tolerance = 1e-6)

  # Overlap check: s2 shares stim1-5 with s1, stim6-10 with s3.
  expect_equal(
    res@aligned[["s2"]][, 1:5, drop = FALSE],
    X1[, 1:5, drop = FALSE],
    tolerance = 1e-6
  )
  expect_equal(
    res@aligned[["s2"]][, 6:10, drop = FALSE],
    res@aligned[["s3"]],
    tolerance = 1e-6
  )
})

test_that("fit_alignment quality metrics handle per-subject obs_labels", {
  set.seed(10)
  .ensure_procrustes_graph_registered()

  d <- 5
  labels <- paste0("stim", 1:10)
  Z <- matrix(rnorm(d * length(labels)), d, length(labels))

  A1 <- random_rotation(d)
  A2 <- random_rotation(d)
  A3 <- random_rotation(d)

  adat <- AlignmentData(
    data = list(
      s1 = A1 %*% Z[, 1:5, drop = FALSE],
      s2 = A2 %*% Z[, 1:10, drop = FALSE],
      s3 = A3 %*% Z[, 6:10, drop = FALSE]
    ),
    obs_labels = list(
      s1 = labels[1:5],
      s2 = labels[1:10],
      s3 = labels[6:10]
    )
  )

  res <- fit_alignment(
    adat,
    method = "procrustes_graph",
    reference = "s1",
    min_overlap = d,
    compute_quality = TRUE
  )

  expect_true("mean_pairwise_correlation" %in% names(res@quality))
  expect_true(is.finite(res@quality$mean_pairwise_correlation))
  expect_gt(res@quality$mean_pairwise_correlation, 0.99)
})

test_that("procrustes_graph errors on disconnected overlap graph", {
  set.seed(2)
  .ensure_procrustes_graph_registered()

  d <- 4
  labels <- paste0("stim", 1:4)
  Z <- matrix(rnorm(d * length(labels)), d, length(labels))

  A1 <- random_rotation(d)
  A2 <- random_rotation(d)
  A3 <- random_rotation(d)

  adat <- AlignmentData(
    data = list(
      s1 = A1 %*% Z[, 1:2, drop = FALSE],
      s2 = A2 %*% Z[, 1:2, drop = FALSE],
      s3 = A3 %*% Z[, 3:4, drop = FALSE]
    ),
    obs_labels = list(
      s1 = labels[1:2],
      s2 = labels[1:2],
      s3 = labels[3:4]
    )
  )

  expect_error(
    fit_alignment(
      adat,
      method = "procrustes_graph",
      reference = "s1",
      compute_quality = FALSE
    ),
    "disconnected"
  )
})

test_that("procrustes_graph errors when min_overlap yields no edges", {
  set.seed(3)
  .ensure_procrustes_graph_registered()

  d <- 3
  labels <- paste0("stim", 1:8)
  Z <- matrix(rnorm(d * length(labels)), d, length(labels))
  A1 <- random_rotation(d)
  A2 <- random_rotation(d)
  A3 <- random_rotation(d)

  adat <- AlignmentData(
    data = list(
      s1 = A1 %*% Z[, 1:4, drop = FALSE],
      s2 = A2 %*% Z[, 3:6, drop = FALSE],
      s3 = A3 %*% Z[, 5:8, drop = FALSE]
    ),
    obs_labels = list(
      s1 = labels[1:4],
      s2 = labels[3:6],
      s3 = labels[5:8]
    )
  )

  expect_error(
    fit_alignment(
      adat,
      method = "procrustes_graph",
      reference = "s1",
      min_overlap = 3L,
      compute_quality = FALSE
    ),
    "No edges in overlap graph"
  )
})

test_that("procrustes_graph can fit a new subject via apply_alignment()", {
  set.seed(4)
  .ensure_procrustes_graph_registered()

  d <- 5
  labels <- paste0("stim", 1:10)
  Z <- matrix(rnorm(d * length(labels)), d, length(labels))

  A1 <- random_rotation(d)
  A2 <- random_rotation(d)
  A3 <- random_rotation(d)

  adat <- AlignmentData(
    data = list(
      s1 = A1 %*% Z[, 1:10, drop = FALSE],
      s2 = A2 %*% Z[, 1:10, drop = FALSE]
    ),
    obs_labels = list(
      s1 = labels[1:10],
      s2 = labels[1:10]
    )
  )

  fit <- fit_alignment(
    adat,
    method = "procrustes_graph",
    reference = "s1",
    min_overlap = d,
    compute_quality = FALSE
  )

  new_data <- AlignmentData(
    data = list(s3 = A3 %*% Z[, 6:10, drop = FALSE]),
    obs_labels = list(s3 = labels[6:10])
  )

  applied <- apply_alignment(fit, new_data)
  Q3 <- get_transform(get_model(applied), "s3")
  expect_equal(Q3, A1 %*% t(A3), tolerance = 1e-6)
})

test_that("procrustes_graph apply_alignment errors when overlap is too small", {
  set.seed(5)
  .ensure_procrustes_graph_registered()

  d <- 4
  labels <- paste0("stim", 1:8)
  Z <- matrix(rnorm(d * length(labels)), d, length(labels))

  A1 <- random_rotation(d)
  A2 <- random_rotation(d)

  adat <- AlignmentData(
    data = list(
      s1 = A1 %*% Z[, 1:8, drop = FALSE],
      s2 = A2 %*% Z[, 1:8, drop = FALSE]
    ),
    obs_labels = list(
      s1 = labels[1:8],
      s2 = labels[1:8]
    )
  )

  fit <- fit_alignment(
    adat,
    method = "procrustes_graph",
    reference = "s1",
    min_overlap = d,
    compute_quality = FALSE
  )

  new_data <- AlignmentData(
    data = list(s3 = A2 %*% Z[, 1:3, drop = FALSE]),
    obs_labels = list(s3 = labels[1:3])
  )

  expect_error(apply_alignment(fit, new_data), "Not enough shared observation labels")
})
