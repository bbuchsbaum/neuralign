random_rotation <- function(d) {
  q <- qr.Q(qr(matrix(rnorm(d * d), d, d)))
  if (det(q) < 0) q[, 1] <- -q[, 1]
  q
}

.ensure_procrustes_graph_registered <- function() {
  if (!is_aligner_registered("procrustes_graph")) {
    if (!is_aligner_registered("procrustes")) {
      ensure_test_aligner("procrustes")
    }
    ensure_test_aligner("procrustes_graph")
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

test_that("procrustes_graph can fit a new subject without reference overlap via apply_alignment()", {
  set.seed(4)
  .ensure_procrustes_graph_registered()

  d <- 5
  labels <- paste0("stim", 1:10)
  Z <- matrix(rnorm(d * length(labels)), d, length(labels))

  A1 <- random_rotation(d)
  A2 <- random_rotation(d)
  A3 <- random_rotation(d)

  # Reference subject only sees stim1-5. Subject 2 provides stim6-10 to the
  # template (via union-fill), so a new subject that sees only stim6-10 can
  # still be aligned even with zero overlap to the reference.
  adat <- AlignmentData(
    data = list(
      s1 = A1 %*% Z[, 1:5, drop = FALSE],
      s2 = A2 %*% Z[, 1:10, drop = FALSE]
    ),
    obs_labels = list(
      s1 = labels[1:5],
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

test_that("procrustes_graph can fit held-out subjects with train_idx using union template", {
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
    train_idx = 1:2, # fit on s1+s2 only
    min_overlap = d,
    compute_quality = FALSE
  )

  expect_equal(res@model@transforms[["s3"]], A1 %*% t(A3), tolerance = 1e-6)
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


# ---------- Additional procrustes graph tests ----------

test_that("procrustes_graph errors without obs_labels", {
  .ensure_procrustes_graph_registered()

  set.seed(20)
  d <- 3
  adat <- AlignmentData(list(
    s1 = matrix(rnorm(d * 4), d, 4),
    s2 = matrix(rnorm(d * 4), d, 4)
  ))

  expect_error(
    fit_alignment(adat, method = "procrustes_graph", reference = "s1", compute_quality = FALSE),
    "obs_labels"
  )
})

test_that("procrustes_graph with uniform weight works", {
  .ensure_procrustes_graph_registered()

  set.seed(21)
  d <- 4
  labels <- paste0("stim", 1:6)
  Z <- matrix(rnorm(d * length(labels)), d, length(labels))

  A1 <- random_rotation(d)
  A2 <- random_rotation(d)

  adat <- AlignmentData(
    data = list(
      s1 = A1 %*% Z,
      s2 = A2 %*% Z
    ),
    obs_labels = list(
      s1 = labels,
      s2 = labels
    )
  )

  res <- fit_alignment(
    adat,
    method = "procrustes_graph",
    reference = "s1",
    weight = "uniform",
    compute_quality = FALSE
  )

  expect_s4_class(res, "AlignmentResult")
  expect_equal(get_model(res)@method_state$weight, "uniform")
})

test_that("procrustes_graph with shared atomic obs_labels", {
  .ensure_procrustes_graph_registered()

  set.seed(22)
  d <- 3
  n <- 5
  labels <- paste0("s", 1:n)
  Z <- matrix(rnorm(d * n), d, n)

  A1 <- random_rotation(d)
  A2 <- random_rotation(d)

  adat <- AlignmentData(
    data = list(
      s1 = A1 %*% Z,
      s2 = A2 %*% Z
    ),
    obs_labels = labels  # shared atomic labels
  )

  res <- fit_alignment(
    adat,
    method = "procrustes_graph",
    reference = "s1",
    compute_quality = FALSE
  )

  expect_s4_class(res, "AlignmentResult")
  expect_equal(res@model@transforms[["s1"]], diag(d), tolerance = 1e-6)
})

test_that("procrustes_graph subject-axis CV with held-out subject", {
  .ensure_procrustes_graph_registered()

  set.seed(23)
  d <- 4
  n <- 8
  labels <- paste0("stim", 1:n)
  Z <- matrix(rnorm(d * n), d, n)

  A1 <- random_rotation(d)
  A2 <- random_rotation(d)
  A3 <- random_rotation(d)

  adat <- AlignmentData(
    data = list(
      s1 = A1 %*% Z,
      s2 = A2 %*% Z,
      s3 = A3 %*% Z
    ),
    obs_labels = list(
      s1 = labels,
      s2 = labels,
      s3 = labels
    )
  )

  res <- fit_alignment(
    adat,
    method = "procrustes_graph",
    reference = "medoid",
    cv = "loso",
    compute_quality = FALSE
  )

  expect_s4_class(res, "AlignmentResult")
})

test_that("procrustes_graph apply errors on non-AlignmentData input", {
  .ensure_procrustes_graph_registered()

  set.seed(24)
  d <- 3
  Z <- matrix(rnorm(d * 4), d, 4)
  labels <- paste0("s", 1:4)

  adat <- AlignmentData(
    data = list(s1 = Z, s2 = Z),
    obs_labels = list(s1 = labels, s2 = labels)
  )
  fit <- fit_alignment(adat, method = "procrustes_graph", reference = "s1", compute_quality = FALSE)

  expect_error(
    neuralign:::.procrustes_graph_apply(
      list(
        method_state = list(reference = "s1"),
        reference_data = Z
      ),
      "not_adat"
    ),
    "AlignmentData"
  )
})

test_that("procrustes_graph apply errors on multiple new subjects", {
  .ensure_procrustes_graph_registered()

  set.seed(25)
  d <- 3
  Z <- matrix(rnorm(d * 4), d, 4)
  labels <- paste0("s", 1:4)

  adat <- AlignmentData(
    data = list(s1 = Z, s2 = Z),
    obs_labels = list(s1 = labels, s2 = labels)
  )
  fit <- fit_alignment(adat, method = "procrustes_graph", reference = "s1", compute_quality = FALSE)

  new_data <- AlignmentData(
    data = list(s3 = Z, s4 = Z),
    obs_labels = list(s3 = labels, s4 = labels)
  )

  expect_error(
    neuralign:::.procrustes_graph_apply(
      list(
        method_state = get_model(fit)@method_state,
        reference_data = get_model(fit)@reference_data
      ),
      new_data
    ),
    "single new subject"
  )
})

test_that("procrustes_graph respects reflection=FALSE in synchronization", {
  .ensure_procrustes_graph_registered()

  set.seed(26)
  d <- 4
  n <- 6
  labels <- paste0("stim", seq_len(n))
  Z <- matrix(rnorm(d * n), d, n)
  R_reflect <- diag(c(-1, rep(1, d - 1L)))

  adat <- AlignmentData(
    data = list(
      s1 = Z,
      s2 = R_reflect %*% Z
    ),
    obs_labels = labels
  )

  no_reflect <- fit_alignment(
    adat,
    method = "procrustes_graph",
    reference = "s1",
    min_overlap = d,
    reflection = FALSE,
    compute_quality = FALSE
  )
  Q2 <- get_transform(get_model(no_reflect), "s2")
  expect_gt(det(Q2), 0)
  expect_gt(norm(no_reflect@aligned[["s2"]] - Z, "F"), 1e-6)

  allow_reflect <- fit_alignment(
    adat,
    method = "procrustes_graph",
    reference = "s1",
    min_overlap = d,
    reflection = TRUE,
    compute_quality = FALSE
  )
  Q2r <- get_transform(get_model(allow_reflect), "s2")
  expect_lt(det(Q2r), 0)
  expect_true(max(abs(allow_reflect@aligned[["s2"]] - Z)) < 1e-8)
})
