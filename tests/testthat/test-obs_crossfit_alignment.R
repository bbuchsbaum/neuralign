test_that("run_obs_crossfit_alignment fits per fold and aligns held-out data", {
  neuralign:::.register_procrustes()

  set.seed(1)
  d <- 4
  n_obs <- 6

  Z <- matrix(rnorm(d * n_obs), d, n_obs)
  Q <- qr.Q(qr(matrix(rnorm(d * d), d, d)))
  if (det(Q) < 0) Q[, 1] <- -Q[, 1]

  X1 <- Z
  X2 <- Q %*% Z

  folds <- list(
    f1 = list(train = 1:3, test = 4:6),
    f2 = list(train = 4:6, test = 1:3)
  )

  train_data <- lapply(folds, function(f) {
    list(
      s1 = X1[, f$train, drop = FALSE],
      s2 = X2[, f$train, drop = FALSE]
    )
  })
  test_data <- lapply(folds, function(f) {
    list(
      s1 = X1[, f$test, drop = FALSE],
      s2 = X2[, f$test, drop = FALSE]
    )
  })

  res <- run_obs_crossfit_alignment(
    train_data_by_fold = train_data,
    test_data_by_fold = test_data,
    obs_labels_train = paste0("tr", 1:3),
    obs_labels_test = paste0("te", 1:3),
    method = "procrustes",
    reference = "s1",
    compute_quality = FALSE
  )

  expect_s3_class(res, "ObsCrossfitAlignment")
  expect_true(isTRUE(res$anchor_common))
  expect_equal(res$reference_kind, "fixed_subject")
  expect_equal(names(res$models_by_fold), c("f1", "f2"))
  expect_equal(res$models_by_fold$f1@method, "procrustes")

  for (fid in names(folds)) {
    aligned <- res$aligned_test_by_fold[[fid]]
    expect_true(max(abs(aligned$s1 - aligned$s2)) < 1e-8)

    tf <- res$transforms_by_fold[[fid]]
    expect_equal(dim(tf$s2), c(d, d))
  }
})

test_that("run_obs_crossfit_alignment enforces anchor_policy for data-driven references", {
  neuralign:::.register_procrustes()

  set.seed(2)
  d <- 3
  labels <- paste0("t", 1:6)

  Z <- matrix(rnorm(d * length(labels)), d, length(labels))
  train_data <- list(
    f1 = list(s1 = Z[, 1:3, drop = FALSE], s2 = Z[, 1:3, drop = FALSE]),
    f2 = list(s1 = Z[, 4:6, drop = FALSE], s2 = Z[, 4:6, drop = FALSE])
  )

  obs_labels_train <- list(
    f1 = labels[1:3],
    f2 = labels[4:6]
  )

  expect_error(
    run_obs_crossfit_alignment(
      train_data_by_fold = train_data,
      method = "procrustes",
      reference = "consensus",
      obs_labels_train = obs_labels_train,
      anchor_policy = "common_or_error"
    ),
    "common_or_error"
  )

  ok <- run_obs_crossfit_alignment(
    train_data_by_fold = train_data,
    method = "procrustes",
    reference = "consensus",
    obs_labels_train = obs_labels_train,
    anchor_policy = "fold_specific_ok"
  )
  expect_false(isTRUE(ok$anchor_common))
  expect_equal(ok$reference_kind, "data_driven")

  template <- Z
  colnames(template) <- labels
  mapped <- run_obs_crossfit_alignment(
    train_data_by_fold = train_data,
    method = "procrustes",
    reference = "consensus",
    obs_labels_train = obs_labels_train,
    anchor_policy = "map_to_template",
    template = template,
    template_obs_labels = labels,
    min_overlap = 2
  )
  expect_true(isTRUE(mapped$anchor_common))
  expect_equal(mapped$provenance$anchor_policy, "map_to_template")
})

test_that("run_obs_crossfit_alignment accepts per-subject obs_labels lists", {
  neuralign:::.register_procrustes()

  set.seed(3)
  d <- 3
  Z <- matrix(rnorm(d * 4), d, 4)

  train_data <- list(
    f1 = list(s1 = Z[, 1:2, drop = FALSE], s2 = Z[, 1:2, drop = FALSE]),
    f2 = list(s1 = Z[, 3:4, drop = FALSE], s2 = Z[, 3:4, drop = FALSE])
  )

  labs_by_subj <- list(
    s1 = c("a", "b"),
    s2 = c("a", "b")
  )

  res <- run_obs_crossfit_alignment(
    train_data_by_fold = train_data,
    method = "procrustes",
    reference = "s1",
    obs_labels_train = labs_by_subj,
    anchor_policy = "common_or_error"
  )
  expect_true(isTRUE(res$anchor_common))
})
