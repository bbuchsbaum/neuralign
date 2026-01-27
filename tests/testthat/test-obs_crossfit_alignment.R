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

  expect_warning(
    ok <- run_obs_crossfit_alignment(
      train_data_by_fold = train_data,
      method = "procrustes",
      reference = "consensus",
      obs_labels_train = obs_labels_train,
      anchor_policy = "fold_specific_ok"
    ),
    "Fold-specific anchors"
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


# ---------- Additional obs crossfit tests ----------

test_that("run_obs_crossfit_alignment train-only (no test data) returns NULL aligned_test", {
  neuralign:::.register_procrustes()

  set.seed(10)
  d <- 3
  Z <- matrix(rnorm(d * 6), d, 6)

  train_data <- list(
    f1 = list(s1 = Z[, 1:3, drop = FALSE], s2 = Z[, 1:3, drop = FALSE]),
    f2 = list(s1 = Z[, 4:6, drop = FALSE], s2 = Z[, 4:6, drop = FALSE])
  )

  res <- run_obs_crossfit_alignment(
    train_data_by_fold = train_data,
    test_data_by_fold = NULL,
    method = "procrustes",
    reference = "s1"
  )

  expect_s3_class(res, "ObsCrossfitAlignment")
  expect_null(res$aligned_test_by_fold)
  expect_equal(length(res$transforms_by_fold), 2)
})

test_that("run_obs_crossfit_alignment errors on mismatched fold ids", {
  neuralign:::.register_procrustes()

  set.seed(11)
  d <- 3
  Z <- matrix(rnorm(d * 6), d, 6)

  train_data <- list(
    f1 = list(s1 = Z[, 1:3, drop = FALSE], s2 = Z[, 1:3, drop = FALSE])
  )
  test_data <- list(
    f_wrong = list(s1 = Z[, 4:6, drop = FALSE], s2 = Z[, 4:6, drop = FALSE])
  )

  expect_error(
    run_obs_crossfit_alignment(
      train_data_by_fold = train_data,
      test_data_by_fold = test_data,
      method = "procrustes",
      reference = "s1"
    ),
    "same fold ids"
  )
})

test_that("run_obs_crossfit_alignment errors on mismatched subjects in test data", {
  neuralign:::.register_procrustes()

  set.seed(12)
  d <- 3
  Z <- matrix(rnorm(d * 6), d, 6)

  train_data <- list(
    f1 = list(s1 = Z[, 1:3, drop = FALSE], s2 = Z[, 1:3, drop = FALSE])
  )
  test_data <- list(
    f1 = list(s1 = Z[, 4:6, drop = FALSE], s3 = Z[, 4:6, drop = FALSE])
  )

  expect_error(
    run_obs_crossfit_alignment(
      train_data_by_fold = train_data,
      test_data_by_fold = test_data,
      method = "procrustes",
      reference = "s1"
    ),
    "same subject ids"
  )
})

test_that("run_obs_crossfit_alignment errors on invalid min_overlap", {
  expect_error(
    run_obs_crossfit_alignment(
      train_data_by_fold = list(f1 = list(s1 = matrix(1, 2, 2))),
      method = "procrustes",
      reference = "s1",
      min_overlap = 0L
    ),
    "positive integer"
  )
})

test_that("run_obs_crossfit_alignment errors on unnamed fold list", {
  expect_error(
    run_obs_crossfit_alignment(
      train_data_by_fold = list(list(s1 = matrix(1, 2, 2))),
      method = "procrustes",
      reference = "s1"
    ),
    "named list"
  )
})

test_that("run_obs_crossfit_alignment errors on missing subjects in fold", {
  train_data <- list(
    f1 = list(s1 = matrix(1, 2, 3), s2 = matrix(1, 2, 3)),
    f2 = list(s1 = matrix(1, 2, 3))  # missing s2
  )

  expect_error(
    run_obs_crossfit_alignment(
      train_data_by_fold = train_data,
      method = "procrustes",
      reference = "s1"
    ),
    "missing subjects"
  )
})

test_that("run_obs_crossfit_alignment errors on non-matrix fold entry", {
  train_data <- list(
    f1 = list(s1 = "not_a_matrix", s2 = matrix(1, 2, 3))
  )

  expect_error(
    run_obs_crossfit_alignment(
      train_data_by_fold = train_data,
      method = "procrustes",
      reference = "s1"
    ),
    "matrix-like"
  )
})

test_that("run_obs_crossfit_alignment map_to_template errors without template", {
  neuralign:::.register_procrustes()

  set.seed(13)
  d <- 3
  Z <- matrix(rnorm(d * 4), d, 4)

  train_data <- list(
    f1 = list(s1 = Z[, 1:2, drop = FALSE], s2 = Z[, 1:2, drop = FALSE]),
    f2 = list(s1 = Z[, 3:4, drop = FALSE], s2 = Z[, 3:4, drop = FALSE])
  )

  expect_error(
    run_obs_crossfit_alignment(
      train_data_by_fold = train_data,
      method = "procrustes",
      reference = "consensus",
      anchor_policy = "map_to_template"
    ),
    "template"
  )
})

test_that("run_obs_crossfit_alignment compute_quality returns quality_by_fold", {
  neuralign:::.register_procrustes()

  set.seed(14)
  d <- 4
  n <- 6
  Z <- matrix(rnorm(d * n), d, n)
  Q <- qr.Q(qr(matrix(rnorm(d * d), d, d)))
  if (det(Q) < 0) Q[, 1] <- -Q[, 1]

  train_data <- list(
    f1 = list(s1 = Z[, 1:3, drop = FALSE], s2 = (Q %*% Z)[, 1:3, drop = FALSE]),
    f2 = list(s1 = Z[, 4:6, drop = FALSE], s2 = (Q %*% Z)[, 4:6, drop = FALSE])
  )

  res <- run_obs_crossfit_alignment(
    train_data_by_fold = train_data,
    method = "procrustes",
    reference = "s1",
    compute_quality = TRUE
  )

  expect_true(!is.null(res$quality_by_fold))
  expect_equal(length(res$quality_by_fold), 2)
})

test_that("run_obs_crossfit_alignment with per-fold obs_labels", {
  neuralign:::.register_procrustes()

  set.seed(15)
  d <- 3
  Z <- matrix(rnorm(d * 6), d, 6)

  train_data <- list(
    f1 = list(s1 = Z[, 1:3, drop = FALSE], s2 = Z[, 1:3, drop = FALSE]),
    f2 = list(s1 = Z[, 4:6, drop = FALSE], s2 = Z[, 4:6, drop = FALSE])
  )

  per_fold_labels <- list(
    f1 = c("a", "b", "c"),
    f2 = c("d", "e", "f")
  )

  res <- run_obs_crossfit_alignment(
    train_data_by_fold = train_data,
    method = "procrustes",
    reference = "s1",
    obs_labels_train = per_fold_labels
  )

  expect_s3_class(res, "ObsCrossfitAlignment")
  expect_equal(length(res$models_by_fold), 2)
})

test_that("run_obs_crossfit_alignment errors on inconsistent feature dims across subjects in fold", {
  train_data <- list(
    f1 = list(s1 = matrix(1, 3, 2), s2 = matrix(1, 4, 2))  # s1=3 rows, s2=4 rows
  )

  expect_error(
    run_obs_crossfit_alignment(
      train_data_by_fold = train_data,
      method = "procrustes",
      reference = "s1"
    ),
    "inconsistent feature dimensions"
  )
})

test_that("run_obs_crossfit_alignment fold_specific_ok with data-driven ref generates warning", {
  neuralign:::.register_procrustes()

  set.seed(16)
  d <- 3
  Z <- matrix(rnorm(d * 6), d, 6)

  train_data <- list(
    f1 = list(s1 = Z[, 1:3, drop = FALSE], s2 = Z[, 1:3, drop = FALSE]),
    f2 = list(s1 = Z[, 4:6, drop = FALSE], s2 = Z[, 4:6, drop = FALSE])
  )

  expect_warning(
    res <- run_obs_crossfit_alignment(
      train_data_by_fold = train_data,
      method = "procrustes",
      reference = "consensus",
      anchor_policy = "fold_specific_ok"
    ),
    "Fold-specific anchors"
  )

  expect_false(isTRUE(res$anchor_common))
  expect_true(length(res$warnings) > 0)
  expect_true(any(grepl("Fold-specific", res$warnings)))
})

test_that(".normalize_obs_labels_for_folds handles factor input", {
  neuralign:::.register_procrustes()

  set.seed(17)
  d <- 3
  Z <- matrix(rnorm(d * 4), d, 4)

  train_data <- list(
    f1 = list(s1 = Z[, 1:2, drop = FALSE], s2 = Z[, 1:2, drop = FALSE]),
    f2 = list(s1 = Z[, 3:4, drop = FALSE], s2 = Z[, 3:4, drop = FALSE])
  )

  # Factor labels
  labs <- factor(c("a", "b"))

  res <- run_obs_crossfit_alignment(
    train_data_by_fold = train_data,
    method = "procrustes",
    reference = "s1",
    obs_labels_train = labs
  )

  expect_s3_class(res, "ObsCrossfitAlignment")
})

test_that(".normalize_obs_labels_for_folds errors on mismatched names", {
  train_data <- list(
    f1 = list(s1 = matrix(1, 2, 3), s2 = matrix(1, 2, 3))
  )

  bad_labels <- list(wrong_key = c("a", "b", "c"))

  expect_error(
    neuralign:::.normalize_obs_labels_for_folds(bad_labels, train_data, "test_labels"),
    "fold ids.*subject ids"
  )
})


test_that("run_obs_crossfit_from_data slices run folds and aligns held-out data", {
  neuralign:::.register_procrustes()

  set.seed(100)
  d <- 4
  runs <- c("r1", "r1", "r2", "r2", "r3", "r3")

  Z <- matrix(rnorm(d * length(runs)), d, length(runs))
  Q <- qr.Q(qr(matrix(rnorm(d * d), d, d)))
  if (det(Q) < 0) Q[, 1] <- -Q[, 1]

  adat <- AlignmentData(
    list(s1 = Z, s2 = Q %*% Z),
    obs_labels = paste0("t", seq_along(runs))
  )

  spec <- create_obs_folds(runs, method = "run", guard_tr = 0)
  res <- run_obs_crossfit_from_data(
    data = adat,
    obs_folds = spec,
    method = "procrustes",
    reference = "s1",
    anchor_policy = "common_or_error"
  )

  expect_s3_class(res, "ObsCrossfitAlignment")
  expect_true(isTRUE(res$anchor_common))
  expect_equal(res$fold_info$obs_folds$guard_tr, 0L)
  expect_setequal(res$fold_info$obs_folds$fold_ids, spec$fold_ids)

  for (fid in names(spec$folds)) {
    aligned <- res$aligned_test_by_fold[[fid]]
    expect_true(max(abs(aligned$s1 - aligned$s2)) < 1e-8)
  }
})

test_that("run_obs_crossfit_from_data supports per-subject run ids with overlap (union policy)", {
  neuralign:::.register_procrustes_graph()

  set.seed(101)
  d <- 2

  obs_ids <- list(
    s1 = c("A", "A", "B", "B", "C", "C"),
    s2 = c("B", "B", "C", "C", "D", "D")
  )
  obs_labels <- list(
    s1 = c("A-1", "A-2", "B-1", "B-2", "C-1", "C-2"),
    s2 = c("B-1", "B-2", "C-1", "C-2", "D-1", "D-2")
  )

  all_labels <- c(obs_labels$s1, "D-1", "D-2")
  Z <- matrix(rnorm(d * length(all_labels)), d, length(all_labels))
  Q <- qr.Q(qr(matrix(rnorm(d * d), d, d)))
  if (det(Q) < 0) Q[, 1] <- -Q[, 1]

  X1 <- Z[, match(obs_labels$s1, all_labels), drop = FALSE]
  X2 <- Q %*% Z[, match(obs_labels$s2, all_labels), drop = FALSE]

  adat <- AlignmentData(
    list(s1 = X1, s2 = X2),
    obs_labels = obs_labels
  )

  spec <- create_obs_folds(obs_ids, method = "run", id_policy = "union", guard_tr = 0)
  res <- run_obs_crossfit_from_data(
    data = adat,
    obs_folds = spec,
    method = "procrustes_graph",
    reference = "s1",
    anchor_policy = "common_or_error"
  )

  expect_s3_class(res, "ObsCrossfitAlignment")
  expect_equal(res$fold_info$obs_folds$id_policy, "union")
  expect_setequal(res$fold_info$obs_folds$common_ids, c("B", "C"))

  # Folds with missing runs yield empty test sets for some subjects.
  expect_equal(ncol(res$aligned_test_by_fold$A$s2), 0L)
  expect_equal(ncol(res$aligned_test_by_fold$D$s1), 0L)

  # Overlapping folds should align correctly.
  for (fid in c("B", "C")) {
    aligned <- res$aligned_test_by_fold[[fid]]
    expect_true(max(abs(aligned$s1 - aligned$s2)) < 1e-8)
  }
})


# ---------- More obs crossfit tests ----------

test_that("run_obs_crossfit_alignment errors on extra subjects in fold", {
  train_data <- list(
    f1 = list(s1 = matrix(1, 2, 3), s2 = matrix(1, 2, 3)),
    f2 = list(s1 = matrix(1, 2, 3), s2 = matrix(1, 2, 3), s3 = matrix(1, 2, 3))
  )

  expect_error(
    run_obs_crossfit_alignment(
      train_data_by_fold = train_data,
      method = "procrustes",
      reference = "s1"
    ),
    "unexpected subjects"
  )
})

test_that(".normalize_obs_labels_for_folds per-fold per-subject named list", {
  train_data <- list(
    f1 = list(s1 = matrix(1, 2, 3), s2 = matrix(1, 2, 3)),
    f2 = list(s1 = matrix(1, 2, 4), s2 = matrix(1, 2, 4))
  )

  per_fold_per_subj <- list(
    f1 = list(s1 = c("a", "b", "c"), s2 = c("a", "b", "c")),
    f2 = list(s1 = c("d", "e", "f", "g"), s2 = c("d", "e", "f", "g"))
  )

  result <- neuralign:::.normalize_obs_labels_for_folds(
    per_fold_per_subj, train_data, "obs_labels_train"
  )
  expect_equal(length(result), 2)
  expect_equal(result$f1$s1, c("a", "b", "c"))
  expect_equal(result$f2$s2, c("d", "e", "f", "g"))
})

test_that(".normalize_obs_labels_for_folds per-fold errors on bad nested list", {
  train_data <- list(
    f1 = list(s1 = matrix(1, 2, 3), s2 = matrix(1, 2, 3))
  )

  bad_nested <- list(f1 = list(c("a", "b", "c")))  # unnamed inner list

  expect_error(
    neuralign:::.normalize_obs_labels_for_folds(bad_nested, train_data, "obs_labels_train"),
    "named per-subject list"
  )
})

test_that(".normalize_obs_labels_for_folds per-fold errors on missing subjects", {
  train_data <- list(
    f1 = list(s1 = matrix(1, 2, 3), s2 = matrix(1, 2, 3))
  )

  missing_subj <- list(f1 = list(s1 = c("a", "b", "c")))  # missing s2

  expect_error(
    neuralign:::.normalize_obs_labels_for_folds(missing_subj, train_data, "obs_labels_train"),
    "missing subjects"
  )
})

test_that(".normalize_obs_labels_for_folds per-fold errors on length mismatch", {
  train_data <- list(
    f1 = list(s1 = matrix(1, 2, 3), s2 = matrix(1, 2, 3))
  )

  bad_len <- list(f1 = list(s1 = c("a", "b"), s2 = c("a", "b", "c")))  # s1 has 2 not 3

  expect_error(
    neuralign:::.normalize_obs_labels_for_folds(bad_len, train_data, "obs_labels_train"),
    "length mismatch"
  )
})

test_that("run_obs_crossfit_alignment template mapping with obs_labels", {
  neuralign:::.register_procrustes()

  set.seed(30)
  d <- 4
  n <- 6
  Z <- matrix(rnorm(d * n), d, n)
  Q <- qr.Q(qr(matrix(rnorm(d * d), d, d)))
  if (det(Q) < 0) Q[, 1] <- -Q[, 1]

  labels <- paste0("obs", 1:3)
  train_data <- list(
    f1 = list(s1 = Z[, 1:3, drop = FALSE], s2 = (Q %*% Z)[, 1:3, drop = FALSE]),
    f2 = list(s1 = Z[, 4:6, drop = FALSE], s2 = (Q %*% Z)[, 4:6, drop = FALSE])
  )

  template <- Z
  colnames(template) <- paste0("obs", 1:6)

  res <- run_obs_crossfit_alignment(
    train_data_by_fold = train_data,
    method = "procrustes",
    reference = "consensus",
    anchor_policy = "map_to_template",
    template = template,
    template_obs_labels = paste0("obs", 1:6),
    obs_labels_train = labels,
    min_overlap = 2
  )

  expect_true(isTRUE(res$anchor_common))
  expect_equal(res$provenance$anchor_policy, "map_to_template")
})

test_that("run_obs_crossfit_alignment template mapping errors without labels", {
  neuralign:::.register_procrustes()

  set.seed(31)
  d <- 3
  Z <- matrix(rnorm(d * 4), d, 4)

  train_data <- list(
    f1 = list(s1 = Z[, 1:2, drop = FALSE], s2 = Z[, 1:2, drop = FALSE]),
    f2 = list(s1 = Z[, 3:4, drop = FALSE], s2 = Z[, 3:4, drop = FALSE])
  )

  # Template without colnames and without template_obs_labels
  template <- matrix(rnorm(d * 4), d, 4)

  expect_error(
    run_obs_crossfit_alignment(
      train_data_by_fold = train_data,
      method = "procrustes",
      reference = "consensus",
      anchor_policy = "map_to_template",
      template = template
    ),
    "template"
  )
})
