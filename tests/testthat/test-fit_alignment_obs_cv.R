test_that("fit_alignment resolves reference per observation fold (medoid)", {
  ensure_test_aligner("procrustes")

  set.seed(1)
  n_feat <- 6
  n_obs <- 6

  base <- matrix(rnorm(n_feat * n_obs), n_feat, n_obs)
  data_list <- list(
    s1 = base,
    s2 = base + 0.1 * matrix(rnorm(n_feat * n_obs), n_feat, n_obs),
    s3 = base + 0.2 * matrix(rnorm(n_feat * n_obs), n_feat, n_obs)
  )
  adat <- AlignmentData(data_list)

  runs <- rep(c("r1", "r2", "r3"), each = 2)
  spec <- create_obs_folds(runs, method = "run")

  res <- NULL
  expect_warning(
    res <- fit_alignment(
      adat,
      method = "procrustes",
      reference = "medoid",
      cv_folds = spec,
      compute_quality = FALSE,
      return_aligned = TRUE
    ),
    "fold-specific anchors"
  )

  expect_s4_class(res, "AlignmentResult")
  info <- get_cv_info(res)
  expect_equal(info$axis, "observation")
  expect_false(isTRUE(info$anchor_common))
  expect_true(all(names(info$reference_by_fold) %in% spec$fold_ids))
  expect_true(all(info$reference_by_fold %in% names(data_list)))
})

test_that("fit_alignment obs-CV can retain fold transforms", {
  ensure_test_aligner("procrustes")

  set.seed(2)
  n_feat <- 5
  n_obs <- 6
  data_list <- list(
    s1 = matrix(rnorm(n_feat * n_obs), n_feat, n_obs),
    s2 = matrix(rnorm(n_feat * n_obs), n_feat, n_obs),
    s3 = matrix(rnorm(n_feat * n_obs), n_feat, n_obs)
  )
  adat <- AlignmentData(data_list)

  runs <- rep(c("r1", "r2", "r3"), each = 2)
  spec <- create_obs_folds(runs, method = "run")

  res <- fit_alignment(
    adat,
    method = "procrustes",
    reference = "s1",
    cv_folds = spec,
    compute_quality = FALSE,
    return_aligned = FALSE,
    return_fold_transforms = TRUE
  )

  info <- get_cv_info(res)
  expect_true(is.list(info$transforms_by_fold))
  expect_equal(length(info$transforms_by_fold), length(spec$fold_ids))

  for (fid in spec$fold_ids) {
    tf <- info$transforms_by_fold[[fid]]
    expect_true(is.list(tf))
    expect_true(all(names(data_list) %in% names(tf)))
    expect_equal(dim(tf$s1), c(n_feat, n_feat))
  }
})


# ---------- Obs-axis CV: fixed-subject reference (anchor_common = TRUE) ----------

test_that("obs-CV with fixed-subject reference does NOT warn about fold-specific anchors", {

  ensure_test_aligner("procrustes")

  set.seed(10)
  n_feat <- 5
  n_obs <- 6
  data_list <- list(
    s1 = matrix(rnorm(n_feat * n_obs), n_feat, n_obs),
    s2 = matrix(rnorm(n_feat * n_obs), n_feat, n_obs)
  )
  adat <- AlignmentData(data_list)

  runs <- rep(c("r1", "r2", "r3"), each = 2)
  spec <- create_obs_folds(runs, method = "run")

  # Fixed-subject reference => anchor_common = TRUE => no warning
  expect_no_warning(
    res <- fit_alignment(
      adat,
      method = "procrustes",
      reference = "s1",
      cv_folds = spec,
      compute_quality = FALSE
    )
  )

  info <- get_cv_info(res)
  expect_true(isTRUE(info$anchor_common))
  expect_equal(info$axis, "observation")
  # reference_by_fold should all be "s1"
  expect_true(all(info$reference_by_fold == "s1"))
})


# ---------- Obs-axis CV: full-coverage reassembly ----------

test_that("obs-CV with complete folds reassembles aligned data in original observation order", {
  ensure_test_aligner("procrustes")

  set.seed(20)
  n_feat <- 4
  n_obs <- 10
  data_list <- list(
    s1 = matrix(rnorm(n_feat * n_obs), n_feat, n_obs),
    s2 = matrix(rnorm(n_feat * n_obs), n_feat, n_obs)
  )
  adat <- AlignmentData(data_list)

  # Two folds covering all 10 observations exactly once
  spec <- list(
    axis = "observation",
    method = "custom",
    folds = list(
      fold1 = list(train_idx = 1:5, test_idx = 6:10),
      fold2 = list(train_idx = 6:10, test_idx = 1:5)
    ),
    n_folds = 2L,
    fold_ids = c("fold1", "fold2")
  )

  res <- fit_alignment(
    adat,
    method = "procrustes",
    reference = "s1",
    cv_folds = spec,
    compute_quality = FALSE,
    return_aligned = TRUE
  )

  # Aligned data should have same ncol as original
  expect_equal(ncol(res@aligned$s1), n_obs)
  expect_equal(ncol(res@aligned$s2), n_obs)
  expect_equal(nrow(res@aligned$s1), n_feat)
})


# ---------- Obs-axis CV: return_fold_transforms = FALSE (default) ----------

test_that("obs-CV with return_fold_transforms=FALSE has NULL transforms_by_fold", {
  ensure_test_aligner("procrustes")

  set.seed(30)
  n_feat <- 4
  n_obs <- 6
  data_list <- list(
    s1 = matrix(rnorm(n_feat * n_obs), n_feat, n_obs),
    s2 = matrix(rnorm(n_feat * n_obs), n_feat, n_obs)
  )
  adat <- AlignmentData(data_list)

  runs <- rep(c("r1", "r2", "r3"), each = 2)
  spec <- create_obs_folds(runs, method = "run")

  res <- fit_alignment(
    adat,
    method = "procrustes",
    reference = "s1",
    cv_folds = spec,
    compute_quality = FALSE,
    return_fold_transforms = FALSE
  )

  info <- get_cv_info(res)
  expect_null(info$transforms_by_fold)
})


# ---------- Obs-axis CV: consensus reference triggers fold-specific warning ----------

test_that("obs-CV with consensus reference warns about fold-specific anchors", {
  ensure_test_aligner("procrustes")

  set.seed(40)
  n_feat <- 5
  n_obs <- 6
  data_list <- list(
    s1 = matrix(rnorm(n_feat * n_obs), n_feat, n_obs),
    s2 = matrix(rnorm(n_feat * n_obs), n_feat, n_obs),
    s3 = matrix(rnorm(n_feat * n_obs), n_feat, n_obs)
  )
  adat <- AlignmentData(data_list)

  runs <- rep(c("r1", "r2", "r3"), each = 2)
  spec <- create_obs_folds(runs, method = "run")

  expect_warning(
    res <- fit_alignment(
      adat,
      method = "procrustes",
      reference = "consensus",
      cv_folds = spec,
      compute_quality = FALSE
    ),
    "fold-specific anchors"
  )

  info <- get_cv_info(res)
  expect_false(isTRUE(info$anchor_common))
  expect_equal(info$reference_kind, "data_driven")
})


# ---------- Obs-axis CV: return_aligned = FALSE ----------

test_that("obs-CV with return_aligned=FALSE produces empty aligned list", {
  ensure_test_aligner("procrustes")

  set.seed(50)
  n_feat <- 4
  n_obs <- 6
  data_list <- list(
    s1 = matrix(rnorm(n_feat * n_obs), n_feat, n_obs),
    s2 = matrix(rnorm(n_feat * n_obs), n_feat, n_obs)
  )
  adat <- AlignmentData(data_list)

  spec <- list(
    axis = "observation",
    method = "custom",
    folds = list(
      fold1 = list(train_idx = 1:3, test_idx = 4:6),
      fold2 = list(train_idx = 4:6, test_idx = 1:3)
    ),
    n_folds = 2L,
    fold_ids = c("fold1", "fold2")
  )

  res <- fit_alignment(
    adat,
    method = "procrustes",
    reference = "s1",
    cv_folds = spec,
    return_aligned = FALSE
  )

  expect_equal(length(res@aligned), 0)
  # Model should still be present
  expect_equal(length(get_model(res)@transforms), 2)
})


# ---------- Obs-axis CV: template reference is anchor_common ----------

test_that("obs-CV with template reference sets anchor_common=TRUE", {
  ensure_test_aligner("procrustes")

  set.seed(60)
  n_feat <- 5
  n_obs <- 6
  template <- matrix(rnorm(n_feat * n_obs), n_feat, n_obs)
  data_list <- list(
    s1 = template + 0.1 * matrix(rnorm(n_feat * n_obs), n_feat, n_obs),
    s2 = template + 0.2 * matrix(rnorm(n_feat * n_obs), n_feat, n_obs)
  )
  adat <- AlignmentData(data_list)

  # Template reference: provide a matrix
  # Build a fold spec for just 2 folds
  spec <- list(
    axis = "observation",
    method = "custom",
    folds = list(
      fold1 = list(train_idx = 1:3, test_idx = 4:6),
      fold2 = list(train_idx = 4:6, test_idx = 1:3)
    ),
    n_folds = 2L,
    fold_ids = c("fold1", "fold2")
  )

  # Template reference => anchor_common = TRUE => no fold-specific warning
  expect_no_warning(
    res <- fit_alignment(
      adat,
      method = "procrustes",
      reference = template,
      cv_folds = spec,
      compute_quality = FALSE
    )
  )

  info <- get_cv_info(res)
  expect_true(isTRUE(info$anchor_common))
  expect_equal(info$reference_kind, "template")
})
