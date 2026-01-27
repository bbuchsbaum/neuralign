test_that("fit_alignment resolves reference per observation fold (medoid)", {
  neuralign:::.register_procrustes()

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
  neuralign:::.register_procrustes()

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
