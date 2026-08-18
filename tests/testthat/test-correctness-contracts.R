test_that("low-rank maps report target nrow(U) and source nrow(V)", {
  U <- matrix(seq_len(14), 7, 2)
  V <- matrix(seq_len(10), 5, 2)
  lr <- neuralign:::.new_low_rank_transform(U, V)
  expect_equal(neuralign:::.transform_target_dim(lr), 7L)
  expect_equal(neuralign:::.transform_source_dim(lr), 5L)
  expect_equal(unname(neuralign:::.transform_dims(lr)[["target"]]), 7)
  expect_equal(unname(neuralign:::.transform_dims(lr)[["source"]]), 5)
})


test_that("scaled Procrustes auto inverse is scaled transpose, not raw transpose", {
  ensure_test_aligner("procrustes")
  set.seed(11)
  d <- 4
  n <- 8
  Z <- matrix(rnorm(d * n), d, n)
  Q <- qr.Q(qr(matrix(rnorm(d * d), d, d)))
  if (det(Q) < 0) Q[, 1] <- -Q[, 1]
  scales <- c(0.5, 2, 3)
  for (s in scales) {
    adat <- AlignmentData(list(
      s1 = Z,
      s2 = s * (Q %*% Z)
    ))
    res <- fit_alignment(
      adat, method = "procrustes", reference = "s1",
      compute_quality = FALSE, scale = TRUE
    )
    mdl <- get_model(res)
    A <- get_transform(mdl, "s2")
    inv_auto <- inverse_transform(mdl, "s2", method = "auto")
    expect_equal(as.matrix(A %*% inv_auto), diag(d), tolerance = 1e-6)
    expect_equal(as.matrix(inv_auto %*% A), diag(d), tolerance = 1e-6)
    raw_t <- t(as.matrix(A))
    expect_false(isTRUE(all.equal(as.matrix(inv_auto), raw_t, tolerance = 1e-6)))
  }
})


test_that("data-driven CV cannot expose or stack a global aligned matrix", {
  ensure_test_aligner("procrustes")
  set.seed(12)
  adat <- AlignmentData(list(
    s1 = matrix(rnorm(40), 8, 5),
    s2 = matrix(rnorm(40), 8, 5),
    s3 = matrix(rnorm(40), 8, 5)
  ))
  res <- fit_alignment(
    adat, method = "procrustes", reference = "medoid",
    cv = "loso", compute_quality = FALSE
  )
  expect_false(isTRUE(get_cv_info(res)$anchor_common))
  expect_error(get_aligned(res), "fold-specific|common space|anchor")
  expect_error(as_aligned_matrix(res, by = "observation"), "fold-specific|common space|anchor")
  expect_error(alignment_quality(res), "fold-specific|common space|anchor")
  expect_error(compose_alignment(get_model(res), get_model(res)), "fold-specific|common space|anchor")
})


test_that("retained fold models reproduce stored assessment outputs", {
  ensure_test_aligner("procrustes")
  set.seed(13)
  adat <- AlignmentData(list(
    s1 = matrix(rnorm(40), 8, 5),
    s2 = matrix(rnorm(40), 8, 5),
    s3 = matrix(rnorm(40), 8, 5)
  ))
  res <- fit_alignment(
    adat, method = "procrustes", reference = "medoid",
    cv = "loso", compute_quality = FALSE,
    return_resample_artifacts = TRUE
  )
  artifacts <- get_cv_info(res)$artifacts_by_fold
  expect_true(is.list(artifacts) && length(artifacts) > 0)
  for (fold in artifacts) {
    mdl <- fold$model
    for (subj in fold$assessment_subjects) {
      stored <- fold$aligned_assessment[[subj]]
      replayed <- apply_transform(get_transform(mdl, subj), get_subject_data(adat, subj))
      expect_equal(replayed, stored, tolerance = 1e-8)
    }
  }
})


test_that("AlignmentModel rejects malformed operators and mixed codomains", {
  expect_error(
    AlignmentModel(
      transforms = list(s1 = diag(2), s1 = diag(2)),
      reference = "s1",
      method = "procrustes"
    ),
    "unique|duplicate"
  )
  expect_error(
    AlignmentModel(
      transforms = list(s1 = matrix(c(1, NA, 0, 1), 2, 2)),
      reference = "s1",
      method = "procrustes"
    ),
    "finite|NA"
  )
  expect_error(
    AlignmentModel(
      transforms = list(s1 = matrix(c(1, Inf, 0, 1), 2, 2)),
      reference = "s1",
      method = "procrustes"
    ),
    "finite"
  )
  expect_error(
    AlignmentModel(
      transforms = list(s1 = "not-an-operator"),
      reference = "s1",
      method = "procrustes"
    ),
    "operator|unsupported|transform"
  )
  expect_error(
    AlignmentModel(
      transforms = list(s1 = diag(2), s2 = diag(3)),
      reference = "s1",
      method = "procrustes"
    ),
    "target|codomain|dimension"
  )
})


test_that("composition is sequential apply and fails closed on unknown or unequal spaces", {
  set.seed(14)
  Q1 <- qr.Q(qr(matrix(rnorm(9), 3, 3)))
  Q2 <- qr.Q(qr(matrix(rnorm(9), 3, 3)))
  m1 <- AlignmentModel(
    list(s1 = Q1),
    reference = "s1",
    method = "procrustes",
    space_from = "A",
    space_to = "B"
  )
  m2 <- AlignmentModel(
    list(s1 = Q2),
    reference = "s1",
    method = "procrustes",
    space_from = "B",
    space_to = "C"
  )
  X <- matrix(rnorm(12), 3, 4)
  composed <- compose_alignment(m1, m2)
  expect_equal(
    apply_transform(get_transform(composed, "s1"), X),
    apply_transform(Q2, apply_transform(Q1, X)),
    tolerance = 1e-10
  )

  unknown <- AlignmentModel(list(s1 = Q1), reference = "s1", method = "procrustes")
  expect_error(compose_alignment(unknown, m2), "unknown|unverified|space")

  m_bad <- AlignmentModel(
    list(s1 = Q2),
    reference = "s1",
    method = "procrustes",
    space_from = "Z",
    space_to = "C"
  )
  expect_error(compose_alignment(m1, m_bad), "space|mismatch|unequal")
})


test_that("observation-CV evaluation is not attached to the deployment refit", {
  ensure_test_aligner("procrustes")
  set.seed(15)
  n_feat <- 5
  n_obs <- 8
  adat <- AlignmentData(list(
    s1 = matrix(rnorm(n_feat * n_obs), n_feat, n_obs),
    s2 = matrix(rnorm(n_feat * n_obs), n_feat, n_obs),
    s3 = matrix(rnorm(n_feat * n_obs), n_feat, n_obs)
  ))
  spec <- list(
    axis = "observation",
    method = "custom",
    folds = list(
      fold1 = list(train_idx = 1:4, test_idx = 5:8),
      fold2 = list(train_idx = 5:8, test_idx = 1:4)
    ),
    n_folds = 2L,
    fold_ids = c("fold1", "fold2")
  )

  res <- fit_alignment(
    adat, method = "procrustes", reference = "medoid",
    cv_folds = spec, compute_quality = FALSE
  )
  info <- get_cv_info(res)
  expect_false(isTRUE(info$anchor_common))
  expect_error(get_aligned(res), "fold-specific|common space|anchor")
  expect_true(isTRUE(info$deployment_refit))
  expect_length(res@aligned, 0L)
  artifacts <- info$artifacts_by_fold
  expect_true(is.list(artifacts) && length(artifacts) == 2L)

  for (fid in names(artifacts)) {
    fold <- artifacts[[fid]]
    mdl <- fold$model
    test_idx <- fold$test_idx
    for (subj in fold$assessment_subjects) {
      stored <- fold$aligned_assessment[[subj]]
      replayed <- apply_transform(
        get_transform(mdl, subj),
        get_subject_data(adat, subj)[, test_idx, drop = FALSE]
      )
      expect_equal(replayed, stored, tolerance = 1e-8)
    }
  }

  deploy <- info$deployment_model
  expect_s4_class(deploy, "AlignmentModel")
  expect_false(identical(deploy, artifacts$fold1$model))
})


test_that("partial compose requires an explicit subject policy", {
  Q <- diag(2)
  m1 <- AlignmentModel(
    list(s1 = Q, s2 = Q),
    reference = "s1",
    method = "procrustes",
    space_from = "A",
    space_to = "B"
  )
  m2 <- AlignmentModel(
    list(s1 = Q, s3 = Q),
    reference = "s1",
    method = "procrustes",
    space_from = "B",
    space_to = "C"
  )
  expect_error(compose_alignment(m1, m2), "allow_partial|subjects")
  partial <- compose_alignment(m1, m2, allow_partial = TRUE)
  expect_equal(model_subjects(partial), "s1")
  explicit <- compose_alignment(m1, m2, subjects = "s1")
  expect_equal(model_subjects(explicit), "s1")
})
