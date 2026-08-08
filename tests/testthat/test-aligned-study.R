test_that("SharedFeatureSpace digest is stable and content-addressed", {
  a <- SharedFeatureSpace(
    dimension = 3,
    coordinate_id = "fixture-space",
    method = "procrustes",
    train_subjects = c("s2", "s1"),
    anchor = "s1"
  )
  b <- SharedFeatureSpace(
    dimension = 3,
    coordinate_id = "fixture-space",
    method = "procrustes",
    train_subjects = c("s1", "s2"),
    anchor = "s1"
  )
  expect_identical(a$id, b$id)
  expect_true(startsWith(a$id, "shared-"))

  c <- SharedFeatureSpace(
    dimension = 4,
    coordinate_id = "fixture-space",
    method = "procrustes",
    train_subjects = c("s1", "s2"),
    anchor = "s1"
  )
  expect_false(identical(a$id, c$id))
})


test_that("SharedFeatureSpace identity distinguishes fitted coordinates", {
  m1 <- AlignmentModel(
    transforms = list(s1 = diag(2)),
    reference = "s1",
    method = "procrustes",
    train_subjects = "s1"
  )
  m2 <- AlignmentModel(
    transforms = list(s1 = matrix(c(0, 1, 1, 0), 2)),
    reference = "s1",
    method = "procrustes",
    train_subjects = "s1"
  )

  s1 <- shared_feature_space_from_model(m1)
  s1_again <- shared_feature_space_from_model(m1)
  s2 <- shared_feature_space_from_model(m2)

  expect_identical(s1$id, s1_again$id)
  expect_false(identical(s1$id, s2$id))
  expect_false(shared_spaces_compatible(s1, s2))
  expect_error(shared_spaces_compatible(NULL, s1), "unknown|identity")
})


test_that("AlignedBlock rejects values that cannot support numeric analysis", {
  ss <- SharedFeatureSpace(2, coordinate_id = "numeric-block-fixture")
  expect_error(
    AlignedBlock(matrix(letters[1:4], 2), subject_id = "s1",
                 shared_space_id = ss$id),
    "numeric"
  )
  expect_error(
    AlignedBlock(matrix(c(1, 2, 3, Inf), 2), subject_id = "s1",
                 shared_space_id = ss$id),
    "finite"
  )

  sparse <- Matrix::Matrix(diag(2), sparse = TRUE)
  block <- AlignedBlock(sparse, subject_id = "s1", shared_space_id = ss$id)
  expect_s4_class(block$values, "Matrix")
})


test_that("orientation boundary is one explicit transpose", {
  set.seed(42)
  V <- 6
  N <- 5
  K <- 6
  X_alg <- matrix(rnorm(V * N), V, N)
  Q <- qr.Q(qr(matrix(rnorm(K * V), K, V))) # orthogonal KxV (=VxV here)

  Z_alg <- apply_transform(Q, X_alg)
  Z_an <- to_analysis_matrix(Z_alg)

  expect_equal(Z_an, to_analysis_matrix(Z_alg), tolerance = 1e-10)
  expect_equal(Z_an, to_analysis_matrix(X_alg) %*% t(Q), tolerance = 1e-10)
  expect_equal(to_algorithm_matrix(Z_an), Z_alg, tolerance = 1e-10)
})


test_that("as_aligned_study builds analysis-facing AlignedStudy", {
  set.seed(1)
  adat <- AlignmentData(list(
    s1 = matrix(rnorm(40), 8, 5),
    s2 = matrix(rnorm(40), 8, 5)
  ), obs_labels = paste0("t", 1:5))
  res <- fit_alignment(adat, method = "procrustes", reference = "s1",
                       compute_quality = FALSE)

  study <- as_aligned_study(res, source_data = adat, representation = "raw_tr")
  expect_s4_class(study, "AlignedStudy")
  expect_equal(subject_ids(study), c("s1", "s2"))
  expect_equal(dim(aligned_matrix(study, "s1")), c(5L, 8L))
  expect_equal(nrow(observation_data(study, "s1")), 5L)
  expect_equal(nrow(feature_data(study)), 8L)
  expect_true(inherits(shared_space(study), "SharedFeatureSpace"))
  assert_common_shared_space(study)

  stacked <- stack_subjects(study)
  expect_equal(nrow(stacked$matrix), 10L)
  expect_equal(nrow(stacked$index), 10L)
  expect_true(all(c("row_id", "subject_id", "observation_id") %in% names(stacked$index)))
})


test_that("align_study frozen application path works", {
  set.seed(2)
  train <- AlignmentData(list(
    s1 = matrix(rnorm(40), 8, 5),
    s2 = matrix(rnorm(40), 8, 5)
  ))
  fit <- fit_alignment(train, method = "procrustes", reference = "s1",
                       compute_quality = FALSE)

  betas <- list(
    s1 = matrix(rnorm(24), 8, 3),
    s2 = matrix(rnorm(24), 8, 3)
  )
  study <- align_study(
    fit,
    betas,
    mode = "frozen_application",
    representation = "trial_beta",
    warn_leakage = FALSE
  )
  expect_s4_class(study, "AlignedStudy")
  expect_equal(study@safety$mode, "frozen_application")
  expect_equal(dim(aligned_matrix(study, "s1")), c(3L, 8L))
  expect_equal(study@blocks$s1$representation, "trial_beta")

  expect_error(
    assert_analysis_safe(study, purpose = "confirmatory_cross_subject_prediction"),
    "declared|verified"
  )
  expect_warning(
    assert_analysis_safe(
      study,
      purpose = "confirmatory_cross_subject_prediction",
      allow_declared = TRUE
    ),
    "declared"
  )
  expect_error(
    assert_analysis_safe(
      as_aligned_study(fit, mode = "in_sample"),
      purpose = "confirmatory_cross_subject_prediction"
    ),
    "not valid|in_sample"
  )
})


test_that("align_study analysis orientation path", {
  set.seed(3)
  train <- AlignmentData(list(
    s1 = matrix(rnorm(30), 6, 5),
    s2 = matrix(rnorm(30), 6, 5)
  ))
  fit <- fit_alignment(train, method = "procrustes", reference = "s1",
                       compute_quality = FALSE)
  # analysis-facing: obs x features
  X_an <- list(
    s1 = matrix(rnorm(18), 3, 6),
    s2 = matrix(rnorm(18), 3, 6)
  )
  study <- align_study(
    fit, X_an,
    orientation = "analysis",
    mode = "frozen_application",
    warn_leakage = FALSE
  )
  expect_equal(dim(aligned_matrix(study, "s1")), c(3L, 6L))
})


test_that("model application preserves fitted coordinate identity", {
  set.seed(31)
  train <- AlignmentData(list(
    s1 = matrix(rnorm(30), 6, 5),
    s2 = matrix(rnorm(30), 6, 5)
  ))
  fit <- fit_alignment(train, method = "procrustes", reference = "s1",
                       compute_quality = FALSE)
  fitted_space <- shared_feature_space_from_model(get_model(fit))

  applied <- apply_alignment(
    fit,
    AlignmentData(list(s3 = matrix(rnorm(30), 6, 5))),
    warn_leakage = FALSE
  )
  applied_space <- shared_feature_space_from_model(get_model(applied))
  expect_identical(applied_space$id, fitted_space$id)
  expect_true(shared_spaces_compatible(applied_space, fitted_space))
})


test_that("fold-specific AlignmentResult requires retained resample artifacts", {
  set.seed(4)
  adat <- AlignmentData(list(
    s1 = matrix(rnorm(40), 8, 5),
    s2 = matrix(rnorm(40), 8, 5),
    s3 = matrix(rnorm(40), 8, 5)
  ))
  res <- fit_alignment(
    adat,
    method = "procrustes",
    reference = "medoid",
    cv = "loso",
    compute_quality = FALSE
  )
  expect_false(isTRUE(get_cv_info(res)$anchor_common))
  expect_error(as_aligned_study(res), "as_aligned_resample_set")

  expect_error(
    as_aligned_resample_set(res, source_data = adat),
    "return_resample_artifacts"
  )

  retained <- fit_alignment(
    adat,
    method = "procrustes",
    reference = "medoid",
    cv = "loso",
    compute_quality = FALSE,
    return_resample_artifacts = TRUE
  )
  ars <- as_aligned_resample_set(retained, source_data = adat)
  expect_s4_class(ars, "AlignedResampleSet")
  expect_gt(length(ars), 0)
  expect_error(stack_aligned_resamples(ars), "Cannot stack")
  for (split in ars@splits) {
    expect_s4_class(split$model, "AlignmentModel")
    expect_s4_class(split$assessment, "AlignedStudy")
    expect_identical(
      shared_space(split$assessment)$id,
      split$shared_space$id
    )
    expect_setequal(
      subject_ids(split$assessment),
      split$metadata$assessment_subjects
    )
    for (subject in split$metadata$assessment_subjects) {
      expect_equal(
        to_algorithm_matrix(aligned_matrix(split$assessment, subject)),
        apply_transform(
          get_transform(split$model, subject),
          get_subject_data(adat, subject)
        ),
        tolerance = 1e-10
      )
    }
  }
})


test_that("as_hyperdesign constructs a real multidesign hyperdesign", {
  skip_if_not_installed("multidesign")
  set.seed(5)
  adat <- AlignmentData(list(
    s1 = matrix(rnorm(40), 8, 5),
    s2 = matrix(rnorm(40), 8, 5)
  ), obs_labels = paste0("t", 1:5))
  res <- fit_alignment(adat, method = "procrustes", reference = "s1",
                       compute_quality = FALSE)
  study <- as_aligned_study(res, source_data = adat)

  hd <- as_hyperdesign(study)
  expect_true(inherits(hd, "hyperdesign"))
  expect_true(all(vapply(hd, inherits, logical(1), "multidesign")))
  expect_equal(nrow(hd$s1$x), 5L)
  expect_equal(ncol(hd$s1$x), 8L)
  expect_equal(nrow(hd$s1$design), 5L)
  expect_equal(nrow(multidesign::design(hd, block = 1)), 5L)
  expect_equal(nrow(multidesign::xdata(hd, block = 1)), 5L)

  mats_an <- as_subject_matrices(study, orientation = "analysis")
  mats_alg <- as_subject_matrices(study, orientation = "algorithm")
  expect_equal(dim(mats_an$s1), c(5L, 8L))
  expect_equal(dim(mats_alg$s1), c(8L, 5L))
})


test_that("per-aligner orientation conversion matches apply_transform", {
  set.seed(7)
  adat <- AlignmentData(list(
    s1 = matrix(rnorm(50), 10, 5),
    s2 = matrix(rnorm(50), 10, 5)
  ))
  res <- fit_alignment(adat, method = "procrustes", reference = "s1",
                       compute_quality = FALSE)
  mdl <- get_model(res)
  for (s in c("s1", "s2")) {
    X <- get_subject_data(adat, s)
    Q <- get_transform(mdl, s)
    Z_alg <- apply_transform(Q, X)
    Z_an <- to_analysis_matrix(apply_transform(Q, X))
    expect_equal(Z_an, to_analysis_matrix(Z_alg), tolerance = 1e-8)
  }
})


test_that("declared safety cannot pass a strict confirmatory assertion", {
  ss <- SharedFeatureSpace(2, coordinate_id = "declared-safety-fixture")
  block <- AlignedBlock(
    matrix(1:4, 2),
    subject_id = "s1",
    shared_space_id = ss$id
  )
  study <- AlignedStudy(
    list(s1 = block),
    ss,
    safety = analysis_safety_record(
      mode = "frozen_application",
      leakage_status = "declared_frozen_application"
    )
  )

  expect_identical(analysis_safety(study)$status, "declared")
  expect_error(
    assert_analysis_safe(
      study,
      purpose = "confirmatory_cross_subject_prediction",
      allow_declared = FALSE
    ),
    "declared|verified"
  )
})


test_that("observation crossfit converts from retained fold artifacts", {
  set.seed(8)
  adat <- AlignmentData(list(
    s1 = matrix(rnorm(48), 8, 6),
    s2 = matrix(rnorm(48), 8, 6)
  ), obs_labels = paste0("o", 1:6))
  folds <- create_obs_folds(rep(c("r1", "r2"), each = 3), method = "run")
  crossfit <- NULL
  expect_warning(
    crossfit <- run_obs_crossfit_from_data(
      adat,
      folds,
      method = "procrustes",
      reference = "medoid",
      anchor_policy = "fold_specific_ok"
    ),
    "Fold-specific anchors"
  )

  ars <- as_aligned_resample_set(crossfit)
  expect_s4_class(ars, "AlignedResampleSet")
  expect_identical(length(ars), length(folds$folds))
  expect_false("oof_view" %in% names(ars@splits))
  expect_true(all(vapply(
    ars@splits,
    function(split) inherits(split$assessment, "AlignedStudy"),
    logical(1)
  )))
  for (fold_id in names(ars@splits)) {
    split <- ars@splits[[fold_id]]
    fold <- folds$folds[[fold_id]]
    for (subject in subject_ids(split$assessment)) {
      expected <- apply_transform(
        get_transform(split$model, subject),
        get_subject_data(adat, subject)[, fold$test_idx, drop = FALSE]
      )
      expect_equal(
        to_algorithm_matrix(aligned_matrix(split$assessment, subject)),
        expected,
        tolerance = 1e-10
      )
      expect_identical(
        observation_data(split$assessment, subject)$source_index,
        as.integer(fold$test_idx)
      )
    }
  }
})
