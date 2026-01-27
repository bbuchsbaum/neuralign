test_that("preprocess_matrix preserves dimnames and centers/scales rows", {
  x <- matrix(c(1, 2, 3, 4,
                10, 20, 30, 40),
              nrow = 2,
              byrow = TRUE)
  rownames(x) <- c("f1", "f2")
  colnames(x) <- paste0("o", 1:4)

  centered <- preprocess_matrix(x, center = "rows")
  expect_identical(dimnames(centered), dimnames(x))
  expect_equal(unname(rowMeans(centered)), c(0, 0), tolerance = 1e-12)

  scaled <- preprocess_matrix(x, center = "rows", scale = "rows")
  expect_equal(unname(rowMeans(scaled)), c(0, 0), tolerance = 1e-12)
  expect_equal(unname(apply(scaled, 1, stats::sd)), c(1, 1), tolerance = 1e-12)
})

test_that("preprocess_matrix centers columns and errors on non-finite values", {
  x <- matrix(c(1, 2, 3, 4,
                10, 20, 30, 40),
              nrow = 2,
              byrow = TRUE)
  centered <- preprocess_matrix(x, center = "cols")
  expect_equal(colMeans(centered), rep(0, ncol(x)), tolerance = 1e-12)

  x_bad <- x
  x_bad[1, 1] <- NA_real_
  expect_error(preprocess_matrix(x_bad), "non-finite")
})

test_that("preprocess_matrix robust centering uses medians", {
  x <- matrix(c(0, 0, 0, 100,
                1, 1, 1, 1),
              nrow = 2,
              byrow = TRUE)

  centered <- preprocess_matrix(x, center = "rows", robust = TRUE)
  expect_equal(apply(centered, 1, stats::median), c(0, 0), tolerance = 1e-12)
})

test_that("preprocess_alignment_data applies preprocessing per subject", {
  subject_ids <- make_test_subject_ids(2)
  obs_labels <- make_test_obs_labels(4)
  data_list <- make_test_data_list(
    n_subjects = 2,
    n_features = 3,
    n_obs = 4,
    subject_ids = subject_ids
  )
  adat <- AlignmentData(data_list, obs_labels = obs_labels)

  out <- preprocess_alignment_data(adat, center = "rows")
  expect_s4_class(out, "AlignmentData")
  expect_identical(out@subjects, adat@subjects)
  expect_identical(out@obs_labels, adat@obs_labels)
  means <- lapply(get_data_list(out), rowMeans)
  expect_true(all(vapply(means, function(m) all(abs(m) < 1e-12), logical(1))))
})

test_that("preprocess_alignment_data no-op does not change fit_alignment() results", {
  subject_ids <- make_test_subject_ids(3)
  obs_labels <- make_test_obs_labels(6)
  data_list <- make_test_data_list(
    n_subjects = 3,
    n_features = 5,
    n_obs = 6,
    subject_ids = subject_ids
  )
  adat <- AlignmentData(data_list, obs_labels = obs_labels)

  res1 <- fit_alignment(
    adat,
    method = "procrustes",
    reference = subject_ids[[1]],
    compute_quality = FALSE,
    return_aligned = FALSE
  )
  res2 <- fit_alignment(
    preprocess_alignment_data(adat),
    method = "procrustes",
    reference = subject_ids[[1]],
    compute_quality = FALSE,
    return_aligned = FALSE
  )

  expect_equal(get_transforms(res1@model), get_transforms(res2@model), tolerance = 1e-12)
})

test_that("preprocess_feature_blocks preprocesses each block matrix", {
  x1 <- matrix(c(1, 2, 3, 4), nrow = 2, byrow = TRUE)
  x2 <- matrix(c(2, 4, 6, 8), nrow = 2, byrow = TRUE)

  blocks_by_subject <- list(
    "sub-01" = list(
      a = alignment_feature_block(x1, name = "a", feature_names = c("f1", "f2")),
      b = alignment_feature_block(x2, name = "b", feature_names = c("g1", "g2"))
    ),
    "sub-02" = list(
      a = alignment_feature_block(x2, name = "a", feature_names = c("f1", "f2")),
      b = alignment_feature_block(x1, name = "b", feature_names = c("g1", "g2"))
    )
  )

  out <- preprocess_feature_blocks(blocks_by_subject, center = "rows")
  expect_true(all(vapply(out, function(bl) all(vapply(bl, inherits, logical(1), "alignment_feature_block")), logical(1))))

  means <- lapply(out, function(bl) lapply(bl, function(b) rowMeans(b$x)))
  expect_true(all(vapply(unlist(means, recursive = FALSE), function(m) all(abs(m) < 1e-12), logical(1))))
})
