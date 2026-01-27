# Tests for AlignmentResult class and methods

test_that("AlignmentResult constructor validates inputs", {
  transforms <- list("sub-01" = diag(5), "sub-02" = diag(5))
  model <- AlignmentModel(transforms = transforms, reference = "consensus", method = "test")
  aligned <- list("sub-01" = matrix(1:10, 5, 2), "sub-02" = matrix(11:20, 5, 2))

  # Valid result
  result <- AlignmentResult(model = model, aligned = aligned)
  expect_s4_class(result, "AlignmentResult")

  # Invalid: model must be AlignmentModel
  expect_error(
    AlignmentResult(model = "not a model", aligned = aligned),
    "must be an AlignmentModel"
  )

  # Invalid: aligned must be a list
  expect_error(
    AlignmentResult(model = model, aligned = "not a list"),
    "must be a list"
  )
})

test_that("AlignmentResult accessors work correctly", {
  transforms <- list("sub-01" = diag(5), "sub-02" = diag(5))
  model <- AlignmentModel(transforms = transforms, reference = "consensus", method = "test")
  aligned <- list("sub-01" = matrix(1:10, 5, 2), "sub-02" = matrix(11:20, 5, 2))
  quality <- list(mean_pairwise_correlation = 0.85, sd_pairwise_correlation = 0.05)
  cv_info <- list(method = "loso", n_folds = 2)

  result <- AlignmentResult(
    model = model,
    aligned = aligned,
    quality = quality,
    cv_info = cv_info
  )

  # get_model
  expect_s4_class(get_model(result), "AlignmentModel")
  expect_equal(length(get_model(result)), 2)

  # get_aligned - all
  all_aligned <- get_aligned(result)
  expect_equal(names(all_aligned), c("sub-01", "sub-02"))

  # get_aligned - specific subject
  subj_aligned <- get_aligned(result, "sub-01")
  expect_equal(subj_aligned, matrix(1:10, 5, 2))

  # get_aligned - missing subject
  expect_error(get_aligned(result, "sub-99"), "not found")

  # aligned_data (synonym)
  expect_equal(aligned_data(result), aligned)

  # get_quality - all
  expect_equal(get_quality(result), quality)

  # get_quality - specific metric
  expect_equal(get_quality(result, "mean_pairwise_correlation"), 0.85)

  # get_quality - missing metric
  expect_error(get_quality(result, "nonexistent"), "not found")

  # get_cv_info
  expect_equal(get_cv_info(result), cv_info)
  expect_equal(get_cv_info(result)$method, "loso")

  # length
  expect_equal(length(result), 2)
})

test_that("AlignmentResult accessors validate input class", {
  expect_error(get_aligned(matrix(1, 1, 1)), "AlignmentResult")
  expect_error(aligned_data(matrix(1, 1, 1)), "AlignmentResult")
  expect_error(get_quality(matrix(1, 1, 1)), "AlignmentResult")
  expect_error(get_model(matrix(1, 1, 1)), "AlignmentResult")
  expect_error(get_cv_info(matrix(1, 1, 1)), "AlignmentResult")
})

test_that("as_aligned_matrix works", {
  transforms <- list("sub-01" = diag(3), "sub-02" = diag(3))
  model <- AlignmentModel(transforms = transforms, reference = "consensus", method = "test")
  aligned <- list(
    "sub-01" = matrix(1:6, 3, 2),
    "sub-02" = matrix(7:12, 3, 2)
  )

  result <- AlignmentResult(model = model, aligned = aligned)

  # By subject (default)
  by_subj <- as_aligned_matrix(result, by = "subject")
  expect_equal(by_subj, aligned)

  # By observation (cbind)
  by_obs <- as_aligned_matrix(result, by = "observation")
  expect_equal(dim(by_obs), c(3, 4))  # 3 features, 2+2 observations
  expect_equal(by_obs[, 1:2], matrix(1:6, 3, 2))
  expect_equal(by_obs[, 3:4], matrix(7:12, 3, 2))
})

test_that("AlignmentResult subsetting works", {
  transforms <- list(
    "sub-01" = diag(3),
    "sub-02" = diag(3),
    "sub-03" = diag(3)
  )
  model <- AlignmentModel(transforms = transforms, reference = "consensus", method = "test")
  aligned <- list(
    "sub-01" = matrix(1, 3, 2),
    "sub-02" = matrix(2, 3, 2),
    "sub-03" = matrix(3, 3, 2)
  )

  result <- AlignmentResult(model = model, aligned = aligned)

  # Subset by character
  sub_result <- result[c("sub-01", "sub-03")]
  expect_equal(length(sub_result), 2)
  expect_equal(names(get_aligned(sub_result)), c("sub-01", "sub-03"))

  # Subset by index
  sub_result2 <- result[1:2]
  expect_equal(length(sub_result2), 2)

  # Invalid subset
  expect_error(result["sub-99"], "Unknown subjects")
})

test_that("AlignmentResult subsetting validates index bounds and types", {
  transforms <- list(
    "sub-01" = diag(3),
    "sub-02" = diag(3),
    "sub-03" = diag(3)
  )
  model <- AlignmentModel(transforms = transforms, reference = "consensus", method = "test")
  aligned <- list(
    "sub-01" = matrix(1, 3, 2),
    "sub-02" = matrix(2, 3, 2),
    "sub-03" = matrix(3, 3, 2)
  )
  result <- AlignmentResult(model = model, aligned = aligned)

  expect_error(result[4], "out of bounds")
  expect_error(result[-4], "out of bounds")
  expect_error(result[c(-1, 2)], "cannot mix negative and positive")
  expect_error(result[1.5], "integer-valued")

  empty <- result[0]
  expect_s4_class(empty, "AlignmentResult")
  expect_equal(length(empty), 0)
})

test_that("AlignmentResult show method works", {
  transforms <- list("sub-01" = diag(10), "sub-02" = diag(10))
  model <- AlignmentModel(transforms = transforms, reference = "consensus", method = "procrustes")
  aligned <- list(
    "sub-01" = matrix(rnorm(50), 10, 5),
    "sub-02" = matrix(rnorm(50), 10, 5)
  )
  quality <- list(
    mean_pairwise_correlation = 0.85,
    pairwise_correlations = c(0.8, 0.9)
  )
  cv_info <- list(method = "loso", n_folds = 2)

  result <- AlignmentResult(
    model = model,
    aligned = aligned,
    quality = quality,
    cv_info = cv_info
  )

  output <- capture.output(show(result))

  expect_true(any(grepl("AlignmentResult", output)))
  expect_true(any(grepl("procrustes", output)))
  expect_true(any(grepl("Aligned subjects: 2", output)))
  expect_true(any(grepl("10 x 5", output)))
  expect_true(any(grepl("Quality metrics", output)))
  expect_true(any(grepl("mean_pairwise_correlation", output)))
  expect_true(any(grepl("Cross-validation", output)))
  expect_true(any(grepl("loso", output)))
})

test_that("AlignmentResult show handles vector quality metrics", {
  transforms <- list("sub-01" = diag(5))
  model <- AlignmentModel(transforms = transforms, reference = "consensus", method = "test")
  aligned <- list("sub-01" = matrix(1, 5, 3))
  quality <- list(
    per_subject_scores = c(0.7, 0.8, 0.9)  # Vector metric
  )

  result <- AlignmentResult(model = model, aligned = aligned, quality = quality)

  output <- capture.output(show(result))
  expect_true(any(grepl("per_subject_scores", output)))
  expect_true(any(grepl("mean=", output)))
})
