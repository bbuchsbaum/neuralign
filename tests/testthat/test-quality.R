test_that("alignment_quality computes correlation metrics", {
  neuralign:::.register_procrustes()

  set.seed(500)
  data_list <- list(
    "sub-01" = matrix(rnorm(100), 10, 10),
    "sub-02" = matrix(rnorm(100), 10, 10),
    "sub-03" = matrix(rnorm(100), 10, 10)
  )
  adat <- AlignmentData(data_list)

  result <- fit_alignment(adat, method = "procrustes")
  quality <- alignment_quality(result, metrics = "correlation")

  expect_true("mean_pairwise_correlation" %in% names(quality))
  expect_true(is.numeric(quality$mean_pairwise_correlation))
  expect_true(quality$mean_pairwise_correlation >= -1 &&
    quality$mean_pairwise_correlation <= 1)

  expect_true("pairwise_correlations" %in% names(quality))
  expect_equal(length(quality$pairwise_correlations), 3)  # C(3,2) = 3 pairs
})

test_that("alignment_quality handles single subject", {
  aligned <- list(
    "sub-01" = matrix(rnorm(100), 10, 10)
  )

  quality <- alignment_quality(aligned, metrics = "correlation")

  expect_true(is.na(quality$mean_pairwise_correlation))
})

test_that("alignment_quality with isc metrics works", {
  set.seed(501)

  # Create data with shared time-varying pattern
  shared_pattern <- matrix(rnorm(10), 10, 1) %*% matrix(rnorm(10), 1, 10)
  aligned <- list(
    "sub-01" = shared_pattern + matrix(rnorm(100, sd = 0.5), 10, 10),
    "sub-02" = shared_pattern + matrix(rnorm(100, sd = 0.5), 10, 10),
    "sub-03" = shared_pattern + matrix(rnorm(100, sd = 0.5), 10, 10)
  )

  quality <- alignment_quality(aligned, metrics = "isc")

  expect_true("mean_isc" %in% names(quality))
  expect_true(is.numeric(quality$mean_isc))
  # With shared pattern, ISC should be positive
  expect_true(quality$mean_isc > 0)
})

test_that("alignment_quality with reference computes reconstruction", {
  set.seed(502)

  reference <- matrix(rnorm(100), 10, 10)
  aligned <- list(
    "sub-01" = reference + matrix(rnorm(100, sd = 0.1), 10, 10),
    "sub-02" = reference + matrix(rnorm(100, sd = 0.1), 10, 10)
  )

  quality <- alignment_quality(
    aligned,
    metrics = "reconstruction",
    reference = reference
  )

  expect_true("mean_reference_correlation" %in% names(quality))
  # Should have high correlation with reference since noise is small
  expect_true(quality$mean_reference_correlation > 0.9)
})

test_that("alignment_quality computes improvement over original", {
  set.seed(503)

  # Create original data with different patterns
  original_list <- list(
    "sub-01" = matrix(rnorm(100), 10, 10),
    "sub-02" = matrix(rnorm(100), 10, 10),
    "sub-03" = matrix(rnorm(100), 10, 10)
  )
  original <- AlignmentData(original_list)

  # Create aligned data with shared pattern (should have higher correlation)
  shared <- matrix(rnorm(100), 10, 10)
  aligned_list <- list(
    "sub-01" = shared + matrix(rnorm(100, sd = 0.1), 10, 10),
    "sub-02" = shared + matrix(rnorm(100, sd = 0.1), 10, 10),
    "sub-03" = shared + matrix(rnorm(100, sd = 0.1), 10, 10)
  )

  quality <- alignment_quality(
    aligned_list,
    original = original,
    metrics = "correlation"
  )

  expect_true("correlation_improvement" %in% names(quality))
  expect_true("original_mean_correlation" %in% names(quality))
  expect_true("aligned_mean_correlation" %in% names(quality))

  # Aligned should have higher correlation
  expect_true(quality$aligned_mean_correlation > quality$original_mean_correlation)
})

test_that("print_quality_summary works", {
  quality <- list(
    mean_pairwise_correlation = 0.75,
    sd_pairwise_correlation = 0.1,
    min_pairwise_correlation = 0.6,
    max_pairwise_correlation = 0.9
  )

  expect_output(print_quality_summary(quality), "Pairwise Correlation")
  expect_output(print_quality_summary(quality), "0.75")
})

test_that("alignment_quality works with AlignmentResult", {
  neuralign:::.register_procrustes()

  set.seed(504)
  data_list <- list(
    "sub-01" = matrix(rnorm(100), 10, 10),
    "sub-02" = matrix(rnorm(100), 10, 10)
  )
  adat <- AlignmentData(data_list)

  result <- fit_alignment(adat, method = "procrustes")

  quality <- alignment_quality(result)

  expect_type(quality, "list")
  expect_true("mean_pairwise_correlation" %in% names(quality))
})

test_that("alignment_quality variance metrics require abind", {
  skip_if_not_installed("abind")

  set.seed(505)
  aligned <- list(
    "sub-01" = matrix(rnorm(100), 10, 10),
    "sub-02" = matrix(rnorm(100), 10, 10)
  )

  quality <- alignment_quality(aligned, metrics = "variance")

  expect_true("total_variance" %in% names(quality))
  expect_true("between_subject_variance" %in% names(quality))
  expect_true("within_subject_variance" %in% names(quality))
})
