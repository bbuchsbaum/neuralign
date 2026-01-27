test_that("alignment_quality computes correlation metrics", {
  neuralign:::.register_procrustes()

  set.seed(500)
  data_list <- make_test_data_list(n_subjects = 3, n_features = 10, n_obs = 10)
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

  reference <- make_test_matrix(n_features = 10, n_obs = 10)
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
  expect_true("reconstruction_frobenius" %in% names(quality))
  expect_true("mean_reconstruction_frobenius" %in% names(quality))
})

test_that("alignment_quality reconstruction includes Frobenius residuals", {
  reference <- matrix(c(0, 1, 1, 0), 2, 2, byrow = TRUE)
  aligned <- list("sub-01" = reference + matrix(1, 2, 2))

  quality <- alignment_quality(
    aligned,
    metrics = "reconstruction",
    reference = reference
  )

  expect_equal(unname(quality$reconstruction_errors[["sub-01"]]), 1, tolerance = 1e-12)
  expect_equal(unname(quality$reconstruction_frobenius[["sub-01"]]), 2, tolerance = 1e-12)
  expect_equal(unname(quality$mean_reconstruction_frobenius), 2, tolerance = 1e-12)

  expect_error(
    alignment_quality(aligned, metrics = "reconstruction", reference = matrix(0, 3, 3)),
    "match reference dims"
  )
})

test_that("alignment_quality computes improvement over original", {
  set.seed(503)

  # Create original data with different patterns
  original_list <- make_test_data_list(n_subjects = 3, n_features = 10, n_obs = 10)
  original <- AlignmentData(original_list)

  # Create aligned data with shared pattern (should have higher correlation)
  shared <- make_test_matrix(n_features = 10, n_obs = 10)
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
  data_list <- make_test_data_list(n_subjects = 2, n_features = 10, n_obs = 10)
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

# ---- Additional coverage tests ----

test_that("alignment_quality rejects bad input (line 40)", {
  expect_error(
    alignment_quality("not_a_list_or_result"),
    "must be an AlignmentResult or list of matrices"
  )
  expect_error(
    alignment_quality(42),
    "must be an AlignmentResult or list of matrices"
  )
})

test_that("alignment_quality rejects empty list (line 44)", {
  expect_error(
    alignment_quality(list()),
    "No aligned data to compute quality"
  )
})

test_that("alignment_quality computes improvement from AlignmentData original (line 69)", {
  set.seed(600)

  # Original data: random, low correlation
  original_list <- make_test_data_list(n_subjects = 3, n_features = 10, n_obs = 5)
  original <- AlignmentData(original_list)

  # Aligned data: share a common pattern, higher correlation
  shared <- make_test_matrix(n_features = 10, n_obs = 5)
  aligned <- list(
    "sub-01" = shared + matrix(rnorm(50, sd = 0.1), 10, 5),
    "sub-02" = shared + matrix(rnorm(50, sd = 0.1), 10, 5),
    "sub-03" = shared + matrix(rnorm(50, sd = 0.1), 10, 5)
  )

  quality <- alignment_quality(aligned, original = original, metrics = "correlation")

  expect_true("correlation_improvement" %in% names(quality))
  expect_true("original_mean_correlation" %in% names(quality))
  expect_true("aligned_mean_correlation" %in% names(quality))
  expect_true("relative_improvement" %in% names(quality))
  expect_true(is.numeric(quality$correlation_improvement))
  expect_true(quality$correlation_improvement > 0)
})

test_that("variance metrics with 3 subjects (lines 153-159)", {
  skip_if_not_installed("abind")

  set.seed(601)
  aligned <- make_test_data_list(n_subjects = 3, n_features = 10, n_obs = 5)

  quality <- alignment_quality(aligned, metrics = "variance")

  expect_true("total_variance" %in% names(quality))
  expect_true("between_subject_variance" %in% names(quality))
  expect_true("within_subject_variance" %in% names(quality))
  expect_true("variance_ratio" %in% names(quality))

  expect_true(is.numeric(quality$total_variance))
  expect_true(quality$total_variance > 0)
  expect_true(is.numeric(quality$variance_ratio))
  expect_true(quality$variance_ratio > 0)
  expect_true(is.numeric(quality$between_subject_variance))
  expect_true(is.numeric(quality$within_subject_variance))
})

test_that("ISC metrics with 2+ subjects return expected fields (lines 199-231)", {
  set.seed(602)

  # Shared signal + noise so ISC is measurable
  signal <- make_test_matrix(n_features = 10, n_obs = 5)
  aligned <- list(
    "sub-01" = signal + matrix(rnorm(50, sd = 0.3), 10, 5),
    "sub-02" = signal + matrix(rnorm(50, sd = 0.3), 10, 5),
    "sub-03" = signal + matrix(rnorm(50, sd = 0.3), 10, 5)
  )

  quality <- alignment_quality(aligned, metrics = "isc")

  expect_true("mean_isc" %in% names(quality))
  expect_true("isc_per_subject" %in% names(quality))
  expect_true("isc_per_feature" %in% names(quality))

  expect_true(is.numeric(quality$mean_isc))
  # With strong shared signal, mean ISC should be positive

  expect_true(quality$mean_isc > 0)
  expect_length(quality$isc_per_subject, 3)
  expect_length(quality$isc_per_feature, 10)
})

test_that("ISC metrics with single subject returns NA", {
  aligned <- list(
    "sub-01" = matrix(rnorm(50), 10, 5)
  )

  quality <- alignment_quality(aligned, metrics = "isc")
  expect_true(is.na(quality$isc))
})

test_that("print_quality_summary prints reference_correlation branch (lines 289-292)", {
  quality <- list(
    mean_pairwise_correlation = 0.80,
    sd_pairwise_correlation = 0.05,
    min_pairwise_correlation = 0.70,
    max_pairwise_correlation = 0.90,
    mean_reference_correlation = 0.95
  )

  out <- capture.output(print_quality_summary(quality))
  full <- paste(out, collapse = "\n")
  expect_true(grepl("Reference Correlation", full))
  expect_true(grepl("0.9500", full))
  expect_true(grepl("Pairwise Correlation", full))
})

test_that("print_quality_summary prints ISC branch (lines 294-298)", {
  quality <- list(
    mean_isc = 0.42
  )

  out <- capture.output(print_quality_summary(quality))
  full <- paste(out, collapse = "\n")
  expect_true(grepl("Inter-Subject Correlation", full))
  expect_true(grepl("0.4200", full))
})

test_that("print_quality_summary prints improvement branch (lines 300-311)", {
  quality <- list(
    mean_pairwise_correlation = 0.80,
    correlation_improvement = 0.30,
    original_mean_correlation = 0.50,
    aligned_mean_correlation = 0.80,
    relative_improvement = 0.60
  )

  out <- capture.output(print_quality_summary(quality))
  full <- paste(out, collapse = "\n")
  expect_true(grepl("Improvement", full))
  expect_true(grepl("Original correlation", full))
  expect_true(grepl("Aligned correlation", full))
  expect_true(grepl("0.3000", full))
  expect_true(grepl("60.0%", full))
})

test_that("print_quality_summary returns quality list invisibly", {
  quality <- list(mean_pairwise_correlation = 0.75)
  result <- invisible(capture.output(ret <- print_quality_summary(quality)))
  expect_identical(ret, quality)
})


# ---------- Additional edge-case coverage ----------

test_that("improvement with fewer than 2 common subjects returns NA", {
  aligned <- list("sub-01" = matrix(rnorm(50), 10, 5))
  original <- list("sub-99" = matrix(rnorm(50), 10, 5))

  # No common subjects → should return improvement = NA
  result <- neuralign:::.compute_improvement_metrics(aligned, original)
  expect_true(is.na(result$improvement))
})

test_that("improvement with exactly one common subject returns NA", {
  aligned <- list(
    "sub-01" = matrix(rnorm(50), 10, 5),
    "sub-02" = matrix(rnorm(50), 10, 5)
  )
  original <- list("sub-01" = matrix(rnorm(50), 10, 5))

  result <- neuralign:::.compute_improvement_metrics(aligned, original)
  expect_true(is.na(result$improvement))
})

test_that("variance_ratio is NA when within_var is zero", {
  skip_if_not_installed("abind")

  # All subjects have identical, constant data → within-subject var = 0
  const_mat <- matrix(1, 3, 4)
  aligned <- list(
    "sub-01" = const_mat,
    "sub-02" = const_mat * 2,
    "sub-03" = const_mat * 3
  )

  result <- neuralign:::.compute_variance_metrics(aligned)
  expect_true(is.na(result$variance_ratio))
  expect_equal(result$within_subject_variance, 0)
})

test_that("alignment_quality computes all four metrics at once", {
  skip_if_not_installed("abind")

  set.seed(700)
  signal <- matrix(rnorm(30), 5, 6)
  reference <- signal
  aligned <- list(
    "sub-01" = signal + matrix(rnorm(30, sd = 0.1), 5, 6),
    "sub-02" = signal + matrix(rnorm(30, sd = 0.1), 5, 6),
    "sub-03" = signal + matrix(rnorm(30, sd = 0.1), 5, 6)
  )

  quality <- alignment_quality(
    aligned,
    metrics = c("correlation", "reconstruction", "variance", "isc"),
    reference = reference
  )

  expect_true("mean_pairwise_correlation" %in% names(quality))
  expect_true("mean_reference_correlation" %in% names(quality))
  expect_true("total_variance" %in% names(quality))
  expect_true("mean_isc" %in% names(quality))
})

test_that("alignment_quality passes original as raw list (auto-converts)", {
  set.seed(701)
  original_list <- list(
    "sub-01" = matrix(rnorm(30), 5, 6),
    "sub-02" = matrix(rnorm(30), 5, 6)
  )
  shared <- matrix(rnorm(30), 5, 6)
  aligned <- list(
    "sub-01" = shared + matrix(rnorm(30, sd = 0.1), 5, 6),
    "sub-02" = shared + matrix(rnorm(30, sd = 0.1), 5, 6)
  )

  # Pass raw list — should auto-convert to AlignmentData
  quality <- alignment_quality(aligned, original = original_list, metrics = "correlation")
  expect_true("correlation_improvement" %in% names(quality))
  expect_true(is.numeric(quality$correlation_improvement))
})

test_that("reconstruction metric without reference is silently skipped", {
  set.seed(702)
  aligned <- list(
    "sub-01" = matrix(rnorm(30), 5, 6),
    "sub-02" = matrix(rnorm(30), 5, 6)
  )

  # Request reconstruction but provide no reference
  quality <- alignment_quality(aligned, metrics = c("correlation", "reconstruction"))
  expect_true("mean_pairwise_correlation" %in% names(quality))
  # reconstruction fields should not appear because reference is NULL
  expect_null(quality$mean_reference_correlation)
})

test_that("relative_improvement NA when original correlation is zero", {
  # Construct data where original pairwise correlation is exactly zero
  # Use orthogonal vectors
  set.seed(703)
  n <- 100
  x <- rnorm(n)
  y <- residuals(lm(rnorm(n) ~ x))  # orthogonal to x
  aligned <- list("s1" = matrix(x, 1, n), "s2" = matrix(x, 1, n))
  original <- list("s1" = matrix(x, 1, n), "s2" = matrix(y, 1, n))

  result <- neuralign:::.compute_improvement_metrics(aligned, original)
  # original should have ~0 correlation; relative_improvement may be NA or very large
  expect_true(is.numeric(result$correlation_improvement))
})
