test_that("fit_alignment works with procrustes", {
  # Register procrustes for testing
  neuralign:::.register_procrustes()

  # Create test data
  set.seed(123)
  n_features <- 20
  n_obs <- 30

  # Create data with shared structure
  reference_pattern <- matrix(rnorm(n_features * n_obs), n_features, n_obs)

  data_list <- lapply(1:4, function(i) {
    # Add subject-specific noise (procrustes will align in feature space)
    reference_pattern + matrix(rnorm(n_features * n_obs, sd = 0.5), n_features, n_obs)
  })
  names(data_list) <- paste0("sub-0", 1:4)

  adat <- AlignmentData(data_list)

  # Fit alignment
  result <- fit_alignment(adat, method = "procrustes", reference = "consensus")

  expect_s4_class(result, "AlignmentResult")
  expect_equal(length(result@aligned), 4)
  expect_equal(names(result@aligned), names(data_list))

  # Check model
  model <- get_model(result)
  expect_s4_class(model, "AlignmentModel")
  expect_equal(model@method, "procrustes")
})

test_that("fit_alignment with medoid reference works", {
  neuralign:::.register_procrustes()

  set.seed(456)
  data_list <- list(
    "sub-01" = matrix(rnorm(100), 10, 10),
    "sub-02" = matrix(rnorm(100), 10, 10),
    "sub-03" = matrix(rnorm(100), 10, 10)
  )
  adat <- AlignmentData(data_list)

  result <- fit_alignment(adat, method = "procrustes", reference = "medoid")

  expect_s4_class(result, "AlignmentResult")
  expect_true(result@model@reference %in% adat@subjects)
})

test_that("fit_alignment with specific reference subject works", {
  neuralign:::.register_procrustes()

  set.seed(789)
  data_list <- list(
    "sub-01" = matrix(rnorm(100), 10, 10),
    "sub-02" = matrix(rnorm(100), 10, 10),
    "sub-03" = matrix(rnorm(100), 10, 10)
  )
  adat <- AlignmentData(data_list)

  result <- fit_alignment(adat, method = "procrustes", reference = "sub-02")

  expect_equal(result@model@reference, "sub-02")
})

test_that("fit_alignment errors on unknown method", {
  data_list <- list(
    "sub-01" = matrix(1, 10, 10),
    "sub-02" = matrix(1, 10, 10)
  )
  adat <- AlignmentData(data_list)

  expect_error(
    fit_alignment(adat, method = "nonexistent_method"),
    "Unknown method"
  )
})

test_that("fit_alignment produces transforms with correct dimensions", {
  neuralign:::.register_procrustes()

  n_features <- 15
  data_list <- list(
    "sub-01" = matrix(rnorm(n_features * 20), n_features, 20),
    "sub-02" = matrix(rnorm(n_features * 20), n_features, 20),
    "sub-03" = matrix(rnorm(n_features * 20), n_features, 20)
  )
  adat <- AlignmentData(data_list)

  result <- fit_alignment(adat, method = "procrustes")
  model <- get_model(result)

  for (subj in names(model@transforms)) {
    transform <- model@transforms[[subj]]
    # Transforms should be (target x source) = (n_features x n_features)
    expect_equal(dim(transform), c(n_features, n_features))
  }
})

test_that("fit_alignment with train_idx works", {
  neuralign:::.register_procrustes()

  set.seed(111)
  data_list <- list(
    "sub-01" = matrix(rnorm(100), 10, 10),
    "sub-02" = matrix(rnorm(100), 10, 10),
    "sub-03" = matrix(rnorm(100), 10, 10),
    "sub-04" = matrix(rnorm(100), 10, 10)
  )
  adat <- AlignmentData(data_list)

  # Fit on first 3 subjects only
  result <- fit_alignment(
    adat,
    method = "procrustes",
    train_idx = 1:3
  )

  # Should still produce transforms for all subjects
  expect_equal(length(result@model@transforms), 4)

  # Train subjects should be recorded
  expect_equal(result@model@train_subjects, c("sub-01", "sub-02", "sub-03"))
})

test_that("fit_alignment computes quality metrics", {
  neuralign:::.register_procrustes()

  set.seed(222)
  data_list <- list(
    "sub-01" = matrix(rnorm(100), 10, 10),
    "sub-02" = matrix(rnorm(100), 10, 10),
    "sub-03" = matrix(rnorm(100), 10, 10)
  )
  adat <- AlignmentData(data_list)

  result <- fit_alignment(adat, method = "procrustes", compute_quality = TRUE)

  expect_true(length(result@quality) > 0)
  expect_true("mean_pairwise_correlation" %in% names(result@quality))
})
