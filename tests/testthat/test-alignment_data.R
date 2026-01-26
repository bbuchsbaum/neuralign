test_that("AlignmentData creation works", {
  data_list <- list(
    "sub-01" = matrix(rnorm(100), 10, 10),
    "sub-02" = matrix(rnorm(100), 10, 10),
    "sub-03" = matrix(rnorm(100), 10, 10)
  )

  adat <- AlignmentData(data_list)

  expect_s4_class(adat, "AlignmentData")
  expect_equal(length(adat), 3)
  expect_equal(adat@subjects, c("sub-01", "sub-02", "sub-03"))
})

test_that("AlignmentData auto-generates subject IDs", {
  data_list <- list(
    matrix(rnorm(100), 10, 10),
    matrix(rnorm(100), 10, 10)
  )

  adat <- AlignmentData(data_list)

  expect_equal(length(adat), 2)
  expect_equal(adat@subjects, c("sub-01", "sub-02"))
})

test_that("AlignmentData rejects duplicate subjects", {
  data_list <- list(
    "sub-01" = matrix(1, 10, 10),
    "sub-02" = matrix(1, 10, 10)
  )

  expect_error(
    AlignmentData(data_list, subjects = c("sub-01", "sub-01")),
    "unique"
  )
})

test_that("AlignmentData subsetting works", {
  data_list <- list(
    "sub-01" = matrix(1, 10, 10),
    "sub-02" = matrix(2, 10, 10),
    "sub-03" = matrix(3, 10, 10)
  )
  adat <- AlignmentData(data_list)

  # Subset by index
  subset1 <- adat[1:2]
  expect_equal(length(subset1), 2)
  expect_equal(subset1@subjects, c("sub-01", "sub-02"))

  # Subset by name
  subset2 <- adat[c("sub-02", "sub-03")]
  expect_equal(length(subset2), 2)
  expect_equal(subset2@subjects, c("sub-02", "sub-03"))

  # Invalid subset
  expect_error(adat["sub-99"], "Unknown")
})

test_that("get_subject_data works", {
  data_list <- list(
    "sub-01" = matrix(1, 10, 10),
    "sub-02" = matrix(2, 10, 10)
  )
  adat <- AlignmentData(data_list)

  d1 <- get_subject_data(adat, "sub-01")
  expect_equal(d1[1, 1], 1)

  d2 <- get_subject_data(adat, 2)
  expect_equal(d2[1, 1], 2)

  expect_error(get_subject_data(adat, "sub-99"), "Unknown")
})

test_that("validate_alignment_data works", {
  # Valid data with same features
  data_list <- list(
    "sub-01" = matrix(1, 10, 5),
    "sub-02" = matrix(1, 10, 5)
  )
  adat <- AlignmentData(data_list)
  expect_true(validate_alignment_data(adat))

  # Invalid: different features
  data_list2 <- list(
    "sub-01" = matrix(1, 10, 5),
    "sub-02" = matrix(1, 8, 5)
  )
  adat2 <- AlignmentData(data_list2)
  expect_error(
    validate_alignment_data(adat2, check_features = TRUE),
    "different"
  )
})

test_that("obs_labels are stored and validated when provided", {
  data_list <- list(
    "sub-01" = matrix(1, 10, 5),
    "sub-02" = matrix(1, 10, 5)
  )

  labs <- c("A", "B", "C", "D", "E")
  adat <- AlignmentData(data_list, obs_labels = labs)

  expect_equal(get_obs_labels(adat), labs)
  expect_true(validate_alignment_data(adat))

  # Length mismatch triggers error
  adat_bad <- AlignmentData(data_list, obs_labels = c("A", "B"))
  expect_error(validate_alignment_data(adat_bad), "length mismatch")
})

test_that("obs_labels can be validated against column names", {
  m1 <- matrix(1, 10, 3)
  m2 <- matrix(1, 10, 3)
  colnames(m1) <- c("x", "y", "z")
  colnames(m2) <- c("x", "y", "z")

  adat <- AlignmentData(list("sub-01" = m1, "sub-02" = m2), obs_labels = c("x", "y", "z"))
  expect_true(validate_alignment_data(adat))

  # Mismatched colnames vs obs_labels
  colnames(m2) <- c("x", "y", "ZZ")
  adat_bad <- AlignmentData(list("sub-01" = m1, "sub-02" = m2), obs_labels = c("x", "y", "z"))
  expect_error(validate_alignment_data(adat_bad), "colnames do not match")
})

test_that("as_alignment_data coercion works", {
  data_list <- list(
    "sub-01" = matrix(1, 5, 5),
    "sub-02" = matrix(2, 5, 5)
  )

  adat <- as_alignment_data(data_list)
  expect_s4_class(adat, "AlignmentData")

  # Identity for AlignmentData
  adat2 <- as_alignment_data(adat)
  expect_identical(adat, adat2)
})
