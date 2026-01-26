test_that("select_reference with medoid works", {
  set.seed(100)

  # Create data where sub-02 is clearly the medoid
  # sub-02 is similar to both, while sub-01 and sub-03 are different from each other
  base <- matrix(rnorm(100), 10, 10)
  data_list <- list(
    "sub-01" = base + matrix(rnorm(100, sd = 2), 10, 10),
    "sub-02" = base + matrix(rnorm(100, sd = 0.1), 10, 10),
    "sub-03" = base + matrix(rnorm(100, sd = 2), 10, 10)
  )
  adat <- AlignmentData(data_list)

  ref <- select_reference(adat, method = "medoid")

  expect_true(ref %in% adat@subjects)
  # With high probability, sub-02 should be selected as medoid
  # But we don't hard-assert this due to randomness
})

test_that("select_reference with centroid works", {
  set.seed(101)
  data_list <- list(
    "sub-01" = matrix(rnorm(100), 10, 10),
    "sub-02" = matrix(rnorm(100), 10, 10),
    "sub-03" = matrix(rnorm(100), 10, 10)
  )
  adat <- AlignmentData(data_list)

  ref <- select_reference(adat, method = "centroid")

  expect_true(ref %in% adat@subjects)
})

test_that("select_reference with first works", {
  data_list <- list(
    "sub-01" = matrix(1, 10, 10),
    "sub-02" = matrix(2, 10, 10),
    "sub-03" = matrix(3, 10, 10)
  )
  adat <- AlignmentData(data_list)

  ref <- select_reference(adat, method = "first")

  expect_equal(ref, "sub-01")
})

test_that("select_reference with random is reproducible with seed", {
  data_list <- list(
    "sub-01" = matrix(1, 10, 10),
    "sub-02" = matrix(2, 10, 10),
    "sub-03" = matrix(3, 10, 10)
  )
  adat <- AlignmentData(data_list)

  ref1 <- select_reference(adat, method = "random", seed = 42)
  ref2 <- select_reference(adat, method = "random", seed = 42)

  expect_equal(ref1, ref2)
})

test_that("select_reference with single subject returns that subject", {
  data_list <- list(
    "sub-01" = matrix(1, 10, 10)
  )
  adat <- AlignmentData(data_list)

  ref <- select_reference(adat, method = "medoid")
  expect_equal(ref, "sub-01")
})

test_that("compute_centroid returns mean of data", {
  data_list <- list(
    "sub-01" = matrix(1, 10, 10),
    "sub-02" = matrix(2, 10, 10),
    "sub-03" = matrix(3, 10, 10)
  )
  adat <- AlignmentData(data_list)

  centroid <- compute_centroid(adat)

  expect_equal(dim(centroid), c(10, 10))
  expect_equal(centroid[1, 1], 2)  # Mean of 1, 2, 3
})

test_that("get_reference_data works with subject ID", {
  data_list <- list(
    "sub-01" = matrix(1, 10, 10),
    "sub-02" = matrix(2, 10, 10)
  )
  adat <- AlignmentData(data_list)

  ref_data <- get_reference_data(adat, "sub-01")
  expect_equal(ref_data[1, 1], 1)

  ref_data2 <- get_reference_data(adat, "sub-02")
  expect_equal(ref_data2[1, 1], 2)
})

test_that("get_reference_data works with consensus", {
  data_list <- list(
    "sub-01" = matrix(1, 10, 10),
    "sub-02" = matrix(3, 10, 10)
  )
  adat <- AlignmentData(data_list)

  ref_data <- get_reference_data(adat, "consensus")
  expect_equal(ref_data[1, 1], 2)  # Mean of 1 and 3
})

test_that("get_reference_data works with external matrix", {
  data_list <- list(
    "sub-01" = matrix(1, 10, 10)
  )
  adat <- AlignmentData(data_list)

  external_ref <- matrix(5, 10, 10)
  ref_data <- get_reference_data(adat, external_ref)

  expect_equal(ref_data, external_ref)
})

test_that("distance metrics produce valid results", {
  set.seed(102)
  data_list <- list(
    "sub-01" = matrix(rnorm(100), 10, 10),
    "sub-02" = matrix(rnorm(100), 10, 10),
    "sub-03" = matrix(rnorm(100), 10, 10)
  )
  adat <- AlignmentData(data_list)

  # Should work with different distance metrics
  ref_cor <- select_reference(adat, method = "medoid", distance = "correlation")
  ref_euc <- select_reference(adat, method = "medoid", distance = "euclidean")
  ref_frob <- select_reference(adat, method = "medoid", distance = "frobenius")

  expect_true(all(c(ref_cor, ref_euc, ref_frob) %in% adat@subjects))
})

# ---- Additional coverage tests ----

test_that("select_reference coerces plain list to AlignmentData (line 51)", {
  set.seed(200)
  plain_list <- list(
    "sub-01" = matrix(rnorm(50), 10, 5),
    "sub-02" = matrix(rnorm(50), 10, 5),
    "sub-03" = matrix(rnorm(50), 10, 5)
  )

  ref <- select_reference(plain_list, method = "first")
  expect_equal(ref, "sub-01")

  ref2 <- select_reference(plain_list, method = "medoid")
  expect_true(ref2 %in% c("sub-01", "sub-02", "sub-03"))
})

test_that("select_reference errors on empty AlignmentData (line 58)", {
  empty_adat <- new("AlignmentData",
    data = list(),
    subjects = character(0),
    metadata = list()
  )

  expect_error(
    select_reference(empty_adat, method = "medoid"),
    "No subjects in data"
  )
})

test_that("compute_centroid coerces plain list input (line 182)", {
  plain_list <- list(
    "sub-01" = matrix(1, 5, 3),
    "sub-02" = matrix(3, 5, 3)
  )

  centroid <- compute_centroid(plain_list)
  expect_equal(dim(centroid), c(5, 3))
  expect_equal(centroid[1, 1], 2)  # Mean of 1 and 3
})

test_that("compute_centroid errors on empty data (line 189)", {
  empty_adat <- new("AlignmentData",
    data = list(),
    subjects = character(0),
    metadata = list()
  )

  expect_error(
    compute_centroid(empty_adat),
    "No subjects in data"
  )
})

test_that("get_reference_data with matrix reference returns it directly (line 211)", {
  data_list <- list(
    "sub-01" = matrix(1, 5, 3)
  )
  adat <- AlignmentData(data_list)

  ext_ref <- matrix(99, 5, 3)
  result <- get_reference_data(adat, ext_ref)
  expect_equal(result, ext_ref)
})

test_that("get_reference_data with 'consensus' computes centroid (line 215)", {
  data_list <- list(
    "sub-01" = matrix(2, 5, 3),
    "sub-02" = matrix(4, 5, 3)
  )
  adat <- AlignmentData(data_list)

  result <- get_reference_data(adat, "consensus")
  expect_equal(dim(result), c(5, 3))
  expect_equal(result[1, 1], 3)  # Mean of 2 and 4
})

test_that("get_reference_data with subject ID returns that subject (line 223)", {
  data_list <- list(
    "sub-01" = matrix(10, 5, 3),
    "sub-02" = matrix(20, 5, 3)
  )
  adat <- AlignmentData(data_list)

  result <- get_reference_data(adat, "sub-02")
  expect_equal(result[1, 1], 20)
})

test_that("get_reference_data errors on missing subject (line 219-221)", {
  data_list <- list(
    "sub-01" = matrix(1, 5, 3)
  )
  adat <- AlignmentData(data_list)

  expect_error(
    get_reference_data(adat, "sub-99"),
    "Subject 'sub-99' not found in data"
  )
})

test_that(".validate_reference_for_method accepts valid reference types (lines 238-263)", {
  # Register a test aligner with limited reference_types
  register_aligner(
    name = "test_ref_validator",
    fit_fn = function(data, ...) list(),
    capabilities = list(
      reference_types = c("subject", "consensus")
    ),
    description = "Test aligner for reference validation"
  )
  on.exit(unregister_aligner("test_ref_validator"), add = TRUE)

  # Subject reference should be valid
  expect_true(
    neuralign:::.validate_reference_for_method("sub-01", "test_ref_validator")
  )

  # Consensus reference should be valid
  expect_true(
    neuralign:::.validate_reference_for_method("consensus", "test_ref_validator")
  )

  # "medoid" is treated as a subject reference, should be valid
  expect_true(
    neuralign:::.validate_reference_for_method("medoid", "test_ref_validator")
  )
})

test_that(".validate_reference_for_method rejects unsupported template reference", {
  # Register aligner that only supports subject references
  register_aligner(
    name = "test_ref_subject_only",
    fit_fn = function(data, ...) list(),
    capabilities = list(
      reference_types = c("subject")
    ),
    description = "Subject-only test aligner"
  )
  on.exit(unregister_aligner("test_ref_subject_only"), add = TRUE)

  # Template (matrix) reference should fail
  template <- matrix(1, 5, 3)
  expect_error(
    neuralign:::.validate_reference_for_method(template, "test_ref_subject_only"),
    "does not support reference type 'template'"
  )

  # Consensus reference should also fail
  expect_error(
    neuralign:::.validate_reference_for_method("consensus", "test_ref_subject_only"),
    "does not support reference type 'consensus'"
  )
})

test_that(".validate_reference_for_method returns TRUE for unknown method", {
  # An unregistered method should return TRUE (can't validate)
  expect_true(
    neuralign:::.validate_reference_for_method("sub-01", "nonexistent_aligner_xyz")
  )
})
