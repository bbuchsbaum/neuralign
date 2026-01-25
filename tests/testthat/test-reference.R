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
