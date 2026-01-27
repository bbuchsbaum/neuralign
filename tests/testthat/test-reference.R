test_that("select_reference with medoid works", {
  set.seed(100)

  # Create data where sub-02 is clearly the medoid
  # sub-02 is similar to both, while sub-01 and sub-03 are different from each other
  base <- make_test_matrix(n_features = 10, n_obs = 10)
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
  data_list <- make_test_data_list(n_subjects = 3, n_features = 10, n_obs = 10)
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
  data_list <- make_test_data_list(n_subjects = 3, n_features = 10, n_obs = 10)
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
  plain_list <- make_test_data_list(n_subjects = 3, n_features = 10, n_obs = 5)

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


# ---------- More reference coverage tests ----------

test_that("select_reference with procrustes distance metric", {
  ensure_test_aligner("procrustes")

  set.seed(200)
  d <- 5
  n_obs <- 8
  Z <- matrix(rnorm(d * n_obs), d, n_obs)
  Q <- qr.Q(qr(matrix(rnorm(d * d), d, d)))
  if (det(Q) < 0) Q[, 1] <- -Q[, 1]

  data_list <- list(
    "sub-01" = Z,
    "sub-02" = Q %*% Z + matrix(rnorm(d * n_obs, sd = 0.01), d, n_obs),
    "sub-03" = matrix(rnorm(d * n_obs), d, n_obs)
  )
  adat <- AlignmentData(data_list)

  ref <- select_reference(adat, method = "medoid", distance = "procrustes")
  expect_true(ref %in% adat@subjects)
})

test_that("select_reference centroid errors with per-subject obs_labels", {
  set.seed(201)
  data_list <- list(
    "sub-01" = matrix(rnorm(30), 5, 6),
    "sub-02" = matrix(rnorm(25), 5, 5)
  )
  adat <- AlignmentData(
    data_list,
    obs_labels = list(
      "sub-01" = paste0("obs", 1:6),
      "sub-02" = paste0("obs", 1:5)
    )
  )

  expect_error(
    select_reference(adat, method = "centroid"),
    "obs_labels differ across subjects"
  )
})

test_that("select_reference random without seed is nondeterministic", {
  data_list <- list(
    "sub-01" = matrix(1, 5, 3),
    "sub-02" = matrix(2, 5, 3),
    "sub-03" = matrix(3, 5, 3)
  )
  adat <- AlignmentData(data_list)

  # Without seed, just verify it returns a valid subject
  ref <- select_reference(adat, method = "random")
  expect_true(ref %in% adat@subjects)
})

test_that("get_reference_data coerces plain list to AlignmentData", {
  plain_list <- list("sub-01" = matrix(10, 5, 3))
  ref_data <- get_reference_data(plain_list, "sub-01")
  expect_equal(ref_data[1, 1], 10)
})

test_that(".resolve_obs_labels_by_subject handles unnamed list", {
  data_list <- list(
    "sub-01" = matrix(1, 5, 3),
    "sub-02" = matrix(2, 5, 3)
  )
  adat <- AlignmentData(data_list)

  # Manually set obs_labels as an unnamed list
  adat@obs_labels <- list(c("a", "b", "c"), c("d", "e", "f"))

  result <- neuralign:::.resolve_obs_labels_by_subject(adat)
  expect_equal(names(result), c("sub-01", "sub-02"))
  expect_equal(result[["sub-01"]], c("a", "b", "c"))
})

test_that(".resolve_obs_labels_by_subject errors on unnamed list length mismatch", {
  data_list <- list(
    "sub-01" = matrix(1, 5, 3),
    "sub-02" = matrix(2, 5, 3)
  )
  adat <- AlignmentData(data_list)

  # Set obs_labels as unnamed list with wrong length
  adat@obs_labels <- list(c("a", "b", "c"))

  expect_error(
    neuralign:::.resolve_obs_labels_by_subject(adat),
    "length must match"
  )
})

test_that(".resolve_obs_labels_by_subject errors on missing subjects in named list", {
  data_list <- list(
    "sub-01" = matrix(1, 5, 3),
    "sub-02" = matrix(2, 5, 3)
  )
  adat <- AlignmentData(data_list)

  adat@obs_labels <- list("sub-01" = c("a", "b", "c"))

  expect_error(
    neuralign:::.resolve_obs_labels_by_subject(adat),
    "missing subjects"
  )
})

test_that(".resolve_obs_labels_by_subject errors on non-list non-atomic", {
  data_list <- list("sub-01" = matrix(1, 5, 3))
  adat <- AlignmentData(data_list)

  adat@obs_labels <- environment()

  expect_error(
    neuralign:::.resolve_obs_labels_by_subject(adat),
    "atomic vector.*list"
  )
})

test_that(".match_obs_indices errors on insufficient overlap", {
  expect_error(
    neuralign:::.match_obs_indices(c("a", "b"), c("c", "d"), min_overlap = 1L),
    "Not enough shared"
  )
})

test_that(".subset_to_overlap errors on incompatible dimensions without labels", {
  expect_error(
    neuralign:::.subset_to_overlap(matrix(1, 2, 3), matrix(1, 3, 4)),
    "incompatible dimensions"
  )
})

test_that(".compute_pairwise_distance returns NA when overlap fails", {
  x <- matrix(rnorm(10), 2, 5)
  y <- matrix(rnorm(6), 2, 3)

  # Different obs_labels with no overlap
  d <- neuralign:::.compute_pairwise_distance(
    x, y, "correlation",
    obs_labels_x = paste0("a", 1:5),
    obs_labels_y = paste0("b", 1:3),
    min_overlap = 1L
  )
  expect_true(is.na(d))
})


# ---------- select_reference edge cases ----------

test_that("select_reference medoid errors when all pairwise distances are NA", {
  # When subjects have disjoint obs_labels, all off-diagonal distances are NA
  # and there is no meaningful medoid.
  adat <- AlignmentData(
    list(
      s1 = matrix(rnorm(6), 2, 3),
      s2 = matrix(rnorm(6), 2, 3)
    ),
    obs_labels = list(
      s1 = c("a1", "a2", "a3"),
      s2 = c("b1", "b2", "b3")
    )
  )

  expect_error(
    select_reference(adat, method = "medoid", distance = "procrustes"),
    "Cannot select medoid: no finite pairwise distances"
  )
})

test_that(".compute_pairwise_distance with euclidean metric", {
  set.seed(70)
  x <- matrix(rnorm(10), 2, 5)
  y <- matrix(rnorm(10), 2, 5)

  d <- neuralign:::.compute_pairwise_distance(x, y, "euclidean")
  expect_true(is.numeric(d))
  expect_true(d >= 0)
})

test_that(".compute_pairwise_distance with frobenius metric", {
  set.seed(71)
  x <- matrix(rnorm(10), 2, 5)
  y <- matrix(rnorm(10), 2, 5)

  d <- neuralign:::.compute_pairwise_distance(x, y, "frobenius")
  expect_true(is.numeric(d))
  expect_true(d >= 0)
  # Frobenius should equal sqrt of sum of squared differences
  expected <- sqrt(sum((x - y)^2))
  expect_equal(d, expected, tolerance = 1e-10)
})

test_that(".compute_pairwise_distance errors on unknown metric", {
  x <- matrix(1, 2, 3)
  y <- matrix(1, 2, 3)
  expect_error(
    neuralign:::.compute_pairwise_distance(x, y, "unknown_metric"),
    "Unknown distance"
  )
})
