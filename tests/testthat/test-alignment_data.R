test_that("AlignmentData creation works", {
  data_list <- make_test_data_list(n_subjects = 3, n_features = 10, n_obs = 10)

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

# ---- Constructor error / warning paths ----

test_that("AlignmentData errors when data is not a list", {
  expect_error(
    AlignmentData(matrix(1, 5, 5)),
    "must be a list"
  )
  expect_error(
    AlignmentData("not a list"),
    "must be a list"
  )
})

test_that("AlignmentData errors when subjects length mismatches data length", {
  data_list <- list(matrix(1, 5, 5), matrix(2, 5, 5))
  expect_error(
    AlignmentData(data_list, subjects = c("a", "b", "c")),
    "must match length"
  )
  expect_error(
    AlignmentData(data_list, subjects = "a"),
    "must match length"
  )
})

test_that("AlignmentData warns for non-matrix, non-NeuroVec elements", {
  data_list <- list(
    "sub-01" = data.frame(x = 1:5, y = 6:10),
    "sub-02" = matrix(1, 5, 2)
  )
  expect_warning(
    AlignmentData(data_list),
    "not a matrix or NeuroVec"
  )
})

# ---- show() method branches ----

test_that("show() displays truncated IDs when >6 subjects", {
  mats <- setNames(
    lapply(1:8, function(i) matrix(rnorm(20), 5, 4)),
    paste0("sub-", sprintf("%02d", 1:8))
  )
  adat <- AlignmentData(mats)
  out <- capture.output(show(adat))
  combined <- paste(out, collapse = "\n")
  expect_match(combined, "5 more")
  expect_match(combined, "sub-01, sub-02, sub-03")
})

test_that("show() prints data class for non-matrix data", {
  data_list <- list("sub-01" = data.frame(x = 1:5))
  suppressWarnings(adat <- AlignmentData(data_list))
  out <- capture.output(show(adat))
  combined <- paste(out, collapse = "\n")
  expect_match(combined, "Data class: data.frame")
})

test_that("show() prints design when present", {
  data_list <- list("sub-01" = matrix(1, 5, 3), "sub-02" = matrix(2, 5, 3))
  adat <- AlignmentData(data_list, design = list(conditions = c("A", "B")))
  out <- capture.output(show(adat))
  combined <- paste(out, collapse = "\n")
  expect_match(combined, "Design: present")
})

test_that("show() prints geometry when present", {
  data_list <- list("sub-01" = matrix(1, 5, 3), "sub-02" = matrix(2, 5, 3))
  adj <- matrix(1, 5, 5)
  adat <- AlignmentData(data_list, geometry = adj)
  out <- capture.output(show(adat))
  combined <- paste(out, collapse = "\n")
  expect_match(combined, "Geometry: present")
})

test_that("show() prints space when present", {
  data_list <- list("sub-01" = matrix(1, 5, 3))
  fake_space <- structure(list(), class = "FakeSpace")
  adat <- AlignmentData(data_list, space = fake_space)
  out <- capture.output(show(adat))
  combined <- paste(out, collapse = "\n")
  expect_match(combined, "Space: FakeSpace")
})

test_that("show() prints obs_labels (<=6 labels)", {
  data_list <- list("sub-01" = matrix(1, 5, 3), "sub-02" = matrix(2, 5, 3))
  labs <- c("A", "B", "C")
  adat <- AlignmentData(data_list, obs_labels = labs)
  out <- capture.output(show(adat))
  combined <- paste(out, collapse = "\n")
  expect_match(combined, "Observation labels: 3")
  expect_match(combined, "Labels: A, B, C")
})

test_that("show() truncates obs_labels when >6 labels", {
  data_list <- list("sub-01" = matrix(1, 5, 8), "sub-02" = matrix(2, 5, 8))
  labs <- paste0("lab", 1:8)
  adat <- AlignmentData(data_list, obs_labels = labs)
  out <- capture.output(show(adat))
  combined <- paste(out, collapse = "\n")
  expect_match(combined, "Observation labels: 8")
  expect_match(combined, "lab1, lab2, lab3")
  expect_match(combined, "5 more")
})

# ---- validate_alignment_data() edge cases ----

test_that("validate_alignment_data returns NA dims for non-matrix, non-NeuroVec data", {
  data_list <- list("sub-01" = list(1, 2, 3), "sub-02" = list(4, 5, 6))
  suppressWarnings(adat <- AlignmentData(data_list))
  # Features check should error because dims are NA (unique NA counts as 1 but

  # check_observations should still pass because unique(NA) has length 1)
  # The key thing: it should not crash, and with check_features=FALSE we get TRUE
  expect_true(validate_alignment_data(adat, check_features = FALSE))
})

test_that("validate_alignment_data check_obs_labels with NULL obs_labels errors", {
  data_list <- list("sub-01" = matrix(1, 5, 3), "sub-02" = matrix(1, 5, 3))
  adat <- AlignmentData(data_list)
  expect_error(
    validate_alignment_data(adat, check_obs_labels = TRUE),
    "obs_labels is NULL"
  )
})

test_that("validate_alignment_data errors for non-atomic obs_labels", {
  data_list <- list("sub-01" = matrix(1, 5, 3), "sub-02" = matrix(1, 5, 3))
  adat <- AlignmentData(data_list, obs_labels = list("a", "b", "c"))
  expect_error(
    validate_alignment_data(adat),
    "named list keyed by subject"
  )
})

test_that("validate_alignment_data supports per-subject obs_labels lists", {
  data_list <- list("sub-01" = matrix(1, 5, 3), "sub-02" = matrix(1, 5, 4))
  adat <- AlignmentData(
    data_list,
    obs_labels = list(
      "sub-01" = c("a", "b", "c"),
      "sub-02" = c("a", "b", "c", "d")
    )
  )
  expect_true(validate_alignment_data(adat, check_features = TRUE, check_observations = FALSE))
})

test_that("validate_alignment_data errors for obs_labels length mismatch", {
  data_list <- list("sub-01" = matrix(1, 5, 4), "sub-02" = matrix(1, 5, 4))
  adat <- AlignmentData(data_list, obs_labels = c("a", "b"))
  expect_error(
    validate_alignment_data(adat),
    "length mismatch"
  )
})

test_that("validate_alignment_data errors for obs_labels with NAs", {
  data_list <- list("sub-01" = matrix(1, 5, 3), "sub-02" = matrix(1, 5, 3))
  adat <- AlignmentData(data_list, obs_labels = c("a", NA, "c"))
  expect_error(
    validate_alignment_data(adat),
    "contains NA"
  )
})

test_that("validate_alignment_data errors when some subjects have colnames and others don't", {
  m1 <- matrix(1, 5, 3)
  colnames(m1) <- c("x", "y", "z")
  m2 <- matrix(1, 5, 3)
  # m2 has no colnames
  data_list <- list("sub-01" = m1, "sub-02" = m2)
  adat <- AlignmentData(data_list, obs_labels = c("x", "y", "z"))
  expect_error(
    validate_alignment_data(adat),
    "Some subjects have colnames but others do not"
  )
})

test_that("validate_alignment_data errors when all colnames present but don't match obs_labels", {
  m1 <- matrix(1, 5, 3)
  colnames(m1) <- c("a", "b", "c")
  m2 <- matrix(1, 5, 3)
  colnames(m2) <- c("a", "b", "c")
  data_list <- list("sub-01" = m1, "sub-02" = m2)
  adat <- AlignmentData(data_list, obs_labels = c("x", "y", "z"))
  expect_error(
    validate_alignment_data(adat),
    "colnames do not match"
  )
})

test_that("validate_alignment_data errors for different observation counts", {
  data_list <- list(
    "sub-01" = matrix(1, 5, 3),
    "sub-02" = matrix(1, 5, 4)
  )
  adat <- AlignmentData(data_list)
  expect_error(
    validate_alignment_data(adat, check_observations = TRUE),
    "different numbers of observations"
  )
})

test_that("validate_alignment_data errors for differing obs when obs_labels set", {
  # When obs_labels are set, check_observations is auto-enabled
  data_list <- list(
    "sub-01" = matrix(1, 5, 3),
    "sub-02" = matrix(1, 5, 4)
  )
  adat <- AlignmentData(data_list, obs_labels = c("a", "b", "c"))
  expect_error(
    validate_alignment_data(adat),
    "different numbers of observations"
  )
})


# ---------- list obs_labels validation ----------

test_that("validate_alignment_data per-subject obs_labels with non-atomic entry errors", {
  data_list <- list("sub-01" = matrix(1, 5, 3), "sub-02" = matrix(1, 5, 3))
  adat <- AlignmentData(data_list, obs_labels = list(
    "sub-01" = c("a", "b", "c"),
    "sub-02" = list("a", "b", "c")
  ))
  expect_error(validate_alignment_data(adat), "must be an atomic vector")
})

test_that("validate_alignment_data per-subject obs_labels with length mismatch errors", {
  data_list <- list("sub-01" = matrix(1, 5, 3), "sub-02" = matrix(1, 5, 3))
  adat <- AlignmentData(data_list, obs_labels = list(
    "sub-01" = c("a", "b", "c"),
    "sub-02" = c("a", "b")
  ))
  expect_error(validate_alignment_data(adat), "length mismatch")
})

test_that("validate_alignment_data per-subject obs_labels with NA errors", {
  data_list <- list("sub-01" = matrix(1, 5, 3), "sub-02" = matrix(1, 5, 3))
  adat <- AlignmentData(data_list, obs_labels = list(
    "sub-01" = c("a", "b", "c"),
    "sub-02" = c("a", NA, "c")
  ))
  expect_error(validate_alignment_data(adat), "NA values")
})

test_that("validate_alignment_data per-subject obs_labels with colnames mismatch errors", {
  m1 <- matrix(1, 5, 3)
  colnames(m1) <- c("a", "b", "c")
  m2 <- matrix(1, 5, 3)
  colnames(m2) <- c("x", "y", "z")
  data_list <- list("sub-01" = m1, "sub-02" = m2)
  adat <- AlignmentData(data_list, obs_labels = list(
    "sub-01" = c("a", "b", "c"),
    "sub-02" = c("a", "b", "c")
  ))
  expect_error(validate_alignment_data(adat), "colnames do not match")
})

test_that("validate_alignment_data per-subject obs_labels missing subject errors", {
  data_list <- list("sub-01" = matrix(1, 5, 3), "sub-02" = matrix(1, 5, 3))
  adat <- AlignmentData(data_list, obs_labels = list(
    "sub-01" = c("a", "b", "c")
  ))
  expect_error(validate_alignment_data(adat), "missing subjects")
})

test_that("validate_alignment_data with non-list non-atomic obs_labels errors", {
  data_list <- list("sub-01" = matrix(1, 5, 3))
  adat <- AlignmentData(data_list, obs_labels = environment())
  expect_error(
    validate_alignment_data(adat, check_obs_labels = TRUE),
    "atomic vector.*named list"
  )
})
