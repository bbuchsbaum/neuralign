# Tests for fmrigds compatibility functions

.is_pkg_installed <- function(pkg) {
  nzchar(system.file(package = pkg))
}

.load_pkg_quietly_or_skip <- function(pkg) {
  if (!.is_pkg_installed(pkg)) {
    skip(paste0("{", pkg, "} is not installed"))
  }
  suppressWarnings(requireNamespace(pkg, quietly = TRUE))
  invisible(TRUE)
}

test_that("as_map_family requires fmrigds", {
  skip_if(.is_pkg_installed("fmrigds"), "fmrigds is available")

  transforms <- list("sub-01" = diag(5))
  model <- AlignmentModel(
    transforms = transforms,
    reference = "consensus",
    method = "procrustes",
    space_from = "native",
    space_to = "MNI"
  )

  expect_error(as_map_family(model), "fmrigds")
})

test_that("as_map_family requires space_from and space_to", {
  .load_pkg_quietly_or_skip("fmrigds")

  transforms <- list("sub-01" = diag(5))

  # Missing space_from
  model_no_from <- AlignmentModel(
    transforms = transforms,
    reference = "consensus",
    method = "procrustes",
    space_to = "MNI"
  )
  expect_error(as_map_family(model_no_from), "space_from")

  # Missing space_to
  model_no_to <- AlignmentModel(
    transforms = transforms,
    reference = "consensus",
    method = "procrustes",
    space_from = "native"
  )
  expect_error(as_map_family(model_no_to), "space_to")
})

test_that("as_map_family validates input type", {
  .load_pkg_quietly_or_skip("fmrigds")

  expect_error(as_map_family("not a model"), "must be an AlignmentModel")
})

test_that("as_map_family accepts AlignmentResult", {
  .load_pkg_quietly_or_skip("fmrigds")

  transforms <- list("sub-01" = diag(5))
  model <- AlignmentModel(
    transforms = transforms,
    reference = "consensus",
    method = "procrustes",
    space_from = "native",
    space_to = "MNI"
  )
  result <- AlignmentResult(
    model = model,
    aligned = list("sub-01" = matrix(1, 5, 3))
  )

  # Should extract model from result
  map_fam <- as_map_family(result)
  expect_true(!is.null(map_fam))
})

test_that("from_map_family creates AlignmentModel", {
  # Create a mock MapFamily-like object
  mock_map_family <- list(
    by_subject = list(
      "sub-01" = diag(10),
      "sub-02" = diag(10) * 2
    ),
    from = "native",
    to = "MNI"
  )

  model <- from_map_family(mock_map_family)

  expect_s4_class(model, "AlignmentModel")
  expect_equal(length(model), 2)
  expect_equal(model@method, "fmrigds_imported")
  expect_equal(model@space_from, "native")
  expect_equal(model@space_to, "MNI")
  expect_equal(model@train_subjects, c("sub-01", "sub-02"))
})

test_that("from_map_family accepts custom method name", {
  mock_map_family <- list(
    by_subject = list("sub-01" = diag(5)),
    from = "A",
    to = "B"
  )

  model <- from_map_family(mock_map_family, method = "custom_method")
  expect_equal(model@method, "custom_method")
})

test_that("from_map_family validates input", {
  # Not a list
  expect_error(from_map_family("not a list"), "MapFamily-like")

  # Missing by_subject
  expect_error(from_map_family(list(from = "A")), "by_subject")
})

test_that("apply_to_gds requires fmrigds", {
  skip_if(.is_pkg_installed("fmrigds"), "fmrigds is available")

  transforms <- list("sub-01" = diag(5))
  model <- AlignmentModel(
    transforms = transforms,
    reference = "consensus",
    method = "procrustes"
  )

  expect_error(apply_to_gds(model, NULL), "fmrigds")
})

test_that("alignment_data_from_gds requires fmrigds", {
  skip_if(.is_pkg_installed("fmrigds"), "fmrigds is available")

  expect_error(alignment_data_from_gds(NULL), "fmrigds")
})

# Integration tests with fmrigds
test_that("full fmrigds round-trip works", {
  .load_pkg_quietly_or_skip("fmrigds")

  # Create alignment model
  transforms <- list(
    "sub-01" = diag(10),
    "sub-02" = matrix(rnorm(100), 10, 10)
  )
  model <- AlignmentModel(
    transforms = transforms,
    reference = "consensus",
    method = "procrustes",
    space_from = "native",
    space_to = "MNI"
  )

  # Convert to MapFamily
  map_fam <- as_map_family(model, name = "test_alignment")

  # Convert back
  model2 <- from_map_family(list(
    by_subject = map_fam$by_subject,
    from = map_fam$from_space,
    to = map_fam$to_space
  ))

  # Check transforms preserved
  expect_equal(get_transform(model2, "sub-01"), diag(10))
})
