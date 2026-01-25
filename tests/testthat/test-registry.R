test_that("register_aligner works", {
  # Clean registry first
  neuralign:::.clear_registry()

  # Register a test aligner
  test_fit <- function(data, reference, train_idx = NULL, ...) {
    n <- length(data@subjects)
    transforms <- lapply(data@subjects, function(s) diag(10))
    names(transforms) <- data@subjects
    list(
      transforms = transforms,
      reference_data = matrix(0, 10, 10),
      space_from = NULL,
      space_to = NULL
    )
  }

  register_aligner(
    name = "test_method",
    fit_fn = test_fit,
    capabilities = list(supports_cv = TRUE),
    package = "testpkg",
    description = "Test method"
  )

  expect_true(is_aligner_registered("test_method"))
  expect_false(is_aligner_registered("nonexistent"))
})

test_that("available_aligners returns registered methods", {
  neuralign:::.clear_registry()

  # Register some aligners
  dummy_fit <- function(data, reference, train_idx = NULL, ...) {
    list(transforms = list(), reference_data = NULL, space_from = NULL, space_to = NULL)
  }

  register_aligner("method_a", dummy_fit, package = "pkg_a")
  register_aligner("method_b", dummy_fit, package = "pkg_b")

  aligners <- available_aligners()
  expect_true("method_a" %in% aligners)
  expect_true("method_b" %in% aligners)

  # With details
  details <- available_aligners(details = TRUE)
  expect_s3_class(details, "data.frame")
  expect_true("name" %in% names(details))
  expect_true("package" %in% names(details))
})

test_that("get_aligner returns aligner info", {
  neuralign:::.clear_registry()

  dummy_fit <- function(data, reference, train_idx = NULL, ...) {
    list(transforms = list(), reference_data = NULL, space_from = NULL, space_to = NULL)
  }

  register_aligner(
    "my_aligner",
    dummy_fit,
    capabilities = list(transform_type = "orthogonal"),
    description = "My aligner"
  )

  info <- get_aligner("my_aligner")
  expect_type(info, "list")
  expect_equal(info$name, "my_aligner")
  expect_equal(info$description, "My aligner")

  # Non-existent
  expect_null(get_aligner("nonexistent"))
})

test_that("aligner_capabilities returns capability info", {
  neuralign:::.clear_registry()

  dummy_fit <- function(data, reference, train_idx = NULL, ...) {
    list(transforms = list(), reference_data = NULL, space_from = NULL, space_to = NULL)
  }

  register_aligner(
    "cap_test",
    dummy_fit,
    capabilities = list(
      supports_cv = TRUE,
      needs_geometry = TRUE,
      transform_type = "linear"
    )
  )

  caps <- aligner_capabilities("cap_test")
  expect_type(caps, "list")
  expect_true(caps$supports_cv)
  expect_true(caps$needs_geometry)
  expect_equal(caps$transform_type, "linear")
})

test_that("unregister_aligner removes aligner", {
  neuralign:::.clear_registry()

  dummy_fit <- function(data, reference, train_idx = NULL, ...) {
    list(transforms = list(), reference_data = NULL, space_from = NULL, space_to = NULL)
  }

  register_aligner("to_remove", dummy_fit)
  expect_true(is_aligner_registered("to_remove"))

  result <- unregister_aligner("to_remove")
  expect_true(result)
  expect_false(is_aligner_registered("to_remove"))

  # Removing non-existent returns FALSE
  result2 <- unregister_aligner("nonexistent")
  expect_false(result2)
})

test_that("default capabilities are set correctly", {
  neuralign:::.clear_registry()

  dummy_fit <- function(data, reference, train_idx = NULL, ...) {
    list(transforms = list(), reference_data = NULL, space_from = NULL, space_to = NULL)
  }

  # Register with minimal capabilities
  register_aligner("minimal", dummy_fit)

  caps <- aligner_capabilities("minimal")

  # Check defaults
  expect_false(caps$supports_cv)  # Default is FALSE
  expect_false(caps$needs_geometry)
  expect_equal(caps$transform_type, "linear")  # Default type
})

test_that("embedding-returning aligners are rejected with guidance", {
  neuralign:::.clear_registry()

  dummy_fit <- function(data, reference, train_idx = NULL, ...) {
    list(transforms = list(), reference_data = NULL, space_from = NULL, space_to = NULL)
  }

  expect_error(
    register_aligner(
      "emb",
      dummy_fit,
      capabilities = list(returns = "embedding")
    ),
    "operator-returning"
  )
})
