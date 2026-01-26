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


# ---------- New tests appended below ----------

test_that("register_aligner errors when name is not a single character string", {
  dummy_fit <- function(data, reference, train_idx = NULL, ...) {
    list(transforms = list(), reference_data = NULL, space_from = NULL, space_to = NULL)
  }

  # Numeric name

  expect_error(
    register_aligner(name = 42, fit_fn = dummy_fit),
    "must be a single character string"
  )

  # Character vector of length > 1
  expect_error(
    register_aligner(name = c("a", "b"), fit_fn = dummy_fit),
    "must be a single character string"
  )
})

test_that("register_aligner errors when fit_fn is not a function", {
  expect_error(
    register_aligner(name = "bad_fit", fit_fn = "not_a_function"),
    "must be a function"
  )
})

test_that("register_aligner errors when apply_fn is not a function or NULL", {
  dummy_fit <- function(data, reference, train_idx = NULL, ...) {
    list(transforms = list(), reference_data = NULL, space_from = NULL, space_to = NULL)
  }

  expect_error(
    register_aligner(name = "bad_apply", fit_fn = dummy_fit, apply_fn = 123),
    "must be a function or NULL"
  )
})

test_that("register_aligner errors when capabilities$returns is not a string", {
  dummy_fit <- function(data, reference, train_idx = NULL, ...) {
    list(transforms = list(), reference_data = NULL, space_from = NULL, space_to = NULL)
  }

  expect_error(
    register_aligner(
      "bad_returns_type",
      dummy_fit,
      capabilities = list(returns = 123)
    ),
    "must be a single string"
  )
})

test_that("register_aligner errors when capabilities$returns is an invalid string", {
  dummy_fit <- function(data, reference, train_idx = NULL, ...) {
    list(transforms = list(), reference_data = NULL, space_from = NULL, space_to = NULL)
  }

  expect_error(
    register_aligner(
      "bad_returns_value",
      dummy_fit,
      capabilities = list(returns = "xyz")
    ),
    "must be one of"
  )
})

test_that("available_aligners(details=TRUE) returns empty data.frame when registry is empty", {
  # Save current registry state
  saved <- as.list(neuralign:::.aligner_registry)
  neuralign:::.clear_registry()
  on.exit({
    # Restore registry
    neuralign:::.clear_registry()
    for (nm in names(saved)) {
      assign(nm, saved[[nm]], envir = neuralign:::.aligner_registry)
    }
  })

  result_details <- available_aligners(details = TRUE)
  expect_s3_class(result_details, "data.frame")
  expect_equal(nrow(result_details), 0)
  expect_true("name" %in% names(result_details))
  expect_true("package" %in% names(result_details))
  expect_true("description" %in% names(result_details))

  result_names <- available_aligners(details = FALSE)
  expect_equal(length(result_names), 0)
  expect_type(result_names, "character")
})

test_that(".try_autoload_aligner returns FALSE for known method whose package is not installed", {
  skip_if(
    requireNamespace("topofmri", quietly = TRUE),
    "topofmri is installed, cannot test missing-package path"
  )
  # "fugw" maps to "topofmri" which should not be installed
  result <- neuralign:::.try_autoload_aligner("fugw")
  expect_false(result)
})

test_that(".validate_aligner_requirements errors for unknown aligner", {
  set.seed(42)
  data_list <- list(
    "sub-01" = matrix(rnorm(50), 10, 5),
    "sub-02" = matrix(rnorm(50), 10, 5)
  )
  adat <- AlignmentData(data_list)

  expect_error(
    neuralign:::.validate_aligner_requirements("totally_unknown_aligner_xyz", adat),
    "Unknown aligner"
  )
})

test_that(".validate_aligner_requirements errors when needs_geometry but data has no geometry", {
  neuralign:::.clear_registry()

  dummy_fit <- function(data, reference, train_idx = NULL, ...) {
    list(transforms = list(), reference_data = NULL, space_from = NULL, space_to = NULL)
  }

  register_aligner(
    "geom_test",
    dummy_fit,
    capabilities = list(needs_geometry = TRUE)
  )
  on.exit(unregister_aligner("geom_test"))

  data_list <- list(
    "sub-01" = matrix(rnorm(50), 10, 5),
    "sub-02" = matrix(rnorm(50), 10, 5)
  )
  adat <- AlignmentData(data_list)  # no geometry

  expect_error(
    neuralign:::.validate_aligner_requirements("geom_test", adat),
    "requires geometry"
  )
})

test_that(".validate_aligner_requirements errors when needs_design but data has no design", {
  neuralign:::.clear_registry()

  dummy_fit <- function(data, reference, train_idx = NULL, ...) {
    list(transforms = list(), reference_data = NULL, space_from = NULL, space_to = NULL)
  }

  register_aligner(
    "design_test",
    dummy_fit,
    capabilities = list(needs_design = TRUE)
  )
  on.exit(unregister_aligner("design_test"))

  data_list <- list(
    "sub-01" = matrix(rnorm(50), 10, 5),
    "sub-02" = matrix(rnorm(50), 10, 5)
  )
  adat <- AlignmentData(data_list)  # no design

  expect_error(
    neuralign:::.validate_aligner_requirements("design_test", adat),
    "requires design"
  )
})
