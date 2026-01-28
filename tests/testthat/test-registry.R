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

test_that("list_aligners is an alias for available_aligners", {
  neuralign:::.clear_registry()

  dummy_fit <- function(data, reference, train_idx = NULL, ...) {
    list(transforms = list(), reference_data = NULL, space_from = NULL, space_to = NULL)
  }

  register_aligner("method_a", dummy_fit, package = "pkg_a")

  expect_identical(list_aligners(), available_aligners())
  expect_identical(list_aligners(details = TRUE), available_aligners(details = TRUE))
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

test_that("register_aligner enforces supports_new_data=FALSE for embedding-returning aligners", {
  neuralign:::.clear_registry()

  dummy_fit <- function(data, reference, train_idx = NULL, ...) {
    k <- 3
    n_obs <- ncol(get_subject_data(data, data@subjects[[1L]]))
    aligned <- lapply(data@subjects, function(s) matrix(0, k, n_obs))
    names(aligned) <- data@subjects
    list(aligned = aligned, reference_data = NULL, space_from = NULL, space_to = NULL)
  }

  expect_error(
    register_aligner(
      "emb_bad",
      dummy_fit,
      capabilities = list(returns = "embedding")
    ),
    "supports_new_data"
  )

  register_aligner(
    "emb_ok",
    dummy_fit,
    capabilities = list(returns = "embedding", supports_new_data = FALSE)
  )
  expect_true(is_aligner_registered("emb_ok"))
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
    "must be .operator. or .embedding"
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
    "must be .operator. or .embedding"
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


# ---------- validate_aligner_contract() tests ----------

test_that("validate_aligner_contract passes for valid aligner", {
  fit_fn <- function(data, reference, train_idx = NULL, ...) {
    list(transforms = list(), reference_data = NULL)
  }
  apply_fn <- function(fit_result, new_data, ...) {
    list(transforms = list())
  }

  expect_true(
    validate_aligner_contract("valid", fit_fn, apply_fn,
                              capabilities = list(returns = "operator"))
  )
})

test_that("validate_aligner_contract errors when fit_fn is not a function", {
  expect_error(
    validate_aligner_contract("bad", "not_a_function"),
    "fit_fn must be a function"
  )
})

test_that("validate_aligner_contract errors when fit_fn lacks required formals", {
  bad_fit <- function(x, y) {}
  expect_error(
    validate_aligner_contract("bad", bad_fit),
    "missing required formals"
  )
})

test_that("validate_aligner_contract accepts fit_fn with ... in lieu of formals", {
  fit_dots <- function(...) {
    list(transforms = list())
  }
  expect_true(validate_aligner_contract("dots", fit_dots))
})

test_that("validate_aligner_contract errors when apply_fn lacks required formals", {
  good_fit <- function(data, reference, ...) {}
  bad_apply <- function(x) {}

  expect_error(
    validate_aligner_contract("bad_apply", good_fit, bad_apply),
    "apply_fn missing required formals"
  )
})

test_that("validate_aligner_contract errors for invalid capabilities$returns", {
  good_fit <- function(data, reference, ...) {}

  expect_error(
    validate_aligner_contract("bad", good_fit,
                              capabilities = list(returns = "xyz")),
    "must be .operator. or .embedding"
  )
})

test_that("validate_aligner_contract errors when capabilities is not a list", {
  good_fit <- function(data, reference, ...) {}
  expect_error(
    validate_aligner_contract("bad", good_fit, capabilities = "not_list"),
    "capabilities must be a list"
  )
})

test_that("validate_aligner_contract errors for future api_version", {
  good_fit <- function(data, reference, ...) {}
  expect_error(
    validate_aligner_contract("future", good_fit, api_version = 999L),
    "supports up to"
  )
})

test_that("register_aligner stores api_version in registry entry", {
  neuralign:::.clear_registry()

  dummy_fit <- function(data, reference, train_idx = NULL, ...) {
    list(transforms = list(), reference_data = NULL)
  }

  register_aligner("api_test", dummy_fit, api_version = 1L)

  entry <- get_aligner("api_test")
  expect_equal(entry$api_version, 1L)
  unregister_aligner("api_test")
})

test_that("NEURALIGN_ALIGNER_API_VERSION is exported and integer", {
  expect_true(is.integer(NEURALIGN_ALIGNER_API_VERSION))
  expect_equal(NEURALIGN_ALIGNER_API_VERSION, 1L)
})


# ---------- validate_aligner_contract capability field validation ----------

test_that("validate_aligner_contract errors on non-character transform_type", {
  good_fit <- function(data, reference, ...) {}
  expect_error(
    validate_aligner_contract("bad_tt", good_fit,
      capabilities = list(transform_type = 42)),
    "transform_type must be a character"
  )
})

test_that("validate_aligner_contract errors on non-character cv_axes", {
  good_fit <- function(data, reference, ...) {}
  expect_error(
    validate_aligner_contract("bad_cva", good_fit,
      capabilities = list(cv_axes = 123)),
    "cv_axes must be a character"
  )
})

test_that("validate_aligner_contract errors on non-character reference_types", {
  good_fit <- function(data, reference, ...) {}
  expect_error(
    validate_aligner_contract("bad_ref", good_fit,
      capabilities = list(reference_types = TRUE)),
    "reference_types must be a character"
  )
})

test_that("validate_aligner_contract errors on invalid returns value", {
  good_fit <- function(data, reference, ...) {}
  expect_error(
    validate_aligner_contract("bad_ret", good_fit,
      capabilities = list(returns = "bogus")),
    "must be.*operator.*embedding"
  )
})

test_that("validate_aligner_contract passes with valid embedding returns", {
  good_fit <- function(data, reference, ...) {}
  expect_invisible(
    validate_aligner_contract("emb_ok", good_fit,
      capabilities = list(returns = "embedding"))
  )
})


# ---------- .validate_aligner_requirements guidance ----------

test_that(".validate_aligner_requirements errors when guidance is needed but absent", {
  neuralign:::.clear_registry()

  dummy_fit <- function(data, reference, ...) {
    list(transforms = list(), reference_data = NULL)
  }
  register_aligner("needs_guide", dummy_fit,
    capabilities = list(needs_guidance = TRUE))

  adat <- AlignmentData(list(s1 = matrix(1, 5, 3), s2 = matrix(1, 5, 3)))
  expect_error(
    neuralign:::.validate_aligner_requirements("needs_guide", adat),
    "requires guidance channels"
  )
  neuralign:::.clear_registry()
})

test_that(".validate_aligner_requirements errors when specific guidance types missing", {
  neuralign:::.clear_registry()

  dummy_fit <- function(data, reference, ...) {
    list(transforms = list(), reference_data = NULL)
  }
  register_aligner("needs_type", dummy_fit,
    capabilities = list(needs_guidance = TRUE,
                        guidance_types = c("intrinsic_geometry")))

  adat <- AlignmentData(list(s1 = matrix(1, 5, 3), s2 = matrix(1, 5, 3)))
  wrong_guidance <- list(
    s1 = list(list(type = "roi_anchor", value = diag(3))),
    s2 = list(list(type = "roi_anchor", value = diag(3)))
  )
  adat <- set_guidance(adat, wrong_guidance)

  expect_error(
    neuralign:::.validate_aligner_requirements("needs_type", adat),
    "requires guidance types.*intrinsic_geometry"
  )
  neuralign:::.clear_registry()
})

test_that(".validate_aligner_requirements passes with correct guidance", {
  neuralign:::.clear_registry()

  dummy_fit <- function(data, reference, ...) {
    list(transforms = list(), reference_data = NULL)
  }
  register_aligner("guide_ok", dummy_fit,
    capabilities = list(needs_guidance = TRUE,
                        guidance_types = c("intrinsic_geometry")))

  adat <- AlignmentData(list(s1 = matrix(1, 5, 3), s2 = matrix(1, 5, 3)))
  good_guidance <- list(
    s1 = list(list(type = "intrinsic_geometry", value = diag(3))),
    s2 = list(list(type = "intrinsic_geometry", value = diag(3)))
  )
  adat <- set_guidance(adat, good_guidance)

  expect_true(
    neuralign:::.validate_aligner_requirements("guide_ok", adat)
  )
  neuralign:::.clear_registry()
})

test_that(".validate_aligner_requirements errors on unknown aligner", {
  neuralign:::.clear_registry()
  expect_error(
    neuralign:::.validate_aligner_requirements("nonexistent", NULL),
    "Unknown aligner"
  )
})
