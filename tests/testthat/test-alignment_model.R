# Tests for AlignmentModel class and methods

test_that("AlignmentModel constructor validates inputs", {
  # Valid model
  transforms <- list(
    "sub-01" = diag(10),
    "sub-02" = diag(10)
  )

  model <- AlignmentModel(
    transforms = transforms,
    reference = "sub-01",
    method = "test"
  )

  expect_s4_class(model, "AlignmentModel")
  expect_equal(length(model), 2)

  # Invalid: transforms must be a list

  expect_error(
    AlignmentModel(transforms = diag(10), reference = "x", method = "test"),
    "must be a named list"
  )

  # Invalid: transforms must be named
  expect_error(
    AlignmentModel(transforms = list(diag(10)), reference = "x", method = "test"),
    "must be named"
  )
})

test_that("AlignmentModel accessors work correctly", {
  transforms <- list(
    "sub-01" = matrix(1:4, 2, 2),
    "sub-02" = matrix(5:8, 2, 2),
    "sub-03" = matrix(9:12, 2, 2)
  )
  ref_data <- matrix(rnorm(4), 2, 2)

  model <- AlignmentModel(
    transforms = transforms,
    reference = "sub-01",
    reference_data = ref_data,
    method = "procrustes",
    space_from = "native",
    space_to = "MNI",
    train_subjects = c("sub-01", "sub-02", "sub-03")
  )

  # get_transform
  expect_equal(get_transform(model, "sub-01"), matrix(1:4, 2, 2))
  expect_equal(get_transform(model, "sub-02"), matrix(5:8, 2, 2))
  expect_error(get_transform(model, "sub-99"), "not found")

  # get_transforms
  expect_equal(get_transforms(model), transforms)

  # model_subjects
  expect_equal(model_subjects(model), c("sub-01", "sub-02", "sub-03"))

  # has_transform
  expect_true(has_transform(model, "sub-01"))
  expect_false(has_transform(model, "sub-99"))

  # get_reference
  expect_equal(get_reference(model), ref_data)

  # length
  expect_equal(length(model), 3)
})

test_that("AlignmentModel subsetting works", {
  transforms <- list(
    "sub-01" = diag(5),
    "sub-02" = diag(5) * 2,
    "sub-03" = diag(5) * 3
  )

  model <- AlignmentModel(
    transforms = transforms,
    reference = "consensus",
    method = "test"
  )

  # Subset by character
  sub_model <- model[c("sub-01", "sub-03")]
  expect_equal(length(sub_model), 2)
  expect_equal(model_subjects(sub_model), c("sub-01", "sub-03"))

  # Subset by index
  sub_model2 <- model[1:2]
  expect_equal(length(sub_model2), 2)
  expect_equal(model_subjects(sub_model2), c("sub-01", "sub-02"))

  # Invalid subset
  expect_error(model["sub-99"], "Unknown subjects")
})

test_that("AlignmentModel subsetting validates index bounds and types", {
  transforms <- list(
    "sub-01" = diag(3),
    "sub-02" = diag(3),
    "sub-03" = diag(3)
  )
  model <- AlignmentModel(transforms = transforms, reference = "consensus", method = "test")

  expect_error(model[4], "out of bounds")
  expect_error(model[-4], "out of bounds")
  expect_error(model[c(-1, 2)], "cannot mix negative and positive")
  expect_error(model[1.5], "integer-valued")

  empty <- model[0]
  expect_s4_class(empty, "AlignmentModel")
  expect_equal(length(empty), 0)
})

test_that("add_transform creates new model", {
  transforms <- list("sub-01" = diag(3))

  model <- AlignmentModel(
    transforms = transforms,
    reference = "consensus",
    method = "test"
  )

  new_transform <- matrix(1:9, 3, 3)
  new_model <- add_transform(model, "sub-02", new_transform)

  # Original unchanged
  expect_equal(length(model), 1)
  expect_false(has_transform(model, "sub-02"))

  # New model has both
  expect_equal(length(new_model), 2)
  expect_true(has_transform(new_model, "sub-02"))
  expect_equal(get_transform(new_model, "sub-02"), new_transform)
})

test_that("AlignmentModel show method works", {
  transforms <- list(
    "sub-01" = diag(10),
    "sub-02" = diag(10)
  )

  model <- AlignmentModel(
    transforms = transforms,
    reference = "consensus",
    method = "procrustes",
    space_from = "native",
    space_to = "MNI",
    train_subjects = c("sub-01", "sub-02")
  )

  # Capture output
  output <- capture.output(show(model))

  expect_true(any(grepl("AlignmentModel", output)))
  expect_true(any(grepl("procrustes", output)))
  expect_true(any(grepl("Subjects: 2", output)))
  expect_true(any(grepl("10 x 10", output)))
  expect_true(any(grepl("consensus", output)))
})

test_that(".format_space handles different space types", {
  # Character
  expect_equal(neuralign:::.format_space("MNI152"), "MNI152")

  # Other class
  x <- structure(list(), class = "my_space")
  expect_equal(neuralign:::.format_space(x), "my_space")
})

test_that("AlignmentModel provenance is recorded", {
  transforms <- list("sub-01" = diag(5))

  model <- AlignmentModel(
    transforms = transforms,
    reference = "consensus",
    method = "procrustes",
    params = list(scale = TRUE)
  )

  expect_true(!is.null(model@provenance$fitted_at))
  expect_true(!is.null(model@provenance$neuralign_version))
  expect_equal(model@provenance$params$scale, TRUE)
})

# ---------- NEW TESTS ----------

test_that("show() displays fold_specific reference correctly", {
  transforms <- list(
    "sub-01" = diag(5),
    "sub-02" = diag(5)
  )

  model <- AlignmentModel(
    transforms = transforms,
    reference = "fold_specific",
    method = "test",
    train_subjects = c("sub-01", "sub-02")
  )

  output <- capture.output(show(model))

  expect_true(any(grepl("fold-specific", output)))
  expect_true(any(grepl("no common anchor", output)))
})

test_that("show() displays subject ID reference correctly", {
  transforms <- list(
    "sub-01" = diag(5),
    "sub-02" = diag(5)
  )

  model <- AlignmentModel(
    transforms = transforms,
    reference = "sub-01",
    method = "test",
    train_subjects = c("sub-01", "sub-02")
  )

  output <- capture.output(show(model))

  expect_true(any(grepl("Reference: subject 'sub-01'", output, fixed = TRUE)))
})

test_that("show() displays template matrix reference correctly", {
  transforms <- list(
    "sub-01" = diag(5),
    "sub-02" = diag(5)
  )

  ref_matrix <- matrix(rnorm(25), 5, 5)
  model <- AlignmentModel(
    transforms = transforms,
    reference = ref_matrix,
    method = "test",
    train_subjects = c("sub-01", "sub-02")
  )

  output <- capture.output(show(model))

  expect_true(any(grepl("Reference: template matrix", output, fixed = TRUE)))
})

test_that(".format_space with gds_space-like object", {
  space <- structure(list(name = "MNI"), class = "gds_space")
  result <- neuralign:::.format_space(space)
  expect_true(grepl("gds_space", result))
  expect_true(grepl("MNI", result))
})

test_that(".format_space with gds_space without name", {
  space <- structure(list(), class = "gds_space")
  result <- neuralign:::.format_space(space)
  expect_true(grepl("gds_space", result))
  expect_true(grepl("unnamed", result))
})

test_that(".format_space with non-character, non-gds_space object", {
  obj <- structure(list(), class = "custom_space_type")
  result <- neuralign:::.format_space(obj)
  expect_equal(result, "custom_space_type")
})

test_that(".format_space with numeric object falls back to class name", {
  obj <- structure(42, class = "special_space")
  result <- neuralign:::.format_space(obj)
  expect_equal(result, "special_space")
})


# ---------- spaces_compatible() tests ----------

test_that("spaces_compatible returns TRUE when both are NULL", {
  expect_true(spaces_compatible(NULL, NULL))
})

test_that("spaces_compatible returns TRUE when one is NULL", {
  expect_true(spaces_compatible(NULL, "MNI"))
  expect_true(spaces_compatible("MNI", NULL))
})

test_that("spaces_compatible returns TRUE for identical strings", {
  expect_true(spaces_compatible("MNI152", "MNI152"))
})

test_that("spaces_compatible returns FALSE for different strings", {
  expect_false(spaces_compatible("MNI152", "native"))
})

test_that("spaces_compatible compares gds_space by name", {
  s1 <- structure(list(name = "MNI152"), class = "gds_space")
  s2 <- structure(list(name = "MNI152"), class = "gds_space")
  s3 <- structure(list(name = "native"), class = "gds_space")

  expect_true(spaces_compatible(s1, s2))
  expect_false(spaces_compatible(s1, s3))
})

test_that("spaces_compatible handles gds_space with missing name", {
  s1 <- structure(list(), class = "gds_space")
  s2 <- structure(list(), class = "gds_space")
  expect_true(spaces_compatible(s1, s2))
})

test_that("spaces_compatible uses all.equal fallback", {
  expect_true(spaces_compatible(1:5, 1:5))
  expect_false(spaces_compatible(1:5, 1:4))
})


# ---------- AlignmentModel provenance parameter ----------

test_that("AlignmentModel constructor accepts pre-built provenance", {
  custom_prov <- list(
    composed_from = list(method1 = "m1", method2 = "m2"),
    composed_at = Sys.time()
  )
  model <- AlignmentModel(
    transforms = list("sub-01" = diag(3)),
    reference = "consensus",
    method = "composed",
    provenance = custom_prov
  )
  expect_equal(model@provenance$composed_from$method1, "m1")
  expect_true(!is.null(model@provenance$composed_at))
  # Should NOT have auto-generated fields

  expect_null(model@provenance$fitted_at)

  out <- capture.output(show(model))
  expect_false(any(grepl("Fitted at:", out, fixed = TRUE)))
})

test_that("AlignmentModel with provenance=NULL auto-generates provenance", {
  model <- AlignmentModel(
    transforms = list("sub-01" = diag(3)),
    reference = "consensus",
    method = "test",
    params = list(scale = TRUE)
  )
  expect_true(!is.null(model@provenance$fitted_at))
  expect_true(!is.null(model@provenance$neuralign_version))
  expect_equal(model@provenance$params$scale, TRUE)
})


# ---------- Template reference with sparse Matrix ----------

test_that("show() works with sparse Matrix reference", {
  skip_if_not_installed("Matrix")

  sparse_ref <- Matrix::sparseMatrix(i = 1:5, j = 1:5, x = 1, dims = c(5, 5))
  model <- AlignmentModel(
    transforms = list("sub-01" = diag(5)),
    reference = sparse_ref,
    method = "test"
  )

  output <- capture.output(show(model))
  expect_true(any(grepl("template matrix", output)))
})
