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
