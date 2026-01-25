# Tests for transform validation helpers

test_that(".is_matrixish identifies matrix types", {
  expect_true(neuralign:::.is_matrixish(matrix(1:4, 2, 2)))
  expect_true(neuralign:::.is_matrixish(Matrix::Matrix(1:4, 2, 2)))
  expect_false(neuralign:::.is_matrixish(1:4))
  expect_false(neuralign:::.is_matrixish(data.frame(a = 1:2)))
  expect_false(neuralign:::.is_matrixish(list()))
})

test_that(".validate_operator_transforms accepts valid transforms", {
  data_list <- list(
    "sub-01" = matrix(rnorm(20), 10, 2),
    "sub-02" = matrix(rnorm(20), 10, 2)
  )

  transforms <- list(
    "sub-01" = diag(10),  # 10x10 for 10-feature data
    "sub-02" = diag(10)
  )

  # Should not error
  expect_true(
    neuralign:::.validate_operator_transforms(transforms, data_list, "test")
  )
})

test_that(".validate_operator_transforms rejects invalid transforms list", {
  data_list <- list("sub-01" = matrix(1, 5, 2))

  # Not a list
  expect_error(
    neuralign:::.validate_operator_transforms(diag(5), data_list, "ctx"),
    "named list"
  )

  # Empty list
  expect_error(
    neuralign:::.validate_operator_transforms(list(), data_list, "ctx"),
    "named list"
  )

  # Unnamed list
  expect_error(
    neuralign:::.validate_operator_transforms(list(diag(5)), data_list, "ctx"),
    "named list"
  )
})

test_that(".validate_operator_transforms catches missing subjects", {
  data_list <- list(
    "sub-01" = matrix(1, 5, 2),
    "sub-02" = matrix(1, 5, 2)
  )

  transforms <- list(
    "sub-01" = diag(5)
    # Missing sub-02
  )

  expect_error(
    neuralign:::.validate_operator_transforms(transforms, data_list, "ctx"),
    "missing transforms.*sub-02"
  )
})

test_that(".validate_operator_transforms catches non-matrix transforms", {
  data_list <- list("sub-01" = matrix(1, 5, 2))

  transforms <- list(
    "sub-01" = 1:25  # Not a matrix
  )

  expect_error(
    neuralign:::.validate_operator_transforms(transforms, data_list, "ctx"),
    "not a matrix"
  )
})

test_that(".validate_operator_transforms catches dimension mismatches", {
  data_list <- list(
    "sub-01" = matrix(1, 10, 5)  # 10 features
  )

  transforms <- list(
    "sub-01" = diag(5)  # 5x5, but data has 10 features
  )

  expect_error(
    neuralign:::.validate_operator_transforms(transforms, data_list, "ctx"),
    "dimension mismatch"
  )
})

test_that(".validate_operator_transforms works with sparse matrices", {
  data_list <- list(
    "sub-01" = matrix(rnorm(20), 10, 2)
  )

  transforms <- list(
    "sub-01" = Matrix::Diagonal(10)  # Sparse identity
  )

  expect_true(
    neuralign:::.validate_operator_transforms(transforms, data_list, "test")
  )
})

test_that(".validate_operator_transforms coerces non-matrix data", {
  # Data that can be coerced to matrix
  data_list <- list(
    "sub-01" = 1:10  # Vector, becomes 10x1 matrix
  )

  transforms <- list(
    "sub-01" = diag(10)
  )

  expect_true(
    neuralign:::.validate_operator_transforms(transforms, data_list, "test")
  )
})

test_that(".validate_operator_transforms provides informative error messages", {
  data_list <- list(
    "sub-01" = matrix(1, 8, 3)  # 8 features, 3 obs
  )

  transforms <- list(
    "sub-01" = matrix(1, 10, 5)  # 10x5, wrong dims
  )

  # Check error message includes dimensions
  expect_error(
    neuralign:::.validate_operator_transforms(transforms, data_list, "my_context"),
    "my_context.*sub-01.*10 x 5.*8 x 3"
  )
})
