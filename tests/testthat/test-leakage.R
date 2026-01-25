# Tests for leakage detection functions

test_that(".check_leakage detects overlap correctly", {
  train <- c("sub-01", "sub-02", "sub-03")
  apply_no_overlap <- c("sub-04", "sub-05")
  apply_with_overlap <- c("sub-02", "sub-04")

  # No overlap
  result_no <- neuralign:::.check_leakage(apply_no_overlap, train, action = "silent")
  expect_false(result_no$has_leakage)
  expect_equal(length(result_no$overlap), 0)

  # With overlap
  result_yes <- neuralign:::.check_leakage(apply_with_overlap, train, action = "silent")
  expect_true(result_yes$has_leakage)
  expect_equal(result_yes$overlap, "sub-02")
  expect_equal(result_yes$n_overlap, 1)
  expect_equal(result_yes$n_apply, 2)
})

test_that(".check_leakage warns correctly", {
  train <- c("sub-01", "sub-02")
  apply_overlap <- c("sub-01", "sub-03")

  expect_warning(
    neuralign:::.check_leakage(apply_overlap, train, context = "test", action = "warn"),
    "Potential data leakage"
  )
})

test_that(".check_leakage errors correctly", {
  train <- c("sub-01", "sub-02")
  apply_overlap <- c("sub-01", "sub-02")

  expect_error(
    neuralign:::.check_leakage(apply_overlap, train, context = "test", action = "error"),
    "Potential data leakage"
  )
})

test_that(".check_leakage handles many overlapping subjects", {
  train <- paste0("sub-", sprintf("%02d", 1:10))
  apply_overlap <- paste0("sub-", sprintf("%02d", 1:8))

  # Should abbreviate the message
  expect_warning(
    result <- neuralign:::.check_leakage(apply_overlap, train, action = "warn"),
    "\\+\\d+ more"
  )
  expect_equal(result$n_overlap, 8)
})

test_that("validate_cv_setup accepts valid setups", {
  # Create mock CV folds
  cv_folds <- list(
    assignments = c("sub-01" = 1, "sub-02" = 1, "sub-03" = 2, "sub-04" = 2),
    folds = list(
      fold_1 = list(train = c(3, 4), test = c(1, 2)),
      fold_2 = list(train = c(1, 2), test = c(3, 4))
    )
  )

  # Medoid reference (computed from training) - always OK

  expect_true(validate_cv_setup(cv_folds, reference = "medoid"))
  expect_true(validate_cv_setup(cv_folds, reference = "centroid"))
  expect_true(validate_cv_setup(cv_folds, reference = "consensus"))
})

test_that("validate_cv_setup catches fixed reference leakage", {
  # Create CV folds where sub-01 is in test set for fold 1
  cv_folds <- list(
    assignments = c("sub-01" = 1, "sub-02" = 2, "sub-03" = 2),
    folds = list(
      fold_1 = list(train = c(2, 3), test = c(1)),  # sub-01 is test
      fold_2 = list(train = c(1, 3), test = c(2))
    )
  )

  # Using sub-01 as fixed reference causes leakage in fold 1
  expect_error(
    validate_cv_setup(cv_folds, reference = "sub-01"),
    "in test set"
  )

  # Using sub-03 as reference is OK (always in training)
  expect_true(validate_cv_setup(cv_folds, reference = "sub-03"))
})

test_that("assess_leakage_risk identifies no-CV risk", {
  transforms <- list("sub-01" = diag(5), "sub-02" = diag(5))
  model <- AlignmentModel(
    transforms = transforms,
    reference = "consensus",
    method = "procrustes",
    train_subjects = c("sub-01", "sub-02")
  )

  # No CV used
  risk <- assess_leakage_risk(model)
  expect_equal(risk$overall_risk, "medium")
  expect_true(any(grepl("No cross-validation", risk$issues)))
  expect_true(any(grepl("cv='loso'", risk$recommendations)))
})

test_that("assess_leakage_risk identifies test/train overlap", {
  transforms <- list("sub-01" = diag(5), "sub-02" = diag(5), "sub-03" = diag(5))
  model <- AlignmentModel(
    transforms = transforms,
    reference = "consensus",
    method = "procrustes",
    train_subjects = c("sub-01", "sub-02", "sub-03")
  )

  # Test subjects overlap with training
  risk <- assess_leakage_risk(model, test_subjects = c("sub-01", "sub-04"))
  expect_equal(risk$overall_risk, "high")
  expect_true(any(grepl("test subjects were in training", risk$issues)))
})

test_that("assess_leakage_risk handles AlignmentResult input", {
  transforms <- list("sub-01" = diag(5), "sub-02" = diag(5))
  model <- AlignmentModel(
    transforms = transforms,
    reference = "consensus",
    method = "procrustes",
    train_subjects = c("sub-01", "sub-02")
  )
  aligned <- list("sub-01" = matrix(1, 5, 3), "sub-02" = matrix(1, 5, 3))
  cv_info <- list(method = "loso", n_folds = 2)

  result <- AlignmentResult(model = model, aligned = aligned, cv_info = cv_info)

  # With CV, risk should be low
  risk <- assess_leakage_risk(result)
  expect_equal(risk$overall_risk, "low")
})

test_that("assess_leakage_risk checks method CV support", {
  # Register a method that doesn't support CV
  register_aligner(
    name = "no_cv_method",
    fit_fn = function(data, reference, ...) {
      list(transforms = list(), reference_data = NULL)
    },
    capabilities = list(supports_cv = FALSE, returns = "operator"),
    package = "neuralign"
  )

  transforms <- list("sub-01" = diag(5))
  model <- AlignmentModel(
    transforms = transforms,
    reference = "consensus",
    method = "no_cv_method",
    train_subjects = c("sub-01")
  )
  cv_info <- list(method = "loso", n_folds = 1)
  result <- AlignmentResult(
    model = model,
    aligned = list("sub-01" = diag(5)),
    cv_info = cv_info
  )

  risk <- assess_leakage_risk(result)
  expect_true(risk$overall_risk %in% c("medium", "high"))
  expect_true(any(grepl("may not fully support CV", risk$issues)))

  # Cleanup
  unregister_aligner("no_cv_method")
})

test_that("print_leakage_assessment formats output correctly", {
  # Low risk
  risk_low <- list(
    overall_risk = "low",
    issues = character(0),
    recommendations = character(0)
  )
  output_low <- capture.output(print_leakage_assessment(risk_low))
  expect_true(any(grepl("LOW", output_low)))
  expect_true(any(grepl("No leakage issues", output_low)))

  # High risk
  risk_high <- list(
    overall_risk = "high",
    issues = c("Issue 1", "Issue 2"),
    recommendations = c("Fix it")
  )
  output_high <- capture.output(print_leakage_assessment(risk_high))
  expect_true(any(grepl("HIGH", output_high)))
  expect_true(any(grepl("Issues Found", output_high)))
  expect_true(any(grepl("Issue 1", output_high)))
  expect_true(any(grepl("Recommendations", output_high)))
})
