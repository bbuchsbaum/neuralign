test_that("save_alignment and load_alignment work", {
  neuralign:::.register_procrustes()

  set.seed(888)
  data_list <- list(
    "sub-01" = matrix(rnorm(100), 10, 10),
    "sub-02" = matrix(rnorm(100), 10, 10)
  )
  adat <- AlignmentData(data_list)

  result <- fit_alignment(adat, method = "procrustes")

  # Save to temp file
  tmp_file <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp_file))

  save_alignment(result, tmp_file)
  expect_true(file.exists(tmp_file))

  # Load back
  loaded <- load_alignment(tmp_file)

  expect_s4_class(loaded, "AlignmentModel")
  expect_equal(loaded@method, "procrustes")
  expect_equal(names(loaded@transforms), names(result@model@transforms))
})

test_that("save_alignment adds .rds extension if missing", {
  transforms <- list("sub-01" = diag(5))
  model <- AlignmentModel(transforms, reference = NULL, method = "test")

  tmp_path <- tempfile()  # No extension
  on.exit(unlink(paste0(tmp_path, ".rds")))

  result_path <- save_alignment(model, tmp_path)

  expect_true(grepl("\\.rds$", result_path))
  expect_true(file.exists(result_path))
})

test_that("load_alignment verifies integrity", {
  transforms <- list("sub-01" = diag(5))
  model <- AlignmentModel(transforms, reference = NULL, method = "test")

  tmp_file <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp_file))

  save_alignment(model, tmp_file)

  # Load with verification (default)
  loaded <- load_alignment(tmp_file, verify = TRUE)
  expect_s4_class(loaded, "AlignmentModel")
})

test_that("load_alignment errors on missing file", {
  expect_error(
    load_alignment("/nonexistent/path/model.rds"),
    "not found"
  )
})

test_that("save_alignment can include aligned data", {
  neuralign:::.register_procrustes()

  set.seed(887)
  data_list <- list(
    "sub-01" = matrix(rnorm(100), 10, 10),
    "sub-02" = matrix(rnorm(100), 10, 10)
  )
  adat <- AlignmentData(data_list)

  result <- fit_alignment(adat, method = "procrustes")

  tmp_file <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp_file))

  # Save with data
  save_alignment(result, tmp_file, include_data = TRUE)

  loaded <- load_alignment(tmp_file)

  # Should be AlignmentResult with aligned data
  expect_s4_class(loaded, "AlignmentResult")
  expect_true(length(loaded@aligned) > 0)
})

test_that("export_alignment to csv works", {
  transforms <- list(
    "sub-01" = diag(3),
    "sub-02" = matrix(1:9, 3, 3)
  )
  model <- AlignmentModel(transforms, reference = NULL, method = "test")

  tmp_path <- tempfile()
  dir_path <- paste0(tmp_path, "_transforms")
  on.exit(unlink(dir_path, recursive = TRUE))

  files <- export_alignment(model, tmp_path, format = "csv")

  expect_true(dir.exists(dir_path))
  expect_true(file.exists(file.path(dir_path, "sub-01.csv")))
  expect_true(file.exists(file.path(dir_path, "sub-02.csv")))
  expect_true(file.exists(file.path(dir_path, "_metadata.csv")))
})

test_that("export_alignment to json works", {
  skip_if_not_installed("jsonlite")

  transforms <- list(
    "sub-01" = diag(3),
    "sub-02" = matrix(1:9, 3, 3)
  )
  model <- AlignmentModel(transforms, reference = NULL, method = "test")

  tmp_file <- tempfile()
  json_file <- paste0(tmp_file, ".json")
  on.exit(unlink(json_file))

  files <- export_alignment(model, tmp_file, format = "json")

  expect_true(file.exists(json_file))

  # Read back and verify
  json_data <- jsonlite::fromJSON(json_file)
  expect_equal(json_data$method, "test")
  expect_equal(json_data$subjects, c("sub-01", "sub-02"))
})

test_that("import_alignment from json works", {
  skip_if_not_installed("jsonlite")

  # Export first
  transforms <- list(
    "sub-01" = diag(3),
    "sub-02" = matrix(1:9, 3, 3)
  )
  original <- AlignmentModel(transforms, reference = NULL, method = "original")

  tmp_file <- tempfile()
  json_file <- paste0(tmp_file, ".json")
  on.exit(unlink(json_file))

  export_alignment(original, tmp_file, format = "json")

  # Import
  imported <- import_alignment(tmp_file, format = "json")

  expect_s4_class(imported, "AlignmentModel")
  expect_equal(imported@method, "original")
  expect_equal(names(imported@transforms), names(original@transforms))
})

test_that("round-trip save/load preserves transforms exactly", {
  neuralign:::.register_procrustes()

  set.seed(886)
  data_list <- list(
    "sub-01" = matrix(rnorm(100), 10, 10),
    "sub-02" = matrix(rnorm(100), 10, 10)
  )
  adat <- AlignmentData(data_list)

  result <- fit_alignment(adat, method = "procrustes")
  original_model <- get_model(result)

  tmp_file <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp_file))

  save_alignment(original_model, tmp_file)
  loaded_model <- load_alignment(tmp_file)

  # Transforms should be identical
  for (subj in names(original_model@transforms)) {
    expect_equal(
      loaded_model@transforms[[subj]],
      original_model@transforms[[subj]],
      tolerance = 1e-15
    )
  }
})
