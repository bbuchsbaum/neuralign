test_that("save_alignment and load_alignment work", {
  ensure_test_aligner("procrustes")

  set.seed(888)
  data_list <- make_test_data_list(n_subjects = 2, n_features = 10, n_obs = 10)
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
  ensure_test_aligner("procrustes")

  set.seed(887)
  data_list <- make_test_data_list(n_subjects = 2, n_features = 10, n_obs = 10)
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
  ensure_test_aligner("procrustes")

  set.seed(886)
  data_list <- make_test_data_list(n_subjects = 2, n_features = 10, n_obs = 10)
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

# ---------------------------------------------------------------------------
# Additional coverage tests
# ---------------------------------------------------------------------------

test_that("save_alignment errors on invalid model argument", {

  tmp_file <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp_file))

  expect_error(
    save_alignment("not_a_model", tmp_file),
    "must be an AlignmentModel, AlignmentResult, AlignedStudy, or AlignedResampleSet"
  )

  expect_error(
    save_alignment(42, tmp_file),
    "must be an AlignmentModel, AlignmentResult, AlignedStudy, or AlignedResampleSet"
  )

  expect_error(
    save_alignment(data.frame(x = 1), tmp_file),
    "must be an AlignmentModel, AlignmentResult, AlignedStudy, or AlignedResampleSet"
  )
})

test_that("load_alignment warns for raw model saved directly", {
  transforms <- list("sub-01" = diag(5))
  model <- AlignmentModel(transforms, reference = NULL, method = "test")

  tmp_file <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp_file))

  # Save directly with saveRDS (not via save_alignment)
  saveRDS(model, tmp_file)

  expect_warning(
    loaded <- load_alignment(tmp_file),
    "not saved with save_alignment"
  )

  expect_s4_class(loaded, "AlignmentModel")
  expect_equal(loaded@method, "test")
})

test_that("load_alignment warns for raw AlignmentResult saved directly", {
  transforms <- list("sub-01" = diag(5))
  model <- AlignmentModel(transforms, reference = NULL, method = "test")
  result <- AlignmentResult(model, aligned = list("sub-01" = diag(5)))

  tmp_file <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp_file))

  # Save directly with saveRDS
  saveRDS(result, tmp_file)

  expect_warning(
    loaded <- load_alignment(tmp_file),
    "not saved with save_alignment"
  )

  expect_s4_class(loaded, "AlignmentResult")
})

test_that("load_alignment errors on non-alignment RDS file", {
  tmp_file <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp_file))

  # Save a plain data frame -- not an alignment model at all
  saveRDS(data.frame(x = 1:5, y = letters[1:5]), tmp_file)

  expect_error(
    load_alignment(tmp_file),
    "does not contain a supported alignment object"
  )
})

test_that("load_alignment fails closed on integrity check failure", {
  transforms <- list("sub-01" = diag(5))
  model <- AlignmentModel(transforms, reference = NULL, method = "test")

  tmp_file <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp_file))

  # Save properly
  save_alignment(model, tmp_file)

  # Load, tamper with the object, and re-save to corrupt the hash
  save_data <- readRDS(tmp_file)
  save_data$object@transforms[["sub-01"]] <- diag(3)  # change transform
  # Keep original hash so verification will fail
  saveRDS(save_data, tmp_file)

  expect_error(
    load_alignment(tmp_file, verify = TRUE),
    "integrity check failed"
  )
})

test_that("load_alignment shows version mismatch message", {
  transforms <- list("sub-01" = diag(5))
  model <- AlignmentModel(transforms, reference = NULL, method = "test")

  tmp_file <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp_file))

  save_alignment(model, tmp_file)

  # Tamper with the saved version so it differs from the current version

  save_data <- readRDS(tmp_file)
  save_data$package_version <- "0.0.0.9000"
  saveRDS(save_data, tmp_file)

  expect_message(
    loaded <- load_alignment(tmp_file, verify = TRUE),
    "Object was saved with neuralign 0\\.0\\.0\\.9000"
  )

  expect_s4_class(loaded, "AlignmentModel")
})

test_that("versioned serialization round-trips analysis-facing representations", {
  model <- AlignmentModel(
    transforms = list(s1 = diag(2)),
    reference = "s1",
    method = "procrustes",
    train_subjects = "s1"
  )
  space <- shared_feature_space_from_model(model, dimension = 2L)
  block <- AlignedBlock(
    values = matrix(1:6, nrow = 3L, ncol = 2L),
    subject_id = "s1",
    shared_space_id = space$id
  )
  study <- AlignedStudy(
    blocks = list(s1 = block),
    shared_space = space,
    model = model
  )
  resamples <- AlignedResampleSet(splits = list(
    fold_1 = list(
      model = model,
      shared_space = space,
      analysis = study,
      assessment = study,
      metadata = list(fold_id = "fold_1")
    )
  ))
  result <- AlignmentResult(model, aligned = list(s1 = diag(2)))

  objects <- list(
    model = model,
    result = result,
    study = study,
    resamples = resamples
  )
  for (nm in names(objects)) {
    path <- tempfile(fileext = ".rds")
    on.exit(unlink(path), add = TRUE)
    save_alignment(objects[[nm]], path, include_data = TRUE)
    envelope <- readRDS(path)

    expect_identical(envelope$format_id, "neuralign-alignment", info = nm)
    expect_identical(envelope$format_version, 2L, info = nm)
    expect_identical(envelope$integrity$algorithm, "sha256", info = nm)
    expect_s4_class(load_alignment(path), class(objects[[nm]])[[1L]])
  }
})

test_that("load_alignment reads the legacy version-1 envelope", {
  model <- AlignmentModel(list(s1 = diag(2)), reference = "s1", method = "test")
  legacy <- list(
    object = model,
    neuralign_version = as.character(utils::packageVersion("neuralign")),
    hash = digest::digest(model, algo = "md5")
  )
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path))
  saveRDS(legacy, path)

  expect_s4_class(load_alignment(path), "AlignmentModel")
})

test_that("load_alignment rejects a newer envelope version", {
  model <- AlignmentModel(list(s1 = diag(2)), reference = "s1", method = "test")
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path))
  save_alignment(model, path)
  envelope <- readRDS(path)
  envelope$format_version <- 3L
  saveRDS(envelope, path)

  expect_error(load_alignment(path), "newer than supported")
})

test_that("export_alignment handles AlignmentResult by extracting model", {
  transforms <- list(
    "sub-01" = diag(3),
    "sub-02" = matrix(1:9, 3, 3)
  )
  model <- AlignmentModel(transforms, reference = NULL, method = "test_result")
  result <- AlignmentResult(model, aligned = list(
    "sub-01" = matrix(rnorm(9), 3, 3),
    "sub-02" = matrix(rnorm(9), 3, 3)
  ))

  tmp_path <- tempfile()
  dir_path <- paste0(tmp_path, "_transforms")
  on.exit(unlink(dir_path, recursive = TRUE))

  # Pass an AlignmentResult rather than an AlignmentModel
  files <- export_alignment(result, tmp_path, format = "csv")

  expect_true(dir.exists(dir_path))
  expect_true(file.exists(file.path(dir_path, "sub-01.csv")))
  expect_true(file.exists(file.path(dir_path, "sub-02.csv")))
  expect_true(file.exists(file.path(dir_path, "_metadata.csv")))
})

test_that("export_alignment to JSON with matrix reference", {
  skip_if_not_installed("jsonlite")

  ref_matrix <- matrix(rnorm(9), 3, 3)
  transforms <- list(
    "sub-01" = diag(3),
    "sub-02" = matrix(1:9, 3, 3)
  )
  model <- AlignmentModel(
    transforms,
    reference = NULL,
    reference_data = ref_matrix,
    method = "test_json_ref"
  )

  tmp_file <- tempfile()
  json_file <- paste0(tmp_file, ".json")
  on.exit(unlink(json_file))

  files <- export_alignment(model, tmp_file, format = "json")

  expect_true(file.exists(json_file))

  json_data <- jsonlite::fromJSON(json_file)
  expect_equal(json_data$method, "test_json_ref")

  # reference should be a matrix (numeric array) in the JSON
  expect_true(is.matrix(json_data$reference) || is.numeric(json_data$reference))
})

test_that("export_alignment to JSON with character reference", {
  skip_if_not_installed("jsonlite")

  transforms <- list(
    "sub-01" = diag(3),
    "sub-02" = matrix(1:9, 3, 3)
  )
  # reference_data is NULL so reference falls through to as.character path
  model <- AlignmentModel(
    transforms,
    reference = "sub-01",
    reference_data = NULL,
    method = "test_json_char"
  )

  tmp_file <- tempfile()
  json_file <- paste0(tmp_file, ".json")
  on.exit(unlink(json_file))

  files <- export_alignment(model, tmp_file, format = "json")

  expect_true(file.exists(json_file))

  json_data <- jsonlite::fromJSON(json_file)
  expect_equal(json_data$method, "test_json_char")
  # reference should be a character
  expect_true(is.character(json_data$reference))
})

test_that("export_alignment to MATLAB format works or skips", {
  skip_if_not_installed("R.matlab")

  transforms <- list(
    "sub-01" = diag(3),
    "sub-02" = matrix(1:9, 3, 3)
  )
  ref_matrix <- matrix(rnorm(9), 3, 3)
  model <- AlignmentModel(
    transforms,
    reference = NULL,
    reference_data = ref_matrix,
    method = "test_mat"
  )

  tmp_path <- tempfile()
  mat_file <- paste0(tmp_path, ".mat")
  on.exit(unlink(mat_file))

  files <- export_alignment(model, tmp_path, format = "mat")

  expect_true(file.exists(mat_file))
  expect_true(mat_file %in% files)
})

test_that("export_alignment to MATLAB errors without R.matlab", {
  # This test only makes sense if R.matlab is NOT installed
  skip_if(requireNamespace("R.matlab", quietly = TRUE),
          "R.matlab is installed, cannot test missing-package error")

  transforms <- list("sub-01" = diag(3))
  model <- AlignmentModel(transforms, reference = NULL, method = "test_mat")

  tmp_path <- tempfile()
  expect_error(
    export_alignment(model, tmp_path, format = "mat"),
    "R.matlab"
  )
})

test_that("import_alignment from CSV round-trip", {
  transforms <- list(
    "sub-01" = diag(4),
    "sub-02" = matrix(seq_len(16), 4, 4)
  )
  original <- AlignmentModel(transforms, reference = NULL, method = "csv_rt")

  tmp_path <- tempfile()
  dir_path <- paste0(tmp_path, "_transforms")
  on.exit(unlink(dir_path, recursive = TRUE))

  export_alignment(original, tmp_path, format = "csv")

  imported <- import_alignment(tmp_path, format = "csv")

  expect_s4_class(imported, "AlignmentModel")
  expect_equal(names(imported@transforms), names(original@transforms))

  # Verify transforms are close (CSV round-trip may lose some precision;

  # dimnames may be added by read.csv so strip them for comparison)
  for (subj in names(original@transforms)) {
    expect_equal(
      unname(imported@transforms[[subj]]),
      unname(original@transforms[[subj]]),
      tolerance = 1e-10
    )
  }
})

test_that("import_alignment from JSON round-trip", {
  skip_if_not_installed("jsonlite")

  transforms <- list(
    "sub-01" = diag(3),
    "sub-02" = matrix(1:9, 3, 3)
  )
  original <- AlignmentModel(transforms, reference = NULL, method = "json_rt")

  tmp_file <- tempfile()
  json_file <- paste0(tmp_file, ".json")
  on.exit(unlink(json_file))

  export_alignment(original, tmp_file, format = "json")
  imported <- import_alignment(tmp_file, format = "json")

  expect_s4_class(imported, "AlignmentModel")
  expect_equal(imported@method, "json_rt")
  expect_equal(names(imported@transforms), names(original@transforms))

  for (subj in names(original@transforms)) {
    expect_equal(
      unname(imported@transforms[[subj]]),
      unname(original@transforms[[subj]]),
      tolerance = 1e-10
    )
  }
})

test_that("import_alignment from CSV errors on missing directory", {
  tmp_path <- tempfile()  # this directory will not exist
  expect_error(
    import_alignment(tmp_path, format = "csv"),
    "Directory not found"
  )
})

test_that("import_alignment from CSV errors on missing metadata file", {
  tmp_path <- tempfile()
  dir_path <- paste0(tmp_path, "_transforms")
  dir.create(dir_path, recursive = TRUE)
  on.exit(unlink(dir_path, recursive = TRUE))

  # Directory exists but _metadata.csv is missing
  expect_error(
    import_alignment(tmp_path, format = "csv"),
    "Metadata file not found"
  )
})


# ---------- import/export edge cases ----------

test_that("import_alignment CSV warns on missing transform file", {
  tmp_path <- tempfile()
  dir_path <- paste0(tmp_path, "_transforms")
  dir.create(dir_path, recursive = TRUE)
  on.exit(unlink(dir_path, recursive = TRUE))

  meta <- data.frame(
    subject = c("s1", "s_missing"),
    nrow = c(3, 3),
    ncol = c(3, 3),
    stringsAsFactors = FALSE
  )
  utils::write.csv(meta, file.path(dir_path, "_metadata.csv"), row.names = FALSE)
  utils::write.csv(diag(3), file.path(dir_path, "s1.csv"), row.names = FALSE)

  expect_warning(
    model <- import_alignment(tmp_path, format = "csv"),
    "Transform file not found for s_missing"
  )
  expect_true(inherits(model, "AlignmentModel"))
  expect_true("s1" %in% names(model@transforms))
  expect_false("s_missing" %in% names(model@transforms))
})

test_that("import_alignment JSON with .json extension in path works", {
  skip_if_not_installed("jsonlite")
  ensure_test_aligner("procrustes")

  model <- AlignmentModel(
    transforms = list(s1 = diag(3), s2 = diag(3)),
    reference = "s1",
    reference_data = matrix(rnorm(9), 3, 3),
    method = "procrustes",
    train_subjects = c("s1", "s2")
  )

  tmp <- tempfile(fileext = ".json")
  on.exit(unlink(tmp))
  export_alignment(model, sub("\\.json$", "", tmp), format = "json")

  imported <- import_alignment(tmp, format = "json")
  expect_true(inherits(imported, "AlignmentModel"))
  expect_equal(imported@method, "procrustes")
})

test_that("import_alignment JSON errors on missing file", {
  skip_if_not_installed("jsonlite")
  expect_error(
    import_alignment("/nonexistent/path", format = "json"),
    "File not found"
  )
})

test_that("load_alignment with raw model (not save_alignment format)", {
  ensure_test_aligner("procrustes")
  model <- AlignmentModel(
    transforms = list(s1 = diag(5)),
    reference = "s1",
    method = "procrustes",
    train_subjects = "s1"
  )

  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp))
  saveRDS(model, tmp)

  expect_warning(
    loaded <- load_alignment(tmp),
    "not saved with save_alignment"
  )
  expect_true(inherits(loaded, "AlignmentModel"))
})

test_that("load_alignment errors on non-alignment object", {
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp))
  saveRDS(list(foo = "bar"), tmp)

  expect_error(load_alignment(tmp), "does not contain a supported alignment object")
})

test_that("load_alignment errors on nonexistent file", {
  expect_error(load_alignment("/no/such/file.rds"), "File not found")
})

test_that("save_alignment rejects non-model objects", {
  expect_error(
    save_alignment("not_a_model", tempfile()),
    "must be an AlignmentModel"
  )
})

test_that("load_alignment detects hash mismatch", {
  ensure_test_aligner("procrustes")
  model <- AlignmentModel(
    transforms = list(s1 = diag(5)),
    reference = "s1",
    method = "procrustes",
    train_subjects = "s1"
  )

  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp))

  # Save normally
  save_alignment(model, tmp)

  # Corrupt the hash
  save_data <- readRDS(tmp)
  save_data$integrity$hash <- "bogus_hash"
  saveRDS(save_data, tmp)

  expect_error(
    load_alignment(tmp),
    "integrity check failed"
  )
})
