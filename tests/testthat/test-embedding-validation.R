test_that("embedding providers must return complete, conformable named matrices", {
  scenario <- "valid"
  method <- "embedding_validation_shapes"
  on.exit(unregister_aligner(method), add = TRUE)

  fit_embedding <- function(data, reference, train_idx = NULL, ...) {
    aligned <- setNames(
      lapply(data@subjects, function(subject) matrix(seq_len(6L), 2L, 3L)),
      data@subjects
    )
    aligned <- switch(
      scenario,
      unnamed = unname(aligned),
      missing = aligned[1L],
      non_matrix = { aligned[[1L]] <- seq_len(3L); aligned },
      observation_mismatch = {
        aligned[[1L]] <- matrix(seq_len(4L), 2L, 2L)
        aligned
      },
      dimension_mismatch = {
        aligned[[2L]] <- matrix(seq_len(9L), 3L, 3L)
        aligned
      },
      aligned
    )
    list(aligned = aligned)
  }
  register_aligner(
    method,
    fit_embedding,
    capabilities = list(returns = "embedding", supports_new_data = FALSE)
  )
  data <- AlignmentData(list(
    s1 = matrix(seq_len(6L), 2L, 3L),
    s2 = matrix(seq_len(6L) + 10L, 2L, 3L)
  ))

  scenario <- "unnamed"
  expect_error(
    fit_alignment(data, method = method, reference = "s1", compute_quality = FALSE),
    "expected 'aligned' to be a named list"
  )
  scenario <- "missing"
  expect_error(
    fit_alignment(data, method = method, reference = "s1", compute_quality = FALSE),
    "missing aligned embeddings.*s2"
  )
  scenario <- "non_matrix"
  expect_error(
    fit_alignment(data, method = method, reference = "s1", compute_quality = FALSE),
    "not matrix-like"
  )
  scenario <- "observation_mismatch"
  expect_error(
    fit_alignment(data, method = method, reference = "s1", compute_quality = FALSE),
    "observation mismatch"
  )
  scenario <- "dimension_mismatch"
  expect_error(
    fit_alignment(data, method = method, reference = "s1", compute_quality = FALSE),
    "embedding dimension mismatch"
  )
  scenario <- "valid"
  expect_s4_class(
    fit_alignment(data, method = method, reference = "s1", compute_quality = FALSE),
    "AlignmentResult"
  )
})


test_that("embedding providers must return numeric finite values", {
  scenario <- "character"
  method <- "embedding_validation_values"
  on.exit(unregister_aligner(method), add = TRUE)

  fit_embedding <- function(data, reference, train_idx = NULL, ...) {
    aligned <- setNames(
      lapply(data@subjects, function(subject) matrix(seq_len(6L), 2L, 3L)),
      data@subjects
    )
    if (identical(scenario, "character")) {
      aligned[[1L]] <- matrix(letters[seq_len(6L)], 2L, 3L)
    } else if (identical(scenario, "non_finite")) {
      aligned[[1L]][[1L]] <- Inf
    }
    list(aligned = aligned)
  }
  register_aligner(
    method,
    fit_embedding,
    capabilities = list(returns = "embedding", supports_new_data = FALSE)
  )
  data <- AlignmentData(list(
    s1 = matrix(seq_len(6L), 2L, 3L),
    s2 = matrix(seq_len(6L) + 10L, 2L, 3L)
  ))

  expect_error(
    fit_alignment(data, method = method, reference = "s1", compute_quality = FALSE),
    "not numeric"
  )
  scenario <- "non_finite"
  expect_error(
    fit_alignment(data, method = method, reference = "s1", compute_quality = FALSE),
    "non-finite values"
  )
})
