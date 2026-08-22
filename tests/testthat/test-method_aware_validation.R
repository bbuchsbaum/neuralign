test_that("fit_alignment checks shared features only when required", {
  set.seed(123)
  adat <- AlignmentData(list(
    s1 = matrix(rnorm(30), 10, 3),
    s2 = matrix(rnorm(36), 12, 3)
  ))

  if (!is_aligner_registered("procrustes")) {
    ensure_test_aligner("procrustes")
  }

  # procrustes requires shared feature dims
  expect_error(
    fit_alignment(adat, method = "procrustes", cv = "none"),
    "different numbers of features"
  )

  # A test aligner that declares variable feature dims support should not fail
  # early validation when nrow differs across subjects.
  tryCatch(unregister_aligner("dummy_varfeat"), error = function(e) NULL)

  dummy_fit <- function(data, reference = "first", train_idx = NULL, fit_context = NULL, provider_plan = NULL, ...) {
    if (is.null(train_idx)) train_idx <- seq_along(data@subjects)
    train_data <- data[train_idx]
    data_list <- get_data_list(train_data)
    transforms <- lapply(data_list, function(x) {
      diag(nrow = 2L, ncol = nrow(x))
    })
    names(transforms) <- names(data_list)
    list(
      transforms = transforms,
      reference_data = get_subject_data(
        train_data,
        train_data@subjects[[1L]]
      )[seq_len(2L), , drop = FALSE],
      space_from = train_data@space,
      space_to = train_data@space,
      method_state = list()
    )
  }

  register_aligner(
    name = "dummy_varfeat",
    fit_fn = dummy_fit,
    apply_fn = NULL,
    capabilities = list(
      supports_cv = FALSE,
      cv_axes = character(0),
      needs_geometry = FALSE,
      needs_design = FALSE,
      requires_shared_features = FALSE,
      requires_shared_observations = FALSE,
      returns_invertible = TRUE,
      transform_type = "linear",
      mass_preserving = FALSE,
      returns = "operator",
      supports_new_subject = FALSE,
      supports_new_data = FALSE,
      reference_types = c("subject")
    ),
    package = "neuralign",
    description = "test-only aligner for variable feature dims",
    version = "0.0.0"
  )

  on.exit({
    tryCatch(unregister_aligner("dummy_varfeat"), error = function(e) NULL)
  }, add = TRUE)

  res <- fit_alignment(adat,
    method = "dummy_varfeat",
    reference = "s1",
    cv = "none",
    compute_quality = FALSE
  )
  expect_s4_class(res, "AlignmentResult")
})
