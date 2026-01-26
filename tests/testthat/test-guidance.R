test_that("set_guidance/get_guidance store per-subject guidance channels", {
  set.seed(1)
  adat <- AlignmentData(list(
    s1 = matrix(rnorm(30), 10, 3),
    s2 = matrix(rnorm(36), 12, 3)
  ))

  g <- list(
    s1 = list(
      roi = guidance_channel("projector", matrix(1, 2, 10), name = "roi")
    ),
    s2 = list(
      roi = guidance_channel("projector", matrix(1, 2, 12), name = "roi")
    )
  )

  adat2 <- set_guidance(adat, g)
  g_all <- get_guidance(adat2)
  expect_equal(names(g_all), c("s1", "s2"))
  expect_equal(get_guidance(adat2, subject = "s1")$roi$type, "projector")
})

test_that("set_guidance validates dimensional compatibility", {
  adat <- AlignmentData(list(
    s1 = matrix(1, 10, 3),
    s2 = matrix(1, 12, 3)
  ))

  bad <- list(
    s1 = list(roi = guidance_channel("projector", matrix(1, 2, 9))),
    s2 = list(roi = guidance_channel("projector", matrix(1, 2, 12)))
  )
  expect_error(set_guidance(adat, bad), "dimension mismatch")
})

test_that("aligners can require guidance channels via capabilities", {
  adat <- AlignmentData(list(
    s1 = matrix(1, 10, 3),
    s2 = matrix(1, 12, 3)
  ))

  g <- list(
    s1 = list(roi = guidance_channel("projector", matrix(1, 2, 10))),
    s2 = list(roi = guidance_channel("projector", matrix(1, 2, 12)))
  )

  tryCatch(unregister_aligner("dummy_guided"), error = function(e) NULL)

  dummy_fit <- function(data, reference = "first", train_idx = NULL, ...) {
    if (is.null(train_idx)) train_idx <- seq_along(data@subjects)
    train_data <- data[train_idx]
    data_list <- get_data_list(train_data)
    transforms <- lapply(data_list, function(x) diag(nrow(x)))
    names(transforms) <- names(data_list)
    list(
      transforms = transforms,
      reference_data = get_subject_data(train_data, train_data@subjects[[1L]]),
      space_from = train_data@space,
      space_to = train_data@space,
      method_state = list()
    )
  }

  register_aligner(
    name = "dummy_guided",
    fit_fn = dummy_fit,
    apply_fn = NULL,
    capabilities = list(
      supports_cv = FALSE,
      cv_axes = character(0),
      needs_geometry = FALSE,
      needs_design = FALSE,
      needs_guidance = TRUE,
      guidance_types = c("projector"),
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
    description = "test-only guided aligner",
    version = "0.0.0"
  )

  on.exit({
    tryCatch(unregister_aligner("dummy_guided"), error = function(e) NULL)
  }, add = TRUE)

  expect_error(
    fit_alignment(adat, method = "dummy_guided", compute_quality = FALSE),
    "requires guidance"
  )

  adat2 <- set_guidance(adat, g)
  res <- fit_alignment(adat2, method = "dummy_guided", compute_quality = FALSE)
  expect_s4_class(res, "AlignmentResult")
})

