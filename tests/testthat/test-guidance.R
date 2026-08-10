test_that("set_guidance/get_guidance store per-subject guidance channels", {
  set.seed(1)
  adat <- AlignmentData(list(
    s1 = make_test_matrix(n_features = 10, n_obs = 3),
    s2 = make_test_matrix(n_features = 12, n_obs = 3)
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

  dummy_fit <- function(data, reference = "first", train_idx = NULL, fit_context = NULL, provider_plan = NULL, ...) {
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
  res <- fit_alignment(adat2, method = "dummy_guided", reference = "s1", compute_quality = FALSE)
  expect_s4_class(res, "AlignmentResult")
})


# ---------- Additional guidance tests ----------

test_that("guidance_channel validates inputs", {
  expect_error(guidance_channel(type = "", value = 1), "non-empty")
  expect_error(guidance_channel(type = 123, value = 1), "non-empty")
  expect_error(guidance_channel(type = "coords", value = 1, name = ""), "non-empty")
  expect_error(guidance_channel(type = "coords", value = 1, name = 42), "non-empty")

  ch <- guidance_channel("projector", matrix(1, 2, 3), name = "roi")
  expect_equal(ch$type, "projector")
  expect_equal(ch$name, "roi")
  expect_true(is.matrix(ch$value))
})

test_that("guidance_channel passes extra fields through", {
  ch <- guidance_channel("coords", value = matrix(0, 5, 3), name = "geo", k = 3, normalized = TRUE)
  expect_equal(ch$k, 3)
  expect_true(ch$normalized)
})

test_that("get_guidance returns empty channels for no-guidance data", {
  adat <- AlignmentData(list(
    s1 = matrix(1, 4, 3),
    s2 = matrix(1, 5, 3)
  ))
  g <- get_guidance(adat)
  expect_equal(names(g), c("s1", "s2"))
  expect_equal(length(g$s1), 0)
  expect_equal(length(g$s2), 0)
})

test_that("get_guidance filters by type", {
  adat <- AlignmentData(list(
    s1 = matrix(1, 4, 3),
    s2 = matrix(1, 5, 3)
  ))
  g <- list(
    s1 = list(
      roi = guidance_channel("projector", matrix(1, 2, 4)),
      geo = guidance_channel("coords", matrix(1, 4, 2))
    ),
    s2 = list(
      roi = guidance_channel("projector", matrix(1, 2, 5))
    )
  )
  adat2 <- set_guidance(adat, g)

  proj_only <- get_guidance(adat2, type = "projector")
  expect_equal(length(proj_only$s1), 1)
  expect_equal(proj_only$s1$roi$type, "projector")
  expect_equal(length(proj_only$s2), 1)

  coords_only <- get_guidance(adat2, type = "coords")
  expect_equal(length(coords_only$s1), 1)
  expect_equal(coords_only$s1$geo$type, "coords")
  expect_equal(length(coords_only$s2), 0)
})

test_that("get_guidance filters by subject and type", {
  adat <- AlignmentData(list(
    s1 = matrix(1, 4, 3),
    s2 = matrix(1, 5, 3)
  ))
  g <- list(
    s1 = list(
      roi = guidance_channel("projector", matrix(1, 2, 4)),
      geo = guidance_channel("coords", matrix(1, 4, 2))
    ),
    s2 = list(
      roi = guidance_channel("projector", matrix(1, 2, 5))
    )
  )
  adat2 <- set_guidance(adat, g)

  s1_proj <- get_guidance(adat2, subject = "s1", type = "projector")
  expect_equal(length(s1_proj), 1)
  expect_equal(s1_proj$roi$type, "projector")

  s1_coords <- get_guidance(adat2, subject = "s1", type = "coords")
  expect_equal(length(s1_coords), 1)
})

test_that("get_guidance errors on unknown subject", {
  adat <- AlignmentData(list(s1 = matrix(1, 4, 3)))
  expect_error(get_guidance(adat, subject = "s99"), "Unknown subject")
})

test_that("set_guidance errors on non-AlignmentData", {
  expect_error(set_guidance("not_adat", list()), "AlignmentData")
})

test_that("set_guidance errors when guidance is missing subjects", {
  adat <- AlignmentData(list(
    s1 = matrix(1, 4, 3),
    s2 = matrix(1, 5, 3)
  ))
  g <- list(s1 = list())
  expect_error(set_guidance(adat, g), "missing subjects")
})

test_that("set_guidance warns on extra subjects in guidance", {
  adat <- AlignmentData(list(s1 = matrix(1, 4, 3)))
  g <- list(
    s1 = list(geo = guidance_channel("coords", matrix(1, 4, 2))),
    s_extra = list(ch = guidance_channel("coords", matrix(1, 4, 2)))
  )
  expect_warning(set_guidance(adat, g), "unknown subjects")
})

test_that("set_guidance validates coords dimension mismatch", {
  adat <- AlignmentData(list(s1 = matrix(1, 4, 3)))
  bad <- list(
    s1 = list(geo = guidance_channel("coords", matrix(1, 5, 2)))  # nrow=5, but subject has 4 features
  )
  expect_error(set_guidance(adat, bad), "dimension mismatch")
})

test_that("set_guidance validates sparse Matrix guidance without densifying", {
  skip_if_not_installed("Matrix")

  n_feat <- 40000L
  adat <- AlignmentData(list(s1 = matrix(1, n_feat, 1L)))
  P <- Matrix::Diagonal(n_feat)

  g <- list(
    s1 = list(
      proj = guidance_channel("projector", P)
    )
  )

  adat2 <- set_guidance(adat, g, validate = TRUE)
  ch <- get_guidance(adat2, subject = "s1")$proj
  expect_true(inherits(ch$value, "Matrix"))
})

test_that("set_guidance with validate=FALSE skips dimension check", {
  adat <- AlignmentData(list(s1 = matrix(1, 4, 3)))
  bad <- list(
    s1 = list(geo = guidance_channel("coords", matrix(1, 99, 2)))
  )
  adat2 <- set_guidance(adat, bad, validate = FALSE)
  expect_equal(nrow(get_guidance(adat2, subject = "s1")$geo$value), 99)
})

test_that("guidance_channel errors on missing value", {
  expect_error(
    neuralign:::.as_guidance_channel(list(type = "coords")),
    "value"
  )
})

test_that("set_guidance errors on non-list guidance", {
  adat <- AlignmentData(list(s1 = matrix(1, 4, 3)))
  expect_error(set_guidance(adat, "not_a_list"), "list")
})

test_that("set_guidance errors on unnamed guidance list", {
  adat <- AlignmentData(list(s1 = matrix(1, 4, 3)))
  expect_error(set_guidance(adat, list(list())), "named list")
})

test_that("set_guidance handles per-subject non-list channels", {
  adat <- AlignmentData(list(s1 = matrix(1, 4, 3)))
  g <- list(s1 = "not_a_list")
  expect_error(set_guidance(adat, g), "must be a list")
})

test_that("guidance_channel auto-assigns name from .as_guidance_channel", {
  ch <- list(type = "coords", value = 42)
  result <- neuralign:::.as_guidance_channel(ch, name = "my_channel")
  expect_equal(result$name, "my_channel")
  expect_equal(result$type, "coords")
})

test_that("validate_guidance_dims skips non-projector/non-coords types", {
  adat <- AlignmentData(list(s1 = matrix(1, 4, 3)))
  g <- list(
    s1 = list(custom = guidance_channel("geometry", value = list(something = TRUE)))
  )
  # "geometry" type is not validated for dimensions; should not error
  adat2 <- set_guidance(adat, g)
  expect_equal(get_guidance(adat2, subject = "s1")$custom$type, "geometry")
})

test_that("validate_guidance_dims errors on non-matrix projector", {
  adat <- AlignmentData(list(s1 = matrix(1, 4, 3)))
  g <- list(
    s1 = list(roi = guidance_channel("projector", value = c(1, 2, 3)))
  )
  expect_error(set_guidance(adat, g), "matrix-like")
})
