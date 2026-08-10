identity_provider_fit <- function(data, reference, train_idx = NULL, fit_context = NULL, provider_plan = NULL, ...) {
  subjects <- data@subjects[train_idx %||% seq_along(data@subjects)]
  transforms <- lapply(subjects, function(subject) {
    diag(nrow(get_subject_data(data, subject)))
  })
  names(transforms) <- subjects
  list(
    transforms = transforms,
    reference_data = get_subject_data(data, subjects[[1L]]),
    space_from = data@space,
    space_to = data@space,
    method_state = list()
  )
}

test_that("observation-CV overlap fails before any provider callback", {
  with_temp_registry(code = {
    calls <- new.env(parent = emptyenv())
    calls$prepare <- 0L
    calls$fit <- 0L

    prepare <- function(data, reference, resampling_plan, ...) {
      calls$prepare <- calls$prepare + 1L
      list(token = "prepared")
    }
    fit <- function(data, reference, train_idx = NULL, fit_context,
                    provider_plan, ...) {
      calls$fit <- calls$fit + 1L
      identity_provider_fit(data, reference, train_idx)
    }

    register_aligner(
      "provider_overlap_guard",
      fit,
      prepare_fn = prepare,
      api_version = 2L,
      capabilities = list(
        supports_cv = TRUE,
        cv_axes = c("subject", "observation"),
        reference_types = "subject"
      )
    )

    adat <- make_test_alignment_data(
      n_subjects = 2,
      n_features = 4,
      n_obs = 8
    )
    spec <- list(
      axis = "observation",
      method = "custom",
      folds = list(
        fold1 = list(train_idx = 1:4, test_idx = 4:6),
        fold2 = list(train_idx = 5:8, test_idx = 1:3)
      ),
      n_folds = 2L,
      fold_ids = c("fold1", "fold2")
    )

    expect_error(
      fit_alignment(
        adat,
        method = "provider_overlap_guard",
        reference = "sub-01",
        cv_folds = spec,
        compute_quality = FALSE,
        return_aligned = FALSE
      ),
      "Observation-level leakage"
    )
    expect_identical(calls$prepare, 0L)
    expect_identical(calls$fit, 0L)
  })
})

test_that("API v2 preflight and fit callbacks receive exact immutable contexts", {
  with_temp_registry(code = {
    calls <- new.env(parent = emptyenv())
    calls$events <- character(0)
    calls$plan <- NULL
    calls$contexts <- list()

    prepare <- function(data, reference, resampling_plan, ...) {
      calls$events <- c(calls$events, "prepare")
      calls$plan <- resampling_plan
      list(token = "prepared-once")
    }
    fit <- function(data, reference, train_idx = NULL, fit_context,
                    provider_plan, ...) {
      calls$events <- c(calls$events, paste0("fit:", fit_context$role))
      calls$contexts[[length(calls$contexts) + 1L]] <- fit_context
      expect_identical(provider_plan$token, "prepared-once")

      # A provider-local mutation must not alter the engine's plan or a later
      # callback context.
      if (identical(fit_context$fold_id, "fold1")) {
        fit_context$train_observation_indices$`sub-01`[[1L]] <- 999L
      }
      identity_provider_fit(data, reference, train_idx)
    }

    register_aligner(
      "provider_context_v2",
      fit,
      prepare_fn = prepare,
      api_version = 2L,
      capabilities = list(
        supports_cv = TRUE,
        cv_axes = c("subject", "observation"),
        reference_types = "subject"
      )
    )

    adat <- make_test_alignment_data(
      n_subjects = 2,
      n_features = 4,
      n_obs = 8
    )
    spec <- list(
      axis = "observation",
      method = "custom",
      folds = list(
        fold1 = list(train_idx = 1:3, test_idx = 5:6),
        fold2 = list(train_idx = 5:6, test_idx = 1:3)
      ),
      n_folds = 2L,
      fold_ids = c("fold1", "fold2")
    )

    result <- fit_alignment(
      adat,
      method = "provider_context_v2",
      reference = "sub-01",
      cv_folds = spec,
      compute_quality = FALSE,
      return_aligned = FALSE
    )

    expect_s4_class(result, "AlignmentResult")
    expect_identical(
      calls$events,
      c("prepare", "fit:observation_cv_fold", "fit:observation_cv_fold", "fit:full_fit")
    )
    expect_identical(calls$plan$contract, "neuralign_resampling_plan_v2")
    expect_identical(calls$plan$axis, "observation")
    expect_named(calls$plan$folds, c("fold1", "fold2"))

    fold1 <- calls$contexts[[1L]]
    expect_identical(fold1$contract, "neuralign_fit_context_v2")
    expect_identical(fold1$fold_id, "fold1")
    expect_identical(fold1$train_subject_ids, adat@subjects)
    expect_identical(fold1$test_subject_ids, character(0))
    expect_identical(fold1$train_observation_indices$`sub-01`, 1:3)
    expect_identical(fold1$test_observation_indices$`sub-01`, 5:6)
    expect_identical(fold1$guard_observation_indices$`sub-01`, c(4L, 7L, 8L))

    fold2 <- calls$contexts[[2L]]
    expect_identical(fold2$train_observation_indices$`sub-01`, 5:6)
    expect_identical(calls$plan$folds$fold1$train_observation_indices$`sub-01`, 1:3)

    full <- calls$contexts[[3L]]
    expect_identical(full$role, "full_fit")
    expect_identical(full$train_observation_indices$`sub-01`, 1:8)
    expect_identical(full$test_observation_indices$`sub-01`, integer(0))
    expect_identical(full$guard_observation_indices$`sub-01`, integer(0))
  })
})

test_that("API v2 subject CV exposes the generated folds used for execution", {
  with_temp_registry(code = {
    calls <- new.env(parent = emptyenv())
    calls$plan <- NULL
    calls$contexts <- list()

    prepare <- function(data, reference, resampling_plan, ...) {
      calls$plan <- resampling_plan
      list(token = "subject-plan")
    }
    fit <- function(data, reference, train_idx = NULL, fit_context,
                    provider_plan, ...) {
      calls$contexts[[fit_context$fold_id]] <- fit_context
      expect_identical(provider_plan$token, "subject-plan")
      identity_provider_fit(data, reference, train_idx)
    }
    apply <- function(fit_result, new_data, ...) {
      subject <- new_data@subjects[[1L]]
      transform <- diag(nrow(get_subject_data(new_data, subject)))
      list(transforms = setNames(list(transform), subject))
    }

    register_aligner(
      "provider_subject_context_v2",
      fit,
      apply_fn = apply,
      prepare_fn = prepare,
      api_version = 2L,
      capabilities = list(
        supports_cv = TRUE,
        cv_axes = c("subject", "observation"),
        reference_types = "consensus"
      )
    )
    adat <- make_test_alignment_data(
      n_subjects = 3,
      n_features = 4,
      n_obs = 6
    )

    result <- fit_alignment(
      adat,
      method = "provider_subject_context_v2",
      reference = "consensus",
      cv = "loso",
      compute_quality = FALSE,
      return_aligned = FALSE
    )

    expect_s4_class(result, "AlignmentResult")
    expect_identical(calls$plan$axis, "subject")
    expect_length(calls$plan$folds, 3L)
    expect_identical(names(calls$contexts), names(calls$plan$folds))
    for (fold_id in names(calls$plan$folds)) {
      context <- calls$contexts[[fold_id]]
      expect_identical(context$role, "subject_cv_fold")
      expect_length(context$train_subject_ids, 2L)
      expect_length(context$test_subject_ids, 1L)
      expect_length(context$guard_subject_ids, 0L)
      expect_setequal(
        c(context$train_subject_ids, context$test_subject_ids),
        adat@subjects
      )
    }
  })
})

test_that("API v2 fallback new-subject fits receive assessment contexts", {
  with_temp_registry(code = {
    roles <- character(0)
    fit <- function(data, reference, train_idx = NULL, fit_context,
                    provider_plan, ...) {
      roles <<- c(roles, fit_context$role)
      if (identical(fit_context$role, "subject_cv_assessment_fit")) {
        expect_length(fit_context$train_subject_ids, 1L)
        expect_length(fit_context$guard_subject_ids, 2L)
      }
      identity_provider_fit(data, reference, train_idx)
    }

    register_aligner(
      "provider_subject_fallback_v2",
      fit,
      api_version = 2L,
      capabilities = list(
        supports_cv = TRUE,
        cv_axes = "subject",
        reference_types = "consensus"
      )
    )
    adat <- make_test_alignment_data(
      n_subjects = 3,
      n_features = 4,
      n_obs = 6
    )

    result <- fit_alignment(
      adat,
      method = "provider_subject_fallback_v2",
      reference = "consensus",
      cv = "loso",
      compute_quality = FALSE,
      return_aligned = FALSE
    )

    expect_s4_class(result, "AlignmentResult")
    expect_identical(sum(roles == "subject_cv_fold"), 3L)
    expect_identical(sum(roles == "subject_cv_assessment_fit"), 3L)
  })
})

test_that("API v1 provider registrations fail closed", {
  with_temp_registry(code = {
    fit <- function(data, reference, train_idx = NULL, fit_context = NULL, provider_plan = NULL, ...) {
      identity_provider_fit(data, reference, train_idx)
    }

    expect_error(
      register_aligner(
        "provider_v1_compat",
        fit,
        api_version = 1L,
        capabilities = list(reference_types = "subject")
      ),
      "supports only API 2"
    )
    expect_false(is_aligner_registered("provider_v1_compat"))
  })
})

test_that("positional registration uses the current API", {
  with_temp_registry(code = {
    expect_true(register_aligner(
      "provider_positional",
      identity_provider_fit,
      NULL,
      list(reference_types = "subject")
    ))
    entry <- get_aligner("provider_positional")
    expect_identical(entry$api_version, NEURALIGN_ALIGNER_API_VERSION)
    expect_null(entry$prepare_fn)
    expect_identical(entry$capabilities$reference_types, "subject")
  })
})

test_that("the provider API validates lifecycle signatures", {
  fit_v2 <- function(data, reference, train_idx = NULL, fit_context,
                     provider_plan, ...) {
    identity_provider_fit(data, reference, train_idx)
  }
  bad_prepare <- function(x) NULL

  expect_error(
    validate_aligner_contract(
      "bad_prepare_v2",
      fit_v2,
      prepare_fn = bad_prepare,
      api_version = 2L
    ),
    "prepare_fn missing required formals"
  )

  bad_fit <- function(data, reference) NULL
  expect_error(
    validate_aligner_contract("bad_fit_v2", bad_fit, api_version = 2L),
    "fit_fn missing required formals"
  )
})
