#' Internal: Copy a provider lifecycle value
#'
#' Provider callbacks receive detached copies so mutations cannot change the
#' engine-owned resampling plan or a later callback's context.
#' @keywords internal
.copy_provider_value <- function(x) {
  unserialize(serialize(x, connection = NULL, version = 3L))
}

.empty_obs_indices <- function(subjects) {
  setNames(rep(list(integer(0)), length(subjects)), subjects)
}

.full_obs_indices <- function(data) {
  setNames(lapply(data@subjects, function(subject) {
    seq_len(ncol(get_subject_data(data, subject)))
  }), data@subjects)
}

.obs_ids_from_indices <- function(data, indices) {
  labels <- .resolve_obs_labels_by_subject(data)
  out <- lapply(data@subjects, function(subject) {
    idx <- indices[[subject]] %||% integer(0)
    if (is.null(labels)) {
      as.character(idx)
    } else {
      as.character(labels[[subject]][idx])
    }
  })
  names(out) <- data@subjects
  out
}

.subject_obs_indices <- function(data, subject_indices) {
  out <- .empty_obs_indices(data@subjects)
  for (subject in data@subjects[subject_indices]) {
    out[[subject]] <- seq_len(ncol(get_subject_data(data, subject)))
  }
  out
}

.new_fit_context <- function(data,
                             role,
                             axis,
                             fold_id = NULL,
                             train_subject_indices,
                             test_subject_indices = integer(0),
                             guard_subject_indices = integer(0),
                             train_observation_indices,
                             test_observation_indices,
                             guard_observation_indices) {
  subjects <- data@subjects
  full_subject_indices <- seq_along(subjects)
  full_observation_indices <- .full_obs_indices(data)

  list(
    contract = "neuralign_fit_context_v2",
    role = role,
    axis = axis,
    fold_id = fold_id,
    full_subject_indices = as.integer(full_subject_indices),
    train_subject_indices = as.integer(train_subject_indices),
    test_subject_indices = as.integer(test_subject_indices),
    guard_subject_indices = as.integer(guard_subject_indices),
    full_subject_ids = as.character(subjects),
    train_subject_ids = as.character(subjects[train_subject_indices]),
    test_subject_ids = as.character(subjects[test_subject_indices]),
    guard_subject_ids = as.character(subjects[guard_subject_indices]),
    full_observation_indices = full_observation_indices,
    train_observation_indices = train_observation_indices,
    test_observation_indices = test_observation_indices,
    guard_observation_indices = guard_observation_indices,
    full_observation_ids = .obs_ids_from_indices(data, full_observation_indices),
    train_observation_ids = .obs_ids_from_indices(data, train_observation_indices),
    test_observation_ids = .obs_ids_from_indices(data, test_observation_indices),
    guard_observation_ids = .obs_ids_from_indices(data, guard_observation_indices)
  )
}

.validate_obs_indices <- function(value, n_obs, context, allow_empty = FALSE) {
  if (!is.numeric(value)) {
    stop(sprintf("%s must contain integer indices", context), call. = FALSE)
  }
  idx <- as.integer(value)
  if (anyNA(idx) || any(value != idx)) {
    stop(sprintf("%s must contain finite integer indices", context), call. = FALSE)
  }
  if (!allow_empty && length(idx) == 0L) {
    stop(sprintf("%s must not be empty", context), call. = FALSE)
  }
  if (anyDuplicated(idx)) {
    stop(sprintf("%s contains duplicate indices", context), call. = FALSE)
  }
  if (length(idx) > 0L && any(idx < 1L | idx > n_obs)) {
    stop(sprintf("%s contains out-of-range indices", context), call. = FALSE)
  }
  idx
}

.normalize_observation_folds <- function(data, cv_folds) {
  folds <- cv_folds$folds
  fold_names <- names(folds)
  if (length(folds) < 2L) {
    stop("Observation cv_folds$folds must contain >= 2 folds", call. = FALSE)
  }
  if (is.null(fold_names) || anyNA(fold_names) || any(!nzchar(fold_names)) ||
      anyDuplicated(fold_names)) {
    stop("Observation cv_folds$folds must be uniquely named", call. = FALSE)
  }

  subjects <- data@subjects
  first_fold <- folds[[1L]]
  per_subject <- is.list(first_fold) && !is.null(names(first_fold)) &&
    is.null(first_fold$train_idx)

  normalized <- lapply(fold_names, function(fold_name) {
    fold <- folds[[fold_name]]
    train_obs <- .empty_obs_indices(subjects)
    test_obs <- .empty_obs_indices(subjects)

    if (isTRUE(per_subject)) {
      if (!is.list(fold) || is.null(names(fold))) {
        stop(sprintf("Observation fold '%s' must be keyed by subject", fold_name), call. = FALSE)
      }
      missing_subjects <- setdiff(subjects, names(fold))
      if (length(missing_subjects) > 0L) {
        stop(sprintf(
          "Observation fold '%s' is missing subjects: %s",
          fold_name, paste(missing_subjects, collapse = ", ")
        ), call. = FALSE)
      }
      for (subject in subjects) {
        subject_fold <- fold[[subject]]
        if (is.null(subject_fold$train_idx) || is.null(subject_fold$test_idx)) {
          stop(sprintf(
            "Observation fold '%s' subject '%s' must contain train_idx and test_idx",
            fold_name, subject
          ), call. = FALSE)
        }
        n_obs <- ncol(get_subject_data(data, subject))
        train_obs[[subject]] <- .validate_obs_indices(
          subject_fold$train_idx, n_obs,
          sprintf("Observation fold '%s' subject '%s' train_idx", fold_name, subject)
        )
        test_obs[[subject]] <- .validate_obs_indices(
          subject_fold$test_idx, n_obs,
          sprintf("Observation fold '%s' subject '%s' test_idx", fold_name, subject),
          allow_empty = TRUE
        )
      }
    } else {
      if (is.null(fold$train_idx) || is.null(fold$test_idx)) {
        stop(sprintf(
          "Observation fold '%s' must contain train_idx and test_idx",
          fold_name
        ), call. = FALSE)
      }
      for (subject in subjects) {
        n_obs <- ncol(get_subject_data(data, subject))
        train_obs[[subject]] <- .validate_obs_indices(
          fold$train_idx, n_obs,
          sprintf("Observation fold '%s' train_idx", fold_name)
        )
        test_obs[[subject]] <- .validate_obs_indices(
          fold$test_idx, n_obs,
          sprintf("Observation fold '%s' test_idx", fold_name),
          allow_empty = TRUE
        )
      }
    }

    # This is the authoritative observation-axis leakage gate. It runs while
    # the whole plan is still engine-owned and before provider preflight.
    check_obs_leakage(train_obs, test_obs, action = "error")

    guard_obs <- lapply(subjects, function(subject) {
      setdiff(
        seq_len(ncol(get_subject_data(data, subject))),
        union(train_obs[[subject]], test_obs[[subject]])
      )
    })
    names(guard_obs) <- subjects

    .new_fit_context(
      data = data,
      role = "observation_cv_fold",
      axis = "observation",
      fold_id = fold_name,
      train_subject_indices = seq_along(subjects),
      train_observation_indices = train_obs,
      test_observation_indices = test_obs,
      guard_observation_indices = guard_obs
    )
  })
  names(normalized) <- fold_names

  list(folds = normalized, per_subject = per_subject)
}

.build_resampling_plan <- function(data, method, cv, cv_folds, train_idx = NULL) {
  subjects <- data@subjects
  all_subjects <- seq_along(subjects)
  all_obs <- .full_obs_indices(data)
  empty_obs <- .empty_obs_indices(subjects)

  plan <- list(
    contract = "neuralign_resampling_plan_v2",
    api_version = 2L,
    method = method,
    axis = "none",
    strategy = cv,
    full_subject_ids = as.character(subjects),
    full_subject_indices = as.integer(all_subjects),
    full_observation_indices = all_obs,
    full_observation_ids = .obs_ids_from_indices(data, all_obs),
    folds = list(),
    contexts = list()
  )

  if (!.is_cv_folds_spec(cv_folds)) {
    if (is.null(train_idx)) {
      train_idx <- all_subjects
    } else {
      train_idx <- .validate_obs_indices(
        train_idx, length(subjects), "train_idx"
      )
    }
    guard_idx <- setdiff(all_subjects, train_idx)
    plan$contexts$single_fit <- .new_fit_context(
      data = data,
      role = "single_fit",
      axis = "none",
      train_subject_indices = train_idx,
      guard_subject_indices = guard_idx,
      train_observation_indices = .subject_obs_indices(data, train_idx),
      test_observation_indices = empty_obs,
      guard_observation_indices = empty_obs
    )
    return(plan)
  }

  if (identical(cv_folds$axis, "observation")) {
    normalized <- .normalize_observation_folds(data, cv_folds)
    plan$axis <- "observation"
    plan$strategy <- cv_folds$method %||% "custom"
    plan$folds <- normalized$folds
    plan$per_subject_observation_folds <- normalized$per_subject
    plan$contexts$full_fit <- .new_fit_context(
      data = data,
      role = "full_fit",
      axis = "observation",
      train_subject_indices = all_subjects,
      train_observation_indices = all_obs,
      test_observation_indices = empty_obs,
      guard_observation_indices = empty_obs
    )
    return(plan)
  }

  .validate_cv_folds_spec(cv_folds, n_subjects = length(subjects))
  validate_cv_setup(cv_folds, reference = "medoid")
  fold_names <- names(cv_folds$folds)
  if (is.null(fold_names) || anyNA(fold_names) || any(!nzchar(fold_names)) ||
      anyDuplicated(fold_names)) {
    stop("Subject cv_folds$folds must be uniquely named", call. = FALSE)
  }

  plan$axis <- "subject"
  plan$strategy <- cv_folds$method %||% "custom"
  plan$folds <- lapply(fold_names, function(fold_name) {
    fold <- cv_folds$folds[[fold_name]]
    train_idx <- as.integer(fold$train)
    test_idx <- as.integer(fold$test)
    guard_idx <- setdiff(all_subjects, union(train_idx, test_idx))
    .new_fit_context(
      data = data,
      role = "subject_cv_fold",
      axis = "subject",
      fold_id = fold_name,
      train_subject_indices = train_idx,
      test_subject_indices = test_idx,
      guard_subject_indices = guard_idx,
      train_observation_indices = .subject_obs_indices(data, train_idx),
      test_observation_indices = .subject_obs_indices(data, test_idx),
      guard_observation_indices = empty_obs
    )
  })
  names(plan$folds) <- fold_names
  plan$assessment_contexts <- lapply(fold_names, function(fold_name) {
    fold <- cv_folds$folds[[fold_name]]
    test_idx <- as.integer(fold$test)
    contexts <- lapply(test_idx, function(subject_idx) {
      .new_fit_context(
        data = data,
        role = "subject_cv_assessment_fit",
        axis = "subject",
        fold_id = fold_name,
        train_subject_indices = subject_idx,
        guard_subject_indices = setdiff(all_subjects, subject_idx),
        train_observation_indices = .subject_obs_indices(data, subject_idx),
        test_observation_indices = empty_obs,
        guard_observation_indices = empty_obs
      )
    })
    names(contexts) <- subjects[test_idx]
    contexts
  })
  names(plan$assessment_contexts) <- fold_names
  plan$contexts$anchor_fit <- .new_fit_context(
    data = data,
    role = "anchor_fit",
    axis = "subject",
    train_subject_indices = all_subjects,
    train_observation_indices = all_obs,
    test_observation_indices = empty_obs,
    guard_observation_indices = empty_obs
  )
  plan
}

.prepare_aligner_runtime <- function(aligner, data, reference,
                                     resampling_plan, dots) {
  aligner$resampling_plan <- resampling_plan
  aligner$provider_plan <- NULL

  if (!identical(as.integer(aligner$api_version %||% 1L), 2L)) {
    return(aligner)
  }

  if (!is.null(aligner$prepare_fn)) {
    args <- c(list(
      data = data,
      reference = reference,
      resampling_plan = .copy_provider_value(resampling_plan)
    ), dots)
    aligner$provider_plan <- do.call(aligner$prepare_fn, args)
  }
  aligner
}

.provider_fit_context <- function(aligner, fit_role, fold_id = NULL,
                                  subject_id = NULL) {
  plan <- aligner$resampling_plan %||% NULL
  if (is.null(plan)) {
    stop("API-v2 aligner runtime is missing its resampling plan", call. = FALSE)
  }

  context <- if (identical(fit_role, "subject_cv_assessment_fit")) {
    plan$assessment_contexts[[fold_id]][[subject_id]] %||% NULL
  } else if (!is.null(fold_id)) {
    plan$folds[[fold_id]] %||% NULL
  } else {
    plan$contexts[[fit_role]] %||% NULL
  }
  if (is.null(context) || !identical(context$role, fit_role)) {
    stop(sprintf(
      "API-v2 aligner runtime has no '%s' context%s",
      fit_role,
      if (is.null(fold_id)) "" else sprintf(" for fold '%s'", fold_id)
    ), call. = FALSE)
  }
  context
}

.invoke_aligner_fit <- function(aligner,
                                data,
                                reference,
                                train_idx,
                                fit_role,
                                fold_id = NULL,
                                subject_id = NULL,
                                ...) {
  args <- c(list(
    data = data,
    reference = reference,
    train_idx = train_idx
  ), list(...))

  if (identical(as.integer(aligner$api_version %||% 1L), 2L)) {
    context <- .provider_fit_context(
      aligner, fit_role, fold_id, subject_id = subject_id
    )
    args$fit_context <- .copy_provider_value(context)
    args$provider_plan <- .copy_provider_value(aligner$provider_plan)
  }

  do.call(aligner$fit_fn, args)
}
