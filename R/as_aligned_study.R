#' Coerce AlignmentResult to AlignedStudy
#'
#' Converts algorithm-facing aligned matrices (features × observations) into
#' an analysis-facing [AlignedStudy] (observations × shared features).
#'
#' If `result` is a CV result without a common anchor
#' (`cv_info$anchor_common` is not `TRUE`), this errors and directs you to
#' [as_aligned_resample_set()].
#'
#' @param result An [AlignmentResult].
#' @param source_data Optional [AlignmentData] used to recover observation
#'   labels / metadata.
#' @param observation_data Optional named list of per-subject `data.frame`s,
#'   or a single `data.frame` shared across subjects with equal \(N\).
#' @param subject_data Optional subject-level `data.frame`.
#' @param representation Optional representation label for all blocks.
#' @param mode Alignment-use mode (see [analysis_safety_record()]). When
#'   `NULL`, a common-anchor CV result is treated as `"cross_fitted"` and a
#'   non-CV result as `"in_sample"`.
#' @param shared_space Optional [SharedFeatureSpace]; inferred from model if
#'   `NULL`.
#' @param alignment_role Role stamped on blocks.
#' @param ... Unused; for future extensions.
#'
#' @return An [AlignedStudy].
#'
#' @family aligned_study
#' @export
as_aligned_study <- function(result,
                             source_data = NULL,
                             observation_data = NULL,
                             subject_data = NULL,
                             representation = NA_character_,
                             mode = NULL,
                             shared_space = NULL,
                             alignment_role = "analysis",
                             ...) {
  result <- .ensure_result(result, what = "result")
  cv <- result@cv_info
  if (.cv_has_fold_specific_spaces(cv)) {
    stop(
      "Result has fold-specific shared spaces (cv_info$anchor_common is not TRUE). ",
      "Use as_aligned_resample_set(result, ...) instead of as_aligned_study().",
      call. = FALSE
    )
  }

  model <- result@model
  aligned <- result@aligned
  if (!length(aligned)) {
    if (identical(cv$axis, "observation") ||
        isTRUE(cv$deployment_refit)) {
      stop(
        "Observation-axis CV does not attach evaluation matrices to the deployment refit. ",
        "Use as_aligned_resample_set(result, ...) for fold assessments.",
        call. = FALSE
      )
    }
    stop("AlignmentResult has no aligned matrices; refit with return_aligned=TRUE",
         call. = FALSE)
  }

  if (is.null(shared_space)) {
    k <- nrow(as.matrix(aligned[[1L]]))
    shared_space <- shared_feature_space_from_model(model, dimension = k)
  }
  if (!inherits(shared_space, "SharedFeatureSpace")) {
    stop("'shared_space' must be a SharedFeatureSpace", call. = FALSE)
  }
  expected_coordinate <- .model_coordinate_id(model)
  if (!identical(shared_space$coordinate_id, expected_coordinate)) {
    stop(
      "'shared_space' does not identify the AlignmentResult model's coordinate artifact",
      call. = FALSE
    )
  }

  inferred_mode <- if (.is_cross_fitted_cv(cv)) "cross_fitted" else "in_sample"
  if (is.null(mode)) {
    mode <- inferred_mode
  }

  if (!is.null(source_data) && !inherits(source_data, "AlignmentData")) {
    source_data <- as_alignment_data(source_data)
  }

  blocks <- list()
  for (subj in names(aligned)) {
    Z_alg <- aligned[[subj]]
    values <- to_analysis_matrix(Z_alg)
    obs <- .resolve_observation_data(
      subject = subj,
      n_obs = nrow(values),
      observation_data = observation_data,
      source_data = source_data
    )
    blocks[[subj]] <- AlignedBlock(
      values = values,
      observation_data = obs,
      subject_id = subj,
      representation = representation,
      alignment_role = alignment_role,
      shared_space_id = shared_space$id,
      source_ref = if (!is.null(source_data)) {
        list(space = source_data@space, subject = subj)
      } else {
        NULL
      }
    )
  }

  verification <- if (identical(mode, "cross_fitted")) {
    .verify_common_anchor_crossfit(result)
  } else {
    NULL
  }
  safety_args <- list(
    mode = mode,
    model_fit_subjects = model@train_subjects,
    application_subjects = names(aligned),
    cross_fitted = identical(mode, "cross_fitted"),
    inductive = identical(mode, "inductive_calibration"),
    leakage_status = .leakage_status_for_mode(mode)
  )
  if (!is.null(verification)) {
    safety <- do.call(
      .verified_analysis_safety_record,
      c(safety_args, list(
        evidence = verification$evidence,
        checks = verification$checks
      ))
    )
  } else {
    safety <- do.call(analysis_safety_record, safety_args)
  }

  lineage <- list(
    created_by = "as_aligned_study",
    method = model@method,
    cv = cv[intersect(names(cv), c("method", "axis", "n_folds", "anchor_common"))],
    representation = representation,
    orientation = list(
      algorithm = "features_x_observations",
      analysis = "observations_x_shared_features",
      boundary = "to_analysis_matrix()"
    )
  )

  AlignedStudy(
    blocks = blocks,
    shared_space = shared_space,
    model = model,
    subject_data = subject_data,
    lineage = lineage,
    safety = safety,
    storage = list(mode = "eager"),
    metadata = list(quality = result@quality)
  )
}


.is_cross_fitted_cv <- function(cv) {
  is.list(cv) &&
    !identical(cv$method %||% "none", "none") &&
    !identical(cv$method %||% "none", "applied") &&
    isTRUE(cv$anchor_common) &&
    (identical(cv$axis, "subject") || identical(cv$axis, "observation"))
}


.verify_common_anchor_crossfit <- function(result) {
  cv <- result@cv_info
  if (!.is_cross_fitted_cv(cv) || !is.list(cv$folds) || !length(cv$folds)) {
    return(NULL)
  }
  aligned <- result@aligned
  if (!length(aligned) || is.null(names(aligned))) return(NULL)

  axis <- cv$axis
  if (identical(axis, "subject")) {
    subjects <- result@model@train_subjects
    counts <- setNames(integer(length(subjects)), subjects)
    disjoint <- TRUE
    valid_indices <- TRUE
    for (fold in cv$folds) {
      train <- fold$train %||% integer(0)
      test <- fold$test %||% integer(0)
      valid_indices <- valid_indices &&
        all(train %in% seq_along(subjects)) &&
        all(test %in% seq_along(subjects))
      disjoint <- disjoint && !length(intersect(train, test))
      if (all(test %in% seq_along(subjects))) {
        counts[subjects[test]] <- counts[subjects[test]] + 1L
      }
    }
    checks <- list(
      common_anchor = isTRUE(cv$anchor_common),
      valid_fold_indices = isTRUE(valid_indices),
      train_test_disjoint = isTRUE(disjoint),
      each_output_assessed_once = identical(
        unname(counts[names(aligned)]),
        rep.int(1L, length(aligned))
      )
    )
  } else {
    subjects <- names(aligned)
    coverage <- setNames(vector("list", length(subjects)), subjects)
    disjoint <- TRUE
    valid_indices <- TRUE
    for (subject in subjects) coverage[[subject]] <- integer(0)
    for (fold in cv$folds) {
      per_subject <- is.list(fold) && !is.null(names(fold)) &&
        is.null(fold$train_idx)
      for (subject in subjects) {
        part <- if (per_subject) fold[[subject]] else fold
        train <- part$train_idx %||% integer(0)
        test <- part$test_idx %||% integer(0)
        n_obs <- ncol(aligned[[subject]])
        valid_indices <- valid_indices &&
          all(train %in% seq_len(n_obs)) && all(test %in% seq_len(n_obs))
        disjoint <- disjoint && !length(intersect(train, test))
        coverage[[subject]] <- c(coverage[[subject]], test)
      }
    }
    exact_coverage <- vapply(subjects, function(subject) {
      idx <- coverage[[subject]]
      n_obs <- ncol(aligned[[subject]])
      length(idx) == n_obs && !anyDuplicated(idx) &&
        setequal(idx, seq_len(n_obs))
    }, logical(1))
    checks <- list(
      common_anchor = isTRUE(cv$anchor_common),
      valid_fold_indices = isTRUE(valid_indices),
      train_test_disjoint = isTRUE(disjoint),
      each_output_assessed_once = all(exact_coverage)
    )
  }

  if (!all(vapply(checks, isTRUE, logical(1)))) return(NULL)
  list(
    evidence = list(
      source = "AlignmentResult@cv_info$folds",
      method = cv$method,
      axis = cv$axis,
      n_folds = length(cv$folds)
    ),
    checks = checks
  )
}


.resolve_observation_data <- function(subject, n_obs, observation_data, source_data) {
  if (!is.null(observation_data)) {
    if (is.data.frame(observation_data)) {
      if (nrow(observation_data) != n_obs) {
        stop(sprintf(
          "observation_data has %d rows but subject '%s' has %d observations",
          nrow(observation_data), subject, n_obs
        ), call. = FALSE)
      }
      return(observation_data)
    }
    if (is.list(observation_data)) {
      if (!subject %in% names(observation_data)) {
        stop(sprintf("observation_data has no entry for subject '%s'", subject),
             call. = FALSE)
      }
      obs <- observation_data[[subject]]
      if (!is.data.frame(obs)) {
        stop("Each observation_data entry must be a data.frame", call. = FALSE)
      }
      if (nrow(obs) != n_obs) {
        stop(sprintf(
          "observation_data[['%s']] has %d rows but aligned data has %d",
          subject, nrow(obs), n_obs
        ), call. = FALSE)
      }
      return(obs)
    }
    stop("'observation_data' must be a data.frame or named list of data.frames",
         call. = FALSE)
  }

  # Recover from AlignmentData obs_labels when possible
  if (!is.null(source_data)) {
    labs <- .subject_obs_labels(source_data, subject)
    if (!is.null(labs)) {
      if (length(labs) != n_obs) {
        stop(sprintf(
          "obs_labels length (%d) does not match aligned observations (%d) for '%s'",
          length(labs), n_obs, subject
        ), call. = FALSE)
      }
      return(data.frame(
        observation_id = as.character(labs),
        stringsAsFactors = FALSE
      ))
    }
  }

  data.frame(
    observation_id = seq_len(n_obs),
    stringsAsFactors = FALSE
  )
}


.subject_obs_labels <- function(data, subject) {
  labs <- data@obs_labels
  if (is.null(labs)) return(NULL)
  if (is.atomic(labs) || is.factor(labs)) {
    return(as.character(labs))
  }
  if (is.list(labs)) {
    if (subject %in% names(labs)) {
      return(as.character(labs[[subject]]))
    }
  }
  NULL
}
