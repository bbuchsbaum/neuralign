#' AlignedResampleSet Class
#'
#' Collection of retained fold artifacts. Use this when each fold may define a
#' different aligned coordinate system. Downstream analyses train and assess
#' within a split, then aggregate metrics or predictions rather than raw shared
#' coordinates.
#'
#' @slot resampling_spec Description of the resampling scheme.
#' @slot splits Named list of split records. Each split contains the exact
#'   `model` used for that fold, its `shared_space`, an optional analysis
#'   [AlignedStudy], an assessment [AlignedStudy], and metadata.
#' @slot aggregation_contract Declares what may be aggregated across splits.
#' @slot metadata Extra metadata.
#' @param object,x An `AlignedResampleSet` object.
#'
#' @family aligned_study
#' @export
setClass("AlignedResampleSet",
  slots = c(
    resampling_spec = "list",
    splits = "list",
    aggregation_contract = "list",
    metadata = "list"
  ),
  prototype = list(
    resampling_spec = list(),
    splits = list(),
    aggregation_contract = list(
      aggregate = c("metrics", "predictions"),
      forbid = c("raw_coordinates", "stack_shared_features")
    ),
    metadata = list()
  )
)

setValidity("AlignedResampleSet", function(object) {
  errors <- character()
  split_names <- names(object@splits)
  if (!is.list(object@splits) ||
      (length(object@splits) &&
       (is.null(split_names) || anyNA(split_names) ||
        any(!nzchar(split_names)) || anyDuplicated(split_names)))) {
    errors <- c(errors, "'splits' must be a named list with unique names")
  }
  for (nm in split_names) {
    split <- object@splits[[nm]]
    if (!is.list(split)) {
      errors <- c(errors, sprintf("splits[['%s']] must be a list", nm))
      next
    }
    if (!inherits(split$model, "AlignmentModel")) {
      errors <- c(errors, sprintf(
        "splits[['%s']]$model must be an AlignmentModel", nm
      ))
    }
    if (!inherits(split$shared_space, "SharedFeatureSpace")) {
      errors <- c(errors, sprintf(
        "splits[['%s']]$shared_space must be a SharedFeatureSpace", nm
      ))
    }
    if (!is.null(split$analysis) &&
        !inherits(split$analysis, "AlignedStudy")) {
      errors <- c(errors, sprintf(
        "splits[['%s']]$analysis must be an AlignedStudy or NULL", nm
      ))
    }
    if (!inherits(split$assessment, "AlignedStudy")) {
      errors <- c(errors, sprintf(
        "splits[['%s']]$assessment must be an AlignedStudy", nm
      ))
    }
    if (!is.list(split$metadata)) {
      errors <- c(errors, sprintf(
        "splits[['%s']]$metadata must be a list", nm
      ))
    }
    if (inherits(split$shared_space, "SharedFeatureSpace")) {
      sid <- split$shared_space$id
      if (inherits(split$model, "AlignmentModel") &&
          !identical(
            split$shared_space$coordinate_id,
            .model_coordinate_id(split$model)
          )) {
        errors <- c(errors, sprintf(
          "splits[['%s']] model and shared space identify different coordinates",
          nm
        ))
      }
      for (role in c("analysis", "assessment")) {
        study <- split[[role]]
        if (inherits(study, "AlignedStudy") &&
            !identical(study@shared_space$id, sid)) {
          errors <- c(errors, sprintf(
            "splits[['%s']]$%s does not use its split shared space", nm, role
          ))
        }
      }
    }
  }
  if (length(errors) == 0L) TRUE else errors
})


#' Create an AlignedResampleSet
#'
#' @param splits Named list of retained split records.
#' @param resampling_spec Resampling description.
#' @param aggregation_contract Aggregation policy list.
#' @param metadata Extra metadata.
#' @return An [AlignedResampleSet].
#' @family aligned_study
#' @export
AlignedResampleSet <- function(splits,
                               resampling_spec = list(),
                               aggregation_contract = list(
                                 aggregate = c("metrics", "predictions"),
                                 forbid = c("raw_coordinates", "stack_shared_features")
                               ),
                               metadata = list()) {
  if (!is.list(splits) || (length(splits) && is.null(names(splits)))) {
    stop("'splits' must be a named list", call. = FALSE)
  }
  if (!is.list(resampling_spec) || !is.list(aggregation_contract) ||
      !is.list(metadata)) {
    stop(
      "resampling_spec, aggregation_contract, and metadata must be lists",
      call. = FALSE
    )
  }
  new("AlignedResampleSet",
    resampling_spec = resampling_spec,
    splits = splits,
    aggregation_contract = aggregation_contract,
    metadata = metadata
  )
}


#' @rdname AlignedResampleSet-class
#' @export
setMethod("show", "AlignedResampleSet", function(object) {
  cat("AlignedResampleSet\n")
  cat(sprintf("  splits: %d\n", length(object@splits)))
  if (!is.null(object@resampling_spec$method)) {
    cat(sprintf("  method: %s\n", object@resampling_spec$method))
  }
  if (!is.null(object@resampling_spec$axis)) {
    cat(sprintf("  axis: %s\n", object@resampling_spec$axis))
  }
  invisible(object)
})


#' @rdname AlignedResampleSet-class
#' @export
setMethod("length", "AlignedResampleSet", function(x) length(x@splits))


#' Convert Retained Cross-Validation Artifacts
#'
#' An [AlignmentResult] must have been fitted with
#' `return_resample_artifacts = TRUE`. Observation-axis conversion accepts the
#' retained result from [run_obs_crossfit_from_data()]. No fold models or
#' aligned matrices are reconstructed from summary metadata.
#'
#' @param result An [AlignmentResult] or `ObsCrossfitAlignment`.
#' @param source_data Optional [AlignmentData] used for observation labels.
#' @param observation_data Optional observation tables.
#' @param representation Optional client-defined representation label.
#' @param mode Safety mode stamped on assessment blocks. The current contract
#'   requires `"cross_fitted"`.
#' @param ... Unused.
#' @return An [AlignedResampleSet].
#' @family aligned_study
#' @export
as_aligned_resample_set <- function(result,
                                    source_data = NULL,
                                    observation_data = NULL,
                                    representation = NA_character_,
                                    mode = "cross_fitted",
                                    ...) {
  if (!identical(mode, "cross_fitted")) {
    stop(
      "AlignedResampleSet assessment safety is derived as 'cross_fitted'; ",
      "caller-selected safety modes are not accepted",
      call. = FALSE
    )
  }
  if (inherits(result, "ObsCrossfitAlignment")) {
    return(.resample_set_from_obs_crossfit(
      result,
      representation = representation,
      mode = mode
    ))
  }

  result <- .ensure_result(result, what = "result")
  cv <- result@cv_info
  artifacts <- cv$artifacts_by_fold %||% NULL
  if (!is.list(artifacts) || !length(artifacts)) {
    stop(
      "AlignmentResult does not retain fold artifacts; refit subject CV with ",
      "return_resample_artifacts=TRUE",
      call. = FALSE
    )
  }
  if (!is.null(source_data) && !inherits(source_data, "AlignmentData")) {
    source_data <- as_alignment_data(source_data)
  }

  splits <- lapply(names(artifacts), function(fold_id) {
    artifact <- artifacts[[fold_id]]
    .validate_subject_fold_artifact(artifact, fold_id)
    model <- artifact$model
    assessment_subjects <- artifact$assessment_subjects
    analysis_subjects <- artifact$analysis_subjects
    k <- nrow(artifact$aligned_assessment[[1L]])
    space <- shared_feature_space_from_model(
      model,
      dimension = k,
      extras = list(fold_id = fold_id, axis = "subject")
    )

    checks <- list(
      retained_fold_model = inherits(model, "AlignmentModel"),
      analysis_assessment_disjoint =
        !length(intersect(analysis_subjects, assessment_subjects)),
      model_training_membership =
        setequal(model@train_subjects, analysis_subjects),
      assessment_membership =
        setequal(names(artifact$aligned_assessment), assessment_subjects)
    )
    assessment_safety <- .verified_analysis_safety_record(
      mode = mode,
      model_fit_subjects = analysis_subjects,
      application_subjects = assessment_subjects,
      cross_fitted = TRUE,
      leakage_status = "verified_subject_crossfit",
      evidence = list(source = "AlignmentResult@cv_info$artifacts_by_fold",
                      fold_id = fold_id),
      checks = checks
    )

    analysis <- .aligned_study_from_matrices(
      matrices = artifact$aligned_analysis,
      model = model,
      space = space,
      role = "analysis",
      representation = representation,
      observation_data = observation_data,
      source_data = source_data,
      safety = analysis_safety_record(
        mode = "in_sample",
        model_fit_subjects = analysis_subjects,
        application_subjects = analysis_subjects,
        leakage_status = "declared_in_sample"
      ),
      fold_id = fold_id
    )
    assessment <- .aligned_study_from_matrices(
      matrices = artifact$aligned_assessment,
      model = model,
      space = space,
      role = "assessment",
      representation = representation,
      observation_data = observation_data,
      source_data = source_data,
      safety = assessment_safety,
      fold_id = fold_id
    )

    list(
      model = model,
      shared_space = space,
      analysis = analysis,
      assessment = assessment,
      metadata = list(
        fold_id = fold_id,
        analysis_subjects = analysis_subjects,
        assessment_subjects = assessment_subjects,
        train_indices = artifact$train_indices,
        test_indices = artifact$test_indices
      )
    )
  })
  names(splits) <- names(artifacts)

  AlignedResampleSet(
    splits = splits,
    resampling_spec = list(
      method = cv$method,
      axis = cv$axis %||% "subject",
      n_folds = length(splits),
      anchor_common = isTRUE(cv$anchor_common)
    ),
    metadata = list(parent_method = result@model@method)
  )
}


#' Forbid Raw Stacking Across Fold Spaces
#'
#' @param x An [AlignedResampleSet].
#' @return Never returns successfully; always errors.
#' @family aligned_study
#' @export
stack_aligned_resamples <- function(x) {
  if (!inherits(x, "AlignedResampleSet")) {
    stop("'x' must be an AlignedResampleSet", call. = FALSE)
  }
  stop(
    "Cannot stack raw aligned coordinates from AlignedResampleSet; ",
    "aggregate fold metrics or predictions instead",
    call. = FALSE
  )
}


#' @export
subject_ids.AlignedResampleSet <- function(x) {
  unique(unlist(lapply(x@splits, function(split) {
    c(
      if (!is.null(split$analysis)) subject_ids(split$analysis),
      subject_ids(split$assessment)
    )
  }), use.names = FALSE))
}


#' @export
shared_space.AlignedResampleSet <- function(x) {
  stop(
    "AlignedResampleSet has per-split shared spaces; inspect a split explicitly",
    call. = FALSE
  )
}


.validate_subject_fold_artifact <- function(artifact, fold_id) {
  required <- c(
    "model", "aligned_analysis", "aligned_assessment",
    "analysis_subjects", "assessment_subjects"
  )
  missing <- setdiff(required, names(artifact))
  if (length(missing)) {
    stop(sprintf(
      "Retained fold '%s' is missing: %s",
      fold_id, paste(missing, collapse = ", ")
    ), call. = FALSE)
  }
  if (!inherits(artifact$model, "AlignmentModel") ||
      !is.list(artifact$aligned_analysis) ||
      !is.list(artifact$aligned_assessment) ||
      !length(artifact$aligned_assessment)) {
    stop(sprintf("Retained fold '%s' is invalid", fold_id), call. = FALSE)
  }
  invisible(artifact)
}


.aligned_study_from_matrices <- function(matrices,
                                         model,
                                         space,
                                         role,
                                         representation,
                                         observation_data,
                                         source_data,
                                         safety,
                                         fold_id) {
  blocks <- lapply(names(matrices), function(subject) {
    values <- to_analysis_matrix(matrices[[subject]])
    obs <- .resolve_observation_data(
      subject = subject,
      n_obs = nrow(values),
      observation_data = observation_data,
      source_data = source_data
    )
    AlignedBlock(
      values = values,
      observation_data = obs,
      subject_id = subject,
      representation = representation,
      alignment_role = role,
      shared_space_id = space$id
    )
  })
  names(blocks) <- names(matrices)
  AlignedStudy(
    blocks = blocks,
    shared_space = space,
    model = model,
    lineage = list(
      created_by = "as_aligned_resample_set",
      fold_id = fold_id,
      role = role
    ),
    safety = safety
  )
}


.resample_set_from_obs_crossfit <- function(result,
                                            representation,
                                            mode) {
  models <- result$models_by_fold
  aligned <- result$aligned_test_by_fold
  if (!is.list(models) || !length(models) ||
      !is.list(aligned) || !length(aligned) ||
      !identical(names(models), names(aligned))) {
    stop(
      "ObsCrossfitAlignment must retain matching models_by_fold and ",
      "aligned_test_by_fold artifacts",
      call. = FALSE
    )
  }
  test_labels <- result$obs_labels_test_by_fold %||% NULL
  train_labels <- result$obs_labels_train_by_fold %||% NULL
  fold_specs <- result$fold_info$obs_folds$folds %||% NULL

  splits <- lapply(names(models), function(fold_id) {
    model <- models[[fold_id]]
    matrices <- aligned[[fold_id]]
    if (!inherits(model, "AlignmentModel") || !is.list(matrices) ||
        !length(matrices)) {
      stop(sprintf("Observation fold '%s' is invalid", fold_id), call. = FALSE)
    }
    k <- nrow(matrices[[1L]])
    space <- shared_feature_space_from_model(
      model,
      dimension = k,
      extras = list(fold_id = fold_id, axis = "observation")
    )
    obs <- NULL
    fold <- fold_specs[[fold_id]] %||% NULL
    if (!is.null(fold)) {
      per_subject <- is.list(fold) && !is.null(names(fold)) &&
        is.null(fold$train_idx)
      obs <- lapply(names(matrices), function(subject) {
        test_idx <- if (per_subject) {
          fold[[subject]]$test_idx
        } else {
          fold$test_idx
        }
        labels <- test_labels[[fold_id]][[subject]] %||% test_idx
        data.frame(
          observation_id = as.character(labels),
          source_index = as.integer(test_idx),
          stringsAsFactors = FALSE
        )
      })
      names(obs) <- names(matrices)
    } else if (!is.null(test_labels) && !is.null(test_labels[[fold_id]])) {
      obs <- lapply(test_labels[[fold_id]], function(labels) {
        data.frame(
          observation_id = as.character(labels),
          stringsAsFactors = FALSE
        )
      })
    }

    checks <- NULL
    if (!is.null(train_labels) && !is.null(test_labels) &&
        !is.null(train_labels[[fold_id]]) &&
        !is.null(test_labels[[fold_id]])) {
      subjects <- names(matrices)
      checks <- list(
        retained_fold_model = TRUE,
        train_test_observations_disjoint = all(vapply(subjects, function(s) {
          !length(intersect(
            as.character(train_labels[[fold_id]][[s]]),
            as.character(test_labels[[fold_id]][[s]])
          ))
        }, logical(1)))
      )
    }
    safety <- if (!is.null(checks) && all(vapply(checks, isTRUE, logical(1)))) {
      .verified_analysis_safety_record(
        mode = mode,
        model_fit_subjects = model@train_subjects,
        application_subjects = names(matrices),
        cross_fitted = TRUE,
        leakage_status = "verified_observation_crossfit",
        evidence = list(source = "ObsCrossfitAlignment", fold_id = fold_id),
        checks = checks
      )
    } else {
      analysis_safety_record(
        mode = mode,
        model_fit_subjects = model@train_subjects,
        application_subjects = names(matrices),
        cross_fitted = TRUE,
        leakage_status = "declared_observation_crossfit"
      )
    }
    assessment <- .aligned_study_from_matrices(
      matrices = matrices,
      model = model,
      space = space,
      role = "assessment",
      representation = representation,
      observation_data = obs,
      source_data = NULL,
      safety = safety,
      fold_id = fold_id
    )
    list(
      model = model,
      shared_space = space,
      analysis = NULL,
      assessment = assessment,
      metadata = list(
        fold_id = fold_id,
        assessment_subjects = names(matrices)
      )
    )
  })
  names(splits) <- names(models)

  AlignedResampleSet(
    splits = splits,
    resampling_spec = list(
      method = result$fold_info$obs_folds$method %||% "custom",
      axis = "observation",
      n_folds = length(splits),
      anchor_common = isTRUE(result$anchor_common)
    ),
    metadata = list(parent_method = result$provenance$method)
  )
}
