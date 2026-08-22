#' Align a Study Into an AlignedStudy
#'
#' Apply an [AlignmentModel] (or extract one from an [AlignmentResult]) to
#' analysis-ready or algorithm-facing data and return an [AlignedStudy].
#'
#' @param model An [AlignmentModel] or [AlignmentResult].
#' @param data An [AlignmentData], named list of algorithm-facing matrices
#'   (features × observations), or named list of analysis-facing matrices
#'   when `orientation = "analysis"`.
#' @param mode Alignment-use mode: `"frozen_application"`, `"in_sample"`,
#'   `"cross_fitted"`, `"inductive_calibration"`, or `"outer_assessment"`.
#' @param orientation `"algorithm"` (features × observations, default) or
#'   `"analysis"` (observations × features).
#' @param representation Optional client-defined representation label.
#' @param observation_data Optional per-subject observation tables.
#' @param subject_data Optional subject table.
#' @param application_source_id Optional stable identifier for the observations
#'   being transformed. When omitted, `data@metadata$source_id` is used if
#'   present. Together with `fit_source_id` retained by [fit_alignment()], this
#'   can verify a truly frozen application for confirmatory group inference.
#' @param materialize Currently must be `TRUE` (lazy backends are reserved).
#' @param fit_new Passed to [apply_alignment()] when new subjects appear.
#' @param warn_leakage Passed to [apply_alignment()].
#' @param shared_space Optional [SharedFeatureSpace] override.
#' @param ... Passed to [apply_alignment()].
#'
#' @return An [AlignedStudy].
#'
#' @examples
#' set.seed(1)
#' train <- AlignmentData(list(
#'   s1 = matrix(rnorm(40), 8, 5),
#'   s2 = matrix(rnorm(40), 8, 5)
#' ))
#' fit <- fit_alignment(train, method = "procrustes", reference = "s1",
#'                      compute_quality = FALSE)
#' held_out <- list(
#'   s1 = matrix(rnorm(24), 8, 3),
#'   s2 = matrix(rnorm(24), 8, 3)
#' )
#' study <- align_study(fit, held_out, mode = "frozen_application",
#'                      representation = "derived_features")
#' aligned_matrix(study, "s1")
#'
#' @family aligned_study
#' @export
align_study <- function(model,
                        data,
                        mode = c(
                          "frozen_application",
                          "in_sample",
                          "cross_fitted",
                          "inductive_calibration",
                          "outer_assessment"
                        ),
                        orientation = c("algorithm", "analysis"),
                        representation = NA_character_,
                        observation_data = NULL,
                        subject_data = NULL,
                        application_source_id = NULL,
                        materialize = TRUE,
                        fit_new = TRUE,
                        warn_leakage = TRUE,
                        shared_space = NULL,
                        ...) {
  mode <- match.arg(mode)
  orientation <- match.arg(orientation)
  if (!isTRUE(materialize)) {
    stop("materialize=FALSE is reserved for a future lazy backend", call. = FALSE)
  }
  if (identical(mode, "cross_fitted")) {
    stop(
      "align_study() applies one fitted model and cannot verify cross-fitting; ",
      "convert a retained CV result with as_aligned_study() or ",
      "as_aligned_resample_set() instead",
      call. = FALSE
    )
  }

  if (inherits(model, "AlignmentResult")) {
    model <- get_model(model)
  }
  model <- .ensure_model(model, what = "model")

  if (identical(orientation, "analysis")) {
    data <- .analysis_list_to_alignment_data(data)
  } else if (!inherits(data, "AlignmentData")) {
    data <- as_alignment_data(data)
  }
  application_source_id <- .resolve_alignment_source_id(
    application_source_id,
    data,
    argument = "application_source_id"
  )

  # Space check against model source space (soft warning via apply_alignment)
  result <- apply_alignment(
    model,
    data,
    fit_new = fit_new,
    warn_leakage = warn_leakage,
    ...
  )

  # Applying a model may add subject transforms, but it must not define a new
  # coordinate system. Derive identity from the fitted artifact supplied by the
  # caller, before application-specific transforms are added.
  if (is.null(shared_space)) {
    k <- nrow(as.matrix(result@aligned[[1L]]))
    shared_space <- shared_feature_space_from_model(model, dimension = k)
  } else if (!inherits(shared_space, "SharedFeatureSpace") ||
             !identical(shared_space$coordinate_id, .model_coordinate_id(model))) {
    stop(
      "'shared_space' must identify the fitted model's coordinate artifact",
      call. = FALSE
    )
  }

  study <- as_aligned_study(
    result,
    source_data = data,
    observation_data = observation_data,
    subject_data = subject_data,
    representation = representation,
    mode = mode,
    shared_space = shared_space,
    alignment_role = if (identical(mode, "outer_assessment")) "assessment" else "analysis"
  )

  # Strengthen lineage for align_study path
  study@lineage$created_by <- "align_study"
  study@lineage$mode <- mode
  study@lineage$representation <- representation
  study@safety <- if (identical(mode, "frozen_application")) {
    .frozen_application_safety_record(
      model = model,
      application_subjects = data@subjects,
      application_source_id = application_source_id
    )
  } else {
    analysis_safety_record(
      mode = mode,
      model_fit_subjects = model@train_subjects,
      application_subjects = data@subjects,
      cross_fitted = identical(mode, "cross_fitted"),
      inductive = identical(mode, "inductive_calibration"),
      leakage_status = .leakage_status_for_mode(mode)
    )
  }
  study
}


#' Build an Analysis Safety Record
#'
#' @param mode Alignment-use mode.
#' @param model_fit_subjects Subjects used to fit the model.
#' @param model_fit_observations Optional observation ids used in fitting.
#' @param calibration_observations Optional calibration observation ids.
#' @param application_subjects Subjects transformed in this view.
#' @param application_observations Optional application observation ids.
#' @param outer_assessment_observations Optional held-out observation ids.
#' @param preprocessing_scope Optional preprocessing declaration.
#' @param cross_fitted Logical.
#' @param inductive Logical.
#' @param leakage_status Character status token.
#' @param extras Extra named fields.
#'
#' @return A list suitable for `AlignedStudy@safety`.
#'
#' @family aligned_study
#' @export
analysis_safety_record <- function(mode,
                                   model_fit_subjects = character(0),
                                   model_fit_observations = NULL,
                                   calibration_observations = NULL,
                                   application_subjects = character(0),
                                   application_observations = NULL,
                                   outer_assessment_observations = NULL,
                                   preprocessing_scope = NULL,
                                   cross_fitted = FALSE,
                                   inductive = FALSE,
                                   leakage_status = "unknown",
                                   extras = list()) {
  allowed <- c(
    "in_sample",
    "cross_fitted",
    "frozen_application",
    "inductive_calibration",
    "outer_assessment"
  )
  if (!is.character(mode) || length(mode) != 1L || is.na(mode) ||
      !mode %in% allowed) {
    stop(sprintf(
      "'mode' must be one of %s",
      paste(allowed, collapse = ", ")
    ), call. = FALSE)
  }
  if (!is.character(leakage_status) || length(leakage_status) != 1L ||
      is.na(leakage_status) || !nzchar(leakage_status)) {
    stop("'leakage_status' must be a non-empty character string",
         call. = FALSE)
  }
  if (!is.list(extras)) {
    stop("'extras' must be a list", call. = FALSE)
  }
  list(
    status = if (identical(mode, "in_sample") ||
                 identical(leakage_status, "declared_in_sample")) {
      "unsafe"
    } else if (identical(leakage_status, "unknown")) {
      "unknown"
    } else {
      "declared"
    },
    mode = mode,
    model_fit_subjects = as.character(model_fit_subjects),
    model_fit_observations = model_fit_observations,
    calibration_observations = calibration_observations,
    application_subjects = as.character(application_subjects),
    application_observations = application_observations,
    outer_assessment_observations = outer_assessment_observations,
    preprocessing_scope = preprocessing_scope,
    cross_fitted = isTRUE(cross_fitted),
    inductive = isTRUE(inductive),
    leakage_status = leakage_status,
    extras = extras
  )
}


#' Assert That an AlignedStudy Is Safe for a Purpose
#'
#' Machine-checkable subset of the safety contract. Declarative modes that
#' cannot be verified from available ids are accepted only when
#' `allow_declared = TRUE`.
#'
#' @param aligned An [AlignedStudy].
#' @param purpose One of `"exploratory"`,
#'   `"confirmatory_cross_subject_prediction"`,
#'   `"confirmatory_group_inference"`, or `"within_subject"`.
#' @param allow_declared Accept declared modes without source-id proof.
#'
#' @return Invisibly `TRUE`.
#'
#' @family aligned_study
#' @export
assert_analysis_safe <- function(aligned,
                                 purpose = c(
                                   "exploratory",
                                   "confirmatory_cross_subject_prediction",
                                   "confirmatory_group_inference",
                                   "within_subject"
                                 ),
                                 allow_declared = FALSE) {
  purpose <- match.arg(purpose)
  if (!inherits(aligned, "AlignedStudy")) {
    stop("'aligned' must be an AlignedStudy", call. = FALSE)
  }
  assert_common_shared_space(aligned)
  safety <- aligned@safety
  mode <- safety$mode %||% NA_character_
  status <- safety$status %||% "unknown"

  if (identical(purpose, "exploratory") || identical(purpose, "within_subject")) {
    return(invisible(TRUE))
  }

  # A mode label is a declaration, not proof. Confirmatory use requires a
  # record created by an internal verifier from retained split artifacts.
  if (identical(status, "verified_safe")) {
    return(invisible(TRUE))
  }
  if (identical(status, "declared") && isTRUE(allow_declared)) {
    warning(
      sprintf(
        "Safety mode '%s' is declared but not verified from retained artifacts",
        mode
      ),
      call. = FALSE
    )
    return(invisible(TRUE))
  }
  if (identical(status, "declared")) {
    stop(
      "Analysis safety is declared, not verified; strict confirmatory use ",
      "requires retained split evidence",
      call. = FALSE
    )
  }
  if (identical(status, "unsafe")) {
    purpose_label <- switch(
      purpose,
      confirmatory_group_inference = "confirmatory group inference",
      confirmatory_cross_subject_prediction = "confirmatory cross-subject prediction",
      purpose
    )
    stop(sprintf(
      "AlignedStudy mode '%s' is unsafe for %s",
      mode, purpose_label
    ), call. = FALSE)
  }
  if (identical(status, "unknown")) {
    stop(
      "Analysis safety is unknown; retained split evidence is required",
      call. = FALSE
    )
  }
  stop(sprintf("Unknown analysis safety status '%s'", status), call. = FALSE)
}


.frozen_application_safety_record <- function(model,
                                               application_subjects,
                                               application_source_id) {
  fit_source_id <- model@provenance[["fit_source_id"]] %||% NULL
  base <- analysis_safety_record(
    mode = "frozen_application",
    model_fit_subjects = model@train_subjects,
    application_subjects = application_subjects,
    leakage_status = "declared_frozen_application",
    extras = list(
      fit_source_id = fit_source_id,
      application_source_id = application_source_id
    )
  )

  if (is.null(fit_source_id) || is.null(application_source_id)) {
    return(base)
  }

  checks <- list(
    distinct_source_ids = !identical(fit_source_id, application_source_id),
    pretrained_subject_transforms = all(
      application_subjects %in% names(model@transforms)
    )
  )
  evidence <- list(
    kind = "frozen_application",
    fit_source_id = fit_source_id,
    application_source_id = application_source_id,
    model_subjects = names(model@transforms),
    application_subjects = application_subjects
  )

  if (all(vapply(checks, isTRUE, logical(1)))) {
    return(.verified_analysis_safety_record(
      mode = "frozen_application",
      model_fit_subjects = model@train_subjects,
      application_subjects = application_subjects,
      leakage_status = "verified_frozen_application",
      extras = base$extras,
      evidence = evidence,
      checks = checks
    ))
  }

  base$status <- "unsafe"
  base$leakage_status <- if (!isTRUE(checks$distinct_source_ids)) {
    "verified_source_overlap"
  } else {
    "application_fitted_subject_transform"
  }
  base$verification <- list(evidence = evidence, checks = checks)
  base
}


.verified_analysis_safety_record <- function(..., evidence, checks) {
  if (missing(evidence) || is.null(evidence) ||
      missing(checks) || !is.list(checks) || !length(checks) ||
      !all(vapply(checks, isTRUE, logical(1)))) {
    stop("Verified safety requires retained evidence and passing checks",
         call. = FALSE)
  }
  record <- analysis_safety_record(...)
  record$status <- "verified_safe"
  record$verification <- list(evidence = evidence, checks = checks)
  record
}


.leakage_status_for_mode <- function(mode) {
  switch(
    mode,
    outer_assessment = "declared_outer_assessment",
    frozen_application = "declared_frozen_application",
    cross_fitted = "declared_cross_fitted",
    inductive_calibration = "declared_inductive_calibration",
    in_sample = "declared_in_sample",
    "unknown"
  )
}


.analysis_list_to_alignment_data <- function(data) {
  if (inherits(data, "AlignmentData")) {
    # Interpret existing AlignmentData as already algorithm-facing
    return(data)
  }
  if (!is.list(data) || is.null(names(data))) {
    stop("For orientation='analysis', 'data' must be a named list of matrices",
         call. = FALSE)
  }
  alg <- lapply(data, function(x) {
    if (!.is_matrixish(x)) {
      stop("Each data element must be a matrix", call. = FALSE)
    }
    to_algorithm_matrix(x)
  })
  as_alignment_data(alg)
}
