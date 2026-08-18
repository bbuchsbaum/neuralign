#' Internal: Single Fit (No CV)
#' @keywords internal
.fit_single <- function(data, aligner, reference, train_idx,
                        compute_quality, return_aligned = TRUE,
                        restrict_to_identified = FALSE,
                        restrict_tol = sqrt(.Machine$double.eps),
                        ...) {
  # Determine training subjects
  if (is.null(train_idx)) {
    train_idx <- seq_along(data@subjects)
  } else {
    n <- length(data@subjects)
    if (length(train_idx) == 0L) {
      stop("'train_idx' must not be empty", call. = FALSE)
    }
    if (any(train_idx < 1L | train_idx > n)) {
      stop(sprintf("'train_idx' contains out-of-bounds indices (valid range: 1 to %d)", n),
        call. = FALSE)
    }
  }
  train_subjects <- data@subjects[train_idx]
  returns <- aligner$capabilities$returns %||% "operator"

  # Resolve reference
  ref_resolved <- .resolve_reference_spec(data, reference, train_idx)

  # Call the fit function
  fit_result <- .invoke_aligner_fit(
    aligner = aligner,
    data = data,
    reference = ref_resolved$reference,
    train_idx = train_idx,
    fit_role = "single_fit",
    ...
  )

  if (identical(returns, "embedding")) {
    aligned_full <- fit_result$aligned %||% NULL
    if (is.null(aligned_full)) {
      stop(
        sprintf("fit_alignment(%s): embedding-returning methods must return 'aligned'", aligner$name),
        call. = FALSE
      )
    }

    .validate_embedding_aligned(
      aligned = aligned_full,
      data_list = get_data_list(data),
      context = sprintf("fit_alignment(%s)", aligner$name)
    )

    transforms <- lapply(names(aligned_full), function(subj) {
      z <- if (isTRUE(return_aligned)) aligned_full[[subj]] else matrix(numeric(0), 0L, 0L)
      .new_embedding_transform(z, subject = subj)
    })
    names(transforms) <- names(aligned_full)

    model <- AlignmentModel(
      transforms = transforms,
      reference = ref_resolved$reference_spec,
      reference_data = fit_result$reference_data %||% NULL,
      method = aligner$name,
      space_from = fit_result$space_from %||% data@space,
      space_to = fit_result$space_to %||% data@space,
      params = list(...),
      method_state = fit_result$method_state %||% list(),
      train_subjects = train_subjects
    )

    aligned <- if (return_aligned) aligned_full else list()

    quality <- list()
    if (compute_quality && return_aligned) {
      quality <- .compute_basic_quality(data, aligned, model)
    }

    return(AlignmentResult(
      model = model,
      aligned = aligned,
      quality = quality,
      cv_info = list(method = "none")
    ))
  }

  .validate_operator_transforms(
    transforms = fit_result$transforms,
    data_list = get_data_list(data),
    context = sprintf("fit_alignment(%s)", aligner$name)
  )

  transforms <- fit_result$transforms
  if (isTRUE(restrict_to_identified)) {
    data_list <- get_data_list(data)
    transforms <- lapply(names(transforms), function(subj) {
      canonicalize_orthogonal_operator(
        Q = transforms[[subj]],
        x = data_list[[subj]],
        convention = "left",
        tol = restrict_tol
      )
    })
    names(transforms) <- names(fit_result$transforms)
  }

  # Build AlignmentModel
  params <- c(list(...), list(
    restrict_to_identified = isTRUE(restrict_to_identified),
    restrict_tol = restrict_tol
  ))
  method_state <- fit_result$method_state %||% list()
  method_state$restrict_to_identified <- isTRUE(restrict_to_identified)
  method_state$restrict_tol <- restrict_tol

  model <- AlignmentModel(
    transforms = transforms,
    reference = ref_resolved$reference_spec,
    reference_data = fit_result$reference_data,
    method = aligner$name,
    space_from = fit_result$space_from,
    space_to = fit_result$space_to,
    params = params,
    method_state = method_state,
    train_subjects = train_subjects
  )

  # Apply alignment to get aligned data
  aligned <- if (return_aligned) .apply_transforms(model, data) else list()

  # Compute quality metrics if requested
  quality <- list()
  if (compute_quality && return_aligned) {
    quality <- .compute_basic_quality(data, aligned, model)
  }

  AlignmentResult(
    model = model,
    aligned = aligned,
    quality = quality,
    cv_info = list(method = "none")
  )
}


#' Internal: Leave-One-Subject-Out CV
#' @keywords internal
.fit_cv_loso <- function(data, aligner, reference, compute_quality,
                         return_aligned = TRUE,
                         return_resample_artifacts = FALSE,
                         restrict_to_identified = FALSE,
                         restrict_tol = sqrt(.Machine$double.eps),
                         ...) {
  cv_folds <- create_cv_folds(data, method = "loso")
  validate_cv_setup(cv_folds, reference = reference)
  returns <- aligner$capabilities$returns %||% "operator"
  if (identical(returns, "embedding")) {
    if (isTRUE(return_resample_artifacts)) {
      stop(
        "return_resample_artifacts is not yet supported for embedding-returning methods",
        call. = FALSE
      )
    }
    .fit_cv_folds_embedding(data, aligner, reference, cv_folds,
      compute_quality = compute_quality,
      return_aligned = return_aligned,
      ...
    )
  } else {
    .fit_cv_folds(data, aligner, reference, cv_folds,
      compute_quality = compute_quality,
      return_aligned = return_aligned,
      return_resample_artifacts = return_resample_artifacts,
      restrict_to_identified = restrict_to_identified,
      restrict_tol = restrict_tol,
      ...
    )
  }
}


#' Internal: K-Fold CV
#' @keywords internal
.fit_cv_kfold <- function(data, aligner, reference, k, compute_quality,
                          return_aligned = TRUE,
                          return_resample_artifacts = FALSE,
                          restrict_to_identified = FALSE,
                          restrict_tol = sqrt(.Machine$double.eps),
                          ...) {
  cv_folds <- create_cv_folds(data, method = "kfold", k = k)
  # create_cv_folds will fall back to LOSO if k > n.
  validate_cv_setup(cv_folds, reference = reference)
  returns <- aligner$capabilities$returns %||% "operator"
  if (identical(returns, "embedding")) {
    if (isTRUE(return_resample_artifacts)) {
      stop(
        "return_resample_artifacts is not yet supported for embedding-returning methods",
        call. = FALSE
      )
    }
    .fit_cv_folds_embedding(data, aligner, reference, cv_folds,
      compute_quality = compute_quality,
      return_aligned = return_aligned,
      ...
    )
  } else {
    .fit_cv_folds(data, aligner, reference, cv_folds,
      compute_quality = compute_quality,
      return_aligned = return_aligned,
      return_resample_artifacts = return_resample_artifacts,
      restrict_to_identified = restrict_to_identified,
      restrict_tol = restrict_tol,
      ...
    )
  }
}


#' Internal: Generic CV with Provided Folds
#' @keywords internal
.resolve_unevaluated_subjects <- function(cv_folds, subjects) {
  unevaluated_subjects <- as.character(cv_folds$unevaluated_subjects %||% character(0))
  if (length(unevaluated_subjects) == 0) return(character(0))

  if (anyNA(unevaluated_subjects) || any(!nzchar(unevaluated_subjects))) {
    stop("'cv_folds$unevaluated_subjects' must not contain NA/empty values", call. = FALSE)
  }
  unknown <- setdiff(unevaluated_subjects, subjects)
  if (length(unknown) > 0) {
    stop(
      sprintf(
        "'cv_folds$unevaluated_subjects' contains unknown subjects: %s",
        paste(unknown, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  unique(unevaluated_subjects)
}

.fit_cv_anchor_common <- function(data, aligner, reference, n_subjects, ...) {
  # For a fixed/external anchor, it's safe to store the corresponding reference_data.
  ref_resolved <- .resolve_reference_spec(data, reference, seq_len(n_subjects))
  fit_result_all <- .invoke_aligner_fit(
    aligner = aligner,
    data = data,
    reference = ref_resolved$reference,
    train_idx = seq_len(n_subjects),
    fit_role = "anchor_fit",
    ...
  )

  .validate_operator_transforms(
    transforms = fit_result_all$transforms,
    data_list = get_data_list(data),
    context = sprintf("fit_alignment(%s) [cv anchor]", aligner$name)
  )

  list(
    reference = ref_resolved$reference_spec,
    reference_data = fit_result_all$reference_data,
    method_state = fit_result_all$method_state %||% list(),
    space_from = fit_result_all$space_from,
    space_to = fit_result_all$space_to
  )
}

.fit_cv_anchor_common_embedding <- function(data, aligner, reference, n_subjects, ...) {
  ref_resolved <- .resolve_reference_spec(data, reference, seq_len(n_subjects))
  fit_result_all <- .invoke_aligner_fit(
    aligner = aligner,
    data = data,
    reference = ref_resolved$reference,
    train_idx = seq_len(n_subjects),
    fit_role = "anchor_fit",
    ...
  )

  aligned_full <- fit_result_all$aligned %||% NULL
  if (is.null(aligned_full)) {
    stop(
      sprintf("fit_alignment(%s) [cv anchor]: embedding-returning methods must return 'aligned'", aligner$name),
      call. = FALSE
    )
  }

  .validate_embedding_aligned(
    aligned = aligned_full,
    data_list = get_data_list(data),
    context = sprintf("fit_alignment(%s) [cv anchor]", aligner$name)
  )

  list(
    reference = ref_resolved$reference_spec,
    reference_data = fit_result_all$reference_data %||% NULL,
    method_state = fit_result_all$method_state %||% list(),
    space_from = fit_result_all$space_from %||% data@space,
    space_to = fit_result_all$space_to %||% data@space,
    aligned = aligned_full
  )
}

.cv_info_note_fold_specific <- function(anchor_common) {
  if (isTRUE(anchor_common)) return(NULL)
  "Aligned outputs are in fold-specific anchor spaces; do not use for group-level comparisons without mapping to a common anchor."
}

.fit_cv_folds <- function(data, aligner, reference, cv_folds,
                          compute_quality, return_aligned = TRUE,
                          return_resample_artifacts = FALSE,
                          restrict_to_identified = FALSE,
                          restrict_tol = sqrt(.Machine$double.eps),
                          ...) {
  .validate_cv_folds_spec(cv_folds, n_subjects = length(data@subjects))

  subjects <- data@subjects
  n_subjects <- length(subjects)
  unevaluated_subjects <- .resolve_unevaluated_subjects(cv_folds, subjects)

  # Check CV support
  if (!isTRUE(aligner$capabilities$supports_cv)) {
    warning(sprintf(
      "Method '%s' may not fully support CV; results may have leakage",
      aligner$name
    ), call. = FALSE)
  }

  reference_kind <- .reference_kind(reference)
  anchor_common <- reference_kind %in% c("fixed_subject", "template")
  retain_artifacts <- isTRUE(return_resample_artifacts) || !isTRUE(anchor_common)
  store_global_aligned <- isTRUE(return_aligned) && isTRUE(anchor_common)

  all_transforms <- list()
  all_aligned <- list()
  artifacts_by_fold <- if (isTRUE(retain_artifacts)) list() else NULL
  anchor_by_subject <- setNames(rep(NA_character_, n_subjects), subjects)

  # Fit/apply for each fold
  for (fold_name in names(cv_folds$folds)) {
    fold <- cv_folds$folds[[fold_name]]
    train_idx <- fold$train
    test_idx <- fold$test
    test_subjects <- subjects[test_idx]

    # Resolve reference using only training subjects
    ref_resolved <- .resolve_reference_spec(data, reference, train_idx)

    # Fit on training subjects
    fit_result <- .invoke_aligner_fit(
      aligner = aligner,
      data = data,
      reference = ref_resolved$reference,
      train_idx = train_idx,
      fit_role = "subject_cv_fold",
      fold_id = fold_name,
      ...
    )

    .validate_operator_transforms(
      transforms = fit_result$transforms,
      data_list = get_data_list(data)[subjects[train_idx]],
      context = sprintf("fit_alignment(%s) [cv train]", aligner$name)
    )

    if (isTRUE(restrict_to_identified)) {
      fit_result$transforms <- lapply(subjects[train_idx], function(subject) {
        canonicalize_orthogonal_operator(
          Q = fit_result$transforms[[subject]],
          x = get_subject_data(data, subject),
          convention = "left",
          tol = restrict_tol
        )
      })
      names(fit_result$transforms) <- subjects[train_idx]
    }

    fold_transforms <- fit_result$transforms
    fold_assessment_aligned <- list()

    # Apply to held-out subjects
    for (test_subj in test_subjects) {
      if (test_subj %in% names(all_transforms)) {
        stop(
          sprintf(
            "CV fold spec assigns subject '%s' to multiple test folds; expected each subject exactly once",
            test_subj
          ),
          call. = FALSE
        )
      }

      test_i <- match(test_subj, subjects)
      test_transform <- .fit_new_subject(
        aligner, fit_result, data, test_i, ref_resolved$reference,
        fold_id = fold_name
      )
      test_data <- get_subject_data(data, test_subj)
      if (isTRUE(restrict_to_identified)) {
        test_transform <- canonicalize_orthogonal_operator(
          Q = test_transform,
          x = test_data,
          convention = "left",
          tol = restrict_tol
        )
      }
      all_transforms[[test_subj]] <- test_transform
      fold_transforms[[test_subj]] <- test_transform

      if (store_global_aligned || isTRUE(retain_artifacts)) {
        test_aligned <- apply_transform(test_transform, test_data)
        fold_assessment_aligned[[test_subj]] <- test_aligned
        if (store_global_aligned) all_aligned[[test_subj]] <- test_aligned
      }
      anchor_by_subject[[test_subj]] <- as.character(ref_resolved$reference_spec)
    }

    if (isTRUE(retain_artifacts)) {
      train_subjects <- subjects[train_idx]
      fold_analysis_aligned <- lapply(train_subjects, function(subject) {
        apply_transform(
          fit_result$transforms[[subject]],
          get_subject_data(data, subject)
        )
      })
      names(fold_analysis_aligned) <- train_subjects
      fold_model <- AlignmentModel(
        transforms = fold_transforms,
        reference = ref_resolved$reference_spec,
        reference_data = fit_result$reference_data %||% NULL,
        method = aligner$name,
        space_from = fit_result$space_from %||% data@space,
        space_to = fit_result$space_to %||% data@space,
        params = c(list(...), list(
          restrict_to_identified = isTRUE(restrict_to_identified),
          restrict_tol = restrict_tol
        )),
        method_state = fit_result$method_state %||% list(),
        train_subjects = train_subjects
      )
      artifacts_by_fold[[fold_name]] <- list(
        model = fold_model,
        aligned_analysis = fold_analysis_aligned,
        aligned_assessment = fold_assessment_aligned,
        analysis_subjects = train_subjects,
        assessment_subjects = test_subjects,
        train_indices = train_idx,
        test_indices = test_idx,
        reference = ref_resolved$reference_spec
      )
    }
  }

  missing_subjects <- setdiff(subjects, names(all_transforms))
  remaining_missing <- setdiff(missing_subjects, unevaluated_subjects)
  if (length(remaining_missing) > 0) {
    stop(
      sprintf(
        "CV fold spec never assigns these subjects to a test fold: %s",
        paste(remaining_missing, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  if (length(missing_subjects) > 0) {
    for (subj in missing_subjects) {
      X <- get_subject_data(data, subj)
      all_transforms[[subj]] <- diag(nrow(X))
      if (store_global_aligned) {
        all_aligned[[subj]] <- X
      }
      anchor_by_subject[[subj]] <- as.character(reference)
    }
  }

  model_reference <- reference
  model_reference_data <- NULL
  model_method_state <- list()
  model_space_from <- data@space
  model_space_to <- data@space

  if (anchor_common) {
    anchor_fit <- .fit_cv_anchor_common(data, aligner, reference, n_subjects, ...)
    model_reference <- anchor_fit$reference
    model_reference_data <- anchor_fit$reference_data
    model_method_state <- anchor_fit$method_state
    model_space_from <- anchor_fit$space_from
    model_space_to <- anchor_fit$space_to
  } else {
    # Fold-specific anchors: do not pretend there is a single shared reference.
    model_reference <- "fold_specific"
  }

  params <- c(list(...), list(
    restrict_to_identified = isTRUE(restrict_to_identified),
    restrict_tol = restrict_tol
  ))
  model_method_state <- model_method_state %||% list()
  model_method_state$restrict_to_identified <- isTRUE(restrict_to_identified)
  model_method_state$restrict_tol <- restrict_tol

  model <- AlignmentModel(
    transforms = all_transforms,
    reference = model_reference,
    reference_data = model_reference_data,
    method = aligner$name,
    space_from = model_space_from,
    space_to = model_space_to,
    params = params,
    method_state = model_method_state,
    train_subjects = subjects
  )

  # Quality metrics
  quality <- list()
  if (compute_quality && return_aligned && isTRUE(anchor_common)) {
    quality <- .compute_basic_quality(data, all_aligned, model)
  }

  cv_method <- cv_folds$method %||% "custom"
  cv_axis <- cv_folds$axis %||% "subject"

  AlignmentResult(
    model = model,
    aligned = if (isTRUE(anchor_common) && isTRUE(return_aligned)) all_aligned else list(),
    quality = quality,
    cv_info = list(
      method = cv_method,
      axis = cv_axis,
      n_folds = cv_folds$n_folds %||% length(cv_folds$folds),
      fold_assignments = cv_folds$assignments %||% NULL,
      folds = cv_folds$folds,
      reference_kind = reference_kind,
      anchor_common = anchor_common,
      anchor_by_subject = anchor_by_subject,
      artifacts_by_fold = artifacts_by_fold,
      aligned_origin = "fold_assessment",
      deployment_refit = FALSE,
      artifact_kind = "subject_cv",
      anchor_note = .cv_info_note_fold_specific(anchor_common)
    )
  )
}

.fit_cv_folds_embedding <- function(data, aligner, reference, cv_folds,
                                    compute_quality, return_aligned = TRUE,
                                    ...) {
  .validate_cv_folds_spec(cv_folds, n_subjects = length(data@subjects))

  caps <- aligner$capabilities %||% list()
  returns <- caps$returns %||% "operator"
  if (!identical(returns, "embedding")) {
    stop(
      sprintf(".fit_cv_folds_embedding called for method '%s' with returns='%s'", aligner$name, returns),
      call. = FALSE
    )
  }

  if (is.null(aligner$apply_fn)) {
    stop(
      sprintf(
        "Method '%s' returns embeddings; CV requires an apply_fn to embed held-out subjects",
        aligner$name
      ),
      call. = FALSE
    )
  }

  if (isFALSE(caps$supports_new_subject %||% TRUE)) {
    stop(
      sprintf(
        "Method '%s' returns embeddings but does not support embedding held-out subjects (supports_new_subject=FALSE)",
        aligner$name
      ),
      call. = FALSE
    )
  }

  subjects <- data@subjects
  n_subjects <- length(subjects)
  unevaluated_subjects <- .resolve_unevaluated_subjects(cv_folds, subjects)

  if (!isTRUE(caps$supports_cv)) {
    warning(sprintf(
      "Method '%s' may not fully support CV; results may have leakage",
      aligner$name
    ), call. = FALSE)
  }

  all_transforms <- list()
  all_aligned <- list()
  anchor_by_subject <- setNames(rep(NA_character_, n_subjects), subjects)

  space_from_first <- NULL
  space_to_first <- NULL

  for (fold_name in names(cv_folds$folds)) {
    fold <- cv_folds$folds[[fold_name]]
    train_idx <- fold$train
    test_idx <- fold$test
    test_subjects <- subjects[test_idx]

    ref_resolved <- .resolve_reference_spec(data, reference, train_idx)

    fit_result <- .invoke_aligner_fit(
      aligner = aligner,
      data = data,
      reference = ref_resolved$reference,
      train_idx = train_idx,
      fit_role = "subject_cv_fold",
      fold_id = fold_name,
      ...
    )

    if (is.null(space_from_first)) {
      space_from_first <- fit_result$space_from %||% data@space
    }
    if (is.null(space_to_first)) {
      space_to_first <- fit_result$space_to %||% data@space
    }

    for (test_subj in test_subjects) {
      if (test_subj %in% names(all_transforms)) {
        stop(
          sprintf(
            "CV fold spec assigns subject '%s' to multiple test folds; expected each subject exactly once",
            test_subj
          ),
          call. = FALSE
        )
      }

      test_i <- match(test_subj, subjects)
      test_transform <- .fit_new_subject(
        aligner, fit_result, data, test_i, ref_resolved$reference,
        fold_id = fold_name
      )

      if (!.is_embedding_transform(test_transform)) {
        stop(
          sprintf(
            "Method '%s' did not produce an embedding transform for subject '%s' during CV",
            aligner$name, test_subj
          ),
          call. = FALSE
        )
      }

      all_transforms[[test_subj]] <- test_transform

      if (return_aligned) {
        test_data <- get_subject_data(data, test_subj)
        all_aligned[[test_subj]] <- apply_transform(test_transform, test_data)
      }

      anchor_by_subject[[test_subj]] <- as.character(ref_resolved$reference_spec)
    }
  }

  missing_subjects <- setdiff(subjects, names(all_transforms))
  remaining_missing <- setdiff(missing_subjects, unevaluated_subjects)
  if (length(remaining_missing) > 0) {
    stop(
      sprintf(
        "CV fold spec never assigns these subjects to a test fold: %s",
        paste(remaining_missing, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  reference_kind <- .reference_kind(reference)
  anchor_common <- reference_kind %in% c("fixed_subject", "template")

  anchor_fit <- NULL
  if (isTRUE(anchor_common) || length(missing_subjects) > 0) {
    anchor_fit <- .fit_cv_anchor_common_embedding(data, aligner, reference, n_subjects, ...)
  }

  if (length(missing_subjects) > 0) {
    aligned_full <- anchor_fit$aligned %||% NULL
    if (is.null(aligned_full)) {
      stop(
        "Cannot populate unevaluated subjects for embedding CV: full-fit 'aligned' is missing",
        call. = FALSE
      )
    }

    for (subj in missing_subjects) {
      z_full <- aligned_full[[subj]]
      z <- if (isTRUE(return_aligned)) z_full else matrix(numeric(0), 0L, 0L)
      all_transforms[[subj]] <- .new_embedding_transform(z, subject = subj)
      if (isTRUE(return_aligned)) {
        all_aligned[[subj]] <- z_full
      }
      anchor_by_subject[[subj]] <- as.character(reference)
    }
  }

  model_reference <- "fold_specific"
  model_reference_data <- NULL
  model_method_state <- list()
  model_space_from <- space_from_first %||% data@space
  model_space_to <- space_to_first %||% data@space

  if (isTRUE(anchor_common)) {
    model_reference <- anchor_fit$reference
    model_reference_data <- anchor_fit$reference_data
    model_method_state <- anchor_fit$method_state
    model_space_from <- anchor_fit$space_from
    model_space_to <- anchor_fit$space_to
  } else {
    all_aligned <- list()
  }

  model <- AlignmentModel(
    transforms = all_transforms,
    reference = model_reference,
    reference_data = model_reference_data,
    method = aligner$name,
    space_from = model_space_from,
    space_to = model_space_to,
    params = list(...),
    method_state = model_method_state,
    train_subjects = subjects
  )

  quality <- list()
  if (compute_quality && return_aligned && isTRUE(anchor_common)) {
    quality <- .compute_basic_quality(data, all_aligned, model)
  }

  cv_method <- cv_folds$method %||% "custom"
  cv_axis <- cv_folds$axis %||% "subject"

  AlignmentResult(
    model = model,
    aligned = if (isTRUE(anchor_common) && isTRUE(return_aligned)) all_aligned else list(),
    quality = quality,
    cv_info = list(
      method = cv_method,
      axis = cv_axis,
      n_folds = cv_folds$n_folds %||% length(cv_folds$folds),
      fold_assignments = cv_folds$assignments %||% NULL,
      folds = cv_folds$folds,
      reference_kind = reference_kind,
      anchor_common = anchor_common,
      anchor_by_subject = anchor_by_subject,
      anchor_note = .cv_info_note_fold_specific(anchor_common)
    )
  )
}


#' Internal: Identify Fold Spec
#' @keywords internal
.is_cv_folds_spec <- function(x) {
  is.list(x) && !is.null(x$folds) && is.list(x$folds)
}


#' Internal: Validate Fold Spec
#' @keywords internal
.validate_cv_folds_spec <- function(cv_folds, n_subjects) {
  if (!is.list(cv_folds) || is.null(cv_folds$folds) || !is.list(cv_folds$folds)) {
    stop("cv_folds must be a fold spec list with a $folds list", call. = FALSE)
  }

  folds <- cv_folds$folds
  if (length(folds) < 2) {
    stop("cv_folds$folds must contain >= 2 folds", call. = FALSE)
  }

  for (fold_name in names(folds)) {
    fold <- folds[[fold_name]]
    if (is.null(fold$train) || is.null(fold$test)) {
      stop(
        sprintf("cv_folds$folds[['%s']] must contain $train and $test", fold_name),
        call. = FALSE
      )
    }
    train_idx <- as.integer(fold$train)
    test_idx <- as.integer(fold$test)
    if (any(is.na(train_idx)) || any(is.na(test_idx))) {
      stop(sprintf("cv_folds fold '%s' has NA indices", fold_name), call. = FALSE)
    }
    if (any(train_idx < 1 | train_idx > n_subjects) ||
        any(test_idx < 1 | test_idx > n_subjects)) {
      stop(sprintf("cv_folds fold '%s' has out-of-range indices", fold_name), call. = FALSE)
    }
    if (length(intersect(train_idx, test_idx)) > 0) {
      stop(sprintf("cv_folds fold '%s' has overlapping train/test indices", fold_name), call. = FALSE)
    }
  }

  invisible(TRUE)
}


#' Internal: Classify Reference Kind
#' @keywords internal
.reference_kind <- function(reference) {
  if (.is_matrixish(reference)) {
    return("template")
  }

  if (is.character(reference) && length(reference) == 1) {
    if (reference %in% c("medoid", "centroid", "consensus")) {
      return("data_driven")
    }
    return("fixed_subject")
  }

  "unknown"
}

.reference_type_for_validation <- function(reference, data) {
  if (.is_matrixish(reference)) {
    return("template")
  }

  if (!is.character(reference) || length(reference) != 1L) {
    return("unknown")
  }

  if (reference %in% c("medoid", "centroid")) {
    return("subject")
  }
  if (reference %in% c("consensus")) {
    return("consensus")
  }

  if (reference %in% data@subjects) {
    return("subject")
  }

  "unknown"
}

.validate_reference_for_aligner <- function(reference, data, method, capabilities) {
  ref_types <- capabilities$reference_types %||% NULL
  if (is.null(ref_types)) return(invisible(TRUE))
  ref_types <- as.character(ref_types)

  ref_type <- .reference_type_for_validation(reference, data)

  if (identical(ref_type, "unknown")) {
    if (is.character(reference) && length(reference) == 1L) {
      stop(
        sprintf(
          "Unknown reference '%s'. Expected a subject id in the data, one of {medoid, centroid, consensus}, or a template matrix.",
          reference
        ),
        call. = FALSE
      )
    }
    stop(
      "Invalid reference specification. Expected a subject id, a supported token, or a template matrix.",
      call. = FALSE
    )
  }

  if (!ref_type %in% ref_types) {
    stop(
      sprintf(
        "Method '%s' does not support reference type '%s'. Supported reference types: %s",
        method,
        ref_type,
        paste(ref_types, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}


#' Internal: Shared Aligner Preamble
#'
#' Many aligners implement the same boilerplate: default train_idx, subset the
#' training data, and extract a data_list. This helper centralizes that logic.
#'
#' @param data AlignmentData.
#' @param train_idx Optional integer indices of training subjects.
#' @return List with train_idx, train_data, train_subjects, data_list.
#' @keywords internal
.aligner_preamble <- function(data, train_idx = NULL) {
  if (is.null(train_idx)) {
    train_idx <- seq_along(data@subjects)
  }

  train_data <- data[train_idx]
  list(
    train_idx = train_idx,
    train_data = train_data,
    train_subjects = train_data@subjects,
    data_list = get_data_list(train_data)
  )
}


#' Internal: Resolve Reference Specification
#' @keywords internal
.resolve_reference_spec <- function(data, reference, train_idx) {
  subjects <- data@subjects
  train_subjects <- subjects[train_idx]

  if (is.character(reference) && length(reference) == 1) {
    if (reference == "medoid") {
      # Select medoid from training subjects only
      ref_subj <- select_reference(data[train_idx], method = "medoid")
      return(list(
        reference = ref_subj,
        reference_spec = ref_subj
      ))
    } else if (reference == "consensus") {
      return(list(
        reference = "consensus",
        reference_spec = "consensus"
      ))
    } else if (reference == "centroid") {
      ref_subj <- select_reference(data[train_idx], method = "centroid")
      return(list(
        reference = ref_subj,
        reference_spec = ref_subj
      ))
    } else {
      # Treat as a subject id when it matches the data; otherwise pass through.
      if (reference %in% subjects && !reference %in% train_subjects) {
        warning(
          sprintf("Reference subject '%s' not in training set; potential leakage", reference),
          call. = FALSE
        )
      }
      return(list(
        reference = reference,
        reference_spec = reference
      ))
    }
  } else if (.is_matrixish(reference)) {
    # External template (supports both base matrix and sparse Matrix)
    return(list(
      reference = reference,
      reference_spec = "template"
    ))
  } else {
    stop("Invalid reference specification", call. = FALSE)
  }
}


#' Internal: Fit Transform for New Subject
#' @keywords internal
.fit_new_subject <- function(aligner, fit_result, data, subject_idx,
                             reference, fold_id = NULL) {
  subject <- data@subjects[subject_idx]
  returns <- aligner$capabilities$returns %||% "operator"

  # If aligner has an apply_fn, use it
  if (!is.null(aligner$apply_fn)) {
    apply_result <- aligner$apply_fn(
      fit_result = fit_result,
      new_data = data[subject_idx]
    )
    if (is.list(apply_result$transforms) && length(apply_result$transforms) > 0L) {
      return(apply_result$transforms[[1L]])
    }
    if (identical(returns, "embedding")) {
      z <- apply_result$aligned %||% NULL
      if (is.list(z)) {
        if (!is.null(names(z)) && subject %in% names(z)) {
          z <- z[[subject]]
        } else {
          z <- z[[1L]]
        }
      }
      if (!is.null(z)) {
        return(.new_embedding_transform(z, subject = subject))
      }
    }
    stop(
      sprintf("apply_fn for method '%s' did not return a transform for subject '%s'", aligner$name, subject),
      call. = FALSE
    )
  }

  # Otherwise, use the method's fit_fn with single subject
  # Use the reference_data (actual matrix) from fit_result, not the reference spec
  # This avoids issues when reference is a subject ID not in the single-subject data
  ref_to_use <- fit_result$reference_data
  if (is.null(ref_to_use)) {
    ref_to_use <- reference
  }

  single_fit <- .invoke_aligner_fit(
    aligner = aligner,
    data = data,
    reference = ref_to_use,
    train_idx = subject_idx,
    fit_role = "subject_cv_assessment_fit",
    fold_id = fold_id,
    subject_id = subject
  )

  if (is.list(single_fit$transforms) && subject %in% names(single_fit$transforms)) {
    return(single_fit$transforms[[subject]])
  }
  if (identical(returns, "embedding")) {
    z <- single_fit$aligned %||% NULL
    if (is.list(z) && subject %in% names(z)) {
      return(.new_embedding_transform(z[[subject]], subject = subject))
    }
  }
  stop(
    sprintf("Method '%s' did not produce a transform for subject '%s'", aligner$name, subject),
    call. = FALSE
  )
}


#' Internal: Subset AlignmentData by observation indices
#' @param data AlignmentData object.
#' @param obs_idx Integer vector of column indices to keep (shared), or a named
#'   list of per-subject integer vectors.
#' @return A new AlignmentData with subsetted matrices and updated obs_labels.
#' @keywords internal
.subset_obs <- function(data, obs_idx) {
  subjects <- data@subjects
  data_list <- get_data_list(data)
  obs_labels <- data@obs_labels
  obs_labels_by_subj <- .resolve_obs_labels_by_subject(data)

  new_data <- list()
  new_labels <- NULL

  if (is.list(obs_idx)) {
    # Per-subject observation indices
    for (subj in subjects) {
      idx <- obs_idx[[subj]]
      if (is.null(idx) || length(idx) == 0) {
        stop(sprintf("Empty observation index for subject '%s'", subj), call. = FALSE)
      }
      new_data[[subj]] <- data_list[[subj]][, idx, drop = FALSE]
    }
    if (!is.null(obs_labels_by_subj)) {
      new_labels <- lapply(subjects, function(s) obs_labels_by_subj[[s]][obs_idx[[s]]])
      names(new_labels) <- subjects
    }
  } else {
    # Shared observation indices
    for (subj in subjects) {
      new_data[[subj]] <- data_list[[subj]][, obs_idx, drop = FALSE]
    }
    if (!is.null(obs_labels) && is.atomic(obs_labels)) {
      new_labels <- obs_labels[obs_idx]
    } else if (!is.null(obs_labels_by_subj)) {
      new_labels <- lapply(subjects, function(s) obs_labels_by_subj[[s]][obs_idx])
      names(new_labels) <- subjects
    }
  }

  as_alignment_data(new_data, obs_labels = new_labels, space = data@space)
}


.reassemble_obs_cv_aligned <- function(data, subjects, aligned_test_by_fold, test_obs_by_fold) {
  has_full_coverage <- TRUE
  for (subj in subjects) {
    n_obs <- ncol(get_subject_data(data, subj))
    idx <- unlist(test_obs_by_fold[[subj]] %||% list(), use.names = FALSE)
    if (length(idx) == 0 || length(idx) != length(unique(idx)) ||
        !setequal(idx, seq_len(n_obs))) {
      has_full_coverage <- FALSE
      break
    }
  }

  aligned <- list()
  for (subj in subjects) {
    fold_parts <- aligned_test_by_fold[[subj]]
    if (is.null(fold_parts) || !length(fold_parts)) next
    if (isTRUE(has_full_coverage)) {
      n_obs <- ncol(get_subject_data(data, subj))
      n_feat <- nrow(fold_parts[[1L]])
      out <- matrix(NA_real_, n_feat, n_obs)
      for (fname in names(fold_parts)) {
        out[, test_obs_by_fold[[subj]][[fname]]] <- fold_parts[[fname]]
      }
      aligned[[subj]] <- out
    } else {
      aligned[[subj]] <- do.call(cbind, fold_parts)
    }
  }
  list(aligned = aligned, has_full_coverage = has_full_coverage)
}


#' Internal: Observation-Axis Cross-Validation
#'
#' Fits alignment using all subjects but partitions observations into
#' training and test sets within each fold.
#'
#' @keywords internal
.fit_cv_obs_folds <- function(data, aligner, reference, cv_folds,
                               compute_quality, return_aligned = TRUE,
                               return_fold_transforms = FALSE,
                               restrict_to_identified = FALSE,
                               restrict_tol = sqrt(.Machine$double.eps),
                               ...) {
  subjects <- data@subjects
  n_subjects <- length(subjects)
  folds <- cv_folds$folds
  n_folds <- length(folds)

  # Check method supports observation-axis CV
  caps <- aligner$capabilities %||% list()
  cv_axes <- caps$cv_axes %||% c("subject")
  if (!"observation" %in% cv_axes) {
    warning(sprintf(
      "Method '%s' declares cv_axes = [%s]; observation-axis CV may not be meaningful",
      aligner$name, paste(cv_axes, collapse = ", ")
    ), call. = FALSE)
  }

  # Determine if folds are per-subject or shared
  first_fold <- folds[[1]]
  per_subject_folds <- is.list(first_fold) && !is.null(names(first_fold)) &&
    is.null(first_fold$train_idx)

  aligned_test_by_fold <- list() # subject -> fold -> aligned test
  test_obs_by_fold <- list()     # subject -> fold -> test indices
  artifacts_by_fold <- list()
  transforms_by_fold <- if (isTRUE(return_fold_transforms)) list() else NULL
  reference_by_fold <- setNames(rep(NA_character_, length(folds)), names(folds))

  for (fold_name in names(folds)) {
    fold <- folds[[fold_name]]

    # Extract train/test observation indices
    if (per_subject_folds) {
      train_obs <- lapply(subjects, function(s) fold[[s]]$train_idx)
      test_obs <- lapply(subjects, function(s) fold[[s]]$test_idx)
      names(train_obs) <- names(test_obs) <- subjects
    } else {
      train_obs <- fold$train_idx
      test_obs <- fold$test_idx
    }

    # Subset data to training observations
    train_data <- .subset_obs(data, train_obs)

    # Resolve reference within this fold using only training observations.
    # This makes reference="medoid"/"centroid"/"consensus" behave correctly in obs-CV.
    ref_resolved_fold <- .resolve_reference_spec(train_data, reference, seq_len(n_subjects))
    reference_by_fold[[fold_name]] <- as.character(ref_resolved_fold$reference_spec)
    ref_for_fold <- ref_resolved_fold$reference

    # If the reference is a template matrix, subset it to the fold's training
    # observations so dimensions match and only training observations are used.
    if (.is_matrixish(reference)) {
      if (isTRUE(per_subject_folds)) {
        stop(
          "Template references are not supported for per-subject observation folds (provide shared folds or use obs_labels)",
          call. = FALSE
        )
      }
      ref_for_fold <- as.matrix(reference)[, train_obs, drop = FALSE]
    }

    # Fit alignment on training observations (all subjects)
    fit_result <- .invoke_aligner_fit(
      aligner = aligner,
      data = train_data,
      reference = ref_for_fold,
      train_idx = seq_len(n_subjects),
      fit_role = "observation_cv_fold",
      fold_id = fold_name,
      ...
    )

    if (isTRUE(restrict_to_identified)) {
      train_list <- get_data_list(train_data)
      fit_result$transforms <- lapply(subjects, function(subj) {
        canonicalize_orthogonal_operator(
          Q = fit_result$transforms[[subj]],
          x = train_list[[subj]],
          convention = "left",
          tol = restrict_tol
        )
      })
      names(fit_result$transforms) <- subjects
    }

    if (isTRUE(return_fold_transforms)) {
      transforms_by_fold[[fold_name]] <- fit_result$transforms
    }

    fold_assessment <- list()
    for (subj in subjects) {
      transform <- fit_result$transforms[[subj]]
      if (is.null(transform)) next
      test_idx <- if (per_subject_folds) test_obs[[subj]] else test_obs
      if (length(test_idx) == 0) next
      subj_test_data <- get_subject_data(data, subj)[, test_idx, drop = FALSE]
      aligned_test <- apply_transform(transform, subj_test_data)
      fold_assessment[[subj]] <- aligned_test
      if (is.null(aligned_test_by_fold[[subj]])) {
        aligned_test_by_fold[[subj]] <- list()
        test_obs_by_fold[[subj]] <- list()
      }
      aligned_test_by_fold[[subj]][[fold_name]] <- aligned_test
      test_obs_by_fold[[subj]][[fold_name]] <- test_idx
    }

    fold_model <- AlignmentModel(
      transforms = fit_result$transforms,
      reference = ref_resolved_fold$reference_spec,
      reference_data = fit_result$reference_data %||% NULL,
      method = aligner$name,
      space_from = fit_result$space_from %||% data@space,
      space_to = fit_result$space_to %||% data@space,
      params = c(list(...), list(
        restrict_to_identified = isTRUE(restrict_to_identified),
        restrict_tol = restrict_tol
      )),
      method_state = fit_result$method_state %||% list(),
      train_subjects = subjects
    )
    artifacts_by_fold[[fold_name]] <- list(
      model = fold_model,
      aligned_analysis = list(),
      aligned_assessment = fold_assessment,
      analysis_subjects = subjects,
      assessment_subjects = names(fold_assessment),
      axis = "observation",
      train_idx = train_obs,
      test_idx = test_obs,
      reference = ref_resolved_fold$reference_spec
    )
  }

  # Do full fit on all data for the model
  ref_resolved <- .resolve_reference_spec(data, reference, seq_len(n_subjects))
  fit_result_all <- .invoke_aligner_fit(
    aligner = aligner,
    data = data,
    reference = ref_resolved$reference,
    train_idx = seq_len(n_subjects),
    fit_role = "full_fit",
    ...
  )

  .validate_operator_transforms(
    transforms = fit_result_all$transforms,
    data_list = get_data_list(data),
    context = sprintf("fit_alignment(%s) [obs-cv full fit]", aligner$name)
  )

  transforms_all <- fit_result_all$transforms
  if (isTRUE(restrict_to_identified)) {
    data_list <- get_data_list(data)
    transforms_all <- lapply(subjects, function(subj) {
      canonicalize_orthogonal_operator(
        Q = transforms_all[[subj]],
        x = data_list[[subj]],
        convention = "left",
        tol = restrict_tol
      )
    })
    names(transforms_all) <- subjects
  }

  params <- c(list(...), list(
    restrict_to_identified = isTRUE(restrict_to_identified),
    restrict_tol = restrict_tol
  ))
  method_state <- fit_result_all$method_state %||% list()
  method_state$restrict_to_identified <- isTRUE(restrict_to_identified)
  method_state$restrict_tol <- restrict_tol

  model <- AlignmentModel(
    transforms = transforms_all,
    reference = ref_resolved$reference_spec,
    reference_data = fit_result_all$reference_data,
    method = aligner$name,
    space_from = fit_result_all$space_from %||% data@space,
    space_to = fit_result_all$space_to %||% data@space,
    params = params,
    method_state = method_state,
    train_subjects = subjects
  )

  reference_kind <- .reference_kind(reference)
  anchor_common <- reference_kind %in% c("fixed_subject", "template")

  reassembled <- .reassemble_obs_cv_aligned(
    data = data,
    subjects = subjects,
    aligned_test_by_fold = aligned_test_by_fold,
    test_obs_by_fold = test_obs_by_fold
  )
  evaluation_aligned <- if (isTRUE(return_aligned) && isTRUE(anchor_common)) {
    reassembled$aligned
  } else {
    list()
  }

  quality <- list()
  if (compute_quality && isTRUE(return_aligned) && isTRUE(anchor_common) &&
      length(evaluation_aligned) >= 2) {
    obs_labels_by_subject <- .resolve_obs_labels_by_subject(data)
    if (!is.null(obs_labels_by_subject) && !isTRUE(reassembled$has_full_coverage)) {
      warning(
        paste0(
          "Skipping quality computation for observation-axis CV: folds do not cover all observations, ",
          "so aligned matrices cannot be matched back to obs_labels safely."
        ),
        call. = FALSE
      )
    } else {
      quality <- .compute_basic_quality(data, evaluation_aligned, model)
    }
  }

  cv_axis <- cv_folds$axis %||% "observation"

  AlignmentResult(
    model = model,
    aligned = list(),
    quality = list(),
    cv_info = list(
      method = cv_folds$method %||% "custom",
      axis = cv_axis,
      n_folds = n_folds,
      fold_ids = cv_folds$fold_ids %||% names(folds),
      folds = folds,
      reference_kind = reference_kind,
      anchor_common = anchor_common,
      reference_by_fold = reference_by_fold,
      transforms_by_fold = transforms_by_fold,
      artifacts_by_fold = artifacts_by_fold,
      evaluation_aligned = evaluation_aligned,
      aligned_origin = "fold_assessment",
      deployment_refit = TRUE,
      artifact_kind = "observation_cv",
      model_role = "deployment_refit",
      deployment_model = model,
      anchor_note = .cv_info_note_fold_specific(anchor_common)
    )
  )
}

.fit_cv_obs_folds_embedding <- function(data, aligner, reference, cv_folds,
                                       compute_quality, return_aligned = TRUE,
                                       ...) {
  subjects <- data@subjects
  n_subjects <- length(subjects)
  folds <- cv_folds$folds
  n_folds <- length(folds)

  caps <- aligner$capabilities %||% list()
  returns <- caps$returns %||% "operator"
  if (!identical(returns, "embedding")) {
    stop(
      sprintf(".fit_cv_obs_folds_embedding called for method '%s' with returns='%s'", aligner$name, returns),
      call. = FALSE
    )
  }

  if (is.null(aligner$apply_fn)) {
    stop(
      sprintf(
        "Method '%s' returns embeddings; observation-axis CV requires an apply_fn to embed held-out observations",
        aligner$name
      ),
      call. = FALSE
    )
  }

  cv_axes <- caps$cv_axes %||% c("subject")
  if (!"observation" %in% cv_axes) {
    warning(sprintf(
      "Method '%s' declares cv_axes = [%s]; observation-axis CV may not be meaningful",
      aligner$name, paste(cv_axes, collapse = ", ")
    ), call. = FALSE)
  }

  first_fold <- folds[[1]]
  per_subject_folds <- is.list(first_fold) && !is.null(names(first_fold)) &&
    is.null(first_fold$train_idx)

  aligned_test_by_fold <- list() # subject -> fold -> aligned test
  test_obs_by_fold <- list()     # subject -> fold -> test indices
  reference_by_fold <- setNames(rep(NA_character_, length(folds)), names(folds))

  for (fold_name in names(folds)) {
    fold <- folds[[fold_name]]

    if (per_subject_folds) {
      train_obs <- lapply(subjects, function(s) fold[[s]]$train_idx)
      test_obs <- lapply(subjects, function(s) fold[[s]]$test_idx)
      names(train_obs) <- names(test_obs) <- subjects
    } else {
      train_obs <- fold$train_idx
      test_obs <- fold$test_idx
    }

    train_data <- .subset_obs(data, train_obs)

    ref_resolved_fold <- .resolve_reference_spec(train_data, reference, seq_len(n_subjects))
    reference_by_fold[[fold_name]] <- as.character(ref_resolved_fold$reference_spec)
    ref_for_fold <- ref_resolved_fold$reference

    if (.is_matrixish(reference)) {
      if (isTRUE(per_subject_folds)) {
        stop(
          "Template references are not supported for per-subject observation folds (provide shared folds or use obs_labels)",
          call. = FALSE
        )
      }
      ref_for_fold <- as.matrix(reference)[, train_obs, drop = FALSE]
    }

    fit_result <- .invoke_aligner_fit(
      aligner = aligner,
      data = train_data,
      reference = ref_for_fold,
      train_idx = seq_len(n_subjects),
      fit_role = "observation_cv_fold",
      fold_id = fold_name,
      ...
    )

    if (isTRUE(return_aligned)) {
      for (subj in subjects) {
        test_idx <- if (per_subject_folds) test_obs[[subj]] else test_obs
        if (length(test_idx) == 0) next

        subj_i <- match(subj, subjects)
        subj_test_data <- .subset_obs(data[subj_i], test_idx)

        apply_result <- aligner$apply_fn(
          fit_result = fit_result,
          new_data = subj_test_data,
          ...
        )

        z <- apply_result$aligned %||% NULL
        if (is.null(z)) {
          stop(
            sprintf(
              "apply_fn for method '%s' did not return aligned embeddings for subject '%s' in fold '%s'",
              aligner$name, subj, fold_name
            ),
            call. = FALSE
          )
        }
        if (is.list(z)) {
          if (!is.null(names(z)) && subj %in% names(z)) {
            z <- z[[subj]]
          } else {
            z <- z[[1L]]
          }
        }

        .validate_embedding_aligned(
          aligned = setNames(list(z), subj),
          data_list = setNames(get_data_list(subj_test_data), subj),
          context = sprintf("fit_alignment(%s) [obs-cv apply]", aligner$name)
        )

        if (is.null(aligned_test_by_fold[[subj]])) {
          aligned_test_by_fold[[subj]] <- list()
          test_obs_by_fold[[subj]] <- list()
        }
        aligned_test_by_fold[[subj]][[fold_name]] <- as.matrix(z)
        test_obs_by_fold[[subj]][[fold_name]] <- test_idx
      }
    }
  }

  ref_resolved <- .resolve_reference_spec(data, reference, seq_len(n_subjects))
  fit_result_all <- .invoke_aligner_fit(
    aligner = aligner,
    data = data,
    reference = ref_resolved$reference,
    train_idx = seq_len(n_subjects),
    fit_role = "full_fit",
    ...
  )

  aligned_full <- fit_result_all$aligned %||% NULL
  if (is.null(aligned_full)) {
    stop(
      sprintf("fit_alignment(%s) [obs-cv full fit]: embedding-returning methods must return 'aligned'", aligner$name),
      call. = FALSE
    )
  }

  .validate_embedding_aligned(
    aligned = aligned_full,
    data_list = get_data_list(data),
    context = sprintf("fit_alignment(%s) [obs-cv full fit]", aligner$name)
  )

  transforms_all <- lapply(subjects, function(subj) {
    z <- if (isTRUE(return_aligned)) aligned_full[[subj]] else matrix(numeric(0), 0L, 0L)
    .new_embedding_transform(z, subject = subj)
  })
  names(transforms_all) <- subjects

  model <- AlignmentModel(
    transforms = transforms_all,
    reference = ref_resolved$reference_spec,
    reference_data = fit_result_all$reference_data %||% NULL,
    method = aligner$name,
    space_from = fit_result_all$space_from %||% data@space,
    space_to = fit_result_all$space_to %||% data@space,
    params = list(...),
    method_state = fit_result_all$method_state %||% list(),
    train_subjects = subjects
  )

  reference_kind <- .reference_kind(reference)
  anchor_common <- reference_kind %in% c("fixed_subject", "template")

  reassembled <- .reassemble_obs_cv_aligned(
    data = data,
    subjects = subjects,
    aligned_test_by_fold = aligned_test_by_fold,
    test_obs_by_fold = test_obs_by_fold
  )
  evaluation_aligned <- if (isTRUE(return_aligned) && isTRUE(anchor_common)) {
    reassembled$aligned
  } else {
    list()
  }

  quality <- list()
  if (compute_quality && isTRUE(return_aligned) && isTRUE(anchor_common) &&
      length(evaluation_aligned) >= 2) {
    quality <- .compute_basic_quality(data, evaluation_aligned, model)
  }

  cv_axis <- cv_folds$axis %||% "observation"

  AlignmentResult(
    model = model,
    aligned = list(),
    quality = list(),
    cv_info = list(
      method = cv_folds$method %||% "custom",
      axis = cv_axis,
      n_folds = n_folds,
      fold_ids = cv_folds$fold_ids %||% names(folds),
      folds = folds,
      reference_kind = reference_kind,
      anchor_common = anchor_common,
      reference_by_fold = reference_by_fold,
      evaluation_aligned = evaluation_aligned,
      aligned_origin = "fold_assessment",
      deployment_refit = TRUE,
      artifact_kind = "observation_cv_embedding",
      model_role = "deployment_refit",
      deployment_model = model,
      anchor_note = .cv_info_note_fold_specific(anchor_common)
    )
  )
}


#' Internal: Apply Transforms to Data
#' @keywords internal
.apply_transforms <- function(model, data) {
  aligned <- list()

  for (subj in names(model@transforms)) {
    if (subj %in% data@subjects) {
      transform <- model@transforms[[subj]]
      subj_data <- get_subject_data(data, subj)

      # Left-multiply: Y = A %*% X
      aligned[[subj]] <- apply_transform(transform, subj_data)
    }
  }

  aligned
}


#' Internal: Compute Basic Quality Metrics
#' @keywords internal
.compute_basic_quality <- function(data, aligned, model) {
  metrics <- list()

  # Mean pairwise correlation of aligned data
  if (length(aligned) >= 2) {
    obs_labels_by_subject <- .resolve_obs_labels_by_subject(data)
    pairs <- combn(names(aligned), 2)
    cors <- apply(pairs, 2, function(p) {
      x <- as.matrix(aligned[[p[1]]])
      y <- as.matrix(aligned[[p[2]]])

      if (!is.null(obs_labels_by_subject)) {
        labs_x <- obs_labels_by_subject[[p[1]]]
        labs_y <- obs_labels_by_subject[[p[2]]]
        common <- intersect(labs_x, labs_y)
        if (length(common) < 2L) return(NA_real_)
        ix <- match(common, labs_x)
        iy <- match(common, labs_y)
        x <- x[, ix, drop = FALSE]
        y <- y[, iy, drop = FALSE]
      } else {
        if (ncol(x) != ncol(y)) return(NA_real_)
      }

      if (nrow(x) != nrow(y) || ncol(x) != ncol(y)) {
        return(NA_real_)
      }

      n_obs <- ncol(x)
      sx <- rowSums(x)
      sy <- rowSums(y)
      sxx <- rowSums(x * x)
      syy <- rowSums(y * y)
      sxy <- rowSums(x * y)

      cov_xy <- sxy - (sx * sy) / n_obs
      var_x <- sxx - (sx * sx) / n_obs
      var_y <- syy - (sy * sy) / n_obs
      denom <- sqrt(var_x * var_y)
      row_cors <- cov_xy / denom
      row_cors[denom == 0] <- NA_real_

      mean(row_cors, na.rm = TRUE)
    })
    metrics$mean_pairwise_correlation <- mean(cors, na.rm = TRUE)
    metrics$pairwise_correlations <- setNames(
      cors,
      apply(pairs, 2, paste, collapse = "-")
    )
  }

  metrics
}
