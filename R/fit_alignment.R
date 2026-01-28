#' Fit an Alignment Model
#'
#' Fit an alignment model to multi-subject data. This is the main entry point
#' for alignment in neuralign. The fitted model can be saved and applied to
#' new data using \code{\link{apply_alignment}}.
#'
#' @details
#' \code{fit_alignment()} does not perform any centering or scaling of inputs by
#' default. If you want preprocessing, apply it explicitly to your matrices or
#' feature blocks (e.g., \code{\link{preprocess_alignment_data}} or
#' \code{\link{preprocess_matrix}}) before fitting.
#'
#' @param data AlignmentData object, or a named list of matrices that will
#'   be coerced to AlignmentData.
#' @param method Character string specifying the alignment method. Use
#'   \code{\link{available_aligners}} to see registered methods.
#' @param reference How to select the reference for alignment:
#'   \itemize{
#'     \item "medoid" - Select the medoid subject (default)
#'     \item "consensus" - Compute a consensus/mean reference
#'     \item Subject ID - Use a specific subject as reference
#'     \item Matrix - Use an external template matrix
#'   }
#' @param cv Cross-validation strategy:
#'   \itemize{
#'     \item "none" - No cross-validation (default)
#'     \item "loso" - Leave-one-subject-out
#'     \item "kfold" - K-fold cross-validation (specify k via cv_folds)
#'   }
#' @param cv_folds Either an integer number of folds for k-fold CV (used when
#'   \code{cv = "kfold"}), or a fold specification list from
#'   \code{\link{create_cv_folds}} (or compatible) to control the exact train/test
#'   splits.
#' @param train_idx Optional integer indices specifying which subjects to use
#'   for fitting. If NULL, all subjects are used. This enables manual CV schemes.
#' @param obs_labels Optional shared observation labels. If \code{data} is a
#'   list, this is passed through to \code{\link{AlignmentData}} via
#'   \code{\link{as_alignment_data}}. If \code{data} is already an
#'   \code{AlignmentData} and has no observation labels, they will be set.
#' @param compute_quality Logical; if TRUE, compute quality metrics.
#' @param return_aligned Logical; if TRUE (default), apply transforms and store
#'   aligned data in the result. Set to FALSE to return only the model (saves
#'   memory in pipeline workflows where aligned data is not needed immediately).
#' @param return_fold_transforms Logical; if TRUE and using observation-axis
#'   fold specs (see \code{\link{create_obs_folds}}), retain per-fold transforms
#'   in \code{result@cv_info$transforms_by_fold}. Default FALSE.
#' @param restrict_to_identified Logical; if TRUE, canonicalize orthogonal
#'   operators so their action on the orthogonal complement of the identified
#'   subspace is deterministic. This is useful when correspondence signals are
#'   rank-deficient (e.g., transform_dim > numeric_rank). Defaults to FALSE.
#'   Only supported for methods that return orthogonal operators.
#' @param restrict_tol Relative tolerance for numeric rank determination used by
#'   \code{\link{canonicalize_orthogonal_operator}}.
#' @param ... Additional arguments passed to the method's fit function.
#'
#' @return An AlignmentResult object containing the fitted model and
#'   (optionally) aligned data and quality metrics.
#'
#' @examples
#' \dontrun{
#' # Create example data
#' data_list <- list(
#'   "sub-01" = matrix(rnorm(100*50), 100, 50),
#'   "sub-02" = matrix(rnorm(100*50), 100, 50),
#'   "sub-03" = matrix(rnorm(100*50), 100, 50)
#' )
#' adat <- AlignmentData(data_list)
#'
#' # Fit alignment with default settings
#' result <- fit_alignment(adat, method = "procrustes")
#'
#' # Fit with specific reference subject
#' result <- fit_alignment(adat, method = "procrustes", reference = "sub-01")
#'
#' # Fit with leave-one-subject-out CV
#' result <- fit_alignment(adat, method = "procrustes", cv = "loso")
#' }
#'
#' @seealso \code{\link{apply_alignment}}, \code{\link{AlignmentData}},
#'   \code{\link{available_aligners}}
#'
#' @export
fit_alignment <- function(data,
                          method = "procrustes",
                          reference = "medoid",
                          cv = c("none", "loso", "kfold"),
                          cv_folds = 5,
                          train_idx = NULL,
                          obs_labels = NULL,
                          compute_quality = TRUE,
                          return_aligned = TRUE,
                          return_fold_transforms = FALSE,
                          restrict_to_identified = FALSE,
                          restrict_tol = sqrt(.Machine$double.eps),
                          ...) {
  cv <- match.arg(cv)

  # Coerce data to AlignmentData if needed
  if (inherits(data, "AlignmentData")) {
    data <- .set_obs_labels_if_missing(data, obs_labels)
  } else {
    data <- as_alignment_data(data, obs_labels = obs_labels)
  }

  # Try to load aligner if not registered
  if (!is_aligner_registered(method)) {
    if (!.try_autoload_aligner(method)) {
      stop(sprintf(
        "Unknown method '%s'. Available methods: %s",
        method,
        paste(available_aligners(), collapse = ", ")
      ))
    }
  }

  # Validate method requirements
  .validate_aligner_requirements(method, data)

  # Get aligner
  aligner <- get_aligner(method)
  caps <- aligner$capabilities %||% list()
  returns <- caps$returns %||% "operator"

  restrict_to_identified <- isTRUE(restrict_to_identified)
  if (restrict_to_identified) {
    if (!identical(returns, "operator")) {
      stop(
        sprintf(
          "restrict_to_identified=TRUE requires operator-returning methods; method '%s' returns '%s'",
          method, returns
        ),
        call. = FALSE
      )
    }
    transform_type <- caps$transform_type %||% NA_character_
    if (!is.character(transform_type) || length(transform_type) != 1L ||
        !transform_type %in% c("orthogonal")) {
      stop(
        sprintf(
          "restrict_to_identified=TRUE is only supported for orthogonal operators; method '%s' has transform_type='%s'",
          method, as.character(transform_type)
        ),
        call. = FALSE
      )
    }
    if (!is.numeric(restrict_tol) || length(restrict_tol) != 1L || !is.finite(restrict_tol) ||
        restrict_tol <= 0 || restrict_tol >= 1) {
      stop("'restrict_tol' must be a single number in (0, 1)", call. = FALSE)
    }
  }

  if (identical(returns, "embedding")) {
    if (!identical(cv, "none") || .is_cv_folds_spec(cv_folds)) {
      stop(
        sprintf(
          "Method '%s' returns embeddings; cross-validation is not yet supported for embedding-returning methods. Use cv='none'.",
          method
        ),
        call. = FALSE
      )
    }
  }

  .validate_reference_for_aligner(reference, data, method, caps)

  # Validate data dimensions based on method capabilities
  validate_alignment_data(
    data,
    check_features = isTRUE(caps$requires_shared_features %||% TRUE),
    check_observations = isTRUE(caps$requires_shared_observations %||% FALSE)
  )

  # Route based on CV strategy
  # Check for explicit fold spec first (overrides cv="none" default)
  if (.is_cv_folds_spec(cv_folds)) {
    if (identical(cv_folds$axis, "observation")) {
      result <- .fit_cv_obs_folds(
        data = data,
        aligner = aligner,
        reference = reference,
        cv_folds = cv_folds,
        compute_quality = compute_quality,
        return_aligned = return_aligned,
        return_fold_transforms = return_fold_transforms,
        restrict_to_identified = restrict_to_identified,
        restrict_tol = restrict_tol,
        ...
      )
    } else {
      validate_cv_setup(cv_folds, reference = reference)
      result <- .fit_cv_folds(
        data = data,
        aligner = aligner,
        reference = reference,
        cv_folds = cv_folds,
        compute_quality = compute_quality,
        return_aligned = return_aligned,
        restrict_to_identified = restrict_to_identified,
        restrict_tol = restrict_tol,
        ...
      )
    }
  } else if (cv == "none") {
    result <- .fit_single(data, aligner, reference, train_idx,
      compute_quality = compute_quality,
      return_aligned = return_aligned,
      restrict_to_identified = restrict_to_identified,
      restrict_tol = restrict_tol,
      ...
    )
  } else if (cv == "loso") {
    result <- .fit_cv_loso(
      data = data,
      aligner = aligner,
      reference = reference,
      compute_quality = compute_quality,
      return_aligned = return_aligned,
      restrict_to_identified = restrict_to_identified,
      restrict_tol = restrict_tol,
      ...
    )
  } else if (cv == "kfold") {
    result <- .fit_cv_kfold(
      data = data,
      aligner = aligner,
      reference = reference,
      k = cv_folds,
      compute_quality = compute_quality,
      return_aligned = return_aligned,
      restrict_to_identified = restrict_to_identified,
      restrict_tol = restrict_tol,
      ...
    )
  }

  result
}


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
  fit_result <- aligner$fit_fn(
    data = data,
    reference = ref_resolved$reference,
    train_idx = train_idx,
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
                         restrict_to_identified = FALSE,
                         restrict_tol = sqrt(.Machine$double.eps),
                         ...) {
  cv_folds <- create_cv_folds(data, method = "loso")
  validate_cv_setup(cv_folds, reference = reference)
  .fit_cv_folds(data, aligner, reference, cv_folds,
    compute_quality = compute_quality,
    return_aligned = return_aligned,
    restrict_to_identified = restrict_to_identified,
    restrict_tol = restrict_tol,
    ...
  )
}


#' Internal: K-Fold CV
#' @keywords internal
.fit_cv_kfold <- function(data, aligner, reference, k, compute_quality,
                          return_aligned = TRUE,
                          restrict_to_identified = FALSE,
                          restrict_tol = sqrt(.Machine$double.eps),
                          ...) {
  cv_folds <- create_cv_folds(data, method = "kfold", k = k)
  # create_cv_folds will fall back to LOSO if k > n.
  validate_cv_setup(cv_folds, reference = reference)
  .fit_cv_folds(data, aligner, reference, cv_folds,
    compute_quality = compute_quality,
    return_aligned = return_aligned,
    restrict_to_identified = restrict_to_identified,
    restrict_tol = restrict_tol,
    ...
  )
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
  fit_result_all <- aligner$fit_fn(
    data = data,
    reference = ref_resolved$reference,
    train_idx = seq_len(n_subjects),
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

.cv_info_note_fold_specific <- function(anchor_common) {
  if (isTRUE(anchor_common)) return(NULL)
  "Aligned outputs are in fold-specific anchor spaces; do not use for group-level comparisons without mapping to a common anchor."
}

.fit_cv_folds <- function(data, aligner, reference, cv_folds,
                          compute_quality, return_aligned = TRUE,
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

  all_transforms <- list()
  all_aligned <- list()
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
    fit_result <- aligner$fit_fn(
      data = data,
      reference = ref_resolved$reference,
      train_idx = train_idx,
      ...
    )

    .validate_operator_transforms(
      transforms = fit_result$transforms,
      data_list = get_data_list(data)[subjects[train_idx]],
      context = sprintf("fit_alignment(%s) [cv train]", aligner$name)
    )

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
        aligner, fit_result, data, test_i, ref_resolved$reference
      )
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
  if (length(missing_subjects) > 0) {
    for (subj in missing_subjects) {
      X <- get_subject_data(data, subj)
      all_transforms[[subj]] <- diag(nrow(X))
      if (return_aligned) {
        all_aligned[[subj]] <- X
      }
      anchor_by_subject[[subj]] <- as.character(reference)
    }
  }

  if (isTRUE(restrict_to_identified)) {
    data_list <- get_data_list(data)
    all_transforms <- lapply(subjects, function(subj) {
      canonicalize_orthogonal_operator(
        Q = all_transforms[[subj]],
        x = data_list[[subj]],
        convention = "left",
        tol = restrict_tol
      )
    })
    names(all_transforms) <- subjects
  }

  reference_kind <- .reference_kind(reference)
  anchor_common <- reference_kind %in% c("fixed_subject", "template")

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
  if (compute_quality && return_aligned) {
    quality <- .compute_basic_quality(data, all_aligned, model)
  }

  cv_method <- cv_folds$method %||% "custom"
  cv_axis <- cv_folds$axis %||% "subject"

  AlignmentResult(
    model = model,
    aligned = all_aligned,
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
                             reference) {
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

  single_fit <- aligner$fit_fn(
    data = data,
    reference = ref_to_use,
    train_idx = subject_idx
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
    fit_result <- aligner$fit_fn(
      data = train_data,
      reference = ref_for_fold,
      train_idx = seq_len(n_subjects),
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

    # Apply transforms to held-out observations
    if (return_aligned) {
      for (subj in subjects) {
        transform <- fit_result$transforms[[subj]]
        if (is.null(transform)) next
        test_idx <- if (per_subject_folds) test_obs[[subj]] else test_obs
        if (length(test_idx) == 0) next
        subj_test_data <- get_subject_data(data, subj)[, test_idx, drop = FALSE]
        aligned_test <- apply_transform(transform, subj_test_data)
        if (is.null(aligned_test_by_fold[[subj]])) {
          aligned_test_by_fold[[subj]] <- list()
          test_obs_by_fold[[subj]] <- list()
        }
        aligned_test_by_fold[[subj]][[fold_name]] <- aligned_test
        test_obs_by_fold[[subj]][[fold_name]] <- test_idx
      }
    }
  }

  # Do full fit on all data for the model
  ref_resolved <- .resolve_reference_spec(data, reference, seq_len(n_subjects))
  fit_result_all <- aligner$fit_fn(
    data = data,
    reference = ref_resolved$reference,
    train_idx = seq_len(n_subjects),
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

  if (!anchor_common) {
    warning(
      paste0(
        "Observation-axis CV with reference kind '", reference_kind, "' uses fold-specific anchors. ",
        "Aligned outputs across folds may not be directly comparable unless you map them to a common anchor."
      ),
      call. = FALSE
    )
  }

  # Reassemble aligned test data. If folds cover all observations (per subject),
  # return a full features×observations matrix in original order; otherwise fall
  # back to concatenating held-out segments.
  aligned <- list()
  has_full_coverage <- FALSE
  if (return_aligned) {
    # Determine whether test folds cover every observation exactly once.
    has_full_coverage <- TRUE
    for (subj in subjects) {
      n_obs <- ncol(get_subject_data(data, subj))
      idx <- unlist(test_obs_by_fold[[subj]] %||% list(), use.names = FALSE)
      if (length(idx) == 0) {
        has_full_coverage <- FALSE
        break
      }
      if (length(idx) != length(unique(idx)) || !setequal(idx, seq_len(n_obs))) {
        has_full_coverage <- FALSE
        break
      }
    }

    for (subj in subjects) {
      fold_parts <- aligned_test_by_fold[[subj]]
      if (is.null(fold_parts) || length(fold_parts) == 0) next

      if (has_full_coverage) {
        n_obs <- ncol(get_subject_data(data, subj))
        n_feat <- nrow(fold_parts[[1L]])
        out <- matrix(NA_real_, n_feat, n_obs)
        for (fname in names(fold_parts)) {
          out[, test_obs_by_fold[[subj]][[fname]]] <- fold_parts[[fname]]
        }
        aligned[[subj]] <- out
      } else {
        # Partial coverage: return only held-out segments (fold order).
        aligned[[subj]] <- do.call(cbind, fold_parts)
      }
    }
  }

  quality <- list()
  if (compute_quality && return_aligned && length(aligned) >= 2) {
    obs_labels_by_subject <- .resolve_obs_labels_by_subject(data)
    if (!is.null(obs_labels_by_subject) && !isTRUE(has_full_coverage)) {
      warning(
        paste0(
          "Skipping quality computation for observation-axis CV: folds do not cover all observations, ",
          "so aligned matrices cannot be matched back to obs_labels safely."
        ),
        call. = FALSE
      )
    } else {
      quality <- .compute_basic_quality(data, aligned, model)
    }
  }

  cv_axis <- cv_folds$axis %||% "observation"

  AlignmentResult(
    model = model,
    aligned = aligned,
    quality = quality,
    cv_info = list(
      method = cv_folds$method %||% "custom",
      axis = cv_axis,
      n_folds = n_folds,
      fold_ids = cv_folds$fold_ids %||% names(folds),
      folds = folds,
      reference_kind = reference_kind,
      anchor_common = anchor_common,
      reference_by_fold = reference_by_fold,
      transforms_by_fold = transforms_by_fold
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
