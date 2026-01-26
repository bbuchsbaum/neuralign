#' Fit an Alignment Model
#'
#' Fit an alignment model to multi-subject data. This is the main entry point
#' for alignment in neuralign. The fitted model can be saved and applied to
#' new data using \code{\link{apply_alignment}}.
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
                          ...) {
  cv <- match.arg(cv)

  # Coerce data to AlignmentData if needed
  if (inherits(data, "AlignmentData")) {
    if (!is.null(obs_labels)) {
      if (!is.null(data@obs_labels) &&
          !identical(as.character(data@obs_labels), as.character(obs_labels))) {
        stop("obs_labels supplied but AlignmentData already has different obs_labels", call. = FALSE)
      }
      if (is.null(data@obs_labels)) {
        data@obs_labels <- obs_labels
      }
    }
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

  # Validate data dimensions based on method capabilities
  validate_alignment_data(
    data,
    check_features = isTRUE(caps$requires_shared_features %||% TRUE),
    check_observations = isTRUE(caps$requires_shared_observations %||% FALSE)
  )

  # Route based on CV strategy
  if (cv == "none") {
    result <- .fit_single(data, aligner, reference, train_idx,
      compute_quality = compute_quality, ...
    )
  } else if (.is_cv_folds_spec(cv_folds)) {
    validate_cv_setup(cv_folds, reference = reference)
    result <- .fit_cv_folds(
      data = data,
      aligner = aligner,
      reference = reference,
      cv_folds = cv_folds,
      compute_quality = compute_quality,
      ...
    )
  } else if (cv == "loso") {
    result <- .fit_cv_loso(
      data = data,
      aligner = aligner,
      reference = reference,
      compute_quality = compute_quality,
      ...
    )
  } else if (cv == "kfold") {
    result <- .fit_cv_kfold(
      data = data,
      aligner = aligner,
      reference = reference,
      k = cv_folds,
      compute_quality = compute_quality,
      ...
    )
  }

  result
}


#' Internal: Single Fit (No CV)
#' @keywords internal
.fit_single <- function(data, aligner, reference, train_idx,
                        compute_quality, ...) {
  # Determine training subjects
  if (is.null(train_idx)) {
    train_idx <- seq_along(data@subjects)
  }
  train_subjects <- data@subjects[train_idx]

  # Resolve reference
  ref_resolved <- .resolve_reference(data, reference, train_idx)

  # Call the fit function
  fit_result <- aligner$fit_fn(
    data = data,
    reference = ref_resolved$reference,
    train_idx = train_idx,
    ...
  )

  .validate_operator_transforms(
    transforms = fit_result$transforms,
    data_list = get_data_list(data),
    context = sprintf("fit_alignment(%s)", aligner$name)
  )

  # Build AlignmentModel
  model <- AlignmentModel(
    transforms = fit_result$transforms,
    reference = ref_resolved$reference_spec,
    reference_data = fit_result$reference_data,
    method = aligner$name,
    space_from = fit_result$space_from,
    space_to = fit_result$space_to,
    params = list(...),
    method_state = fit_result$method_state %||% list(),
    train_subjects = train_subjects
  )

  # Apply alignment to get aligned data
  aligned <- .apply_transforms(model, data)

  # Compute quality metrics if requested
  quality <- list()
  if (compute_quality) {
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
.fit_cv_loso <- function(data, aligner, reference, compute_quality, ...) {
  cv_folds <- create_cv_folds(data, method = "loso")
  validate_cv_setup(cv_folds, reference = reference)
  .fit_cv_folds(data, aligner, reference, cv_folds,
    compute_quality = compute_quality, ...
  )
}


#' Internal: K-Fold CV
#' @keywords internal
.fit_cv_kfold <- function(data, aligner, reference, k, compute_quality, ...) {
  cv_folds <- create_cv_folds(data, method = "kfold", k = k)
  # create_cv_folds will fall back to LOSO if k > n.
  validate_cv_setup(cv_folds, reference = reference)
  .fit_cv_folds(data, aligner, reference, cv_folds,
    compute_quality = compute_quality, ...
  )
}


#' Internal: Generic CV with Provided Folds
#' @keywords internal
.fit_cv_folds <- function(data, aligner, reference, cv_folds,
                          compute_quality, ...) {
  .validate_cv_folds_spec(cv_folds, n_subjects = length(data@subjects))

  subjects <- data@subjects
  n_subjects <- length(subjects)

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
    ref_resolved <- .resolve_reference(data, reference, train_idx)

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

      test_data <- get_subject_data(data, test_subj)
      all_aligned[[test_subj]] <- test_transform %*% test_data
      anchor_by_subject[[test_subj]] <- as.character(ref_resolved$reference_spec)
    }
  }

  missing_subjects <- setdiff(subjects, names(all_transforms))
  if (length(missing_subjects) > 0) {
    stop(
      sprintf(
        "CV fold spec never assigns these subjects to a test fold: %s",
        paste(missing_subjects, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  reference_kind <- .reference_kind(reference)
  anchor_common <- reference_kind %in% c("fixed_subject", "template")

  model_reference <- reference
  model_reference_data <- NULL
  model_method_state <- list()
  model_space_from <- data@space
  model_space_to <- data@space

  if (anchor_common) {
    # For a fixed/external anchor, it's safe to store the corresponding reference_data.
    ref_resolved <- .resolve_reference(data, reference, seq_len(n_subjects))
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

    model_reference <- ref_resolved$reference_spec
    model_reference_data <- fit_result_all$reference_data
    model_method_state <- fit_result_all$method_state %||% list()
    model_space_from <- fit_result_all$space_from
    model_space_to <- fit_result_all$space_to
  } else {
    # Fold-specific anchors: do not pretend there is a single shared reference.
    model_reference <- "fold_specific"
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

  # Quality metrics
  quality <- list()
  if (compute_quality) {
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
      anchor_note = if (!anchor_common) {
        "Aligned outputs are in fold-specific anchor spaces; do not use for group-level comparisons without mapping to a common anchor."
      } else {
        NULL
      }
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
  if (is.matrix(reference)) {
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


#' Internal: Resolve Reference Specification
#' @keywords internal
.resolve_reference <- function(data, reference, train_idx) {
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
      # Assume it's a subject ID
      if (!reference %in% train_subjects) {
        warning(sprintf(
          "Reference subject '%s' not in training set; potential leakage",
          reference
        ))
      }
      return(list(
        reference = reference,
        reference_spec = reference
      ))
    }
  } else if (is.matrix(reference)) {
    # External template
    return(list(
      reference = reference,
      reference_spec = "template"
    ))
  } else {
    stop("Invalid reference specification")
  }
}


#' Internal: Fit Transform for New Subject
#' @keywords internal
.fit_new_subject <- function(aligner, fit_result, data, subject_idx,
                             reference) {
  subject <- data@subjects[subject_idx]

  # If aligner has an apply_fn, use it
  if (!is.null(aligner$apply_fn)) {
    apply_result <- aligner$apply_fn(
      fit_result = fit_result,
      new_data = data[subject_idx]
    )
    return(apply_result$transforms[[1]])
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

  single_fit$transforms[[subject]]
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
      aligned[[subj]] <- transform %*% subj_data
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
