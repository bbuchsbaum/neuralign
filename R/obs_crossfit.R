#' Observation-Axis Crossfit Alignment
#'
#' Utilities for running crossfitting along the observation axis (e.g., run-wise
#' or blocked-time splits) while aligning across subjects within each fold.
#'
#' The core pattern is:
#'   1. Fit transforms on train observations only (per fold).
#'   2. Apply those transforms to held-out observations.
#'
#' This is domain-agnostic and can be used wherever a workflow can provide
#' fold-wise train/test matrices in `(features × observations)` form.
#'
#' @name obs_crossfit
NULL

.is_named_list <- function(x) {
  is.list(x) && !is.null(names(x)) && all(nzchar(names(x)))
}

.validate_folded_data <- function(x, arg_name) {
  if (!is.list(x) || length(x) < 1L) {
    stop(sprintf("'%s' must be a non-empty list of folds", arg_name), call. = FALSE)
  }
  if (is.null(names(x)) || any(!nzchar(names(x)))) {
    stop(sprintf("'%s' must be a named list keyed by fold id", arg_name), call. = FALSE)
  }

  fold_ids <- names(x)
  first <- x[[1L]]
  if (!.is_named_list(first)) {
    stop(
      sprintf("Each fold in '%s' must be a named list keyed by subject", arg_name),
      call. = FALSE
    )
  }
  subjects <- names(first)
  if (length(subjects) < 1L) {
    stop(sprintf("'%s' folds must contain >= 1 subject", arg_name), call. = FALSE)
  }

  # Validate structure and feature dimensions
  for (fid in fold_ids) {
    fold <- x[[fid]]
    if (!.is_named_list(fold)) {
      stop(
        sprintf("Fold '%s' in '%s' must be a named list keyed by subject", fid, arg_name),
        call. = FALSE
      )
    }
    missing <- setdiff(subjects, names(fold))
    extra <- setdiff(names(fold), subjects)
    if (length(missing) > 0) {
      stop(
        sprintf("Fold '%s' in '%s' is missing subjects: %s", fid, arg_name, paste(missing, collapse = ", ")),
        call. = FALSE
      )
    }
    if (length(extra) > 0) {
      stop(
        sprintf("Fold '%s' in '%s' has unexpected subjects: %s", fid, arg_name, paste(extra, collapse = ", ")),
        call. = FALSE
      )
    }

    n_feat <- NULL
    for (subj in subjects) {
      m <- fold[[subj]]
      if (!.is_matrixish(m)) {
        stop(
          sprintf("'%s' fold '%s' subject '%s' must be matrix-like", arg_name, fid, subj),
          call. = FALSE
        )
      }
      m <- as.matrix(m)
      if (is.null(n_feat)) n_feat <- nrow(m)
      if (nrow(m) != n_feat) {
        stop(
          sprintf(
            "'%s' fold '%s' has inconsistent feature dimensions across subjects",
            arg_name, fid
          ),
          call. = FALSE
        )
      }
    }
  }

  list(fold_ids = fold_ids, subjects = subjects)
}

.normalize_obs_labels_for_folds <- function(obs_labels, folded_data, label_name) {
  fold_ids <- names(folded_data)
  subjects <- names(folded_data[[1L]])

  # NULL -> all folds NULL
  if (is.null(obs_labels)) {
    out <- replicate(length(fold_ids), NULL, simplify = FALSE)
    names(out) <- fold_ids
    return(out)
  }

  # Atomic -> shared across folds/subjects
  if (is.atomic(obs_labels) || is.factor(obs_labels)) {
    v <- as.character(obs_labels)
    if (anyNA(v)) stop(sprintf("'%s' contains NA values", label_name), call. = FALSE)
    out <- lapply(fold_ids, function(fid) {
      fold <- folded_data[[fid]]
      for (subj in subjects) {
        if (length(v) != ncol(as.matrix(fold[[subj]]))) {
          stop(
            sprintf(
              "%s length mismatch for fold '%s' subject '%s': expected %d, got %d",
              label_name, fid, subj, ncol(as.matrix(fold[[subj]])), length(v)
            ),
            call. = FALSE
          )
        }
      }
      setNames(rep(list(v), length(subjects)), subjects)
    })
    names(out) <- fold_ids
    return(out)
  }

  if (!is.list(obs_labels) || is.null(names(obs_labels)) || any(!nzchar(names(obs_labels)))) {
    stop(
      sprintf("'%s' must be NULL, an atomic vector, or a named list (by fold or by subject)", label_name),
      call. = FALSE
    )
  }

  # Per-fold labels: names match fold ids
  if (all(fold_ids %in% names(obs_labels))) {
    out <- lapply(fold_ids, function(fid) {
      lab_f <- obs_labels[[fid]]
      fold <- folded_data[[fid]]

      if (is.atomic(lab_f) || is.factor(lab_f)) {
        v <- as.character(lab_f)
        if (anyNA(v)) stop(sprintf("'%s' contains NA values", label_name), call. = FALSE)
        for (subj in subjects) {
          if (length(v) != ncol(as.matrix(fold[[subj]]))) {
            stop(
              sprintf(
                "%s length mismatch for fold '%s' subject '%s': expected %d, got %d",
                label_name, fid, subj, ncol(as.matrix(fold[[subj]])), length(v)
              ),
              call. = FALSE
            )
          }
        }
        return(setNames(rep(list(v), length(subjects)), subjects))
      }

      if (!is.list(lab_f) || is.null(names(lab_f)) || any(!nzchar(names(lab_f)))) {
        stop(
          sprintf(
            "%s for fold '%s' must be an atomic vector or a named per-subject list",
            label_name, fid
          ),
          call. = FALSE
        )
      }
      missing <- setdiff(subjects, names(lab_f))
      if (length(missing) > 0) {
        stop(
          sprintf(
            "%s for fold '%s' missing subjects: %s",
            label_name, fid, paste(missing, collapse = ", ")
          ),
          call. = FALSE
        )
      }
      out_subj <- lapply(subjects, function(subj) {
        v <- as.character(lab_f[[subj]])
        if (anyNA(v)) stop(sprintf("'%s' contains NA values", label_name), call. = FALSE)
        if (length(v) != ncol(as.matrix(fold[[subj]]))) {
          stop(
            sprintf(
              "%s length mismatch for fold '%s' subject '%s': expected %d, got %d",
              label_name, fid, subj, ncol(as.matrix(fold[[subj]])), length(v)
            ),
            call. = FALSE
          )
        }
        v
      })
      names(out_subj) <- subjects
      out_subj
    })
    names(out) <- fold_ids
    return(out)
  }

  # Per-subject labels: names match subjects
  if (all(subjects %in% names(obs_labels))) {
    out <- lapply(fold_ids, function(fid) {
      fold <- folded_data[[fid]]
      out_subj <- lapply(subjects, function(subj) {
        v <- as.character(obs_labels[[subj]])
        if (anyNA(v)) stop(sprintf("'%s' contains NA values", label_name), call. = FALSE)
        if (length(v) != ncol(as.matrix(fold[[subj]]))) {
          stop(
            sprintf(
              "%s length mismatch for fold '%s' subject '%s': expected %d, got %d",
              label_name, fid, subj, ncol(as.matrix(fold[[subj]])), length(v)
            ),
            call. = FALSE
          )
        }
        v
      })
      names(out_subj) <- subjects
      out_subj
    })
    names(out) <- fold_ids
    return(out)
  }

  stop(
    sprintf(
      "'%s' names must match fold ids (%s) or subject ids (%s)",
      label_name,
      paste(fold_ids, collapse = ", "),
      paste(subjects, collapse = ", ")
    ),
    call. = FALSE
  )
}

.template_obs_labels <- function(template, template_obs_labels) {
  if (!is.null(template_obs_labels)) {
    v <- as.character(template_obs_labels)
    if (anyNA(v)) stop("'template_obs_labels' contains NA values", call. = FALSE)
    return(v)
  }
  nm <- colnames(template)
  if (!is.null(nm)) return(as.character(nm))
  NULL
}

#' Run Observation-Axis Crossfit Alignment
#'
#' Fit alignment transforms on per-fold training matrices (train observations
#' only), then apply those transforms to held-out matrices. This helper is
#' intentionally domain-agnostic: downstream packages are responsible for
#' producing fold-wise train/test matrices.
#'
#' @param train_data_by_fold Named list of folds. Each fold is a named list of
#'   subject matrices in `(features × observations)` orientation.
#' @param test_data_by_fold Optional fold list matching `train_data_by_fold`.
#' @param obs_labels_train Optional observation labels for train matrices.
#'   Accepts a shared vector, a per-subject named list, or a per-fold named list
#'   (whose entries are vectors or per-subject lists).
#' @param obs_labels_test Optional observation labels for test matrices (same
#'   forms as `obs_labels_train`).
#' @param method Alignment method passed to \code{\link{fit_alignment}}.
#' @param reference Reference specification passed to \code{\link{fit_alignment}}.
#' @param anchor_policy Policy for handling fold-specific anchors when
#'   \code{reference} is data-driven:
#'   \itemize{
#'     \item \code{"common_or_error"}: require a common anchor (fixed subject or template).
#'     \item \code{"fold_specific_ok"}: allow fold-specific anchors but mark
#'       \code{anchor_common=FALSE}.
#'     \item \code{"map_to_template"}: map fold-wise transforms into a provided
#'       template anchor (requires \code{template}).
#'   }
#' @param template Optional template matrix used when
#'   \code{anchor_policy="map_to_template"}.
#' @param template_obs_labels Optional observation labels for the template (or
#'   use \code{colnames(template)} when present).
#' @param min_overlap Minimum overlap used when mapping fold anchors to a
#'   template via \code{\link{procrustes_rotation}}.
#' @param compute_quality Logical; if TRUE, compute quality metrics for each
#'   fold model on its training data.
#' @param data_template Optional \code{\link{AlignmentData}}. When provided, its
#'   \code{space}, \code{design}, \code{geometry}, and \code{metadata} are copied
#'   into each per-fold training \code{\link{AlignmentData}} object, enabling
#'   design/geometry-dependent methods to work in observation-axis crossfit.
#' @param ... Additional arguments passed to \code{\link{fit_alignment}}.
#'
#' @return An object of class \code{"ObsCrossfitAlignment"} with fields:
#'   \itemize{
#'     \item \code{models_by_fold}: fold id -> \code{\link{AlignmentModel}}
#'     \item \code{transforms_by_fold}: fold id -> named list subject -> operator
#'     \item \code{aligned_test_by_fold}: fold id -> subject -> aligned matrix (if test data provided)
#'     \item \code{fold_info}: provenance and anchor semantics
#'   }
#'
#' @export
run_obs_crossfit_alignment <- function(train_data_by_fold,
                                      test_data_by_fold = NULL,
                                      obs_labels_train = NULL,
                                      obs_labels_test = NULL,
                                      method = "procrustes",
                                      reference = "medoid",
                                      anchor_policy = c(
                                        "common_or_error",
                                        "fold_specific_ok",
                                        "map_to_template"
                                      ),
                                      template = NULL,
                                      template_obs_labels = NULL,
                                      min_overlap = 2L,
                                      compute_quality = FALSE,
                                      data_template = NULL,
                                      ...) {
  anchor_policy <- match.arg(anchor_policy)
  min_overlap <- as.integer(min_overlap)
  if (!is.finite(min_overlap) || min_overlap < 1L) {
    stop("'min_overlap' must be a positive integer", call. = FALSE)
  }

  template_space <- NULL
  template_design <- NULL
  template_geometry <- NULL
  template_metadata <- list()
  if (!is.null(data_template)) {
    if (!inherits(data_template, "AlignmentData")) {
      stop("'data_template' must be an AlignmentData object (or NULL)", call. = FALSE)
    }
    template_space <- data_template@space
    template_design <- data_template@design
    template_geometry <- data_template@geometry
    template_metadata <- data_template@metadata %||% list()
  }

  train_meta <- .validate_folded_data(train_data_by_fold, "train_data_by_fold")
  fold_ids <- train_meta$fold_ids
  subjects <- train_meta$subjects

  if (!is.null(test_data_by_fold)) {
    test_meta <- .validate_folded_data(test_data_by_fold, "test_data_by_fold")
    if (!setequal(test_meta$fold_ids, fold_ids)) {
      stop(
        "'test_data_by_fold' must have the same fold ids as 'train_data_by_fold'",
        call. = FALSE
      )
    }
    if (!identical(test_meta$subjects, subjects)) {
      stop(
        "'test_data_by_fold' must have the same subject ids as 'train_data_by_fold'",
        call. = FALSE
      )
    }
  }

  train_labels <- .normalize_obs_labels_for_folds(
    obs_labels_train,
    folded_data = train_data_by_fold,
    label_name = "obs_labels_train"
  )
  test_labels <- if (!is.null(test_data_by_fold)) {
    .normalize_obs_labels_for_folds(
      obs_labels_test,
      folded_data = test_data_by_fold,
      label_name = "obs_labels_test"
    )
  } else {
    NULL
  }

  reference_kind <- .reference_kind(reference)
  anchor_common <- reference_kind %in% c("fixed_subject", "template")
  warnings <- character(0)

  if (anchor_policy == "common_or_error" && !anchor_common) {
    stop(
      paste0(
        "anchor_policy='common_or_error' requires a fixed subject or template reference. ",
        "Got reference kind '", reference_kind, "'."
      ),
      call. = FALSE
    )
  }

  if (anchor_policy == "fold_specific_ok" && !anchor_common) {
    warnings <- c(
      warnings,
      paste0(
        "Fold-specific anchors (reference kind '", reference_kind, "'): do not aggregate across folds ",
        "unless you map results to a common anchor."
      )
    )
  }

  if (anchor_policy == "map_to_template" && reference_kind != "template") {
    if (is.null(template) || !.is_matrixish(template)) {
      stop(
        "anchor_policy='map_to_template' requires a matrix 'template' argument",
        call. = FALSE
      )
    }
  }

  template_mat <- if (!is.null(template)) as.matrix(template) else NULL
  template_labels <- if (!is.null(template_mat)) {
    .template_obs_labels(template_mat, template_obs_labels)
  } else {
    NULL
  }

  models_by_fold <- list()
  transforms_by_fold <- list()
  aligned_test_by_fold <- if (!is.null(test_data_by_fold)) list() else NULL
  reference_by_fold <- setNames(rep(NA_character_, length(fold_ids)), fold_ids)
  quality_by_fold <- if (isTRUE(compute_quality)) list() else NULL

  for (fid in fold_ids) {
    train_fold <- train_data_by_fold[[fid]]
    train_adat <- AlignmentData(
      data = train_fold,
      obs_labels = train_labels[[fid]],
      space = template_space,
      design = template_design,
      geometry = template_geometry,
      metadata = template_metadata
    )

    # Fit fold model on training observations only.
    fit_res <- fit_alignment(
      data = train_adat,
      method = method,
      reference = reference,
      cv = "none",
      compute_quality = isTRUE(compute_quality),
      return_aligned = FALSE,
      ...
    )
    model <- get_model(fit_res)
    models_by_fold[[fid]] <- model
    transforms <- get_transforms(model)
    transforms_by_fold[[fid]] <- transforms
    reference_by_fold[[fid]] <- as.character(model@reference)

    if (isTRUE(compute_quality)) {
      quality_by_fold[[fid]] <- get_quality(fit_res)
    }

    # Optional anchor mapping: map fold transforms into template anchor.
    if (anchor_policy == "map_to_template" && reference_kind != "template") {
      ref_data <- model@reference_data
      if (is.null(ref_data) || !.is_matrixish(ref_data)) {
        stop(
          sprintf("Fold '%s' does not provide reference_data needed for template mapping", fid),
          call. = FALSE
        )
      }
      if (is.null(template_labels)) {
        stop(
          "Template mapping requires template_obs_labels or colnames(template)",
          call. = FALSE
        )
      }

      # Source labels: prefer the fold's labels for the reference subject if present.
      ref_subj <- as.character(model@reference)
      src_labels <- NULL
      if (!is.null(train_labels[[fid]]) && is.list(train_labels[[fid]])) {
        if (ref_subj %in% names(train_labels[[fid]])) {
          src_labels <- train_labels[[fid]][[ref_subj]]
        } else {
          src_labels <- train_labels[[fid]][[subjects[[1L]]]]
        }
      }

      map_res <- procrustes_rotation(
        source = ref_data,
        target = template_mat,
        convention = "left",
        obs_labels_source = src_labels,
        obs_labels_target = template_labels,
        min_overlap = min_overlap
      )
      Q_map <- map_res$Q

      transforms <- lapply(transforms, function(Q_subj) Q_map %*% as.matrix(Q_subj))
      transforms_by_fold[[fid]] <- transforms
      # Do not mutate the stored model object; record mapping in fold_info.
    }

    # Apply transforms to held-out data
    if (!is.null(test_data_by_fold)) {
      test_fold <- test_data_by_fold[[fid]]
      aligned_fold <- lapply(subjects, function(subj) {
        Q_subj <- transforms_by_fold[[fid]][[subj]]
        X <- as.matrix(test_fold[[subj]])
        Q_subj %*% X
      })
      names(aligned_fold) <- subjects
      aligned_test_by_fold[[fid]] <- aligned_fold
    }
  }

  if (anchor_policy == "map_to_template") {
    anchor_common <- TRUE
  }

  if (length(warnings) > 0) {
    warning(paste(unique(warnings), collapse = "\n"), call. = FALSE)
  }

  out <- list(
    models_by_fold = models_by_fold,
    transforms_by_fold = transforms_by_fold,
    aligned_test_by_fold = aligned_test_by_fold,
    fold_info = list(
      fold_ids = fold_ids,
      subjects = subjects,
      reference_by_fold = reference_by_fold
    ),
    reference_kind = reference_kind,
    anchor_common = anchor_common,
    provenance = list(
      method = method,
      reference = reference,
      anchor_policy = anchor_policy,
      min_overlap = min_overlap,
      params = list(...)
    ),
    quality_by_fold = quality_by_fold,
    warnings = warnings
  )
  class(out) <- "ObsCrossfitAlignment"
  out
}


#' Run Observation-Axis Crossfit from AlignmentData
#'
#' Convenience wrapper around \code{\link{run_obs_crossfit_alignment}} that
#' slices an \code{\link{AlignmentData}} object into per-fold train/test matrices
#' using a fold specification produced by \code{\link{create_obs_folds}}.
#'
#' @param data An \code{\link{AlignmentData}} object (or input coercible via
#'   \code{\link{as_alignment_data}}).
#' @param obs_folds Observation-axis fold spec as returned by
#'   \code{\link{create_obs_folds}}.
#' @param method Alignment method passed to \code{\link{fit_alignment}}.
#' @param reference Reference specification passed to \code{\link{fit_alignment}}.
#' @param anchor_policy Passed to \code{\link{run_obs_crossfit_alignment}}.
#' @param compute_quality Logical; if TRUE, compute per-fold training quality.
#' @param ... Additional arguments passed to \code{\link{fit_alignment}}.
#'
#' @return An object of class \code{"ObsCrossfitAlignment"}.
#'
#' @export
run_obs_crossfit_from_data <- function(data,
                                       obs_folds,
                                       method = "procrustes",
                                       reference = "medoid",
                                       anchor_policy = c(
                                         "common_or_error",
                                         "fold_specific_ok",
                                         "map_to_template"
                                       ),
                                       compute_quality = FALSE,
                                       ...) {
  anchor_policy <- match.arg(anchor_policy)

  if (!inherits(data, "AlignmentData")) {
    data <- as_alignment_data(data)
  }

  if (!is.list(obs_folds) || !identical(obs_folds$axis %||% NULL, "observation") ||
      !is.list(obs_folds$folds)) {
    stop("'obs_folds' must be a create_obs_folds() spec with axis='observation' and a $folds list", call. = FALSE)
  }

  folds <- obs_folds$folds
  if (length(folds) < 1L || is.null(names(folds)) || any(!nzchar(names(folds)))) {
    stop("'obs_folds$folds' must be a non-empty named list keyed by fold id", call. = FALSE)
  }

  subjects <- data@subjects
  data_list <- get_data_list(data)

  if (is.null(names(data_list)) || any(!nzchar(names(data_list))) || !setequal(names(data_list), subjects)) {
    stop("AlignmentData must contain a named data list keyed by subject id", call. = FALSE)
  }
  data_list <- data_list[subjects]

  is_matrix_subj <- vapply(data_list, .is_matrixish, logical(1))
  if (!all(is_matrix_subj)) {
    bad <- names(is_matrix_subj)[!is_matrix_subj]
    stop(
      sprintf(
        "run_obs_crossfit_from_data requires matrix-like subject data; non-matrix subjects: %s",
        paste(bad, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  fold_ids <- obs_folds$fold_ids %||% names(folds)
  if (!setequal(fold_ids, names(folds))) {
    stop("'obs_folds$fold_ids' must match names(obs_folds$folds)", call. = FALSE)
  }
  fold_ids <- fold_ids[fold_ids %in% names(folds)]

  first_fold <- folds[[fold_ids[[1L]]]]
  per_subject_folds <- is.list(first_fold) && !is.null(names(first_fold)) &&
    is.null(first_fold$train_idx)

  labels_by_subject <- .resolve_obs_labels_by_subject(data)

  check_idx <- function(idx, n_obs, context, allow_empty = FALSE) {
    idx <- as.integer(idx)
    if (anyNA(idx)) stop(sprintf("%s contains NA indices", context), call. = FALSE)
    if (!allow_empty && length(idx) < 1L) stop(sprintf("%s is empty", context), call. = FALSE)
    if (length(idx) > 0 && any(idx < 1L | idx > n_obs)) {
      stop(sprintf("%s contains out-of-range indices", context), call. = FALSE)
    }
    idx
  }

  train_data_by_fold <- list()
  test_data_by_fold <- list()
  obs_labels_train <- if (!is.null(labels_by_subject)) list() else NULL
  obs_labels_test <- if (!is.null(labels_by_subject)) list() else NULL

  for (fid in fold_ids) {
    fold <- folds[[fid]]

    train_fold <- vector("list", length(subjects))
    test_fold <- vector("list", length(subjects))
    names(train_fold) <- names(test_fold) <- subjects

    if (isTRUE(per_subject_folds)) {
      if (!is.list(fold) || is.null(names(fold))) {
        stop(sprintf("Fold '%s' must be a named list keyed by subject", fid), call. = FALSE)
      }
      missing <- setdiff(subjects, names(fold))
      if (length(missing) > 0) {
        stop(sprintf("Fold '%s' is missing subjects: %s", fid, paste(missing, collapse = ", ")), call. = FALSE)
      }

      for (subj in subjects) {
        X <- data_list[[subj]]
        n_obs <- ncol(X)
        train_idx <- check_idx(fold[[subj]]$train_idx, n_obs, sprintf("Fold '%s' subject '%s' train_idx", fid, subj))
        test_idx <- check_idx(fold[[subj]]$test_idx, n_obs, sprintf("Fold '%s' subject '%s' test_idx", fid, subj), allow_empty = TRUE)

        train_fold[[subj]] <- X[, train_idx, drop = FALSE]
        test_fold[[subj]] <- X[, test_idx, drop = FALSE]
      }

      if (!is.null(labels_by_subject)) {
        obs_labels_train[[fid]] <- setNames(lapply(subjects, function(s) labels_by_subject[[s]][fold[[s]]$train_idx]), subjects)
        obs_labels_test[[fid]] <- setNames(lapply(subjects, function(s) labels_by_subject[[s]][fold[[s]]$test_idx]), subjects)
      }
    } else {
      train_idx <- fold$train_idx %||% NULL
      test_idx <- fold$test_idx %||% NULL
      if (is.null(train_idx) || is.null(test_idx)) {
        stop(sprintf("Fold '%s' must contain $train_idx and $test_idx", fid), call. = FALSE)
      }

      for (subj in subjects) {
        X <- data_list[[subj]]
        n_obs <- ncol(X)
        tr <- check_idx(train_idx, n_obs, sprintf("Fold '%s' train_idx", fid))
        te <- check_idx(test_idx, n_obs, sprintf("Fold '%s' test_idx", fid), allow_empty = TRUE)
        train_fold[[subj]] <- X[, tr, drop = FALSE]
        test_fold[[subj]] <- X[, te, drop = FALSE]
      }

      if (!is.null(labels_by_subject)) {
        obs_labels_train[[fid]] <- setNames(lapply(subjects, function(s) labels_by_subject[[s]][train_idx]), subjects)
        obs_labels_test[[fid]] <- setNames(lapply(subjects, function(s) labels_by_subject[[s]][test_idx]), subjects)
      }
    }

    train_data_by_fold[[fid]] <- train_fold
    test_data_by_fold[[fid]] <- test_fold
  }

  res <- run_obs_crossfit_alignment(
    train_data_by_fold = train_data_by_fold,
    test_data_by_fold = test_data_by_fold,
    obs_labels_train = obs_labels_train,
    obs_labels_test = obs_labels_test,
    method = method,
    reference = reference,
    anchor_policy = anchor_policy,
    compute_quality = compute_quality,
    data_template = data,
    ...
  )

  res$fold_info$obs_folds <- list(
    method = obs_folds$method %||% "custom",
    axis = obs_folds$axis %||% "observation",
    fold_ids = fold_ids,
    n_folds = obs_folds$n_folds %||% length(fold_ids),
    guard_tr = obs_folds$guard_tr %||% NA_integer_,
    id_policy = obs_folds$id_policy %||% NA_character_,
    common_ids = obs_folds$common_ids %||% NA
  )

  res
}
