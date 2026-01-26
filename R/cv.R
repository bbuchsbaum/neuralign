#' Cross-Validation Utilities for Alignment
#'
#' Functions for setting up and managing cross-validation in alignment workflows.
#'
#' @name cv
NULL


#' Create CV Fold Assignments
#'
#' Create fold assignments for cross-validation. This is useful for
#' implementing custom CV schemes or for reproducible experiments.
#'
#' @param data AlignmentData object, or just the number of subjects.
#' @param method CV method:
#'   \itemize{
#'     \item "loso" - Leave-one-subject-out
#'     \item "kfold" - K-fold cross-validation
#'     \item "stratified" - Stratified k-fold (requires group info)
#'   }
#' @param k Number of folds for k-fold CV.
#' @param groups Optional factor for stratified CV.
#' @param seed Random seed for reproducibility.
#'
#' @return List with fold information:
#'   \itemize{
#'     \item assignments - Named vector of fold numbers
#'     \item folds - List of train/test indices for each fold
#'     \item n_folds - Number of folds
#'   }
#'
#' @examples
#' \dontrun{
#' data <- AlignmentData(list(
#'   "sub-01" = matrix(1, 10, 5),
#'   "sub-02" = matrix(1, 10, 5),
#'   "sub-03" = matrix(1, 10, 5),
#'   "sub-04" = matrix(1, 10, 5)
#' ))
#'
#' # Leave-one-subject-out
#' cv_loso <- create_cv_folds(data, method = "loso")
#'
#' # 2-fold CV
#' cv_kfold <- create_cv_folds(data, method = "kfold", k = 2)
#' }
#'
#' @export
create_cv_folds <- function(data,
                            method = c("loso", "kfold", "stratified"),
                            k = 5,
                            groups = NULL,
                            seed = NULL) {
  method <- match.arg(method)

  # Get subjects
  if (inherits(data, "AlignmentData")) {
    subjects <- data@subjects
  } else if (is.numeric(data) && length(data) == 1) {
    subjects <- paste0("sub-", sprintf("%02d", seq_len(data)))
  } else if (is.character(data)) {
    subjects <- data
  } else {
    stop("'data' must be AlignmentData, number of subjects, or subject vector")
  }

  n <- length(subjects)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (method == "loso") {
    # Leave-one-subject-out
    assignments <- setNames(seq_len(n), subjects)
    folds <- lapply(seq_len(n), function(i) {
      list(
        train = setdiff(seq_len(n), i),
        test = i
      )
    })
    names(folds) <- subjects

    return(list(
      method = "loso",
      axis = "subject",
      assignments = assignments,
      folds = folds,
      n_folds = n
    ))
  }

  if (method == "kfold") {
    if (k > n) {
      warning(sprintf("k=%d > n=%d; using LOSO instead", k, n))
      return(create_cv_folds(data, method = "loso"))
    }

    # Random fold assignment
    assignments <- sample(rep(seq_len(k), length.out = n))
    names(assignments) <- subjects

    folds <- lapply(seq_len(k), function(fold) {
      list(
        train = which(assignments != fold),
        test = which(assignments == fold)
      )
    })
    names(folds) <- paste0("fold-", seq_len(k))

    return(list(
      method = "kfold",
      axis = "subject",
      assignments = assignments,
      folds = folds,
      n_folds = k
    ))
  }

  if (method == "stratified") {
    if (is.null(groups)) {
      stop("'groups' required for stratified CV")
    }

    if (length(groups) != n) {
      stop("'groups' must have same length as number of subjects")
    }

    groups <- as.factor(groups)

    # Stratified sampling within each group
    assignments <- integer(n)
    names(assignments) <- subjects

    for (g in levels(groups)) {
      g_idx <- which(groups == g)
      n_g <- length(g_idx)
      g_folds <- sample(rep(seq_len(k), length.out = n_g))
      assignments[g_idx] <- g_folds
    }

    folds <- lapply(seq_len(k), function(fold) {
      list(
        train = which(assignments != fold),
        test = which(assignments == fold)
      )
    })
    names(folds) <- paste0("fold-", seq_len(k))

    return(list(
      method = "stratified",
      axis = "subject",
      assignments = assignments,
      folds = folds,
      n_folds = k,
      groups = groups
    ))
  }
}

.as_obs_ids_list <- function(obs_ids) {
  if (is.atomic(obs_ids)) {
    return(NULL)
  }
  if (!is.list(obs_ids)) {
    stop("'obs_ids' must be an atomic vector or a list", call. = FALSE)
  }
  if (is.null(names(obs_ids))) {
    stop("When 'obs_ids' is a list, it must be a named list (one entry per subject)", call. = FALSE)
  }
  lapply(obs_ids, as.character)
}

.blocked_time_folds <- function(n_obs, k, guard_tr) {
  if (!is.numeric(n_obs) || length(n_obs) != 1L || n_obs < 2) {
    stop("blocked_time folds require n_obs >= 2", call. = FALSE)
  }
  k <- as.integer(k)
  guard_tr <- as.integer(guard_tr)
  if (k < 2L) stop("'k' must be >= 2 for blocked_time folds", call. = FALSE)
  if (guard_tr < 0L) stop("'guard_tr' must be >= 0", call. = FALSE)
  if (k > n_obs) stop("'k' cannot exceed number of observations", call. = FALSE)

  # Split indices into k contiguous blocks as evenly as possible
  sizes <- rep(floor(n_obs / k), k)
  sizes[seq_len(n_obs %% k)] <- sizes[seq_len(n_obs %% k)] + 1L
  ends <- cumsum(sizes)
  starts <- c(1L, utils::head(ends, -1L) + 1L)

  folds <- vector("list", k)
  for (i in seq_len(k)) {
    test_idx <- seq.int(starts[i], ends[i])
    train_idx <- setdiff(seq_len(n_obs), test_idx)
    if (guard_tr > 0L) {
      guard_lo <- max(1L, starts[i] - guard_tr)
      guard_hi <- min(n_obs, ends[i] + guard_tr)
      guard_idx <- seq.int(guard_lo, guard_hi)
      train_idx <- setdiff(train_idx, guard_idx)
    }
    if (length(train_idx) < 1L) {
      stop("guard_tr too large: training set becomes empty", call. = FALSE)
    }
    folds[[i]] <- list(train_idx = train_idx, test_idx = test_idx)
  }
  names(folds) <- paste0("block-", seq_len(k))
  folds
}

#' Create Observation-Axis Fold Specifications
#'
#' Construct cross-validation folds along the observation axis (e.g., runs or
#' contiguous time blocks). Unlike subject-axis CV, these folds return per-subject
#' index sets suitable for timepoint/run-based cross-fitting workflows.
#'
#' @param obs_ids Observation identifiers. Either a single vector (shared across
#'   subjects), or a named list of per-subject vectors.
#' @param method Fold construction method:
#'   \itemize{
#'     \item `"run"`: leave-one-run-out (each unique run id defines a fold)
#'     \item `"blocked_time"`: split a single time series into contiguous blocks
#'   }
#' @param k Number of blocks for `method="blocked_time"`.
#' @param guard_tr Number of TRs to exclude around each held-out block/run from
#'   the training set (to reduce temporal leakage).
#' @param seed Optional random seed (currently only used for future extensions).
#' @param id_policy For per-subject `obs_ids` lists and `method="run"`, choose how
#'   fold ids are defined: `"intersection"` uses only run ids present for all
#'   subjects; `"union"` uses the union (subjects without a given run contribute
#'   an empty test set for that fold).
#'
#' @return A list with elements `method`, `axis="observation"`, `folds`,
#'   `n_folds`, and `fold_ids`. When `obs_ids` is a list, each fold contains a
#'   named list of per-subject `train_idx`/`test_idx` vectors.
#'
#' @export
create_obs_folds <- function(obs_ids,
                             method = c("run", "blocked_time"),
                             k = 5,
                             guard_tr = 0L,
                             seed = NULL,
                             id_policy = c("intersection", "union")) {
  method <- match.arg(method)
  id_policy <- match.arg(id_policy)
  guard_tr <- as.integer(guard_tr)
  if (!is.finite(guard_tr) || guard_tr < 0L) {
    stop("'guard_tr' must be a non-negative integer", call. = FALSE)
  }
  if (!is.null(seed)) set.seed(seed)

  obs_list <- .as_obs_ids_list(obs_ids)

  if (method == "blocked_time") {
    if (is.null(obs_list)) {
      folds <- .blocked_time_folds(length(obs_ids), k = k, guard_tr = guard_tr)
      return(list(
        method = "blocked_time",
        axis = "observation",
        folds = folds,
        n_folds = length(folds),
        fold_ids = names(folds)
      ))
    }
    folds <- lapply(obs_list, function(v) .blocked_time_folds(length(v), k = k, guard_tr = guard_tr))
    fold_ids <- unique(unlist(lapply(folds, names)))
    out_folds <- lapply(fold_ids, function(fid) {
      per_subject <- lapply(names(folds), function(subj) folds[[subj]][[fid]])
      names(per_subject) <- names(folds)
      per_subject
    })
    names(out_folds) <- fold_ids
    return(list(
      method = "blocked_time",
      axis = "observation",
      folds = out_folds,
      n_folds = length(out_folds),
      fold_ids = fold_ids
    ))
  }

  # method == "run"
  if (is.null(obs_list)) {
    ids <- as.character(obs_ids)
    if (anyNA(ids)) stop("'obs_ids' contains NA", call. = FALSE)
    run_ids <- unique(ids)
    if (length(run_ids) < 2) stop("Need >= 2 unique run ids for method='run'", call. = FALSE)
    folds <- lapply(run_ids, function(run) {
      test_idx <- which(ids == run)
      train_idx <- which(ids != run)
      if (guard_tr > 0L) {
        guard <- unique(unlist(lapply(test_idx, function(i) seq.int(max(1L, i - guard_tr), min(length(ids), i + guard_tr)))))
        train_idx <- setdiff(train_idx, guard)
      }
      list(train_idx = train_idx, test_idx = test_idx)
    })
    names(folds) <- as.character(run_ids)
    return(list(
      method = "run",
      axis = "observation",
      folds = folds,
      n_folds = length(folds),
      fold_ids = names(folds)
    ))
  }

  # Per-subject run ids
  obs_list <- lapply(obs_list, function(v) {
    if (anyNA(v)) stop("Per-subject 'obs_ids' contains NA", call. = FALSE)
    v
  })
  run_sets <- lapply(obs_list, unique)
  fold_ids <- if (id_policy == "intersection") {
    Reduce(intersect, run_sets)
  } else {
    Reduce(union, run_sets)
  }
  if (length(fold_ids) < 2) {
    stop(
      sprintf("Need >= 2 fold ids for method='run' (id_policy='%s')", id_policy),
      call. = FALSE
    )
  }

  folds <- lapply(fold_ids, function(fid) {
    per_subject <- lapply(names(obs_list), function(subj) {
      ids <- obs_list[[subj]]
      if (!fid %in% ids) {
        if (id_policy == "intersection") {
          # Should not happen by construction
          return(list(train_idx = integer(0), test_idx = integer(0)))
        }
        return(list(train_idx = seq_along(ids), test_idx = integer(0)))
      }
      test_idx <- which(ids == fid)
      train_idx <- which(ids != fid)
      if (guard_tr > 0L) {
        guard <- unique(unlist(lapply(test_idx, function(i) seq.int(max(1L, i - guard_tr), min(length(ids), i + guard_tr)))))
        train_idx <- setdiff(train_idx, guard)
      }
      if (length(train_idx) < 1L) {
        stop(
          sprintf("guard_tr too large for subject '%s' fold '%s': training set empty", subj, fid),
          call. = FALSE
        )
      }
      list(train_idx = train_idx, test_idx = test_idx)
    })
    names(per_subject) <- names(obs_list)
    per_subject
  })
  names(folds) <- as.character(fold_ids)

  if (id_policy == "union") {
    missing_any <- vapply(folds, function(f) any(vapply(f, function(s) length(s$test_idx) == 0L, logical(1))), logical(1))
    if (any(missing_any)) {
      warning(
        "Some folds have empty test sets for one or more subjects (id_policy='union'); downstream evaluation should account for this.",
        call. = FALSE
      )
    }
  }

  list(
    method = "run",
    axis = "observation",
    folds = folds,
    n_folds = length(folds),
    fold_ids = names(folds),
    id_policy = id_policy
  )
}


#' Run Cross-Validated Alignment
#'
#' Convenience function to run alignment with cross-validation and
#' collect per-fold results.
#'
#' @param data AlignmentData object.
#' @param method Alignment method.
#' @param cv_folds CV fold information from create_cv_folds, or a string
#'   ("loso", "kfold") to create folds automatically.
#' @param k Number of folds if cv_folds is "kfold".
#' @param reference Reference selection method.
#' @param obs_labels Optional shared observation labels; passed through to
#'   \code{\link{as_alignment_data}} when \code{data} is not already an
#'   \code{\link{AlignmentData}}.
#' @param ... Additional arguments passed to fit_alignment.
#'
#' @return List with:
#'   \itemize{
#'     \item result - Combined AlignmentResult with all CV transforms
#'     \item fold_results - List of per-fold results
#'     \item cv_info - CV configuration info
#'   }
#'
#' @export
run_cv_alignment <- function(data,
                             method = "procrustes",
                             cv_folds = "loso",
                             k = 5,
                             reference = "medoid",
                             obs_labels = NULL,
                             ...) {
  if (!inherits(data, "AlignmentData")) {
    data <- as_alignment_data(data, obs_labels = obs_labels)
  } else {
    if (!is.null(obs_labels)) {
      if (!is.null(data@obs_labels) &&
          !identical(as.character(data@obs_labels), as.character(obs_labels))) {
        stop("obs_labels supplied but AlignmentData already has different obs_labels", call. = FALSE)
      }
      if (is.null(data@obs_labels)) {
        data@obs_labels <- obs_labels
      }
    }
  }

  # Create folds if needed
  if (is.character(cv_folds)) {
    cv_folds <- create_cv_folds(data, method = cv_folds, k = k)
  }

  subjects <- data@subjects
  n_folds <- cv_folds$n_folds

  # Storage
  all_transforms <- list()
  all_aligned <- list()
  fold_results <- list()

  # Run each fold
  for (fold_name in names(cv_folds$folds)) {
    fold <- cv_folds$folds[[fold_name]]
    train_idx <- fold$train
    test_idx <- fold$test

    # Fit on training data
    result <- fit_alignment(
      data = data,
      method = method,
      reference = reference,
      cv = "none",
      train_idx = train_idx,
      ...
    )

    # Apply to test subjects
    test_subjects <- subjects[test_idx]
    for (subj in test_subjects) {
      # Get transform for test subject
      if (subj %in% names(result@model@transforms)) {
        transform <- result@model@transforms[[subj]]
      } else {
        # Need to fit for this subject
        test_result <- apply_alignment(
          result@model,
          data[subj],
          warn_leakage = FALSE
        )
        transform <- test_result@model@transforms[[subj]]
      }

      all_transforms[[subj]] <- transform
      all_aligned[[subj]] <- transform %*% get_subject_data(data, subj)
    }

    fold_results[[fold_name]] <- result
  }

  # Build combined model
  # Fit on all data for reference
  full_fit <- fit_alignment(
    data = data,
    method = method,
    reference = reference,
    cv = "none",
    compute_quality = FALSE,
    ...
  )

  reference_kind <- if (.is_matrixish(reference)) {
    "template"
  } else if (is.character(reference) && length(reference) == 1) {
    if (reference %in% c("medoid", "centroid", "consensus")) "data_driven" else "fixed_subject"
  } else {
    "unknown"
  }

  anchor_common <- reference_kind %in% c("fixed_subject", "template")

  combined_model <- AlignmentModel(
    transforms = all_transforms,
    reference = if (anchor_common) full_fit@model@reference else "fold_specific",
    reference_data = if (anchor_common) full_fit@model@reference_data else NULL,
    method = method,
    space_from = full_fit@model@space_from,
    space_to = full_fit@model@space_to,
    params = list(...),
    method_state = if (anchor_common) full_fit@model@method_state else list(),
    train_subjects = subjects
  )

  quality <- .compute_basic_quality(data, all_aligned, combined_model)

  combined_result <- AlignmentResult(
    model = combined_model,
    aligned = all_aligned,
    quality = quality,
    cv_info = list(
      method = cv_folds$method,
      axis = cv_folds$axis %||% "subject",
      n_folds = n_folds,
      fold_assignments = cv_folds$assignments,
      folds = cv_folds$folds,
      reference_kind = reference_kind,
      anchor_common = anchor_common,
      anchor_note = if (!anchor_common) {
        "Aligned outputs are in fold-specific anchor spaces; do not use for group-level comparisons without mapping to a common anchor."
      } else {
        NULL
      }
    )
  )

  list(
    result = combined_result,
    fold_results = fold_results,
    cv_info = cv_folds
  )
}


#' Check if Result is Cross-Validated
#'
#' @param result An AlignmentResult.
#'
#' @return Logical indicating if result used CV.
#'
#' @export
is_cv_result <- function(result) {
  cv_info <- get_cv_info(result)
  !is.null(cv_info$method) && cv_info$method != "none"
}


#' Get Fold Assignments from Result
#'
#' @param result An AlignmentResult.
#'
#' @return Named vector of fold assignments, or NULL if not CV.
#'
#' @export
get_fold_assignments <- function(result) {
  cv_info <- get_cv_info(result)
  cv_info$fold_assignments
}


#' Check if Result/Model Has a Common Anchor Space
#'
#' In cross-validation, each fold can define a different reference/anchor
#' (e.g., medoid/consensus computed from training data only). In that case,
#' the aligned outputs are expressed in fold-specific spaces and should not be
#' used for group-level comparisons without a leak-free mapping into a single
#' common anchor.
#'
#' @param x An \code{AlignmentResult} or \code{AlignmentModel}.
#'
#' @return Logical; TRUE if outputs are expressed in a single anchor space.
#'
#' @export
has_common_anchor <- function(x) {
  if (inherits(x, "AlignmentResult")) {
    cv_info <- get_cv_info(x)
    # No CV info implies a single fit (common anchor).
    if (length(cv_info) == 0 || is.null(cv_info$method) || cv_info$method == "none") {
      return(TRUE)
    }
    if (!is.null(cv_info$anchor_common)) {
      return(isTRUE(cv_info$anchor_common))
    }
    # Fallback: treat fold-specific sentinel as non-common.
    model <- get_model(x)
    return(!isTRUE(identical(model@reference, "fold_specific")))
  }

  if (inherits(x, "AlignmentModel")) {
    return(!isTRUE(identical(x@reference, "fold_specific")))
  }

  stop("'x' must be an AlignmentResult or AlignmentModel", call. = FALSE)
}


#' Validate That Outputs Are in a Single Anchor Space
#'
#' Convenience guard for downstream group analysis. If the result/model does
#' not have a common anchor, you typically need a fixed/external anchor or an
#' outer split plus a leak-free mapping from fold-specific spaces into a common
#' analysis space.
#'
#' @param x An \code{AlignmentResult} or \code{AlignmentModel}.
#' @param action What to do if no common anchor is detected:
#'   \itemize{
#'     \item "warn" - Issue a warning (default)
#'     \item "error" - Throw an error
#'     \item "silent" - Do nothing
#'   }
#' @param context Short context string used in messages.
#'
#' @return Invisibly returns TRUE if valid; otherwise warns/errors.
#'
#' @export
validate_common_anchor <- function(x,
                                   action = c("warn", "error", "silent"),
                                   context = "group analysis") {
  action <- match.arg(action)

  if (has_common_anchor(x)) {
    return(invisible(TRUE))
  }

  msg <- paste0(
    "No common anchor detected for ", context, ": aligned outputs are in fold-specific spaces. ",
    "Use a fixed/external reference, or perform an outer split and map fold-specific anchors ",
    "into a single analysis space before computing group statistics."
  )

  if (action == "error") {
    stop(msg, call. = FALSE)
  } else if (action == "warn") {
    warning(msg, call. = FALSE)
  }

  invisible(FALSE)
}
