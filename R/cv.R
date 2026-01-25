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
      assignments = assignments,
      folds = folds,
      n_folds = k,
      groups = groups
    ))
  }
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
                             ...) {
  if (!inherits(data, "AlignmentData")) {
    data <- as_alignment_data(data)
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

  combined_model <- AlignmentModel(
    transforms = all_transforms,
    reference = full_fit@model@reference,
    reference_data = full_fit@model@reference_data,
    method = method,
    space_from = full_fit@model@space_from,
    space_to = full_fit@model@space_to,
    params = list(...),
    method_state = full_fit@model@method_state,
    train_subjects = subjects
  )

  quality <- .compute_basic_quality(data, all_aligned, combined_model)

  combined_result <- AlignmentResult(
    model = combined_model,
    aligned = all_aligned,
    quality = quality,
    cv_info = list(
      method = cv_folds$method,
      n_folds = n_folds,
      fold_assignments = cv_folds$assignments
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
