#' Data Leakage Detection
#'
#' Functions for detecting and warning about potential data leakage in
#' alignment workflows. Leakage occurs when test data influences the
#' alignment transformation, leading to overly optimistic results.
#'
#' @name leakage
NULL


#' Check for Leakage in Alignment Application
#'
#' Check if applying alignment to data might involve leakage because
#' the subjects overlap with training data.
#'
#' @param apply_subjects Character vector of subjects being aligned.
#' @param train_subjects Character vector of subjects used for training.
#' @param context Character string describing the context (for warning message).
#' @param action What to do on leakage detection:
#'   \itemize{
#'     \item "warn" - Issue a warning (default)
#'     \item "error" - Throw an error
#'     \item "silent" - Do nothing
#'   }
#'
#' @return Invisibly returns list with leakage info.
#'
#' @examples
#' \dontrun{
#' # This will warn about leakage
#' .check_leakage(c("sub-01", "sub-02"), c("sub-01", "sub-03"), "test")
#' }
#'
#' @keywords internal
.check_leakage <- function(apply_subjects, train_subjects,
                           context = "alignment",
                           action = c("warn", "error", "silent")) {
  action <- match.arg(action)

  overlap <- intersect(apply_subjects, train_subjects)

  if (length(overlap) == 0) {
    return(invisible(list(
      has_leakage = FALSE,
      overlap = character(0)
    )))
  }

  msg <- sprintf(
    paste0(
      "Potential data leakage in %s: %d of %d subjects (%s) ",
      "were in the training set. Use cross-validation or separate ",
      "train/test sets to avoid inflated accuracy."
    ),
    context,
    length(overlap),
    length(apply_subjects),
    if (length(overlap) <= 3) {
      paste(overlap, collapse = ", ")
    } else {
      paste(c(overlap[1:2], sprintf("... +%d more", length(overlap) - 2)),
        collapse = ", "
      )
    }
  )

  if (action == "error") {
    stop(msg, call. = FALSE)
  } else if (action == "warn") {
    warning(msg, call. = FALSE)
  }

  invisible(list(
    has_leakage = TRUE,
    overlap = overlap,
    n_overlap = length(overlap),
    n_apply = length(apply_subjects),
    n_train = length(train_subjects)
  ))
}


#' Validate CV Setup for Leakage
#'
#' Validate that a cross-validation setup doesn't have leakage.
#'
#' @param cv_folds CV fold info from create_cv_folds.
#' @param reference Reference selection method or subject.
#'
#' @return TRUE if valid, otherwise throws error.
#'
#' @export
validate_cv_setup <- function(cv_folds, reference = "medoid") {
  # Check train/test disjointness for every fold
  for (fold_name in names(cv_folds$folds)) {
    fold <- cv_folds$folds[[fold_name]]
    overlap <- intersect(fold$train, fold$test)
    if (length(overlap) > 0) {
      stop(sprintf(
        "Fold '%s' has %d overlapping indices between train and test sets",
        fold_name, length(overlap)
      ), call. = FALSE)
    }
  }

  # Check that reference selection doesn't leak into test
  if (is.character(reference) && length(reference) == 1) {
    if (reference %in% c("medoid", "centroid", "consensus")) {
      # These are computed from training data only - OK
      return(TRUE)
    }

    # Fixed subject reference - check it's always in training
    for (fold_name in names(cv_folds$folds)) {
      fold <- cv_folds$folds[[fold_name]]
      fold_subjects <- names(cv_folds$assignments)[fold$train]

      if (!reference %in% fold_subjects) {
        stop(sprintf(
          "Fixed reference subject '%s' is in test set for fold '%s'. This causes leakage. Use reference='medoid' instead.",
          reference, fold_name
        ), call. = FALSE)
      }
    }
  }

  TRUE
}


#' Summarize Leakage Risk
#'
#' Analyze an alignment workflow for potential leakage risks.
#'
#' @param model AlignmentModel or AlignmentResult.
#' @param test_subjects Optional character vector of intended test subjects.
#'
#' @return List with leakage risk assessment.
#'
#' @export
assess_leakage_risk <- function(model, test_subjects = NULL) {
  if (inherits(model, "AlignmentResult")) {
    cv_info <- get_cv_info(model)
    model <- get_model(model)
  } else {
    cv_info <- list()
  }

  risk <- list(
    overall_risk = "low",
    issues = character(0),
    recommendations = character(0)
  )

  train_subjects <- model@train_subjects

  # Check if CV was used
  if (length(cv_info) == 0 || is.null(cv_info$method) ||
      cv_info$method == "none") {
    risk$issues <- c(risk$issues,
      "No cross-validation was used"
    )
    risk$recommendations <- c(risk$recommendations,
      "Use cv='loso' or cv='kfold' in fit_alignment for held-out evaluation"
    )
    risk$overall_risk <- "medium"
  }

  # Check test/train overlap if test subjects provided
  if (!is.null(test_subjects)) {
    overlap <- intersect(test_subjects, train_subjects)
    if (length(overlap) > 0) {
      risk$issues <- c(risk$issues,
        sprintf("%d test subjects were in training set", length(overlap))
      )
      risk$recommendations <- c(risk$recommendations,
        "Ensure test subjects are held out during fitting"
      )
      risk$overall_risk <- "high"
    }
  }

  # Check aligner CV support
  caps <- aligner_capabilities(model@method)
  if (!is.null(caps) && !isTRUE(caps$supports_cv)) {
    risk$issues <- c(risk$issues,
      sprintf("Method '%s' may not fully support CV", model@method)
    )
    risk$recommendations <- c(risk$recommendations,
      "Verify method properly handles train/test splits"
    )
    if (risk$overall_risk == "low") {
      risk$overall_risk <- "medium"
    }
  }

  risk
}


#' Check for Observation-Level Leakage
#'
#' Check if observation-level train/test splits have overlapping indices,
#' which would constitute data leakage in observation-axis cross-validation.
#'
#' @param train_obs Named list (keyed by subject) of training observation
#'   labels or indices.
#' @param test_obs Named list (keyed by subject) of test observation
#'   labels or indices.
#' @param action What to do on leakage detection: "warn", "error", or "silent".
#'
#' @return Invisibly returns list with leakage info per subject.
#'
#' @export
check_obs_leakage <- function(train_obs, test_obs,
                              action = c("warn", "error", "silent")) {
  action <- match.arg(action)

  common_subjects <- intersect(names(train_obs), names(test_obs))
  leaky_subjects <- character(0)

  for (subj in common_subjects) {
    overlap <- intersect(train_obs[[subj]], test_obs[[subj]])
    if (length(overlap) > 0) {
      leaky_subjects <- c(leaky_subjects, subj)
    }
  }

  if (length(leaky_subjects) == 0) {
    return(invisible(list(has_leakage = FALSE, leaky_subjects = character(0))))
  }

  msg <- sprintf(
    "Observation-level leakage: %d of %d subjects have overlapping train/test observations (%s).",
    length(leaky_subjects),
    length(common_subjects),
    if (length(leaky_subjects) <= 3) {
      paste(leaky_subjects, collapse = ", ")
    } else {
      paste(c(leaky_subjects[1:2], sprintf("... +%d more", length(leaky_subjects) - 2)),
        collapse = ", "
      )
    }
  )

  if (action == "error") {
    stop(msg, call. = FALSE)
  } else if (action == "warn") {
    warning(msg, call. = FALSE)
  }

  invisible(list(
    has_leakage = TRUE,
    leaky_subjects = leaky_subjects,
    n_leaky = length(leaky_subjects),
    n_subjects = length(common_subjects)
  ))
}


#' Print Leakage Risk Assessment
#'
#' @param risk Risk assessment from assess_leakage_risk.
#'
#' @export
print_leakage_assessment <- function(risk) {
  cat("Leakage Risk Assessment\n")
  cat("=======================\n\n")

  color <- switch(risk$overall_risk,
    "low" = "green",
    "medium" = "yellow",
    "high" = "red"
  )

  cat(sprintf("Overall Risk: %s\n\n", toupper(risk$overall_risk)))

  if (length(risk$issues) > 0) {
    cat("Issues Found:\n")
    for (issue in risk$issues) {
      cat(sprintf("  - %s\n", issue))
    }
    cat("\n")
  }

  if (length(risk$recommendations) > 0) {
    cat("Recommendations:\n")
    for (rec in risk$recommendations) {
      cat(sprintf("  - %s\n", rec))
    }
  }

  if (length(risk$issues) == 0) {
    cat("No leakage issues detected.\n")
  }

  invisible(risk)
}
