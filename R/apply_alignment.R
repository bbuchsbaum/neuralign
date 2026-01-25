#' Apply an Alignment Model to New Data
#'
#' Apply a previously fitted alignment model to new data. This enables
#' out-of-sample alignment for cross-validation or applying trained
#' alignments to new subjects.
#'
#' @param model An AlignmentModel object from \code{\link{fit_alignment}},
#'   or an AlignmentResult from which the model will be extracted.
#' @param new_data AlignmentData object, or a named list of matrices.
#'   For subjects already in the model, their existing transforms are used.
#'   For new subjects, transforms are computed using the model's reference.
#' @param fit_new Logical; if TRUE (default), fit transforms for subjects
#'   not in the model. If FALSE, only apply to subjects already in model.
#' @param warn_leakage Logical; if TRUE (default), warn when applying to
#'   subjects that were used in training.
#' @param ... Additional arguments passed to the aligner's apply function.
#'
#' @return An AlignmentResult object with aligned data for the new subjects.
#'
#' @examples
#' \dontrun{
#' # Fit on training subjects
#' train_data <- AlignmentData(list(
#'   "sub-01" = matrix(rnorm(100*50), 100, 50),
#'   "sub-02" = matrix(rnorm(100*50), 100, 50)
#' ))
#' result <- fit_alignment(train_data, method = "procrustes")
#'
#' # Apply to new subject
#' new_data <- AlignmentData(list(
#'   "sub-03" = matrix(rnorm(100*50), 100, 50)
#' ))
#' new_result <- apply_alignment(result, new_data)
#' }
#'
#' @seealso \code{\link{fit_alignment}}, \code{\link{AlignmentModel}}
#'
#' @export
apply_alignment <- function(model,
                            new_data,
                            fit_new = TRUE,
                            warn_leakage = TRUE,
                            ...) {
  # Extract model from result if needed
  if (inherits(model, "AlignmentResult")) {
    model <- get_model(model)
  }

  if (!inherits(model, "AlignmentModel")) {
    stop("'model' must be an AlignmentModel or AlignmentResult")
  }

  # Coerce data to AlignmentData if needed
  if (!inherits(new_data, "AlignmentData")) {
    new_data <- as_alignment_data(new_data)
  }

  # Check for leakage
  if (warn_leakage) {
    train_subjects <- model@train_subjects
    overlap <- intersect(new_data@subjects, train_subjects)
    if (length(overlap) > 0) {
      .check_leakage(new_data@subjects, train_subjects, "apply_alignment")
    }
  }

  # Separate subjects into those with existing transforms and new subjects
  existing_subjects <- intersect(new_data@subjects, names(model@transforms))
  new_subjects <- setdiff(new_data@subjects, names(model@transforms))

  aligned <- list()
  new_transforms <- list()

  # Apply existing transforms
  for (subj in existing_subjects) {
    transform <- model@transforms[[subj]]
    subj_data <- get_subject_data(new_data, subj)
    aligned[[subj]] <- transform %*% subj_data
  }

  # Fit and apply for new subjects
  if (length(new_subjects) > 0) {
    if (!fit_new) {
      warning(sprintf(
        "No transforms for subjects: %s (fit_new=FALSE)",
        paste(new_subjects, collapse = ", ")
      ))
    } else {
      # Get the aligner
      aligner <- get_aligner(model@method)

      if (is.null(aligner)) {
        stop(sprintf(
          "Method '%s' not registered; cannot fit new subjects",
          model@method
        ))
      }

      # Check if method supports new subjects
      caps <- aligner$capabilities
      if (!is.null(caps$supports_new_subject) && !caps$supports_new_subject) {
        stop(sprintf(
          "Method '%s' does not support fitting transforms for new subjects. ",
          "Use fit_new=FALSE to apply only to subjects with existing transforms.",
          model@method
        ))
      }

      # Fit transforms for new subjects
      for (subj in new_subjects) {
        new_transform <- .fit_transform_for_subject(
          aligner, model, new_data, subj, ...
        )
        new_transforms[[subj]] <- new_transform

        subj_data <- get_subject_data(new_data, subj)
        aligned[[subj]] <- new_transform %*% subj_data
      }
    }
  }

  # Create updated model with new transforms
  all_transforms <- c(model@transforms, new_transforms)
  updated_model <- new("AlignmentModel",
    transforms = all_transforms,
    reference = model@reference,
    reference_data = model@reference_data,
    method = model@method,
    space_from = model@space_from,
    space_to = model@space_to,
    provenance = model@provenance,
    method_state = model@method_state,
    train_subjects = model@train_subjects
  )

  AlignmentResult(
    model = updated_model,
    aligned = aligned,
    quality = list(),
    cv_info = list(method = "applied")
  )
}


#' Internal: Fit Transform for a Single Subject
#' @keywords internal
.fit_transform_for_subject <- function(aligner, model, data, subject, ...) {
  subject_idx <- match(subject, data@subjects)

  # If aligner has custom apply function, use it
  if (!is.null(aligner$apply_fn)) {
    fit_result <- list(
      transforms = model@transforms,
      reference_data = model@reference_data,
      space_from = model@space_from,
      space_to = model@space_to,
      method_state = model@method_state
    )

    apply_result <- aligner$apply_fn(
      fit_result = fit_result,
      new_data = data[subject_idx],
      ...
    )
    return(apply_result$transforms[[1]])
  }

  # Default: use fit_fn with single subject
  # Need to rebuild enough state for the aligner
  fit_result <- aligner$fit_fn(
    data = data,
    reference = model@reference_data %||% model@reference,
    train_idx = subject_idx,
    ...
  )

  fit_result$transforms[[subject]]
}


#' Apply Transform to Data
#'
#' Apply a single transform operator to data using left-multiply convention.
#'
#' @param transform The transform operator (matrix, target x source).
#' @param data Data to transform (matrix, features x observations).
#'
#' @return Transformed data.
#'
#' @export
apply_transform <- function(transform, data) {
  # Validate dimensions
  if (ncol(transform) != nrow(data)) {
    stop(sprintf(
      "Transform dimension mismatch: transform is %d x %d, data is %d x %d",
      nrow(transform), ncol(transform), nrow(data), ncol(data)
    ))
  }

  transform %*% data
}


#' Get Inverse Transform
#'
#' Get the inverse of a transform, if it exists.
#'
#' @param model An AlignmentModel.
#' @param subject Subject ID.
#' @param method How to compute inverse:
#'   \itemize{
#'     \item "transpose" - Use transpose (valid for orthogonal transforms)
#'     \item "solve" - Use matrix solve (valid for invertible transforms)
#'     \item "auto" - Choose based on transform type
#'   }
#'
#' @return The inverse transform operator.
#'
#' @export
inverse_transform <- function(model, subject, method = "auto") {
  transform <- get_transform(model, subject)
  caps <- aligner_capabilities(model@method)

  if (method == "auto") {
    if (!is.null(caps) && isTRUE(caps$returns_invertible)) {
      if (caps$transform_type == "orthogonal") {
        method <- "transpose"
      } else {
        method <- "solve"
      }
    } else {
      method <- "solve"
    }
  }

  if (method == "transpose") {
    return(t(transform))
  } else if (method == "solve") {
    tryCatch(
      solve(transform),
      error = function(e) {
        stop(sprintf(
          "Transform is not invertible: %s",
          conditionMessage(e)
        ))
      }
    )
  } else {
    stop(sprintf("Unknown inverse method: %s", method))
  }
}
