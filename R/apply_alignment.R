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
    stop("'model' must be an AlignmentModel or AlignmentResult", call. = FALSE)
  }

  caps <- aligner_capabilities(model@method)
  if (!is.null(caps) && !identical(caps$returns %||% "operator", "operator")) {
    stop(
      sprintf(
        "Method '%s' does not return operator transforms; apply_alignment() currently supports operators only",
        model@method
      ),
      call. = FALSE
    )
  }

  # Coerce data to AlignmentData if needed
  if (!inherits(new_data, "AlignmentData")) {
    new_data <- as_alignment_data(new_data)
  }

  # Check space compatibility
  if (!spaces_compatible(model@space_from, new_data@space)) {
    warning(sprintf(
      "Space mismatch: model expects '%s' but data is in '%s'",
      .format_space(model@space_from), .format_space(new_data@space)
    ), call. = FALSE)
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
    aligned[[subj]] <- apply_transform(transform, subj_data)
  }

  # Fit and apply for new subjects
  if (length(new_subjects) > 0) {
    if (isTRUE(identical(model@reference, "fold_specific")) &&
        is.null(model@reference_data)) {
      stop(
        paste0(
          "This model uses fold-specific anchors (no common reference); ",
          "cannot fit transforms for new subjects. Fit a non-CV model or use a fixed/external reference, ",
          "or call apply_alignment(..., fit_new=FALSE) to apply only existing transforms."
        ),
        call. = FALSE
      )
    }
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
          "Method '%s' does not support fitting transforms for new subjects. Use fit_new=FALSE to apply only to subjects with existing transforms.",
          model@method
        ), call. = FALSE)
      }

      # Fit transforms for new subjects
      for (subj in new_subjects) {
        new_transform <- .fit_transform_for_subject(
          aligner, model, new_data, subj, ...
        )
        new_transforms[[subj]] <- new_transform

        subj_data <- get_subject_data(new_data, subj)
        aligned[[subj]] <- apply_transform(new_transform, subj_data)
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
  if (!.is_matrixish(transform)) {
    stop("Transform must be a matrix/Matrix operator", call. = FALSE)
  }
  if (!.is_matrixish(data)) {
    data <- as.matrix(data)
  }

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
#'     \item "pinv" - Moore-Penrose pseudoinverse (works for non-square / rank-deficient operators)
#'     \item "ridge" - Tikhonov-regularized inverse (stable approximate inverse; requires \code{lambda > 0})
#'     \item "auto" - Choose based on transform type
#'   }
#' @param tol Tolerance used by \code{method = "pinv"} (relative to the largest singular value).
#' @param lambda Regularization strength used by \code{method = "ridge"}.
#'
#' @return The inverse transform operator.
#'
#' @export
inverse_transform <- function(model,
                              subject,
                              method = "auto",
                              tol = sqrt(.Machine$double.eps),
                              lambda = 1e-6) {
  transform <- get_transform(model, subject)
  caps <- aligner_capabilities(model@method)

  if (!is.null(caps) && identical(caps$transform_type, "ot")) {
    stop(
      sprintf(
        "Method '%s' returns OT-style couplings (non-bijective); inverse_transform() is not defined",
        model@method
      ),
      call. = FALSE
    )
  }

  if (method == "auto") {
    if (!is.null(caps)) {
      if (!isTRUE(caps$returns_invertible)) {
        stop(
          sprintf(
            "Method '%s' does not declare invertible operators; use method='pinv' or method='ridge' for an approximate inverse",
            model@method
          ),
          call. = FALSE
        )
      }
      if (caps$transform_type %in% c("orthogonal", "permutation")) {
        method <- "transpose"
      } else {
        method <- "solve"
      }
    } else {
      if (nrow(transform) == ncol(transform)) {
        method <- "solve"
      } else {
        stop(
          "Transform is not square; no exact inverse. Use method='pinv' or method='ridge' for an approximate inverse.",
          call. = FALSE
        )
      }
    }
  }

  if (method == "transpose") {
    return(t(transform))
  } else if (method == "solve") {
    if (nrow(transform) != ncol(transform)) {
      stop(
        "Transform is not square; no exact inverse. Use method='pinv' or method='ridge' for an approximate inverse.",
        call. = FALSE
      )
    }
    tryCatch(
      solve(transform),
      error = function(e) {
        stop(sprintf(
          "Transform is not invertible: %s. Consider method='pinv' or method='ridge'.",
          conditionMessage(e)
        ), call. = FALSE)
      }
    )
  } else if (method == "pinv") {
    return(.pseudoinverse(transform, tol = tol))
  } else if (method == "ridge") {
    return(.ridge_inverse(transform, lambda = lambda))
  } else {
    stop(sprintf("Unknown inverse method: %s", method), call. = FALSE)
  }
}

.as_dense_matrix <- function(x) {
  if (inherits(x, "Matrix")) return(as.matrix(x))
  if (is.matrix(x)) return(x)
  as.matrix(x)
}

.pseudoinverse <- function(x, tol = sqrt(.Machine$double.eps)) {
  x <- .as_dense_matrix(x)
  svd_result <- svd(x)
  d <- svd_result$d

  if (!length(d)) {
    return(matrix(numeric(0), ncol(x), nrow(x)))
  }

  cutoff <- max(dim(x)) * max(d) * tol
  d_inv <- ifelse(d > cutoff, 1 / d, 0)

  svd_result$v %*% (t(svd_result$u) * d_inv)
}

.ridge_inverse <- function(x, lambda = 1e-6) {
  if (!is.numeric(lambda) || length(lambda) != 1L || !is.finite(lambda) || lambda <= 0) {
    stop("'lambda' must be a single positive number for method='ridge'", call. = FALSE)
  }

  x <- .as_dense_matrix(x)
  n_out <- nrow(x)
  n_in <- ncol(x)

  if (n_out >= n_in) {
    gram <- crossprod(x) + diag(lambda, n_in)
    solve(gram, t(x))
  } else {
    gram <- tcrossprod(x) + diag(lambda, n_out)
    t(x) %*% solve(gram)
  }
}
