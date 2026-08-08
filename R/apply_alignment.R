#' Apply an Alignment Model to New Data
#'
#' Apply a previously fitted alignment model to new data. This enables
#' out-of-sample alignment for cross-validation or applying trained
#' alignments to new subjects.
#'
#' @family alignment_workflow
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
#' # Fit on training subjects
#' set.seed(1)
#' train_data <- AlignmentData(list(
#'   "sub-01" = matrix(rnorm(30), 6, 5),
#'   "sub-02" = matrix(rnorm(30), 6, 5)
#' ))
#' result <- fit_alignment(train_data, method = "procrustes")
#'
#' # Apply to new subject
#' new_data <- AlignmentData(list(
#'   "sub-03" = matrix(rnorm(30), 6, 5)
#' ))
#' new_result <- apply_alignment(result, new_data)
#'
#' @seealso \code{\link{fit_alignment}}, \code{\link{AlignmentModel}}
#'
#' @export
apply_alignment <- function(model,
                            new_data,
                            fit_new = TRUE,
                            warn_leakage = TRUE,
                            ...) {
  model <- .ensure_model(model, what = "model")

  caps <- aligner_capabilities(model@method)

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

  if (!is.null(caps) && isFALSE(caps$supports_new_data %||% TRUE) && length(existing_subjects) > 0) {
    stop(
      sprintf(
        "Method '%s' does not support applying existing transforms to new data. ",
        model@method
      ),
      "Pass only new subjects (not already in the model), or use a method that declares capabilities$supports_new_data = TRUE.",
      call. = FALSE
    )
  }

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
      ), call. = FALSE)
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

      # Validate method requirements for fitting transforms on the new subjects
      # (design/geometry/guidance). This keeps errors user-friendly and avoids
      # cryptic failures inside apply_fn/fit_fn.
      new_subset <- new_data[match(new_subjects, new_data@subjects)]
      .validate_aligner_requirements(model@method, new_subset)

      # Fit transforms for new subjects
      for (subj in new_subjects) {
        new_transform <- .fit_transform_for_subject(
          aligner, model, new_data, subj, ...
        )
        subj_data <- get_subject_data(new_data, subj)

        if (isTRUE(model@method_state$restrict_to_identified %||% FALSE)) {
          tol <- model@method_state$restrict_tol %||% sqrt(.Machine$double.eps)
          new_transform <- canonicalize_orthogonal_operator(
            Q = new_transform,
            x = subj_data,
            convention = "left",
            tol = tol
          )
        }
        new_transforms[[subj]] <- new_transform

        aligned[[subj]] <- apply_transform(new_transform, subj_data)
      }
    }
  }

  # Create updated model with new transforms
  all_transforms <- c(model@transforms, new_transforms)
  updated_provenance <- model@provenance
  updated_provenance$shared_space_coordinate_id <- .model_coordinate_id(model)
  updated_model <- new("AlignmentModel",
    transforms = all_transforms,
    reference = model@reference,
    reference_data = model@reference_data,
    method = model@method,
    space_from = model@space_from,
    space_to = model@space_to,
    provenance = updated_provenance,
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
  returns <- aligner$capabilities$returns %||% "operator"

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

    if (is.list(apply_result$transforms) && length(apply_result$transforms) > 0) {
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
      sprintf(
        "apply_fn for method '%s' did not return transforms (and no aligned embedding for returns='embedding')",
        aligner$name %||% "<unknown>"
      ),
      call. = FALSE
    )
  }

  # Default: use fit_fn with single subject
  # Need to rebuild enough state for the aligner
  fit_result <- aligner$fit_fn(
    data = data,
    reference = model@reference_data %||% model@reference,
    train_idx = subject_idx,
    ...
  )

  if (is.list(fit_result$transforms) && subject %in% names(fit_result$transforms)) {
    return(fit_result$transforms[[subject]])
  }
  if (identical(returns, "embedding")) {
    z <- fit_result$aligned %||% NULL
    if (is.list(z) && subject %in% names(z)) {
      return(.new_embedding_transform(z[[subject]], subject = subject))
    }
  }

  stop(
    sprintf("Method '%s' did not produce a transform for subject '%s'", aligner$name, subject),
    call. = FALSE
  )
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
#' @examples
#' X <- matrix(1:6, 3, 2)
#' Q <- diag(3)
#' apply_transform(Q, X)
#'
#' @export
apply_transform <- function(transform, data) {
  if (.is_embedding_transform(transform)) {
    if (!.is_matrixish(data)) data <- as.matrix(data)
    z <- transform$aligned
    if (ncol(z) != ncol(data)) {
      stop(sprintf(
        "Embedding transform cannot be applied: embedding has %d observations but data has %d",
        ncol(z), ncol(data)
      ), call. = FALSE)
    }
    return(z)
  }

  if (.is_low_rank_transform(transform)) {
    if (!.is_matrixish(data)) data <- as.matrix(data)
    if (nrow(data) != nrow(transform$V)) {
      stop(sprintf(
        "Transform dimension mismatch: low-rank source=%d, data has %d rows",
        nrow(transform$V), nrow(data)
      ), call. = FALSE)
    }
    return(transform$U %*% (t(transform$V) %*% data))
  }

  if (!.is_matrixish(transform)) {
    stop("Transform must be an operator matrix/Matrix or a supported transform object", call. = FALSE)
  }
  if (!.is_matrixish(data)) data <- as.matrix(data)

  # Validate dimensions
  if (ncol(transform) != nrow(data)) {
    stop(sprintf(
      "Transform dimension mismatch: transform is %d x %d, data is %d x %d",
      nrow(transform), ncol(transform), nrow(data), ncol(data)
    ), call. = FALSE)
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
#' @examples
#' set.seed(1)
#' adat <- AlignmentData(list(
#'   s1 = matrix(rnorm(20), 4, 5),
#'   s2 = matrix(rnorm(20), 4, 5)
#' ))
#' res <- fit_alignment(adat, method = "procrustes", reference = "s1", compute_quality = FALSE)
#' mdl <- get_model(res)
#' inv <- inverse_transform(mdl, "s2", method = "transpose")
#' dim(inv)
#'
#' @export
inverse_transform <- function(model,
                              subject,
                              method = "auto",
                              tol = sqrt(.Machine$double.eps),
                              lambda = 1e-6) {
  transform <- get_transform(model, subject)
  caps <- aligner_capabilities(model@method)

  if (.is_embedding_transform(transform)) {
    stop(
      sprintf(
        "Method '%s' returns embedding transforms; inverse_transform() is not defined",
        model@method
      ),
      call. = FALSE
    )
  }
  if (.is_low_rank_transform(transform)) {
    if (identical(method, "auto")) {
      method <- "pinv"
    }
    if (method == "pinv") return(.pseudoinverse_low_rank(transform, tol = tol))
    if (method == "ridge") return(.ridge_inverse_low_rank(transform, lambda = lambda))
    stop(
      "inverse_transform() for low-rank operator transforms supports method='pinv' or method='ridge'",
      call. = FALSE
    )
  }

  if (!is.null(caps) && identical(caps$transform_type, "ot")) {
    stop(
      sprintf(
        "Method '%s' returns OT-style couplings (non-bijective); inverse_transform() is not defined",
        model@method
      ),
      call. = FALSE
    )
  }

  method <- .resolve_inverse_method(method, caps, transform, model@method)

  if (method == "transpose") return(t(transform))

  if (method == "pinv") return(.pseudoinverse(transform, tol = tol))
  if (method == "ridge") return(.ridge_inverse(transform, lambda = lambda))

  if (method != "solve") {
    stop(sprintf("Unknown inverse method: %s", method), call. = FALSE)
  }

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
}

.resolve_inverse_method <- function(method, caps, transform, method_name) {
  if (!identical(method, "auto")) return(method)

  if (!is.null(caps)) {
    if (!isTRUE(caps$returns_invertible)) {
      stop(
        sprintf(
          "Method '%s' does not declare invertible operators; use method='pinv' or method='ridge' for an approximate inverse",
          method_name
        ),
        call. = FALSE
      )
    }
    if (caps$transform_type %in% c("orthogonal", "permutation")) {
      return("transpose")
    }
    return("solve")
  }

  if (nrow(transform) == ncol(transform)) return("solve")
  stop(
    "Transform is not square; no exact inverse. Use method='pinv' or method='ridge' for an approximate inverse.",
    call. = FALSE
  )
}

.pseudoinverse_low_rank <- function(transform, tol = sqrt(.Machine$double.eps)) {
  sv <- .low_rank_operator_svd(transform, context = "low-rank pseudoinverse")
  d <- sv$d
  if (!length(d)) {
    return(.new_low_rank_transform(sv$V, sv$U))
  }
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol <= 0) {
    stop("'tol' must be a single positive number for method='pinv'", call. = FALSE)
  }
  cutoff <- max(nrow(sv$U), nrow(sv$V)) * max(d) * tol
  d_inv <- ifelse(d > cutoff, 1 / d, 0)
  U_new <- sweep(sv$V, 2L, d_inv, "*")
  .new_low_rank_transform(U_new, sv$U)
}

.ridge_inverse_low_rank <- function(transform, lambda = 1e-6) {
  if (!is.numeric(lambda) || length(lambda) != 1L || !is.finite(lambda) || lambda <= 0) {
    stop("'lambda' must be a single positive number for method='ridge'", call. = FALSE)
  }
  sv <- .low_rank_operator_svd(transform, context = "low-rank ridge inverse")
  d <- sv$d
  if (!length(d)) {
    return(.new_low_rank_transform(sv$V, sv$U))
  }
  w <- d / (d^2 + lambda)
  U_new <- sweep(sv$V, 2L, w, "*")
  .new_low_rank_transform(U_new, sv$U)
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

#' Lift Aligned Data Into an Ambient Feature Space
#'
#' Many aligners (e.g., manifoldalign projection methods) map high-dimensional
#' feature data `X (p x n)` into a shared bottleneck space `Z (k x n)` using
#' non-square operators `A (k x p)`. This helper maps aligned data back into an
#' ambient feature space using an approximate inverse (pseudoinverse or ridge).
#'
#' Typical uses:
#' - Lift aligned latent data into the reference subject's voxel space.
#' - "Paint" a shared latent representation into a particular subject's voxel
#'   space for visualization or decoding.
#'
#' @param x An \code{\link{AlignmentResult}} (preferred) or \code{\link{AlignmentModel}}.
#' @param to Target ambient space to lift into:
#'   \itemize{
#'     \item \code{"reference"}: use the model's reference subject space.
#'     \item \code{"subject"}: use a specific \code{subject}.
#'     \item \code{"each"}: lift each aligned matrix into its own subject space.
#'   }
#' @param subject Subject ID used when \code{to="subject"}.
#' @param aligned Optional aligned data to lift. If \code{x} is an
#'   \code{AlignmentResult} and \code{aligned} is \code{NULL}, uses
#'   \code{\link{get_aligned}(x)}.
#'   Can be a single matrix or a named list of matrices.
#' @param inverse How to compute the inverse operator; forwarded to
#'   \code{\link{inverse_transform}}. Use \code{"pinv"} or \code{"ridge"} for
#'   bottleneck (non-square) operators.
#' @param tol Tolerance for \code{inverse="pinv"} (relative to largest singular value).
#' @param lambda Regularization strength for \code{inverse="ridge"}.
#'
#' @return A matrix (if \code{aligned} is a matrix) or a named list of matrices.
#'
#' @examples
#' set.seed(1)
#' adat <- AlignmentData(list(
#'   s1 = matrix(rnorm(20), 4, 5),
#'   s2 = matrix(rnorm(20), 4, 5)
#' ))
#' res <- fit_alignment(adat, method = "procrustes", reference = "s1", compute_quality = FALSE)
#' lifted <- lift_aligned(res, to = "reference", inverse = "transpose")
#' lapply(lifted, dim)
#'
#' @export
lift_aligned <- function(x,
                         to = c("reference", "subject", "each"),
                         subject = NULL,
                         aligned = NULL,
                         inverse = c("pinv", "ridge", "transpose", "solve", "auto"),
                         tol = sqrt(.Machine$double.eps),
                         lambda = 1e-6) {
  to <- match.arg(to)
  inverse <- match.arg(inverse)

  if (inherits(x, "AlignmentResult") && is.null(aligned)) {
    aligned <- get_aligned(x)
  }
  model <- .ensure_model(x, what = "x")

  if (is.null(aligned)) {
    stop(
      "No aligned data provided. Pass an AlignmentResult, or supply 'aligned=' explicitly.",
      call. = FALSE
    )
  }

  single <- FALSE
  if (!is.list(aligned)) {
    single <- TRUE
    aligned <- list(.aligned = aligned)
  } else if (length(aligned) == 0L) {
    stop("'aligned' is empty", call. = FALSE)
  }

  if (to == "reference") {
    target <- model@reference
    if (!(is.character(target) && length(target) == 1L && target %in% names(model@transforms))) {
      stop(
        "Model does not have a subject-id reference; use to='subject' and supply 'subject=' explicitly.",
        call. = FALSE
      )
    }
    inv <- inverse_transform(model, target, method = inverse, tol = tol, lambda = lambda)
    out <- lapply(aligned, function(z) apply_transform(inv, z))
  } else if (to == "subject") {
    if (is.null(subject) || !is.character(subject) || length(subject) != 1L) {
      stop("to='subject' requires a single character 'subject'", call. = FALSE)
    }
    if (!subject %in% names(model@transforms)) {
      stop(sprintf("Unknown subject '%s' in model transforms", subject), call. = FALSE)
    }
    inv <- inverse_transform(model, subject, method = inverse, tol = tol, lambda = lambda)
    out <- lapply(aligned, function(z) apply_transform(inv, z))
  } else {
    nm <- names(aligned)
    if (is.null(nm) || any(nm == "")) {
      stop("to='each' requires 'aligned' to be a named list keyed by subject", call. = FALSE)
    }
    missing <- setdiff(nm, names(model@transforms))
    if (length(missing) > 0) {
      stop(
        sprintf("Missing transforms for subject(s): %s", paste(missing, collapse = ", ")),
        call. = FALSE
      )
    }
    out <- lapply(nm, function(s) {
      inv <- inverse_transform(model, s, method = inverse, tol = tol, lambda = lambda)
      apply_transform(inv, aligned[[s]])
    })
    names(out) <- nm
  }

  if (single) out[[1L]] else out
}
