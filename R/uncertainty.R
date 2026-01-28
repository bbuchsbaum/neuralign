#' Uncertainty Propagation Helpers
#'
#' Domain-agnostic helpers for propagating simple uncertainty summaries through
#' alignment transforms.
#'
#' @name uncertainty
NULL

#' Diagonal Covariance Propagation Under an Orthogonal Transform
#'
#' Compute the diagonal of a covariance matrix after applying a (typically
#' orthogonal) change of basis.
#'
#' If `x` has diagonal covariance `diag(variances)` and `y = Q %*% x`, then:
#' \deqn{\mathrm{diag}(\mathrm{Cov}(y)) = \mathrm{diag}(Q \, \mathrm{diag}(v) \, Q^\top).}
#' This can be computed without forming the full covariance via element-wise
#' squares:
#' \deqn{\mathrm{diag}(Q \, \mathrm{diag}(v) \, Q^\top) = (Q \odot Q)\, v.}
#'
#' `convention="right"` computes the analogous quantity for `y = t(Q) %*% x`
#' (equivalently, for right-multiplying row vectors by `Q`).
#'
#' @param variances Numeric vector of length `d` (or `d x k` matrix) giving the
#'   diagonal variances in the source basis.
#' @param Q Square `d x d` matrix/Matrix operator, or a list of such operators.
#' @param convention Convention for how `Q` is applied:
#'   `"left"` computes `diag(Q %*% diag(v) %*% t(Q))`;
#'   `"right"` computes `diag(t(Q) %*% diag(v) %*% Q)`.
#'
#' @return If `Q` is a single matrix, returns a numeric vector (or `d x k`
#'   matrix if `variances` is a matrix). If `Q` is a list and `variances` is a
#'   vector, returns a `d x length(Q)` matrix whose columns correspond to each
#'   transform.
#'
#' @examples
#' v <- c(a = 1, b = 2)
#' Q <- matrix(c(0, 1, -1, 0), 2, 2)  # 90-degree rotation (orthogonal)
#' diag_cov_orthogonal(v, Q, convention = "left")
#'
#' @export
diag_cov_orthogonal <- function(variances, Q, convention = c("left", "right")) {
  convention <- match.arg(convention)

  diag_one <- function(v, Qi) {
    if (!is.numeric(v)) {
      stop("'variances' must be numeric", call. = FALSE)
    }
    if (!.is_matrixish(Qi)) {
      stop("'Q' must be a matrix/Matrix operator", call. = FALSE)
    }
    Qi <- as.matrix(Qi)
    if (nrow(Qi) != ncol(Qi)) {
      stop("'Q' must be square", call. = FALSE)
    }

    v_dim <- dim(v)
    if (is.null(v_dim)) {
      if (length(v) != nrow(Qi)) {
        stop("Dimension mismatch: length(variances) must equal nrow(Q)", call. = FALSE)
      }
      if (any(!is.finite(v))) stop("'variances' must be finite", call. = FALSE)
      if (any(v < 0)) stop("'variances' must be non-negative", call. = FALSE)
      nm <- names(v)
      v <- as.numeric(v)
      if (!is.null(nm)) names(v) <- nm
      w <- (Qi * Qi)
      out <- if (convention == "left") w %*% v else t(w) %*% v
      out <- as.numeric(out)
      if (!is.null(nm)) names(out) <- nm
      out
    } else {
      if (nrow(v) != nrow(Qi)) {
        stop("Dimension mismatch: nrow(variances) must equal nrow(Q)", call. = FALSE)
      }
      if (any(!is.finite(v))) stop("'variances' must be finite", call. = FALSE)
      if (any(v < 0)) stop("'variances' must be non-negative", call. = FALSE)
      rn <- rownames(v)
      w <- (Qi * Qi)
      out <- if (convention == "left") w %*% v else t(w) %*% v
      if (!is.null(rn)) rownames(out) <- rn
      out
    }
  }

  if (is.list(Q)) {
    if (length(Q) < 1L) stop("'Q' must be non-empty", call. = FALSE)
    if (!is.null(dim(variances))) {
      stop(
        "When 'Q' is a list, 'variances' must be a numeric vector (not a matrix).",
        call. = FALSE
      )
    }
    out <- vapply(Q, function(Qi) diag_one(variances, Qi), numeric(length(variances)))
    rownames(out) <- names(variances)
    if (!is.null(names(Q)) && all(nzchar(names(Q)))) {
      colnames(out) <- names(Q)
    }
    out
  } else {
    diag_one(variances, Q)
  }
}

#' Diagonal Variance Propagation Through a Transform Chain
#'
#' Convenience wrapper to propagate diagonal variances through a sequence of
#' transforms without materializing dense covariances.
#'
#' This function tracks **diagonal variances only**. When chaining multiple
#' transforms, it applies `diag_cov_orthogonal()` stepwise, discarding any
#' off-diagonal covariances induced by intermediate transforms. This is a useful
#' approximation for uncertainty summaries but is not, in general, equivalent to
#' computing the diagonal of the fully composed covariance.
#'
#' @param variances Numeric vector of length `d` giving diagonal variances in the
#'   source basis.
#' @param transform A transform specification:
#'   \itemize{
#'     \item A single square matrix/Matrix `Q`.
#'     \item A list of square matrix/Matrix operators, applied in order
#'       (first element applied first).
#'     \item An [AlignmentModel] (or [AlignmentResult]) containing per-subject
#'       operator transforms.
#'     \item A list of [AlignmentModel] (or [AlignmentResult]) objects,
#'       interpreted as a transform chain to apply in order (first model applied
#'       first). Subjects are matched by ID across models.
#'   }
#' @param convention Convention for how `Q` is applied at each step:
#'   `"left"` computes `diag(Q %*% diag(v) %*% t(Q))`;
#'   `"right"` computes `diag(t(Q) %*% diag(v) %*% Q)`.
#' @param subject Optional subject id(s)/index when `transform` is an
#'   [AlignmentModel] or list of models. If `NULL`, returns results for all
#'   subjects (or all subjects common to all models in the chain).
#'
#' @return If `transform` is a single matrix or list of matrices, returns a
#'   numeric vector. If `transform` is a model (or list of models) and
#'   `subject=NULL`, returns a `d x n_subjects` matrix with one column per
#'   subject; if `subject` selects a single subject, returns a numeric vector.
#'
#' @examples
#' v <- c(a = 1, b = 2)
#' Q <- matrix(c(0, 1, -1, 0), 2, 2)
#' diag_cov_orthogonal_chain(v, list(Q, Q))
#'
#' @export
diag_cov_orthogonal_chain <- function(variances,
                                      transform,
                                      convention = c("left", "right"),
                                      subject = NULL) {
  convention <- match.arg(convention)

  if (!is.numeric(variances) || !is.null(dim(variances))) {
    stop("'variances' must be a numeric vector", call. = FALSE)
  }
  if (any(!is.finite(variances))) stop("'variances' must be finite", call. = FALSE)
  if (any(variances < 0)) stop("'variances' must be non-negative", call. = FALSE)
  nm <- names(variances)

  apply_one <- function(v, Q) {
    diag_cov_orthogonal(v, Q, convention = convention)
  }

  apply_chain_Q <- function(v, Qs) {
    out <- v
    for (Qi in Qs) {
      out <- apply_one(out, Qi)
    }
    out
  }

  apply_model_subject <- function(v, model, subj) {
    Q <- model@transforms[[subj]]
    apply_one(v, Q)
  }

  if (.is_matrixish(transform)) {
    return(apply_one(variances, transform))
  }

  if (is.list(transform) && length(transform) > 0 && all(vapply(transform, .is_matrixish, logical(1)))) {
    return(apply_chain_Q(variances, transform))
  }

  if (inherits(transform, c("AlignmentModel", "AlignmentResult"))) {
    model <- .ensure_model(transform, what = "transform")
    subj_names <- names(model@transforms)
    subj_sel <- if (is.null(subject)) subj_names else .resolve_subject_subset(subject, subj_names, what = "subject")
    if (length(subj_sel) == 0L) {
      out <- matrix(numeric(0), nrow = length(variances), ncol = 0)
      rownames(out) <- nm
      return(out)
    }
    if (length(subj_sel) == 1L) {
      return(apply_model_subject(variances, model, subj_sel))
    }
    out <- vapply(
      subj_sel,
      function(s) unname(apply_model_subject(variances, model, s)),
      numeric(length(variances))
    )
    if (!is.null(nm)) rownames(out) <- nm
    colnames(out) <- subj_sel
    return(out)
  } else if (is.list(transform) && length(transform) > 0 &&
             all(vapply(transform, function(x) inherits(x, c("AlignmentModel", "AlignmentResult")), logical(1)))) {
    models <- lapply(transform, .ensure_model, what = "transform")
    subj_common <- Reduce(intersect, lapply(models, function(m) names(m@transforms)))
    if (length(subj_common) == 0L) {
      stop("Transform chain has no subjects in common", call. = FALSE)
    }
    subj_sel <- if (is.null(subject)) subj_common else .resolve_subject_subset(subject, subj_common, what = "subject")
    if (length(subj_sel) == 0L) {
      out <- matrix(numeric(0), nrow = length(variances), ncol = 0)
      rownames(out) <- nm
      return(out)
    }

    chain_one_subject <- function(subj) {
      v <- variances
      for (m in models) {
        v <- apply_model_subject(v, m, subj)
      }
      v
    }

    if (length(subj_sel) == 1L) {
      return(chain_one_subject(subj_sel))
    }

    out <- vapply(subj_sel, function(s) unname(chain_one_subject(s)), numeric(length(variances)))
    if (!is.null(nm)) rownames(out) <- nm
    colnames(out) <- subj_sel
    return(out)
  }

  stop(
    "'transform' must be a matrix/Matrix, a list of matrices, an AlignmentModel/AlignmentResult, or a list of models",
    call. = FALSE
  )
}
