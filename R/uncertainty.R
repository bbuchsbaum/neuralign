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
