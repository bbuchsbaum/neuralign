#' Identified Subspace Helpers
#'
#' Domain-agnostic utilities for working in a lower-dimensional subspace of a
#' transform space (e.g., when correspondence signals are rank-deficient).
#'
#' @name subspace
NULL

.svd_rank <- function(d, tol = sqrt(.Machine$double.eps)) {
  d <- as.numeric(d)
  if (length(d) == 0) return(0L)
  if (all(d == 0)) return(0L)
  as.integer(sum(d > (max(d) * tol)))
}

.is_orthogonal_operator <- function(Q, tol = 1e-6) {
  if (!.is_matrixish(Q)) return(FALSE)
  Q <- as.matrix(Q)
  if (nrow(Q) != ncol(Q)) return(FALSE)
  if (any(!is.finite(Q))) return(FALSE)
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol <= 0) {
    stop("'tol' must be a single positive number", call. = FALSE)
  }
  # Check Q'Q ≈ I. Use max absolute deviation as a simple diagnostic.
  I <- diag(nrow(Q))
  dev <- max(abs(crossprod(Q) - I))
  is.finite(dev) && dev <= tol
}

.validate_subspace_basis <- function(basis) {
  if (!.is_matrixish(basis)) {
    stop("'basis' must be matrix-like", call. = FALSE)
  }
  basis <- as.matrix(basis)
  if (ncol(basis) < 0L) stop("'basis' must have non-negative column count", call. = FALSE)
  if (any(!is.finite(basis))) stop("'basis' must be finite", call. = FALSE)
  basis
}

#' Compute an Identified Subspace Basis
#'
#' Given a correspondence matrix that drives alignment, compute an orthonormal
#' basis for the identified subspace of the transform space.
#'
#' If `convention="right"` (transforms act on columns), the identified subspace
#' is the span of the right singular vectors corresponding to non-zero singular
#' values. If `convention="left"` (transforms act on rows), the identified
#' subspace is the span of the corresponding left singular vectors.
#'
#' @param x Numeric matrix (correspondence features x transform dimension, or
#'   vice versa depending on convention). Must be finite.
#' @param convention Convention for interpreting the transform dimension:
#'   `"right"` means transforms act on columns; `"left"` means transforms act on
#'   rows.
#' @param tol Relative tolerance for numeric rank determination (fraction of
#'   the largest singular value).
#'
#' @return An object of class `"identified_subspace"` with components:
#' \describe{
#'   \item{basis}{`d x r` orthonormal basis matrix (columns).}
#'   \item{singular_values}{Singular values of `x`.}
#'   \item{numeric_rank}{Computed numeric rank `r`.}
#'   \item{transform_dim}{Transform dimension `d`.}
#'   \item{convention}{The convention used.}
#'   \item{tol}{The tolerance used.}
#' }
#'
#' @export
identified_subspace_basis <- function(x,
                                      convention = c("right", "left"),
                                      tol = sqrt(.Machine$double.eps)) {
  convention <- match.arg(convention)
  if (!.is_matrixish(x)) stop("'x' must be matrix-like", call. = FALSE)
  x <- as.matrix(x)
  if (any(!is.finite(x))) stop("'x' must be finite", call. = FALSE)
  if (!is.numeric(x)) storage.mode(x) <- "double"

  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol <= 0 || tol >= 1) {
    stop("'tol' must be a single number in (0, 1)", call. = FALSE)
  }

  d <- if (convention == "left") nrow(x) else ncol(x)
  # Use an economy SVD; we only need singular vectors in transform space.
  sv <- if (convention == "left") {
    svd(x, nu = min(nrow(x), ncol(x)), nv = 0)
  } else {
    svd(x, nu = 0, nv = min(nrow(x), ncol(x)))
  }

  r <- .svd_rank(sv$d, tol = tol)
  basis <- if (r == 0L) {
    matrix(0, nrow = d, ncol = 0L)
  } else if (convention == "left") {
    sv$u[, seq_len(r), drop = FALSE]
  } else {
    sv$v[, seq_len(r), drop = FALSE]
  }

  out <- list(
    basis = basis,
    singular_values = sv$d,
    numeric_rank = as.integer(r),
    transform_dim = as.integer(d),
    convention = convention,
    tol = tol
  )
  class(out) <- "identified_subspace"
  out
}

#' Project Data Onto a Subspace
#'
#' @param x Numeric vector or matrix.
#' @param basis Orthonormal basis matrix with `d` rows and `r` columns.
#' @param side If `x` is a matrix, which side carries the `d`-dimensional space:
#'   `"left"` means `nrow(x) == d` and projects rows; `"right"` means
#'   `ncol(x) == d` and projects columns.
#'
#' @return `x` projected onto the column span of `basis`, with the same shape as
#'   the input.
#'
#' @export
project_to_subspace <- function(x, basis, side = c("left", "right")) {
  side <- match.arg(side)
  basis <- .validate_subspace_basis(basis)

  if (!.is_matrixish(x) && !is.numeric(x)) {
    stop("'x' must be numeric or matrix-like", call. = FALSE)
  }

  if (is.null(dim(x))) {
    x <- as.numeric(x)
    if (length(x) != nrow(basis)) {
      stop("Dimension mismatch: length(x) must equal nrow(basis)", call. = FALSE)
    }
    out <- basis %*% (crossprod(basis, x))
    out <- as.numeric(out)
    names(out) <- names(x)
    return(out)
  }

  x <- as.matrix(x)
  if (any(!is.finite(x))) stop("'x' must be finite", call. = FALSE)

  if (side == "left") {
    if (nrow(x) != nrow(basis)) stop("Dimension mismatch: nrow(x) must equal nrow(basis)", call. = FALSE)
    out <- basis %*% (crossprod(basis, x))
  } else {
    if (ncol(x) != nrow(basis)) stop("Dimension mismatch: ncol(x) must equal nrow(basis)", call. = FALSE)
    out <- (x %*% basis) %*% t(basis)
  }

  dimnames(out) <- dimnames(x)
  out
}

#' Subspace Coordinates
#'
#' Convert vectors/matrices in the ambient space into coordinates in a subspace
#' basis.
#'
#' @param x Numeric vector or matrix.
#' @param basis Orthonormal basis matrix with `d` rows and `r` columns.
#' @param side If `x` is a matrix, which side carries the `d`-dimensional space:
#'   `"left"` means `nrow(x) == d` and returns an r-by-`ncol(x)` matrix;
#'   `"right"` means `ncol(x) == d` and returns an `nrow(x)`-by-r matrix.
#'
#' @return Coordinates of `x` in the subspace.
#'
#' @export
subspace_coordinates <- function(x, basis, side = c("left", "right")) {
  side <- match.arg(side)
  basis <- .validate_subspace_basis(basis)

  if (!.is_matrixish(x) && !is.numeric(x)) {
    stop("'x' must be numeric or matrix-like", call. = FALSE)
  }

  if (is.null(dim(x))) {
    x <- as.numeric(x)
    if (length(x) != nrow(basis)) {
      stop("Dimension mismatch: length(x) must equal nrow(basis)", call. = FALSE)
    }
    return(as.numeric(crossprod(basis, x)))
  }

  x <- as.matrix(x)
  if (any(!is.finite(x))) stop("'x' must be finite", call. = FALSE)

  if (side == "left") {
    if (nrow(x) != nrow(basis)) stop("Dimension mismatch: nrow(x) must equal nrow(basis)", call. = FALSE)
    crossprod(basis, x)
  } else {
    if (ncol(x) != nrow(basis)) stop("Dimension mismatch: ncol(x) must equal nrow(basis)", call. = FALSE)
    x %*% basis
  }
}

#' Restrict an Operator to a Subspace
#'
#' Compute the r-by-r operator representing `Q` in the coordinates of `basis`:
#' \deqn{Q_r = B^\top Q B.}
#'
#' For any `z` in `R^r`, defining `x = B z` (a vector in the subspace), we have:
#' \deqn{B^\top (Q x) = Q_r z.}
#'
#' @param Q Square `d x d` operator.
#' @param basis Orthonormal basis matrix `d x r`.
#'
#' @return r-by-r matrix.
#'
#' @export
restrict_operator_to_subspace <- function(Q, basis) {
  basis <- .validate_subspace_basis(basis)
  if (!.is_matrixish(Q)) stop("'Q' must be matrix-like", call. = FALSE)
  Q <- as.matrix(Q)
  if (nrow(Q) != ncol(Q)) stop("'Q' must be square", call. = FALSE)
  if (nrow(Q) != nrow(basis)) stop("Dimension mismatch: nrow(Q) must equal nrow(basis)", call. = FALSE)
  if (any(!is.finite(Q))) stop("'Q' must be finite", call. = FALSE)
  crossprod(basis, Q %*% basis)
}

#' Lift a Subspace Operator to Full Space
#'
#' @param Q_sub r-by-r operator in subspace coordinates.
#' @param basis Orthonormal basis matrix `d x r`.
#' @param fill How to act on the orthogonal complement: `"zero"` leaves it as
#'   zero; `"identity"` uses identity on the complement.
#'
#' @return `d x d` matrix operator.
#'
#' @export
lift_operator_from_subspace <- function(Q_sub, basis, fill = c("zero", "identity")) {
  fill <- match.arg(fill)
  basis <- .validate_subspace_basis(basis)
  if (!.is_matrixish(Q_sub)) stop("'Q_sub' must be matrix-like", call. = FALSE)
  Q_sub <- as.matrix(Q_sub)
  r <- ncol(basis)
  if (!identical(dim(Q_sub), c(r, r))) {
    stop("Dimension mismatch: Q_sub must be r x r where r=ncol(basis)", call. = FALSE)
  }
  if (any(!is.finite(Q_sub))) stop("'Q_sub' must be finite", call. = FALSE)

  full <- basis %*% Q_sub %*% t(basis)
  if (fill == "identity") {
    P <- basis %*% t(basis)
    full <- full + (diag(nrow(P)) - P)
  }
  full
}

#' Canonicalize an Orthogonal Operator in the Identified Subspace
#'
#' When a correspondence matrix is rank-deficient, orthogonal operators are
#' under-identified: their action on the orthogonal complement of the identified
#' subspace does not affect the fit. This helper replaces that arbitrary part
#' with a deterministic extension, while preserving the action of `Q` on the
#' identified subspace of `x`.
#'
#' In particular, if `S` is the identified subspace of the transform space
#' implied by `x`, the returned operator `Q_canon` satisfies:
#' \deqn{Q_{\mathrm{canon}} v = Q v \quad \forall v \in S.}
#'
#' @param Q Square orthogonal operator (matrix or Matrix).
#' @param x Correspondence matrix used to determine the identified subspace.
#' @param convention Convention for interpreting the transform dimension:
#'   `"left"` means transforms act on rows (`d = nrow(x)`); `"right"` means
#'   transforms act on columns (`d = ncol(x)`).
#' @param tol Relative tolerance for numeric rank determination (fraction of the
#'   largest singular value).
#' @param check_orthogonal Logical; if TRUE (default), check that `Q` is
#'   orthogonal (within `orthogonal_tol`) before canonicalizing.
#' @param orthogonal_tol Tolerance for the orthogonality check.
#'
#' @return A `d x d` orthogonal matrix.
#' @export
canonicalize_orthogonal_operator <- function(Q,
                                             x,
                                             convention = c("left", "right"),
                                             tol = sqrt(.Machine$double.eps),
                                             check_orthogonal = TRUE,
                                             orthogonal_tol = 1e-6) {
  convention <- match.arg(convention)

  if (!.is_matrixish(Q)) stop("'Q' must be matrix-like", call. = FALSE)
  Q <- as.matrix(Q)
  if (nrow(Q) != ncol(Q)) stop("'Q' must be square", call. = FALSE)
  if (any(!is.finite(Q))) stop("'Q' must be finite", call. = FALSE)

  if (isTRUE(check_orthogonal) && !.is_orthogonal_operator(Q, tol = orthogonal_tol)) {
    stop("'Q' must be orthogonal (within 'orthogonal_tol')", call. = FALSE)
  }

  subspace <- identified_subspace_basis(x, convention = convention, tol = tol)
  B <- subspace$basis
  d <- subspace$transform_dim

  if (nrow(Q) != d) {
    stop(
      sprintf("Dimension mismatch: nrow(Q)=%d but transform_dim=%d", nrow(Q), d),
      call. = FALSE
    )
  }

  # Compute deterministic full bases for the identified subspace and its image.
  # Note: .complete_orthonormal_basis() lives in the package namespace.
  B_full <- .complete_orthonormal_basis(B)
  B_img <- Q %*% B
  B_img_full <- .complete_orthonormal_basis(B_img)

  # Canonical orthogonal operator that matches Q on span(B).
  B_img_full %*% t(B_full)
}
