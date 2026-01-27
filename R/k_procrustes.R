#' K-metric Procrustes Utilities
#'
#' Utilities for K-metric (kernel-weighted) orthonormalization and Procrustes
#' alignment. These are general-purpose helpers for aligners that operate in a
#' metric defined by a positive (semi-)definite kernel `K`.
#'
#' @name k_procrustes
NULL

.k_assert_square <- function(K, name = "K") {
  if (!.is_matrixish(K)) {
    stop(sprintf("'%s' must be matrix-like", name), call. = FALSE)
  }
  K <- as.matrix(K)
  if (nrow(K) != ncol(K)) {
    stop(sprintf("'%s' must be square", name), call. = FALSE)
  }
  K
}

#' Kernel Roots (K^(1/2) and K^(-1/2))
#'
#' Compute numerically stable square-root and inverse square-root matrices for a
#' (symmetric) kernel `K`, using a small eigenvalue floor (`jitter`).
#'
#' @param K Square kernel matrix (PSD recommended).
#' @param jitter Non-negative eigenvalue floor (default `1e-10`).
#'
#' @return A list with elements:
#' \describe{
#'   \item{Khalf}{Matrix `K^{1/2}`.}
#'   \item{Kihalf}{Matrix `K^{-1/2}`.}
#'   \item{evals}{Eigenvalues after flooring.}
#' }
#'
#' @examples
#' K <- diag(3)
#' roots <- k_kernel_roots(K)
#' roots$Khalf
#' roots$Kihalf
#'
#' @family k_procrustes
#' @export
k_kernel_roots <- function(K, jitter = 1e-10) {
  K <- .k_assert_square(K, name = "K")
  jitter <- as.numeric(jitter)
  if (!is.finite(jitter) || jitter < 0) {
    stop("'jitter' must be a non-negative numeric scalar", call. = FALSE)
  }

  Ksym <- (K + t(K)) / 2
  ee <- eigen(Ksym, symmetric = TRUE)
  vals <- pmax(ee$values, jitter)
  V <- ee$vectors
  Khalf <- V %*% diag(sqrt(vals), length(vals)) %*% t(V)
  Kihalf <- V %*% diag(1 / sqrt(vals), length(vals)) %*% t(V)
  list(Khalf = Khalf, Kihalf = Kihalf, evals = vals)
}

#' K-orthonormalize Columns
#'
#' Produce a matrix `U` with K-orthonormal columns: `t(U) %*% K %*% U = I`.
#' Uses Cholesky factorization of the Gram matrix `t(W) K W` to preserve column
#' order and sign determinism.
#' Falls back to Euclidean QR on `K^{1/2} W` (mapped back with `K^{-1/2}`)
#' when the Gram matrix is not positive definite.
#'
#' @param W Matrix whose columns span the subspace to orthonormalize.
#' @param K Square kernel matrix defining the metric.
#' @param Kroots Optional precomputed kernel roots from `k_kernel_roots()`.
#'
#' @return Matrix `U` with the same number of columns as `W`.
#'
#' @examples
#' set.seed(1)
#' K <- diag(5)
#' W <- matrix(rnorm(10), 5, 2)
#' U <- k_orthonormalize(W, K)
#' round(t(U) %*% K %*% U, 10)
#'
#' @family k_procrustes
#' @export
k_orthonormalize <- function(W, K, Kroots = NULL) {
  if (!.is_matrixish(W)) {
    stop("'W' must be matrix-like", call. = FALSE)
  }
  W <- as.matrix(W)
  K <- .k_assert_square(K, name = "K")
  if (nrow(K) != nrow(W)) {
    stop("Row dimension mismatch: nrow(K) must equal nrow(W)", call. = FALSE)
  }

  r <- ncol(W)
  KW <- K %*% W
  G <- crossprod(W, KW)            # t(W) K W, r x r
  G <- (G + t(G)) / 2              # enforce exact symmetry

  R_chol <- tryCatch(chol(G), error = function(e) NULL)
  if (!is.null(R_chol)) {
    U <- W %*% backsolve(R_chol, diag(r))
    # Verify K-orthonormality; fall back to QR if poor
    check <- crossprod(U, K %*% U)
    if (max(abs(check - diag(r))) < 1e-8) {
      return(U)
    }
  }

  # Fallback: QR in transformed space (handles rank-deficient W)
  Kroots <- Kroots %||% k_kernel_roots(K)
  qrw <- qr(Kroots$Khalf %*% W)
  Q <- qr.Q(qrw)
  R <- qr.R(qrw)
  sgn <- sign(diag(R))
  sgn[!is.finite(sgn) | sgn == 0] <- 1
  Q <- sweep(Q, 2, sgn, `*`)
  Kroots$Kihalf %*% Q
}

#' K-metric Procrustes Alignment
#'
#' Align basis `U` to reference `Uref` in a kernel-weighted metric by solving
#' an orthogonal Procrustes problem:
#' \deqn{U R \approx U_{ref}}
#' where orthogonality is Euclidean in the reduced (column) space and the match
#' is measured in the `K`-metric.
#'
#' @param Uref Reference basis (q x r).
#' @param U Basis to align (q x r).
#' @param K Square kernel matrix (q x q).
#' @param allow_reflection Logical; if `FALSE`, forces `det(R) = +1`.
#' @param reflection Alias for `allow_reflection`.
#'
#' @return A list with elements:
#' \describe{
#'   \item{U_aligned}{Aligned basis `U %*% R`.}
#'   \item{R}{Orthogonal r x r rotation/reflection matrix.}
#'   \item{d}{Sum of singular values of the cross-gram matrix.}
#'   \item{cosines}{Principal cosines (clipped to `[0,1]`).}
#'   \item{determinant}{`det(R)`.}
#' }
#'
#' @examples
#' set.seed(1)
#' q <- 6; r <- 3
#' A <- matrix(rnorm(q * q), q, q)
#' K <- t(A) %*% A + diag(0.1, q)
#' Uref <- k_orthonormalize(matrix(rnorm(q * r), q, r), K)
#' Rtrue <- qr.Q(qr(matrix(rnorm(r * r), r, r)))
#' U <- Uref %*% Rtrue
#' pr <- k_procrustes(Uref, U, K)
#' norm(pr$U_aligned - Uref, "F")
#'
#' @family k_procrustes
#' @export
k_procrustes <- function(Uref, U, K, allow_reflection = TRUE, reflection = NULL) {
  if (!is.null(reflection)) {
    if (!missing(allow_reflection) && isTRUE(allow_reflection) != isTRUE(reflection)) {
      stop("Provide only one of 'allow_reflection' and 'reflection' (and they must agree).", call. = FALSE)
    }
    allow_reflection <- reflection
  }

  if (!.is_matrixish(Uref) || !.is_matrixish(U)) {
    stop("'Uref' and 'U' must be matrix-like", call. = FALSE)
  }
  Uref <- as.matrix(Uref)
  U <- as.matrix(U)
  if (!identical(dim(Uref), dim(U))) {
    stop("'Uref' and 'U' must have identical dimensions", call. = FALSE)
  }

  K <- .k_assert_square(K, name = "K")
  if (nrow(K) != nrow(U)) {
    stop("Dimension mismatch: nrow(K) must equal nrow(Uref/U)", call. = FALSE)
  }

  allow_reflection <- isTRUE(allow_reflection)

  # Cross-gram for the weighted Procrustes problem: (U^T K Uref).
  C <- t(U) %*% K %*% Uref
  sv <- svd(C)
  R <- sv$u %*% t(sv$v)

  if (!allow_reflection && det(R) < 0) {
    Uu <- sv$u
    Uu[, ncol(Uu)] <- -Uu[, ncol(Uu)]
    R <- Uu %*% t(sv$v)
  }

  list(
    U_aligned = U %*% R,
    R = R,
    d = sum(sv$d),
    cosines = pmin(sv$d, 1),
    determinant = det(R)
  )
}

#' Align Multiple Bases to a Reference
#'
#' Convenience wrapper aligning a list of bases to a single reference (either a
#' matrix or an index into `U_list`).
#'
#' @param U_list List of q x r bases.
#' @param K Square kernel matrix.
#' @param ref Reference basis or index (default 1).
#' @param allow_reflection Passed to `k_procrustes()`.
#' @param weights Optional numeric weights (stored in return value).
#' @param reflection Alias for `allow_reflection`.
#'
#' @return A list with elements `U_aligned`, `R`, `Uref`, `score`, `weights`.
#'
#' @family k_procrustes
#' @export
k_align_bases <- function(U_list, K, ref = 1L,
                          allow_reflection = TRUE,
                          weights = NULL,
                          reflection = NULL) {
  if (!is.null(reflection)) {
    if (!missing(allow_reflection) && isTRUE(allow_reflection) != isTRUE(reflection)) {
      stop("Provide only one of 'allow_reflection' and 'reflection' (and they must agree).", call. = FALSE)
    }
    allow_reflection <- reflection
  }

  if (!is.list(U_list) || length(U_list) < 1L) {
    stop("'U_list' must be a non-empty list of bases", call. = FALSE)
  }

  Uref <- if (.is_matrixish(ref)) {
    as.matrix(ref)
  } else {
    U_list[[as.integer(ref)]]
  }

  aligned <- vector("list", length(U_list))
  Rlist <- vector("list", length(U_list))
  scores <- numeric(length(U_list))

  for (i in seq_along(U_list)) {
    pr <- k_procrustes(Uref, U_list[[i]], K, allow_reflection = allow_reflection)
    aligned[[i]] <- pr$U_aligned
    Rlist[[i]] <- pr$R
    scores[[i]] <- pr$d
  }

  list(
    U_aligned = aligned,
    R = Rlist,
    Uref = Uref,
    score = scores,
    weights = weights
  )
}

#' Consensus K-orthonormal Basis (K-Procrustes Mean)
#'
#' Iteratively aligns bases to a current estimate, averages (optionally
#' weighted), and retracts back to the Stiefel manifold via K-orthonormalization.
#'
#' @param U_list List of q x r K-orthonormal bases.
#' @param K Square kernel matrix.
#' @param weights Optional numeric weights (default equal).
#' @param Kroots Optional precomputed kernel roots.
#' @param max_iter Maximum iterations.
#' @param tol Convergence tolerance on `1 - mean(principal cosines)`.
#' @param allow_reflection Passed to alignment step.
#' @param reflection Alias for `allow_reflection`.
#'
#' @return A list with elements `U`, `iters`, `converged`, `gaps`, `scores`.
#'
#' @family k_procrustes
#' @export
k_consensus_basis <- function(U_list, K,
                              weights = NULL,
                              Kroots = NULL,
                              max_iter = 50,
                              tol = 1e-6,
                              allow_reflection = TRUE,
                              reflection = NULL) {
  if (!is.null(reflection)) {
    if (!missing(allow_reflection) && isTRUE(allow_reflection) != isTRUE(reflection)) {
      stop("Provide only one of 'allow_reflection' and 'reflection' (and they must agree).", call. = FALSE)
    }
    allow_reflection <- reflection
  }

  if (!is.list(U_list) || length(U_list) < 1L) {
    stop("'U_list' must be a non-empty list of bases", call. = FALSE)
  }

  K <- .k_assert_square(K, name = "K")
  Kroots <- Kroots %||% k_kernel_roots(K)

  Uc <- as.matrix(U_list[[1L]])
  r <- ncol(Uc)
  S <- length(U_list)

  if (is.null(weights)) {
    weights <- rep(1, S)
  }
  weights <- as.numeric(weights)
  if (length(weights) != S || any(!is.finite(weights)) || any(weights < 0)) {
    stop("'weights' must be a non-negative numeric vector of length length(U_list)", call. = FALSE)
  }
  if (sum(weights) == 0) {
    stop("'weights' must have positive sum", call. = FALSE)
  }
  weights <- weights / sum(weights)

  gaps <- numeric(0)
  scores_hist <- numeric(0)

  for (it in seq_len(as.integer(max_iter))) {
    aligned <- k_align_bases(
      U_list = U_list,
      K = K,
      ref = Uc,
      allow_reflection = allow_reflection,
      weights = weights
    )

    Abar <- Reduce(`+`, Map(function(Ui, w) w * Ui, aligned$U_aligned, weights))
    Unew <- k_orthonormalize(Abar, K, Kroots = Kroots)

    svals <- svd(t(Uc) %*% K %*% Unew)$d
    gap <- 1 - mean(pmin(svals, 1))
    gaps <- c(gaps, gap)
    scores_hist <- c(scores_hist, mean(aligned$score))

    Uc <- Unew
    if (gap < tol) {
      return(list(
        U = Uc,
        iters = it,
        converged = TRUE,
        gaps = gaps,
        scores = scores_hist
      ))
    }
  }

  list(
    U = Uc,
    iters = as.integer(max_iter),
    converged = FALSE,
    gaps = gaps,
    scores = scores_hist
  )
}
