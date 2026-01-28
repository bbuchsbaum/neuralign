#' Intrinsic Geometry Guidance (Laplacian Eigenmodes)
#'
#' Compute template-free intrinsic geometry coordinates (eigenmodes) from an
#' adjacency graph and store them as `"coords"` guidance channels for alignment.
#'
#' These coordinates are *guidance* only: they provide a smooth positional
#' encoding that methods can use as priors/costs, not the alignment target.
#'
#' @name intrinsic_geometry
#' @family guidance
NULL

.fix_eigenmode_signs <- function(Z) {
  Z <- as.matrix(Z)
  for (k in seq_len(ncol(Z))) {
    v <- Z[, k]
    idx <- which(abs(v) > 0)[1]
    if (!is.na(idx) && v[[idx]] < 0) {
      Z[, k] <- -v
    }
  }
  Z
}

.as_symmetric_adjacency <- function(adj) {
  if (!.is_matrixish(adj)) {
    stop("'adj' must be a matrix-like adjacency", call. = FALSE)
  }
  A <- if (inherits(adj, "Matrix")) adj else as.matrix(adj)
  if (nrow(A) != ncol(A)) {
    stop("'adj' must be square", call. = FALSE)
  }

  # Basic symmetry check (cheap): compare a few entries, then full if dense
  if (!inherits(A, "Matrix")) {
    if (!isTRUE(all.equal(A, t(A)))) {
      stop("'adj' must be symmetric", call. = FALSE)
    }
  } else {
    if (!Matrix::isSymmetric(A)) {
      stop("'adj' must be symmetric", call. = FALSE)
    }
  }

  A
}

#' Laplacian Eigenmodes
#'
#' Compute intrinsic geometry coordinates from a symmetric adjacency graph via
#' Laplacian eigenmodes.
#'
#' @family guidance
#'
#' @param adj Symmetric adjacency matrix (dense matrix or sparse `Matrix`).
#' @param k Number of non-trivial eigenmodes to return.
#' @param normalized Logical; if `TRUE`, use the symmetric normalized Laplacian.
#' @param drop_zero Logical; if `TRUE` (default), drop near-zero eigenmodes
#'   (constant components) and return the next `k`.
#' @param tol Tolerance for considering an eigenvalue "zero".
#' @param method Eigen solver: `"auto"` uses `RSpectra` if available (for sparse),
#'   otherwise base `eigen()`.
#' @param return_values Logical; if `TRUE`, return a list with eigenvalues.
#'
#' @return Matrix of size `(n_nodes x k)` or a list with `$vectors` and `$values`.
#' @examples
#' # Path graph with 6 nodes
#' n <- 6
#' A <- matrix(0, n, n)
#' A[cbind(1:(n - 1), 2:n)] <- 1
#' A <- A + t(A)
#'
#' Z <- laplacian_eigenmodes(A, k = 2)
#' dim(Z)
#' @export
laplacian_eigenmodes <- function(adj,
                                 k = 10,
                                 normalized = TRUE,
                                 drop_zero = TRUE,
                                 tol = 1e-10,
                                 method = c("auto", "eigen", "rspectra"),
                                 return_values = FALSE) {
  method <- match.arg(method)
  k <- as.integer(k)
  if (!is.finite(k) || k < 1L) {
    stop("'k' must be a positive integer", call. = FALSE)
  }
  tol <- as.numeric(tol)
  if (!is.finite(tol) || tol < 0) {
    stop("'tol' must be a non-negative number", call. = FALSE)
  }

  A <- .as_symmetric_adjacency(adj)
  n <- nrow(A)

  d <- Matrix::rowSums(A)
  if (any(!is.finite(d)) || any(d < 0)) {
    stop("Invalid adjacency: non-finite or negative degrees", call. = FALSE)
  }

  if (isTRUE(normalized)) {
    inv_sqrt <- rep(0, n)
    ok <- d > 0
    inv_sqrt[ok] <- 1 / sqrt(d[ok])
    D_inv_sqrt <- Matrix::Diagonal(x = inv_sqrt)
    L <- Matrix::Diagonal(n) - (D_inv_sqrt %*% A %*% D_inv_sqrt)
  } else {
    L <- Matrix::Diagonal(x = d) - A
  }

  # Determine solver
  use_rspectra <- FALSE
  if (method == "rspectra") {
    use_rspectra <- TRUE
  } else if (method == "auto") {
    use_rspectra <- inherits(L, "Matrix") && requireNamespace("RSpectra", quietly = TRUE)
  }

  if (use_rspectra) {
    if (!requireNamespace("RSpectra", quietly = TRUE)) {
      stop("method='rspectra' requires the RSpectra package", call. = FALSE)
    }

    # Request extra modes to allow dropping near-zero components.
    extra_modes <- 5L
    k_req <- min(n, k + extra_modes)
    eig <- RSpectra::eigs_sym(L, k = k_req, which = "SM")
    vals <- eig$values
    vecs <- eig$vectors
    ord <- order(vals, decreasing = FALSE)
    vals <- vals[ord]
    vecs <- vecs[, ord, drop = FALSE]
  } else {
    eigen_warn_threshold <- 4000L
    if (n > eigen_warn_threshold) {
      warning(
        sprintf(
          "Computing %dx%d eigen-decomposition with base eigen(); consider installing RSpectra for scalability",
          n, n
        ),
        call. = FALSE
      )
    }
    eig <- eigen(as.matrix(L), symmetric = TRUE)
    vals <- eig$values
    vecs <- eig$vectors
    ord <- order(vals, decreasing = FALSE)
    vals <- vals[ord]
    vecs <- vecs[, ord, drop = FALSE]
  }

  start <- 1L
  if (isTRUE(drop_zero)) {
    start <- which(vals > tol)[1]
    if (is.na(start)) {
      stop("No non-trivial eigenmodes found (graph may be too disconnected)", call. = FALSE)
    }
  }

  end <- start + k - 1L
  if (end > ncol(vecs)) {
    stop(
      sprintf("Requested k=%d modes but only %d available", k, ncol(vecs) - start + 1L),
      call. = FALSE
    )
  }

  Z <- vecs[, start:end, drop = FALSE]
  Z <- .fix_eigenmode_signs(Z)

  if (isTRUE(return_values)) {
    return(list(vectors = Z, values = vals[start:end]))
  }
  Z
}

#' Add Intrinsic Geometry Guidance to AlignmentData
#'
#' Compute Laplacian eigenmodes per subject and store them as `"coords"` guidance
#' channels on an `AlignmentData` object.
#'
#' @param data An `AlignmentData` object.
#' @param adjacency_by_subject Named list of symmetric adjacency matrices, one
#'   per subject.
#' @param channel_name Guidance channel name (default `"intrinsic"`).
#' @param k Number of eigenmodes to compute.
#' @param normalized Logical; use normalized Laplacian.
#' @param drop_zero Logical; drop near-zero eigenmodes.
#' @param tol Tolerance for zero eigenvalues.
#' @param method Eigen solver (see `laplacian_eigenmodes()`).
#' @param validate Logical; if `TRUE`, validate dimensional compatibility.
#'
#' @return The modified `AlignmentData`.
#'
#' @examples
#' set.seed(1)
#' # Line-graph adjacency for 5 nodes
#' A <- matrix(0, 5, 5)
#' A[cbind(1:4, 2:5)] <- 1
#' A[cbind(2:5, 1:4)] <- 1
#' adat <- AlignmentData(list(
#'   s1 = matrix(rnorm(20), 5, 4),
#'   s2 = matrix(rnorm(20), 5, 4)
#' ))
#' adat2 <- set_intrinsic_geometry_guidance(adat, list(s1 = A, s2 = A), k = 2, method = "eigen")
#' Z1 <- get_guidance(adat2, subject = "s1")$intrinsic$value
#' dim(Z1)
#' @export
set_intrinsic_geometry_guidance <- function(data,
                                            adjacency_by_subject,
                                            channel_name = "intrinsic",
                                            k = 10,
                                            normalized = TRUE,
                                            drop_zero = TRUE,
                                            tol = 1e-10,
                                            method = c("auto", "eigen", "rspectra"),
                                            validate = TRUE) {
  if (!inherits(data, "AlignmentData")) {
    stop("'data' must be an AlignmentData object", call. = FALSE)
  }
  if (!is.list(adjacency_by_subject) || is.null(names(adjacency_by_subject))) {
    stop("'adjacency_by_subject' must be a named list keyed by subject", call. = FALSE)
  }
  if (!is.character(channel_name) || length(channel_name) != 1L || !nzchar(channel_name)) {
    stop("'channel_name' must be a non-empty character scalar", call. = FALSE)
  }
  method <- match.arg(method)

  subjects <- data@subjects
  missing <- setdiff(subjects, names(adjacency_by_subject))
  if (length(missing) > 0) {
    stop(
      sprintf("adjacency_by_subject is missing subjects: %s", paste(missing, collapse = ", ")),
      call. = FALSE
    )
  }

  coords <- lapply(subjects, function(subj) {
    A <- adjacency_by_subject[[subj]]
    laplacian_eigenmodes(
      adj = A,
      k = k,
      normalized = normalized,
      drop_zero = drop_zero,
      tol = tol,
      method = method,
      return_values = FALSE
    )
  })
  names(coords) <- subjects

  g <- get_guidance(data)
  for (subj in subjects) {
    g_sub <- g[[subj]] %||% list()
    Z <- coords[[subj]]
    g_sub[[channel_name]] <- guidance_channel(
      type = "coords",
      value = Z,
      name = channel_name,
      k = ncol(Z),
      normalized = normalized
    )
    g[[subj]] <- g_sub
  }

  set_guidance(data, g, validate = validate)
}
