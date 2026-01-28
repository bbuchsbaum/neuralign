#' Procrustes Alignment
#'
#' Procrustes alignment via manifoldalign. Finds orthogonal transformations
#' that minimize the sum of squared distances between subjects' data.
#'
#' @name aligner_procrustes
NULL


#' Procrustes Fit Function
#'
#' Internal fit function for Procrustes alignment.
#'
#' @param data AlignmentData object.
#' @param reference Reference specification.
#' @param train_idx Indices of subjects to use for fitting.
#' @param scale Logical; if TRUE, also optimize scaling.
#' @param reflection Logical; if TRUE, allow reflections.
#' @param tol Convergence tolerance.
#' @param max_iter Maximum iterations for GPA.
#' @param ... Additional arguments.
#'
#' @return List with transforms, reference_data, etc.
#'
#' @keywords internal
.procrustes_fit <- function(data,
                            reference = "consensus",
                            train_idx = NULL,
                            scale = FALSE,
                            reflection = FALSE,
                            rank = NULL,
                            tol = 1e-6,
                            max_iter = 100,
                            allow_reflection = NULL,
                            ...) {
  if (!is.null(allow_reflection)) {
    if (!missing(reflection) && isTRUE(reflection) != isTRUE(allow_reflection)) {
      stop("Provide only one of 'reflection' and 'allow_reflection' (and they must agree).", call. = FALSE)
    }
    reflection <- isTRUE(allow_reflection)
  }

  pre <- .aligner_preamble(data, train_idx = train_idx)
  train_data <- pre$train_data
  data_list <- pre$data_list

  # Use built-in implementation. (manifoldalign's GPA API has changed over time,
  # and neuralign avoids relying on unexported or unstable symbols.)
  result <- .procrustes_fit_builtin(
    data_list, reference, train_data,
    scale = scale, reflection = reflection, rank = rank, tol = tol, max_iter = max_iter
  )

  # Ensure all original subjects have transforms
  all_subjects <- data@subjects
  transforms <- result$transforms

  for (subj in all_subjects) {
    if (!subj %in% names(transforms)) {
      # Fit transform for held-out subject
      subj_data <- get_subject_data(data, subj)
      ref_data <- result$reference_data

      # .procrustes_single returns left-multiply directly
      Q <- .procrustes_single(subj_data, ref_data,
        scale = scale,
        reflection = reflection,
        rank = rank
      )
      attr(Q, "scale_factor") <- NULL
      transforms[[subj]] <- Q
    }
  }

  list(
    transforms = transforms,
    reference_data = result$reference_data,
    space_from = train_data@space,
    space_to = train_data@space,
    method_state = list(
      scale = scale,
      reflection = reflection,
      rank = rank,
      consensus = result$reference_data
    )
  )
}


#' Built-in Procrustes Implementation
#' @keywords internal
.procrustes_fit_builtin <- function(data_list, reference, train_data,
                                    scale, reflection, rank, tol, max_iter) {
  if (is.character(reference) && reference == "consensus") {
    # Generalized Procrustes Analysis
    result <- .gpa_builtin(data_list, scale, reflection, rank, tol, max_iter)
    return(result)
  }

  # Fixed reference
  if (.is_matrixish(reference)) {
    reference_data <- as.matrix(reference)
  } else {
    reference_data <- get_subject_data(train_data, reference)
  }

  # Align each subject to reference
  # .procrustes_single returns left-multiply transform directly
  transforms <- lapply(names(data_list), function(subj) {
    if (is.character(reference) && subj == reference) {
      diag(nrow(data_list[[subj]]))
    } else {
      Q <- .procrustes_single(data_list[[subj]], reference_data,
        scale = scale,
        reflection = reflection,
        rank = rank
      )
      attr(Q, "scale_factor") <- NULL
      Q
    }
  })
  names(transforms) <- names(data_list)

  list(
    transforms = transforms,
    reference_data = reference_data
  )
}


#' Single Procrustes Alignment
#' @keywords internal
.procrustes_single <- function(X, Y, scale = FALSE, reflection = FALSE, rank = NULL) {
  # X: source data (n_features x n_obs)
  # Y: target/reference data (n_features x n_obs)
  # Find Q (n_features x n_features) such that Q %*% X ≈ Y (left-multiply)
  #
  # This is equivalent to finding Q such that t(X) %*% t(Q) ≈ t(Y)
  # Or: t(X) %*% R ≈ t(Y) where R = t(Q)
  # Standard Procrustes: H = X %*% t(Y), SVD gives rotation in feature space

  if (!is.null(rank)) {
    rank <- as.integer(rank)
    if (length(rank) != 1L || is.na(rank) || rank < 1L) {
      stop("'rank' must be a positive integer", call. = FALSE)
    }
    rmax <- min(nrow(X), ncol(X))
    if (rank > rmax) {
      stop(
        sprintf("'rank' must be <= min(n_features, n_obs) = %d", as.integer(rmax)),
        call. = FALSE
      )
    }
  }

  if (is.null(rank) || rank >= nrow(X)) {
    # Full SVD of cross-covariance in feature space (existing behavior)
    H <- X %*% t(Y)
    sv <- svd(H)
    U <- sv$u
    V <- sv$v
  } else {
    sv <- .procrustes_truncated_svd(X, Y, k = rank)
    U <- sv$u
    V <- sv$v
  }

  # Optimal rotation Q (left-multiply convention)
  Q <- V %*% t(U)

  # Handle reflection
  if (!reflection && det(Q) < 0) {
    U[, ncol(U)] <- -U[, ncol(U)]
    Q <- V %*% t(U)
  }

  if (isTRUE(scale)) {
    # Optimal scale factor s minimizing || s Q X - Y ||_F^2
    # s = <Y, QX> / ||X||_F^2
    QX <- Q %*% X
    denom <- sum(X * X)
    if (is.finite(denom) && denom > 0) {
      s <- sum(Y * QX) / denom
    } else {
      s <- 1
    }
    Q <- s * Q
    attr(Q, "scale_factor") <- s
  }

  Q
}

#' Truncated-SVD Procrustes Helper
#' @keywords internal
.procrustes_truncated_svd <- function(X, Y, k) {
  # Compute leading singular vectors of H = X %*% t(Y) without forming H,
  # then deterministically complete them to full orthonormal bases so that the
  # resulting Procrustes operator is a full orthogonal matrix.

  k <- as.integer(k)
  if (length(k) != 1L || is.na(k) || k < 1L) {
    stop("'k' must be a positive integer", call. = FALSE)
  }

  r0 <- min(nrow(X), ncol(X))
  if (k > r0) {
    stop(
      sprintf("'k' must be <= min(n_features, n_obs) = %d", as.integer(r0)),
      call. = FALSE
    )
  }

  # Economy QR: X = Qx Rx, Y = Qy Ry, with Qx/Qy having r0 columns.
  qx <- qr(X)
  Qx <- qr.Q(qx)
  Rx <- qr.R(qx)

  qy <- qr(Y)
  Qy <- qr.Q(qy)
  Ry <- qr.R(qy)

  # Small cross-covariance in the reduced space.
  C <- Rx %*% t(Ry)

  sv <- NULL
  if (k < nrow(C) && requireNamespace("RSpectra", quietly = TRUE)) {
    sv <- RSpectra::svds(C, k = k, nu = k, nv = k)
  }
  if (is.null(sv)) {
    sv <- svd(C, nu = k, nv = k)
  }

  U_k <- Qx %*% sv$u
  V_k <- Qy %*% sv$v

  list(
    u = .complete_orthonormal_basis(U_k),
    v = .complete_orthonormal_basis(V_k),
    d = sv$d
  )
}

.complete_orthonormal_basis <- function(U, tol = 1e-12) {
  if (!.is_matrixish(U)) {
    stop("'U' must be matrix-like", call. = FALSE)
  }
  U <- as.matrix(U)

  p <- nrow(U)
  k <- ncol(U)

  if (k == 0L) {
    return(diag(p))
  }
  if (k >= p) {
    return(U)
  }

  # Deterministic Gram-Schmidt completion using standard basis vectors.
  B <- matrix(0, p, p - k)
  nB <- 0L

  for (j in seq_len(p)) {
    v <- rep(0, p)
    v[[j]] <- 1

    # Remove projection onto U
    v <- v - U %*% drop(crossprod(U, v))

    # Remove projection onto existing B columns
    if (nB > 0L) {
      Bj <- B[, seq_len(nB), drop = FALSE]
      v <- v - Bj %*% drop(crossprod(Bj, v))
    }

    nv <- sqrt(sum(v * v))
    if (is.finite(nv) && nv > tol) {
      nB <- nB + 1L
      B[, nB] <- v / nv
      if (nB == (p - k)) break
    }
  }

  if (nB != (p - k)) {
    stop(
      sprintf(
        "Failed to complete orthonormal basis: needed %d columns, got %d",
        as.integer(p - k), as.integer(nB)
      ),
      call. = FALSE
    )
  }

  cbind(U, B)
}

#' Procrustes Rotation (Convention-Safe)
#'
#' Compute the (optionally scaled) Procrustes operator aligning a source matrix
#' to a target matrix, with explicit left/right multiplication conventions.
#'
#' - `convention="left"` solves for `Q` such that `Q %*% source ≈ target`
#'   where matrices are `(features x observations)`.
#' - `convention="right"` solves for `Q` such that `source %*% Q ≈ target`
#'   where matrices are `(observations x features)`.
#'
#' When per-matrix observation labels are supplied, alignment is computed on
#' the intersection of labels (in the order they appear in `obs_labels_source`).
#'
#' @family procrustes
#'
#' @param source Source matrix.
#' @param target Target/reference matrix.
#' @param convention Multiplication convention (`"left"` or `"right"`).
#' @param scale Logical; if `TRUE`, include optimal scale factor.
#' @param reflection Logical; if `TRUE`, allow reflections (det(Q) may be -1).
#' @param allow_reflection Alias for `reflection`.
#' @param obs_labels_source Optional observation labels for `source`.
#' @param obs_labels_target Optional observation labels for `target`.
#' @param min_overlap Minimum number of shared labels when labels are supplied.
#' @param rank Optional truncated-SVD rank for Procrustes fitting. If `NULL`
#'   (default), uses the full SVD. If provided, must be a positive integer
#'   `<= min(n_features, n_obs)`.
#'
#' @return A list with elements `Q`, `scale_factor`, `residual`,
#'   `convention`, and `matched_labels`.
#'
#' @examples
#' set.seed(1)
#' Q_true <- qr.Q(qr(matrix(rnorm(9), 3, 3)))
#' if (det(Q_true) < 0) Q_true[, 1] <- -Q_true[, 1]
#' X <- matrix(rnorm(15), 3, 5)
#' Y <- Q_true %*% X
#'
#' fit <- procrustes_rotation(X, Y, convention = "left")
#' max(abs(fit$Q - Q_true)) < 1e-6
#'
#' @export
procrustes_rotation <- function(source,
                                target,
                                convention = c("left", "right"),
                                scale = FALSE,
                                reflection = FALSE,
                                obs_labels_source = NULL,
                                obs_labels_target = NULL,
                                min_overlap = 2L,
                                allow_reflection = NULL,
                                rank = NULL) {
  if (!is.null(allow_reflection)) {
    if (!missing(reflection) && isTRUE(reflection) != isTRUE(allow_reflection)) {
      stop("Provide only one of 'reflection' and 'allow_reflection' (and they must agree).", call. = FALSE)
    }
    reflection <- isTRUE(allow_reflection)
  }

  convention <- match.arg(convention)
  min_overlap <- as.integer(min_overlap)
  if (!is.finite(min_overlap) || min_overlap < 1L) {
    stop("'min_overlap' must be a positive integer", call. = FALSE)
  }

  if (!.is_matrixish(source) || !.is_matrixish(target)) {
    stop("'source' and 'target' must be matrix-like", call. = FALSE)
  }

  if (convention == "left") {
    X <- as.matrix(source)
    Y <- as.matrix(target)
    if (!is.null(obs_labels_source) || !is.null(obs_labels_target)) {
      idx <- .match_obs_indices(obs_labels_source, obs_labels_target, min_overlap = min_overlap)
      X <- X[, idx$source, drop = FALSE]
      Y <- Y[, idx$target, drop = FALSE]
      matched_labels <- idx$labels
    } else {
      matched_labels <- NULL
    }
    if (ncol(X) != ncol(Y)) {
      stop(
        sprintf(
          "Left-convention requires matching observations: ncol(source)=%d, ncol(target)=%d",
          ncol(X), ncol(Y)
        ),
        call. = FALSE
      )
    }
    if (nrow(X) != nrow(Y)) {
      stop(
        sprintf(
          "Left-convention requires matching feature dimensions: nrow(source)=%d, nrow(target)=%d",
          nrow(X), nrow(Y)
        ),
        call. = FALSE
      )
    }

    Q <- .procrustes_single(X, Y,
      scale = scale,
      reflection = reflection,
      rank = rank
    )
    aligned <- Q %*% X
    resid <- norm(aligned - Y, "F")
    scale_factor <- if (isTRUE(scale)) attr(Q, "scale_factor") %||% 1 else 1
    attr(Q, "scale_factor") <- NULL

    return(list(
      Q = Q,
      scale_factor = scale_factor,
      residual = resid,
      convention = "left",
      matched_labels = matched_labels
    ))
  }

  Xr <- as.matrix(source)
  Yr <- as.matrix(target)
  if (!is.null(obs_labels_source) || !is.null(obs_labels_target)) {
    idx <- .match_obs_indices(obs_labels_source, obs_labels_target, min_overlap = min_overlap)
    Xr <- Xr[idx$source, , drop = FALSE]
    Yr <- Yr[idx$target, , drop = FALSE]
    matched_labels <- idx$labels
  } else {
    matched_labels <- NULL
  }

  if (nrow(Xr) != nrow(Yr)) {
    stop(
      sprintf(
        "Right-convention requires matching observations: nrow(source)=%d, nrow(target)=%d",
        nrow(Xr), nrow(Yr)
      ),
      call. = FALSE
    )
  }
  if (ncol(Xr) != ncol(Yr)) {
    stop(
      sprintf(
        "Right-convention requires matching feature dimensions: ncol(source)=%d, ncol(target)=%d",
        ncol(Xr), ncol(Yr)
      ),
      call. = FALSE
    )
  }

  left_res <- procrustes_rotation(
    t(Xr),
    t(Yr),
    convention = "left",
    scale = scale,
    reflection = reflection,
    obs_labels_source = NULL,
    obs_labels_target = NULL,
    min_overlap = min_overlap,
    rank = rank
  )
  Qr <- t(left_res$Q)
  aligned <- Xr %*% Qr
  resid <- norm(aligned - Yr, "F")

  list(
    Q = Qr,
    scale_factor = left_res$scale_factor,
    residual = resid,
    convention = "right",
    matched_labels = matched_labels
  )
}

#' Procrustes Distance
#'
#' Convenience wrapper returning the Procrustes residual (Frobenius norm of
#' the alignment error) between two matrices.
#'
#' @family procrustes
#'
#' @param x First matrix.
#' @param y Second matrix.
#' @param convention Multiplication convention (`"left"` or `"right"`).
#' @param obs_labels_x Optional observation labels for `x`.
#' @param obs_labels_y Optional observation labels for `y`.
#' @param min_overlap Minimum number of shared labels when labels are supplied.
#' @return Numeric scalar residual.
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(15), 3, 5)
#' Y <- X + 0.01 * matrix(rnorm(15), 3, 5)
#' procrustes_distance(X, Y)
#' @export
procrustes_distance <- function(x,
                                y,
                                convention = c("left", "right"),
                                obs_labels_x = NULL,
                                obs_labels_y = NULL,
                                min_overlap = 2L) {
  res <- procrustes_rotation(
    source = x,
    target = y,
    convention = convention,
    scale = FALSE,
    reflection = FALSE,
    obs_labels_source = obs_labels_x,
    obs_labels_target = obs_labels_y,
    min_overlap = min_overlap
  )
  res$residual
}


#' Built-in GPA Implementation
#' @keywords internal
.gpa_builtin <- function(data_list, scale, reflection, rank, tol, max_iter) {
  n_subjects <- length(data_list)
  subjects <- names(data_list)

  # Initialize with mean
  consensus <- compute_centroid(data_list)

  transforms <- vector("list", n_subjects)
  names(transforms) <- subjects

  for (iter in seq_len(max_iter)) {
    old_consensus <- consensus

    # Align each subject to current consensus
    aligned <- vector("list", n_subjects)
    for (i in seq_len(n_subjects)) {
      # Get left-multiply transform Q such that Q %*% X ≈ consensus
      Q <- .procrustes_single(data_list[[i]], consensus,
        scale = scale,
        reflection = reflection,
        rank = rank
      )
      attr(Q, "scale_factor") <- NULL
      transforms[[i]] <- Q
      # Apply left-multiply: Q %*% X
      aligned[[i]] <- Q %*% data_list[[i]]
    }

    # Update consensus
    consensus <- compute_centroid(aligned)

    # Check convergence
    diff <- norm(consensus - old_consensus, "F") / norm(old_consensus, "F")
    if (diff < tol) {
      break
    }
  }

  list(
    transforms = transforms,
    reference_data = consensus
  )
}


#' Convert List to Array
#' @keywords internal
.list_to_array <- function(data_list) {
  n <- length(data_list)
  dims <- dim(data_list[[1]])
  if (is.null(dims)) {
    dims <- c(nrow(data_list[[1]]), ncol(data_list[[1]]))
  }

  arr <- array(NA, dim = c(dims, n))
  for (i in seq_len(n)) {
    arr[, , i] <- as.matrix(data_list[[i]])
  }
  arr
}


#' Procrustes Capabilities
#' @keywords internal
.procrustes_capabilities <- list(
  supports_cv = TRUE,
  cv_axes = c("subject", "observation"),
  needs_geometry = FALSE,
  needs_design = FALSE,
  requires_shared_features = TRUE,
  requires_shared_observations = TRUE,
  returns_invertible = TRUE,
  transform_type = "orthogonal",
  mass_preserving = TRUE,
  returns = "operator",
  supports_new_subject = TRUE,
  supports_new_data = TRUE,
  reference_types = c("subject", "consensus", "template")
)


#' Register Procrustes Aligner
#'
#' Called during package load to register the Procrustes aligner.
#'
#' @keywords internal
.register_procrustes <- function() {
  register_aligner(
    name = "procrustes",
    fit_fn = .procrustes_fit,
    apply_fn = NULL,  # Uses default left-multiply
    capabilities = .procrustes_capabilities,
    package = "neuralign",
    description = "Orthogonal Procrustes via GPA",
    version = as.character(utils::packageVersion("neuralign"))
  )
}
