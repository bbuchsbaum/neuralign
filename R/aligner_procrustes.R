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
                            tol = 1e-6,
                            max_iter = 100,
                            ...) {
  # Handle train_idx
  if (is.null(train_idx)) {
    train_idx <- seq_along(data@subjects)
  }

  train_data <- data[train_idx]
  train_subjects <- train_data@subjects
  data_list <- get_data_list(train_data)

  # Check if manifoldalign is available
  use_manifoldalign <- requireNamespace("manifoldalign", quietly = TRUE)

  if (use_manifoldalign) {
    # Use manifoldalign for fitting
    result <- .procrustes_fit_manifoldalign(
      data_list, reference, train_data,
      scale = scale, tol = tol, max_iter = max_iter, ...
    )
  } else {
    # Use built-in implementation
    result <- .procrustes_fit_builtin(
      data_list, reference, train_data,
      scale = scale, reflection = reflection, tol = tol, max_iter = max_iter
    )
  }

  # Ensure all original subjects have transforms
  all_subjects <- data@subjects
  transforms <- result$transforms

  for (subj in all_subjects) {
    if (!subj %in% names(transforms)) {
      # Fit transform for held-out subject
      subj_data <- get_subject_data(data, subj)
      ref_data <- result$reference_data

      if (use_manifoldalign) {
        # manifoldalign returns right-multiply, need to transpose
        Q <- manifoldalign::procrustes(subj_data, ref_data)$Q
        transforms[[subj]] <- t(Q)
      } else {
        # .procrustes_single now returns left-multiply directly
        Q <- .procrustes_single(subj_data, ref_data, scale, reflection)
        attr(Q, "scale_factor") <- NULL
        transforms[[subj]] <- Q
      }
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
      consensus = result$reference_data
    )
  )
}


#' Procrustes Fit via manifoldalign
#' @keywords internal
.procrustes_fit_manifoldalign <- function(data_list, reference, train_data,
                                          scale, tol, max_iter, ...) {
  # Determine reference
  if (is.character(reference) && reference == "consensus") {
    # Use GPA to find consensus
    # Convert to format expected by manifoldalign
    X_array <- .list_to_array(data_list)

    gpa_result <- manifoldalign::generalized_procrustes(
      X_array,
      scale = scale,
      tol = tol,
      maxiter = max_iter
    )

    # Extract transforms (need to transpose for left-multiply)
    transforms <- lapply(seq_along(data_list), function(i) {
      t(gpa_result$rotations[[i]])  # Transpose: (target x source)
    })
    names(transforms) <- names(data_list)

    reference_data <- gpa_result$consensus
  } else {
    # Fixed reference
    if (.is_matrixish(reference)) {
      reference_data <- as.matrix(reference)
    } else {
      reference_data <- get_subject_data(train_data, reference)
    }

    # Align each subject to reference
    transforms <- lapply(names(data_list), function(subj) {
      if (is.character(reference) && subj == reference) {
        # Reference subject gets identity
        diag(nrow(data_list[[subj]]))
      } else {
        Q <- manifoldalign::procrustes(data_list[[subj]], reference_data)$Q
        t(Q)  # Transpose for left-multiply
      }
    })
    names(transforms) <- names(data_list)
  }

  list(
    transforms = transforms,
    reference_data = reference_data
  )
}


#' Built-in Procrustes Implementation
#' @keywords internal
.procrustes_fit_builtin <- function(data_list, reference, train_data,
                                    scale, reflection, tol, max_iter) {
  if (is.character(reference) && reference == "consensus") {
    # Generalized Procrustes Analysis
    result <- .gpa_builtin(data_list, scale, reflection, tol, max_iter)
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
      Q <- .procrustes_single(data_list[[subj]], reference_data, scale, reflection)
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
.procrustes_single <- function(X, Y, scale = FALSE, reflection = FALSE) {
  # X: source data (n_features x n_obs)
  # Y: target/reference data (n_features x n_obs)
  # Find Q (n_features x n_features) such that Q %*% X ≈ Y (left-multiply)
  #
  # This is equivalent to finding Q such that t(X) %*% t(Q) ≈ t(Y)
  # Or: t(X) %*% R ≈ t(Y) where R = t(Q)
  # Standard Procrustes: H = X %*% t(Y), SVD gives rotation in feature space

  # SVD of cross-covariance in feature space
  H <- X %*% t(Y)
  svd_result <- svd(H)

  # Optimal rotation Q (left-multiply convention)
  Q <- svd_result$v %*% t(svd_result$u)

  # Handle reflection
  if (!reflection && det(Q) < 0) {
    svd_result$u[, ncol(svd_result$u)] <- -svd_result$u[, ncol(svd_result$u)]
    Q <- svd_result$v %*% t(svd_result$u)
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
#' @param source Source matrix.
#' @param target Target/reference matrix.
#' @param convention Multiplication convention (`"left"` or `"right"`).
#' @param scale Logical; if `TRUE`, include optimal scale factor.
#' @param reflection Logical; if `TRUE`, allow reflections (det(Q) may be -1).
#' @param obs_labels_source Optional observation labels for `source`.
#' @param obs_labels_target Optional observation labels for `target`.
#' @param min_overlap Minimum number of shared labels when labels are supplied.
#'
#' @return A list with elements `Q`, `scale_factor`, `residual`,
#'   `convention`, and `matched_labels`.
#'
#' @export
procrustes_rotation <- function(source,
                                target,
                                convention = c("left", "right"),
                                scale = FALSE,
                                reflection = FALSE,
                                obs_labels_source = NULL,
                                obs_labels_target = NULL,
                                min_overlap = 2L) {
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

    Q <- .procrustes_single(X, Y, scale = scale, reflection = reflection)
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
    min_overlap = min_overlap
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
#' @param x First matrix.
#' @param y Second matrix.
#' @param convention Multiplication convention (`"left"` or `"right"`).
#' @param obs_labels_x Optional observation labels for `x`.
#' @param obs_labels_y Optional observation labels for `y`.
#' @param min_overlap Minimum number of shared labels when labels are supplied.
#' @return Numeric scalar residual.
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
.gpa_builtin <- function(data_list, scale, reflection, tol, max_iter) {
  n_subjects <- length(data_list)
  subjects <- names(data_list)

  # Initialize with mean
  consensus <- Reduce(`+`, data_list) / n_subjects

  transforms <- vector("list", n_subjects)
  names(transforms) <- subjects

  for (iter in seq_len(max_iter)) {
    old_consensus <- consensus

    # Align each subject to current consensus
    aligned <- vector("list", n_subjects)
    for (i in seq_len(n_subjects)) {
      # Get left-multiply transform Q such that Q %*% X ≈ consensus
      Q <- .procrustes_single(data_list[[i]], consensus, scale, reflection)
      attr(Q, "scale_factor") <- NULL
      transforms[[i]] <- Q
      # Apply left-multiply: Q %*% X
      aligned[[i]] <- Q %*% data_list[[i]]
    }

    # Update consensus
    consensus <- Reduce(`+`, aligned) / n_subjects

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
