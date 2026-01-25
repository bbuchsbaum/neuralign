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
    if (is.matrix(reference)) {
      reference_data <- reference
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
  if (is.matrix(reference)) {
    reference_data <- reference
  } else {
    reference_data <- get_subject_data(train_data, reference)
  }

  # Align each subject to reference
  # .procrustes_single returns left-multiply transform directly
  transforms <- lapply(names(data_list), function(subj) {
    if (is.character(reference) && subj == reference) {
      diag(nrow(data_list[[subj]]))
    } else {
      .procrustes_single(data_list[[subj]], reference_data, scale, reflection)
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

  Q
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
  cv_axes = c("subject"),
  needs_geometry = FALSE,
  needs_design = FALSE,
  returns_invertible = TRUE,
  transform_type = "orthogonal",
  mass_preserving = TRUE,
  returns = "operator",
  supports_new_subject = TRUE,
  supports_new_data = TRUE,
  reference_types = c("subject", "consensus")
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
