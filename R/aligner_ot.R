#' Optimal Transport Alignment
#'
#' Optimal transport-based alignment methods via manifoldalign.
#' Includes Gromov-Wasserstein (GW) and related methods.
#'
#' @name aligner_ot
NULL


#' Gromov-Wasserstein Fit Function
#'
#' Internal fit function for Gromov-Wasserstein alignment.
#'
#' @param data AlignmentData object.
#' @param reference Reference specification.
#' @param train_idx Indices of subjects to use for fitting.
#' @param epsilon Entropic regularization parameter.
#' @param max_iter Maximum iterations.
#' @param tol Convergence tolerance.
#' @param ... Additional arguments.
#'
#' @return List with transforms, reference_data, etc.
#'
#' @keywords internal
.gw_fit <- function(data,
                    reference = "medoid",
                    train_idx = NULL,
                    epsilon = 0.01,
                    max_iter = 100,
                    tol = 1e-6,
                    ...) {
  if (!requireNamespace("manifoldalign", quietly = TRUE)) {
    stop("Package 'manifoldalign' required for GW alignment")
  }

  # Handle train_idx
  if (is.null(train_idx)) {
    train_idx <- seq_along(data@subjects)
  }

  train_data <- data[train_idx]
  train_subjects <- train_data@subjects
  data_list <- get_data_list(train_data)

  # Resolve reference
  if (is.character(reference) && reference %in% c("medoid", "centroid")) {
    reference <- select_reference(train_data, method = reference)
  }

  if (is.character(reference) && reference != "barycenter") {
    reference_data <- get_subject_data(train_data, reference)
  } else if (is.matrix(reference)) {
    reference_data <- reference
  } else {
    # Barycenter reference - compute iteratively
    reference_data <- .compute_gw_barycenter(data_list, epsilon, max_iter, tol)
  }

  # Compute GW couplings to reference
  transforms <- lapply(names(data_list), function(subj) {
    X <- data_list[[subj]]

    if (is.character(reference) && subj == reference) {
      # Reference subject gets identity-like coupling
      n <- nrow(X)
      diag(n) / n * n  # Scaled identity
    } else {
      # Compute GW coupling
      coupling <- manifoldalign::gromov_wasserstein(
        X, reference_data,
        epsilon = epsilon,
        max_iter = max_iter,
        tol = tol
      )

      # Convert coupling to operator
      # Coupling P is (n_source x n_target)
      # For left-multiply convention, we need (n_target x n_source)
      .coupling_to_operator(coupling$P)
    }
  })
  names(transforms) <- names(data_list)

  # Add transforms for held-out subjects
  all_subjects <- data@subjects
  for (subj in all_subjects) {
    if (!subj %in% names(transforms)) {
      subj_data <- get_subject_data(data, subj)

      coupling <- manifoldalign::gromov_wasserstein(
        subj_data, reference_data,
        epsilon = epsilon,
        max_iter = max_iter,
        tol = tol
      )

      transforms[[subj]] <- .coupling_to_operator(coupling$P)
    }
  }

  list(
    transforms = transforms,
    reference_data = reference_data,
    space_from = train_data@space,
    space_to = train_data@space,
    method_state = list(
      epsilon = epsilon,
      reference = reference
    )
  )
}


#' Convert OT Coupling to Operator
#'
#' Convert an optimal transport coupling matrix to an operator
#' suitable for left-multiply application.
#'
#' @param P Coupling matrix (n_source x n_target).
#'
#' @return Operator matrix (n_target x n_source).
#'
#' @keywords internal
.coupling_to_operator <- function(P) {
  # Normalize rows to create a stochastic operator
  # P is (source x target), we want (target x source)
  P_t <- t(P)

  # Row-normalize for transport operator
  row_sums <- rowSums(P_t)
  row_sums[row_sums == 0] <- 1  # Avoid division by zero

  P_t / row_sums
}


#' Compute GW Barycenter
#' @keywords internal
.compute_gw_barycenter <- function(data_list, epsilon, max_iter, tol) {
  if (!requireNamespace("manifoldalign", quietly = TRUE)) {
    stop("Package 'manifoldalign' required for GW barycenter")
  }

  # Use manifoldalign's barycenter function if available
  tryCatch({
    manifoldalign::gw_barycenter(
      data_list,
      epsilon = epsilon,
      max_iter = max_iter,
      tol = tol
    )
  }, error = function(e) {
    # Fallback: use mean as approximate barycenter
    warning("GW barycenter not available; using arithmetic mean")
    Reduce(`+`, data_list) / length(data_list)
  })
}


#' GW Capabilities
#' @keywords internal
.gw_capabilities <- list(
  supports_cv = TRUE,
  cv_axes = c("subject"),
  needs_geometry = FALSE,
  needs_design = FALSE,
  returns_invertible = FALSE,  # OT couplings are not bijections
  transform_type = "ot",
  mass_preserving = TRUE,
  returns = "operator",  # Couplings converted to operators via .coupling_to_operator
  supports_new_subject = TRUE,
  supports_new_data = TRUE,
  reference_types = c("subject", "barycenter", "template")
)


#' Register GW Aligner
#' @keywords internal
.register_gw <- function() {
  register_aligner(
    name = "gw",
    fit_fn = .gw_fit,
    apply_fn = NULL,
    capabilities = .gw_capabilities,
    package = "manifoldalign",
    description = "Gromov-Wasserstein alignment",
    version = "0.1.0"
  )
}


#' Fused Procrustes-GW Fit Function
#'
#' Combines Procrustes rotation with GW optimal transport.
#'
#' @keywords internal
.fpgw_fit <- function(data,
                      reference = "medoid",
                      train_idx = NULL,
                      alpha = 0.5,
                      epsilon = 0.01,
                      max_iter = 100,
                      tol = 1e-6,
                      ...) {
  if (!requireNamespace("manifoldalign", quietly = TRUE)) {
    stop("Package 'manifoldalign' required for FPGW alignment")
  }

  # Handle train_idx
  if (is.null(train_idx)) {
    train_idx <- seq_along(data@subjects)
  }

  train_data <- data[train_idx]
  train_subjects <- train_data@subjects
  data_list <- get_data_list(train_data)

  # Resolve reference
  if (is.character(reference) && reference %in% c("medoid", "centroid")) {
    reference <- select_reference(train_data, method = reference)
  }

  if (is.character(reference)) {
    reference_data <- get_subject_data(train_data, reference)
  } else {
    reference_data <- reference
  }

  # Compute FPGW alignments
  transforms <- lapply(names(data_list), function(subj) {
    X <- data_list[[subj]]

    if (is.character(reference) && subj == reference) {
      diag(nrow(X))
    } else {
      result <- manifoldalign::fused_procrustes_gw(
        X, reference_data,
        alpha = alpha,
        epsilon = epsilon,
        max_iter = max_iter,
        tol = tol
      )

      # FPGW returns both rotation Q and coupling P
      # Combined operator: Q then P
      if (!is.null(result$Q)) {
        P_op <- .coupling_to_operator(result$P)
        P_op %*% t(result$Q)
      } else {
        .coupling_to_operator(result$P)
      }
    }
  })
  names(transforms) <- names(data_list)

  list(
    transforms = transforms,
    reference_data = reference_data,
    space_from = train_data@space,
    space_to = train_data@space,
    method_state = list(
      alpha = alpha,
      epsilon = epsilon
    )
  )
}


#' FPGW Capabilities
#' @keywords internal
.fpgw_capabilities <- list(
  supports_cv = TRUE,
  cv_axes = c("subject"),
  needs_geometry = FALSE,
  needs_design = FALSE,
  returns_invertible = FALSE,
  transform_type = "ot",
  mass_preserving = TRUE,
  returns = "operator",  # Couplings converted to operators via .coupling_to_operator
  supports_new_subject = TRUE,
  supports_new_data = TRUE,
  reference_types = c("subject", "template")
)


#' Register FPGW Aligner
#' @keywords internal
.register_fpgw <- function() {
  register_aligner(
    name = "fpgw",
    fit_fn = .fpgw_fit,
    apply_fn = NULL,
    capabilities = .fpgw_capabilities,
    package = "manifoldalign",
    description = "Fused Procrustes Gromov-Wasserstein",
    version = "0.1.0"
  )
}
