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

  dots <- list(...)

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

  # Determine reference matrix + name (for provenance)
  reference_data <- NULL
  reference_name <- reference
  if (is.character(reference) && reference != "barycenter") {
    reference_data <- get_subject_data(train_data, reference)
  } else if (.is_matrixish(reference)) {
    reference_data <- as.matrix(reference)
    reference_name <- "template"
  } else {
    reference_data <- .compute_gw_barycenter(data_list, epsilon, max_iter, tol)
    reference_name <- "barycenter"
  }

  # Build a minimal hyperdesign (list of domains with $x) for manifoldalign.
  # In neuralign, each subject matrix is (features x observations).
  # For fMRI alignment use-cases, treat rows (features/voxels) as samples and
  # columns (maps/contrasts) as features, which matches manifoldalign's expectations.
  domain_names <- names(data_list)
  domains <- lapply(data_list, function(X) list(x = X))
  names(domains) <- domain_names

  # Ensure reference is present as a domain
  ref_idx <- match(reference, domain_names)
  if (is.na(ref_idx)) {
    domains[[reference_name]] <- list(x = reference_data)
    ref_idx <- length(domains)
  }

  hd <- structure(domains, class = "hyperdesign")

  gw_args <- utils::modifyList(
    list(epsilon = epsilon, max_iter = max_iter, tol = tol),
    dots
  )
  gw <- do.call(manifoldalign::gromov_wasserstein, c(list(data = hd), gw_args))

  # Extract a transport plan for each subject -> reference and convert to operator.
  n_domains <- length(domains)
  transforms <- list()
  for (i in seq_along(domain_names)) {
    subj <- domain_names[[i]]
    Xi <- data_list[[subj]]
    if (i == ref_idx) {
      transforms[[subj]] <- diag(nrow(Xi))
      next
    }
    P <- .extract_pair_plan(gw$transport_plans, i, ref_idx, n_domains)
    # Ensure coupling is (source=subj x target=ref)
    P_subj_ref <- if (i < ref_idx) P else t(P)
    transforms[[subj]] <- .coupling_to_operator(P_subj_ref)
  }

  # Add transforms for held-out subjects (fit to reference only)
  all_subjects <- data@subjects
  heldout <- setdiff(all_subjects, names(transforms))
  if (length(heldout) > 0) {
    for (subj in heldout) {
      subj_data <- get_subject_data(data, subj)
      hd_new <- structure(list(
        subj = list(x = subj_data),
        ref = list(x = reference_data)
      ), class = "hyperdesign")
      gw_new <- do.call(manifoldalign::gromov_wasserstein, c(list(data = hd_new), gw_args))
      P_subj_ref <- gw_new$transport_plans[[1L]]
      transforms[[subj]] <- .coupling_to_operator(P_subj_ref)
    }
  }

  list(
    transforms = transforms,
    reference_data = reference_data,
    space_from = train_data@space,
    space_to = train_data@space,
    method_state = list(
      epsilon = epsilon,
      reference = reference_name,
      gw_args = gw_args
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
  # Fallback: use mean as approximate barycenter (requires matched feature dims).
  warning("GW barycenter not implemented; using arithmetic mean", call. = FALSE)
  Reduce(`+`, data_list) / length(data_list)
}


#' GW Capabilities
#' @keywords internal
.gw_capabilities <- list(
  supports_cv = TRUE,
  cv_axes = c("subject"),
  needs_geometry = FALSE,
  needs_design = FALSE,
  requires_shared_features = FALSE,
  requires_shared_observations = FALSE,
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
                      omega1 = 0.001,
                      epsilon = 0.01,
                      max_iter = 100,
                      tol = 1e-6,
                      ...) {
  if (!requireNamespace("manifoldalign", quietly = TRUE)) {
    stop("Package 'manifoldalign' required for FPGW alignment")
  }

  dots <- list(...)

  # Handle train_idx
  if (is.null(train_idx)) {
    train_idx <- seq_along(data@subjects)
  }

  train_data <- data[train_idx]
  data_list <- get_data_list(train_data)

  # Resolve reference
  if (is.character(reference) && reference %in% c("medoid", "centroid")) {
    reference <- select_reference(train_data, method = reference)
  }

  reference_data <- NULL
  reference_name <- reference
  if (is.character(reference) && reference != "barycenter") {
    reference_data <- get_subject_data(train_data, reference)
  } else if (.is_matrixish(reference)) {
    reference_data <- as.matrix(reference)
    reference_name <- "template"
  } else {
    reference_data <- Reduce(`+`, data_list) / length(data_list)
    reference_name <- "barycenter"
  }

  domain_names <- names(data_list)
  domains <- lapply(data_list, function(X) list(x = X))
  names(domains) <- domain_names
  ref_idx <- match(reference, domain_names)
  if (is.na(ref_idx)) {
    domains[[reference_name]] <- list(x = reference_data)
    ref_idx <- length(domains)
  }
  hd <- structure(domains, class = "hyperdesign")

  fpgw_args <- utils::modifyList(
    list(omega1 = omega1, epsilon = epsilon, max_iter = max_iter, tol = tol),
    dots
  )
  fpgw <- do.call(manifoldalign::fpgw, c(list(data = hd), fpgw_args))

  n_domains <- length(domains)
  transforms <- list()
  for (i in seq_along(domain_names)) {
    subj <- domain_names[[i]]
    Xi <- data_list[[subj]]
    if (i == ref_idx) {
      transforms[[subj]] <- diag(nrow(Xi))
      next
    }
    P <- .extract_pair_plan(fpgw$transport_plans, i, ref_idx, n_domains)
    P_subj_ref <- if (i < ref_idx) P else t(P)
    transforms[[subj]] <- .coupling_to_operator(P_subj_ref)
  }

  list(
    transforms = transforms,
    reference_data = reference_data,
    space_from = train_data@space,
    space_to = train_data@space,
    method_state = list(
      omega1 = omega1,
      epsilon = epsilon,
      reference = reference_name,
      fpgw_args = fpgw_args
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
  requires_shared_features = FALSE,
  requires_shared_observations = FALSE,
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
    description = "Fused-Partial Gromov-Wasserstein alignment",
    version = "0.1.0"
  )
}

#' Extract a pairwise transport plan from a packed (i,j) list
#'
#' manifoldalign stores pairwise plans in the order:
#' (1,2), (1,3), ..., (1,n), (2,3), (2,4), ..., (n-1,n)
#'
#' @param plans List of plans in packed order.
#' @param i First domain index.
#' @param j Second domain index.
#' @param n Total number of domains.
#'
#' @return Transport plan matrix for the (i,j) pair.
#' @keywords internal
.extract_pair_plan <- function(plans, i, j, n) {
  if (!is.numeric(i) || !is.numeric(j) || i == j) {
    stop("Invalid pair indices", call. = FALSE)
  }
  if (i > j) {
    tmp <- i
    i <- j
    j <- tmp
  }
  # Number of pairs before i: sum_{a=1}^{i-1} (n-a)
  before_i <- (i - 1) * (2 * n - i) / 2
  idx <- as.integer(before_i + (j - i))
  plans[[idx]]
}
