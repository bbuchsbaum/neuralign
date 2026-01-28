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
  dots <- list(...)

  gw_args <- utils::modifyList(
    list(epsilon = epsilon, max_iter = max_iter, tol = tol),
    dots
  )
  .ot_fit_common(
    data = data,
    reference = reference,
    train_idx = train_idx,
    solver_name = "gromov_wasserstein",
    solver_args = gw_args,
    method_label = "GW",
    args_name = "gw_args",
    method_state_extra = list(epsilon = gw_args$epsilon)
  )
}

.gw_apply <- function(fit_result, new_data, ...) {
  .ma_require_manifoldalign("GW")

  if (!inherits(new_data, "AlignmentData") || length(new_data@subjects) != 1L) {
    stop("gw apply_fn expects new_data to contain exactly one subject", call. = FALSE)
  }
  subj <- new_data@subjects[[1L]]
  X <- get_subject_data(new_data, subj)

  ref_data <- fit_result$reference_data %||% NULL
  if (is.null(ref_data)) stop("gw apply_fn requires reference_data in model", call. = FALSE)

  st <- fit_result$method_state %||% list()
  gw_args <- st$gw_args %||% list()
  gw_args <- utils::modifyList(gw_args, list(...))

  hd_new <- structure(list(
    subj = list(x = X),
    ref = list(x = ref_data)
  ), class = "hyperdesign")

  ot_new <- do.call(manifoldalign::gromov_wasserstein, c(list(data = hd_new), gw_args))
  P_subj_ref <- ot_new$transport_plans[[1L]]
  A_new <- .coupling_to_operator(P_subj_ref)

  list(transforms = setNames(list(A_new), subj))
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

.assert_finite_ot_matrix <- function(x, name = "matrix") {
  if (!.is_matrixish(x)) {
    stop(sprintf("OT alignment expects '%s' to be matrix-like", name), call. = FALSE)
  }

  if (inherits(x, "Matrix")) {
    vals <- x@x
    if (!is.null(vals) && any(!is.finite(vals))) {
      stop(sprintf("OT alignment expects '%s' to contain only finite values (NA/NaN/Inf)", name), call. = FALSE)
    }
    return(invisible(TRUE))
  }

  x <- as.matrix(x)
  if (any(!is.finite(x))) {
    stop(sprintf("OT alignment expects '%s' to contain only finite values (NA/NaN/Inf)", name), call. = FALSE)
  }
  invisible(TRUE)
}


.ot_fit_common <- function(data,
                           reference,
                           train_idx = NULL,
                           solver_name,
                           solver_args,
                           method_label,
                           args_name,
                           method_state_extra = list()) {
  .ma_require_manifoldalign(method_label)

  pre <- .aligner_preamble(data, train_idx = train_idx)
  train_data <- pre$train_data
  data_list <- pre$data_list

  for (subj in names(data_list)) {
    .assert_finite_ot_matrix(data_list[[subj]], name = sprintf("data for subject '%s'", subj))
  }

  if (.is_matrixish(reference)) {
    reference_name <- "template"
    reference_data <- as.matrix(reference)
  } else if (is.character(reference) && length(reference) == 1L) {
    if (reference %in% c("medoid", "centroid")) {
      ref_subj <- select_reference(train_data, method = reference)
      reference_name <- ref_subj
      reference_data <- get_subject_data(train_data, ref_subj)
    } else if (identical(reference, "consensus")) {
      reference_name <- "consensus"
      reference_data <- compute_centroid(train_data)
    } else if (identical(reference, "barycenter")) {
      stop(
        "GW barycenter reference is not implemented. Use reference='consensus' (arithmetic mean), a subject id, or a template matrix.",
        call. = FALSE
      )
    } else {
      reference_name <- reference
      reference_data <- get_subject_data(train_data, reference)
    }
  } else {
    stop("Invalid reference specification", call. = FALSE)
  }

  .assert_finite_ot_matrix(reference_data, name = "reference_data")

  # Build a minimal hyperdesign (list of domains with $x) for manifoldalign.
  # In neuralign, each subject matrix is (features x observations). For voxel-level
  # OT, treat rows (features) as samples and columns (observations/maps) as features.
  domain_names <- names(data_list)
  domains <- lapply(domain_names, function(nm) list(x = data_list[[nm]]))
  names(domains) <- domain_names

  # Ensure reference is present as a domain
  ref_idx <- match(reference_name, domain_names)
  if (is.na(ref_idx)) {
    domains[[reference_name]] <- list(x = reference_data)
    ref_idx <- length(domains)
  }

  hd <- structure(domains, class = "hyperdesign")

  solver <- get(solver_name, envir = asNamespace("manifoldalign"))
  ot <- do.call(solver, c(list(data = hd), solver_args))

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
    P <- .extract_pair_plan(ot$transport_plans, i, ref_idx, n_domains)
    # Ensure coupling is (source=subj x target=ref)
    P_subj_ref <- if (i < ref_idx) P else t(P)
    transforms[[subj]] <- .coupling_to_operator(P_subj_ref)
  }

  # Add transforms for held-out subjects (fit to reference only)
  heldout <- setdiff(data@subjects, names(transforms))
  if (length(heldout) > 0) {
    for (subj in heldout) {
      subj_data <- get_subject_data(data, subj)
      hd_new <- structure(list(
        subj = list(x = subj_data),
        ref = list(x = reference_data)
      ), class = "hyperdesign")
      ot_new <- do.call(solver, c(list(data = hd_new), solver_args))
      P_subj_ref <- ot_new$transport_plans[[1L]]
      transforms[[subj]] <- .coupling_to_operator(P_subj_ref)
    }
  }

  method_state <- method_state_extra
  method_state$reference <- reference_name
  method_state[[args_name]] <- solver_args

  list(
    transforms = transforms,
    reference_data = reference_data,
    space_from = train_data@space,
    space_to = train_data@space,
    method_state = method_state
  )
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
  # We return a row-normalized operator derived from the coupling (barycentric
  # projection). This preserves constant fields but does not preserve mass
  # (marginals) in the OT sense.
  mass_preserving = FALSE,
  returns = "operator",  # Couplings converted to operators via .coupling_to_operator
  supports_new_subject = TRUE,
  supports_new_data = TRUE,
  reference_types = c("subject", "consensus", "template")
)


#' Register GW Aligner
#' @keywords internal
.register_gw <- function() {
  register_aligner(
    name = "gw",
    fit_fn = .gw_fit,
    apply_fn = .gw_apply,
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
  dots <- list(...)

  fpgw_args <- utils::modifyList(
    list(omega1 = omega1, epsilon = epsilon, max_iter = max_iter, tol = tol),
    dots
  )
  .ot_fit_common(
    data = data,
    reference = reference,
    train_idx = train_idx,
    solver_name = "fpgw",
    solver_args = fpgw_args,
    method_label = "FPGW",
    args_name = "fpgw_args",
    method_state_extra = list(
      omega1 = fpgw_args$omega1,
      epsilon = fpgw_args$epsilon
    )
  )
}

.fpgw_apply <- function(fit_result, new_data, ...) {
  .ma_require_manifoldalign("FPGW")

  if (!inherits(new_data, "AlignmentData") || length(new_data@subjects) != 1L) {
    stop("fpgw apply_fn expects new_data to contain exactly one subject", call. = FALSE)
  }
  subj <- new_data@subjects[[1L]]
  X <- get_subject_data(new_data, subj)

  ref_data <- fit_result$reference_data %||% NULL
  if (is.null(ref_data)) stop("fpgw apply_fn requires reference_data in model", call. = FALSE)

  st <- fit_result$method_state %||% list()
  fpgw_args <- st$fpgw_args %||% list()
  fpgw_args <- utils::modifyList(fpgw_args, list(...))

  hd_new <- structure(list(
    subj = list(x = X),
    ref = list(x = ref_data)
  ), class = "hyperdesign")

  ot_new <- do.call(manifoldalign::fpgw, c(list(data = hd_new), fpgw_args))
  P_subj_ref <- ot_new$transport_plans[[1L]]
  A_new <- .coupling_to_operator(P_subj_ref)

  list(transforms = setNames(list(A_new), subj))
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
  # See .gw_capabilities: the returned operator is a barycentric projection
  # derived from the coupling, not a mass-preserving transport plan.
  mass_preserving = FALSE,
  returns = "operator",  # Couplings converted to operators via .coupling_to_operator
  supports_new_subject = TRUE,
  supports_new_data = TRUE,
  reference_types = c("subject", "consensus", "template")
)


#' Register FPGW Aligner
#' @keywords internal
.register_fpgw <- function() {
  register_aligner(
    name = "fpgw",
    fit_fn = .fpgw_fit,
    apply_fn = .fpgw_apply,
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
