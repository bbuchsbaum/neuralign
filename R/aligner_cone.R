#' CONE-Align
#'
#' Consensus Optimisation for Node Embedding Alignment via manifoldalign.
#' Iteratively alternates between Procrustes rotation of spectral
#' embeddings and Hungarian assignment of node correspondences.
#'
#' @name aligner_cone
#' @keywords internal
NULL


#' CONE-Align Fit Function
#'
#' @param data AlignmentData object.
#' @param reference Reference specification.
#' @param train_idx Indices of subjects to use for fitting.
#' @param ncomp Number of spectral components.
#' @param sigma Bandwidth for kNN graph.
#' @param lambda Regularisation.
#' @param max_iter Maximum iterations.
#' @param tol Convergence tolerance.
#' @param ... Additional arguments forwarded to
#'   manifoldalign::cone_align() or manifoldalign::cone_align_multiple().
#'
#' @return List with transforms, reference_data, etc.
#' @keywords internal
.cone_fit <- function(data,
                      reference = "medoid",
                      train_idx = NULL,
                      ncomp = 10L,
                      sigma = 0.73,
                      lambda = 0.1,
                      max_iter = 30L,
                      tol = 0.01,
                      ...) {
  .ma_require_manifoldalign("CONE-Align")

  dots <- list(...)

  if (is.null(train_idx)) {
    train_idx <- seq_along(data@subjects)
  }
  train_data <- data[train_idx]
  ref <- .ma_resolve_reference_spec(train_data, reference, method = "cone", allow_template = TRUE)

  cone_args <- utils::modifyList(
    list(ncomp = ncomp, sigma = sigma, lambda = lambda, max_iter = max_iter, tol = tol),
    dots
  )

  # Pairwise align each domain to the reference, and build a sparse permutation
  # operator from the returned assignment. This is the primary output of CONE-Align.
  data_list_all <- get_data_list(data)
  n_target <- nrow(ref$data)

  transforms <- list()
  for (subj in names(data_list_all)) {
    Xi <- data_list_all[[subj]]
    if (identical(ref$name, subj)) {
      transforms[[subj]] <- Matrix::Diagonal(n_target)
      next
    }

    pair_hd <- structure(
      list(
        target = list(x = ref$data),
        source = list(x = Xi)
      ),
      class = c("hyperdesign", "list")
    )
    fit <- do.call(manifoldalign::cone_align, c(list(data = pair_hd), cone_args))
    P <- fit$assignment
    transforms[[subj]] <- .ma_assignment_to_operator(P, n_target = n_target, n_source = nrow(Xi))
  }

  list(
    transforms     = transforms,
    reference_data = ref$data,
    space_from     = train_data@space,
    space_to       = train_data@space,
    method_state   = list(reference = ref$name, cone_args = cone_args)
  )
}

.cone_apply <- function(fit_result, new_data, ...) {
  .ma_require_manifoldalign("CONE-Align")
  if (!inherits(new_data, "AlignmentData") || length(new_data@subjects) != 1L) {
    stop("cone apply_fn expects new_data to contain exactly one subject", call. = FALSE)
  }
  subj <- new_data@subjects[[1L]]
  X <- get_subject_data(new_data, subj)

  ref_data <- fit_result$reference_data %||% NULL
  st <- fit_result$method_state %||% list()
  cone_args <- st$cone_args %||% list()
  if (is.null(ref_data)) stop("cone apply_fn requires reference_data in model", call. = FALSE)

  n_target <- nrow(ref_data)
  pair_hd <- structure(
    list(target = list(x = ref_data), source = list(x = X)),
    class = c("hyperdesign", "list")
  )
  fit <- do.call(manifoldalign::cone_align, c(list(data = pair_hd), cone_args))
  P <- fit$assignment
  A_new <- .ma_assignment_to_operator(P, n_target = n_target, n_source = nrow(X))
  list(transforms = setNames(list(A_new), subj))
}


.cone_capabilities <- list(
  supports_cv                  = TRUE,
  cv_axes                      = c("subject"),
  needs_geometry               = FALSE,
  needs_design                 = FALSE,
  requires_shared_features     = FALSE,
  requires_shared_observations = FALSE,
  returns_invertible           = FALSE,
  transform_type               = "permutation",
  mass_preserving              = FALSE,
  returns                      = "operator",
  supports_new_subject         = TRUE,
  supports_new_data            = TRUE,
  reference_types              = c("subject", "template")
)


#' @keywords internal
.register_cone <- function() {
  register_aligner(
    name         = "cone",
    fit_fn       = .cone_fit,
    apply_fn     = .cone_apply,
    capabilities = .cone_capabilities,
    package      = "manifoldalign",
    description  = "CONE-Align graph alignment",
    version      = "0.1.0"
  )
}
