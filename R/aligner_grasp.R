#' GRASP Alignment
#'
#' Graph Alignment by Spectral Corresponding Functions via manifoldalign.
#' Uses spectral embeddings and multi-scale descriptors to find node
#' correspondences between domains.
#'
#' @name aligner_grasp
#' @keywords internal
NULL


#' GRASP Fit Function
#'
#' @param data AlignmentData object.
#' @param reference Reference specification.
#' @param train_idx Indices of subjects to use for fitting.
#' @param ncomp Number of spectral components.
#' @param q_descriptors Number of diffusion-time descriptors.
#' @param sigma Diffusion bandwidth.
#' @param lambda Regularisation for base alignment.
#' @param solver Assignment solver: "linear" or "hungarian".
#' @param ... Additional arguments forwarded to manifoldalign::grasp()
#'   or manifoldalign::grasp_multiset().
#'
#' @return List with transforms, reference_data, etc.
#' @keywords internal
.grasp_fit <- function(data,
                       reference = "medoid",
                       train_idx = NULL,
                       ncomp = 30L,
                       q_descriptors = 100L,
                       sigma = 0.73,
                       lambda = 0.1,
                       solver = "linear",
                       fit_context = NULL,
                       provider_plan = NULL,
                       ...) {
  .ma_require_manifoldalign("GRASP")

  dots <- list(...)

  if (is.null(train_idx)) {
    train_idx <- seq_along(data@subjects)
  }
  train_data <- data[train_idx]
  ref <- .ma_resolve_reference_spec(train_data, reference, method = "grasp", allow_template = TRUE)

  grasp_args <- utils::modifyList(
    list(
      ncomp         = ncomp,
      q_descriptors = q_descriptors,
      sigma         = sigma,
      lambda        = lambda,
      solver        = solver
    ),
    dots
  )
  data_list_all <- get_data_list(data)
  n_target <- nrow(ref$data)

  # Pairwise align each domain to the reference, returning a sparse
  # permutation-style operator derived from the assignment.
  transforms <- list()
  for (subj in names(data_list_all)) {
    Xi <- data_list_all[[subj]]
    if (identical(ref$name, subj)) {
      transforms[[subj]] <- Matrix::Diagonal(n_target)
      next
    }

    pair_hd <- .ma_build_pair_hyperdesign_features(
      ref$data, Xi, target_name = "target", source_name = "source"
    )
    fit <- do.call(manifoldalign::grasp, c(list(data = pair_hd), grasp_args))
    P <- fit$assignment
    transforms[[subj]] <- .ma_assignment_to_operator(P, n_target = n_target, n_source = nrow(Xi))
  }

  list(
    transforms     = transforms,
    reference_data = ref$data,
    space_from     = train_data@space,
    space_to       = train_data@space,
    method_state   = list(reference = ref$name, grasp_args = grasp_args)
  )
}

.grasp_apply <- function(fit_result, new_data, ...) {
  .ma_require_manifoldalign("GRASP")
  if (!inherits(new_data, "AlignmentData") || length(new_data@subjects) != 1L) {
    stop("grasp apply_fn expects new_data to contain exactly one subject", call. = FALSE)
  }
  subj <- new_data@subjects[[1L]]
  X <- get_subject_data(new_data, subj)

  ref_data <- fit_result$reference_data %||% NULL
  st <- fit_result$method_state %||% list()
  grasp_args <- st$grasp_args %||% list()
  if (is.null(ref_data)) stop("grasp apply_fn requires reference_data in model", call. = FALSE)

  n_target <- nrow(ref_data)
  pair_hd <- .ma_build_pair_hyperdesign_features(
    ref_data, X, target_name = "target", source_name = "source"
  )
  fit <- do.call(manifoldalign::grasp, c(list(data = pair_hd), grasp_args))
  P <- fit$assignment
  A_new <- .ma_assignment_to_operator(P, n_target = n_target, n_source = nrow(X))
  list(transforms = setNames(list(A_new), subj))
}


.grasp_capabilities <- list(
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
.register_grasp <- function() {
  register_aligner(
    name         = "grasp",
    fit_fn       = .grasp_fit,
    apply_fn     = .grasp_apply,
    capabilities = .grasp_capabilities,
    package      = "manifoldalign",
    description  = "GRASP spectral graph alignment",
    version      = "0.1.0"
  )
}
