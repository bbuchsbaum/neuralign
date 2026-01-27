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
  .require_manifoldalign("CONE-Align")

  dots <- list(...)

  if (is.null(train_idx)) {
    train_idx <- seq_along(data@subjects)
  }
  train_data <- data[train_idx]
  data_list <- get_data_list(train_data)
  ref <- .resolve_reference(train_data, reference)

  hd <- .build_hyperdesign(train_data, transpose = FALSE)

  cone_args <- utils::modifyList(
    list(
      ncomp    = ncomp,
      sigma    = sigma,
      lambda   = lambda,
      max_iter = max_iter,
      tol      = tol
    ),
    dots
  )

  n_domains <- length(data_list)

  if (n_domains > 2) {
    cone_result <- do.call(
      manifoldalign::cone_align_multiple,
      c(list(data = hd), cone_args)
    )
  } else {
    cone_result <- do.call(
      manifoldalign::cone_align,
      c(list(data = hd), cone_args)
    )
  }

  transforms <- .extract_operators_from_scores(
    cone_result, data_list, ref$name, ref$data
  )

  # Held-out subjects
  heldout <- setdiff(data@subjects, names(transforms))
  if (length(heldout) > 0) {
    for (subj in heldout) {
      subj_data <- get_subject_data(data, subj)
      small_hd <- structure(
        list(subj = list(x = subj_data), ref = list(x = ref$data)),
        class = "hyperdesign"
      )
      small_result <- do.call(
        manifoldalign::cone_align,
        c(list(data = small_hd), cone_args)
      )
      small_dl <- list(subj = subj_data, ref = ref$data)
      transforms[[subj]] <- .extract_operators_from_scores(
        small_result, small_dl, "ref", ref$data
      )[["subj"]]
    }
  }

  list(
    transforms     = transforms,
    reference_data = ref$data,
    space_from     = train_data@space,
    space_to       = train_data@space,
    method_state   = list(reference = ref$name, cone_args = cone_args)
  )
}


.cone_capabilities <- list(
  supports_cv                  = TRUE,
  cv_axes                      = c("subject"),
  needs_geometry               = FALSE,
  needs_design                 = FALSE,
  requires_shared_features     = FALSE,
  requires_shared_observations = FALSE,
  returns_invertible           = FALSE,
  transform_type               = "linear",
  mass_preserving              = FALSE,
  returns                      = "operator",
  supports_new_subject         = TRUE,
  supports_new_data            = TRUE,
  reference_types              = c("subject", "medoid", "template")
)


#' @keywords internal
.register_cone <- function() {
  register_aligner(
    name         = "cone",
    fit_fn       = .cone_fit,
    apply_fn     = NULL,
    capabilities = .cone_capabilities,
    package      = "manifoldalign",
    description  = "CONE-Align graph alignment",
    version      = "0.1.0"
  )
}
