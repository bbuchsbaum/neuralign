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
                       ...) {
  .require_manifoldalign("GRASP")

  dots <- list(...)

  if (is.null(train_idx)) {
    train_idx <- seq_along(data@subjects)
  }
  train_data <- data[train_idx]
  data_list <- get_data_list(train_data)
  ref <- .resolve_reference(train_data, reference)

  hd <- .build_hyperdesign(train_data, transpose = FALSE)

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

  n_domains <- length(data_list)

  # Use multiset variant for 3+ domains, pairwise for 2
  if (n_domains > 2) {
    anchor_idx <- match(ref$name, names(data_list))
    if (is.na(anchor_idx)) anchor_idx <- 1L
    grasp_args$anchor <- anchor_idx
    grasp_result <- do.call(
      manifoldalign::grasp_multiset,
      c(list(data = hd), grasp_args)
    )
  } else {
    grasp_result <- do.call(
      manifoldalign::grasp,
      c(list(data = hd), grasp_args)
    )
  }

  # Extract operators from spectral embeddings
  transforms <- .extract_operators_from_scores(
    grasp_result, data_list, ref$name, ref$data
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
        manifoldalign::grasp,
        c(list(data = small_hd), grasp_args)
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
    method_state   = list(reference = ref$name, grasp_args = grasp_args)
  )
}


.grasp_capabilities <- list(
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
.register_grasp <- function() {
  register_aligner(
    name         = "grasp",
    fit_fn       = .grasp_fit,
    apply_fn     = NULL,
    capabilities = .grasp_capabilities,
    package      = "manifoldalign",
    description  = "GRASP spectral graph alignment",
    version      = "0.1.0"
  )
}
