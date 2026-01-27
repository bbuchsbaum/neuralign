#' Coupled Diagonalization Alignment
#'
#' Coupled matrix diagonalization via manifoldalign. Finds joint spectral
#' bases across domains that diagonalise their respective graph Laplacians
#' while staying coupled through a shared latent subspace.
#'
#' @name aligner_coupled_diag
#' @keywords internal
NULL


#' Coupled Diagonalization Fit Function
#'
#' @param data AlignmentData object.
#' @param reference Reference specification.
#' @param train_idx Indices of subjects to use for fitting.
#' @param ncomp Number of coupled components.
#' @param ncomp_per_domain Number of per-domain eigenvectors.
#' @param mu_coupling Coupling strength.
#' @param knn Neighbours for graph construction.
#' @param sigma Kernel bandwidth (NULL for auto).
#' @param max_iter Maximum optimisation iterations.
#' @param tol Convergence tolerance.
#' @param ... Additional arguments forwarded to
#'   manifoldalign::coupled_diagonalization().
#'
#' @return List with transforms, reference_data, etc.
#'
#' @keywords internal
.coupled_diag_fit <- function(data,
                              reference = "medoid",
                              train_idx = NULL,
                              ncomp = 10L,
                              ncomp_per_domain = 20L,
                              mu_coupling = 1,
                              knn = 10L,
                              sigma = NULL,
                              max_iter = 200L,
                              tol = 1e-6,
                              ...) {
  .require_manifoldalign("coupled diagonalization")

  dots <- list(...)

  if (is.null(train_idx)) {
    train_idx <- seq_along(data@subjects)
  }
  train_data <- data[train_idx]
  data_list <- get_data_list(train_data)

  ref <- .resolve_reference(train_data, reference)

  # Build hyperdesign (no transpose: features-as-points)
  hd <- .build_hyperdesign(train_data, transpose = FALSE)

  cd_args <- utils::modifyList(
    list(
      ncomp            = ncomp,
      ncomp_per_domain = ncomp_per_domain,
      mu_coupling      = mu_coupling,
      knn              = knn,
      sigma            = sigma,
      max_iter         = max_iter,
      tol              = tol
    ),
    dots
  )

  cd_result <- do.call(
    manifoldalign::coupled_diagonalization,
    c(list(data = hd), cd_args)
  )

  # Extract per-subject operators from aligned scores
  transforms <- .extract_operators_from_scores(
    cd_result, data_list, ref$name, ref$data
  )

  # Handle held-out subjects
  all_subjects <- data@subjects
  heldout <- setdiff(all_subjects, names(transforms))
  if (length(heldout) > 0) {
    for (subj in heldout) {
      subj_data <- get_subject_data(data, subj)
      small_hd <- structure(
        list(
          subj = list(x = subj_data),
          ref  = list(x = ref$data)
        ),
        class = "hyperdesign"
      )
      small_result <- do.call(
        manifoldalign::coupled_diagonalization,
        c(list(data = small_hd), cd_args)
      )
      small_dl <- list(subj = subj_data, ref = ref$data)
      small_transforms <- .extract_operators_from_scores(
        small_result, small_dl, "ref", ref$data
      )
      transforms[[subj]] <- small_transforms[["subj"]]
    }
  }

  list(
    transforms = transforms,
    reference_data = ref$data,
    space_from = train_data@space,
    space_to   = train_data@space,
    method_state = list(
      ncomp       = ncomp,
      mu_coupling = mu_coupling,
      reference   = ref$name,
      cd_args     = cd_args
    )
  )
}


#' Coupled Diagonalization Capabilities
#' @keywords internal
.coupled_diag_capabilities <- list(
  supports_cv                = TRUE,
  cv_axes                    = c("subject"),
  needs_geometry             = FALSE,
  needs_design               = FALSE,
  requires_shared_features   = TRUE,
  requires_shared_observations = FALSE,
  returns_invertible         = FALSE,
  transform_type             = "linear",
  mass_preserving            = FALSE,
  returns                    = "operator",
  supports_new_subject       = TRUE,
  supports_new_data          = TRUE,
  reference_types            = c("subject", "medoid", "template")
)


#' Register Coupled Diagonalization Aligner
#' @keywords internal
.register_coupled_diag <- function() {
  register_aligner(
    name         = "coupled_diag",
    fit_fn       = .coupled_diag_fit,
    apply_fn     = NULL,
    capabilities = .coupled_diag_capabilities,
    package      = "manifoldalign",
    description  = "Coupled Diagonalization alignment",
    version      = "0.1.0"
  )
}
