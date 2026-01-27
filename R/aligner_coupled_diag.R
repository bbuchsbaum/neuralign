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
                              lambda_oos = 1e-2,
                              ...) {
  .ma_require_manifoldalign("coupled diagonalization")

  dots <- list(...)

  if (is.null(train_idx)) {
    train_idx <- seq_along(data@subjects)
  }
  train_data <- data[train_idx]
  labels <- .ma_require_shared_obs_labels(train_data, method = "coupled_diag")
  ref <- .ma_resolve_reference_spec(train_data, reference, method = "coupled_diag", allow_template = FALSE)

  # observations-as-samples: x = t(X)
  hd <- .ma_build_hyperdesign_obs(train_data, labels = labels, label_name = "label")

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

  data_list_train <- get_data_list(train_data)
  feature_counts <- vapply(data_list_train, nrow, integer(1))
  v_blocks <- .ma_mbp_split_loadings(cd_result, names(data_list_train), feature_counts)
  transforms_train <- .ma_projection_transforms_from_loadings(v_blocks)

  A_ref <- transforms_train[[ref$name]]
  Z_ref <- .ma_reference_scores(A_ref, ref$data)

  data_list_all <- get_data_list(data)
  transforms <- list()
  for (subj in names(data_list_all)) {
    if (!is.null(transforms_train[[subj]])) {
      transforms[[subj]] <- transforms_train[[subj]]
    } else {
      transforms[[subj]] <- .ma_ridge_map_to_reference_scores(data_list_all[[subj]], Z_ref, lambda = lambda_oos)
    }
  }

  list(
    transforms = transforms,
    reference_data = Z_ref,
    space_from = train_data@space,
    space_to   = NULL,
    method_state = list(
      ncomp       = ncomp,
      mu_coupling = mu_coupling,
      reference   = ref$name,
      obs_labels_ref = labels,
      X_ref = ref$data,
      Z_ref = Z_ref,
      lambda_oos = lambda_oos,
      cd_args     = cd_args
    )
  )
}

.coupled_diag_apply <- function(fit_result, new_data, ...) {
  if (!inherits(new_data, "AlignmentData") || length(new_data@subjects) != 1L) {
    stop("coupled_diag apply_fn expects new_data to contain exactly one subject", call. = FALSE)
  }
  subj <- new_data@subjects[[1L]]
  X <- get_subject_data(new_data, subj)

  st <- fit_result$method_state %||% list()
  Z_ref <- st$Z_ref %||% NULL
  obs_ref <- st$obs_labels_ref %||% NULL
  lambda <- st$lambda_oos %||% 1e-2
  if (is.null(Z_ref) || is.null(obs_ref)) {
    stop("coupled_diag apply_fn missing reference latent scores/labels in method_state", call. = FALSE)
  }
  obs_new <- .ma_get_single_subject_obs_labels(new_data)
  if (is.null(obs_new)) stop("coupled_diag apply_fn requires obs_labels on new_data", call. = FALSE)
  idx <- .ma_match_obs_indices(obs_ref, obs_new)
  A_new <- .ma_ridge_map_to_reference_scores(
    X[, idx$new, drop = FALSE],
    Z_ref[, idx$ref, drop = FALSE],
    lambda = lambda
  )
  list(transforms = setNames(list(A_new), subj))
}


#' Coupled Diagonalization Capabilities
#' @keywords internal
.coupled_diag_capabilities <- list(
  supports_cv                = TRUE,
  cv_axes                    = c("subject"),
  needs_geometry             = FALSE,
  needs_design               = FALSE,
  requires_shared_features   = TRUE,
  requires_shared_observations = TRUE,
  returns_invertible         = FALSE,
  transform_type             = "linear",
  mass_preserving            = FALSE,
  returns                    = "operator",
  supports_new_subject       = TRUE,
  supports_new_data          = TRUE,
  reference_types            = c("subject")
)


#' Register Coupled Diagonalization Aligner
#' @keywords internal
.register_coupled_diag <- function() {
  register_aligner(
    name         = "coupled_diag",
    fit_fn       = .coupled_diag_fit,
    apply_fn     = .coupled_diag_apply,
    capabilities = .coupled_diag_capabilities,
    package      = "manifoldalign",
    description  = "Coupled Diagonalization alignment",
    version      = "0.1.0"
  )
}
