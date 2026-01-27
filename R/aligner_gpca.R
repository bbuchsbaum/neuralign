#' GPCA Alignment
#'
#' Generalized PCA alignment via manifoldalign. Combines within- and
#' between-domain similarity constraints in a structured PCA framework.
#'
#' @name aligner_gpca
#' @keywords internal
NULL


#' GPCA Fit Function
#'
#' @param data AlignmentData object.
#' @param reference Reference specification.
#' @param train_idx Indices of subjects to use for fitting.
#' @param ncomp Number of components.
#' @param u Trade-off parameter (0-1).
#' @param lambda Regularisation.
#' @param ... Additional arguments forwarded to manifoldalign::gpca_align().
#'
#' @return List with transforms, reference_data, etc.
#' @keywords internal
.gpca_fit <- function(data,
                      reference = "medoid",
                      train_idx = NULL,
                      ncomp = 10L,
                      u = 0.5,
                      lambda = 1e-2,
                      ...) {
  .ma_require_manifoldalign("GPCA")

  dots <- list(...)

  if (is.null(train_idx)) {
    train_idx <- seq_along(data@subjects)
  }
  train_data <- data[train_idx]
  labels <- .ma_require_shared_obs_labels(train_data, method = "gpca")
  ref <- .ma_resolve_reference_spec(train_data, reference, method = "gpca", allow_template = FALSE)

  label_name <- "label"
  y_sym <- as.name(label_name)
  hd <- .ma_build_hyperdesign_obs(train_data, labels = labels, label_name = label_name)

  gpca_args <- utils::modifyList(
    list(y = y_sym, ncomp = ncomp, u = u, lambda = lambda),
    dots
  )

  gpca_result <- do.call(
    manifoldalign::gpca_align.hyperdesign,
    c(list(data = hd), gpca_args)
  )

  data_list_train <- get_data_list(train_data)
  feature_counts <- vapply(data_list_train, nrow, integer(1))
  v_blocks <- .ma_mbp_split_loadings(gpca_result, names(data_list_train), feature_counts)
  transforms_train <- .ma_projection_transforms_from_loadings(v_blocks)

  A_ref <- transforms_train[[ref$name]]
  Z_ref <- .ma_reference_scores(A_ref, ref$data)

  data_list_all <- get_data_list(data)
  transforms <- list()
  for (subj in names(data_list_all)) {
    if (!is.null(transforms_train[[subj]])) {
      transforms[[subj]] <- transforms_train[[subj]]
    } else {
      transforms[[subj]] <- .ma_ridge_map_to_reference_scores(data_list_all[[subj]], Z_ref, lambda = lambda)
    }
  }

  list(
    transforms     = transforms,
    reference_data = ref$data,
    space_from     = train_data@space,
    space_to       = NULL,
    method_state   = list(
      reference = ref$name,
      obs_labels_ref = labels,
      Z_ref = Z_ref,
      lambda = lambda,
      gpca_args = gpca_args
    )
  )
}

.gpca_apply <- function(fit_result, new_data, ...) {
  if (!inherits(new_data, "AlignmentData") || length(new_data@subjects) != 1L) {
    stop("gpca apply_fn expects new_data to contain exactly one subject", call. = FALSE)
  }
  subj <- new_data@subjects[[1L]]
  X <- get_subject_data(new_data, subj)

  st <- fit_result$method_state %||% list()
  Z_ref <- st$Z_ref %||% NULL
  obs_ref <- st$obs_labels_ref %||% NULL
  lambda <- st$lambda %||% 1e-2
  if (is.null(Z_ref) || is.null(obs_ref)) {
    stop("gpca apply_fn missing reference latent scores/labels in method_state", call. = FALSE)
  }

  obs_new <- .ma_get_single_subject_obs_labels(new_data)
  if (is.null(obs_new)) stop("gpca apply_fn requires obs_labels on new_data", call. = FALSE)
  idx <- .ma_match_obs_indices(obs_ref, obs_new)
  A_new <- .ma_ridge_map_to_reference_scores(
    X[, idx$new, drop = FALSE],
    Z_ref[, idx$ref, drop = FALSE],
    lambda = lambda
  )
  list(transforms = setNames(list(A_new), subj))
}


.gpca_capabilities <- list(
  supports_cv                  = TRUE,
  cv_axes                      = c("subject"),
  needs_geometry               = FALSE,
  needs_design                 = FALSE,
  requires_shared_features     = TRUE,
  requires_shared_observations = TRUE,
  returns_invertible           = FALSE,
  transform_type               = "linear",
  mass_preserving              = FALSE,
  returns                      = "operator",
  supports_new_subject         = TRUE,
  supports_new_data            = TRUE,
  reference_types              = c("subject")
)


#' @keywords internal
.register_gpca <- function() {
  register_aligner(
    name         = "gpca",
    fit_fn       = .gpca_fit,
    apply_fn     = .gpca_apply,
    capabilities = .gpca_capabilities,
    package      = "manifoldalign",
    description  = "Generalized PCA alignment",
    version      = "0.1.0"
  )
}
