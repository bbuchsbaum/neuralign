#' KEMA Alignment
#'
#' Kernel Manifold Alignment (KEMA) via manifoldalign. Produces per-subject
#' projection operators into a shared latent space.
#'
#' @name aligner_kema
#' @keywords internal
NULL


#' KEMA Fit Function
#'
#' @param data AlignmentData object.
#' @param reference Reference specification (subject id, "medoid", or "centroid").
#' @param train_idx Indices of subjects to use for fitting.
#' @param ncomp Number of latent components.
#' @param knn Number of nearest neighbours for graph construction.
#' @param sigma Kernel bandwidth (NULL for auto).
#' @param u Trade-off between manifold and class alignment (0-1).
#' @param solver KEMA solver: "regression" or "exact".
#' @param lambda Regularisation parameter.
#' @param ... Additional arguments forwarded to manifoldalign::kema()
#'   (except \code{preproc}; neuralign forces \code{preproc=multivarious::pass()}
#'   to keep transforms linear).
#'
#' @return List with transforms, reference_data, etc.
#'
#' @keywords internal
.kema_fit <- function(data,
                      reference = "medoid",
                      train_idx = NULL,
                      ncomp = 10L,
                      knn = 5L,
                      sigma = NULL,
                      u = 0.5,
                      solver = "regression",
                      lambda = 1e-2,
                      ...) {
  .ma_require_manifoldalign("KEMA")

  dots <- list(...)
  if ("preproc" %in% names(dots)) {
    stop(
      "manifoldalign 'preproc' is not supported via neuralign adapters. ",
      "neuralign forces preproc=multivarious::pass() to keep transforms linear. ",
      "Preprocess your inputs explicitly (e.g. preprocess_alignment_data(center='rows')) before fit_alignment().",
      call. = FALSE
    )
  }

  if (is.null(train_idx)) {
    train_idx <- seq_along(data@subjects)
  }
  train_data <- data[train_idx]
  labels <- .ma_require_shared_obs_labels(train_data, method = "kema")
  ref <- .ma_resolve_reference_spec(train_data, reference, method = "kema", allow_template = FALSE)

  # observations-as-samples: x = t(X), y = obs_labels
  label_name <- "label"
  y_sym <- as.name(label_name)
  hd <- .ma_build_hyperdesign_obs(train_data, labels = labels, label_name = label_name)

  # --- Assemble arguments ---
  kema_args <- utils::modifyList(
    list(
      y     = y_sym,
      preproc = multivarious::pass(),
      ncomp = ncomp,
      knn   = knn,
      sigma = sigma,
      u     = u,
      solver = solver,
      lambda = lambda
    ),
    dots
  )

  kema_result <- do.call(
    manifoldalign::kema,
    c(list(data = hd), kema_args)
  )

  data_list_train <- get_data_list(train_data)
  feature_counts <- vapply(data_list_train, nrow, integer(1))
  v_blocks <- .ma_mbp_split_loadings(kema_result, names(data_list_train), feature_counts)
  transforms_train <- .ma_projection_transforms_from_loadings(v_blocks)

  # Reference latent scores (k x n_obs)
  A_ref <- transforms_train[[ref$name]]
  Z_ref <- .ma_reference_scores(A_ref, ref$data)

  # Produce transforms for all subjects (train + heldout/new)
  data_list_all <- get_data_list(data)
  transforms <- list()
  for (subj in names(data_list_all)) {
    if (!is.null(transforms_train[[subj]])) {
      transforms[[subj]] <- transforms_train[[subj]]
    } else {
      X_subj <- data_list_all[[subj]]
      transforms[[subj]] <- .ma_ridge_map_to_reference_scores(X_subj, Z_ref, lambda = lambda)
    }
  }

  list(
    transforms = transforms,
    reference_data = Z_ref,
    space_from = train_data@space,
    space_to   = NULL,
    method_state = list(
      reference = ref$name,
      obs_labels_ref = labels,
      X_ref = ref$data,
      Z_ref = Z_ref,
      lambda = lambda,
      kema_args = kema_args
    )
  )
}

.kema_apply <- function(fit_result, new_data, ...) {
  if (!inherits(new_data, "AlignmentData") || length(new_data@subjects) != 1L) {
    stop("kema apply_fn expects new_data to contain exactly one subject", call. = FALSE)
  }
  subj <- new_data@subjects[[1L]]
  X <- get_subject_data(new_data, subj)

  st <- fit_result$method_state %||% list()
  Z_ref <- st$Z_ref %||% NULL
  obs_ref <- st$obs_labels_ref %||% NULL
  lambda <- st$lambda %||% 1e-2
  if (is.null(Z_ref) || is.null(obs_ref)) {
    stop("kema apply_fn missing reference latent scores/labels in method_state", call. = FALSE)
  }

  obs_new <- .ma_get_single_subject_obs_labels(new_data)
  if (is.null(obs_new)) {
    stop("kema apply_fn requires obs_labels on new_data", call. = FALSE)
  }
  idx <- .ma_match_obs_indices(obs_ref, obs_new)
  Xc <- X[, idx$new, drop = FALSE]
  Zc <- Z_ref[, idx$ref, drop = FALSE]
  A_new <- .ma_ridge_map_to_reference_scores(Xc, Zc, lambda = lambda)
  list(transforms = setNames(list(A_new), subj))
}


#' KEMA Capabilities
#' @keywords internal
.kema_capabilities <- list(
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


#' Register KEMA Aligner
#' @keywords internal
.register_kema <- function() {
  register_aligner(
    name        = "kema",
    fit_fn      = .kema_fit,
    apply_fn    = .kema_apply,
    capabilities = .kema_capabilities,
    package     = "manifoldalign",
    description = "Kernel Manifold Alignment (KEMA)",
    version     = "0.1.0"
  )
}
