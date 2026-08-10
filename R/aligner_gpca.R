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
#' @param ... Additional arguments forwarded to manifoldalign::gpca_align()
#'   (except \code{preproc}; neuralign forces \code{preproc=multivarious::pass()}
#'   to keep transforms linear).
#'
#' @return List with transforms, reference_data, etc.
#' @keywords internal
.gpca_fit <- function(data,
                      reference = "medoid",
                      train_idx = NULL,
                      ncomp = 10L,
                      u = 0.5,
                      lambda = 1e-2,
                      target_space = c("latent", "reference"),
                      fit_context = NULL,
                      provider_plan = NULL,
                      ...) {
  .ma_require_manifoldalign("GPCA")

  target_space <- match.arg(target_space)

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
  labels <- .ma_require_shared_obs_labels(train_data, method = "gpca")
  ref <- .ma_resolve_reference_spec(train_data, reference, method = "gpca", allow_template = FALSE)

  label_name <- "label"
  y_sym <- as.name(label_name)
  hd <- .ma_build_hyperdesign_obs(train_data, labels = labels, label_name = label_name)

  gpca_args <- utils::modifyList(
    list(y = y_sym, preproc = multivarious::pass(), ncomp = ncomp, u = u, lambda = lambda),
    dots
  )

  gpca_result <- do.call(
    manifoldalign::gpca_align,
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

  U_ref <- NULL
  if (identical(target_space, "reference")) {
    lifted <- .ma_lift_latent_transforms_to_reference(transforms, ref$name, context = "gpca")
    transforms <- lifted$transforms
    U_ref <- lifted$U_ref
  }

  list(
    transforms     = transforms,
    reference_data = Z_ref,
    space_from     = train_data@space,
    space_to       = NULL,
    method_state   = list(
      reference = ref$name,
      target_space = target_space,
      U_ref = U_ref,
      obs_labels_ref = labels,
      X_ref = ref$data,
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

  target_space <- st$target_space %||% "latent"
  if (identical(target_space, "reference")) {
    U_ref <- st$U_ref %||% NULL
    if (is.null(U_ref)) stop("gpca apply_fn missing U_ref in method_state", call. = FALSE)
    T_new <- .new_low_rank_transform(U_ref, t(A_new))
    return(list(transforms = setNames(list(T_new), subj)))
  }

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
