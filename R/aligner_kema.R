#' KEMA Alignment
#'
#' Kernel Manifold Alignment (KEMA) via manifoldalign. Produces the fitted
#' nonlinear training embeddings in a shared latent space. It does not pretend
#' that KEMA's primal-vector diagnostic is an exact out-of-sample operator.
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
#' @param solver KEMA solver. Only the truthful upstream contract, "exact", is
#'   supported.
#' @param lambda Regularisation parameter.
#' @param ... Additional arguments forwarded to manifoldalign::kema()
#'   (except \code{preproc}; neuralign owns the input orientation and forces
#'   \code{preproc=multivarious::pass()}).
#'
#' @return List with named training embeddings, reference data, spaces, and
#'   method state.
#'
#' @keywords internal
.kema_fit <- function(data,
                      reference = "medoid",
                      train_idx = NULL,
                      ncomp = 10L,
                      knn = 5L,
                      sigma = NULL,
                      u = 0.5,
                      solver = "exact",
                      lambda = 1e-2,
                      target_space = "latent",
                      fit_context = NULL,
                      provider_plan = NULL,
                      ...) {
  .ma_require_manifoldalign("KEMA")

  if (!identical(target_space, "latent")) {
    stop(
      "KEMA returns a nonlinear latent embedding; target_space='reference' ",
      "would require an exact decoder and is not supported.",
      call. = FALSE
    )
  }
  if (!identical(solver, "exact")) {
    stop(
      "KEMA currently supports only solver='exact'; other solvers do not ",
      "provide neuralign's verified training-score contract.",
      call. = FALSE
    )
  }

  dots <- list(...)
  if ("preproc" %in% names(dots)) {
    stop(
      "manifoldalign 'preproc' is not supported via neuralign adapters. ",
      "neuralign forces preproc=multivarious::pass() to preserve the declared ",
      "input orientation. ",
      "Preprocess your inputs explicitly (e.g. preprocess_alignment_data(center='rows')) before fit_alignment().",
      call. = FALSE
    )
  }

  if (is.null(train_idx)) {
    train_idx <- seq_along(data@subjects)
  }
  if (!setequal(as.integer(train_idx), seq_along(data@subjects))) {
    stop(
      "KEMA currently returns training embeddings only; subject-subset fitting ",
      "and subject cross-validation require a truthful out-of-sample kernel ",
      "extension and are not supported.",
      call. = FALSE
    )
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
  observation_counts <- vapply(data_list_train, ncol, integer(1))
  aligned <- .ma_mbp_split_scores(
    kema_result,
    domain_names = names(data_list_train),
    observation_counts = observation_counts
  )
  Z_ref <- aligned[[ref$name]]

  list(
    aligned = aligned,
    reference_data = Z_ref,
    space_from = train_data@space,
    space_to   = NULL,
    method_state = list(
      reference = ref$name,
      target_space = target_space,
      embedding_contract = "training_scores",
      out_of_sample = "unsupported",
      upstream_backend = kema_result$backend %||% NULL,
      upstream_fidelity = kema_result$fidelity %||% NULL,
      kema_args = kema_args
    )
  )
}


#' KEMA Capabilities
#' @keywords internal
.kema_capabilities <- list(
  supports_cv                = FALSE,
  cv_axes                    = character(0),
  needs_geometry             = FALSE,
  needs_design               = FALSE,
  requires_shared_features   = FALSE,
  requires_shared_observations = TRUE,
  returns_invertible         = FALSE,
  transform_type             = "embedding",
  mass_preserving            = FALSE,
  returns                    = "embedding",
  supports_new_subject       = FALSE,
  supports_new_data          = FALSE,
  reference_types            = c("subject")
)


#' Register KEMA Aligner
#' @keywords internal
.register_kema <- function() {
  register_aligner(
    name        = "kema",
    fit_fn      = .kema_fit,
    apply_fn    = NULL,
    capabilities = .kema_capabilities,
    package     = "manifoldalign",
    description = "Kernel Manifold Alignment training embeddings (KEMA)",
    version     = "0.1.0"
  )
}
