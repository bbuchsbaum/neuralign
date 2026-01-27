#' KEMA Alignment
#'
#' Kernel Manifold Alignment (KEMA) via manifoldalign. Produces per-subject
#' operators by extracting score-based embeddings and constructing linear
#' maps to a reference subject's feature space.
#'
#' @name aligner_kema
#' @keywords internal
NULL


#' KEMA Fit Function
#'
#' @param data AlignmentData object.
#' @param reference Reference specification (subject id, "medoid", "consensus").
#' @param train_idx Indices of subjects to use for fitting.
#' @param ncomp Number of latent components.
#' @param knn Number of nearest neighbours for graph construction.
#' @param sigma Kernel bandwidth (NULL for auto).
#' @param u Trade-off between manifold and class alignment (0-1).
#' @param solver KEMA solver: "regression" or "exact".
#' @param lambda Regularisation parameter.
#' @param ... Additional arguments forwarded to manifoldalign::kema().
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
                      lambda = 1e-4,
                      ...) {
  .require_manifoldalign("KEMA")

  dots <- list(...)

  # --- Subset to training subjects ---
  if (is.null(train_idx)) {
    train_idx <- seq_along(data@subjects)
  }
  train_data <- data[train_idx]
  train_subjects <- train_data@subjects
  data_list <- get_data_list(train_data)

  # --- Resolve reference ---
  ref <- .resolve_reference(train_data, reference)

  # --- Build hyperdesign (no transpose: features-as-points) ---
  hd <- .build_hyperdesign(train_data, transpose = FALSE)

  # --- Resolve observation labels for the `y` argument ---
  # KEMA needs a label vector, one per "sample" (here: per feature/voxel).
  # If obs_labels are available they apply to observations (columns), not
  # features (rows).  For the no-transpose convention we rely on KEMA's
  # unsupervised capability (u=1 fully manifold-based) when no feature-level
  # labels are supplied.
  obs_labels <- .resolve_obs_labels(data)

  # --- Assemble arguments ---
  kema_args <- utils::modifyList(
    list(
      ncomp = ncomp,
      knn   = knn,
      sigma = sigma,
      u     = u,
      solver = solver,
      lambda = lambda
    ),
    dots
  )

  # Add y only if labels are present and correspond to observations
  # (i.e. the hyperdesign rows).  With no-transpose, rows are features;
  # only include y if feature-level labels are supplied via metadata.
  feat_labels <- train_data@metadata[["feature_labels"]]
  if (!is.null(feat_labels)) {
    kema_args$y <- feat_labels
  }

  kema_result <- do.call(
    manifoldalign::kema,
    c(list(data = hd), kema_args)
  )

  # --- Extract per-subject operators ---
  transforms <- .extract_operators_from_scores(
    kema_result, data_list, ref$name, ref$data
  )

  # --- Handle held-out subjects ---
  all_subjects <- data@subjects
  heldout <- setdiff(all_subjects, names(transforms))
  if (length(heldout) > 0) {
    for (subj in heldout) {
      subj_data <- get_subject_data(data, subj)
      # Re-fit a small 2-domain problem (subject + reference)
      small_hd <- structure(
        list(
          subj = list(x = subj_data),
          ref  = list(x = ref$data)
        ),
        class = "hyperdesign"
      )
      small_kema_args <- kema_args
      small_kema_args$y <- NULL
      if (!is.null(feat_labels)) {
        small_kema_args$y <- feat_labels
      }
      small_result <- do.call(
        manifoldalign::kema,
        c(list(data = small_hd), small_kema_args)
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
      ncomp     = ncomp,
      knn       = knn,
      sigma     = sigma,
      u         = u,
      solver    = solver,
      lambda    = lambda,
      reference = ref$name,
      kema_args = kema_args
    )
  )
}


#' KEMA Capabilities
#' @keywords internal
.kema_capabilities <- list(
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


#' Register KEMA Aligner
#' @keywords internal
.register_kema <- function() {
  register_aligner(
    name        = "kema",
    fit_fn      = .kema_fit,
    apply_fn    = NULL,
    capabilities = .kema_capabilities,
    package     = "manifoldalign",
    description = "Kernel Manifold Alignment (KEMA)",
    version     = "0.1.0"
  )
}
