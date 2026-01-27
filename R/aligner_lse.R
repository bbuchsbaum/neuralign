#' Linear Similarity Embedding Alignment
#'
#' Linear Similarity Embedding (LSE) via manifoldalign.  Learns a linear
#' projection that preserves a target similarity structure.
#'
#' LSE is a single-domain method.  The adapter runs it independently on
#' each subject (or on the concatenated data) and uses the resulting
#' projection weights to build per-subject operators.
#'
#' @name aligner_lse
#' @keywords internal
NULL


#' LSE Fit Function
#'
#' @param data AlignmentData object.
#' @param reference Reference specification.
#' @param train_idx Indices of subjects to use for fitting.
#' @param ncomp Number of embedding components.
#' @param sigma_P Gaussian bandwidth for similarity ("auto" for auto).
#' @param maxit Maximum optimisation iterations.
#' @param tol Convergence tolerance.
#' @param ... Additional arguments forwarded to
#'   manifoldalign::linear_sim_embed().
#'
#' @return List with transforms, reference_data, etc.
#' @keywords internal
.lse_fit <- function(data,
                     reference = "medoid",
                     train_idx = NULL,
                     ncomp = 10L,
                     sigma_P = "auto",
                     maxit = 500L,
                     tol = 1e-6,
                     ...) {
  .require_manifoldalign("LSE")

  dots <- list(...)

  if (is.null(train_idx)) {
    train_idx <- seq_along(data@subjects)
  }
  train_data <- data[train_idx]
  data_list <- get_data_list(train_data)
  ref <- .resolve_reference(train_data, reference)

  lse_args <- utils::modifyList(
    list(ncomp = ncomp, sigma_P = sigma_P, maxit = maxit, tol = tol),
    dots
  )

  # LSE works on a single matrix.  Concatenate features (rows) across
  # subjects, compute a shared projection, then derive per-subject
  # operators from the score subspace.
  X_concat <- do.call(rbind, data_list)

  lse_result <- do.call(
    manifoldalign::linear_sim_embed,
    c(list(X = X_concat), lse_args)
  )

  # The result has $scores (n_total x ncomp).  Use score-based extraction.
  mbp_like <- list(s = lse_result$scores)
  transforms <- .extract_operators_from_scores(
    mbp_like, data_list, ref$name, ref$data
  )

  # Held-out subjects: project using the learned weights
  heldout <- setdiff(data@subjects, names(transforms))
  if (length(heldout) > 0) {
    for (subj in heldout) {
      subj_data <- get_subject_data(data, subj)
      new_scores <- stats::predict(lse_result, subj_data)
      # Use Procrustes from new scores to reference scores
      ref_idx <- match(ref$name, names(data_list))
      n_feat <- vapply(data_list, nrow, integer(1))
      cum <- cumsum(c(0L, n_feat))
      S_ref <- lse_result$scores[(cum[ref_idx] + 1L):cum[ref_idx + 1L], ,
                                  drop = FALSE]
      nk <- ncol(new_scores)
      gram <- crossprod(new_scores) + 1e-6 * diag(nk)
      transforms[[subj]] <- S_ref %*% solve(gram, t(new_scores))
    }
  }

  list(
    transforms     = transforms,
    reference_data = ref$data,
    space_from     = train_data@space,
    space_to       = train_data@space,
    method_state   = list(reference = ref$name, lse_args = lse_args)
  )
}


.lse_capabilities <- list(
  supports_cv                  = TRUE,
  cv_axes                      = c("subject"),
  needs_geometry               = FALSE,
  needs_design                 = FALSE,
  requires_shared_features     = TRUE,
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
.register_lse <- function() {
  register_aligner(
    name         = "lse",
    fit_fn       = .lse_fit,
    apply_fn     = NULL,
    capabilities = .lse_capabilities,
    package      = "manifoldalign",
    description  = "Linear Similarity Embedding",
    version      = "0.1.0"
  )
}
