#' PARROT Alignment
#'
#' Position-Aware Random-walk Regularised Optimal Transport via manifoldalign.
#' Computes a soft transport plan between domains using anchor correspondences
#' and random-walk graph structure.
#'
#' @name aligner_parrot
#' @keywords internal
NULL


#' PARROT Fit Function
#'
#' @param data AlignmentData object.
#' @param reference Reference specification.
#' @param train_idx Indices of subjects to use for fitting.
#' @param ncomp Number of embedding components (NULL for auto).
#' @param sigma Random-walk restart probability.
#' @param lambda Consistency regularisation.
#' @param tau Sinkhorn entropy regularisation.
#' @param alpha Feature vs graph weight.
#' @param gamma Cross-graph mixing.
#' @param solver "sinkhorn" or "exact".
#' @param ... Additional arguments forwarded to manifoldalign::parrot().
#'
#' @return List with transforms, reference_data, etc.
#' @keywords internal
.parrot_fit <- function(data,
                        reference = "medoid",
                        train_idx = NULL,
                        ncomp = NULL,
                        sigma = 0.15,
                        lambda = 0.1,
                        tau = 0.05,
                        alpha = 0.2,
                        gamma = 0.1,
                        solver = "sinkhorn",
                        ...) {
  .require_manifoldalign("PARROT")

  dots <- list(...)

  if (is.null(train_idx)) {
    train_idx <- seq_along(data@subjects)
  }
  train_data <- data[train_idx]
  data_list <- get_data_list(train_data)
  train_subjects <- train_data@subjects

  ref <- .resolve_reference(train_data, reference)

  # PARROT requires anchors.  obs_labels provide observation-level anchors.
  # In the no-transpose convention (features-as-rows), we need feature-level
  # anchor correspondences.  Use obs_labels as the anchor column if both
  # subjects share the same observations — manifoldalign's parrot expects
  # an "anchors" argument identifying which column in the design holds
  # shared identifiers.
  obs_labels <- .resolve_obs_labels(data)
  if (is.null(obs_labels)) {
    stop(
      "PARROT requires anchor correspondences; ",
      "set obs_labels on AlignmentData to provide them.",
      call. = FALSE
    )
  }

  hd <- .build_hyperdesign(train_data, transpose = FALSE)

  parrot_args <- utils::modifyList(
    list(
      ncomp   = ncomp,
      sigma   = sigma,
      lambda  = lambda,
      tau     = tau,
      alpha   = alpha,
      gamma   = gamma,
      solver  = solver
    ),
    dots
  )
  # Anchors: the design column (populated by .build_hyperdesign when labels
  # are present) or a column name string.  For the no-transpose case we
  # populate anchors via column names matching across subjects.
  parrot_args$anchors <- obs_labels

  # PARROT is pairwise.  Align each subject to reference.
  transforms <- list()

  for (subj in train_subjects) {
    Xi <- data_list[[subj]]
    if (subj == ref$name) {
      transforms[[subj]] <- diag(nrow(Xi))
      next
    }

    pair_hd <- structure(
      list(
        subj = list(x = Xi),
        ref  = list(x = ref$data)
      ),
      class = "hyperdesign"
    )
    pair_result <- do.call(
      manifoldalign::parrot,
      c(list(data = pair_hd), parrot_args)
    )

    # Extract transport plan as operator (same pattern as GW)
    P <- pair_result$transport_plan %||%
      tryCatch(pair_result$alignment_matrix, error = function(e) NULL)
    if (!is.null(P)) {
      transforms[[subj]] <- .coupling_to_operator(P)
    } else {
      # Fallback to score-based extraction
      small_dl <- list(subj = Xi, ref = ref$data)
      transforms[[subj]] <- .extract_operators_from_scores(
        pair_result, small_dl, "ref", ref$data
      )[["subj"]]
    }
  }

  # Held-out subjects
  heldout <- setdiff(data@subjects, names(transforms))
  if (length(heldout) > 0) {
    for (subj in heldout) {
      subj_data <- get_subject_data(data, subj)
      pair_hd <- structure(
        list(subj = list(x = subj_data), ref = list(x = ref$data)),
        class = "hyperdesign"
      )
      pair_result <- do.call(
        manifoldalign::parrot,
        c(list(data = pair_hd), parrot_args)
      )
      P <- pair_result$transport_plan %||%
        tryCatch(pair_result$alignment_matrix, error = function(e) NULL)
      if (!is.null(P)) {
        transforms[[subj]] <- .coupling_to_operator(P)
      } else {
        small_dl <- list(subj = subj_data, ref = ref$data)
        transforms[[subj]] <- .extract_operators_from_scores(
          pair_result, small_dl, "ref", ref$data
        )[["subj"]]
      }
    }
  }

  list(
    transforms     = transforms,
    reference_data = ref$data,
    space_from     = train_data@space,
    space_to       = train_data@space,
    method_state   = list(reference = ref$name, parrot_args = parrot_args)
  )
}


.parrot_capabilities <- list(
  supports_cv                  = TRUE,
  cv_axes                      = c("subject"),
  needs_geometry               = FALSE,
  needs_design                 = FALSE,
  requires_shared_features     = FALSE,
  requires_shared_observations = TRUE,
  returns_invertible           = FALSE,
  transform_type               = "ot",
  mass_preserving              = FALSE,
  returns                      = "operator",
  supports_new_subject         = TRUE,
  supports_new_data            = TRUE,
  reference_types              = c("subject", "medoid", "template")
)


#' @keywords internal
.register_parrot <- function() {
  register_aligner(
    name         = "parrot",
    fit_fn       = .parrot_fit,
    apply_fn     = NULL,
    capabilities = .parrot_capabilities,
    package      = "manifoldalign",
    description  = "PARROT optimal-transport alignment",
    version      = "0.1.0"
  )
}
