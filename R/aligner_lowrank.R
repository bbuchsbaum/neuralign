#' Low-Rank Alignment
#'
#' Low-rank structure plus similarity constraints via manifoldalign.
#'
#' @name aligner_lowrank
#' @keywords internal
NULL


#' Low-rank Fit Function
#'
#' @param data AlignmentData object.
#' @param reference Reference specification.
#' @param train_idx Indices of subjects to use for fitting.
#' @param ncomp Number of components.
#' @param mu Balance between low-rank and similarity (0-1).
#' @param ... Additional arguments forwarded to manifoldalign::lowrank_align().
#'
#' @return List with transforms, reference_data, etc.
#' @keywords internal
.lowrank_fit <- function(data,
                         reference = "medoid",
                         train_idx = NULL,
                         ncomp = 10L,
                         mu = 0.5,
                         ...) {
  .require_manifoldalign("low-rank")

  dots <- list(...)

  if (is.null(train_idx)) {
    train_idx <- seq_along(data@subjects)
  }
  train_data <- data[train_idx]
  data_list <- get_data_list(train_data)
  ref <- .resolve_reference(train_data, reference)

  hd <- .build_hyperdesign(train_data, transpose = FALSE)

  lr_args <- utils::modifyList(
    list(ncomp = ncomp, mu = mu),
    dots
  )

  feat_labels <- train_data@metadata[["feature_labels"]]
  if (!is.null(feat_labels)) {
    lr_args$y <- feat_labels
  }

  lr_result <- do.call(
    manifoldalign::lowrank_align,
    c(list(data = hd), lr_args)
  )

  transforms <- .extract_operators_from_scores(
    lr_result, data_list, ref$name, ref$data
  )

  heldout <- setdiff(data@subjects, names(transforms))
  if (length(heldout) > 0) {
    for (subj in heldout) {
      subj_data <- get_subject_data(data, subj)
      small_hd <- structure(
        list(subj = list(x = subj_data), ref = list(x = ref$data)),
        class = "hyperdesign"
      )
      small_args <- lr_args
      if (!is.null(feat_labels)) small_args$y <- feat_labels
      small_result <- do.call(
        manifoldalign::lowrank_align,
        c(list(data = small_hd), small_args)
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
    method_state   = list(reference = ref$name, lr_args = lr_args)
  )
}


.lowrank_capabilities <- list(
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
.register_lowrank <- function() {
  register_aligner(
    name         = "lowrank",
    fit_fn       = .lowrank_fit,
    apply_fn     = NULL,
    capabilities = .lowrank_capabilities,
    package      = "manifoldalign",
    description  = "Low-rank alignment",
    version      = "0.1.0"
  )
}
