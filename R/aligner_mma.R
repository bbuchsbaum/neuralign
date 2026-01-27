#' MMA Alignment
#'
#' Multiset Manifold Alignment via manifoldalign.  Embeds multiple
#' domains into a shared spectral space and aligns them through
#' rotation synchronisation or consensus matching.
#'
#' @name aligner_mma
#' @keywords internal
NULL


#' MMA Fit Function
#'
#' @param data AlignmentData object.
#' @param reference Reference specification; "medoid" chooses the medoid
#'   subject as anchor; a subject id selects that subject.
#' @param train_idx Indices of subjects to use for fitting.
#' @param ncomp Number of embedding components.
#' @param sigma Bandwidth for nearest-neighbour graph.
#' @param match_to "reference" or "consensus".
#' @param ... Additional arguments forwarded to
#'   manifoldalign::mma_align_multiple().
#'
#' @return List with transforms, reference_data, etc.
#' @keywords internal
.mma_fit <- function(data,
                     reference = "medoid",
                     train_idx = NULL,
                     ncomp = 10L,
                     sigma = 0.73,
                     match_to = "reference",
                     ...) {
  .require_manifoldalign("MMA")

  dots <- list(...)

  if (is.null(train_idx)) {
    train_idx <- seq_along(data@subjects)
  }
  train_data <- data[train_idx]
  data_list <- get_data_list(train_data)
  train_subjects <- train_data@subjects

  ref <- .resolve_reference(train_data, reference)

  # Determine ref_idx in the training set
  ref_idx <- match(ref$name, train_subjects)
  if (is.na(ref_idx)) ref_idx <- 1L

  hd <- .build_hyperdesign(train_data, transpose = FALSE)

  mma_args <- utils::modifyList(
    list(
      ref_idx  = ref_idx,
      ncomp    = ncomp,
      sigma    = sigma,
      match_to = match_to
    ),
    dots
  )

  mma_result <- do.call(
    manifoldalign::mma_align_multiple,
    c(list(data = hd), mma_args)
  )

  transforms <- .extract_operators_from_scores(
    mma_result, data_list, ref$name, ref$data
  )

  # Held-out subjects
  heldout <- setdiff(data@subjects, names(transforms))
  if (length(heldout) > 0) {
    for (subj in heldout) {
      subj_data <- get_subject_data(data, subj)
      small_hd <- structure(
        list(subj = list(x = subj_data), ref = list(x = ref$data)),
        class = "hyperdesign"
      )
      small_args <- mma_args
      small_args$ref_idx <- 2L
      small_result <- do.call(
        manifoldalign::mma_align_multiple,
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
    method_state   = list(
      reference = ref$name,
      match_to  = match_to,
      mma_args  = mma_args
    )
  )
}


.mma_capabilities <- list(
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
.register_mma <- function() {
  register_aligner(
    name         = "mma",
    fit_fn       = .mma_fit,
    apply_fn     = NULL,
    capabilities = .mma_capabilities,
    package      = "manifoldalign",
    description  = "Multiset Manifold Alignment",
    version      = "0.1.0"
  )
}
