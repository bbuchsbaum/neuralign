#' Convert AlignedStudy to hyperdesign
#'
#' Analysis adapter: each subject/session block becomes a domain with
#' `x` (observations × shared features) and `design` (observation metadata).
#'
#' This adapter requires the suggested `multidesign` package and uses its public
#' constructors. The result therefore satisfies the actual `hyperdesign`
#' accessor contract rather than merely carrying a compatible class label.
#'
#' Conversion is lossless with respect to block matrices and observation
#' tables, but drops alignment-specific slots (model, safety, lineage). Keep
#' the [AlignedStudy] as the durable scientific object.
#'
#' @param x An [AlignedStudy].
#' @param ... Unused.
#'
#' @return A `multidesign::hyperdesign` object.
#'
#' @family aligned_study
#' @export
as_hyperdesign <- function(x, ...) {
  UseMethod("as_hyperdesign")
}

#' @export
as_hyperdesign.AlignedStudy <- function(x, ...) {
  assert_common_shared_space(x)
  if (!length(x@blocks)) {
    stop("AlignedStudy has no blocks", call. = FALSE)
  }
  if (!requireNamespace("multidesign", quietly = TRUE)) {
    stop(
      "as_hyperdesign() requires the suggested package 'multidesign'",
      call. = FALSE
    )
  }

  domains <- lapply(x@blocks, function(block) {
    multidesign::multidesign(
      x = block$values,
      y = block$observation_data
    )
  })
  names(domains) <- names(x@blocks)
  out <- multidesign::hyperdesign(domains, block_names = names(domains))
  attr(out, "shared_space_id") <- x@shared_space$id
  out
}

#' @export
as_hyperdesign.default <- function(x, ...) {
  stop("as_hyperdesign() not implemented for class ",
       paste(class(x), collapse = "/"), call. = FALSE)
}


#' Extract Subject Matrices From AlignedStudy
#'
#' @param x An [AlignedStudy].
#' @param orientation `"analysis"` (default, observations × features) or
#'   `"algorithm"` (features × observations).
#'
#' @return Named list of matrices.
#' @family aligned_study
#' @export
as_subject_matrices <- function(x, orientation = c("analysis", "algorithm")) {
  orientation <- match.arg(orientation)
  if (!inherits(x, "AlignedStudy")) {
    stop("'x' must be an AlignedStudy", call. = FALSE)
  }
  out <- lapply(x@blocks, function(b) b$values)
  if (identical(orientation, "algorithm")) {
    out <- lapply(out, to_algorithm_matrix)
  }
  out
}


#' Coerce AlignedStudy Back Toward AlignmentData
#'
#' Returns algorithm-facing [AlignmentData] (features × observations).
#' Alignment metadata (safety/lineage/shared space) is not preserved in the
#' `AlignmentData` object; use the original [AlignedStudy] for that.
#'
#' @param x An [AlignedStudy].
#' @param space Optional space label; defaults to the shared aligned space.
#' @param ... Unused.
#'
#' @return An [AlignmentData].
#' @family aligned_study
#' @export
as_alignment_data.AlignedStudy <- function(x, space = NULL, ...) {
  mats <- as_subject_matrices(x, orientation = "algorithm")
  if (is.null(space)) {
    space <- x@shared_space$id
  }
  # Build per-subject obs labels from observation_id when present
  obs_labels <- lapply(x@blocks, function(b) {
    od <- b$observation_data
    if ("observation_id" %in% names(od)) as.character(od$observation_id) else NULL
  })
  if (all(vapply(obs_labels, is.null, logical(1)))) {
    obs_labels <- NULL
  }
  AlignmentData(mats, space = space, obs_labels = obs_labels)
}
