#' Block-driven Alignment Diagnostics Report
#'
#' Produce a unified diagnostics report for block-driven alignment workflows.
#' This combines:
#' \itemize{
#'   \item feature-name harmonization coverage (per block, per subject),
#'   \item rank/effective-rank diagnostics for each block and the stacked matrix,
#'   \item optional alignment quality metrics for aligned outputs.
#' }
#'
#' The report is domain-agnostic: it treats blocks as arbitrary named feature
#' matrices with metadata, and aligned matrices as generic numeric matrices.
#'
#' @param blocks_by_subject Named list keyed by subject, where each entry is a
#'   named list of `"alignment_feature_block"` objects.
#' @param aligned Optional. Either an `AlignmentResult` or a named list of
#'   aligned matrices (as accepted by `alignment_quality()`).
#' @param original Optional `AlignmentData` (or coercible) passed to
#'   `alignment_quality()` for improvement metrics.
#' @param reference Optional reference matrix used for reconstruction metrics
#'   when `aligned` is a plain list of matrices.
#' @param quality_metrics Character vector of `alignment_quality()` metrics to
#'   compute when `aligned` is provided.
#' @param ... Passed to `build_alignment_features()`.
#'
#' @return An object of class `"block_alignment_report"` with components:
#' \describe{
#'   \item{features}{Output of `build_alignment_features()`.}
#'   \item{diagnostics}{Output of `feature_block_diagnostics()` on harmonized blocks.}
#'   \item{quality}{Output of `alignment_quality()` (or NULL if `aligned` is NULL).}
#' }
#'
#' @export
block_alignment_report <- function(blocks_by_subject,
                                   aligned = NULL,
                                   original = NULL,
                                   reference = NULL,
                                   quality_metrics = c("correlation", "reconstruction"),
                                   ...) {
  features <- build_alignment_features(blocks_by_subject, ...)

  diagnostics <- feature_block_diagnostics(
    blocks = features$blocks,
    convention = features$params$convention,
    block_weights = features$params$block_weights,
    tol = features$params$rank_tol,
    include_singular_values = FALSE
  )

  quality <- NULL
  if (!is.null(aligned)) {
    quality <- alignment_quality(
      result = aligned,
      original = original,
      metrics = quality_metrics,
      reference = reference
    )
  }

  out <- list(
    features = features,
    diagnostics = diagnostics,
    quality = quality
  )
  class(out) <- "block_alignment_report"
  out
}

#' @export
print.block_alignment_report <- function(x, ...) {
  stopifnot(inherits(x, "block_alignment_report"))

  subjects <- names(x$features$matrices)
  n_subjects <- length(subjects)
  n_blocks <- if (is.data.frame(x$features$per_block)) nrow(x$features$per_block) else NA_integer_

  cat("<block_alignment_report>\n")
  cat(sprintf("Subjects: %d\n", n_subjects))
  cat(sprintf("Blocks: %s\n", if (is.na(n_blocks)) "NA" else as.character(n_blocks)))

  if (length(x$features$dropped_blocks) > 0) {
    cat(sprintf("Dropped blocks: %s\n", paste(x$features$dropped_blocks, collapse = ", ")))
  }

  if (is.data.frame(x$features$coverage) && nrow(x$features$coverage) > 0) {
    min_cov <- suppressWarnings(min(x$features$coverage$fraction_observed, na.rm = TRUE))
    if (is.finite(min_cov)) {
      cat(sprintf("Min block coverage: %.2f\n", min_cov))
    }
  }

  if (is.data.frame(x$features$identifiability) && nrow(x$features$identifiability) > 0) {
    min_ident <- suppressWarnings(min(x$features$identifiability$fraction_identified_rank, na.rm = TRUE))
    if (is.finite(min_ident)) {
      cat(sprintf("Min identified fraction (rank/transform_dim): %.2f\n", min_ident))
    }
  }

  if (is.list(x$quality) && length(x$quality) > 0) {
    if (!is.null(x$quality$mean_pairwise_correlation)) {
      cat(sprintf("Mean pairwise correlation: %.3f\n", x$quality$mean_pairwise_correlation))
    }
    if (!is.null(x$quality$mean_reconstruction_error)) {
      cat(sprintf("Mean reconstruction error: %.3f\n", x$quality$mean_reconstruction_error))
    }
  }

  invisible(x)
}

