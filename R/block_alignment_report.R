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
#' @param per_block_quality Logical; if TRUE, compute per-block reconstruction
#'   metrics by slicing aligned and reference matrices using
#'   `features$block_row_ranges`. Requires `aligned` and a matrix `reference`.
#' @param ... Passed to `build_alignment_features()`.
#'
#' @return An object of class `"block_alignment_report"` with components:
#' \describe{
#'   \item{features}{Output of `build_alignment_features()`.}
#'   \item{diagnostics}{Output of `feature_block_diagnostics()` on harmonized blocks.}
#'   \item{quality}{Output of `alignment_quality()` (or NULL if `aligned` is NULL).}
#'   \item{per_block_quality}{Optional per-block reconstruction summaries (or NULL).}
#' }
#'
#' @examples
#' set.seed(1)
#' b1_s1 <- alignment_feature_block(matrix(rnorm(6), 3, 2), "b1",
#'   feature_names = c("a", "b", "c")
#' )
#' b1_s2 <- alignment_feature_block(matrix(rnorm(4), 2, 2), "b1",
#'   feature_names = c("b", "c")
#' )
#' rpt <- block_alignment_report(
#'   list(s1 = list(b1 = b1_s1), s2 = list(b1 = b1_s2)),
#'   harmonize = "intersection",
#'   min_features = 2
#' )
#' names(rpt)
#'
#' @export
block_alignment_report <- function(blocks_by_subject,
                                   aligned = NULL,
                                   original = NULL,
                                   reference = NULL,
                                   quality_metrics = c("correlation", "reconstruction"),
                                   per_block_quality = FALSE,
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
  per_block <- NULL
  if (!is.null(aligned)) {
    quality <- alignment_quality(
      result = aligned,
      original = original,
      metrics = quality_metrics,
      reference = reference
    )

    if (isTRUE(per_block_quality)) {
      aligned_mats <- if (inherits(aligned, "AlignmentResult")) {
        get_aligned(aligned)
      } else if (is.list(aligned)) {
        aligned
      } else {
        stop("'aligned' must be an AlignmentResult or list of matrices", call. = FALSE)
      }

      if (length(aligned_mats) == 0L) {
        stop("No aligned matrices available for per-block quality metrics", call. = FALSE)
      }
      if (is.null(names(aligned_mats)) || any(!nzchar(names(aligned_mats)))) {
        names(aligned_mats) <- paste0("sub-", sprintf("%02d", seq_along(aligned_mats)))
      }

      ref_mat <- reference
      if (is.null(ref_mat) && inherits(aligned, "AlignmentResult")) {
        ref_mat <- get_reference(get_model(aligned))
      }
      if (is.null(ref_mat) || !.is_matrixish(ref_mat)) {
        stop("per_block_quality requires a matrix 'reference' (or a model with matrix reference_data)", call. = FALSE)
      }
      ref_mat <- as.matrix(ref_mat)

      ranges <- features$block_row_ranges
      if (!is.data.frame(ranges) || nrow(ranges) < 1L) {
        stop("Internal error: features$block_row_ranges missing or invalid", call. = FALSE)
      }

      per_block <- .compute_per_block_quality(
        aligned = aligned_mats,
        reference = ref_mat,
        ranges = ranges
      )
    }
  }

  out <- list(
    features = features,
    diagnostics = diagnostics,
    quality = quality,
    per_block_quality = per_block
  )
  class(out) <- "block_alignment_report"
  out
}

.compute_per_block_quality <- function(aligned, reference, ranges) {
  ref_dim <- dim(reference)

  by_subject <- do.call(rbind, lapply(names(aligned), function(subj) {
    x <- aligned[[subj]]
    if (!.is_matrixish(x)) x <- as.matrix(x)
    if (!identical(dim(x), ref_dim)) {
      stop(sprintf(
        "Per-block quality requires aligned matrices to match reference dims; '%s' is %d x %d but reference is %d x %d",
        subj, nrow(x), ncol(x), ref_dim[[1]], ref_dim[[2]]
      ), call. = FALSE)
    }

    do.call(rbind, lapply(seq_len(nrow(ranges)), function(i) {
      r0 <- ranges$row_start[[i]]
      r1 <- ranges$row_end[[i]]
      xb <- x[r0:r1, , drop = FALSE]
      rb <- reference[r0:r1, , drop = FALSE]

      rmse <- sqrt(mean((xb - rb)^2))
      frob <- sqrt(sum((xb - rb)^2))

      # Mean row-wise correlation with reference; rows with zero variance become NA.
      row_cors <- .row_wise_correlation(xb, rb)

      data.frame(
        subject = subj,
        block = ranges$block[[i]],
        rmse = as.numeric(rmse),
        frobenius = as.numeric(frob),
        reference_correlation = as.numeric(mean(row_cors, na.rm = TRUE)),
        stringsAsFactors = FALSE
      )
    }))
  }))

  summary <- stats::aggregate(
    by_subject[, c("rmse", "frobenius", "reference_correlation")],
    by = list(block = by_subject$block),
    FUN = function(x) mean(x, na.rm = TRUE)
  )
  names(summary) <- c("block", "mean_rmse", "mean_frobenius", "mean_reference_correlation")

  list(by_subject = by_subject, summary = summary)
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
