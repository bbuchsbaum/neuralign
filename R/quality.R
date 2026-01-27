#' Compute Alignment Quality Metrics
#'
#' Compute various quality metrics for an alignment result.
#'
#' @param result An AlignmentResult object, or a list of aligned matrices.
#' @param original Optional AlignmentData with original (unaligned) data
#'   for computing improvement metrics.
#' @param metrics Character vector of metrics to compute:
#'   \itemize{
#'     \item "correlation" - Mean pairwise correlation
#'     \item "reconstruction" - Reconstruction vs reference (row-wise correlation,
#'       RMSE, and Frobenius residual)
#'     \item "variance" - Explained variance by alignment
#'     \item "isc" - Inter-subject correlation (time-locked)
#'   }
#' @param reference Optional reference matrix for reconstruction metrics.
#'
#' @return Named list of quality metrics.
#'
#' @examples
#' \dontrun{
#' result <- fit_alignment(data, method = "procrustes")
#' quality <- alignment_quality(result)
#' print(quality$mean_pairwise_correlation)
#' }
#'
#' @export
alignment_quality <- function(result,
                              original = NULL,
                              metrics = c("correlation"),
                              reference = NULL) {
  # Extract aligned data
  if (inherits(result, "AlignmentResult")) {
    aligned <- get_aligned(result)
    if (is.null(reference)) {
      reference <- get_reference(get_model(result))
    }
  } else if (is.list(result)) {
    aligned <- result
  } else {
    stop("'result' must be an AlignmentResult or list of matrices", call. = FALSE)
  }

  if (length(aligned) == 0) {
    stop("No aligned data to compute quality", call. = FALSE)
  }

  if (is.null(names(aligned)) || any(!nzchar(names(aligned)))) {
    names(aligned) <- paste0("sub-", sprintf("%02d", seq_along(aligned)))
  }

  quality <- list()

  # Compute requested metrics
  if ("correlation" %in% metrics) {
    quality <- c(quality, .compute_correlation_metrics(aligned))
  }

  if ("reconstruction" %in% metrics && !is.null(reference)) {
    quality <- c(quality, .compute_reconstruction_metrics(aligned, reference))
  }

  if ("variance" %in% metrics) {
    quality <- c(quality, .compute_variance_metrics(aligned))
  }

  if ("isc" %in% metrics) {
    quality <- c(quality, .compute_isc_metrics(aligned))
  }

  # Compute improvement over original if provided
  if (!is.null(original)) {
    if (!inherits(original, "AlignmentData")) {
      original <- as_alignment_data(original)
    }
    quality <- c(quality, .compute_improvement_metrics(
      aligned, get_data_list(original)
    ))
  }

  quality
}


#' Compute Correlation Metrics
#' @keywords internal
.compute_correlation_metrics <- function(aligned) {
  n <- length(aligned)

  if (n < 2) {
    return(list(
      mean_pairwise_correlation = NA_real_,
      pairwise_correlations = numeric(0)
    ))
  }

  # Compute pairwise correlations
  pairs <- combn(names(aligned), 2)
  cors <- apply(pairs, 2, function(p) {
    x <- as.matrix(aligned[[p[1]]])
    y <- as.matrix(aligned[[p[2]]])

    # Mean row-wise correlation
    mean(.row_wise_correlation(x, y), na.rm = TRUE)
  })

  names(cors) <- apply(pairs, 2, paste, collapse = "-")

  list(
    mean_pairwise_correlation = mean(cors, na.rm = TRUE),
    sd_pairwise_correlation = sd(cors, na.rm = TRUE),
    min_pairwise_correlation = min(cors, na.rm = TRUE),
    max_pairwise_correlation = max(cors, na.rm = TRUE),
    pairwise_correlations = cors
  )
}


.row_wise_correlation <- function(x, y) {
  x <- as.matrix(x)
  y <- as.matrix(y)

  if (!identical(dim(x), dim(y))) {
    stop(
      sprintf(
        "Row-wise correlation requires identical dims; x is %d x %d but y is %d x %d",
        nrow(x), ncol(x), nrow(y), ncol(y)
      ),
      call. = FALSE
    )
  }

  if (ncol(x) < 2L) {
    return(rep(NA_real_, nrow(x)))
  }

  mx <- rowMeans(x)
  my <- rowMeans(y)
  xc <- x - mx
  yc <- y - my

  num <- rowSums(xc * yc)
  den <- sqrt(rowSums(xc^2) * rowSums(yc^2))
  out <- num / den
  out[!is.finite(out)] <- NA_real_
  out
}


#' Compute Reconstruction Metrics
#' @keywords internal
.compute_reconstruction_metrics <- function(aligned, reference) {
  reference <- as.matrix(reference)
  ref_dim <- dim(reference)

  for (subj in names(aligned)) {
    x <- as.matrix(aligned[[subj]])
    if (!identical(dim(x), ref_dim)) {
      stop(sprintf(
        "Reconstruction metrics require aligned matrices to match reference dims; '%s' is %d x %d but reference is %d x %d",
        subj, nrow(x), ncol(x), ref_dim[[1]], ref_dim[[2]]
      ), call. = FALSE)
    }
  }

  # For each subject, compute correlation with reference
  recon_cors <- vapply(aligned, function(x) {
    mean(.row_wise_correlation(x, reference), na.rm = TRUE)
  }, numeric(1))

  # Also compute reconstruction error
  recon_errors <- vapply(aligned, function(x) {
    x <- as.matrix(x)
    sqrt(mean((x - reference)^2))
  }, numeric(1))

  # Frobenius-norm reconstruction residuals (sqrt(sum((x - ref)^2))).
  recon_frobenius <- vapply(aligned, function(x) {
    x <- as.matrix(x)
    sqrt(sum((x - reference)^2))
  }, numeric(1))

  list(
    mean_reference_correlation = mean(recon_cors, na.rm = TRUE),
    reference_correlations = recon_cors,
    mean_reconstruction_error = mean(recon_errors, na.rm = TRUE),
    reconstruction_errors = recon_errors,
    mean_reconstruction_frobenius = mean(recon_frobenius, na.rm = TRUE),
    reconstruction_frobenius = recon_frobenius
  )
}


#' Compute Variance Metrics
#' @keywords internal
.compute_variance_metrics <- function(aligned) {
  if (!requireNamespace("abind", quietly = TRUE)) {
    warning("Package 'abind' required for variance metrics")
    return(list(
      total_variance = NA_real_,
      between_subject_variance = NA_real_,
      within_subject_variance = NA_real_,
      variance_ratio = NA_real_
    ))
  }

  # Stack all aligned data
  stacked <- do.call(abind::abind, c(
    lapply(aligned, function(x) as.matrix(x)),
    list(along = 3)
  ))

  # Total variance across subjects
  total_var <- var(as.vector(stacked))

  # Mean within each feature across subjects
  n_features <- dim(stacked)[1]
  n_obs <- dim(stacked)[2]
  n_subjects <- dim(stacked)[3]

  # Between-subject variance (for each feature-timepoint)
  between_var <- mean(apply(stacked, c(1, 2), var))

  # Within-subject variance (time series variance averaged)
  within_var <- mean(vapply(
    seq_len(n_subjects),
    function(s) var(as.vector(stacked[, , s])),
    numeric(1)
  ))

  list(
    total_variance = total_var,
    between_subject_variance = between_var,
    within_subject_variance = within_var,
    variance_ratio = if (within_var > 0) between_var / within_var else NA_real_
  )
}


#' Compute Inter-Subject Correlation
#' @keywords internal
.compute_isc_metrics <- function(aligned) {
  n <- length(aligned)

  if (n < 2) {
    return(list(isc = NA_real_))
  }

  # Stack aligned data
  data_list <- lapply(aligned, as.matrix)
  n_features <- nrow(data_list[[1]])
  n_obs <- ncol(data_list[[1]])

  # ISC: correlate each subject's time series with mean of others
  isc_values <- matrix(NA, n, n_features)

  for (i in seq_len(n)) {
    # Mean of all others
    others <- data_list[-i]
    mean_others <- compute_centroid(others)

    # Correlate each feature
    for (j in seq_len(n_features)) {
      isc_values[i, j] <- cor(data_list[[i]][j, ], mean_others[j, ])
    }
  }

  # Average ISC across subjects and features
  mean_isc <- mean(isc_values, na.rm = TRUE)

  # ISC per feature
  feature_isc <- colMeans(isc_values, na.rm = TRUE)

  list(
    mean_isc = mean_isc,
    isc_per_subject = rowMeans(isc_values, na.rm = TRUE),
    isc_per_feature = feature_isc
  )
}


#' Compute Improvement Over Original
#' @keywords internal
.compute_improvement_metrics <- function(aligned, original) {
  # Get metrics for original data
  common_subjects <- intersect(names(aligned), names(original))

  if (length(common_subjects) < 2) {
    return(list(improvement = NA_real_))
  }

  aligned_sub <- aligned[common_subjects]
  original_sub <- original[common_subjects]

  # Pairwise correlation improvement
  aligned_cor <- .compute_correlation_metrics(aligned_sub)
  original_cor <- .compute_correlation_metrics(original_sub)

  improvement <- aligned_cor$mean_pairwise_correlation -
    original_cor$mean_pairwise_correlation

  list(
    original_mean_correlation = original_cor$mean_pairwise_correlation,
    aligned_mean_correlation = aligned_cor$mean_pairwise_correlation,
    correlation_improvement = improvement,
    relative_improvement = if (abs(original_cor$mean_pairwise_correlation) > 0)
      improvement / abs(original_cor$mean_pairwise_correlation) else NA_real_
  )
}


#' Print Quality Summary
#'
#' Print a formatted summary of alignment quality metrics.
#'
#' @param quality Quality metrics list from alignment_quality.
#'
#' @export
print_quality_summary <- function(quality) {
  cat("Alignment Quality Summary\n")
  cat("=========================\n\n")

  if (!is.null(quality$mean_pairwise_correlation)) {
    cat("Pairwise Correlation:\n")
    cat(sprintf("  Mean: %.4f\n", quality$mean_pairwise_correlation))
    if (!is.null(quality$sd_pairwise_correlation)) {
      cat(sprintf("  SD:   %.4f\n", quality$sd_pairwise_correlation))
    }
    if (!is.null(quality$min_pairwise_correlation)) {
      cat(sprintf("  Min:  %.4f\n", quality$min_pairwise_correlation))
      cat(sprintf("  Max:  %.4f\n", quality$max_pairwise_correlation))
    }
    cat("\n")
  }

  if (!is.null(quality$mean_reference_correlation)) {
    cat("Reference Correlation:\n")
    cat(sprintf("  Mean: %.4f\n", quality$mean_reference_correlation))
    cat("\n")
  }

  if (!is.null(quality$mean_isc)) {
    cat("Inter-Subject Correlation:\n")
    cat(sprintf("  Mean ISC: %.4f\n", quality$mean_isc))
    cat("\n")
  }

  if (!is.null(quality$correlation_improvement)) {
    cat("Improvement:\n")
    cat(sprintf("  Original correlation: %.4f\n",
      quality$original_mean_correlation
    ))
    cat(sprintf("  Aligned correlation:  %.4f\n",
      quality$aligned_mean_correlation
    ))
    cat(sprintf("  Improvement:          %.4f (%.1f%%)\n",
      quality$correlation_improvement,
      quality$relative_improvement * 100
    ))
  }

  invisible(quality)
}
