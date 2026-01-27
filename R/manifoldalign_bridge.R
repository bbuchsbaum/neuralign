#' Bridge Utilities for manifoldalign Integration
#'
#' Data conversion and operator extraction for manifoldalign methods.
#'
#' @name manifoldalign_bridge
#' @keywords internal
NULL


#' Build a Minimal Hyperdesign from AlignmentData
#'
#' Creates a \code{hyperdesign}-like structure that manifoldalign methods
#' can consume.  neuralign stores matrices as \code{(features x observations)};
#' manifoldalign expects \code{(observations x features)} in its design slots.
#'
#' When \code{transpose = FALSE} (the default) the matrices are passed as-is,
#' treating rows (features / voxels) as samples and columns (observations /
#' conditions) as features.
#' This matches the convention used by the existing GW / FPGW adapters
#' and is appropriate for correspondence-based methods that build
#' intra-domain distance or similarity graphs over features.
#'
#' When \code{transpose = TRUE} the matrices are transposed so that rows are
#' observations and columns are features.  Use this when the manifoldalign
#' method needs a \code{y} (label) argument that corresponds to observation
#' labels (stored in \code{data@@obs_labels}).
#'
#' @param data AlignmentData object.
#' @param train_idx Integer indices of subjects to include (or NULL for all).
#' @param transpose Logical; transpose matrices before wrapping?
#' @param labels Optional labels vector to place in the design.
#'   Ignored when \code{transpose = FALSE}.
#' @param label_name Column name for labels in the design data frame.
#'
#' @return A list of class \code{"hyperdesign"} consumable by manifoldalign
#'   methods.
#'
#' @keywords internal
.build_hyperdesign <- function(data,
                               train_idx = NULL,
                               transpose = FALSE,
                               labels = NULL,
                               label_name = "label") {
  if (is.null(train_idx)) {
    train_idx <- seq_along(data@subjects)
  }
  train_data <- data[train_idx]
  data_list <- get_data_list(train_data)
  domain_names <- names(data_list)

  domains <- lapply(domain_names, function(nm) {
    X <- as.matrix(data_list[[nm]])
    if (transpose) X <- t(X)

    dom <- list(x = X)

    # Attach a design data-frame when labels are available
    if (transpose && !is.null(labels)) {
      dom$design <- data.frame(label = labels)
      names(dom$design)[1] <- label_name
    }
    dom
  })
  names(domains) <- domain_names
  structure(domains, class = "hyperdesign")
}


#' Resolve Observation Labels for manifoldalign
#'
#' Extracts a character/factor vector of observation labels from
#' \code{AlignmentData}.  Returns NULL if none are available.
#'
#' @param data AlignmentData object.
#'
#' @return A vector of labels or NULL.
#'
#' @keywords internal
.resolve_obs_labels <- function(data) {
  labs <- data@obs_labels
  if (is.null(labs)) return(NULL)


  # Atomic or factor → shared across subjects

  if (is.atomic(labs) || is.factor(labs)) {
    return(labs)
  }


  # Per-subject list → take the first subject's labels (all should match

  # for label-supervised methods that require shared observations).
  if (is.list(labs) && length(labs) > 0) {
    return(labs[[1]])
  }
  NULL
}


#' Extract Per-Subject Operators from a multiblock_biprojector
#'
#' Given a \code{multiblock_biprojector} result from manifoldalign and
#' the source \code{AlignmentData}, derive per-subject \code{(target x source)}
#' linear operators suitable for neuralign.
#'
#' The strategy depends on the method output:
#' \itemize{
#'   \item \strong{Scores-based (default)}: Extract per-subject score blocks
#'     from \code{$s}.  Compute \code{T_i = S_ref \%*\% pinv(S_i)} so that
#'     \code{T_i \%*\% X_i \approx X_ref} in the subspace.
#'   \item \strong{Transport/coupling}: Similar to GW adapter — normalise
#'     coupling rows to build a stochastic operator.
#' }
#'
#' @param mbp A multiblock_biprojector (or list with $s, $block_indices).
#' @param data_list Named list of per-subject matrices (features x obs).
#' @param ref_name Name of the reference subject.
#' @param ref_data Reference matrix (features x obs), or NULL to extract
#'   from \code{data_list}.
#' @param lambda Ridge regularisation for pseudo-inverse (default 1e-6).
#'
#' @return Named list of per-subject operator matrices.
#'
#' @keywords internal
.extract_operators_from_scores <- function(mbp,
                                           data_list,
                                           ref_name,
                                           ref_data = NULL,
                                           lambda = 1e-6) {
  S <- mbp$s  # total_points x ncomp
  n_subjects <- length(data_list)
  subject_names <- names(data_list)

  # Determine per-subject score ranges.
  # block_indices / feature_blocks stores column ranges in the concatenated

  # feature dimension, but $s is concatenated along *rows* (one per point).
  # With the no-transpose convention, each subject contributes nrow (= n_features)
  # points.
  n_features <- vapply(data_list, nrow, integer(1))
  cum <- cumsum(c(0L, n_features))

  # Extract per-subject score blocks
  score_blocks <- setNames(lapply(seq_along(subject_names), function(i) {
    rows <- (cum[i] + 1L):cum[i + 1L]
    S[rows, , drop = FALSE]
  }), subject_names)

  # Reference scores
  if (is.null(ref_data)) {
    ref_data <- data_list[[ref_name]]
  }
  S_ref <- score_blocks[[ref_name]]
  ncomp <- ncol(S_ref)

  # Build per-subject operators:  T_i = S_ref %*% pinv(S_i)
  # pinv(S_i) = solve(t(S_i) %*% S_i + lambda * I) %*% t(S_i)
  transforms <- setNames(lapply(subject_names, function(nm) {
    if (nm == ref_name) {
      return(diag(nrow(data_list[[nm]])))
    }
    S_i <- score_blocks[[nm]]
    gram <- crossprod(S_i) + lambda * diag(ncomp)
    # T_i = S_ref %*% solve(gram) %*% t(S_i)
    S_ref %*% solve(gram, t(S_i))
  }), subject_names)

  transforms
}


#' Require manifoldalign at Runtime
#'
#' Utility that stops with a clear message when manifoldalign is not installed.
#'
#' @param method Character string naming the method (for error message).
#' @keywords internal
.require_manifoldalign <- function(method) {
  if (!requireNamespace("manifoldalign", quietly = TRUE)) {
    stop(
      sprintf("Package 'manifoldalign' is required for %s alignment.", method),
      call. = FALSE
    )
  }
}


#' Standard Reference Resolution
#'
#' Resolve the reference argument to a subject name and matrix.
#'
#' @param data AlignmentData object (training subset).
#' @param reference Character or matrix.
#' @return List with \code{name} and \code{data}.
#' @keywords internal
.resolve_reference <- function(data, reference) {
  if (is.character(reference) && reference %in% c("medoid", "centroid")) {
    ref_name <- select_reference(data, method = reference)
    ref_data <- get_subject_data(data, ref_name)
  } else if (.is_matrixish(reference)) {
    ref_name <- "template"
    ref_data <- as.matrix(reference)
  } else if (is.character(reference)) {
    ref_name <- reference
    ref_data <- get_subject_data(data, reference)
  } else {
    stop("Unsupported reference type", call. = FALSE)
  }
  list(name = ref_name, data = ref_data)
}
