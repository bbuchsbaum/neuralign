#' ROI/Parcel Anchor Projectors
#'
#' Build per-subject anchor projectors (anchors x native features) from ROI or
#' parcellation definitions. These projectors can be used as anatomy/parcel
#' *guidance* for functional alignment (e.g., to define shared anchor
#' coordinates when native feature counts differ).
#'
#' @name roi_anchors
NULL

.as_character_vector <- function(x) {
  if (is.null(x)) return(NULL)
  if (is.factor(x)) return(as.character(x))
  if (!is.atomic(x)) {
    stop("Expected an atomic vector", call. = FALSE)
  }
  as.character(x)
}

.normalize_projector_rows <- function(P) {
  rs <- rowSums(P)
  keep <- rs != 0
  if (any(keep)) {
    P[keep, ] <- P[keep, , drop = FALSE] / rs[keep]
  }
  P
}

#' ROI Anchor Projector
#'
#' Construct an (anchors x features) projector `P` for a single subject.
#' If `X` is a subject data matrix `(features x observations)`, then
#' `P %*% X` yields `(anchors x observations)` anchor summaries.
#'
#' `roi` can be:
#' - an atomic vector of length `n_features` giving an ROI label per feature
#' - a named list mapping ROI names to integer feature indices
#' - an (anchors x features) weight matrix (used as-is, optionally normalized)
#'
#' @param roi ROI specification (see details).
#' @param anchors Optional character vector of anchor names to include and order.
#' @param n_features Optional feature count (useful for list-of-indices input).
#' @param normalize Logical; if `TRUE`, normalize each anchor row to sum to 1.
#' @param sparse Logical; if `TRUE`, return a sparse `Matrix` when possible.
#'
#' @return An (anchors x features) projector matrix/Matrix.
#' @export
roi_anchor_projector <- function(roi,
                                 anchors = NULL,
                                 n_features = NULL,
                                 normalize = TRUE,
                                 sparse = TRUE) {
  anchors <- .as_character_vector(anchors)

  if (.is_matrixish(roi)) {
    P <- if (inherits(roi, "Matrix")) roi else as.matrix(roi)
    if (!is.null(anchors)) {
      if (is.null(rownames(P))) {
        stop("roi weight matrix must have rownames when 'anchors' is supplied", call. = FALSE)
      }
      missing <- setdiff(anchors, rownames(P))
      if (length(missing) > 0) {
        stop(
          sprintf("roi weight matrix is missing anchors: %s", paste(missing, collapse = ", ")),
          call. = FALSE
        )
      }
      P <- P[anchors, , drop = FALSE]
    } else if (!is.null(rownames(P))) {
      anchors <- rownames(P)
    } else {
      anchors <- paste0("anchor-", seq_len(nrow(P)))
      rownames(P) <- anchors
    }

    if (isTRUE(normalize)) {
      P <- .normalize_projector_rows(P)
    }
    return(P)
  }

  if (is.list(roi)) {
    if (is.null(names(roi)) || any(!nzchar(names(roi)))) {
      stop("roi list must be named by anchor", call. = FALSE)
    }
    roi <- roi[names(roi) != ""]
    roi_names <- names(roi)
    if (is.null(anchors)) {
      anchors <- roi_names
    } else {
      missing <- setdiff(anchors, roi_names)
      if (length(missing) > 0) {
        stop(
          sprintf("roi list is missing anchors: %s", paste(missing, collapse = ", ")),
          call. = FALSE
        )
      }
    }

    if (is.null(n_features)) {
      idx_all <- unlist(roi, use.names = FALSE)
      if (length(idx_all) == 0) {
        stop("Cannot infer n_features from empty roi index list; supply n_features", call. = FALSE)
      }
      n_features <- max(as.integer(idx_all), na.rm = TRUE)
    }
    n_features <- as.integer(n_features)
    if (!is.finite(n_features) || n_features < 1L) {
      stop("'n_features' must be a positive integer", call. = FALSE)
    }

    rows <- match(roi_names, anchors)
    keep <- !is.na(rows)
    roi <- roi[keep]
    roi_names <- names(roi)
    rows <- match(roi_names, anchors)

    i <- integer(0)
    j <- integer(0)
    x <- numeric(0)

    for (k in seq_along(roi_names)) {
      nm <- roi_names[[k]]
      idx <- as.integer(roi[[nm]])
      idx <- idx[is.finite(idx)]
      idx <- idx[idx >= 1L & idx <= n_features]
      idx <- unique(idx)
      if (length(idx) == 0) next
      w <- if (isTRUE(normalize)) rep(1 / length(idx), length(idx)) else rep(1, length(idx))
      i <- c(i, rep(rows[[k]], length(idx)))
      j <- c(j, idx)
      x <- c(x, w)
    }

    if (length(anchors) == 0) {
      stop("No anchors available to build projector", call. = FALSE)
    }

    if (isTRUE(sparse)) {
      P <- Matrix::sparseMatrix(
        i = i, j = j, x = x,
        dims = c(length(anchors), n_features),
        dimnames = list(anchors, NULL),
        repr = "C"
      )
      return(P)
    }

    P <- matrix(0, nrow = length(anchors), ncol = n_features, dimnames = list(anchors, NULL))
    if (length(i) > 0) {
      P[cbind(i, j)] <- x
    }
    return(P)
  }

  roi_vec <- .as_character_vector(roi)
  if (is.null(roi_vec)) {
    stop("'roi' must be a vector, list, or matrix-like object", call. = FALSE)
  }

  n_features <- length(roi_vec)
  if (is.null(anchors)) {
    anchors <- sort(unique(roi_vec[!is.na(roi_vec) & nzchar(roi_vec)]))
  }

  i <- integer(0)
  j <- integer(0)
  x <- numeric(0)

  for (k in seq_along(anchors)) {
    nm <- anchors[[k]]
    idx <- which(roi_vec == nm)
    if (length(idx) == 0) next
    w <- if (isTRUE(normalize)) rep(1 / length(idx), length(idx)) else rep(1, length(idx))
    i <- c(i, rep(k, length(idx)))
    j <- c(j, idx)
    x <- c(x, w)
  }

  if (isTRUE(sparse)) {
    return(Matrix::sparseMatrix(
      i = i, j = j, x = x,
      dims = c(length(anchors), n_features),
      dimnames = list(anchors, NULL),
      repr = "C"
    ))
  }

  P <- matrix(0, nrow = length(anchors), ncol = n_features, dimnames = list(anchors, NULL))
  if (length(i) > 0) {
    P[cbind(i, j)] <- x
  }
  P
}

#' ROI Anchor Projectors (Per Subject)
#'
#' Build a named list of per-subject ROI anchor projectors with a harmonized
#' anchor vocabulary.
#'
#' @param roi_by_subject Named list of ROI specifications (one per subject).
#' @param anchors Optional anchor vocabulary to enforce (character vector).
#' @param anchor_policy How to harmonize anchors when `anchors` is NULL:
#'   `"intersection"` keeps anchors present in all subjects; `"union"` keeps
#'   the union and fills missing anchors with zero rows.
#' @param normalize Logical; if `TRUE`, normalize each anchor row to sum to 1.
#' @param sparse Logical; if `TRUE`, return sparse matrices where possible.
#'
#' @return Named list of (anchors x features) projector matrices/Matrix objects.
#' @export
roi_anchor_projectors <- function(roi_by_subject,
                                  anchors = NULL,
                                  anchor_policy = c("intersection", "union"),
                                  normalize = TRUE,
                                  sparse = TRUE) {
  if (!is.list(roi_by_subject) || is.null(names(roi_by_subject))) {
    stop("'roi_by_subject' must be a named list keyed by subject", call. = FALSE)
  }

  anchor_policy <- match.arg(anchor_policy)
  anchors <- .as_character_vector(anchors)

  # First pass: build projectors (may have differing anchor sets)
  proj <- lapply(roi_by_subject, function(roi) {
    roi_anchor_projector(roi, anchors = anchors, normalize = normalize, sparse = sparse)
  })

  if (!is.null(anchors)) {
    # Already harmonized
    return(proj)
  }

  anchor_sets <- lapply(proj, rownames)
  if (any(vapply(anchor_sets, is.null, logical(1)))) {
    stop("All projectors must have rownames to harmonize anchors", call. = FALSE)
  }

  anchors_harmonized <- if (anchor_policy == "intersection") {
    Reduce(intersect, anchor_sets)
  } else {
    sort(unique(unlist(anchor_sets, use.names = FALSE)))
  }

  if (length(anchors_harmonized) == 0) {
    stop("No anchors remain after harmonization", call. = FALSE)
  }

  if (anchor_policy == "union") {
    res <- harmonize_union_fill(
      mats = proj,
      axis = "rows",
      union_ids = anchors_harmonized,
      fill = 0,
      min_coverage = 0L,
      warn_sparse_below = 0
    )
    out <- res$mats
    if (isTRUE(sparse)) {
      out <- lapply(out, function(P) if (inherits(P, "Matrix")) P else Matrix::Matrix(P, sparse = TRUE))
    } else {
      out <- lapply(out, as.matrix)
    }
    names(out) <- names(proj)
    return(out)
  }

  out <- lapply(proj, function(P) {
    P <- if (inherits(P, "Matrix")) P else as.matrix(P)
    have <- rownames(P)
    if (anchor_policy == "intersection") {
      return(P[anchors_harmonized, , drop = FALSE])
    }
  })
  names(out) <- names(proj)
  out
}

#' Add ROI Guidance to AlignmentData
#'
#' Convenience wrapper that builds ROI anchor projectors and stores them as
#' `"projector"` guidance channels on an `AlignmentData` object.
#'
#' @param data An `AlignmentData` object.
#' @param roi_by_subject Named list of ROI specifications (one per subject).
#' @param channel_name Name for the guidance channel (default `"roi"`).
#' @param anchors Optional anchor vocabulary to enforce (character vector).
#' @param anchor_policy How to harmonize anchors when `anchors` is NULL.
#' @param normalize Logical; if `TRUE`, normalize projector rows to sum to 1.
#' @param sparse Logical; if `TRUE`, return sparse matrices where possible.
#' @param validate Logical; if `TRUE`, validate dimensional compatibility.
#'
#' @return The modified `AlignmentData`.
#' @export
set_roi_guidance <- function(data,
                             roi_by_subject,
                             channel_name = "roi",
                             anchors = NULL,
                             anchor_policy = c("intersection", "union"),
                             normalize = TRUE,
                             sparse = TRUE,
                             validate = TRUE) {
  if (!inherits(data, "AlignmentData")) {
    stop("'data' must be an AlignmentData object", call. = FALSE)
  }
  if (!is.character(channel_name) || length(channel_name) != 1L || !nzchar(channel_name)) {
    stop("'channel_name' must be a non-empty character scalar", call. = FALSE)
  }

  anchor_policy <- match.arg(anchor_policy)
  proj <- roi_anchor_projectors(
    roi_by_subject = roi_by_subject,
    anchors = anchors,
    anchor_policy = anchor_policy,
    normalize = normalize,
    sparse = sparse
  )

  subjects <- data@subjects
  missing <- setdiff(subjects, names(proj))
  if (length(missing) > 0) {
    stop(
      sprintf("roi_by_subject is missing subjects: %s", paste(missing, collapse = ", ")),
      call. = FALSE
    )
  }

  g <- get_guidance(data)
  for (subj in subjects) {
    g_sub <- g[[subj]] %||% list()
    P <- proj[[subj]]
    g_sub[[channel_name]] <- guidance_channel(
      type = "projector",
      value = P,
      name = channel_name,
      anchors = rownames(P) %||% NULL
    )
    g[[subj]] <- g_sub
  }

  set_guidance(data, g, validate = validate)
}
