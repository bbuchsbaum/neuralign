#' Harmonize Matrices to a Union of IDs
#'
#' Harmonize a named list of matrices to a shared union of row/column IDs,
#' filling missing entries with `fill` (default 0). This is useful when subjects
#' have partially overlapping feature/effect vocabularies but a downstream
#' method expects all matrices to live in a common space.
#'
#' The function also returns provenance metadata describing which IDs were
#' observed in each subject, along with pairwise co-observation counts.
#'
#' @param mats Named list of matrices (or `Matrix` objects).
#' @param ids Optional named list of ID vectors corresponding to the chosen
#'   axis for each matrix. If `NULL`, IDs are inferred from `rownames()` or
#'   `colnames()` depending on `axis`.
#' @param axis Which margin to harmonize: `"rows"` or `"cols"`.
#' @param union_ids Optional character vector specifying the union vocabulary.
#'   If provided, order is used as-is (and `union_order` is ignored).
#' @param union_order Ordering used when computing a union (ignored when
#'   `union_ids` is provided): `"first_seen"` preserves first appearance order;
#'   `"sorted"` sorts IDs.
#' @param fill Fill value for missing entries (default 0).
#' @param min_coverage Minimum number of observed IDs required per subject. Any
#'   subject with fewer observed IDs is dropped.
#' @param warn_sparse_below Warn when a subject observes fewer than this
#'   fraction of union IDs (after filtering).
#'
#' @return A list with elements:
#' \describe{
#'   \item{mats}{Named list of harmonized matrices.}
#'   \item{ids}{Character vector of union IDs.}
#'   \item{obs_mask}{Named list of logical vectors (TRUE = originally present).}
#'   \item{pair_counts}{Integer matrix of co-observation counts for ID pairs.}
#'   \item{coverage}{Data frame with columns `id` and `n_subjects`.}
#'   \item{dropped}{Character vector of dropped subject IDs.}
#' }
#'
#' @examples
#' x1 <- matrix(1:6, 3, 2, dimnames = list(c("a", "b", "c"), c("v1", "v2")))
#' x2 <- matrix(1:4, 2, 2, dimnames = list(c("b", "d"), c("v1", "v2")))
#' res <- harmonize_union_fill(list(s1 = x1, s2 = x2), axis = "rows", fill = 0)
#' res$ids
#' rownames(res$mats$s1)
#' rownames(res$mats$s2)
#'
#' @family harmonize
#' @export
harmonize_union_fill <- function(mats,
                                 ids = NULL,
                                 axis = c("rows", "cols"),
                                 union_ids = NULL,
                                 union_order = c("first_seen", "sorted"),
                                 fill = 0,
                                 min_coverage = 1L,
                                 warn_sparse_below = 0.5) {
  axis <- match.arg(axis)
  union_order <- match.arg(union_order)
  min_coverage <- as.integer(min_coverage)
  if (!is.finite(min_coverage) || min_coverage < 0L) {
    stop("'min_coverage' must be a non-negative integer", call. = FALSE)
  }
  warn_sparse_below <- as.numeric(warn_sparse_below)
  if (!is.finite(warn_sparse_below) || warn_sparse_below < 0 || warn_sparse_below > 1) {
    stop("'warn_sparse_below' must be in [0, 1]", call. = FALSE)
  }

  if (!is.list(mats) || length(mats) == 0L || is.null(names(mats))) {
    stop("'mats' must be a non-empty named list", call. = FALSE)
  }

  is_matrixish <- function(x) {
    is.matrix(x) || inherits(x, "Matrix")
  }
  bad <- names(mats)[!vapply(mats, is_matrixish, logical(1))]
  if (length(bad) > 0) {
    stop("All entries of 'mats' must be matrices/Matrix; invalid: ", paste(bad, collapse = ", "), call. = FALSE)
  }

  if (is.null(ids)) {
    ids <- lapply(mats, function(m) {
      if (axis == "rows") rownames(m) else colnames(m)
    })
  }
  if (!is.list(ids)) {
    stop("'ids' must be a list of the same length as 'mats' (or NULL)", call. = FALSE)
  }

  # Named list: check for missing subjects first (more informative than length error).
  if (!is.null(names(ids))) {
    missing <- setdiff(names(mats), names(ids))
    if (length(missing) > 0) {
      stop("'ids' is missing subjects: ", paste(missing, collapse = ", "), call. = FALSE)
    }
    ids <- ids[names(mats)]
  } else {
    if (length(ids) != length(mats)) {
      stop("'ids' must be a list of the same length as 'mats' (or NULL)", call. = FALSE)
    }
    names(ids) <- names(mats)
  }

  ids <- lapply(ids, function(x) as.character(x %||% character(0)))

  # Validate ids length matches matrix margin.
  for (subj in names(mats)) {
    m <- mats[[subj]]
    idv <- ids[[subj]]
    expected <- if (axis == "rows") nrow(m) else ncol(m)
    if (length(idv) != expected) {
      stop(sprintf(
        "ID length mismatch for '%s': expected %d %s ids but got %d",
        subj, expected, axis, length(idv)
      ), call. = FALSE)
    }
    if (anyDuplicated(idv)) {
      stop(sprintf("Duplicate ids for '%s' along %s", subj, axis), call. = FALSE)
    }
  }

  if (!is.null(union_ids)) {
    union_ids <- as.character(union_ids)
    if (anyDuplicated(union_ids)) {
      stop("'union_ids' must not contain duplicates", call. = FALSE)
    }
    # Ensure provided union covers all observed ids.
    missing <- setdiff(unique(unlist(ids, use.names = FALSE)), union_ids)
    if (length(missing) > 0) {
      stop(
        "Observed ids not present in 'union_ids': ",
        paste(missing, collapse = ", "),
        call. = FALSE
      )
    }
  } else {
    if (union_order == "sorted") {
      union_ids <- sort(unique(unlist(ids, use.names = FALSE)))
    } else {
      union_ids <- character(0)
      for (idv in ids) {
        new_ids <- idv[!idv %in% union_ids]
        union_ids <- c(union_ids, new_ids)
      }
    }
  }

  union_ids <- as.character(union_ids)
  if (length(union_ids) == 0L) {
    stop("No ids remain after union construction", call. = FALSE)
  }

  n_union <- length(union_ids)

  obs_mask <- lapply(ids, function(idv) {
    mask <- union_ids %in% idv
    names(mask) <- union_ids
    mask
  })

  keep <- vapply(obs_mask, function(mask) sum(mask) >= min_coverage, logical(1))
  dropped <- names(obs_mask)[!keep]
  if (!any(keep)) {
    stop("All subjects were dropped by 'min_coverage'", call. = FALSE)
  }

  mats <- mats[keep]
  ids <- ids[keep]
  obs_mask <- obs_mask[keep]

  # Coverage warnings (after filtering).
  for (subj in names(obs_mask)) {
    frac <- sum(obs_mask[[subj]]) / n_union
    if (is.finite(frac) && frac < warn_sparse_below) {
      warning(sprintf(
        "Subject '%s': sparse coverage (%.1f%% of ids present)",
        subj, 100 * frac
      ), call. = FALSE)
    }
  }

  # Build harmonized matrices.
  out_mats <- vector("list", length(mats))
  names(out_mats) <- names(mats)

  for (subj in names(mats)) {
    m <- mats[[subj]]
    idv <- ids[[subj]]
    idx <- match(idv, union_ids)
    if (anyNA(idx)) {
      stop("Internal error: id mapping produced NA indices", call. = FALSE)
    }

    if (axis == "rows") {
      if (inherits(m, "Matrix") && identical(fill, 0)) {
        out <- Matrix::sparseMatrix(
          i = integer(0), j = integer(0), x = numeric(0),
          dims = c(n_union, ncol(m)),
          dimnames = list(union_ids, colnames(m)),
          giveCsparse = TRUE
        )
        out[idx, ] <- m
      } else {
        out <- matrix(fill,
          nrow = n_union,
          ncol = ncol(m),
          dimnames = list(union_ids, colnames(m))
        )
        out[idx, ] <- as.matrix(m)
      }
      out_mats[[subj]] <- out
    } else {
      if (inherits(m, "Matrix") && identical(fill, 0)) {
        out <- Matrix::sparseMatrix(
          i = integer(0), j = integer(0), x = numeric(0),
          dims = c(nrow(m), n_union),
          dimnames = list(rownames(m), union_ids),
          giveCsparse = TRUE
        )
        out[, idx] <- m
      } else {
        out <- matrix(fill,
          nrow = nrow(m),
          ncol = n_union,
          dimnames = list(rownames(m), union_ids)
        )
        out[, idx] <- as.matrix(m)
      }
      out_mats[[subj]] <- out
    }
  }

  # pair_counts via crossprod of mask matrix.
  mask_mat <- do.call(rbind, lapply(obs_mask, function(mask) as.integer(mask)))
  pair_counts <- crossprod(mask_mat)
  pair_counts <- matrix(
    as.integer(pair_counts),
    nrow = n_union,
    ncol = n_union,
    dimnames = list(union_ids, union_ids)
  )

  coverage <- data.frame(
    id = union_ids,
    n_subjects = as.integer(diag(pair_counts)),
    stringsAsFactors = FALSE
  )

  list(
    mats = out_mats,
    ids = union_ids,
    obs_mask = obs_mask,
    pair_counts = pair_counts,
    coverage = coverage,
    dropped = dropped
  )
}

