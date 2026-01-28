#' Feature Block Helpers
#'
#' Lightweight helpers for constructing, harmonizing, and stacking alignment
#' feature matrices in a way that can be shared across downstream packages.
#'
#' @name feature_blocks
#' @family feature_blocks
NULL

.is_scalar_number <- function(x) {
  is.numeric(x) && length(x) == 1L && is.finite(x)
}

.validate_feature_names <- function(feature_names, n_rows, block_name) {
  if (is.null(feature_names)) return(NULL)
  feature_names <- as.character(feature_names)
  if (length(feature_names) != n_rows) {
    stop(
      sprintf(
        "feature_names length mismatch for block '%s': got %d, need %d",
        block_name, length(feature_names), n_rows
      ),
      call. = FALSE
    )
  }
  dups <- feature_names[duplicated(feature_names)]
  if (length(dups) > 0) {
    stop(
      sprintf(
        "feature_names for block '%s' contain duplicates: %s",
        block_name, paste(unique(dups), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  feature_names
}

#' Create an Alignment Feature Block
#'
#' @family feature_blocks
#'
#' @param x Matrix of features for a block (features x observations).
#' @param name Block name (character scalar).
#' @param weight Non-negative scalar weight (applied as `sqrt(weight)` in stacking).
#' @param feature_names Optional character vector naming the rows of `x`.
#' @param meta Optional free-form list of metadata carried with the block.
#'   neuralign does not interpret these fields, but downstream packages may use
#'   them for provenance (e.g., \code{meta$source_type}, \code{meta$requires_independence}).
#'
#' @return An object of class `"alignment_feature_block"`.
#' @export
alignment_feature_block <- function(x,
                                    name,
                                    weight = 1,
                                    feature_names = NULL,
                                    meta = NULL) {
  if (!.is_matrixish(x)) {
    stop("'x' must be matrix-like", call. = FALSE)
  }
  x <- as.matrix(x)
  if (!is.character(name) || length(name) != 1L || !nzchar(name)) {
    stop("'name' must be a non-empty character scalar", call. = FALSE)
  }
  if (!.is_scalar_number(weight) || weight < 0) {
    stop("'weight' must be a single non-negative number", call. = FALSE)
  }
  if (!is.null(meta) && !is.list(meta)) {
    stop("'meta' must be NULL or a list", call. = FALSE)
  }

  feature_names <- .validate_feature_names(feature_names, nrow(x), name)

  structure(
    list(
      x = x,
      name = name,
      weight = as.numeric(weight),
      feature_names = feature_names,
      meta = meta
    ),
    class = "alignment_feature_block"
  )
}

.is_feature_block <- function(x) {
  inherits(x, "alignment_feature_block") &&
    is.list(x) &&
    !is.null(x$x) &&
    !is.null(x$name) &&
    !is.null(x$weight)
}

.validate_blocks_by_subject <- function(blocks_by_subject) {
  if (!is.list(blocks_by_subject) || length(blocks_by_subject) < 1L) {
    stop("'blocks_by_subject' must be a non-empty named list", call. = FALSE)
  }
  if (is.null(names(blocks_by_subject)) || any(!nzchar(names(blocks_by_subject)))) {
    stop("'blocks_by_subject' must be a named list keyed by subject", call. = FALSE)
  }

  subj_blocks <- lapply(blocks_by_subject, function(bl) {
    if (!is.list(bl) || length(bl) < 1L) {
      stop("Each subject must provide a non-empty list of blocks", call. = FALSE)
    }
    if (is.null(names(bl))) {
      stop("Each subject's block list must be named by block name", call. = FALSE)
    }
    if (!all(vapply(bl, .is_feature_block, logical(1)))) {
      stop("All blocks must be alignment_feature_block objects", call. = FALSE)
    }
    bl
  })

  common_block_names <- Reduce(intersect, lapply(subj_blocks, names))
  dropped_missing <- setdiff(unique(unlist(lapply(subj_blocks, names))), common_block_names)

  list(
    subjects = names(blocks_by_subject),
    common_block_names = common_block_names,
    dropped_missing = dropped_missing
  )
}

.is_named_numeric_vector <- function(x) {
  is.numeric(x) && length(x) > 0L && !is.null(names(x)) && all(nzchar(names(x)))
}

.normalize_weights <- function(w, normalize = c("median", "mean", "none")) {
  normalize <- match.arg(normalize)
  if (normalize == "none") return(w)
  ok <- is.finite(w) & w > 0
  if (!any(ok)) return(w)

  center <- if (normalize == "median") stats::median(w[ok]) else mean(w[ok])
  if (!is.finite(center) || center <= 0) return(w)
  w / center
}

.validate_block_weights <- function(block_weights,
                                   subjects = NULL,
                                   arg = "block_weights") {
  if (is.null(block_weights)) return(invisible(TRUE))

  if (.is_named_numeric_vector(block_weights)) {
    if (any(!is.finite(block_weights)) || any(block_weights < 0)) {
      stop(sprintf("'%s' must contain finite non-negative values", arg), call. = FALSE)
    }
    return(invisible(TRUE))
  }

  if (!is.list(block_weights) || length(block_weights) < 1L) {
    stop(
      sprintf("'%s' must be a named numeric vector, a named list of such vectors, or NULL", arg),
      call. = FALSE
    )
  }
  if (is.null(names(block_weights)) || any(!nzchar(names(block_weights)))) {
    stop(sprintf("'%s' list must be named by subject", arg), call. = FALSE)
  }
  if (!is.null(subjects)) {
    missing <- setdiff(subjects, names(block_weights))
    if (length(missing) > 0) {
      stop(
        sprintf("'%s' is missing subjects: %s", arg, paste(missing, collapse = ", ")),
        call. = FALSE
      )
    }
    block_weights <- block_weights[subjects]
  }

  bad <- names(block_weights)[!vapply(block_weights, .is_named_numeric_vector, logical(1))]
  if (length(bad) > 0) {
    stop(
      sprintf("'%s' list entries must be named numeric vectors; invalid subjects: %s", arg, paste(bad, collapse = ", ")),
      call. = FALSE
    )
  }
  for (subj in names(block_weights)) {
    w <- block_weights[[subj]]
    if (any(!is.finite(w)) || any(w < 0)) {
      stop(
        sprintf("'%s' for subject '%s' must contain finite non-negative values", arg, subj),
        call. = FALSE
      )
    }
  }

  invisible(TRUE)
}

.combine_weight_vectors <- function(a, b) {
  if (is.null(a)) return(b)
  if (is.null(b)) return(a)

  out_names <- union(names(a), names(b))
  out <- setNames(rep(1, length(out_names)), out_names)
  out[names(a)] <- out[names(a)] * a
  out[names(b)] <- out[names(b)] * b
  out
}

.resolve_block_weights_for_subject <- function(block_weights, subject) {
  if (is.null(block_weights)) return(NULL)
  if (is.numeric(block_weights)) return(block_weights)
  block_weights[[subject]]
}

#' Suggest Feature Block Weights
#'
#' @family feature_blocks
#'
#' Suggest additional block multipliers to mitigate dominance by large blocks
#' when stacking heterogeneous correspondence signals via \code{\link{stack_feature_blocks}}.
#'
#' Note that \code{stack_feature_blocks()} scales each block by \code{sqrt(weight)}.
#' Therefore, weights correspond to objective weights in a sum-of-squares sense.
#'
#' @param blocks_by_subject Either:
#'   \itemize{
#'     \item A list of \code{"alignment_feature_block"} objects (single subject), or
#'     \item A named list of subjects, each containing such a list.
#'   }
#' @param method Weighting heuristic:
#'   \itemize{
#'     \item \code{"equalize_rows"}: weight \eqn{\propto 1 / n_{\mathrm{rows}}}
#'     \item \code{"equalize_frobenius"}: weight \eqn{\propto 1 / ||X||_F^2}
#'     \item \code{"equalize_rms"}: weight \eqn{\propto 1 / \mathrm{mean}(X^2)}
#'   }
#' @param per Whether to suggest weights globally across subjects (\code{"global"})
#'   or separately per subject (\code{"per_subject"}).
#' @param normalize How to normalize the returned weights for readability:
#'   \code{"median"} (default) rescales so the median positive weight is 1;
#'   \code{"mean"} rescales so the mean positive weight is 1; \code{"none"}
#'   leaves raw proportional weights.
#'
#' @return If \code{per = "global"}, a named numeric vector. If
#'   \code{per = "per_subject"}, a named list keyed by subject, each containing
#'   a named numeric vector.
#' @export
suggest_block_weights <- function(blocks_by_subject,
                                  method = c("equalize_rows", "equalize_frobenius", "equalize_rms"),
                                  per = c("global", "per_subject"),
                                  normalize = c("median", "mean", "none")) {
  method <- match.arg(method)
  per <- match.arg(per)
  normalize <- match.arg(normalize)

  if (!is.list(blocks_by_subject) || length(blocks_by_subject) < 1L) {
    stop("'blocks_by_subject' must be a non-empty list", call. = FALSE)
  }

  if (all(vapply(blocks_by_subject, .is_feature_block, logical(1)))) {
    # Single subject: ensure blocks are named.
    bl <- blocks_by_subject
    if (is.null(names(bl)) || any(!nzchar(names(bl)))) {
      bl <- setNames(bl, vapply(bl, function(b) b$name, character(1)))
    }
    blocks_by_subject <- list(subject = bl)
  }

  meta <- .validate_blocks_by_subject(blocks_by_subject)
  subjects <- meta$subjects
  block_names <- meta$common_block_names
  if (length(block_names) < 1L) {
    stop("No common blocks available for weight suggestion", call. = FALSE)
  }

  compute_one <- function(blocks) {
    denom <- vapply(block_names, function(bname) {
      b <- blocks[[bname]]
      x <- as.matrix(b$x)
      if (any(!is.finite(x))) {
        stop(sprintf("Block '%s' contains non-finite values", bname), call. = FALSE)
      }
      if (nrow(x) < 1L || ncol(x) < 1L) {
        stop(sprintf("Block '%s' has invalid dimensions", bname), call. = FALSE)
      }

      if (method == "equalize_rows") {
        as.numeric(nrow(x))
      } else if (method == "equalize_frobenius") {
        sum(x^2)
      } else {
        mean(x^2)
      }
    }, numeric(1))

    w <- rep(0, length(block_names))
    names(w) <- block_names
    ok <- is.finite(denom) & denom > 0
    w[ok] <- 1 / denom[ok]
    w <- .normalize_weights(w, normalize = normalize)
    w
  }

  if (per == "per_subject") {
    out <- lapply(subjects, function(subj) compute_one(blocks_by_subject[[subj]]))
    names(out) <- subjects
    return(out)
  }

  # Global: use a robust typical denominator across subjects.
  denom_by_subj <- lapply(subjects, function(subj) {
    blocks <- blocks_by_subject[[subj]]
    vapply(block_names, function(bname) {
      x <- as.matrix(blocks[[bname]]$x)
      if (any(!is.finite(x))) {
        stop(sprintf("Block '%s' contains non-finite values", bname), call. = FALSE)
      }
      if (nrow(x) < 1L || ncol(x) < 1L) {
        stop(sprintf("Block '%s' has invalid dimensions", bname), call. = FALSE)
      }

      if (method == "equalize_rows") {
        as.numeric(nrow(x))
      } else if (method == "equalize_frobenius") {
        sum(x^2)
      } else {
        mean(x^2)
      }
    }, numeric(1))
  })
  denom_mat <- do.call(rbind, denom_by_subj)
  colnames(denom_mat) <- block_names

  denom <- apply(denom_mat, 2L, stats::median)
  w <- rep(0, length(block_names))
  names(w) <- block_names
  ok <- is.finite(denom) & denom > 0
  w[ok] <- 1 / denom[ok]
  .normalize_weights(w, normalize = normalize)
}

#' Stack Feature Blocks
#'
#' @family feature_blocks
#'
#' Vertically stack (rbind) multiple feature blocks after scaling each by
#' `sqrt(weight)`. This produces a single matrix ready for alignment methods
#' such as Procrustes.
#'
#' @param blocks List of `"alignment_feature_block"` objects.
#' @param block_weights Optional named numeric vector providing additional
#'   multipliers per block name (applied multiplicatively to each block's
#'   internal weight).
#'
#' @return A matrix (stacked features x observations).
#' @export
stack_feature_blocks <- function(blocks, block_weights = NULL) {
  if (!is.list(blocks) || length(blocks) < 1L) {
    stop("'blocks' must be a non-empty list of feature blocks", call. = FALSE)
  }
  if (!all(vapply(blocks, .is_feature_block, logical(1)))) {
    stop("'blocks' must contain only alignment_feature_block objects", call. = FALSE)
  }

  if (!is.null(block_weights)) {
    if (!is.numeric(block_weights) || is.null(names(block_weights))) {
      stop("'block_weights' must be a named numeric vector", call. = FALSE)
    }
  }

  scaled <- lapply(blocks, function(b) {
    w_extra <- if (!is.null(block_weights) && b$name %in% names(block_weights)) {
      block_weights[[b$name]]
    } else {
      1
    }
    if (!.is_scalar_number(w_extra) || w_extra < 0) {
      stop(sprintf("Invalid block_weights for block '%s'", b$name), call. = FALSE)
    }
    w <- sqrt(b$weight * w_extra)
    x <- b$x
    if (w != 1) x <- w * x

    if (!is.null(b$feature_names)) {
      rownames(x) <- paste0(b$name, ":", b$feature_names)
    } else if (!is.null(rownames(x))) {
      rownames(x) <- paste0(b$name, ":", rownames(x))
    }
    x
  })

  Reduce(rbind, scaled)
}

.block_feature_names <- function(block) {
  if (!is.null(block$feature_names)) return(as.character(block$feature_names))
  if (!is.null(rownames(block$x))) return(as.character(rownames(block$x)))
  NULL
}

#' Harmonize Feature Blocks Across Subjects
#'
#' @family feature_blocks
#'
#' For each block name, intersect feature names across subjects, subset and
#' reorder rows so that each subject has a consistent feature ordering. Blocks
#' missing from any subject are dropped. Blocks with fewer than `min_features`
#' common features are dropped with a warning.
#'
#' @param blocks_by_subject Named list of subjects, each containing a named list
#'   of `"alignment_feature_block"` objects.
#' @param min_features Minimum number of shared features required to keep a
#'   block.
#'
#' @return A harmonized `blocks_by_subject` list with the same structure.
#' @export
harmonize_feature_blocks <- function(blocks_by_subject, min_features = 2L) {
  if (!is.list(blocks_by_subject) || length(blocks_by_subject) < 1L) {
    stop("'blocks_by_subject' must be a non-empty named list", call. = FALSE)
  }
  if (is.null(names(blocks_by_subject)) || any(!nzchar(names(blocks_by_subject)))) {
    stop("'blocks_by_subject' must be a named list keyed by subject", call. = FALSE)
  }
  min_features <- as.integer(min_features)
  if (!is.finite(min_features) || min_features < 1L) {
    stop("'min_features' must be a positive integer", call. = FALSE)
  }

  # Validate blocks and collect common block names
  subj_blocks <- lapply(blocks_by_subject, function(bl) {
    if (!is.list(bl) || length(bl) < 1L) {
      stop("Each subject must provide a non-empty list of blocks", call. = FALSE)
    }
    if (is.null(names(bl))) {
      stop("Each subject's block list must be named by block name", call. = FALSE)
    }
    if (!all(vapply(bl, .is_feature_block, logical(1)))) {
      stop("All blocks must be alignment_feature_block objects", call. = FALSE)
    }
    bl
  })

  common_block_names <- Reduce(intersect, lapply(subj_blocks, names))
  dropped_missing <- setdiff(unique(unlist(lapply(subj_blocks, names))), common_block_names)
  if (length(dropped_missing) > 0) {
    warning(
      sprintf(
        "Dropping blocks not present for all subjects: %s",
        paste(dropped_missing, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  out <- blocks_by_subject
  # Drop blocks not present for all subjects
  for (subj in names(out)) {
    out[[subj]] <- out[[subj]][common_block_names]
  }

  for (bname in common_block_names) {
    feat_lists <- lapply(subj_blocks, function(bl) .block_feature_names(bl[[bname]]))
    if (any(vapply(feat_lists, is.null, logical(1)))) {
      stop(
        sprintf("Cannot harmonize block '%s': feature names are missing", bname),
        call. = FALSE
      )
    }
    common_feats <- Reduce(intersect, feat_lists)
    if (length(common_feats) < min_features) {
      warning(
        sprintf(
          "Dropping block '%s': only %d common features (< %d)",
          bname, length(common_feats), min_features
        ),
        call. = FALSE
      )
      for (subj in names(out)) out[[subj]][[bname]] <- NULL
      next
    }

    for (subj in names(out)) {
      block <- out[[subj]][[bname]]
      feats <- .block_feature_names(block)
      idx <- match(common_feats, feats)
      block$x <- as.matrix(block$x)[idx, , drop = FALSE]
      block$feature_names <- common_feats
      out[[subj]][[bname]] <- block
    }
  }

  out
}

.validate_stackable_blocks <- function(blocks, subject) {
  ncols <- vapply(blocks, function(b) ncol(as.matrix(b$x)), integer(1))
  if (length(unique(ncols)) > 1L) {
    stop(
      sprintf(
        "Subject '%s' has blocks with differing column counts: %s",
        subject,
        paste(sprintf("%s=%d", names(ncols), ncols), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Build Alignment Matrices from Feature Blocks
#'
#' @family feature_blocks
#'
#' Construct per-subject alignment matrices by harmonizing and stacking feature
#' blocks. This is useful for block-driven alignment workflows where different
#' correspondence signals (blocks) are combined into a single matrix per
#' subject.
#'
#' @param blocks_by_subject Named list keyed by subject, where each entry is a
#'   named list of `"alignment_feature_block"` objects.
#' @param harmonize Harmonization strategy for feature names within each block:
#'   `"intersection"` keeps only features common to all subjects (via
#'   `harmonize_feature_blocks()`); `"union_fill"` harmonizes to the union of
#'   feature names within each block, filling missing rows with `fill`.
#'   Blocks missing from any subject are dropped in both modes.
#' @param min_features Minimum number of features required to keep a block. For
#'   `"intersection"`, this is the size of the common intersection. For
#'   `"union_fill"`, blocks are dropped when any subject contributes fewer than
#'   `min_features` observed features.
#' @param block_weights Optional additional multipliers per block name (as in
#'   `stack_feature_blocks()`). Either a named numeric vector (applied to all
#'   subjects) or a named list keyed by subject, each containing a named numeric
#'   vector.
#' @param suggest_weights Optional weighting heuristic to apply before stacking.
#'   If provided, this computes suggested weights from the harmonized blocks and
#'   combines them multiplicatively with any explicit `block_weights`. One of
#'   `"equalize_rows"`, `"equalize_frobenius"`, or `"equalize_rms"` (see
#'   `suggest_block_weights()`).
#' @param obs_crossfit Logical; set TRUE when blocks are constructed within an
#'   observation-axis crossfit workflow (suppresses independence warnings).
#' @param check_independence Logical; if TRUE (default), warn when any block has
#'   `meta$requires_independence=TRUE` and `obs_crossfit=FALSE`.
#' @param check_identifiability Logical; if TRUE (default), compute stacked
#'   numeric/effective rank per subject and warn when the stacked matrix does
#'   not fully identify the transform (numeric_rank < transform_dim).
#' @param convention Convention for interpreting the transform dimension in the
#'   identifiability check: `"right"` means transforms act on columns; `"left"`
#'   means transforms act on rows.
#' @param rank_tol Relative tolerance for numeric rank computation as a fraction
#'   of the largest singular value.
#' @param fill Fill value used for `"union_fill"` (default 0).
#' @param union_order Ordering for union construction when `harmonize="union_fill"`
#'   and no explicit union is provided: `"first_seen"` or `"sorted"`.
#' @param warn_sparse_below Warn when a subject observes fewer than this
#'   fraction of union features within a block (only used for `"union_fill"`).
#' @param emit_warnings Logical; if FALSE, suppress warnings emitted during
#'   harmonization while still returning them in the output.
#'
#' @return A list with elements:
#' \describe{
#'   \item{matrices}{Named list of per-subject stacked matrices.}
#'   \item{blocks}{Harmonized per-subject block lists used to build matrices.}
#'   \item{per_block}{Data frame summarizing retained blocks and feature counts.}
#'   \item{block_row_ranges}{Data frame giving stacked row ranges for each block
#'     (block, row_start, row_end, n_rows).}
#'   \item{coverage}{Data frame of per-subject block coverage counts
#'     (n_total/n_observed/fraction_observed) based on the harmonized feature
#'     vocabulary.}
#'   \item{identifiability}{Data frame summarizing stacked transform dimension
#'     and numeric/effective rank per subject (NULL if `check_identifiability=FALSE`).}
#'   \item{params}{List of parameters used to construct the output.}
#'   \item{dropped_blocks}{Character vector of dropped block names.}
#'   \item{warnings}{Character vector of warning messages encountered.}
#' }
#'
#' @export
build_alignment_features <- function(blocks_by_subject,
                                     harmonize = c("intersection", "union_fill"),
                                     min_features = 2L,
                                     block_weights = NULL,
                                     suggest_weights = NULL,
                                     obs_crossfit = FALSE,
                                     check_independence = TRUE,
                                     check_identifiability = TRUE,
                                     convention = c("right", "left"),
                                     rank_tol = sqrt(.Machine$double.eps),
                                     fill = 0,
                                     union_order = c("first_seen", "sorted"),
                                     warn_sparse_below = 0.5,
                                     emit_warnings = TRUE) {
  harmonize <- match.arg(harmonize)
  union_order <- match.arg(union_order)
  convention <- match.arg(convention)

  min_features <- as.integer(min_features)
  if (!is.finite(min_features) || min_features < 1L) {
    stop("'min_features' must be a positive integer", call. = FALSE)
  }

  .validate_block_weights(block_weights, arg = "block_weights")

  check_independence <- isTRUE(check_independence)
  obs_crossfit <- isTRUE(obs_crossfit)
  check_identifiability <- isTRUE(check_identifiability)
  if (check_identifiability) {
    if (!.is_scalar_number(rank_tol) || rank_tol <= 0 || rank_tol >= 1) {
      stop("'rank_tol' must be a single number in (0, 1)", call. = FALSE)
    }
  }

  warnings <- character(0)
  capture_warnings <- function(expr) {
    withCallingHandlers(
      expr,
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        if (!isTRUE(emit_warnings)) {
          invokeRestart("muffleWarning")
        }
      }
    )
  }

  if (harmonize == "intersection") {
    blocks_h <- capture_warnings(harmonize_feature_blocks(blocks_by_subject, min_features = min_features))
    if (length(blocks_h) < 1L || length(blocks_h[[1L]]) < 1L) {
      stop("No feature blocks remain after harmonization", call. = FALSE)
    }
  } else {
    meta <- .validate_blocks_by_subject(blocks_by_subject)
    subjects <- meta$subjects
    common_block_names <- meta$common_block_names
    dropped_missing <- meta$dropped_missing

    if (length(dropped_missing) > 0) {
      capture_warnings(warning(
        sprintf(
          "Dropping blocks not present for all subjects: %s",
          paste(dropped_missing, collapse = ", ")
        ),
        call. = FALSE
      ))
    }

    blocks_h <- blocks_by_subject
    for (subj in subjects) {
      blocks_h[[subj]] <- blocks_h[[subj]][common_block_names]
    }

    # Union-fill per block.
    for (bname in common_block_names) {
      ids <- lapply(subjects, function(subj) .block_feature_names(blocks_h[[subj]][[bname]]))
      names(ids) <- subjects
      if (any(vapply(ids, is.null, logical(1)))) {
        stop(
          sprintf("Cannot harmonize block '%s': feature names are missing", bname),
          call. = FALSE
        )
      }
      mats <- lapply(subjects, function(subj) {
        b <- blocks_h[[subj]][[bname]]
        m <- as.matrix(b$x)
        rownames(m) <- ids[[subj]]
        m
      })
      names(mats) <- subjects

      harm <- capture_warnings(harmonize_union_fill(
        mats = mats,
        ids = ids,
        axis = "rows",
        union_order = union_order,
        fill = fill,
        min_coverage = 0L,
        warn_sparse_below = warn_sparse_below
      ))

      n_obs <- vapply(harm$obs_mask, function(mask) sum(mask), integer(1))
      if (any(n_obs < min_features)) {
        capture_warnings(warning(
          sprintf(
            "Dropping block '%s': one or more subjects have < %d observed features",
            bname, min_features
          ),
          call. = FALSE
        ))
        for (subj in subjects) blocks_h[[subj]][[bname]] <- NULL
        next
      }

      for (subj in subjects) {
        block <- blocks_h[[subj]][[bname]]
        block$x <- as.matrix(harm$mats[[subj]])
        block$feature_names <- harm$ids
        blocks_h[[subj]][[bname]] <- block
      }
    }

    if (length(blocks_h) < 1L || length(blocks_h[[1L]]) < 1L) {
      stop("No feature blocks remain after harmonization", call. = FALSE)
    }
  }

  # Stack per subject.
  subjects <- names(blocks_h)
  for (subj in subjects) {
    .validate_stackable_blocks(blocks_h[[subj]], subject = subj)
  }

  .validate_block_weights(block_weights, subjects = subjects, arg = "block_weights")
  if (!is.null(suggest_weights)) {
    if (!is.character(suggest_weights) || length(suggest_weights) != 1L || is.na(suggest_weights)) {
      stop("'suggest_weights' must be a non-NA character scalar", call. = FALSE)
    }
    suggested <- suggest_block_weights(blocks_h, method = suggest_weights, per = "global")
    if (is.null(block_weights)) {
      block_weights <- suggested
    } else if (is.numeric(block_weights)) {
      block_weights <- .combine_weight_vectors(block_weights, suggested)
    } else {
      block_weights <- lapply(subjects, function(subj) {
        .combine_weight_vectors(block_weights[[subj]], suggested)
      })
      names(block_weights) <- subjects
    }
    .validate_block_weights(block_weights, subjects = subjects, arg = "block_weights")
  }

  mats_by_subj <- lapply(subjects, function(subj) {
    w_subj <- .resolve_block_weights_for_subject(block_weights, subj)
    stack_feature_blocks(blocks_h[[subj]], block_weights = w_subj)
  })
  names(mats_by_subj) <- subjects

  dims <- vapply(mats_by_subj, function(m) paste0(nrow(m), "x", ncol(m)), character(1))
  if (length(unique(dims)) != 1L) {
    pieces <- paste0(names(dims), "=", dims)
    stop(
      sprintf(
        "Subjects have differing stacked matrix dimensions: %s",
        paste(pieces, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  block_names <- names(blocks_h[[1L]])
  requires_independence <- feature_block_requires_independence(blocks_h)
  requires_independence <- requires_independence[block_names]
  per_block <- data.frame(
    block = block_names,
    requires_independence = as.logical(requires_independence),
    n_features = vapply(block_names, function(bname) {
      length(.block_feature_names(blocks_h[[subjects[[1L]]]][[bname]]))
    }, integer(1)),
    stringsAsFactors = FALSE
  )

  row_end <- cumsum(per_block$n_features)
  row_start <- c(1L, utils::head(row_end, -1L) + 1L)
  block_row_ranges <- data.frame(
    block = per_block$block,
    row_start = as.integer(row_start),
    row_end = as.integer(row_end),
    n_rows = as.integer(per_block$n_features),
    stringsAsFactors = FALSE
  )

  dropped_blocks <- setdiff(
    unique(unlist(lapply(blocks_by_subject, names))),
    block_names
  )

  coverage <- do.call(rbind, lapply(subjects, function(subj) {
    do.call(rbind, lapply(block_names, function(bname) {
      vocab <- .block_feature_names(blocks_h[[subj]][[bname]])
      orig_feats <- .block_feature_names(blocks_by_subject[[subj]][[bname]])
      n_total <- length(vocab)
      n_observed <- if (is.null(orig_feats)) NA_integer_ else sum(vocab %in% orig_feats)
      data.frame(
        subject = subj,
        block = bname,
        n_total = as.integer(n_total),
        n_observed = as.integer(n_observed),
        fraction_observed = if (is.finite(n_observed) && n_total > 0L) {
          n_observed / n_total
        } else {
          NA_real_
        },
        stringsAsFactors = FALSE
      )
    }))
  }))
  rownames(coverage) <- NULL

  params <- list(
    harmonize = harmonize,
    min_features = min_features,
    block_weights = block_weights,
    obs_crossfit = obs_crossfit,
    check_independence = check_independence,
    check_identifiability = check_identifiability,
    convention = convention,
    rank_tol = rank_tol,
    fill = fill,
    union_order = union_order,
    warn_sparse_below = warn_sparse_below
  )

  if (check_independence && !obs_crossfit && any(requires_independence, na.rm = TRUE)) {
    flagged <- names(requires_independence)[which(requires_independence)]
    capture_warnings(warning(
      paste0(
        "Blocks marked meta$requires_independence=TRUE: ",
        paste(flagged, collapse = ", "),
        ". Ensure these signals were estimated on data independent of downstream evaluation ",
        "(e.g., via observation-axis crossfit). Set obs_crossfit=TRUE when appropriate, ",
        "or check_independence=FALSE to disable this warning."
      ),
      call. = FALSE
    ))
  }

  identifiability <- NULL
  if (check_identifiability) {
    nrow_by_subj <- vapply(subjects, function(subj) nrow(mats_by_subj[[subj]]), integer(1))
    ncol_by_subj <- vapply(subjects, function(subj) ncol(mats_by_subj[[subj]]), integer(1))
    transform_dim <- if (convention == "left") nrow_by_subj else ncol_by_subj

    numeric_rank <- integer(length(subjects))
    effective_rank <- numeric(length(subjects))
    for (i in seq_along(subjects)) {
      d <- .matrix_singular_values(mats_by_subj[[subjects[[i]]]])
      numeric_rank[[i]] <- .numeric_rank_from_singular_values(d, tol = rank_tol)
      effective_rank[[i]] <- .effective_rank_from_singular_values(d)
    }

    identifiability <- data.frame(
      subject = subjects,
      nrow = nrow_by_subj,
      ncol = ncol_by_subj,
      transform_dim = as.integer(transform_dim),
      numeric_rank = as.integer(numeric_rank),
      effective_rank = as.numeric(effective_rank),
      stringsAsFactors = FALSE
    )
    identifiability$fraction_identified_rank <- ifelse(
      identifiability$transform_dim > 0L,
      identifiability$numeric_rank / identifiability$transform_dim,
      NA_real_
    )
    identifiability$fraction_identified_effective <- ifelse(
      identifiability$transform_dim > 0L,
      identifiability$effective_rank / identifiability$transform_dim,
      NA_real_
    )
    identifiability$convention <- convention

    if (any(identifiability$numeric_rank < identifiability$transform_dim, na.rm = TRUE)) {
      bad <- identifiability[identifiability$numeric_rank < identifiability$transform_dim, , drop = FALSE]
      pieces <- paste0(bad$subject, " (", bad$numeric_rank, "/", bad$transform_dim, ")")
      capture_warnings(warning(
        paste0(
          "Under-identified transform: numeric_rank < transform_dim (convention='",
          convention,
          "') for ",
          paste(pieces, collapse = ", "),
          ". Downstream inference should be restricted to the identified subspace; ",
          "use feature_block_diagnostics() to inspect singular values."
        ),
        call. = FALSE
      ))
    }
  }

  list(
    matrices = mats_by_subj,
    blocks = blocks_h,
    per_block = per_block,
    block_row_ranges = block_row_ranges,
    coverage = coverage,
    identifiability = identifiability,
    params = params,
    dropped_blocks = dropped_blocks,
    warnings = unique(warnings)
  )
}

#' Feature Block Independence Requirements
#'
#' @family feature_blocks
#'
#' Summarize which feature blocks declare `meta$requires_independence = TRUE`.
#' This is a domain-agnostic mechanism for downstream packages to tag blocks
#' that must be computed on data independent of downstream evaluation (e.g.,
#' via observation-axis crossfitting).
#'
#' neuralign does not enforce a particular crossfitting policy; this helper
#' only reads the metadata flags.
#'
#' @param blocks Either:
#'   \itemize{
#'     \item A named list of `"alignment_feature_block"` objects, or
#'     \item A `blocks_by_subject` list (subject -> named list of blocks).
#'   }
#'
#' @return A named logical vector keyed by block name.
#' @export
feature_block_requires_independence <- function(blocks) {
  if (!is.list(blocks) || length(blocks) < 1L) {
    stop("'blocks' must be a non-empty list", call. = FALSE)
  }

  if (all(vapply(blocks, .is_feature_block, logical(1)))) {
    if (is.null(names(blocks)) || any(!nzchar(names(blocks)))) {
      blocks <- setNames(blocks, vapply(blocks, function(b) b$name, character(1)))
    }
    out <- vapply(blocks, function(b) {
      isTRUE(b$meta$requires_independence %||% FALSE)
    }, logical(1))
    return(out)
  }

  if (is.null(names(blocks)) || any(!nzchar(names(blocks)))) {
    stop(
      "'blocks' must be either a list of feature blocks, or a named list keyed by subject",
      call. = FALSE
    )
  }

  subjects <- names(blocks)
  block_names <- unique(unlist(lapply(blocks, names)))
  out <- setNames(rep(FALSE, length(block_names)), block_names)

  for (subj in subjects) {
    bl <- blocks[[subj]]
    if (!is.list(bl) || length(bl) < 1L) {
      stop(sprintf("Subject '%s' must provide a non-empty list of blocks", subj), call. = FALSE)
    }
    if (is.null(names(bl)) || any(!nzchar(names(bl)))) {
      stop(sprintf("Subject '%s' block list must be named by block", subj), call. = FALSE)
    }
    if (!all(vapply(bl, .is_feature_block, logical(1)))) {
      stop(sprintf("Subject '%s' contains non-feature-block entries", subj), call. = FALSE)
    }

    for (bname in names(bl)) {
      b <- bl[[bname]]
      if (isTRUE(b$meta$requires_independence %||% FALSE)) {
        out[[bname]] <- TRUE
      }
    }
  }

  out
}

.effective_rank_from_singular_values <- function(d) {
  d <- as.numeric(d)
  if (length(d) == 0) return(0)
  if (any(!is.finite(d))) return(NA_real_)
  d2 <- d^2
  total <- sum(d2)
  if (total <= 0) return(0)
  p <- d2 / total
  p <- p[p > 0]
  exp(-sum(p * log(p)))
}

.matrix_singular_values <- function(x) {
  x <- as.matrix(x)
  if (any(!is.finite(x))) {
    stop("Diagnostics require finite numeric values", call. = FALSE)
  }
  svd(x, nu = 0, nv = 0)$d
}

.numeric_rank_from_singular_values <- function(d, tol) {
  d <- as.numeric(d)
  if (length(d) == 0) return(0L)
  if (all(d == 0)) return(0L)
  as.integer(sum(d > (max(d) * tol)))
}

.feature_block_diagnostics_one <- function(blocks,
                                          convention,
                                          block_weights,
                                          tol,
                                          include_singular_values) {
  convention <- match.arg(convention, c("left", "right"))
  if (!.is_scalar_number(tol) || tol <= 0 || tol >= 1) {
    stop("'tol' must be a single number in (0, 1)", call. = FALSE)
  }

  per_block <- lapply(blocks, function(b) {
    w_extra <- if (!is.null(block_weights) && b$name %in% names(block_weights)) {
      block_weights[[b$name]]
    } else {
      1
    }
    if (!.is_scalar_number(w_extra) || w_extra < 0) {
      stop(sprintf("Invalid block_weights for block '%s'", b$name), call. = FALSE)
    }
    w_total <- b$weight * w_extra
    x <- b$x
    if (w_total != 1) x <- sqrt(w_total) * x

    d <- .matrix_singular_values(x)
    numeric_rank <- .numeric_rank_from_singular_values(d, tol = tol)
    effective_rank <- .effective_rank_from_singular_values(d)
    transform_dim <- if (convention == "left") nrow(x) else ncol(x)

    list(
      block = b$name,
      nrow = nrow(x),
      ncol = ncol(x),
      weight = as.numeric(w_total),
      transform_dim = as.integer(transform_dim),
      numeric_rank = as.integer(numeric_rank),
      effective_rank = as.numeric(effective_rank),
      fraction_identified_rank = if (transform_dim > 0) numeric_rank / transform_dim else NA_real_,
      fraction_identified_effective = if (transform_dim > 0) effective_rank / transform_dim else NA_real_,
      singular_values = if (isTRUE(include_singular_values)) d else NULL
    )
  })

  per_block_df <- data.frame(
    block = vapply(per_block, `[[`, character(1), "block"),
    nrow = vapply(per_block, `[[`, integer(1), "nrow"),
    ncol = vapply(per_block, `[[`, integer(1), "ncol"),
    weight = vapply(per_block, `[[`, numeric(1), "weight"),
    transform_dim = vapply(per_block, `[[`, integer(1), "transform_dim"),
    numeric_rank = vapply(per_block, `[[`, integer(1), "numeric_rank"),
    effective_rank = vapply(per_block, `[[`, numeric(1), "effective_rank"),
    fraction_identified_rank = vapply(per_block, `[[`, numeric(1), "fraction_identified_rank"),
    fraction_identified_effective = vapply(per_block, `[[`, numeric(1), "fraction_identified_effective"),
    stringsAsFactors = FALSE
  )

  stacked <- stack_feature_blocks(blocks, block_weights = block_weights)
  d_stacked <- .matrix_singular_values(stacked)
  numeric_rank_stacked <- .numeric_rank_from_singular_values(d_stacked, tol = tol)
  effective_rank_stacked <- .effective_rank_from_singular_values(d_stacked)
  transform_dim_stacked <- if (convention == "left") nrow(stacked) else ncol(stacked)

  out <- list(
    convention = convention,
    per_block = per_block_df,
    stacked = list(
      nrow = nrow(stacked),
      ncol = ncol(stacked),
      transform_dim = as.integer(transform_dim_stacked),
      numeric_rank = as.integer(numeric_rank_stacked),
      effective_rank = as.numeric(effective_rank_stacked),
      fraction_identified_rank = if (transform_dim_stacked > 0) {
        numeric_rank_stacked / transform_dim_stacked
      } else {
        NA_real_
      },
      fraction_identified_effective = if (transform_dim_stacked > 0) {
        effective_rank_stacked / transform_dim_stacked
      } else {
        NA_real_
      },
      singular_values = if (isTRUE(include_singular_values)) d_stacked else NULL
    )
  )

  class(out) <- "feature_block_diagnostics"
  out
}

#' Feature Block Identifiability Diagnostics
#'
#' @family feature_blocks
#'
#' Compute per-block and stacked rank/effective-rank summaries for feature block
#' workflows. This is useful when alignment is driven by correspondence signals
#' (task betas, contrast betas, anchor fingerprints): the singular value spectrum
#' of the (weighted) correspondence matrix indicates how much of the transform
#' space is identified.
#'
#' @param blocks Either:
#'   \itemize{
#'     \item A list of `"alignment_feature_block"` objects (single subject), or
#'     \item A named list of subjects, each containing such a list.
#'   }
#' @param convention Convention for interpreting the transform dimension:
#'   `"left"` means transforms act on rows, `"right"` means transforms act on
#'   columns (see `procrustes_rotation()`).
#' @param block_weights Optional named numeric vector of additional multipliers
#'   per block name (as in `stack_feature_blocks()`).
#' @param tol Relative tolerance for numeric rank computation as a fraction of
#'   the largest singular value.
#' @param include_singular_values Logical; if TRUE, include singular value
#'   vectors in the return object.
#'
#' @return For a single subject, an object of class `"feature_block_diagnostics"`.
#'   For multiple subjects, a named list of `"feature_block_diagnostics"` objects
#'   with class `"feature_block_diagnostics_by_subject"`.
#'
#' @examples
#' set.seed(1)
#' blocks <- list(
#'   alignment_feature_block(matrix(rnorm(12), 3, 4), name = "b1", feature_names = paste0("f", 1:3)),
#'   alignment_feature_block(matrix(rnorm(8), 2, 4), name = "b2", feature_names = paste0("g", 1:2))
#' )
#' diag <- feature_block_diagnostics(blocks, convention = "right", include_singular_values = TRUE)
#' diag$stacked$numeric_rank
#'
#' @export
feature_block_diagnostics <- function(blocks,
                                      convention = c("left", "right"),
                                      block_weights = NULL,
                                      tol = sqrt(.Machine$double.eps),
                                      include_singular_values = FALSE) {
  convention <- match.arg(convention)

  if (!is.list(blocks) || length(blocks) < 1L) {
    stop("'blocks' must be a non-empty list", call. = FALSE)
  }
  if (!is.null(block_weights)) {
    if (!is.numeric(block_weights) || is.null(names(block_weights))) {
      stop("'block_weights' must be a named numeric vector", call. = FALSE)
    }
  }

  if (all(vapply(blocks, .is_feature_block, logical(1)))) {
    return(.feature_block_diagnostics_one(
      blocks = blocks,
      convention = convention,
      block_weights = block_weights,
      tol = tol,
      include_singular_values = include_singular_values
    ))
  }

  if (is.null(names(blocks)) || any(!nzchar(names(blocks)))) {
    stop(
      "'blocks' must be either a list of feature blocks, or a named list keyed by subject",
      call. = FALSE
    )
  }

  out <- lapply(names(blocks), function(subj) {
    bl <- blocks[[subj]]
    if (!is.list(bl) || length(bl) < 1L) {
      stop(sprintf("Subject '%s' must provide a non-empty list of blocks", subj), call. = FALSE)
    }
    if (!all(vapply(bl, .is_feature_block, logical(1)))) {
      stop(
        sprintf("Subject '%s' contains non-feature-block entries", subj),
        call. = FALSE
      )
    }
    .feature_block_diagnostics_one(
      blocks = bl,
      convention = convention,
      block_weights = block_weights,
      tol = tol,
      include_singular_values = include_singular_values
    )
  })
  names(out) <- names(blocks)
  class(out) <- c("feature_block_diagnostics_by_subject", class(out))
  out
}
