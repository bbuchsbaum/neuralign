#' Feature Block Helpers
#'
#' Lightweight helpers for constructing, harmonizing, and stacking alignment
#' feature matrices in a way that can be shared across downstream packages.
#'
#' @name feature_blocks
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
  feature_names
}

#' Create an Alignment Feature Block
#'
#' @param x Matrix of features for a block (features x observations).
#' @param name Block name (character scalar).
#' @param weight Non-negative scalar weight (applied as `sqrt(weight)` in stacking).
#' @param feature_names Optional character vector naming the rows of `x`.
#'
#' @return An object of class `"alignment_feature_block"`.
#' @export
alignment_feature_block <- function(x, name, weight = 1, feature_names = NULL) {
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

  feature_names <- .validate_feature_names(feature_names, nrow(x), name)

  structure(
    list(
      x = x,
      name = name,
      weight = as.numeric(weight),
      feature_names = feature_names
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

#' Stack Feature Blocks
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
