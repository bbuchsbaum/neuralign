#' Preprocessing Helpers
#'
#' Domain-agnostic preprocessing helpers for alignment inputs. These functions
#' are **opt-in**: \code{\link{fit_alignment}} does not perform any centering or
#' scaling unless you explicitly preprocess your matrices or feature blocks.
#'
#' @name preprocess
NULL

.validate_finite_matrix <- function(x, arg = "x") {
  if (!.is_matrixish(x)) {
    stop(sprintf("'%s' must be matrix-like", arg), call. = FALSE)
  }
  x <- as.matrix(x)
  if (length(x) == 0L) {
    stop(sprintf("'%s' must not be empty", arg), call. = FALSE)
  }
  if (any(!is.finite(x))) {
    stop(sprintf("'%s' contains non-finite values (NA/NaN/Inf)", arg), call. = FALSE)
  }
  x
}

.compute_center <- function(x, margin, robust) {
  if (isTRUE(robust)) {
    if (margin == 1L) {
      apply(x, 1L, stats::median)
    } else {
      apply(x, 2L, stats::median)
    }
  } else {
    if (margin == 1L) {
      rowMeans(x)
    } else {
      colMeans(x)
    }
  }
}

.compute_scale <- function(x, margin, robust) {
  if (isTRUE(robust)) {
    if (margin == 1L) {
      apply(x, 1L, stats::mad)
    } else {
      apply(x, 2L, stats::mad)
    }
  } else {
    if (margin == 1L) {
      apply(x, 1L, stats::sd)
    } else {
      apply(x, 2L, stats::sd)
    }
  }
}

.apply_sweep <- function(x, what, center_or_scale, fun) {
  if (identical(what, "rows")) {
    sweep(x, 1L, center_or_scale, fun, check.margin = FALSE)
  } else {
    sweep(x, 2L, center_or_scale, fun, check.margin = FALSE)
  }
}

#' Preprocess a Matrix for Alignment
#'
#' Center and/or scale a numeric matrix along rows or columns. This helper is
#' useful for removing mean offsets and/or equalizing feature scales prior to
#' fitting an alignment model.
#'
#' @param x Matrix-like numeric object.
#' @param center One of \code{"none"}, \code{"rows"}, or \code{"cols"}.
#' @param scale One of \code{"none"}, \code{"rows"}, or \code{"cols"}.
#' @param robust Logical; if TRUE, use median/MAD instead of mean/SD.
#'
#' @return A numeric matrix with the same dimensions and dimnames as \code{x}.
#'
#' @examples
#' X <- matrix(1:6, 3, 2)
#' preprocess_matrix(X, center = "rows", scale = "none")
#' @export
preprocess_matrix <- function(x,
                              center = c("none", "rows", "cols"),
                              scale = c("none", "rows", "cols"),
                              robust = FALSE) {
  x <- .validate_finite_matrix(x, arg = "x")
  center <- match.arg(center)
  scale <- match.arg(scale)

  if (!is.logical(robust) || length(robust) != 1L || is.na(robust)) {
    stop("'robust' must be TRUE or FALSE", call. = FALSE)
  }

  out <- x

  if (center != "none") {
    margin <- if (center == "rows") 1L else 2L
    center_vec <- .compute_center(out, margin = margin, robust = robust)
    out <- .apply_sweep(out, what = center, center_or_scale = center_vec, fun = "-")
  }

  if (scale != "none") {
    margin <- if (scale == "rows") 1L else 2L
    scale_vec <- .compute_scale(out, margin = margin, robust = robust)
    if (any(!is.finite(scale_vec))) {
      stop("Scaling produced non-finite scale factors", call. = FALSE)
    }
    zeros <- which(scale_vec == 0)
    if (length(zeros) > 0L) {
      warning(
        sprintf(
          "Replacing %d zero scale factor(s) with 1 (no scaling for those rows/cols)",
          length(zeros)
        ),
        call. = FALSE
      )
      scale_vec[zeros] <- 1
    }
    out <- .apply_sweep(out, what = scale, center_or_scale = scale_vec, fun = "/")
  }

  out
}

#' Preprocess AlignmentData
#'
#' Apply \code{\link{preprocess_matrix}} to each subject matrix in an
#' \code{\link{AlignmentData}} object (or a list coercible to \code{AlignmentData}).
#'
#' @param data \code{AlignmentData} or a named list of matrices.
#' @param ... Arguments passed to \code{\link{preprocess_matrix}}.
#'
#' @return An \code{AlignmentData} object with preprocessed subject matrices.
#'
#' @examples
#' set.seed(1)
#' adat <- AlignmentData(list(
#'   s1 = matrix(rnorm(20), 4, 5),
#'   s2 = matrix(rnorm(20), 4, 5)
#' ))
#' adat2 <- preprocess_alignment_data(adat, center = "rows")
#' dim(get_subject_data(adat2, "s1"))
#' @export
preprocess_alignment_data <- function(data, ...) {
  data <- as_alignment_data(data)
  data_list <- get_data_list(data)
  if (!all(vapply(data_list, .is_matrixish, logical(1)))) {
    bad <- names(data_list)[!vapply(data_list, .is_matrixish, logical(1))]
    stop(
      sprintf(
        "preprocess_alignment_data() only supports matrix/Matrix data; non-matrix subjects: %s",
        paste(bad, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  out_list <- lapply(data_list, preprocess_matrix, ...)

  AlignmentData(
    data = out_list,
    subjects = data@subjects,
    space = data@space,
    design = data@design,
    geometry = data@geometry,
    obs_labels = data@obs_labels,
    metadata = data@metadata
  )
}

#' Preprocess Feature Blocks
#'
#' Apply \code{\link{preprocess_matrix}} to the \code{$x} matrix of each
#' \code{\link{alignment_feature_block}} in a \code{blocks_by_subject} structure.
#'
#' @param blocks_by_subject Named list of subjects, each containing a named list
#'   of \code{"alignment_feature_block"} objects.
#' @param ... Arguments passed to \code{\link{preprocess_matrix}}.
#'
#' @return A \code{blocks_by_subject} list with the same structure and metadata,
#'   but with preprocessed \code{$x} matrices.
#'
#' @examples
#' set.seed(1)
#' b1 <- alignment_feature_block(matrix(rnorm(6), 3, 2), "b1",
#'   feature_names = c("a", "b", "c")
#' )
#' blocks <- list(s1 = list(b1 = b1))
#' blocks2 <- preprocess_feature_blocks(blocks, center = "rows")
#' dim(blocks2$s1$b1$x)
#' @export
preprocess_feature_blocks <- function(blocks_by_subject, ...) {
  if (!is.list(blocks_by_subject) || length(blocks_by_subject) < 1L) {
    stop("'blocks_by_subject' must be a non-empty named list", call. = FALSE)
  }
  if (is.null(names(blocks_by_subject)) || any(!nzchar(names(blocks_by_subject)))) {
    stop("'blocks_by_subject' must be named by subject", call. = FALSE)
  }

  lapply(blocks_by_subject, function(blocks) {
    if (!is.list(blocks) || length(blocks) < 1L) {
      stop("Each subject must provide a non-empty list of blocks", call. = FALSE)
    }
    if (is.null(names(blocks)) || any(!nzchar(names(blocks)))) {
      stop("Each subject's block list must be named by block name", call. = FALSE)
    }
    if (!all(vapply(blocks, inherits, logical(1), "alignment_feature_block"))) {
      stop("All blocks must be 'alignment_feature_block' objects", call. = FALSE)
    }

    lapply(blocks, function(block) {
      out <- block
      out$x <- preprocess_matrix(block$x, ...)
      out
    })
  })
}
