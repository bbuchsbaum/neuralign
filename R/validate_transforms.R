#' Transform Validation Helpers
#'
#' Internal helpers for validating operator transforms returned by aligners.
#'
#' @name validate_transforms
#' @keywords internal
NULL

.is_matrixish <- function(x) {
  is.matrix(x) || inherits(x, "Matrix")
}

.validate_operator_transforms <- function(transforms, data_list, context) {
  if (!is.list(transforms) || !length(transforms) || is.null(names(transforms))) {
    stop(context, ": expected 'transforms' to be a named list", call. = FALSE)
  }
  missing_subj <- setdiff(names(data_list), names(transforms))
  if (length(missing_subj)) {
    stop(
      context, ": missing transforms for subject(s): ", paste(missing_subj, collapse = ", "),
      call. = FALSE
    )
  }

  for (subj in names(data_list)) {
    tmat <- transforms[[subj]]
    if (!.transform_is_operator(tmat)) {
      stop(
        context, ": transform for '", subj, "' is not an operator transform",
        call. = FALSE
      )
    }

    x <- data_list[[subj]]
    if (!.is_matrixish(x)) x <- as.matrix(x)

    dims <- .transform_dims(tmat)
    if (dims[["source"]] != nrow(x)) {
      stop(sprintf(
        "%s: transform/data dimension mismatch for '%s': transform is %d x %d (target x source), data is %d x %d (features x obs)",
        context, subj, dims[["target"]], dims[["source"]], nrow(x), ncol(x)
      ), call. = FALSE)
    }
  }

  invisible(TRUE)
}

.validate_embedding_aligned <- function(aligned, data_list, context) {
  if (!is.list(aligned) || is.null(names(aligned))) {
    stop(context, ": expected 'aligned' to be a named list", call. = FALSE)
  }
  missing_subj <- setdiff(names(data_list), names(aligned))
  if (length(missing_subj)) {
    stop(
      context, ": missing aligned embeddings for subject(s): ", paste(missing_subj, collapse = ", "),
      call. = FALSE
    )
  }

  emb_rows <- NULL
  for (subj in names(data_list)) {
    z <- aligned[[subj]]
    if (!.is_matrixish(z)) {
      stop(context, ": aligned embedding for '", subj, "' is not matrix-like", call. = FALSE)
    }
    z <- as.matrix(z)
    x <- data_list[[subj]]
    if (!.is_matrixish(x)) x <- as.matrix(x)

    if (ncol(z) != ncol(x)) {
      stop(sprintf(
        "%s: aligned/data observation mismatch for '%s': aligned has %d obs, data has %d obs",
        context, subj, ncol(z), ncol(x)
      ), call. = FALSE)
    }

    if (is.null(emb_rows)) {
      emb_rows <- nrow(z)
    } else if (!identical(nrow(z), emb_rows)) {
      stop(sprintf(
        "%s: embedding dimension mismatch for '%s': expected %d rows, got %d",
        context, subj, emb_rows, nrow(z)
      ), call. = FALSE)
    }
  }

  invisible(TRUE)
}

.validate_embedding_transforms <- function(transforms, data_list, context) {
  if (!is.list(transforms) || is.null(names(transforms))) {
    stop(context, ": expected 'transforms' to be a named list", call. = FALSE)
  }
  missing_subj <- setdiff(names(data_list), names(transforms))
  if (length(missing_subj)) {
    stop(
      context, ": missing transforms for subject(s): ", paste(missing_subj, collapse = ", "),
      call. = FALSE
    )
  }

  for (subj in names(data_list)) {
    tmat <- transforms[[subj]]
    x <- data_list[[subj]]
    if (!.is_matrixish(x)) x <- as.matrix(x)

    if (.transform_is_operator(tmat)) {
      dims <- .transform_dims(tmat)
      if (dims[["source"]] != nrow(x)) {
        stop(sprintf(
          "%s: transform/data dimension mismatch for '%s': transform source=%d, data rows=%d",
          context, subj, dims[["source"]], nrow(x)
        ), call. = FALSE)
      }
      next
    }

    if (.is_embedding_transform(tmat)) {
      z <- tmat$aligned
      if (ncol(z) != ncol(x)) {
        stop(sprintf(
          "%s: embedding transform mismatch for '%s': embedding has %d obs, data has %d obs",
          context, subj, ncol(z), ncol(x)
        ), call. = FALSE)
      }
      next
    }

    stop(context, ": transform for '", subj, "' is not a supported operator/embedding transform", call. = FALSE)
  }

  invisible(TRUE)
}
