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
    if (!.is_matrixish(tmat)) {
      stop(
        context, ": transform for '", subj, "' is not a matrix/Matrix operator",
        call. = FALSE
      )
    }

    x <- data_list[[subj]]
    if (!.is_matrixish(x)) x <- as.matrix(x)

    if (ncol(tmat) != nrow(x)) {
      stop(sprintf(
        "%s: transform/data dimension mismatch for '%s': transform is %d x %d (target x source), data is %d x %d (features x obs)",
        context, subj, nrow(tmat), ncol(tmat), nrow(x), ncol(x)
      ), call. = FALSE)
    }
  }

  invisible(TRUE)
}
