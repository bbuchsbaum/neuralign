#' Infix Helpers
#'
#' Internal infix operators used throughout the package.
#'
#' @keywords internal
NULL

#' Null-coalescing operator
#'
#' Returns `y` when `x` is `NULL`, otherwise `x`.
#'
#' @keywords internal
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

