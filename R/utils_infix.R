#' Infix Helpers
#'
#' Internal infix operators used throughout the package.
#'
#' @name utils_infix
#' @keywords internal
NULL

#' Null-coalescing operator
#'
#' Returns `y` when `x` is `NULL`, otherwise `x`.
#'
#' @rdname utils_infix
#' @keywords internal
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
