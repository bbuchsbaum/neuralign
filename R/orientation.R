#' Orientation Boundary Between Algorithm and Analysis Layers
#'
#' `neuralign` algorithm objects ([AlignmentData], [AlignmentModel] transforms)
#' use **features × observations** with left-multiply `Q %*% X`.
#'
#' Analysis-facing objects ([AlignedStudy]) use **observations × shared features**
#' so observation metadata has one row per observation. The conversion boundary
#' performs exactly one explicit transpose.
#'
#' @name orientation
#' @family aligned_study
NULL


#' Transpose Algorithm Matrix to Analysis Orientation
#'
#' @param x Matrix in features × observations orientation.
#' @return Matrix in observations × features orientation.
#'
#' @keywords internal
to_analysis_matrix <- function(x) {
  if (!.is_matrixish(x)) {
    stop("'x' must be a matrix or Matrix", call. = FALSE)
  }
  t(x)
}


#' Transpose Analysis Matrix to Algorithm Orientation
#'
#' @param x Matrix in observations × features orientation.
#' @return Matrix in features × observations orientation.
#'
#' @keywords internal
to_algorithm_matrix <- function(x) {
  if (!.is_matrixish(x)) {
    stop("'x' must be a matrix or Matrix", call. = FALSE)
  }
  t(x)
}
