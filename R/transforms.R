#' Transform Objects
#'
#' Internal helpers for representing and applying transforms that are not best
#' represented as a dense matrix (e.g., low-rank operators) and for supporting
#' embedding-returning aligners.
#'
#' @name transforms
#' @keywords internal
NULL

.is_low_rank_transform <- function(x) {
  inherits(x, "neuralign_low_rank_transform")
}

.is_embedding_transform <- function(x) {
  inherits(x, "neuralign_embedding_transform")
}

.transform_is_operator <- function(x) {
  .is_matrixish(x) || .is_low_rank_transform(x)
}

.transform_dims <- function(x) {
  if (.is_matrixish(x)) {
    return(c(target = nrow(x), source = ncol(x)))
  }
  if (.is_low_rank_transform(x)) {
    return(c(target = nrow(x$U), source = nrow(x$V)))
  }
  if (.is_embedding_transform(x)) {
    return(c(target = nrow(x$aligned), observations = ncol(x$aligned)))
  }
  stop("Unsupported transform type", call. = FALSE)
}

.new_low_rank_transform <- function(U, V) {
  if (!.is_matrixish(U) || !.is_matrixish(V)) {
    stop("low-rank transform requires matrix-like U and V", call. = FALSE)
  }
  U <- .as_dense_matrix(U)
  V <- .as_dense_matrix(V)
  if (ncol(U) != ncol(V)) {
    stop(
      sprintf(
        "low-rank transform requires matching factor ranks: ncol(U)=%d, ncol(V)=%d",
        ncol(U), ncol(V)
      ),
      call. = FALSE
    )
  }
  structure(
    list(U = U, V = V),
    class = c("neuralign_low_rank_transform", "neuralign_transform")
  )
}

.new_embedding_transform <- function(aligned, subject = NULL, obs_labels = NULL) {
  if (!.is_matrixish(aligned)) {
    stop("embedding transform requires 'aligned' to be a matrix/Matrix", call. = FALSE)
  }
  aligned <- .as_dense_matrix(aligned)
  structure(
    list(
      aligned = aligned,
      subject = subject,
      obs_labels = obs_labels
    ),
    class = c("neuralign_embedding_transform", "neuralign_transform")
  )
}

.as_matrix_transform <- function(x) {
  if (.is_matrixish(x)) return(.as_dense_matrix(x))
  if (.is_low_rank_transform(x)) return(x$U %*% t(x$V))
  stop("Transform cannot be materialized as a matrix", call. = FALSE)
}
