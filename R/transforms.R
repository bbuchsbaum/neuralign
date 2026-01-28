#' Transform Objects
#'
#' Internal helpers for representing and applying transforms that are not best
#' represented as a dense matrix (e.g., low-rank operators) and for supporting
#' embedding-returning aligners.
#'
#' @name transforms
#' @keywords internal
NULL

.as_dense_matrix <- function(x) {
  if (inherits(x, "Matrix")) return(as.matrix(x))
  if (is.matrix(x)) return(x)
  as.matrix(x)
}

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

.compose_operator_transforms <- function(t2, t1, context = "compose transforms") {
  if (!.transform_is_operator(t1) || !.transform_is_operator(t2)) {
    stop(context, ": expected operator transforms", call. = FALSE)
  }

  dims1 <- .transform_dims(t1)
  dims2 <- .transform_dims(t2)
  if (!identical(dims2[["source"]], dims1[["target"]])) {
    stop(
      sprintf(
        "%s: dimension mismatch: right transform target=%d does not match left transform source=%d",
        context, dims1[["target"]], dims2[["source"]]
      ),
      call. = FALSE
    )
  }

  if (.is_low_rank_transform(t2) || .is_low_rank_transform(t1)) {
    # Products involving low-rank transforms remain low-rank. Represent the
    # composition without materializing dense (target x source) operators.
    if (.is_low_rank_transform(t2) && .is_low_rank_transform(t1)) {
      # (U2 V2^T)(U1 V1^T) = U2 (V2^T U1) V1^T
      # Use an SVD of the small cross term to stabilize and cap the rank at
      # min(rank1, rank2).
      C <- crossprod(t2$V, t1$U) # (r2 x r1)
      if (!length(C)) {
        return(.new_low_rank_transform(matrix(0, nrow(t2$U), 0L), matrix(0, nrow(t1$V), 0L)))
      }
      sv <- svd(C)
      U_new <- t2$U %*% sv$u
      U_new <- sweep(U_new, 2L, sv$d, "*")
      V_new <- t1$V %*% sv$v
      return(.new_low_rank_transform(U_new, V_new))
    }

    if (.is_low_rank_transform(t2) && .is_matrixish(t1)) {
      # (U2 V2^T) M = U2 (V2^T M) = U2 ( (t(M) V2)^T )
      V_new <- t(t1) %*% t2$V
      return(.new_low_rank_transform(t2$U, V_new))
    }

    if (.is_matrixish(t2) && .is_low_rank_transform(t1)) {
      # M (U1 V1^T) = (M U1) V1^T
      U_new <- t2 %*% t1$U
      return(.new_low_rank_transform(U_new, t1$V))
    }
  }

  t2 %*% t1
}

.low_rank_operator_svd <- function(transform, context = "low-rank operator") {
  if (!.is_low_rank_transform(transform)) {
    stop(context, ": expected a low-rank operator transform", call. = FALSE)
  }

  U0 <- .as_dense_matrix(transform$U)
  V0 <- .as_dense_matrix(transform$V)

  if (!identical(ncol(U0), ncol(V0))) {
    stop(context, ": low-rank factors must have matching column counts", call. = FALSE)
  }

  r <- ncol(U0)
  if (r == 0L) {
    return(list(
      U = matrix(0, nrow(U0), 0L),
      V = matrix(0, nrow(V0), 0L),
      d = numeric(0)
    ))
  }

  qru <- qr(U0)
  Qu <- qr.Q(qru)
  Ru <- qr.R(qru)

  qrv <- qr(V0)
  Qv <- qr.Q(qrv)
  Rv <- qr.R(qrv)

  B <- Ru %*% t(Rv)
  sv <- svd(B)

  list(
    U = Qu %*% sv$u,
    V = Qv %*% sv$v,
    d = sv$d
  )
}
