.resolve_subject_subset <- function(i, subjects, what = "subjects") {
  if (missing(i)) return(subjects)

  if (is.character(i)) {
    if (anyNA(i) || any(!nzchar(i))) {
      stop(sprintf("Invalid %s: subject ids must be non-empty strings without NA", what), call. = FALSE)
    }
    missing_subj <- setdiff(i, subjects)
    if (length(missing_subj) > 0) {
      stop(
        sprintf("Unknown %s: %s", what, paste(missing_subj, collapse = ", ")),
        call. = FALSE
      )
    }
    return(i)
  }

  if (is.logical(i)) {
    if (anyNA(i)) stop(sprintf("Invalid %s: logical index contains NA", what), call. = FALSE)
    if (length(i) != length(subjects)) {
      stop(sprintf("Invalid %s: logical index must have length %d", what, length(subjects)), call. = FALSE)
    }
    return(subjects[which(i)])
  }

  if (is.numeric(i)) {
    if (anyNA(i) || any(!is.finite(i))) {
      stop(sprintf("Invalid %s: numeric index must be finite and non-NA", what), call. = FALSE)
    }
    i_int <- as.integer(i)
    if (any(i_int != i)) {
      stop(sprintf("Invalid %s: numeric index must be integer-valued", what), call. = FALSE)
    }

    # Drop zeros (base R convention)
    i_int <- i_int[i_int != 0L]
    if (length(i_int) == 0L) return(subjects[integer(0)])

    if (all(i_int < 0L)) {
      drop <- abs(i_int)
      if (any(drop > length(subjects))) {
        stop(
          sprintf("Invalid %s: negative indices out of bounds (max %d)", what, length(subjects)),
          call. = FALSE
        )
      }
      keep <- setdiff(seq_along(subjects), drop)
      return(subjects[keep])
    }

    if (any(i_int < 0L)) {
      stop(sprintf("Invalid %s: cannot mix negative and positive indices", what), call. = FALSE)
    }

    if (any(i_int < 1L | i_int > length(subjects))) {
      stop(
        sprintf("Invalid %s: indices out of bounds (valid range: 1..%d)", what, length(subjects)),
        call. = FALSE
      )
    }
    return(subjects[i_int])
  }

  stop(sprintf("Invalid %s index type: %s", what, class(i)[1]), call. = FALSE)
}

