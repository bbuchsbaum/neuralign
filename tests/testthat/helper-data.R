make_test_subject_ids <- function(n_subjects, prefix = "sub-", width = 2L) {
  n_subjects <- as.integer(n_subjects)
  width <- as.integer(width)
  if (length(n_subjects) != 1L || is.na(n_subjects) || n_subjects < 1L) {
    stop("'n_subjects' must be a positive integer", call. = FALSE)
  }
  if (length(width) != 1L || is.na(width) || width < 1L) {
    stop("'width' must be a positive integer", call. = FALSE)
  }
  if (!is.character(prefix) || length(prefix) != 1L || is.na(prefix)) {
    stop("'prefix' must be a non-NA character scalar", call. = FALSE)
  }
  sprintf("%s%0*d", prefix, width, seq_len(n_subjects))
}

make_test_obs_labels <- function(n_obs, prefix = "obs_") {
  n_obs <- as.integer(n_obs)
  if (length(n_obs) != 1L || is.na(n_obs) || n_obs < 1L) {
    stop("'n_obs' must be a positive integer", call. = FALSE)
  }
  if (!is.character(prefix) || length(prefix) != 1L || is.na(prefix)) {
    stop("'prefix' must be a non-NA character scalar", call. = FALSE)
  }
  paste0(prefix, seq_len(n_obs))
}

make_test_matrix <- function(n_features = 10L, n_obs = 10L, fill = stats::rnorm) {
  n_features <- as.integer(n_features)
  n_obs <- as.integer(n_obs)
  if (length(n_features) != 1L || is.na(n_features) || n_features < 1L) {
    stop("'n_features' must be a positive integer", call. = FALSE)
  }
  if (length(n_obs) != 1L || is.na(n_obs) || n_obs < 1L) {
    stop("'n_obs' must be a positive integer", call. = FALSE)
  }
  if (!is.function(fill)) {
    stop("'fill' must be a function", call. = FALSE)
  }
  matrix(fill(n_features * n_obs), n_features, n_obs)
}

make_test_data_list <- function(n_subjects = 3L,
                                n_features = 10L,
                                n_obs = 10L,
                                subject_ids = NULL,
                                fill = stats::rnorm) {
  if (is.null(subject_ids)) {
    subject_ids <- make_test_subject_ids(n_subjects)
  }
  if (!is.character(subject_ids) || anyNA(subject_ids) || any(!nzchar(subject_ids))) {
    stop("'subject_ids' must be a character vector of non-empty, non-NA values", call. = FALSE)
  }
  mats <- lapply(seq_along(subject_ids), function(i) make_test_matrix(n_features, n_obs, fill = fill))
  names(mats) <- subject_ids
  mats
}

make_test_alignment_data <- function(n_subjects = 3L,
                                     n_features = 10L,
                                     n_obs = 10L,
                                     subject_ids = NULL,
                                     obs_labels = NULL,
                                     fill = stats::rnorm,
                                     ...) {
  data_list <- make_test_data_list(
    n_subjects = n_subjects,
    n_features = n_features,
    n_obs = n_obs,
    subject_ids = subject_ids,
    fill = fill
  )
  AlignmentData(data_list, obs_labels = obs_labels, ...)
}

make_random_orthogonal <- function(d, allow_reflection = FALSE) {
  d <- as.integer(d)
  if (length(d) != 1L || is.na(d) || d < 1L) {
    stop("'d' must be a positive integer", call. = FALSE)
  }
  Q <- qr.Q(qr(matrix(stats::rnorm(d * d), d, d)))
  if (!isTRUE(allow_reflection) && det(Q) < 0) {
    Q[, d] <- -Q[, d]
  }
  Q
}

