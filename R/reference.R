#' Select Reference for Alignment
#'
#' Select a reference subject or compute a reference template for alignment.
#' The reference selection is done using only the provided subjects, making
#' it safe for use within cross-validation folds.
#'
#' @param data AlignmentData object or named list of matrices.
#' @param method Method for reference selection:
#'   \itemize{
#'     \item "medoid" - Subject most similar to all others (default)
#'     \item "centroid" - Subject closest to the mean
#'     \item "first" - First subject in the data
#'     \item "random" - Random subject
#'   }
#' @param distance Distance metric for medoid/centroid:
#'   \itemize{
#'     \item "correlation" - 1 - mean pairwise correlation (default)
#'     \item "euclidean" - Euclidean distance between flattened data
#'     \item "frobenius" - Frobenius norm of difference
#'   }
#' @param seed Random seed for reproducibility (used if method="random").
#'
#' @return Character string with the selected subject ID.
#'
#' @examples
#' \dontrun{
#' data_list <- list(
#'   "sub-01" = matrix(rnorm(100*50), 100, 50),
#'   "sub-02" = matrix(rnorm(100*50), 100, 50),
#'   "sub-03" = matrix(rnorm(100*50), 100, 50)
#' )
#' adat <- AlignmentData(data_list)
#'
#' # Select medoid
#' ref <- select_reference(adat, method = "medoid")
#'
#' # Select centroid
#' ref <- select_reference(adat, method = "centroid")
#' }
#'
#' @export
select_reference <- function(data,
                             method = c("medoid", "centroid", "first", "random"),
                             distance = c("correlation", "euclidean", "frobenius", "procrustes"),
                             seed = NULL) {
  method <- match.arg(method)
  distance <- match.arg(distance)

  # Coerce to AlignmentData if needed
  if (!inherits(data, "AlignmentData")) {
    data <- as_alignment_data(data)
  }

  subjects <- data@subjects
  n <- length(subjects)

  if (n == 0) {
    stop("No subjects in data")
  }

  if (n == 1) {
    return(subjects[1])
  }

  # Simple cases
  if (method == "first") {
    return(subjects[1])
  }

  if (method == "random") {
    if (!is.null(seed)) {
      set.seed(seed)
    }
    return(sample(subjects, 1))
  }

  # Compute distance matrix
  dist_mat <- .compute_distance_matrix(data, distance)

  if (method == "medoid") {
    # Medoid: subject with minimum average distance to all others.
    # When some pairs have no overlap (e.g., disjoint label sets), distances
    # are NA; we ignore those pairs rather than failing outright.
    total_dist <- rowMeans(dist_mat, na.rm = TRUE)
    total_dist[!is.finite(total_dist)] <- Inf
    if (all(is.infinite(total_dist))) {
      stop(
        "Cannot select medoid: no finite pairwise distances (check obs_labels overlap across subjects)",
        call. = FALSE
      )
    }
    idx <- which.min(total_dist)
    return(subjects[idx])
  }

  if (method == "centroid") {
    # Centroid requires a well-defined mean across subjects; if observation
    # labels differ across subjects, this is not meaningful.
    labels_by_subject <- .resolve_obs_labels_by_subject(data)
    if (!is.null(labels_by_subject)) {
      labs <- unique(lapply(labels_by_subject, paste, collapse = "\r"))
      if (length(labs) > 1) {
        stop(
          paste0(
            "Cannot use method='centroid' when obs_labels differ across subjects. ",
            "Use method='medoid' (distance='procrustes' recommended) or harmonize observations first."
          ),
          call. = FALSE
        )
      }
    }

    # Centroid: subject closest to the mean
    # First compute mean data
    data_list <- get_data_list(data)
    mean_data <- Reduce(`+`, data_list) / n

    # Distance from each subject to mean
    dist_to_mean <- sapply(data_list, function(x) {
      .compute_pairwise_distance(x, mean_data, distance)
    })

    idx <- which.min(dist_to_mean)
    return(subjects[idx])
  }

  stop(sprintf("Unknown method: %s", method))
}


#' Compute Distance Matrix Between Subjects
#'
#' @param data AlignmentData object.
#' @param distance Distance metric.
#'
#' @return Square distance matrix.
#'
#' @keywords internal
.compute_distance_matrix <- function(data, distance) {
  subjects <- data@subjects
  n <- length(subjects)
  data_list <- get_data_list(data)

  dist_mat <- matrix(NA_real_, n, n)
  rownames(dist_mat) <- colnames(dist_mat) <- subjects
  diag(dist_mat) <- 0

  labels_by_subject <- .resolve_obs_labels_by_subject(data)

  for (i in seq_len(n - 1)) {
    for (j in (i + 1):n) {
      labs_i <- if (!is.null(labels_by_subject)) labels_by_subject[[i]] else NULL
      labs_j <- if (!is.null(labels_by_subject)) labels_by_subject[[j]] else NULL
      d <- .compute_pairwise_distance(data_list[[i]], data_list[[j]], distance,
        obs_labels_x = labs_i,
        obs_labels_y = labs_j
      )
      dist_mat[i, j] <- dist_mat[j, i] <- d
    }
  }

  dist_mat
}

.match_obs_indices <- function(labels_source, labels_target, min_overlap = 2L) {
  if (is.null(labels_source) || is.null(labels_target)) {
    return(list(source = NULL, target = NULL, labels = NULL))
  }
  labels_source <- as.character(labels_source)
  labels_target <- as.character(labels_target)
  common <- intersect(labels_source, labels_target)
  if (length(common) < min_overlap) {
    stop(
      sprintf(
        "Not enough shared observation labels: %d (need >= %d)",
        length(common), as.integer(min_overlap)
      ),
      call. = FALSE
    )
  }
  list(
    source = match(common, labels_source),
    target = match(common, labels_target),
    labels = common
  )
}

.resolve_obs_labels_by_subject <- function(data) {
  labels <- data@obs_labels
  if (is.null(labels)) return(NULL)

  subjects <- data@subjects
  n <- length(subjects)

  if (is.atomic(labels)) {
    labs <- as.character(labels)
    out <- rep(list(labs), n)
    names(out) <- subjects
    return(out)
  }

  if (is.list(labels)) {
    if (is.null(names(labels))) {
      if (length(labels) != n) {
        stop(
          "AlignmentData@obs_labels is a list without names; its length must match number of subjects",
          call. = FALSE
        )
      }
      names(labels) <- subjects
    }
    missing <- setdiff(subjects, names(labels))
    if (length(missing) > 0) {
      stop(
        sprintf(
          "AlignmentData@obs_labels list is missing subjects: %s",
          paste(missing, collapse = ", ")
        ),
        call. = FALSE
      )
    }
    out <- lapply(subjects, function(s) as.character(labels[[s]]))
    names(out) <- subjects
    return(out)
  }

  stop("AlignmentData@obs_labels must be NULL, an atomic vector, or a (named) list", call. = FALSE)
}

.subset_to_overlap <- function(x, y, obs_labels_x = NULL, obs_labels_y = NULL, min_overlap = 2L) {
  x <- as.matrix(x)
  y <- as.matrix(y)

  if (is.null(obs_labels_x) && is.null(obs_labels_y)) {
    if (!identical(dim(x), dim(y))) {
      stop(
        sprintf(
          "Matrices have incompatible dimensions: x is %d x %d, y is %d x %d",
          nrow(x), ncol(x), nrow(y), ncol(y)
        ),
        call. = FALSE
      )
    }
    return(list(x = x, y = y, labels = NULL))
  }

  idx <- .match_obs_indices(obs_labels_x, obs_labels_y, min_overlap = min_overlap)
  x <- x[, idx$source, drop = FALSE]
  y <- y[, idx$target, drop = FALSE]
  list(x = x, y = y, labels = idx$labels)
}

#' Compute Distance Between Two Matrices
#'
#' @param x First matrix.
#' @param y Second matrix.
#' @param distance Distance metric.
#' @param obs_labels_x Optional observation labels for x.
#' @param obs_labels_y Optional observation labels for y.
#' @param min_overlap Minimum overlap required when labels are supplied.
#'
#' @return Scalar distance.
#'
#' @keywords internal
.compute_pairwise_distance <- function(x,
                                       y,
                                       distance,
                                       obs_labels_x = NULL,
                                       obs_labels_y = NULL,
                                       min_overlap = 2L) {
  min_overlap <- as.integer(min_overlap)
  if (!is.finite(min_overlap) || min_overlap < 1L) {
    stop("'min_overlap' must be a positive integer", call. = FALSE)
  }

  # Subset to matched observations when labels are supplied
  out <- tryCatch(
    .subset_to_overlap(x, y, obs_labels_x = obs_labels_x, obs_labels_y = obs_labels_y, min_overlap = min_overlap),
    error = function(e) NULL
  )
  if (is.null(out)) {
    return(NA_real_)
  }
  x <- out$x
  y <- out$y

  if (distance == "correlation") {
    # Mean correlation across features/rows
    n_features <- nrow(x)
    cors <- sapply(seq_len(n_features), function(i) {
      cor(x[i, ], y[i, ])
    })
    # Distance = 1 - mean correlation
    return(1 - mean(cors, na.rm = TRUE))
  }

  if (distance == "euclidean") {
    # Euclidean distance of flattened matrices
    return(sqrt(sum((x - y)^2)))
  }

  if (distance == "frobenius") {
    # Frobenius norm
    return(norm(x - y, "F"))
  }

  if (distance == "procrustes") {
    return(procrustes_distance(x, y, convention = "left"))
  }

  stop(sprintf("Unknown distance: %s", distance))
}


#' Compute Mean/Centroid Data
#'
#' Compute the element-wise mean across subjects.
#'
#' @param data AlignmentData object.
#'
#' @return Matrix representing the mean across subjects.
#'
#' @export
compute_centroid <- function(data) {
  if (!inherits(data, "AlignmentData")) {
    data <- as_alignment_data(data)
  }

  data_list <- get_data_list(data)
  n <- length(data_list)

  if (n == 0) {
    stop("No subjects in data")
  }

  Reduce(`+`, data_list) / n
}


#' Get Reference Data
#'
#' Get the data for a reference subject, or compute reference if "consensus".
#'
#' @param data AlignmentData object.
#' @param reference Subject ID or "consensus".
#'
#' @return Matrix representing the reference data.
#'
#' @export
get_reference_data <- function(data, reference) {
  if (!inherits(data, "AlignmentData")) {
    data <- as_alignment_data(data)
  }

  if (.is_matrixish(reference)) {
    return(as.matrix(reference))
  }

  if (reference == "consensus") {
    return(compute_centroid(data))
  }

  if (!reference %in% data@subjects) {
    stop(sprintf("Subject '%s' not found in data", reference))
  }

  get_subject_data(data, reference)
}


#' Validate Reference for Aligner
#'
#' Check if a reference specification is valid for a given aligner.
#'
#' @param reference Reference specification.
#' @param method Aligner method name.
#'
#' @return TRUE if valid, throws error otherwise.
#'
#' @keywords internal
.validate_reference_for_method <- function(reference, method) {
  caps <- aligner_capabilities(method)

  if (is.null(caps)) {
    return(TRUE)  # Unknown method, can't validate
  }

  ref_types <- caps$reference_types %||% c("subject", "consensus")

  ref_type <- if (.is_matrixish(reference)) {
    "template"
  } else if (reference == "consensus") {
    "consensus"
  } else if (reference == "medoid" || reference == "centroid") {
    "subject"
  } else {
    "subject"  # Assume it's a subject ID
  }

  if (!ref_type %in% ref_types) {
    stop(sprintf(
      "Method '%s' does not support reference type '%s'. Supported: %s",
      method, ref_type, paste(ref_types, collapse = ", ")
    ))
  }

  TRUE
}
