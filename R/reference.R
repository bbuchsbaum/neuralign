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
                             distance = c("correlation", "euclidean", "frobenius"),
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
    # Medoid: subject with minimum sum of distances to all others
    total_dist <- rowSums(dist_mat)
    idx <- which.min(total_dist)
    return(subjects[idx])
  }

  if (method == "centroid") {
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

  dist_mat <- matrix(0, n, n)
  rownames(dist_mat) <- colnames(dist_mat) <- subjects

  for (i in seq_len(n - 1)) {
    for (j in (i + 1):n) {
      d <- .compute_pairwise_distance(data_list[[i]], data_list[[j]], distance)
      dist_mat[i, j] <- dist_mat[j, i] <- d
    }
  }

  dist_mat
}


#' Compute Distance Between Two Matrices
#'
#' @param x First matrix.
#' @param y Second matrix.
#' @param distance Distance metric.
#'
#' @return Scalar distance.
#'
#' @keywords internal
.compute_pairwise_distance <- function(x, y, distance) {
  # Ensure matrices
  x <- as.matrix(x)
  y <- as.matrix(y)

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

  if (is.matrix(reference)) {
    return(reference)
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

  ref_type <- if (is.matrix(reference)) {
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
