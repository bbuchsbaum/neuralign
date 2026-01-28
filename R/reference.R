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
#'     \item "procrustes" - Procrustes residual after optimal orthogonal alignment
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
    stop("No subjects in data", call. = FALSE)
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
      old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
      } else {
        NULL
      }
      on.exit({
        if (is.null(old_seed)) {
          if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
            rm(".Random.seed", envir = .GlobalEnv)
          }
        } else {
          assign(".Random.seed", old_seed, envir = .GlobalEnv)
        }
      }, add = TRUE)
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
    dm <- dist_mat
    diag(dm) <- NA_real_
    total_dist <- rowMeans(dm, na.rm = TRUE)
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
    mean_data <- compute_centroid(data_list)

    # Distance from each subject to mean
    dist_to_mean <- vapply(
      data_list,
      function(x) .compute_pairwise_distance(x, mean_data, distance),
      numeric(1)
    )

    idx <- which.min(dist_to_mean)
    return(subjects[idx])
  }

  stop(sprintf("Unknown method: %s", method), call. = FALSE)
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

  ok_source <- !is.na(labels_source)
  ok_target <- !is.na(labels_target)
  if (!any(ok_source) || !any(ok_target)) {
    stop(
      sprintf(
        "Not enough shared observation labels: %d (need >= %d)",
        0L, as.integer(min_overlap)
      ),
      call. = FALSE
    )
  }

  src_idx_map <- which(ok_source)
  tgt_idx_map <- which(ok_target)
  src_vals <- labels_source[ok_source]
  tgt_vals <- labels_target[ok_target]

  src_pos_by_label <- split(seq_along(src_vals), src_vals)
  tgt_pos_by_label <- split(seq_along(tgt_vals), tgt_vals)
  common_labels <- intersect(names(src_pos_by_label), names(tgt_pos_by_label))

  if (length(common_labels) < 1L) {
    stop(
      sprintf(
        "Not enough shared observation labels: %d (need >= %d)",
        0L, as.integer(min_overlap)
      ),
      call. = FALSE
    )
  }

  source <- integer(0)
  target <- integer(0)
  labels <- character(0)
  for (lab in common_labels) {
    src_pos <- src_pos_by_label[[lab]]
    tgt_pos <- tgt_pos_by_label[[lab]]
    n_match <- min(length(src_pos), length(tgt_pos))
    if (n_match < 1L) next

    source <- c(source, src_idx_map[src_pos[seq_len(n_match)]])
    target <- c(target, tgt_idx_map[tgt_pos[seq_len(n_match)]])
    labels <- c(labels, rep(lab, n_match))
  }

  if (length(source) < min_overlap) {
    stop(
      sprintf(
        "Not enough shared observation labels: %d (need >= %d)",
        length(source), as.integer(min_overlap)
      ),
      call. = FALSE
    )
  }

  ord <- order(source)
  list(
    source = source[ord],
    target = target[ord],
    labels = labels[ord]
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


.require_obs_labels_by_subject <- function(data, method = NULL, check_length = TRUE) {
  labels_by_subject <- .resolve_obs_labels_by_subject(data)
  if (is.null(labels_by_subject)) {
    label_method <- method %||% "This method"
    stop(
      sprintf(
        "%s requires AlignmentData@obs_labels (atomic shared labels or a per-subject named list)",
        label_method
      ),
      call. = FALSE
    )
  }

  if (isTRUE(check_length)) {
    data_list <- get_data_list(data)
    for (subj in data@subjects) {
      x <- data_list[[subj]]
      if (!.is_matrixish(x)) x <- as.matrix(x)
      if (length(labels_by_subject[[subj]]) != ncol(x)) {
        stop(
          sprintf("obs_labels length mismatch for subject '%s'", subj),
          call. = FALSE
        )
      }
      if (anyNA(labels_by_subject[[subj]])) {
        stop(
          sprintf("obs_labels for subject '%s' contains NA values", subj),
          call. = FALSE
        )
      }
    }
  }

  labels_by_subject
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
    cors <- vapply(seq_len(n_features), function(i) cor(x[i, ], y[i, ]), numeric(1))
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

  stop(sprintf("Unknown distance: %s", distance), call. = FALSE)
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
  if (inherits(data, "AlignmentData")) {
    data_list <- get_data_list(data)
  } else if (is.list(data)) {
    data_list <- data
  } else {
    data <- as_alignment_data(data)
    data_list <- get_data_list(data)
  }
  n <- length(data_list)

  if (n == 0) {
    stop("No subjects in data", call. = FALSE)
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
    stop(sprintf("Subject '%s' not found in data", reference), call. = FALSE)
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
    stop(
      sprintf(
        "Method '%s' does not support reference type '%s'. Supported: %s",
        method, ref_type, paste(ref_types, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  TRUE
}
