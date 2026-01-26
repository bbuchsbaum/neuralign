#' @import methods
#' @importFrom stats cor sd var setNames
#' @importFrom utils combn packageVersion write.csv read.csv
#' @importFrom Matrix Matrix
#' @importFrom digest digest
NULL

#' AlignmentData Class
#'
#' Input container for data to be aligned. Standardizes the representation
#' of multi-subject neuroimaging data for alignment methods.
#'
#' @slot data Named list of matrices or NeuroVec objects (one per subject).
#'   Each matrix should be (features x observations) or (voxels x time).
#' @slot subjects Character vector of subject IDs, matching names in data.
#' @slot space Optional space specification (gds_space or compatible).
#' @slot design Optional design/task structure for methods that need it (e.g., dkge).
#' @slot geometry Optional adjacency/graph structure for graph-based methods.
#' @slot obs_labels Optional vector of length n_obs defining a shared observation
#'   label system across subjects (e.g., design column names, trial ids, or
#'   timepoints for time-locked paradigms). When provided, validation can enforce
#'   that observations are comparable across subjects.
#' @slot metadata Named list for additional method-specific information.
#'
#' @export
setClass("AlignmentData",
  slots = c(
    data = "list",
    subjects = "character",
    space = "ANY",
    design = "ANY",
    geometry = "ANY",
    obs_labels = "ANY",
    metadata = "list"
  ),
  prototype = list(
    data = list(),
    subjects = character(0),
    space = NULL,
    design = NULL,
    geometry = NULL,
    obs_labels = NULL,
    metadata = list()
  )
)


#' Create an AlignmentData Object
#'
#' @param data A named list of matrices or NeuroVec objects, one per subject.
#'   Names should be subject identifiers.
#' @param subjects Optional character vector of subject IDs. If NULL, extracted
#'   from names(data).
#' @param space Optional space specification (e.g., gds_space object).
#' @param design Optional design/task structure for supervised methods.
#' @param geometry Optional adjacency matrix or graph for topology-aware methods.
#' @param obs_labels Optional shared observation labels (see slot description).
#' @param metadata Optional named list of additional metadata.
#'
#' @return An AlignmentData object.
#'
#' @examples
#' # Create from list of matrices
#' data_list <- list(
#'   "sub-01" = matrix(rnorm(100*50), 100, 50),
#'   "sub-02" = matrix(rnorm(100*50), 100, 50),
#'   "sub-03" = matrix(rnorm(100*50), 100, 50)
#' )
#' adat <- AlignmentData(data_list)
#'
#' @export
AlignmentData <- function(data,
                          subjects = NULL,
                          space = NULL,
                          design = NULL,
                          geometry = NULL,
                          obs_labels = NULL,
                          metadata = list()) {
  # Validate data is a list

if (!is.list(data)) {
    stop("'data' must be a list of matrices or NeuroVec objects")
  }

  # Extract subject IDs from names if not provided
  if (is.null(subjects)) {
    subjects <- names(data)
    if (is.null(subjects)) {
      # Auto-generate subject IDs
      subjects <- paste0("sub-", sprintf("%02d", seq_along(data)))
      names(data) <- subjects
    }
  } else {
    # Ensure names match subjects
    if (length(subjects) != length(data)) {
      stop("Length of 'subjects' must match length of 'data'")
    }
    names(data) <- subjects
  }

  # Validate subjects are unique
  if (anyDuplicated(subjects)) {
    stop("Subject IDs must be unique")
  }

  # Basic validation of data elements
  for (i in seq_along(data)) {
    elem <- data[[i]]
    if (!is.matrix(elem) && !inherits(elem, "Matrix") &&
        !inherits(elem, "NeuroVec") && !inherits(elem, "NeuroVol")) {
      warning(sprintf(
        "Element '%s' is not a matrix or NeuroVec; coercion may be needed",
        subjects[i]
      ))
    }
  }

  new("AlignmentData",
    data = data,
    subjects = subjects,
    space = space,
    design = design,
    geometry = geometry,
    obs_labels = obs_labels,
    metadata = metadata
  )
}


#' Subset AlignmentData by Subject Index or ID
#'
#' @name sub-AlignmentData
#' @aliases [,AlignmentData,ANY,ANY,ANY-method
#' @param x An AlignmentData object.
#' @param i Integer indices or character subject IDs.
#'
#' @return A new AlignmentData object containing only the selected subjects.
#'
#' @export
setMethod("[", c("AlignmentData", "ANY"),
  function(x, i) {
    if (is.character(i)) {
      # Convert subject IDs to indices
      idx <- match(i, x@subjects)
      if (any(is.na(idx))) {
        missing <- i[is.na(idx)]
        stop(sprintf("Unknown subjects: %s", paste(missing, collapse = ", ")))
      }
      i <- idx
    }

    AlignmentData(
      data = x@data[i],
      subjects = x@subjects[i],
      space = x@space,
      design = x@design,
      geometry = x@geometry,
      obs_labels = x@obs_labels,
      metadata = x@metadata
    )
  }
)


#' Get Number of Subjects in AlignmentData
#'
#' @param x An AlignmentData object.
#'
#' @return Integer number of subjects.
#'
#' @export
setMethod("length", "AlignmentData", function(x) length(x@subjects))


#' Print Method for AlignmentData
#'
#' @param object An AlignmentData object.
#'
#' @export
setMethod("show", "AlignmentData", function(object) {
  cat("AlignmentData object\n")
  cat(sprintf("  Subjects: %d\n", length(object@subjects)))
  if (length(object@subjects) > 0) {
    if (length(object@subjects) <= 6) {
      cat(sprintf("    IDs: %s\n", paste(object@subjects, collapse = ", ")))
    } else {
      cat(sprintf("    IDs: %s, ... (%d more)\n",
        paste(object@subjects[1:3], collapse = ", "),
        length(object@subjects) - 3
      ))
    }

    # Show data dimensions
    first_data <- object@data[[1]]
    if (is.matrix(first_data) || inherits(first_data, "Matrix")) {
      cat(sprintf("  Data dims: %d x %d (features x observations)\n",
        nrow(first_data), ncol(first_data)
      ))
    } else {
      cat(sprintf("  Data class: %s\n", class(first_data)[1]))
    }
  }

  if (!is.null(object@space)) {
    cat(sprintf("  Space: %s\n", class(object@space)[1]))
  }
  if (!is.null(object@design)) {
    cat("  Design: present\n")
  }
  if (!is.null(object@geometry)) {
    cat("  Geometry: present\n")
  }
  if (!is.null(object@obs_labels)) {
    n <- length(object@obs_labels)
    cat(sprintf("  Observation labels: %d\n", n))
    if (n > 0) {
      if (n <= 6) {
        cat(sprintf("    Labels: %s\n", paste(object@obs_labels, collapse = ", ")))
      } else {
        cat(sprintf("    Labels: %s, ... (%d more)\n",
          paste(object@obs_labels[1:3], collapse = ", "),
          n - 3
        ))
      }
    }
  }
})


#' Coerce List to AlignmentData
#'
#' @param x A named list of matrices.
#' @param ... Additional arguments passed to AlignmentData constructor.
#'
#' @return An AlignmentData object.
#'
#' @export
as_alignment_data <- function(x, ...) {
  UseMethod("as_alignment_data")
}


#' @export
as_alignment_data.list <- function(x, ...) {
  AlignmentData(data = x, ...)
}


#' @export
as_alignment_data.AlignmentData <- function(x, ...) {
  x
}


#' Get Subject Data by ID
#'
#' @param object An AlignmentData object.
#' @param subject Subject ID (character) or index (integer).
#'
#' @return The data for the specified subject.
#'
#' @export
get_subject_data <- function(object, subject) {
  if (is.character(subject)) {
    idx <- match(subject, object@subjects)
    if (is.na(idx)) {
      stop(sprintf("Unknown subject: %s", subject))
    }
  } else {
    idx <- subject
  }
  object@data[[idx]]
}


#' Get All Data as a List
#'
#' @param object An AlignmentData object.
#'
#' @return Named list of subject data.
#'
#' @export
get_data_list <- function(object) {
  object@data
}


#' Get Observation Labels
#'
#' @param object An AlignmentData object.
#'
#' @return Observation labels (or NULL if not set).
#'
#' @export
get_obs_labels <- function(object) {
  object@obs_labels
}


#' Validate AlignmentData for Alignment
#'
#' Check that all subjects have compatible dimensions for alignment.
#'
#' @param object An AlignmentData object.
#' @param check_features Logical; if TRUE, check that all subjects have
#'   the same number of features (rows).
#' @param check_observations Logical; if TRUE, check that all subjects have
#'   the same number of observations (columns).
#'
#' @return TRUE invisibly if valid, otherwise throws an error.
#'
#' @export
validate_alignment_data <- function(object, check_features = TRUE,
                                    check_observations = FALSE,
                                    check_obs_labels = FALSE) {
  if (length(object) == 0) {
    stop("AlignmentData contains no subjects")
  }

  # If obs_labels are present, enforce basic consistency automatically.
  if (!is.null(object@obs_labels)) {
    check_obs_labels <- TRUE
    check_observations <- TRUE
  }

  # Get dimensions for each subject
  dims <- lapply(object@data, function(x) {
    if (is.matrix(x) || inherits(x, "Matrix")) {
      c(nrow(x), ncol(x))
    } else if (inherits(x, "NeuroVec")) {
      # Assume features are spatial, observations are time
      c(prod(dim(x)[1:3]), dim(x)[4])
    } else {
      c(NA, NA)
    }
  })

  dims_mat <- do.call(rbind, dims)

  if (check_features) {
    if (length(unique(dims_mat[, 1])) > 1) {
      stop(sprintf(
        "Subjects have different numbers of features: %s",
        paste(unique(dims_mat[, 1]), collapse = ", ")
      ))
    }
  }

  if (check_observations) {
    if (length(unique(dims_mat[, 2])) > 1) {
      stop(sprintf(
        "Subjects have different numbers of observations: %s",
        paste(unique(dims_mat[, 2]), collapse = ", ")
      ))
    }
  }

  if (check_obs_labels) {
    n_obs <- unique(dims_mat[, 2])
    if (length(n_obs) != 1 || is.na(n_obs)) {
      stop(
        "Cannot validate obs_labels: subjects have differing or unknown observation counts",
        call. = FALSE
      )
    }

    labs <- object@obs_labels
    if (is.null(labs)) {
      stop("check_obs_labels=TRUE but object@obs_labels is NULL", call. = FALSE)
    }

    if (!is.atomic(labs) && !is.factor(labs)) {
      stop("obs_labels must be an atomic vector (character/factor/numeric/etc.)", call. = FALSE)
    }

    if (length(labs) != n_obs) {
      stop(
        sprintf("obs_labels length mismatch: expected %d, got %d", n_obs, length(labs)),
        call. = FALSE
      )
    }

    if (any(is.na(labs))) {
      stop("obs_labels contains NA values", call. = FALSE)
    }

    # If column names are present, require they agree with obs_labels.
    data_list <- get_data_list(object)
    colnames_list <- lapply(data_list, function(x) {
      if (is.matrix(x) || inherits(x, "Matrix")) {
        colnames(x)
      } else {
        NULL
      }
    })
    has_any_names <- any(vapply(colnames_list, function(nm) !is.null(nm), logical(1)))
    has_all_names <- all(vapply(colnames_list, function(nm) !is.null(nm), logical(1)))
    if (has_any_names && !has_all_names) {
      stop("Some subjects have colnames but others do not; cannot validate against obs_labels", call. = FALSE)
    }
    if (has_all_names) {
      for (subj in names(data_list)) {
        nm <- colnames_list[[subj]]
        if (length(nm) != length(labs) || any(nm != as.character(labs))) {
          stop(
            sprintf("Subject '%s' colnames do not match obs_labels", subj),
            call. = FALSE
          )
        }
      }
    }
  }

  invisible(TRUE)
}
