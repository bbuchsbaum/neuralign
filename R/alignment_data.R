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
#' @family alignment_data
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
#'   that observations are comparable across subjects. Duplicate labels are
#'   allowed; when labels are used for matching, duplicates are matched by
#'   occurrence order within each label.
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

.alignment_data_is_supported_element <- function(x) {
  .is_matrixish(x) || inherits(x, "NeuroVec") || inherits(x, "NeuroVol")
}

.alignment_data_element_dims <- function(x) {
  if (.is_matrixish(x)) {
    return(c(nrow(x), ncol(x)))
  }
  if (inherits(x, "NeuroVec") || inherits(x, "NeuroVol")) {
    d <- dim(x)
    if (is.null(d) || !length(d)) {
      return(c(NA_integer_, NA_integer_))
    }
    if (length(d) >= 3L) {
      n_features <- prod(d[seq_len(3L)])
    } else {
      n_features <- d[[1L]]
    }
    n_obs <- if (length(d) >= 4L) d[[4L]] else 1L
    return(c(as.integer(n_features), as.integer(n_obs)))
  }
  c(NA_integer_, NA_integer_)
}

setValidity("AlignmentData", function(object) {
  errors <- character()

  if (!is.list(object@data)) {
    errors <- c(errors, "'data' must be a list")
  }
  if (!is.character(object@subjects)) {
    errors <- c(errors, "'subjects' must be a character vector")
  }

  if (length(object@data) != length(object@subjects)) {
    errors <- c(errors, "length of 'data' must match length of 'subjects'")
  }
  if (anyDuplicated(object@subjects)) {
    errors <- c(errors, "subject IDs must be unique")
  }
  if (length(object@subjects) > 0 && any(!nzchar(object@subjects))) {
    errors <- c(errors, "subject IDs must be non-empty strings")
  }
  if (!is.list(object@metadata)) {
    errors <- c(errors, "'metadata' must be a list")
  }

  if (length(errors) > 0L) {
    return(errors)
  }

  if (length(object@data) > 0L) {
    if (is.null(names(object@data)) || any(!nzchar(names(object@data)))) {
      errors <- c(errors, "'data' must be a named list keyed by subject")
    } else if (!identical(names(object@data), object@subjects)) {
      errors <- c(errors, "names(data) must exactly match 'subjects'")
    }
  }

  dims_mat <- NULL
  if (length(object@data) > 0L && length(errors) == 0L) {
    dims <- vector("list", length(object@data))
    for (i in seq_along(object@data)) {
      subj <- object@subjects[[i]]
      x <- object@data[[i]]

      if (!.alignment_data_is_supported_element(x)) {
        errors <- c(
          errors,
          sprintf(
            "data element for subject '%s' must be matrix-like or a NeuroVec (got class '%s')",
            subj, class(x)[1L]
          )
        )
        next
      }

      d <- .alignment_data_element_dims(x)
      if (any(is.na(d))) {
        errors <- c(
          errors,
          sprintf(
            "data element for subject '%s' must have valid dimensions (got NA dims for class '%s')",
            subj, class(x)[1L]
          )
        )
        next
      }
      if (any(d < 1L)) {
        errors <- c(
          errors,
          sprintf(
            "data element for subject '%s' must have positive dimensions (got %d x %d)",
            subj, d[[1L]], d[[2L]]
          )
        )
      }
      dims[[i]] <- d
    }

    if (length(errors) == 0L) {
      dims_mat <- do.call(rbind, dims)
      rownames(dims_mat) <- object@subjects
    }
  }

  labs <- object@obs_labels
  if (!is.null(labs)) {
    if (length(object@data) == 0L) {
      errors <- c(errors, "obs_labels set but no subjects are present in 'data'")
    } else if (is.null(dims_mat)) {
      errors <- c(errors, "obs_labels set but could not determine data dimensions")
    } else if (is.atomic(labs) || is.factor(labs)) {
      if (any(is.na(labs))) {
        errors <- c(errors, "obs_labels contains NA values")
      }
      n_obs <- dims_mat[, 2L]
      bad <- names(n_obs)[n_obs != length(labs)]
      if (length(bad) > 0L) {
        errors <- c(
          errors,
          sprintf(
            "obs_labels length mismatch: expected %d columns for subjects %s, got %d labels",
            unique(n_obs[bad])[[1L]],
            paste(bad, collapse = ", "),
            length(labs)
          )
        )
      }

      # If column names are present, require they agree with obs_labels and are
      # consistently present across subjects.
      colnames_list <- lapply(object@data, function(x) {
        if (.is_matrixish(x)) colnames(x) else NULL
      })
      has_any_names <- any(vapply(colnames_list, function(nm) !is.null(nm), logical(1)))
      has_all_names <- all(vapply(colnames_list, function(nm) !is.null(nm), logical(1)))
      if (has_any_names && !has_all_names) {
        errors <- c(errors, "Some subjects have colnames but others do not; cannot validate against obs_labels")
      }
      if (has_all_names) {
        labs_chr <- as.character(labs)
        for (subj in object@subjects) {
          nm <- colnames_list[[subj]]
          if (!identical(nm, labs_chr)) {
            errors <- c(errors, sprintf("Subject '%s' colnames do not match obs_labels", subj))
          }
        }
      }
    } else if (is.list(labs)) {
      if (is.null(names(labs)) || any(!nzchar(names(labs)))) {
        errors <- c(errors, "obs_labels list must be a named list keyed by subject")
      } else {
        missing <- setdiff(object@subjects, names(labs))
        if (length(missing) > 0L) {
          errors <- c(
            errors,
            sprintf("obs_labels list is missing subjects: %s", paste(missing, collapse = ", "))
          )
        }
      }

      if (length(errors) == 0L) {
        for (subj in object@subjects) {
          v <- labs[[subj]]
          if (!is.atomic(v) && !is.factor(v)) {
            errors <- c(errors, sprintf("obs_labels for subject '%s' must be an atomic vector", subj))
            next
          }
          v <- as.character(v)
          if (any(is.na(v))) {
            errors <- c(errors, sprintf("obs_labels for subject '%s' contains NA values", subj))
          }
          n_obs_subj <- dims_mat[subj, 2L]
          if (length(v) != n_obs_subj) {
            errors <- c(
              errors,
              sprintf(
                "obs_labels length mismatch for subject '%s': expected %d, got %d",
                subj, n_obs_subj, length(v)
              )
            )
          }
          x <- object@data[[subj]]
          if (.is_matrixish(x)) {
            nm <- colnames(x)
            if (!is.null(nm) && !identical(nm, v)) {
              errors <- c(errors, sprintf("Subject '%s' colnames do not match its obs_labels", subj))
            }
          }
        }
      }
    } else {
      errors <- c(errors, "obs_labels must be NULL, an atomic vector, or a named list")
    }
  }

  if (length(errors) == 0L) TRUE else errors
})

.ensure_alignment_data <- function(x, what = "object") {
  if (!inherits(x, "AlignmentData")) {
    stop(sprintf("'%s' must be an AlignmentData object", what), call. = FALSE)
  }
  x
}

.set_obs_labels_if_missing <- function(data, obs_labels) {
  data <- .ensure_alignment_data(data, what = "data")
  if (is.null(obs_labels)) return(data)

  if (!is.null(data@obs_labels) &&
      !identical(as.character(data@obs_labels), as.character(obs_labels))) {
    stop("obs_labels supplied but AlignmentData already has different obs_labels", call. = FALSE)
  }
  if (is.null(data@obs_labels)) {
    data@obs_labels <- obs_labels
    validObject(data)
  }

  data
}

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
#' @family alignment_data
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
    stop("'data' must be a list of matrices or NeuroVec objects", call. = FALSE)
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
      stop("Length of 'subjects' must match length of 'data'", call. = FALSE)
    }
    names(data) <- subjects
  }

  # Validate subjects are unique
  if (anyDuplicated(subjects)) {
    stop("Subject IDs must be unique", call. = FALSE)
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
    subj_sel <- .resolve_subject_subset(i, x@subjects, what = "subjects")
    i <- match(subj_sel, x@subjects)

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
    max_ids_inline <- 6L
    n_ids_preview <- 3L

    if (length(object@subjects) <= max_ids_inline) {
      cat(sprintf("    IDs: %s\n", paste(object@subjects, collapse = ", ")))
    } else {
      cat(sprintf("    IDs: %s, ... (%d more)\n",
        paste(object@subjects[seq_len(n_ids_preview)], collapse = ", "),
        length(object@subjects) - n_ids_preview
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
#' @family alignment_data
#'
#' @examples
#' set.seed(1)
#' x <- list(
#'   s1 = matrix(rnorm(12), 3, 4),
#'   s2 = matrix(rnorm(12), 3, 4)
#' )
#' adat <- as_alignment_data(x)
#' stopifnot(inherits(adat, "AlignmentData"))
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
#' @family alignment_data
#'
#' @examples
#' set.seed(1)
#' adat <- AlignmentData(list(
#'   s1 = matrix(rnorm(12), 3, 4),
#'   s2 = matrix(rnorm(12), 3, 4)
#' ))
#' x1 <- get_subject_data(adat, "s1")
#' dim(x1)
#'
#' @export
get_subject_data <- function(object, subject) {
  object <- .ensure_alignment_data(object, what = "object")

  subj_ids <- object@subjects
  subj <- .resolve_subject_subset(subject, subj_ids, what = "subject")
  if (length(subj) != 1L) {
    stop("'subject' must select exactly one subject", call. = FALSE)
  }
  idx <- match(subj, subj_ids)
  object@data[[idx]]
}


#' Get All Data as a List
#'
#' @param object An AlignmentData object.
#'
#' @return Named list of subject data.
#'
#' @family alignment_data
#'
#' @examples
#' set.seed(1)
#' adat <- AlignmentData(list(
#'   s1 = matrix(rnorm(12), 3, 4),
#'   s2 = matrix(rnorm(12), 3, 4)
#' ))
#' dl <- get_data_list(adat)
#' names(dl)
#'
#' @export
get_data_list <- function(object) {
  object <- .ensure_alignment_data(object, what = "object")
  object@data
}


#' Get Observation Labels
#'
#' @param object An AlignmentData object.
#'
#' @return Observation labels (or NULL if not set).
#'
#' @family alignment_data
#'
#' @examples
#' adat <- AlignmentData(
#'   list(s1 = matrix(1, 3, 4), s2 = matrix(1, 3, 4)),
#'   obs_labels = paste0("t", 1:4)
#' )
#' get_obs_labels(adat)
#'
#' @export
get_obs_labels <- function(object) {
  object <- .ensure_alignment_data(object, what = "object")
  object@obs_labels
}


.validate_features <- function(dims_mat) {
  if (length(unique(dims_mat[, 1])) > 1) {
    stop(sprintf(
      "Subjects have different numbers of features: %s",
      paste(unique(dims_mat[, 1]), collapse = ", ")
    ), call. = FALSE)
  }
  invisible(TRUE)
}

.validate_observations <- function(dims_mat) {
  if (length(unique(dims_mat[, 2])) > 1) {
    stop(sprintf(
      "Subjects have different numbers of observations: %s",
      paste(unique(dims_mat[, 2]), collapse = ", ")
    ), call. = FALSE)
  }
  invisible(TRUE)
}

.validate_obs_labels_atomic <- function(labels, dims_mat) {
  # Shared observation axis
  n_obs <- dims_mat[1, 2]
  if (length(labels) != n_obs) {
    stop(sprintf(
      "AlignmentData@obs_labels length mismatch: expected %d (n_obs), got %d",
      as.integer(n_obs), length(labels)
    ), call. = FALSE)
  }
  if (any(is.na(labels))) {
    stop("AlignmentData@obs_labels must not contain NA", call. = FALSE)
  }
  invisible(TRUE)
}

.validate_obs_labels_list <- function(labels, subjects, dims_mat, data_names) {
  if (is.null(names(labels))) {
    if (length(labels) != length(subjects)) {
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

  for (subj in subjects) {
    labs_i <- labels[[subj]]
    if (!(is.atomic(labs_i) || is.factor(labs_i))) {
      stop(
        sprintf("AlignmentData@obs_labels[[%s]] must be an atomic vector or factor", subj),
        call. = FALSE
      )
    }
    if (any(is.na(labs_i))) {
      stop(
        sprintf("AlignmentData@obs_labels[[%s]] must not contain NA", subj),
        call. = FALSE
      )
    }

    n_obs_i <- dims_mat[match(subj, data_names), 2]
    if (length(labs_i) != n_obs_i) {
      stop(sprintf(
        "AlignmentData@obs_labels[[%s]] length mismatch: expected %d (n_obs), got %d",
        subj, as.integer(n_obs_i), length(labs_i)
      ), call. = FALSE)
    }
  }

  invisible(TRUE)
}

#' Validate AlignmentData for Alignment
#'
#' Check that all subjects have compatible dimensions for alignment.
#'
#' @param object An [AlignmentData] object.
#' @param check_features Logical; if `TRUE`, check that all subjects have the
#'   same number of features (rows).
#' @param check_observations Logical; if `TRUE`, check that all subjects have
#'   the same number of observations (columns).
#' @param check_obs_labels Logical; if `TRUE`, validate `obs_labels`. For atomic
#'   `obs_labels`, enforces length == n_obs (and requires a shared observation
#'   axis). For list-valued `obs_labels`, validates per-subject label vectors
#'   match each subject's observation count.
#'
#' @return `TRUE` invisibly if valid, otherwise throws an error.
#'
#' @family alignment_data
#'
#' @export
validate_alignment_data <- function(object, check_features = TRUE,
                                    check_observations = FALSE,
                                    check_obs_labels = FALSE) {
  if (!inherits(object, "AlignmentData")) {
    stop("'object' must be an AlignmentData instance", call. = FALSE)
  }

  validObject(object)

  if (length(object) == 0) {
    stop("AlignmentData contains no subjects", call. = FALSE)
  }

  # If obs_labels are present, enforce basic consistency automatically.
  # Atomic obs_labels imply a shared observation axis; list obs_labels are
  # per-subject and should not force matching observation counts.
  if (!is.null(object@obs_labels) && (is.atomic(object@obs_labels) || is.factor(object@obs_labels))) {
    check_obs_labels <- TRUE
    check_observations <- TRUE
  } else if (!is.null(object@obs_labels) && is.list(object@obs_labels)) {
    check_obs_labels <- TRUE
  }

  dims <- lapply(object@data, .alignment_data_element_dims)
  dims_mat <- do.call(rbind, dims)

  if (check_features) {
    .validate_features(dims_mat)
  }

  if (check_observations) {
    .validate_observations(dims_mat)
  }

  if (check_obs_labels) {
    if (is.null(object@obs_labels)) {
      stop("check_obs_labels=TRUE but object@obs_labels is NULL", call. = FALSE)
    }

    labels <- object@obs_labels

    if (is.atomic(labels) || is.factor(labels)) {
      .validate_obs_labels_atomic(labels, dims_mat)
    } else if (is.list(labels)) {
      .validate_obs_labels_list(labels, object@subjects, dims_mat, names(object@data))
    } else {
      stop(
        "AlignmentData@obs_labels must be NULL, an atomic vector/factor, or a per-subject (named) list",
        call. = FALSE
      )
    }
  }

  invisible(TRUE)
}
