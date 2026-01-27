#' AlignmentModel Class
#'
#' Represents a fitted alignment model that can be applied to new data.
#' Separates fit from apply to enable cross-validation workflows.
#'
#' @slot transforms Named list of per-subject operators (target x source dimensions).
#' @slot reference The reference used for alignment (subject ID, matrix, or "consensus").
#' @slot reference_data The actual reference data (matrix or template object).
#' @slot method Character string identifying the alignment method.
#' @slot space_from Source space specification.
#' @slot space_to Target space specification.
#' @slot provenance List containing parameters, package versions, and data hashes.
#' @slot method_state List containing method-specific state needed for apply.
#' @slot train_subjects Character vector of subjects used for training (for CV).
#'
#' @export
setClass("AlignmentModel",
  slots = c(
    transforms = "list",
    reference = "ANY",
    reference_data = "ANY",
    method = "character",
    space_from = "ANY",
    space_to = "ANY",
    provenance = "list",
    method_state = "list",
    train_subjects = "character"
  ),
  prototype = list(
    transforms = list(),
    reference = NULL,
    reference_data = NULL,
    method = character(0),
    space_from = NULL,
    space_to = NULL,
    provenance = list(),
    method_state = list(),
    train_subjects = character(0)
  )
)

setValidity("AlignmentModel", function(object) {
  errors <- character()
  if (!is.list(object@transforms)) {
    errors <- c(errors, "'transforms' must be a list")
  }
  if (length(object@transforms) > 0 && is.null(names(object@transforms))) {
    errors <- c(errors, "'transforms' must be named with subject IDs")
  }
  if (!is.character(object@method)) {
    errors <- c(errors, "'method' must be a character string")
  }
  if (!is.list(object@provenance)) {
    errors <- c(errors, "'provenance' must be a list")
  }
  if (!is.list(object@method_state)) {
    errors <- c(errors, "'method_state' must be a list")
  }
  if (!is.character(object@train_subjects)) {
    errors <- c(errors, "'train_subjects' must be a character vector")
  }
  if (length(errors) == 0L) TRUE else errors
})


#' Create an AlignmentModel Object
#'
#' @param transforms Named list of per-subject operators (target x source).
#' @param reference The reference specification (subject ID, matrix, or "consensus").
#' @param reference_data The actual reference data.
#' @param method Character string identifying the alignment method.
#' @param space_from Source space specification.
#' @param space_to Target space specification.
#' @param params List of parameters used for fitting.
#' @param method_state Method-specific state for applying to new data.
#' @param train_subjects Character vector of subjects used in training.
#' @param provenance Optional pre-built provenance list. If provided, used
#'   directly instead of auto-generating from params. Useful for composed models
#'   or models reconstructed from serialized state.
#'
#' @return An AlignmentModel object.
#'
#' @export
AlignmentModel <- function(transforms,
                           reference,
                           reference_data = NULL,
                           method,
                           space_from = NULL,
                           space_to = NULL,
                           params = list(),
                           method_state = list(),
                           train_subjects = character(0),
                           provenance = NULL) {
  # Validate required parameters
  if (missing(method) || is.null(method) || !is.character(method) || length(method) != 1L) {
    stop("'method' must be a single character string", call. = FALSE)
  }
  if (!is.list(transforms)) {
    stop("'transforms' must be a named list of operators", call. = FALSE)
  }
  if (is.null(names(transforms))) {
    stop("'transforms' must be named with subject IDs", call. = FALSE)
  }

  # Build provenance if not provided
  if (is.null(provenance)) {
    provenance <- list(
      params = params,
      fitted_at = Sys.time(),
      neuralign_version = as.character(utils::packageVersion("neuralign")),
      r_version = R.version.string,
      data_hash = NULL
    )

    # Try to get method package version
    aligner_info <- get_aligner(method)
    if (!is.null(aligner_info)) {
      pkg <- aligner_info$package
      if (!is.null(pkg) && requireNamespace(pkg, quietly = TRUE)) {
        provenance$method_package_version <- as.character(
          utils::packageVersion(pkg)
        )
      }
    }
  }

  new("AlignmentModel",
    transforms = transforms,
    reference = reference,
    reference_data = reference_data,
    method = method,
    space_from = space_from,
    space_to = space_to,
    provenance = provenance,
    method_state = method_state,
    train_subjects = train_subjects
  )
}


#' Print Method for AlignmentModel
#'
#' @param object An AlignmentModel object.
#'
#' @export
setMethod("show", "AlignmentModel", function(object) {
  cat("AlignmentModel object\n")
  cat(sprintf("  Method: %s\n", object@method))
  cat(sprintf("  Subjects: %d\n", length(object@transforms)))

  if (length(object@transforms) > 0) {
    # Show transform dimensions
    first_t <- object@transforms[[1]]
    if (is.matrix(first_t) || inherits(first_t, "Matrix")) {
      cat(sprintf("  Transform dims: %d x %d (target x source)\n",
        nrow(first_t), ncol(first_t)
      ))
    }
  }

  # Reference info
  if (is.character(object@reference) && length(object@reference) == 1) {
    if (object@reference == "consensus") {
      cat("  Reference: consensus\n")
    } else if (object@reference == "fold_specific") {
      cat("  Reference: fold-specific (no common anchor)\n")
    } else {
      cat(sprintf("  Reference: subject '%s'\n", object@reference))
    }
  } else if (.is_matrixish(object@reference)) {
    cat("  Reference: template matrix\n")
  }

  # Space info
  if (!is.null(object@space_from)) {
    cat(sprintf("  Space from: %s\n", .format_space(object@space_from)))
  }
  if (!is.null(object@space_to)) {
    cat(sprintf("  Space to: %s\n", .format_space(object@space_to)))
  }

  # Training info
  if (length(object@train_subjects) > 0) {
    cat(sprintf("  Trained on: %d subjects\n", length(object@train_subjects)))
  }

  # Provenance
  fitted_at <- object@provenance$fitted_at
  if (!is.null(fitted_at) && length(fitted_at) == 1L) {
    cat(sprintf("  Fitted at: %s\n", fitted_at))
  }
})


#' Check Space Compatibility
#'
#' Determine whether two space specifications are compatible.
#' NULL is compatible with anything.
#'
#' @param a First space specification.
#' @param b Second space specification.
#'
#' @return Logical; TRUE if spaces are compatible.
#'
#' @export
spaces_compatible <- function(a, b) {
  if (is.null(a) || is.null(b)) return(TRUE)
  if (identical(a, b)) return(TRUE)

  # Both character: string equality

  if (is.character(a) && is.character(b)) {
    return(identical(a, b))
  }

  # Both gds_space: compare by name
  if (inherits(a, "gds_space") && inherits(b, "gds_space")) {
    a_name <- a$name %||% ""
    b_name <- b$name %||% ""
    return(identical(a_name, b_name))
  }

  # Fallback: all.equal
  isTRUE(all.equal(a, b))
}


#' Helper to Format Space Info
#' @keywords internal
.format_space <- function(space) {
  if (is.character(space)) {
    space
  } else if (inherits(space, "gds_space")) {
    sprintf("gds_space (%s)", space$name %||% "unnamed")
  } else {
    class(space)[1]
  }
}


#' Get Transform for a Subject
#'
#' @param model An AlignmentModel object.
#' @param subject Subject ID.
#'
#' @return The transform operator for the subject.
#'
#' @export
get_transform <- function(model, subject) {
  if (!subject %in% names(model@transforms)) {
    stop(sprintf("Subject '%s' not found in model", subject), call. = FALSE)
  }
  model@transforms[[subject]]
}


#' Get All Transforms
#'
#' @param model An AlignmentModel object.
#'
#' @return Named list of transform operators.
#'
#' @export
get_transforms <- function(model) {
  model@transforms
}


#' Get Subjects in Model
#'
#' @param model An AlignmentModel object.
#'
#' @return Character vector of subject IDs.
#'
#' @export
model_subjects <- function(model) {
  names(model@transforms)
}


#' Check if Model Has Transform for Subject
#'
#' @param model An AlignmentModel object.
#' @param subject Subject ID.
#'
#' @return Logical.
#'
#' @export
has_transform <- function(model, subject) {
  subject %in% names(model@transforms)
}


#' Get Reference Data from Model
#'
#' @param model An AlignmentModel object.
#'
#' @return The reference data (matrix or template object).
#'
#' @export
get_reference <- function(model) {
  model@reference_data
}


#' Add Transform to Model
#'
#' Creates a new model with an additional transform (does not modify in place).
#'
#' @param model An AlignmentModel object.
#' @param subject Subject ID for the new transform.
#' @param transform The transform operator.
#'
#' @return A new AlignmentModel with the added transform.
#'
#' @export
add_transform <- function(model, subject, transform) {
  new_transforms <- model@transforms
  new_transforms[[subject]] <- transform

  new("AlignmentModel",
    transforms = new_transforms,
    reference = model@reference,
    reference_data = model@reference_data,
    method = model@method,
    space_from = model@space_from,
    space_to = model@space_to,
    provenance = model@provenance,
    method_state = model@method_state,
    train_subjects = model@train_subjects
  )
}


#' Subset AlignmentModel by Subject
#'
#' @name sub-AlignmentModel
#' @aliases [,AlignmentModel,ANY,ANY,ANY-method
#' @param x An AlignmentModel object.
#' @param i Character subject IDs or integer indices.
#'
#' @return A new AlignmentModel with only the selected subjects.
#'
#' @export
setMethod("[", c("AlignmentModel", "ANY"),
  function(x, i) {
    subj_names <- names(x@transforms)

    if (is.character(i)) {
      if (!all(i %in% subj_names)) {
        missing <- setdiff(i, subj_names)
        stop(sprintf("Unknown subjects: %s", paste(missing, collapse = ", ")), call. = FALSE)
      }
    } else {
      i <- subj_names[i]
    }

    new("AlignmentModel",
      transforms = x@transforms[i],
      reference = x@reference,
      reference_data = x@reference_data,
      method = x@method,
      space_from = x@space_from,
      space_to = x@space_to,
      provenance = x@provenance,
      method_state = x@method_state,
      train_subjects = x@train_subjects
    )
  }
)


#' Get Number of Subjects in Model
#'
#' @param x An AlignmentModel object.
#'
#' @return Integer number of subjects with transforms.
#'
#' @export
setMethod("length", "AlignmentModel", function(x) length(x@transforms))


.ensure_model <- function(x, what = "x") {
  if (inherits(x, "AlignmentResult")) {
    x <- get_model(x)
  }
  if (!inherits(x, "AlignmentModel")) {
    stop(sprintf("'%s' must be an AlignmentModel or AlignmentResult", what), call. = FALSE)
  }
  x
}


#' NULL-coalescing operator
#' @param a First argument
#' @param b Second argument (returned if a is NULL)
#' @return a if not NULL, otherwise b
#' @noRd
`%||%` <- function(a, b) if (is.null(a)) b else a
