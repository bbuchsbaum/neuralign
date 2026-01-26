#' Aligner Registry
#'
#' Functions for registering and querying alignment methods.
#' External packages can register their methods via \code{register_aligner()}.
#'
#' @name registry
#' @importFrom utils modifyList
NULL


# Internal registry environment
.aligner_registry <- new.env(parent = emptyenv())


#' Register an Alignment Method
#'
#' Register an alignment method with neuralign's registry. This allows external
#' packages to provide alignment implementations that integrate with neuralign's
#' unified interface.
#'
#' @param name Character string identifying the method (e.g., "procrustes", "fugw").
#' @param fit_fn Function implementing the fit operation. See Details for signature.
#' @param apply_fn Optional function for applying to new data. If NULL, default
#'   left-multiply application is used.
#' @param capabilities Named list of capability flags. See Details.
#' @param package Character string identifying the providing package.
#' @param description Brief description of the method.
#' @param version Version string for the method implementation.
#'
#' @details
#' The \code{fit_fn} must have the following signature:
#' \preformatted{
#' fit_fn <- function(data, reference, train_idx = NULL, ...) {
#'   # data: AlignmentData object
#'   # reference: subject_id, matrix, or "consensus"
#'   # train_idx: indices of subjects to use for fitting (for CV)
#'   # ...: method-specific parameters
#'
#'   list(
#'     transforms = list(...),      # Named list: subject_id -> operator
#'     reference_data = ...,        # The actual reference
#'     space_from = ...,            # Source space
#'     space_to = ...,              # Target space
#'     method_state = list()        # Method-specific state for apply
#'   )
#' }
#' }
#'
#' Capability flags:
#' \describe{
#'   \item{supports_cv}{Can handle train/test splits}
#'   \item{cv_axes}{Which axes support CV: "subject", "run", "task"}
#'   \item{needs_geometry}{Requires adjacency/graph in AlignmentData}
#'   \item{needs_design}{Requires task structure in AlignmentData}
#'   \item{returns_invertible}{Transform has exact inverse}
#'   \item{transform_type}{Type: "orthogonal", "linear", "ot", "permutation"}
#'   \item{mass_preserving}{For OT: does transport preserve mass?}
#'   \item{returns}{What the method returns. Currently "operator" only; "embedding" is reserved}
#'   \item{supports_new_subject}{Can compute operator for new subject}
#'   \item{supports_new_data}{Can apply existing operator to new data}
#'   \item{reference_types}{Character vector of supported reference types}
#' }
#'
#' @return Invisibly returns TRUE.
#'
#' @examples
#' \dontrun{
#' register_aligner(
#'   name = "my_method",
#'   fit_fn = my_fit_function,
#'   capabilities = list(
#'     supports_cv = TRUE,
#'     needs_geometry = FALSE,
#'     transform_type = "orthogonal"
#'   ),
#'   package = "mypackage",
#'   description = "My custom alignment method"
#' )
#' }
#'
#' @export
register_aligner <- function(name,
                             fit_fn,
                             apply_fn = NULL,
                             capabilities = list(),
                             package = NA_character_,
                             description = "",
                             version = "0.0.0") {
  # Validate name
  if (!is.character(name) || length(name) != 1) {
    stop("'name' must be a single character string")
  }

  # Validate fit_fn
  if (!is.function(fit_fn)) {
    stop("'fit_fn' must be a function")
  }

  # Validate apply_fn if provided
  if (!is.null(apply_fn) && !is.function(apply_fn)) {
    stop("'apply_fn' must be a function or NULL")
  }

  # Set default capabilities
  default_caps <- list(
    # CV support
    supports_cv = FALSE,
    cv_axes = c("subject"),  # which axes support CV: "subject", "run", "task"

    # Data requirements
    needs_geometry = FALSE,
    needs_design = FALSE,

    # Transform properties
    returns_invertible = FALSE,
    transform_type = "linear",  # "orthogonal", "linear", "ot", "permutation"
    mass_preserving = FALSE,
    returns = "operator",  # "operator" (embedding reserved; see below)

    # Apply semantics
    supports_new_subject = TRUE,   # Can compute operator for new subject
    supports_new_data = TRUE,      # Can apply existing operator to new data

    # Reference options
    reference_types = c("subject", "consensus")
  )
  capabilities <- modifyList(default_caps, capabilities)

  # Validate capabilities
  if (!is.character(capabilities$returns) || length(capabilities$returns) != 1L) {
    stop("capabilities$returns must be a single string", call. = FALSE)
  }
  if (!capabilities$returns %in% c("operator", "embedding")) {
    stop(
      sprintf(
        "capabilities$returns must be one of: operator, embedding (got: %s)",
        capabilities$returns
      ),
      call. = FALSE
    )
  }
  if (identical(capabilities$returns, "embedding")) {
    stop(
      "neuralign currently supports operator-returning aligners only. ",
      "If your method produces embeddings, expose them as (target x source) operators (e.g., projections) ",
      "and store any extra embedding artifacts in method_state. Set capabilities$returns = 'operator'.",
      call. = FALSE
    )
  }

  # Create registration entry
  entry <- list(
    name = name,
    fit_fn = fit_fn,
    apply_fn = apply_fn,
    capabilities = capabilities,
    package = package,
    description = description,
    version = version,
    registered_at = Sys.time()
  )

  # Store in registry
  .aligner_registry[[name]] <- entry

  invisible(TRUE)
}


#' List Available Alignment Methods
#'
#' @param details Logical; if TRUE, return detailed information as a data frame.
#'
#' @return Character vector of method names, or data frame with details.
#'
#' @export
available_aligners <- function(details = FALSE) {
  names_list <- ls(.aligner_registry)

  if (length(names_list) == 0) {
    if (details) {
      return(data.frame(
        name = character(0),
        package = character(0),
        description = character(0),
        stringsAsFactors = FALSE
      ))
    }
    return(character(0))
  }

  if (!details) {
    return(names_list)
  }

  # Build details data frame
  info <- lapply(names_list, function(nm) {
    entry <- .aligner_registry[[nm]]
    data.frame(
      name = nm,
      package = entry$package %||% NA_character_,
      description = entry$description %||% "",
      transform_type = entry$capabilities$transform_type %||% "unknown",
      supports_cv = entry$capabilities$supports_cv %||% FALSE,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, info)
}


#' Get Information About an Aligner
#'
#' @param name Character string identifying the method.
#'
#' @return List with aligner information, or NULL if not found.
#'
#' @export
get_aligner <- function(name) {
  if (!exists(name, envir = .aligner_registry)) {
    return(NULL)
  }
  .aligner_registry[[name]]
}


#' Get Aligner Capabilities
#'
#' @param name Character string identifying the method.
#'
#' @return List of capability flags, or NULL if method not found.
#'
#' @export
aligner_capabilities <- function(name) {
  entry <- get_aligner(name)
  if (is.null(entry)) {
    return(NULL)
  }
  entry$capabilities
}


#' Check if Aligner is Registered
#'
#' @param name Character string identifying the method.
#'
#' @return Logical indicating if the aligner is registered.
#'
#' @export
is_aligner_registered <- function(name) {
  exists(name, envir = .aligner_registry)
}


#' Unregister an Aligner
#'
#' Remove an aligner from the registry. Primarily useful for testing.
#'
#' @param name Character string identifying the method.
#'
#' @return Invisibly returns TRUE if removed, FALSE if not found.
#'
#' @export
unregister_aligner <- function(name) {
  if (!exists(name, envir = .aligner_registry)) {
    return(invisible(FALSE))
  }
  rm(list = name, envir = .aligner_registry)
  invisible(TRUE)
}


#' Attempt to Auto-Load an Aligner
#'
#' When a method is requested but not registered, attempt to load the
#' package that provides it (triggers .onLoad registration).
#'
#' @param name Character string identifying the method.
#'
#' @return Logical indicating if the aligner is now available.
#'
#' @keywords internal
.try_autoload_aligner <- function(name) {
  # Already registered?
  if (is_aligner_registered(name)) {
    return(TRUE)
  }

  # Known method -> package mappings for lazy loading
  known_methods <- list(
    "fugw" = "topofmri",
    "gw" = "manifoldalign",
    "fpgw" = "manifoldalign",
    "kema" = "manifoldalign",
    "dkge" = "dkge",
    "nef" = "fmrireg.gnef"
  )

  pkg <- known_methods[[name]]
  if (!is.null(pkg)) {
    if (requireNamespace(pkg, quietly = TRUE)) {
      # Package loading should trigger registration
      return(is_aligner_registered(name))
    }
  }

  FALSE
}


#' Validate Aligner Requirements
#'
#' Check that data meets the requirements of an aligner.
#'
#' @param name Aligner name.
#' @param data AlignmentData object.
#'
#' @return TRUE if valid, otherwise throws an error.
#'
#' @keywords internal
.validate_aligner_requirements <- function(name, data) {
  caps <- aligner_capabilities(name)
  if (is.null(caps)) {
    stop(sprintf("Unknown aligner: %s", name))
  }

  # Check geometry requirement
  if (isTRUE(caps$needs_geometry)) {
    if (is.null(data@geometry)) {
      stop(sprintf(
        "Method '%s' requires geometry (adjacency matrix) in AlignmentData",
        name
      ))
    }
  }

  # Check design requirement
  if (isTRUE(caps$needs_design)) {
    if (is.null(data@design)) {
      stop(sprintf(
        "Method '%s' requires design (task structure) in AlignmentData",
        name
      ))
    }
  }

  # Check guidance requirement (anatomy/geometry priors)
  if (isTRUE(caps$needs_guidance)) {
    guidance <- tryCatch(get_guidance(data), error = function(e) NULL)
    has_any <- !is.null(guidance) && any(vapply(guidance, function(chs) length(chs %||% list()) > 0L, logical(1)))
    if (!has_any) {
      stop(sprintf(
        "Method '%s' requires guidance channels (set via set_guidance())",
        name
      ), call. = FALSE)
    }

    if (!is.null(caps$guidance_types)) {
      types <- as.character(caps$guidance_types)
      for (subj in data@subjects) {
        chs <- guidance[[subj]] %||% list()
        ok <- any(vapply(chs, function(ch) is.list(ch) && !is.null(ch$type) && ch$type %in% types, logical(1)))
        if (!ok) {
          stop(sprintf(
            "Method '%s' requires guidance types {%s}; subject '%s' has none",
            name, paste(types, collapse = ", "), subj
          ), call. = FALSE)
        }
      }
    }
  }

  TRUE
}


#' Clear All Registered Aligners
#'
#' Remove all aligners from the registry. Primarily useful for testing.
#'
#' @return Invisibly returns TRUE.
#'
#' @keywords internal
.clear_registry <- function() {
  rm(list = ls(.aligner_registry), envir = .aligner_registry)
  invisible(TRUE)
}
