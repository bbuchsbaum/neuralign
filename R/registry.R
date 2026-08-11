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

.disabled_aligner_reasons <- c(
  cone = paste(
    "CONE is disabled because its upstream estimator fails neuralign's",
    "deterministic graph-correspondence accuracy contract (neuralign-tzc).",
    "Use 'grasp' for the supported graph-alignment path."
  ),
  cone_align = paste(
    "CONE is disabled because its upstream estimator fails neuralign's",
    "deterministic graph-correspondence accuracy contract (neuralign-tzc).",
    "Use 'grasp' for the supported graph-alignment path."
  )
)

.disabled_aligner_reason <- function(name) {
  if (!is.character(name) || length(name) != 1L || is.na(name) ||
      !name %in% names(.disabled_aligner_reasons)) {
    return(NULL)
  }
  unname(.disabled_aligner_reasons[[name]])
}

#' Aligner API Version
#'
#' The sole supported version of the aligner registration API. Third-party
#' packages should declare this version so incompatible provider contracts fail
#' during registration rather than during fitting.
#'
#' @examples
#' NEURALIGN_ALIGNER_API_VERSION
#'
#' @export
NEURALIGN_ALIGNER_API_VERSION <- 2L

#' Validate Aligner Contract
#'
#' Check that an aligner registration entry conforms to the expected contract.
#' This validates the fit_fn and apply_fn signatures, capabilities structure,
#' and other requirements.
#'
#' @family registry
#'
#' @param name Method name.
#' @param fit_fn The fit function to validate.
#' @param apply_fn Optional apply function to validate.
#' @param prepare_fn Optional provider preflight function to validate.
#' @param capabilities Capabilities list to validate.
#' @param api_version API version the aligner declares. neuralign accepts only
#'   \code{NEURALIGN_ALIGNER_API_VERSION}.
#'
#' @return TRUE invisibly if valid; otherwise throws an error.
#'
#' @details
#' The \code{fit_fn} must explicitly accept \code{data}, \code{reference},
#' \code{train_idx}, \code{fit_context}, and \code{provider_plan}. Requiring
#' lifecycle arguments as named formals prevents them from leaking through
#' \code{...} into an estimator backend. If supplied, \code{prepare_fn} must
#' explicitly accept \code{data}, \code{reference}, and
#' \code{resampling_plan}.
#'
#' If provided, \code{apply_fn} must accept at minimum: \code{fit_result},
#' \code{new_data}, and \code{...}. It should return a list with a
#' \code{transforms} element (named list of operators).
#'
#' @examples
#' dummy_fit <- function(data, reference, train_idx = NULL,
#'                       fit_context, provider_plan, ...) {
#'   list(
#'     transforms = list(),
#'     reference_data = NULL,
#'     space_from = NULL,
#'     space_to = NULL,
#'     method_state = list()
#'   )
#' }
#' validate_aligner_contract(
#'   name = "dummy",
#'   fit_fn = dummy_fit,
#'   capabilities = list(transform_type = "orthogonal")
#' )
#'
#' @export
validate_aligner_contract <- function(name,
                                      fit_fn,
                                      apply_fn = NULL,
                                      capabilities = list(),
                                      api_version = NEURALIGN_ALIGNER_API_VERSION,
                                      prepare_fn = NULL) {
  errors <- character(0)

  # Check API version compatibility
  api_version_valid <- is.numeric(api_version) && length(api_version) == 1L &&
    !is.na(api_version) && is.finite(api_version) &&
    api_version == as.integer(api_version)
  if (!api_version_valid) {
    errors <- c(errors, "api_version must be one integer")
  } else if (!identical(as.integer(api_version), NEURALIGN_ALIGNER_API_VERSION)) {
    errors <- c(errors, sprintf(
      "Aligner declares api_version=%d but neuralign supports only API %d",
      api_version, NEURALIGN_ALIGNER_API_VERSION
    ))
  }

  # Validate fit_fn signature
  if (!is.function(fit_fn)) {
    errors <- c(errors, "fit_fn must be a function")
  } else {
    fit_args <- names(formals(fit_fn))
    required_fit <- c(
      "data", "reference", "train_idx", "fit_context", "provider_plan"
    )
    missing_fit <- setdiff(required_fit, fit_args)
    if (length(missing_fit) > 0L) {
      errors <- c(errors, sprintf(
        "fit_fn missing required formals: %s",
        paste(missing_fit, collapse = ", ")
      ))
    }
  }

  # Validate prepare_fn signature if provided.
  if (!is.null(prepare_fn)) {
    if (!is.function(prepare_fn)) {
      errors <- c(errors, "prepare_fn must be a function or NULL")
    } else {
      prepare_args <- names(formals(prepare_fn))
      required_prepare <- c("data", "reference", "resampling_plan")
      missing_prepare <- setdiff(required_prepare, prepare_args)
      if (length(missing_prepare) > 0L) {
        errors <- c(errors, sprintf(
          "prepare_fn missing required formals: %s",
          paste(missing_prepare, collapse = ", ")
        ))
      }
    }
  }

  # Validate apply_fn signature if provided
  if (!is.null(apply_fn)) {
    if (!is.function(apply_fn)) {
      errors <- c(errors, "apply_fn must be a function or NULL")
    } else {
      apply_args <- names(formals(apply_fn))
      required_apply <- c("fit_result", "new_data")
      missing_apply <- setdiff(required_apply, apply_args)
      if (length(missing_apply) > 0 && !"..." %in% apply_args) {
        errors <- c(errors, sprintf(
          "apply_fn missing required formals: %s",
          paste(missing_apply, collapse = ", ")
        ))
      }
    }
  }

  # Validate capabilities structure
  if (!is.list(capabilities)) {
    errors <- c(errors, "capabilities must be a list")
  } else {
    # Use [["..."]] to avoid partial matching (e.g., $returns matching $returns_invertible)
    ret <- capabilities[["returns"]]
    if (!is.null(ret) && !ret %in% c("operator", "embedding")) {
      errors <- c(errors, sprintf(
        "capabilities$returns must be 'operator' or 'embedding' (got: '%s')",
        ret
      ))
    }
    tt <- capabilities[["transform_type"]]
    if (!is.null(tt) && !is.character(tt)) {
      errors <- c(errors, "capabilities$transform_type must be a character string")
    }
    cva <- capabilities[["cv_axes"]]
    if (!is.null(cva) && !is.character(cva)) {
      errors <- c(errors, "capabilities$cv_axes must be a character vector")
    }
    reft <- capabilities[["reference_types"]]
    if (!is.null(reft) && !is.character(reft)) {
      errors <- c(errors, "capabilities$reference_types must be a character vector")
    }

    ng <- capabilities[["needs_guidance"]]
    if (!is.null(ng) && (!is.logical(ng) || length(ng) != 1L || is.na(ng))) {
      errors <- c(errors, "capabilities$needs_guidance must be TRUE/FALSE (length 1)")
    }
    gt <- capabilities[["guidance_types"]]
    if (!is.null(gt) && !is.character(gt)) {
      errors <- c(errors, "capabilities$guidance_types must be a character vector or NULL")
    }
    if (!is.null(gt) && !isTRUE(ng)) {
      errors <- c(errors, "capabilities$guidance_types requires capabilities$needs_guidance=TRUE")
    }
  }

  if (length(errors) > 0) {
    stop(sprintf(
      "Aligner contract validation failed for '%s':\n  - %s",
      name, paste(errors, collapse = "\n  - ")
    ), call. = FALSE)
  }

  invisible(TRUE)
}


#' Register an Alignment Method
#'
#' Register an alignment method with neuralign's registry. This allows external
#' packages to provide alignment implementations that integrate with neuralign's
#' unified interface.
#'
#' @family registry
#'
#' @param name Character string identifying the method (e.g., "procrustes", "fugw").
#' @param fit_fn Function implementing the fit operation. See Details for signature.
#' @param apply_fn Optional function for applying to new data. If NULL, default
#'   left-multiply application is used. When provided, must accept
#'   \code{fit_result} (list with transforms, reference_data, space_from,
#'   space_to, method_state) and \code{new_data} (AlignmentData), and return
#'   a list with a \code{transforms} element.
#' @param prepare_fn Optional provider preflight function. It receives the full
#'   \code{data}, unresolved \code{reference}, normalized
#'   \code{resampling_plan}, and method-specific \code{...} once, before any
#'   fit callback. Its return value is passed to fit callbacks as
#'   \code{provider_plan}.
#' @param capabilities Named list of capability flags. See Details.
#' @param package Character string identifying the providing package.
#' @param description Brief description of the method.
#' @param version Version string for the method implementation.
#' @param api_version Integer declaring which neuralign aligner API version
#'   this method implements. The only accepted value is
#'   \code{NEURALIGN_ALIGNER_API_VERSION}.
#'
#' @details
#' The \code{fit_fn} must have the following signature:
#' \preformatted{
#' fit_fn <- function(data, reference, train_idx = NULL,
#'                    fit_context, provider_plan, ...) {
#'   # data: AlignmentData object
#'   # reference: subject_id, matrix, or "consensus"
#'   # train_idx: indices of subjects to use for fitting (for CV)
#'   # fit_context: exact subjects and observations assigned to this fit
#'   # provider_plan: result returned by prepare_fn, or NULL
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
#'   \item{cv_axes}{Which axes support CV: "subject", "observation"}
#'   \item{needs_geometry}{Requires adjacency/graph in AlignmentData}
#'   \item{needs_design}{Requires task structure in AlignmentData}
#'   \item{needs_guidance}{Requires guidance channels (see \code{\link{set_guidance}}).}
#'   \item{guidance_types}{Optional character vector of guidance \code{$type} values
#'     required for each subject when \code{needs_guidance=TRUE} (e.g., \code{"projector"}).}
#'   \item{returns_invertible}{Transform has exact inverse}
#'   \item{transform_type}{Type: "orthogonal", "linear", "ot", "permutation"}
#'   \item{mass_preserving}{For OT: does transport preserve mass?}
#'   \item{returns}{What the method returns: "operator" or "embedding".}
#'   \item{supports_new_subject}{Can compute operator for new subject}
#'   \item{supports_new_data}{Can apply existing operator to new data}
#'   \item{reference_types}{Character vector of supported reference types}
#' }
#'
#' @return Invisibly returns TRUE.
#'
#' @examples
#' # Register a trivial identity aligner
#' if (is_aligner_registered("identity_example")) {
#'   unregister_aligner("identity_example")
#' }
#'
#' identity_fit <- function(data, reference, train_idx = NULL,
#'                          fit_context, provider_plan, ...) {
#'   transforms <- lapply(data@subjects, function(s) {
#'     diag(nrow(get_subject_data(data, s)))
#'   })
#'   names(transforms) <- data@subjects
#'   list(
#'     transforms = transforms,
#'     reference_data = NULL,
#'     space_from = data@space,
#'     space_to = data@space,
#'     method_state = list()
#'   )
#' }
#'
#' register_aligner(
#'   name = "identity_example",
#'   fit_fn = identity_fit,
#'   capabilities = list(transform_type = "orthogonal")
#' )
#'
#' set.seed(1)
#' adat <- AlignmentData(list(
#'   s1 = matrix(rnorm(12), 3, 4),
#'   s2 = matrix(rnorm(12), 3, 4)
#' ))
#' res <- fit_alignment(adat, method = "identity_example", reference = "s1", compute_quality = FALSE)
#' stopifnot(inherits(res, "AlignmentResult"))
#'
#' unregister_aligner("identity_example")
#'
#' @export
register_aligner <- function(name,
                             fit_fn,
                             apply_fn = NULL,
                             capabilities = list(),
                             package = NA_character_,
                             description = "",
                             version = "0.0.0",
                             api_version = NEURALIGN_ALIGNER_API_VERSION,
                             prepare_fn = NULL) {
  # Validate name
  if (!is.character(name) || length(name) != 1) {
    stop("'name' must be a single character string", call. = FALSE)
  }

  # Validate aligner contract
  validate_aligner_contract(
    name = name,
    fit_fn = fit_fn,
    apply_fn = apply_fn,
    prepare_fn = prepare_fn,
    capabilities = capabilities,
    api_version = api_version
  )
  api_version <- as.integer(api_version)

  # Set default capabilities
  default_caps <- list(
    # CV support
    supports_cv = FALSE,
    cv_axes = c("subject"),

    # Data requirements
    needs_geometry = FALSE,
    needs_design = FALSE,
    needs_guidance = FALSE,
    guidance_types = NULL,

    # Transform properties
    returns_invertible = FALSE,
    transform_type = "linear",  # "orthogonal", "linear", "ot", "permutation"
    mass_preserving = FALSE,
    returns = "operator",  # "operator" or "embedding"

    # Apply semantics
    supports_new_subject = TRUE,   # Can compute operator for new subject
    supports_new_data = TRUE,      # Can apply existing operator to new data

    # Reference options
    reference_types = c("subject", "consensus")
  )
  capabilities <- modifyList(default_caps, capabilities)

  if (identical(capabilities[["returns"]], "embedding") &&
      !isFALSE(capabilities[["supports_new_data"]])) {
    stop(
      sprintf(
        "Aligner '%s' declares returns='embedding' but supports_new_data is not FALSE. ",
        name
      ),
      "Embedding-returning aligners must set capabilities$supports_new_data = FALSE.",
      call. = FALSE
    )
  }

  # Create registration entry
  entry <- list(
    name = name,
    fit_fn = fit_fn,
    apply_fn = apply_fn,
    prepare_fn = prepare_fn,
    capabilities = capabilities,
    package = package,
    description = description,
    version = version,
    api_version = api_version,
    registered_at = Sys.time()
  )

  # Warn on overwrite
  if (exists(name, envir = .aligner_registry, inherits = FALSE)) {
    warning(sprintf("Overwriting existing aligner registration '%s'", name),
      call. = FALSE
    )
  }

  # Store in registry
  .aligner_registry[[name]] <- entry

  invisible(TRUE)
}


#' List Available Alignment Methods
#'
#' @family registry
#'
#' @param details Logical; if TRUE, return detailed information as a data frame.
#'
#' @return If `details = FALSE` (default), a character vector of method names.
#'   If `details = TRUE`, a data frame with one row per method and columns:
#'   \describe{
#'     \item{name}{Aligner name.}
#'     \item{package}{Providing package (if known).}
#'     \item{description}{Short description.}
#'     \item{transform_type}{One of `"orthogonal"`, `"linear"`, `"ot"`, ...}
#'     \item{supports_cv}{Logical; whether the aligner declares CV support.}
#'   }
#'
#' @examples
#' available_aligners()
#' available_aligners(details = TRUE)
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

#' @rdname available_aligners
#' @export
list_aligners <- function(details = FALSE) {
  available_aligners(details = details)
}


#' Get Information About an Aligner
#'
#' @family registry
#'
#' @param name Character string identifying the method.
#'
#' @return List with aligner information, or NULL if not found.
#'
#' @examples
#' get_aligner("procrustes")$description
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
#' @family registry
#'
#' @param name Character string identifying the method.
#'
#' @return List of capability flags, or NULL if method not found.
#'
#' @examples
#' aligner_capabilities("procrustes")$transform_type
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
#' @family registry
#'
#' @param name Character string identifying the method.
#'
#' @return Logical indicating if the aligner is registered.
#'
#' @examples
#' is_aligner_registered("procrustes")
#'
#' @export
is_aligner_registered <- function(name) {
  exists(name, envir = .aligner_registry)
}


#' Unregister an Aligner
#'
#' Remove an aligner from the registry. Primarily useful for testing.
#'
#' @family registry
#'
#' @param name Character string identifying the method.
#'
#' @return Invisibly returns TRUE if removed, FALSE if not found.
#'
#' @examples
#' tmp_fit <- function(data, reference, train_idx = NULL,
#'                     fit_context, provider_plan, ...) {
#'   list(
#'     transforms = list(),
#'     reference_data = NULL,
#'     space_from = NULL,
#'     space_to = NULL,
#'     method_state = list()
#'   )
#' }
#' if (is_aligner_registered("tmp_registry_example")) {
#'   unregister_aligner("tmp_registry_example")
#' }
#' register_aligner(
#'   name = "tmp_registry_example",
#'   fit_fn = tmp_fit,
#'   capabilities = list(transform_type = "orthogonal")
#' )
#' unregister_aligner("tmp_registry_example")
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

  # Known method -> provider mappings for lazy loading. Most providers
  # register from .onLoad. A provider may instead expose one explicit,
  # idempotent registration hook; neuralign invokes that hook only when the
  # requested method is still absent after the namespace has loaded.
  known_methods <- list(
    "fugw"          = list(package = "topofmri"),
    "gw"            = list(package = "manifoldalign"),
    "fpgw"          = list(package = "manifoldalign"),
    "kema"          = list(package = "manifoldalign"),
    "coupled_diag"  = list(package = "manifoldalign"),
    "coupled_diagonalization" = list(package = "manifoldalign"),
    "gpca"          = list(package = "manifoldalign"),
    "grasp"         = list(package = "manifoldalign"),
    "lowrank"       = list(package = "manifoldalign"),
    "dkge"          = list(package = "dkge"),
    "nef"           = list(package = "fmrireg.gnef"),
    "a_corrca"      = list(
      package = "xpar",
      register = "register_xpar_neuralign"
    )
  )

  provider <- known_methods[[name]]
  if (is.null(provider) ||
      !requireNamespace(provider$package, quietly = TRUE)) {
    return(FALSE)
  }

  registration <- provider$register %||% NULL
  if (!is_aligner_registered(name) && !is.null(registration)) {
    exports <- getNamespaceExports(provider$package)
    if (!registration %in% exports) {
      stop(sprintf(
        "Provider package '%s' does not export its declared registration hook '%s'",
        provider$package, registration
      ), call. = FALSE)
    }
    getExportedValue(provider$package, registration)()
  }

  is_aligner_registered(name)
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
    stop(sprintf("Unknown aligner: %s", name), call. = FALSE)
  }

  # Check geometry requirement
  if (isTRUE(caps$needs_geometry)) {
    if (is.null(data@geometry)) {
      stop(sprintf(
        "Method '%s' requires geometry (adjacency matrix) in AlignmentData",
        name
      ), call. = FALSE)
    }
  }

  # Check design requirement
  if (isTRUE(caps$needs_design)) {
    if (is.null(data@design)) {
      stop(sprintf(
        "Method '%s' requires design (task structure) in AlignmentData",
        name
      ), call. = FALSE)
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
