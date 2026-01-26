#' fmrigds Compatibility
#'
#' Functions for integrating neuralign with fmrigds workflows.
#'
#' @name compat_fmrigds
NULL


#' Convert AlignmentModel to MapFamily
#'
#' Convert a neuralign AlignmentModel to an fmrigds MapFamily object
#' for use in fmrigds pipelines.
#'
#' @param model An AlignmentModel or AlignmentResult.
#' @param name Optional name for the MapFamily.
#'
#' @return A MapFamily object (if fmrigds is available).
#'
#' @details
#' The MapFamily object represents the collection of per-subject operators
#' that can be applied in fmrigds group analysis workflows. The operator
#' semantics match: left-multiply convention with (target x source) dimensions.
#'
#' @examples
#' \dontrun{
#' result <- fit_alignment(data, method = "procrustes")
#' map_family <- as_map_family(result)
#'
#' # Use in fmrigds pipeline
#' aligned <- map_family$apply(group_data)
#' }
#'
#' @export
as_map_family <- function(model, name = NULL) {
  if (!requireNamespace("fmrigds", quietly = TRUE)) {
    stop("Package 'fmrigds' required for as_map_family()")
  }

  if (inherits(model, "AlignmentResult")) {
    if (is_cv_result(model) && !has_common_anchor(model)) {
      stop(
        paste0(
          "Cannot convert a cross-validated result with fold-specific anchors to fmrigds::MapFamily(). ",
          "Fit without CV or use a fixed/external reference so all subjects share a common anchor space."
        ),
        call. = FALSE
      )
    }
    model <- get_model(model)
  }

  if (!inherits(model, "AlignmentModel")) {
    stop("'model' must be an AlignmentModel")
  }

  # Build name from method if not provided
  if (is.null(name)) {
    name <- sprintf("neuralign_%s", model@method)
  }

  if (is.null(model@space_from) || is.null(model@space_to)) {
    stop("AlignmentModel must have non-NULL space_from/space_to to convert to fmrigds::MapFamily()", call. = FALSE)
  }

  by_subject <- model@transforms

  caps <- aligner_capabilities(model@method) %||% list()
  transform_type <- caps$transform_type %||% "linear"
  fmrigds_type <- transform_type
  if (!transform_type %in% c("linear", "orthogonal", "ot", "affine3d", "deform3d")) {
    # fmrigds has a fixed set; treat unknown/permutation/etc. as linear maps
    fmrigds_type <- "linear"
  }

  traits <- list(
    orthogonal = identical(transform_type, "orthogonal"),
    mass_preserving = isTRUE(caps$mass_preserving)
  )

  fmrigds::MapFamily(
    name = name,
    from_space = model@space_from,
    to_space = model@space_to,
    type = fmrigds_type,
    by_subject = by_subject,
    traits = traits,
    uncertainty = fmrigds::UncertaintyRule("independent")
  )
}


#' Convert MapFamily to AlignmentModel
#'
#' Convert an fmrigds MapFamily to a neuralign AlignmentModel.
#'
#' @param map_family A MapFamily object from fmrigds.
#' @param method Optional method name to assign.
#'
#' @return An AlignmentModel.
#'
#' @export
from_map_family <- function(map_family, method = "fmrigds_imported") {
  if (!is.list(map_family) || is.null(map_family$by_subject)) {
    stop("Expected an fmrigds MapFamily-like object with $by_subject", call. = FALSE)
  }

  maps <- map_family$by_subject
  space_from <- map_family$from
  space_to <- map_family$to

  AlignmentModel(
    transforms = maps,
    reference = NULL,
    method = method,
    space_from = space_from,
    space_to = space_to,
    train_subjects = names(maps)
  )
}


#' Apply AlignmentModel in fmrigds Pipeline
#'
#' Helper for applying neuralign alignment within fmrigds workflows.
#'
#' @param model AlignmentModel.
#' @param gds_data A gds_data object from fmrigds.
#' @param ... Additional arguments.
#'
#' @return Transformed gds_data object.
#'
#' @export
apply_to_gds <- function(model, gds_data, ...) {
  if (!requireNamespace("fmrigds", quietly = TRUE)) {
    stop("Package 'fmrigds' required")
  }

  fam <- as_map_family(model)
  fmrigds::align(gds_data, fam)
}


#' Create AlignmentData from gds_data
#'
#' Extract data from an fmrigds gds_data object into AlignmentData format.
#'
#' @param gds_data A gds_data object from fmrigds.
#' @param ... Additional arguments passed to AlignmentData.
#' @param assay Name of assay to extract (default: "beta"). Extracted as a
#'   per-subject matrix (samples x contrasts).
#'
#' @return An AlignmentData object.
#'
#' @export
alignment_data_from_gds <- function(gds_data, assay = "beta", ...) {
  if (!requireNamespace("fmrigds", quietly = TRUE)) {
    stop("Package 'fmrigds' required")
  }

  subjects <- fmrigds::subjects(gds_data)
  space <- fmrigds::space(gds_data)

  arr <- fmrigds::assay(gds_data, name = assay)
  if (length(dim(arr)) != 3L) {
    stop("fmrigds::assay() must return a 3D array [sample x subject x contrast]", call. = FALSE)
  }
  if (dim(arr)[2] != length(subjects)) {
    stop("Assay subject dimension does not match fmrigds::subjects()", call. = FALSE)
  }

  data_list <- lapply(seq_along(subjects), function(j) {
    mat <- arr[, j, , drop = TRUE]
    if (!is.matrix(mat)) mat <- as.matrix(mat)
    mat
  })
  names(data_list) <- subjects

  AlignmentData(
    data = data_list,
    subjects = subjects,
    space = space,
    ...
  )
}


#' Print Method for neuralign_map_family
#'
#' @param x A neuralign_map_family object.
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns x.
#'
#' @export
print.neuralign_map_family <- function(x, ...) {
  cat("neuralign MapFamily\n")
  cat(sprintf("  Name: %s\n", x$name))
  cat(sprintf("  Subjects: %d\n", x$n_subjects))
  if (x$n_subjects > 0) {
    first_map <- x$maps[[1]]
    cat(sprintf("  Map dims: %d x %d\n", nrow(first_map), ncol(first_map)))
  }
  invisible(x)
}
