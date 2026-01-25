#' Compose Alignment Models
#'
#' Compose two alignment models into a single model that applies both
#' transformations. Useful for combining geometric (spatial) and functional
#' alignments into a single pipeline.
#'
#' @param model1 First AlignmentModel (applied first, closer to source).
#' @param model2 Second AlignmentModel (applied second, closer to target).
#'
#' @return A new AlignmentModel with composed transforms.
#'
#' @details
#' The composition follows standard matrix composition order:
#' if \code{model1} transforms from A to B, and \code{model2} transforms
#' from B to C, then the composed model transforms from A to C.
#'
#' The composed transform for subject s is: T_composed = T2_s \%*\% T1_s
#'
#' @examples
#' \dontrun{
#' # Geometric alignment (native -> MNI)
#' geo_model <- fit_alignment(data, method = "geometric")
#'
#' # Functional alignment (MNI -> reference)
#' func_model <- fit_alignment(geo_aligned, method = "procrustes")
#'
#' # Compose into single transform
#' composed <- compose_alignment(geo_model, func_model)
#'
#' # Or use %*%
#' composed <- func_model %*% geo_model
#' }
#'
#' @export
compose_alignment <- function(model1, model2) {
  # Extract models from results if needed
  if (inherits(model1, "AlignmentResult")) {
    model1 <- get_model(model1)
  }
  if (inherits(model2, "AlignmentResult")) {
    model2 <- get_model(model2)
  }

  if (!inherits(model1, "AlignmentModel")) {
    stop("'model1' must be an AlignmentModel")
  }
  if (!inherits(model2, "AlignmentModel")) {
    stop("'model2' must be an AlignmentModel")
  }

  caps1 <- aligner_capabilities(model1@method)
  if (!is.null(caps1) && !identical(caps1$returns %||% "operator", "operator")) {
    stop(
      sprintf(
        "Method '%s' does not return operator transforms; compose_alignment() currently supports operators only",
        model1@method
      ),
      call. = FALSE
    )
  }
  caps2 <- aligner_capabilities(model2@method)
  if (!is.null(caps2) && !identical(caps2$returns %||% "operator", "operator")) {
    stop(
      sprintf(
        "Method '%s' does not return operator transforms; compose_alignment() currently supports operators only",
        model2@method
      ),
      call. = FALSE
    )
  }

  # Get common subjects
  subjects1 <- names(model1@transforms)
  subjects2 <- names(model2@transforms)
  common_subjects <- intersect(subjects1, subjects2)

  if (length(common_subjects) == 0) {
    stop("Models have no subjects in common")
  }

  if (length(common_subjects) < length(subjects1) ||
      length(common_subjects) < length(subjects2)) {
    warning(sprintf(
      "Only %d subjects in common; %d/%d from model1, %d/%d from model2",
      length(common_subjects),
      length(common_subjects), length(subjects1),
      length(common_subjects), length(subjects2)
    ))
  }

  # Compose transforms: T_composed = T2 %*% T1
  composed_transforms <- lapply(common_subjects, function(subj) {
    t1 <- model1@transforms[[subj]]
    t2 <- model2@transforms[[subj]]

    if (!.is_matrixish(t1) || !.is_matrixish(t2)) {
      stop(
        sprintf("Non-matrix transforms cannot be composed (subject '%s')", subj),
        call. = FALSE
      )
    }

    # Check dimension compatibility
    if (ncol(t2) != nrow(t1)) {
      stop(sprintf(
        "Dimension mismatch for subject '%s': model2 expects %d, model1 provides %d",
        subj, ncol(t2), nrow(t1)
      ))
    }

    t2 %*% t1
  })
  names(composed_transforms) <- common_subjects

  # Build provenance for composed model
  provenance <- list(
    composed_from = list(
      model1 = list(
        method = model1@method,
        provenance = model1@provenance
      ),
      model2 = list(
        method = model2@method,
        provenance = model2@provenance
      )
    ),
    composed_at = Sys.time(),
    neuralign_version = as.character(utils::packageVersion("neuralign"))
  )

  AlignmentModel(
    transforms = composed_transforms,
    reference = model2@reference,
    reference_data = model2@reference_data,
    method = sprintf("%s+%s", model1@method, model2@method),
    space_from = model1@space_from,
    space_to = model2@space_to,
    params = list(),
    method_state = list(
      model1_state = model1@method_state,
      model2_state = model2@method_state
    ),
    train_subjects = common_subjects
  )
}


#' Matrix Multiplication for AlignmentModels
#'
#' Compose two alignment models using matrix multiplication syntax.
#' model2 \%*\% model1 composes model1 followed by model2.
#'
#' @param x Left operand (applied second).
#' @param y Right operand (applied first).
#'
#' @return Composed AlignmentModel.
#'
#' @export
setMethod("%*%", c("AlignmentModel", "AlignmentModel"),
  function(x, y) {
    # x is applied after y, so x = model2, y = model1
    compose_alignment(y, x)
  }
)


#' Compose Transform with Data
#'
#' Apply an alignment model to data using matrix multiplication syntax.
#'
#' @param x An AlignmentModel.
#' @param y Data (matrix or AlignmentData).
#'
#' @return Transformed data.
#'
#' @export
setMethod("%*%", c("AlignmentModel", "matrix"),
  function(x, y) {
    # Assume single subject, apply first transform
    if (length(x@transforms) == 0) {
      stop("Model has no transforms")
    }

    # Use first transform
    transform <- x@transforms[[1]]
    apply_transform(transform, y)
  }
)


#' Check Composition Compatibility
#'
#' Check if two models can be composed.
#'
#' @param model1 First AlignmentModel.
#' @param model2 Second AlignmentModel.
#'
#' @return List with \code{compatible} logical and \code{message} string.
#'
#' @export
check_composition <- function(model1, model2) {
  # Extract from results if needed
  if (inherits(model1, "AlignmentResult")) {
    model1 <- get_model(model1)
  }
  if (inherits(model2, "AlignmentResult")) {
    model2 <- get_model(model2)
  }

  # Check common subjects
  common <- intersect(names(model1@transforms), names(model2@transforms))
  if (length(common) == 0) {
    return(list(
      compatible = FALSE,
      message = "No subjects in common between models"
    ))
  }

  # Check dimensions for first common subject
  subj <- common[1]
  t1 <- model1@transforms[[subj]]
  t2 <- model2@transforms[[subj]]

  if (!.is_matrixish(t1) || !.is_matrixish(t2)) {
    return(list(
      compatible = FALSE,
      message = "Non-matrix transforms cannot be composed"
    ))
  }

  if (ncol(t2) != nrow(t1)) {
    return(list(
      compatible = FALSE,
      message = sprintf(
        "Dimension mismatch: model1 output (%d) != model2 input (%d)",
        nrow(t1), ncol(t2)
      )
    ))
  }

  list(
    compatible = TRUE,
    message = sprintf(
      "Models are compatible (%d common subjects, dims: %d -> %d -> %d)",
      length(common), ncol(t1), nrow(t1), nrow(t2)
    )
  )
}
