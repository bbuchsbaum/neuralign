#' Compose Alignment Models
#'
#' Compose two alignment models into a single model that applies both
#' transformations. Useful for combining geometric (spatial) and functional
#' alignments into a single pipeline.
#'
#' @param model1 First AlignmentModel (applied first, closer to source).
#' @param model2 Second AlignmentModel (applied second, closer to target).
#' @param allow_unverified_spaces Logical; if `TRUE`, compose models whose
#'   intermediate spaces are unknown (`NULL`). The result is stamped unverified.
#' @param allow_partial Logical; if `TRUE`, drop subjects that are not in both
#'   models. Equivalent to passing `subjects = intersect(...)`.
#' @param allow_space_mismatch Logical; expert override that composes models
#'   with known unequal intermediate spaces and stamps the result unverified.
#' @param subjects Optional character vector of subjects to compose. Must be
#'   present in both models.
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
compose_alignment <- function(model1,
                              model2,
                              allow_unverified_spaces = FALSE,
                              allow_partial = FALSE,
                              allow_space_mismatch = FALSE,
                              subjects = NULL) {
  model1 <- .ensure_model(model1, what = "model1")
  model2 <- .ensure_model(model2, what = "model2")

  if (.model_is_fold_specific(model1) || .model_is_fold_specific(model2)) {
    stop(
      "compose_alignment() is not available for fold-specific models: there is no common space/anchor",
      call. = FALSE
    )
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

  space_stamp <- .assert_compose_spaces(
    model1,
    model2,
    allow_unverified_spaces = isTRUE(allow_unverified_spaces),
    allow_space_mismatch = isTRUE(allow_space_mismatch)
  )

  subjects1 <- names(model1@transforms)
  subjects2 <- names(model2@transforms)
  common_subjects <- intersect(subjects1, subjects2)

  if (length(common_subjects) == 0) {
    stop("Models have no subjects in common", call. = FALSE)
  }

  if (!is.null(subjects)) {
    if (!is.character(subjects) || !length(subjects) || any(!nzchar(subjects))) {
      stop("'subjects' must be a non-empty character vector", call. = FALSE)
    }
    missing <- setdiff(subjects, common_subjects)
    if (length(missing)) {
      stop(
        sprintf(
          "Requested subjects are not in both models: %s",
          paste(missing, collapse = ", ")
        ),
        call. = FALSE
      )
    }
    common_subjects <- unique(subjects)
    allow_partial <- TRUE
  } else if (length(common_subjects) < length(subjects1) ||
             length(common_subjects) < length(subjects2)) {
    if (!isTRUE(allow_partial)) {
      stop(
        sprintf(
          "Partial subject drop requires allow_partial=TRUE or subjects=intersect(...); %d in common, %d/%d from model1, %d/%d from model2",
          length(common_subjects),
          length(common_subjects), length(subjects1),
          length(common_subjects), length(subjects2)
        ),
        call. = FALSE
      )
    }
  }

  # Compose transforms: T_composed = T2 %*% T1
  composed_transforms <- lapply(common_subjects, function(subj) {
    t1 <- model1@transforms[[subj]]
    t2 <- model2@transforms[[subj]]

    if (!.transform_is_operator(t1) || !.transform_is_operator(t2)) {
      stop(
        sprintf("Non-operator transforms cannot be composed (subject '%s')", subj),
        call. = FALSE
      )
    }

    .compose_operator_transforms(
      t2,
      t1,
      context = sprintf("compose_alignment: subject '%s'", subj)
    )
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
    neuralign_version = as.character(utils::packageVersion("neuralign")),
    space_verification = space_stamp,
    partial_subjects = isTRUE(allow_partial) && (
      length(common_subjects) < length(subjects1) ||
        length(common_subjects) < length(subjects2)
    )
  )

  AlignmentModel(
    transforms = composed_transforms,
    reference = model2@reference,
    reference_data = model2@reference_data,
    method = sprintf("%s+%s", model1@method, model2@method),
    space_from = model1@space_from,
    space_to = model2@space_to,
    method_state = list(
      model1_state = model1@method_state,
      model2_state = model2@method_state
    ),
    train_subjects = common_subjects,
    provenance = provenance
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
    n_transforms <- length(x@transforms)
    if (n_transforms == 0) {
      stop("Model has no transforms", call. = FALSE)
    }

    if (n_transforms > 1) {
      stop(
        paste0(
          "AlignmentModel has transforms for multiple subjects; cannot apply ",
          "to a bare matrix. Use apply_alignment(model, AlignmentData) or ",
          "apply_transform(get_transforms(model)[[subject_id]], mat)."
        ),
        call. = FALSE
      )
    }

    transform <- x@transforms[[1L]]
    apply_transform(transform, y)
  }
)


.space_is_unknown <- function(space) {
  is.null(space)
}

.assert_compose_spaces <- function(model1,
                                   model2,
                                   allow_unverified_spaces = FALSE,
                                   allow_space_mismatch = FALSE) {
  from <- model2@space_from
  to <- model1@space_to
  unknown <- .space_is_unknown(to) || .space_is_unknown(from)
  stamp <- list(
    unverified_spaces = FALSE,
    space_mismatch_allowed = FALSE,
    model1_space_to = to,
    model2_space_from = from
  )

  if (unknown) {
    if (!isTRUE(allow_unverified_spaces)) {
      stop(
        "Cannot compose models with unknown/unverified intermediate spaces; set allow_unverified_spaces=TRUE to override",
        call. = FALSE
      )
    }
    stamp$unverified_spaces <- TRUE
    return(stamp)
  }

  if (!spaces_compatible(to, from)) {
    if (!isTRUE(allow_space_mismatch)) {
      stop(
        sprintf(
          "Space chain mismatch: model1 maps to '%s' but model2 expects '%s'",
          .format_space(to), .format_space(from)
        ),
        call. = FALSE
      )
    }
    stamp$unverified_spaces <- TRUE
    stamp$space_mismatch_allowed <- TRUE
  }
  stamp
}


#' Check Composition Compatibility
#'
#' Check if two models can be composed.
#'
#' @param model1 First AlignmentModel.
#' @param model2 Second AlignmentModel.
#' @param allow_unverified_spaces,allow_partial,allow_space_mismatch Same
#'   meaning as in [compose_alignment()].
#'
#' @return List with \code{compatible} logical and \code{message} string.
#'
#' @examples
#' Q <- diag(3)
#' m1 <- AlignmentModel(list(s1 = Q), reference = "s1", method = "procrustes")
#' m2 <- AlignmentModel(list(s1 = Q), reference = "s1", method = "procrustes")
#' check_composition(m1, m2, allow_unverified_spaces = TRUE)$compatible
#'
#' @export
check_composition <- function(model1,
                              model2,
                              allow_unverified_spaces = FALSE,
                              allow_partial = FALSE,
                              allow_space_mismatch = FALSE) {
  model1 <- .ensure_model(model1, what = "model1")
  model2 <- .ensure_model(model2, what = "model2")

  if (.model_is_fold_specific(model1) || .model_is_fold_specific(model2)) {
    return(list(
      compatible = FALSE,
      message = "fold-specific models have no common space/anchor"
    ))
  }

  space_ok <- tryCatch(
    {
      .assert_compose_spaces(
        model1,
        model2,
        allow_unverified_spaces = isTRUE(allow_unverified_spaces),
        allow_space_mismatch = isTRUE(allow_space_mismatch)
      )
      TRUE
    },
    error = function(e) conditionMessage(e)
  )
  if (!isTRUE(space_ok)) {
    return(list(compatible = FALSE, message = space_ok))
  }

  # Check common subjects
  common <- intersect(names(model1@transforms), names(model2@transforms))
  if (length(common) == 0) {
    return(list(
      compatible = FALSE,
      message = "No subjects in common between models"
    ))
  }

  n1 <- length(model1@transforms)
  n2 <- length(model2@transforms)
  if ((length(common) < n1 || length(common) < n2) && !isTRUE(allow_partial)) {
    return(list(
      compatible = FALSE,
      message = "Partial subject drop requires allow_partial=TRUE or subjects=intersect(...)"
    ))
  }

  # Check dimensions for all common subjects
  for (subj in common) {
    t1 <- model1@transforms[[subj]]
    t2 <- model2@transforms[[subj]]

    if (!.transform_is_operator(t1) || !.transform_is_operator(t2)) {
      return(list(
        compatible = FALSE,
        message = sprintf(
          "Non-operator transforms cannot be composed (subject '%s')", subj
        )
      ))
    }

    dims1 <- .transform_dims(t1)
    dims2 <- .transform_dims(t2)
    if (!identical(dims2[["source"]], dims1[["target"]])) {
      return(list(
        compatible = FALSE,
        message = sprintf(
          "Dimension mismatch for subject '%s': model1 output (%d) != model2 input (%d)",
          subj, dims1[["target"]], dims2[["source"]]
        )
      ))
    }
  }

  # Use first subject for summary dims
  dims1 <- .transform_dims(model1@transforms[[common[1]]])
  dims2 <- .transform_dims(model2@transforms[[common[1]]])

  list(
    compatible = TRUE,
    message = sprintf(
      "Models are compatible (%d common subjects, dims: %d -> %d -> %d)",
      length(common), dims1[["source"]], dims1[["target"]], dims2[["target"]]
    )
  )
}
