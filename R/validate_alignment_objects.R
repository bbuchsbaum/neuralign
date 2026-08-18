#' Alignment Object Validators
#'
#' Authoritative checks for [AlignmentModel] and [AlignmentResult]. Constructors,
#' S4 validity, `add_transform()`, subsetting, compose, load, and apply all use
#' these helpers so one model means one identified target coordinate system.
#'
#' @name validate_alignment_objects
#' @keywords internal
NULL


.transform_is_finite <- function(transform) {
  if (.is_embedding_transform(transform)) {
    z <- transform$aligned
    return(is.numeric(z) && all(is.finite(z)))
  }
  if (.is_low_rank_transform(transform)) {
    return(is.numeric(transform$U) && is.numeric(transform$V) &&
             all(is.finite(transform$U)) && all(is.finite(transform$V)))
  }
  if (.is_matrixish(transform)) {
    x <- as.matrix(transform)
    return(is.numeric(x) && all(is.finite(x)))
  }
  FALSE
}


validate_alignment_transforms <- function(transforms,
                                          context = "AlignmentModel") {
  errors <- character()
  if (!is.list(transforms)) {
    return(sprintf("%s: 'transforms' must be a list", context))
  }
  if (!length(transforms)) {
    return(character())
  }

  nms <- names(transforms)
  if (is.null(nms) || anyNA(nms) || any(!nzchar(nms))) {
    errors <- c(errors, sprintf(
      "%s: transform names must be non-empty subject IDs", context
    ))
  }
  if (!is.null(nms) && anyDuplicated(nms) > 0L) {
    errors <- c(errors, sprintf(
      "%s: transform names must be unique (duplicate subject IDs)", context
    ))
  }

  target_dims <- integer()
  for (i in seq_along(transforms)) {
    tr <- transforms[[i]]
    lab <- if (!is.null(nms) && i <= length(nms) && !is.na(nms[[i]]) &&
               nzchar(nms[[i]])) {
      nms[[i]]
    } else {
      as.character(i)
    }

    if (!.transform_is_operator(tr) && !.is_embedding_transform(tr)) {
      errors <- c(errors, sprintf(
        "%s: transform '%s' is not a supported operator", context, lab
      ))
      next
    }
    if (!isTRUE(.transform_is_finite(tr))) {
      errors <- c(errors, sprintf(
        "%s: transform '%s' must be finite (no NA/Inf)", context, lab
      ))
    }
    td <- .transform_target_dim(tr)
    if (!is.na(td)) {
      target_dims <- c(target_dims, as.integer(td))
    }
  }

  if (length(unique(target_dims)) > 1L) {
    errors <- c(errors, sprintf(
      "%s: transforms must share one target/codomain dimension", context
    ))
  }
  errors
}


.validate_named_transforms <- function(transforms, context = "AlignmentModel") {
  if (!is.list(transforms)) {
    stop(sprintf("%s: 'transforms' must be a named list of operators", context),
         call. = FALSE)
  }
  if (!length(transforms)) {
    return(invisible(TRUE))
  }
  nms <- names(transforms)
  if (is.null(nms)) {
    stop(sprintf("%s: 'transforms' must be named with subject IDs", context),
         call. = FALSE)
  }
  errors <- validate_alignment_transforms(transforms, context = context)
  if (length(errors)) {
    stop(paste(errors, collapse = "; "), call. = FALSE)
  }
  invisible(TRUE)
}


.validate_alignment_model_object <- function(object) {
  errors <- character()
  if (!is.list(object@transforms)) {
    errors <- c(errors, "'transforms' must be a list")
  }
  if (length(object@transforms) > 0 && is.null(names(object@transforms))) {
    errors <- c(errors, "'transforms' must be named with subject IDs")
  }
  if (!is.character(object@method) || length(object@method) != 1L ||
      is.na(object@method) || !nzchar(object@method)) {
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
  more <- validate_alignment_transforms(object@transforms, context = "AlignmentModel")
  unique(c(errors, more))
}


.validate_alignment_result_object <- function(object) {
  errors <- character()
  if (!inherits(object@model, "AlignmentModel")) {
    errors <- c(errors, "'model' must be an AlignmentModel object")
  }
  if (!is.list(object@aligned)) {
    errors <- c(errors, "'aligned' must be a list")
  }
  if (!is.list(object@quality)) {
    errors <- c(errors, "'quality' must be a list")
  }
  if (!is.list(object@cv_info)) {
    errors <- c(errors, "'cv_info' must be a list")
  }

  model_ok <- tryCatch(
    .validate_alignment_model_object(object@model),
    error = function(e) conditionMessage(e)
  )
  if (!is.character(model_ok) && !isTRUE(model_ok)) {
    errors <- c(errors, as.character(model_ok))
  } else if (is.character(model_ok) && length(model_ok)) {
    errors <- c(errors, model_ok)
  }

  aligned <- object@aligned
  fold_specific <- .cv_has_fold_specific_spaces(object@cv_info) ||
    (inherits(object@model, "AlignmentModel") &&
       identical(object@model@reference, "fold_specific"))
  if (isTRUE(fold_specific) && length(aligned) > 0L) {
    errors <- c(
      errors,
      "fold-specific AlignmentResult cannot store a global aligned matrix"
    )
  }
  if (isTRUE(fold_specific) && length(object@quality) > 0L) {
    errors <- c(
      errors,
      "fold-specific AlignmentResult cannot store cross-fold quality"
    )
  }

  if (length(aligned) > 0L) {
    nms <- names(aligned)
    if (is.null(nms) || anyNA(nms) || any(!nzchar(nms)) ||
        anyDuplicated(nms) > 0L) {
      errors <- c(errors, "aligned names must be unique, non-empty subject IDs")
    } else if (inherits(object@model, "AlignmentModel")) {
      model_subjects <- names(object@model@transforms)
      extra <- setdiff(nms, model_subjects)
      if (length(extra) && length(model_subjects)) {
        errors <- c(errors, "aligned subjects must agree with the attached model")
      }
      target <- NA_integer_
      if (length(object@model@transforms)) {
        target <- .transform_target_dim(object@model@transforms[[1L]])
      }
      for (subj in nms) {
        mat <- aligned[[subj]]
        if (!.is_matrixish(mat)) {
          errors <- c(errors, sprintf("aligned['%s'] must be matrix-like", subj))
          next
        }
        if (!is.na(target) && nrow(as.matrix(mat)) != target) {
          errors <- c(errors, sprintf(
            "aligned['%s'] codomain does not match the model target dimension",
            subj
          ))
        }
      }
    }
  }
  unique(errors)
}


validate_alignment_model <- function(object, action = c("check", "error")) {
  action <- match.arg(action)
  if (!inherits(object, "AlignmentModel")) {
    errors <- "'object' must be an AlignmentModel"
  } else {
    errors <- .validate_alignment_model_object(object)
  }
  if (identical(action, "check")) {
    if (length(errors) == 0L) TRUE else errors
  } else if (length(errors)) {
    stop(paste(errors, collapse = "; "), call. = FALSE)
  } else {
    invisible(TRUE)
  }
}


validate_alignment_result <- function(object, action = c("check", "error")) {
  action <- match.arg(action)
  if (!inherits(object, "AlignmentResult")) {
    errors <- "'object' must be an AlignmentResult"
  } else {
    errors <- .validate_alignment_result_object(object)
  }
  if (identical(action, "check")) {
    if (length(errors) == 0L) TRUE else errors
  } else if (length(errors)) {
    stop(paste(errors, collapse = "; "), call. = FALSE)
  } else {
    invisible(TRUE)
  }
}


.model_is_fold_specific <- function(model) {
  inherits(model, "AlignmentModel") && identical(model@reference, "fold_specific")
}


.abort_if_fold_specific_aligned <- function(result, what) {
  if (inherits(result, "AlignedResampleSet")) {
    stop(
      sprintf(
        "%s cannot use AlignedResampleSet: fold-specific spaces have no common aligned matrix",
        what
      ),
      call. = FALSE
    )
  }
  if (!inherits(result, "AlignmentResult")) {
    return(invisible(TRUE))
  }
  if (.cv_has_fold_specific_spaces(result@cv_info) ||
      .model_is_fold_specific(result@model)) {
    stop(
      sprintf(
        "%s cannot use a fold-specific result; there is no common aligned space/anchor. Use as_aligned_resample_set().",
        what
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}


.cv_has_fold_specific_spaces <- function(cv) {
  if (!is.list(cv) || !length(cv)) return(FALSE)
  method <- cv$method %||% "none"
  if (identical(method, "none") || identical(method, "applied") ||
      is.null(method)) {
    return(FALSE)
  }
  if (!is.null(cv$anchor_common)) {
    return(!isTRUE(cv$anchor_common))
  }
  TRUE
}
