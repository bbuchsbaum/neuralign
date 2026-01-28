#' AlignmentResult Class
#'
#' Represents the result of applying an alignment model to data.
#' Contains the aligned data along with quality metrics and CV info.
#'
#' @family alignment_result
#'
#' @slot model The AlignmentModel used to produce this result.
#' @slot aligned Named list of aligned data (lazy or realized).
#' @slot quality List of quality diagnostics.
#' @slot cv_info Cross-validation information (fold assignments, etc.).
#'
#' @export
setClass("AlignmentResult",
  slots = c(
    model = "AlignmentModel",
    aligned = "list",
    quality = "list",
    cv_info = "list"
  ),
  prototype = list(
    aligned = list(),
    quality = list(),
    cv_info = list()
  )
)

setValidity("AlignmentResult", function(object) {
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
  if (length(errors) == 0L) TRUE else errors
})


#' Create an AlignmentResult Object
#'
#' @param model The AlignmentModel used.
#' @param aligned Named list of aligned data.
#' @param quality Optional list of quality metrics.
#' @param cv_info Optional cross-validation information.
#'
#' @return An AlignmentResult object.
#'
#' @family alignment_result
#'
#' @export
AlignmentResult <- function(model,
                            aligned,
                            quality = list(),
                            cv_info = list()) {
  if (!inherits(model, "AlignmentModel")) {
    stop("'model' must be an AlignmentModel object", call. = FALSE)
  }

  if (!is.list(aligned)) {
    stop("'aligned' must be a list", call. = FALSE)
  }

  new("AlignmentResult",
    model = model,
    aligned = aligned,
    quality = quality,
    cv_info = cv_info
  )
}


#' Print Method for AlignmentResult
#'
#' @param object An AlignmentResult object.
#'
#' @export
setMethod("show", "AlignmentResult", function(object) {
  cat("AlignmentResult object\n")
  cat(sprintf("  Method: %s\n", object@model@method))
  cat(sprintf("  Aligned subjects: %d\n", length(object@aligned)))

  if (length(object@aligned) > 0) {
    # Show data dimensions
    first_data <- object@aligned[[1]]
    if (is.matrix(first_data) || inherits(first_data, "Matrix")) {
      cat(sprintf("  Aligned dims: %d x %d\n",
        nrow(first_data), ncol(first_data)
      ))
    }
  }

  # Quality metrics
  if (length(object@quality) > 0) {
    cat("  Quality metrics:\n")
    for (nm in names(object@quality)) {
      val <- object@quality[[nm]]
      if (is.numeric(val) && length(val) == 1) {
        cat(sprintf("    %s: %.4f\n", nm, val))
      } else if (is.numeric(val)) {
        cat(sprintf("    %s: mean=%.4f, sd=%.4f\n", nm, mean(val), sd(val)))
      }
    }
  }

  # CV info
  if (length(object@cv_info) > 0) {
    cat("  Cross-validation:\n")
    if (!is.null(object@cv_info$method)) {
      cat(sprintf("    Method: %s\n", object@cv_info$method))
    }
    if (!is.null(object@cv_info$axis)) {
      cat(sprintf("    Axis: %s\n", object@cv_info$axis))
    }
    if (!is.null(object@cv_info$n_folds)) {
      cat(sprintf("    Folds: %d\n", object@cv_info$n_folds))
    }
    if (!is.null(object@cv_info$anchor_common)) {
      cat(sprintf(
        "    Common anchor: %s\n",
        if (isTRUE(object@cv_info$anchor_common)) "yes" else "no"
      ))
    }
  }
})

.ensure_result <- function(x, what = "result") {
  if (!inherits(x, "AlignmentResult")) {
    stop(sprintf("'%s' must be an AlignmentResult object", what), call. = FALSE)
  }
  x
}


#' Get aligned data
#'
#' Extract aligned data from an [AlignmentResult] object.
#'
#' @param result An [AlignmentResult] object.
#' @param subject Optional subject ID. If `NULL` (default), returns the full named
#'   list of aligned data for all subjects.
#'
#' @return A numeric matrix (when `subject` is provided) or a named list of
#'   numeric matrices (when `subject` is `NULL`).
#'
#' @family alignment_result
#'
#' @examples
#' set.seed(1)
#' adat <- AlignmentData(list(
#'   s1 = matrix(rnorm(20), 4, 5),
#'   s2 = matrix(rnorm(20), 4, 5)
#' ))
#' res <- fit_alignment(adat, method = "procrustes", reference = "s1", compute_quality = FALSE)
#' all_aligned <- get_aligned(res)
#' s1_aligned <- get_aligned(res, "s1")
#' dim(s1_aligned)
#'
#' @export
get_aligned <- function(result, subject = NULL) {
  result <- .ensure_result(result, what = "result")

  if (is.null(subject)) {
    return(result@aligned)
  }

  if (!subject %in% names(result@aligned)) {
    stop(sprintf("Subject '%s' not found in result", subject), call. = FALSE)
  }
  result@aligned[[subject]]
}


#' @rdname get_aligned
#' @export
aligned_data <- function(result, subject = NULL) {
  get_aligned(result, subject = subject)
}


#' Get Quality Metrics
#'
#' @param result An AlignmentResult object.
#' @param metric Optional specific metric name.
#'
#' @return Quality metrics list or specific metric value.
#'
#' @family alignment_result
#'
#' @examples
#' set.seed(1)
#' adat <- AlignmentData(list(
#'   s1 = matrix(rnorm(20), 4, 5),
#'   s2 = matrix(rnorm(20), 4, 5)
#' ))
#' res <- fit_alignment(adat, method = "procrustes", reference = "s1", compute_quality = TRUE)
#' qm <- get_quality(res)
#' names(qm)
#'
#' @export
get_quality <- function(result, metric = NULL) {
  result <- .ensure_result(result, what = "result")

  if (is.null(metric)) {
    return(result@quality)
  }

  if (!metric %in% names(result@quality)) {
    stop(sprintf("Metric '%s' not found", metric), call. = FALSE)
  }
  result@quality[[metric]]
}


#' Get Model from Result
#'
#' @param result An AlignmentResult object.
#'
#' @return The AlignmentModel.
#'
#' @family alignment_result
#'
#' @examples
#' set.seed(1)
#' adat <- AlignmentData(list(
#'   s1 = matrix(rnorm(20), 4, 5),
#'   s2 = matrix(rnorm(20), 4, 5)
#' ))
#' res <- fit_alignment(adat, method = "procrustes", reference = "s1", compute_quality = FALSE)
#' mdl <- get_model(res)
#' inherits(mdl, "AlignmentModel")
#'
#' @export
get_model <- function(result) {
  result <- .ensure_result(result, what = "result")
  result@model
}


#' Get CV Info from Result
#'
#' @param result An AlignmentResult object.
#'
#' @return A list containing cross-validation metadata. Common fields include:
#' \describe{
#'   \item{method}{CV method identifier (e.g., `"none"`, `"loso"`, `"kfold"`, `"run"`, `"blocked_time"`, or `"custom"`).}
#'   \item{axis}{CV axis (typically `"subject"` or `"observation"`).}
#'   \item{n_folds}{Number of folds.}
#'   \item{folds}{Fold specification used during fitting.}
#'   \item{fold_assignments}{Subject-to-fold assignments (subject-axis CV only; may be NULL).}
#'   \item{fold_ids}{Fold identifiers (observation-axis CV only; may be NULL).}
#'   \item{reference_kind}{One of `"fixed_subject"`, `"template"`, or `"data_driven"`.}
#'   \item{anchor_common}{Logical; whether folds share a common anchor/reference space.}
#'   \item{anchor_by_subject}{Per-subject reference used during evaluation (subject-axis CV only; may be NULL).}
#'   \item{reference_by_fold}{Per-fold reference used during evaluation (observation-axis CV only; may be NULL).}
#'   \item{transforms_by_fold}{Optional per-fold transforms (observation-axis CV for operator-returning methods only; may be NULL).}
#'   \item{anchor_note}{Optional warning/note about fold-specific anchors.}
#' }
#'
#' Additional fields may be present depending on the CV strategy and aligner.
#'
#' @family alignment_result
#'
#' @examples
#' set.seed(1)
#' adat <- AlignmentData(list(
#'   s1 = matrix(rnorm(20), 4, 5),
#'   s2 = matrix(rnorm(20), 4, 5),
#'   s3 = matrix(rnorm(20), 4, 5)
#' ))
#' res <- fit_alignment(adat, method = "procrustes", reference = "s1", cv = "loso",
#'   compute_quality = FALSE, return_aligned = FALSE
#' )
#' info <- get_cv_info(res)
#' info$method
#'
#' @export
get_cv_info <- function(result) {
  result <- .ensure_result(result, what = "result")
  result@cv_info
}


#' Convert Aligned Data to Matrix
#'
#' Stack all aligned subjects into a single matrix.
#'
#' @param result An AlignmentResult object.
#' @param by How to stack: "subject" (list element per subject) or
#'   "observation" (concatenate columns).
#'
#' @return Matrix or list of matrices.
#'
#' @examples
#' set.seed(1)
#' adat <- AlignmentData(list(
#'   s1 = matrix(rnorm(20), 4, 5),
#'   s2 = matrix(rnorm(20), 4, 5)
#' ))
#' res <- fit_alignment(adat, method = "procrustes", reference = "s1", compute_quality = FALSE)
#' mat <- as_aligned_matrix(res, by = "observation")
#' dim(mat)
#'
#' @export
as_aligned_matrix <- function(result, by = c("subject", "observation")) {
  result <- .ensure_result(result, what = "result")
  by <- match.arg(by)

  if (by == "subject") {
    return(result@aligned)
  }

  # Stack by observation (concatenate columns)
  do.call(cbind, result@aligned)
}


#' Subset AlignmentResult by Subject
#'
#' @name sub-AlignmentResult
#' @aliases [,AlignmentResult,ANY,ANY,ANY-method
#' @param x An AlignmentResult object.
#' @param i Subject IDs or indices.
#'
#' @return A new AlignmentResult with only selected subjects.
#'
#' @export
setMethod("[", c("AlignmentResult", "ANY"),
  function(x, i) {
    subj_names <- names(x@aligned)
    i <- .resolve_subject_subset(i, subj_names, what = "subjects")

    AlignmentResult(
      model = x@model[i],
      aligned = x@aligned[i],
      quality = x@quality,
      cv_info = x@cv_info
    )
  }
)


#' Get Number of Aligned Subjects
#'
#' @param x An AlignmentResult object.
#'
#' @return Integer number of aligned subjects.
#'
#' @export
setMethod("length", "AlignmentResult", function(x) length(x@aligned))
