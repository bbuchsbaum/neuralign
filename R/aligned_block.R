#' Aligned Observation Block
#'
#' One subject/session block of analysis-facing aligned data
#' (observations × shared features), with observation metadata and role.
#'
#' Note: this is distinct from [alignment_feature_block()], which describes
#' *pre-alignment* correspondence signals.
#'
#' @param values Numeric matrix: observations × shared features.
#' @param observation_data `data.frame` with one row per observation.
#' @param subject_id Subject identifier.
#' @param session_id Optional session identifier.
#' @param representation Optional client-defined representation label.
#' @param uncertainty Optional uncertainty provider (list or `NULL`).
#' @param source_ref Optional reference to source features / provenance ids.
#' @param alignment_role One of `"fit"`, `"calibration"`, `"analysis"`,
#'   `"assessment"`, `"test"`.
#' @param shared_space_id Optional shared-space id for this block.
#' @param metadata Named list of additional fields.
#'
#' @return An `AlignedBlock` object.
#'
#' @family aligned_study
#' @export
AlignedBlock <- function(values,
                         observation_data = NULL,
                         subject_id,
                         session_id = NULL,
                         representation = NA_character_,
                         uncertainty = NULL,
                         source_ref = NULL,
                         alignment_role = c(
                           "analysis", "fit", "calibration",
                           "assessment", "test"
                         ),
                         shared_space_id = NULL,
                         metadata = list()) {
  alignment_role <- match.arg(alignment_role)

  if (missing(subject_id) || !is.character(subject_id) || length(subject_id) != 1L ||
      is.na(subject_id) || !nzchar(subject_id)) {
    stop("'subject_id' must be a non-empty character string", call. = FALSE)
  }
  if (!.is_matrixish(values)) {
    stop("'values' must be a matrix or Matrix", call. = FALSE)
  }
  numeric_values <- is.numeric(values) ||
    methods::is(values, "dMatrix") || methods::is(values, "iMatrix")
  if (!numeric_values) {
    stop("'values' must be numeric", call. = FALSE)
  }
  if (any(!is.finite(values))) {
    stop("'values' must contain only finite values", call. = FALSE)
  }
  n_obs <- nrow(values)
  n_feat <- ncol(values)

  if (is.null(observation_data)) {
    observation_data <- data.frame(
      observation_id = seq_len(n_obs),
      stringsAsFactors = FALSE
    )
  }
  if (!is.data.frame(observation_data)) {
    stop("'observation_data' must be a data.frame", call. = FALSE)
  }
  if (nrow(observation_data) != n_obs) {
    stop(sprintf(
      "'observation_data' must have %d rows (one per observation), got %d",
      n_obs, nrow(observation_data)
    ), call. = FALSE)
  }
  if (!is.list(metadata)) {
    stop("'metadata' must be a list", call. = FALSE)
  }
  if (!is.null(session_id) &&
      (!is.character(session_id) || length(session_id) != 1L ||
       is.na(session_id) || !nzchar(session_id))) {
    stop("'session_id' must be NULL or a non-empty character string",
         call. = FALSE)
  }
  if (!is.character(representation) || length(representation) != 1L) {
    stop("'representation' must be a single character string", call. = FALSE)
  }
  if (!is.null(shared_space_id) &&
      (!is.character(shared_space_id) || length(shared_space_id) != 1L ||
       is.na(shared_space_id) || !nzchar(shared_space_id))) {
    stop("'shared_space_id' must be NULL or a non-empty character string",
         call. = FALSE)
  }

  structure(
    list(
      values = values,
      observation_data = observation_data,
      subject_id = subject_id,
      session_id = session_id,
      representation = as.character(representation)[[1L]],
      uncertainty = uncertainty,
      source_ref = source_ref,
      alignment_role = alignment_role,
      shared_space_id = shared_space_id,
      metadata = metadata,
      n_obs = n_obs,
      n_features = n_feat
    ),
    class = c("AlignedBlock", "list")
  )
}


#' @export
print.AlignedBlock <- function(x, ...) {
  cat("AlignedBlock\n")
  cat(sprintf("  subject_id: %s\n", x$subject_id))
  if (!is.null(x$session_id)) {
    cat(sprintf("  session_id: %s\n", x$session_id))
  }
  cat(sprintf("  dims: %d obs x %d shared features\n", x$n_obs, x$n_features))
  cat(sprintf("  representation: %s\n", x$representation %||% NA_character_))
  cat(sprintf("  alignment_role: %s\n", x$alignment_role))
  if (!is.null(x$shared_space_id)) {
    cat(sprintf("  shared_space_id: %s\n", x$shared_space_id))
  }
  invisible(x)
}


.validate_aligned_block <- function(block, shared_space = NULL, name = "block") {
  if (!inherits(block, "AlignedBlock")) {
    stop(sprintf("'%s' must be an AlignedBlock", name), call. = FALSE)
  }
  if (!is.null(shared_space)) {
    sid <- shared_space$id
    if (!is.null(block$shared_space_id) && !identical(block$shared_space_id, sid)) {
      stop(sprintf(
        "'%s' shared_space_id '%s' does not match study space '%s'",
        name, block$shared_space_id, sid
      ), call. = FALSE)
    }
    if (ncol(block$values) != shared_space$dimension) {
      stop(sprintf(
        "'%s' has %d shared features but SharedFeatureSpace dimension is %d",
        name, ncol(block$values), shared_space$dimension
      ), call. = FALSE)
    }
  }
  invisible(block)
}
