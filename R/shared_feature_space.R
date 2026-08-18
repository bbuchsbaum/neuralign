#' Shared Feature Space
#'
#' Identifies one fitted or externally defined aligned coordinate system.
#' Columns of an [AlignedStudy] live in this space. Spaces are compatible only
#' when their coordinate identities match.
#'
#' @param dimension Integer shared-space dimensionality.
#' @param coordinate_id Non-empty identifier for the coordinate-defining
#'   artifact. Aligner-produced spaces derive this from the fitted model.
#' @param method Alignment method recorded as lineage.
#' @param feature_data Optional data frame with one row per shared feature.
#' @param metric Declared metric on the shared space.
#' @param anchor Reference or anchor description.
#' @param canonicalization Description of sign, order, or gauge conventions.
#' @param decoder_capabilities Character vector of declared decode targets.
#' @param source_space Optional source-space description.
#' @param target_space Optional target-space description.
#' @param train_subjects Subjects that contributed to the fitted space.
#' @param extras Additional metadata. This is lineage and is not used to decide
#'   coordinate compatibility.
#'
#' @return A `SharedFeatureSpace` object.
#' @family aligned_study
#' @export
SharedFeatureSpace <- function(dimension,
                               coordinate_id,
                               method = NA_character_,
                               feature_data = NULL,
                               metric = "euclidean",
                               anchor = NULL,
                               canonicalization = list(),
                               decoder_capabilities = character(0),
                               source_space = NULL,
                               target_space = NULL,
                               train_subjects = character(0),
                               extras = list()) {
  if (missing(dimension) || length(dimension) != 1L || is.na(dimension) ||
      !is.numeric(dimension) || !is.finite(dimension) || dimension < 0 ||
      dimension > .Machine$integer.max || dimension != floor(dimension)) {
    stop("'dimension' must be a single non-negative integer", call. = FALSE)
  }
  dimension <- as.integer(dimension)

  if (missing(coordinate_id) || !is.character(coordinate_id) ||
      length(coordinate_id) != 1L || is.na(coordinate_id) ||
      !nzchar(coordinate_id)) {
    stop(
      "'coordinate_id' must identify the fitted or external coordinate artifact",
      call. = FALSE
    )
  }
  if (!is.character(method) || length(method) != 1L ||
      (!is.na(method) && !nzchar(method))) {
    stop("'method' must be NA or a non-empty character string", call. = FALSE)
  }
  if (!is.character(metric) || length(metric) != 1L || is.na(metric) ||
      !nzchar(metric)) {
    stop("'metric' must be a non-empty character string", call. = FALSE)
  }
  if (!is.list(canonicalization)) {
    stop("'canonicalization' must be a list", call. = FALSE)
  }
  if (!is.character(decoder_capabilities) || anyNA(decoder_capabilities)) {
    stop("'decoder_capabilities' must be a character vector", call. = FALSE)
  }
  if (!is.character(train_subjects) || anyNA(train_subjects)) {
    stop("'train_subjects' must be a character vector", call. = FALSE)
  }
  if (!is.list(extras)) {
    stop("'extras' must be a list", call. = FALSE)
  }

  if (is.null(feature_data)) {
    feature_data <- data.frame(
      feature_id = if (dimension) {
        sprintf("shared_%03d", seq_len(dimension))
      } else {
        character(0)
      },
      stringsAsFactors = FALSE
    )
  }
  if (!is.data.frame(feature_data)) {
    stop("'feature_data' must be a data.frame or NULL", call. = FALSE)
  }
  if (nrow(feature_data) != dimension) {
    stop(sprintf(
      "'feature_data' must have %d rows, got %d",
      dimension, nrow(feature_data)
    ), call. = FALSE)
  }

  obj <- list(
    dimension = dimension,
    coordinate_id = coordinate_id,
    method = method,
    feature_data = feature_data,
    metric = metric,
    anchor = anchor,
    canonicalization = canonicalization,
    decoder_capabilities = sort(unique(decoder_capabilities)),
    source_space = source_space,
    target_space = target_space,
    train_subjects = sort(unique(train_subjects)),
    extras = extras
  )
  class(obj) <- c("SharedFeatureSpace", "list")
  obj$schema_id <- .shared_space_schema_id(obj)
  obj$id <- shared_space_id(obj)
  obj
}


#' @export
print.SharedFeatureSpace <- function(x, ...) {
  cat("SharedFeatureSpace\n")
  cat(sprintf("  id: %s\n", x$id))
  cat(sprintf("  coordinate_id: %s\n", x$coordinate_id))
  cat(sprintf("  dimension: %d\n", x$dimension))
  cat(sprintf("  method: %s\n", x$method))
  cat(sprintf("  metric: %s\n", x$metric))
  invisible(x)
}


#' Compute a Shared-Space Identifier
#'
#' The identifier combines the coordinate artifact with the feature schema.
#'
#' @param x A [SharedFeatureSpace].
#' @return A character identifier prefixed with `"shared-"`.
#' @family aligned_study
shared_space_id <- function(x) {
  if (!inherits(x, "SharedFeatureSpace")) {
    stop("'x' must be a SharedFeatureSpace", call. = FALSE)
  }
  payload <- list(
    coordinate_id = x$coordinate_id,
    schema_id = x$schema_id %||% .shared_space_schema_id(x)
  )
  paste0("shared-", digest::digest(payload, algo = "xxhash64"))
}


#' Build a SharedFeatureSpace From an AlignmentModel
#'
#' @param model An [AlignmentModel].
#' @param dimension Optional shared dimension. By default it is inferred from
#'   the fitted transforms.
#' @param feature_data Optional shared-feature metadata.
#' @param coordinate_id Optional explicit coordinate artifact identifier.
#' @param extras Additional lineage metadata.
#' @param ... Passed to [SharedFeatureSpace()].
#' @return A [SharedFeatureSpace].
#' @family aligned_study
shared_feature_space_from_model <- function(model,
                                            dimension = NULL,
                                            feature_data = NULL,
                                            coordinate_id = NULL,
                                            extras = list(),
                                            ...) {
  model <- .ensure_model(model, what = "model")
  if (is.null(dimension)) {
    dimension <- .infer_shared_dimension(model)
  }
  if (!is.list(extras)) {
    stop("'extras' must be a list", call. = FALSE)
  }
  if (is.null(coordinate_id)) {
    coordinate_id <- .model_coordinate_id(model)
  }

  SharedFeatureSpace(
    dimension = dimension,
    coordinate_id = coordinate_id,
    method = model@method,
    feature_data = feature_data,
    anchor = model@reference,
    source_space = model@space_from,
    target_space = model@space_to,
    train_subjects = model@train_subjects,
    extras = utils::modifyList(list(
      neuralign_version = model@provenance$neuralign_version %||% NULL
    ), extras),
    ...
  )
}


#' Check Shared-Space Compatibility
#'
#' Unknown identities fail closed. This prevents values with missing
#' coordinate provenance from being combined silently.
#'
#' @param a,b [SharedFeatureSpace] objects or non-empty IDs.
#' @return `TRUE` only when both identities are known and equal.
#' @family aligned_study
#' @export
shared_spaces_compatible <- function(a, b) {
  id_a <- .as_shared_space_id(a)
  id_b <- .as_shared_space_id(b)
  if (is.null(id_a) || is.null(id_b)) {
    stop("Shared-space identity is unknown; compatibility cannot be established",
         call. = FALSE)
  }
  identical(id_a, id_b)
}


.shared_space_schema_id <- function(x) {
  payload <- list(
    dimension = as.integer(x$dimension),
    metric = x$metric,
    feature_data = .canonicalize_for_digest(x$feature_data),
    canonicalization = .canonicalize_for_digest(x$canonicalization),
    target_space = .canonicalize_for_digest(x$target_space)
  )
  paste0("schema-", digest::digest(payload, algo = "xxhash64"))
}


.model_coordinate_id <- function(model) {
  declared <- model@provenance$shared_space_coordinate_id %||%
    model@method_state$shared_space_coordinate_id %||% NULL
  if (!is.null(declared)) {
    if (!is.character(declared) || length(declared) != 1L ||
        is.na(declared) || !nzchar(declared)) {
      stop("Model shared-space coordinate identifier is invalid", call. = FALSE)
    }
    return(declared)
  }

  transforms <- model@transforms
  if (length(transforms)) {
    ord <- order(names(transforms))
    transforms <- transforms[ord]
  }
  payload <- list(
    method = model@method,
    reference = .canonicalize_for_digest(model@reference),
    reference_data = .canonicalize_for_digest(model@reference_data),
    target_space = .canonicalize_for_digest(model@space_to),
    method_state = .canonicalize_for_digest(model@method_state),
    transforms = lapply(transforms, .canonicalize_for_digest)
  )
  paste0("model-", digest::digest(payload, algo = "xxhash64"))
}


.canonicalize_for_digest <- function(x) {
  if (is.null(x)) return(NULL)
  if (is.factor(x)) return(as.character(x))
  if (inherits(x, "Date") || inherits(x, "POSIXt")) return(as.character(x))
  if (is.data.frame(x)) {
    out <- lapply(x, .canonicalize_for_digest)
    names(out) <- names(x)
    return(list(names = names(x), rows = nrow(x), columns = out))
  }
  if (is.matrix(x) || inherits(x, "Matrix")) {
    return(list(
      type = "matrix",
      dim = dim(x),
      hash = digest::digest(as.matrix(x), algo = "xxhash64")
    ))
  }
  if (is.atomic(x) && !is.object(x)) return(unname(x))
  if (inherits(x, "gds_space")) {
    return(list(type = "gds_space", name = x$name %||% ""))
  }
  if (is.list(x)) {
    out <- lapply(x, .canonicalize_for_digest)
    nms <- names(x)
    if (!is.null(nms)) {
      nms[is.na(nms) | !nzchar(nms)] <- paste0("..", which(is.na(nms) | !nzchar(nms)))
      names(out) <- nms
      out <- out[order(names(out))]
    }
    return(out)
  }
  list(type = class(x)[[1L]], hash = digest::digest(x, algo = "xxhash64"))
}


.as_shared_space_id <- function(x) {
  if (is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)) return(x)
  if (inherits(x, "SharedFeatureSpace")) return(x$id)
  NULL
}


.infer_shared_dimension <- function(model) {
  if (!length(model@transforms)) {
    stop("Cannot infer shared dimension: model has no transforms", call. = FALSE)
  }
  for (transform in model@transforms) {
    dimension <- .transform_target_dim(transform)
    if (!is.na(dimension)) return(as.integer(dimension))
  }
  stop("Cannot infer shared dimension from model transforms", call. = FALSE)
}


.transform_target_dim <- function(transform) {
  if (.is_embedding_transform(transform)) return(as.integer(nrow(transform$aligned)))
  if (.is_low_rank_transform(transform)) return(as.integer(nrow(transform$U)))
  if (.is_matrixish(transform)) return(as.integer(nrow(transform)))
  NA_integer_
}


.transform_source_dim <- function(transform) {
  if (.is_embedding_transform(transform)) return(NA_integer_)
  if (.is_low_rank_transform(transform)) return(as.integer(nrow(transform$V)))
  if (.is_matrixish(transform)) return(as.integer(ncol(transform)))
  NA_integer_
}
