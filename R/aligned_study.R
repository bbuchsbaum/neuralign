#' AlignedStudy Class
#'
#' Analysis-ready collection of scientific observations in one common shared
#' space. Distinct from [AlignmentResult], which is an algorithm execution
#' artifact.
#'
#' Orientation: each block's `values` are **observations × shared features**.
#'
#' @slot blocks Named list of [AlignedBlock] objects.
#' @slot shared_space A [SharedFeatureSpace].
#' @slot model The [AlignmentModel] (or empty model) that produced the data.
#' @slot subject_data Optional `data.frame` with one row per subject/session.
#' @slot lineage Extraction → preprocessing → alignment history (list).
#' @slot safety Fit/application/leakage declarations (list).
#' @slot storage Eager storage description. This is not a lazy-backend API.
#' @slot metadata Additional named list.
#' @param object,x An `AlignedStudy` object.
#'
#' @family aligned_study
#' @export
setClass("AlignedStudy",
  slots = c(
    blocks = "list",
    shared_space = "ANY",
    model = "ANY",
    subject_data = "ANY",
    lineage = "list",
    safety = "list",
    storage = "list",
    metadata = "list"
  ),
  prototype = list(
    blocks = list(),
    shared_space = NULL,
    model = NULL,
    subject_data = NULL,
    lineage = list(),
    safety = list(),
    storage = list(mode = "eager"),
    metadata = list()
  )
)

setValidity("AlignedStudy", function(object) {
  errors <- character()
  if (!is.list(object@blocks)) {
    errors <- c(errors, "'blocks' must be a list")
  }
  block_names <- names(object@blocks)
  if (length(object@blocks) > 0 &&
      (is.null(block_names) || anyNA(block_names) || any(!nzchar(block_names)) ||
       anyDuplicated(block_names))) {
    errors <- c(errors, "'blocks' must have unique, non-empty names")
  }
  if (!is.null(object@model) && !inherits(object@model, "AlignmentModel")) {
    errors <- c(errors, "'model' must be an AlignmentModel or NULL")
  } else if (inherits(object@model, "AlignmentModel") &&
             inherits(object@shared_space, "SharedFeatureSpace") &&
             !identical(
               object@shared_space$coordinate_id,
               .model_coordinate_id(object@model)
             )) {
    errors <- c(errors, "'model' and 'shared_space' identify different coordinates")
  }
  if (!is.data.frame(object@subject_data)) {
    errors <- c(errors, "'subject_data' must be a data.frame")
  } else {
    if (nrow(object@subject_data) != length(object@blocks)) {
      errors <- c(errors, "'subject_data' must have one row per block")
    }
    if (!"block_id" %in% names(object@subject_data)) {
      errors <- c(errors, "'subject_data' must contain 'block_id'")
    } else if (!identical(as.character(object@subject_data$block_id), block_names)) {
      errors <- c(errors, "'subject_data$block_id' must match block names in order")
    }
  }
  if (!is.list(object@lineage) || !is.list(object@safety) ||
      !is.list(object@storage) || !is.list(object@metadata)) {
    errors <- c(errors, "lineage, safety, storage, and metadata must be lists")
  }
  if (!identical(object@storage$mode %||% "eager", "eager")) {
    errors <- c(errors, "only storage mode 'eager' is currently supported")
  }
  if (!inherits(object@shared_space, "SharedFeatureSpace")) {
    errors <- c(errors, "'shared_space' must be a SharedFeatureSpace")
  } else {
    sid <- object@shared_space$id
    for (nm in names(object@blocks)) {
      b <- object@blocks[[nm]]
      if (!inherits(b, "AlignedBlock")) {
        errors <- c(errors, sprintf("blocks[['%s']] is not an AlignedBlock", nm))
        next
      }
      if (!is.null(b$shared_space_id) && !identical(b$shared_space_id, sid)) {
        errors <- c(errors, sprintf(
          "blocks[['%s']] space id mismatch", nm
        ))
      }
      if (ncol(b$values) != object@shared_space$dimension) {
        errors <- c(errors, sprintf(
          "blocks[['%s']] feature dim != shared_space$dimension", nm
        ))
      }
    }
  }
  if (length(errors) == 0L) TRUE else errors
})


#' Create an AlignedStudy
#'
#' @param blocks Named list of [AlignedBlock]s.
#' @param shared_space A [SharedFeatureSpace].
#' @param model Optional [AlignmentModel].
#' @param subject_data Optional subject-level `data.frame`.
#' @param lineage Lineage list.
#' @param safety Safety / leakage declarations.
#' @param storage Storage metadata list.
#' @param metadata Extra metadata.
#'
#' @return An [AlignedStudy].
#'
#' @family aligned_study
#' @export
AlignedStudy <- function(blocks,
                         shared_space,
                         model = NULL,
                         subject_data = NULL,
                         lineage = list(),
                         safety = list(),
                         storage = list(mode = "eager"),
                         metadata = list()) {
  if (!inherits(shared_space, "SharedFeatureSpace")) {
    stop("'shared_space' must be a SharedFeatureSpace", call. = FALSE)
  }
  if (!is.list(blocks) || (length(blocks) > 0 && is.null(names(blocks)))) {
    stop("'blocks' must be a named list of AlignedBlock objects", call. = FALSE)
  }
  sid <- shared_space$id
  for (nm in names(blocks)) {
    .validate_aligned_block(blocks[[nm]], shared_space = shared_space, name = nm)
    if (is.null(blocks[[nm]]$shared_space_id)) {
      blocks[[nm]]$shared_space_id <- sid
    }
  }
  if (!is.null(model) && !inherits(model, "AlignmentModel")) {
    stop("'model' must be an AlignmentModel or NULL", call. = FALSE)
  }
  if (inherits(model, "AlignmentModel") &&
      !identical(shared_space$coordinate_id, .model_coordinate_id(model))) {
    stop("'model' and 'shared_space' identify different coordinates",
         call. = FALSE)
  }
  if (!is.null(subject_data) && !is.data.frame(subject_data)) {
    stop("'subject_data' must be a data.frame or NULL", call. = FALSE)
  }
  if (is.null(subject_data) && length(blocks) > 0) {
    subject_data <- data.frame(
      block_id = names(blocks),
      subject_id = vapply(blocks, function(b) b$subject_id, character(1)),
      session_id = vapply(blocks, function(b) b$session_id %||% NA_character_,
                          character(1)),
      stringsAsFactors = FALSE
    )
  } else if (is.null(subject_data)) {
    subject_data <- data.frame(
      block_id = character(0),
      subject_id = character(0),
      session_id = character(0),
      stringsAsFactors = FALSE
    )
  } else {
    if (nrow(subject_data) != length(blocks)) {
      stop("'subject_data' must have one row per block", call. = FALSE)
    }
    if (!"block_id" %in% names(subject_data)) {
      subject_data$block_id <- names(blocks)
    }
    if (!identical(as.character(subject_data$block_id), names(blocks))) {
      stop("'subject_data$block_id' must match block names in order", call. = FALSE)
    }
  }
  if (!is.list(lineage) || !is.list(safety) || !is.list(storage) ||
      !is.list(metadata)) {
    stop("lineage, safety, storage, and metadata must be lists", call. = FALSE)
  }
  if (!identical(storage$mode %||% "eager", "eager")) {
    stop("only storage mode 'eager' is currently supported", call. = FALSE)
  }

  lineage$alignment_artifact_id <- lineage$alignment_artifact_id %||%
    shared_space$coordinate_id

  new("AlignedStudy",
    blocks = blocks,
    shared_space = shared_space,
    model = model,
    subject_data = subject_data,
    lineage = lineage,
    safety = safety,
    storage = storage,
    metadata = metadata
  )
}


#' @rdname AlignedStudy-class
#' @export
setMethod("show", "AlignedStudy", function(object) {
  cat("AlignedStudy\n")
  ss <- object@shared_space
  cat(sprintf("  shared_space: %s\n", ss$id %||% "<unset>"))
  cat(sprintf("  dimension: %d\n", ss$dimension %||% NA_integer_))
  cat(sprintf("  blocks: %d\n", length(object@blocks)))
  if (!is.null(object@model) && inherits(object@model, "AlignmentModel")) {
    cat(sprintf("  model method: %s\n", object@model@method))
  }
  mode <- object@safety$mode %||% object@safety$alignment_use$mode
  if (!is.null(mode)) {
    cat(sprintf("  safety mode: %s\n", mode))
  }
  cat(sprintf("  storage: %s\n", object@storage$mode %||% "eager"))
})


#' @rdname AlignedStudy-class
#' @export
setMethod("length", "AlignedStudy", function(x) length(x@blocks))


# ---- accessors --------------------------------------------------------------

#' Subject IDs in an AlignedStudy
#'
#' @param x An [AlignedStudy] (or compatible).
#' @return Character vector of subject ids (one per block, in block order).
#' @family aligned_study
#' @export
subject_ids <- function(x) {
  UseMethod("subject_ids")
}

#' @export
subject_ids.AlignedStudy <- function(x) {
  unname(vapply(x@blocks, function(b) b$subject_id, character(1)))
}

#' @export
subject_ids.default <- function(x) {
  stop("subject_ids() not implemented for class ", paste(class(x), collapse = "/"),
       call. = FALSE)
}


#' Subject-Level Metadata
#'
#' @param x An [AlignedStudy].
#' @return A `data.frame`.
#' @family aligned_study
#' @export
subject_data <- function(x) {
  UseMethod("subject_data")
}

#' @export
subject_data.AlignedStudy <- function(x) x@subject_data

#' @export
subject_data.default <- function(x) {
  stop("subject_data() not implemented for class ", paste(class(x), collapse = "/"),
       call. = FALSE)
}


#' Observation Metadata for a Subject Block
#'
#' @param x An [AlignedStudy].
#' @param subject Subject id (required when `x` has multiple blocks).
#' @return A `data.frame`.
#' @family aligned_study
#' @export
observation_data <- function(x, subject = NULL) {
  UseMethod("observation_data")
}

#' @export
observation_data.AlignedStudy <- function(x, subject = NULL) {
  block <- .get_study_block(x, subject)
  block$observation_data
}


#' Shared Feature Metadata
#'
#' @param x An [AlignedStudy] or [SharedFeatureSpace].
#' @return A `data.frame` with one row per shared feature.
#' @family aligned_study
#' @export
feature_data <- function(x) {
  UseMethod("feature_data")
}

#' @export
feature_data.AlignedStudy <- function(x) x@shared_space$feature_data

#' @export
feature_data.SharedFeatureSpace <- function(x) x$feature_data

#' @export
feature_data.default <- function(x) {
  stop("feature_data() not implemented for class ", paste(class(x), collapse = "/"),
       call. = FALSE)
}


#' Shared Feature Space Accessor
#'
#' @param x An [AlignedStudy] or [AlignedResampleSet] split.
#' @return A [SharedFeatureSpace].
#' @family aligned_study
#' @export
shared_space <- function(x) {
  UseMethod("shared_space")
}

#' @export
shared_space.AlignedStudy <- function(x) x@shared_space

#' @export
shared_space.default <- function(x) {
  stop("shared_space() not implemented for class ", paste(class(x), collapse = "/"),
       call. = FALSE)
}


#' Alignment Model Accessor
#'
#' @param x An [AlignedStudy].
#' @return An [AlignmentModel] or `NULL`.
#' @family aligned_study
#' @export
alignment_model <- function(x) {
  UseMethod("alignment_model")
}

#' @export
alignment_model.AlignedStudy <- function(x) x@model

#' @export
alignment_model.default <- function(x) {
  stop("alignment_model() not implemented for class ",
       paste(class(x), collapse = "/"), call. = FALSE)
}


#' Analysis Lineage Accessor
#'
#' @param x An [AlignedStudy].
#' @return Lineage list.
#' @family aligned_study
#' @export
analysis_lineage <- function(x) {
  UseMethod("analysis_lineage")
}

#' @export
analysis_lineage.AlignedStudy <- function(x) x@lineage


#' Analysis Safety Accessor
#'
#' @param x An [AlignedStudy].
#' @return Safety list.
#' @family aligned_study
#' @export
analysis_safety <- function(x) {
  UseMethod("analysis_safety")
}

#' @export
analysis_safety.AlignedStudy <- function(x) x@safety


#' Aligned Matrix for One Subject (Analysis Orientation)
#'
#' @param x An [AlignedStudy].
#' @param subject Subject id.
#' @return Observations × shared features matrix.
#' @family aligned_study
#' @export
aligned_matrix <- function(x, subject = NULL) {
  UseMethod("aligned_matrix")
}

#' @export
aligned_matrix.AlignedStudy <- function(x, subject = NULL) {
  .get_study_block(x, subject)$values
}


#' Map a Function Over Subject Blocks
#'
#' @param x An [AlignedStudy].
#' @param .f Function of one [AlignedBlock].
#' @param ... Passed to `.f`.
#' @return Named list of results.
#' @family aligned_study
#' @export
map_subjects <- function(x, .f, ...) {
  UseMethod("map_subjects")
}

#' @export
map_subjects.AlignedStudy <- function(x, .f, ...) {
  lapply(x@blocks, .f, ...)
}


#' Stack Subjects With Explicit Index
#'
#' Unlike silent `rbind`, this returns both the stacked matrix and a complete
#' row index preserving subject / observation identity.
#'
#' @param x An [AlignedStudy].
#' @param include_index Logical; if `FALSE`, return only the matrix.
#' @return If `include_index=TRUE`, a list with `matrix` and `index`. Otherwise
#'   the stacked matrix.
#' @family aligned_study
#' @export
stack_subjects <- function(x, include_index = TRUE) {
  UseMethod("stack_subjects")
}

#' @export
stack_subjects.AlignedStudy <- function(x, include_index = TRUE) {
  if (!length(x@blocks)) {
    mat <- matrix(numeric(0), 0, x@shared_space$dimension)
    idx <- data.frame(
      row_id = integer(0),
      subject_id = character(0),
      session_id = character(0),
      observation_id = character(0),
      stringsAsFactors = FALSE
    )
    if (isTRUE(include_index)) {
      return(list(matrix = mat, index = idx))
    }
    return(mat)
  }

  mats <- lapply(x@blocks, function(b) b$values)
  mat <- do.call(rbind, mats)

  idx_list <- lapply(x@blocks, function(b) {
    obs <- b$observation_data
    obs_id <- if ("observation_id" %in% names(obs)) {
      as.character(obs$observation_id)
    } else {
      as.character(seq_len(nrow(obs)))
    }
    data.frame(
      subject_id = b$subject_id,
      session_id = if (is.null(b$session_id)) NA_character_ else as.character(b$session_id),
      observation_id = obs_id,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
  idx <- do.call(rbind, idx_list)
  rownames(idx) <- NULL
  idx$row_id <- seq_len(nrow(idx))
  idx <- idx[, c("row_id", "subject_id", "session_id", "observation_id"), drop = FALSE]

  # Attach remaining observation columns with subject-disambiguated names if needed
  extra_cols <- setdiff(
    unique(unlist(lapply(x@blocks, function(b) names(b$observation_data)), use.names = FALSE)),
    "observation_id"
  )
  if (length(extra_cols)) {
    extras <- lapply(x@blocks, function(b) {
      obs <- b$observation_data
      out <- data.frame(row.names = NULL)
      for (nm in extra_cols) {
        out[[nm]] <- if (nm %in% names(obs)) obs[[nm]] else rep(NA, nrow(obs))
      }
      out
    })
    extra_df <- do.call(rbind, extras)
    rownames(extra_df) <- NULL
    idx <- cbind(idx, extra_df)
  }

  if (isTRUE(include_index)) {
    list(matrix = mat, index = idx)
  } else {
    mat
  }
}


.get_study_block <- function(x, subject = NULL) {
  blocks <- x@blocks
  if (!length(blocks)) {
    stop("AlignedStudy has no blocks", call. = FALSE)
  }
  if (is.null(subject)) {
    if (length(blocks) != 1L) {
      stop("'subject' is required when AlignedStudy has multiple blocks", call. = FALSE)
    }
    return(blocks[[1L]])
  }
  # Match by block name or subject_id
  if (subject %in% names(blocks)) {
    return(blocks[[subject]])
  }
  ids <- vapply(blocks, function(b) b$subject_id, character(1))
  hit <- which(ids == subject)
  if (!length(hit)) {
    stop(sprintf("Subject '%s' not found in AlignedStudy", subject), call. = FALSE)
  }
  if (length(hit) > 1L) {
    stop(sprintf(
      "Subject '%s' matches multiple blocks; use block names: %s",
      subject, paste(names(blocks)[hit], collapse = ", ")
    ), call. = FALSE)
  }
  blocks[[hit]]
}


#' Assert All Blocks Share One Shared-Space Id
#'
#' @param x An [AlignedStudy].
#' @return Invisibly `TRUE`, or error.
#' @family aligned_study
#' @export
assert_common_shared_space <- function(x) {
  if (!inherits(x, "AlignedStudy")) {
    stop("'x' must be an AlignedStudy", call. = FALSE)
  }
  sid <- x@shared_space$id
  ids <- vapply(x@blocks, function(b) {
    b$shared_space_id %||% sid
  }, character(1))
  if (any(ids != sid)) {
    stop(
      "AlignedStudy blocks do not all share shared_space$id '", sid, "'. ",
      "Use AlignedResampleSet for fold-specific spaces.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}
