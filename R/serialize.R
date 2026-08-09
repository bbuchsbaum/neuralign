#' Save and Load Alignment Artifacts
#'
#' Serialize alignment models, results, studies, and retained resample sets for
#' storage and later reuse.
#'
#' @name serialize
#' @family serialization
NULL

.NEURALIGN_SERIALIZATION_FORMAT <- "neuralign-alignment"
.NEURALIGN_SERIALIZATION_VERSION <- 2L
.NEURALIGN_SERIALIZABLE_CLASSES <- c(
  "AlignmentModel",
  "AlignmentResult",
  "AlignedStudy",
  "AlignedResampleSet"
)

.serialization_object_kind <- function(x) {
  kinds <- .NEURALIGN_SERIALIZABLE_CLASSES[
    vapply(
      .NEURALIGN_SERIALIZABLE_CLASSES,
      function(kind) inherits(x, kind),
      logical(1)
    )
  ]
  if (!length(kinds)) NULL else kinds[[1L]]
}

.validate_serializable_object <- function(x, context = "object") {
  kind <- .serialization_object_kind(x)
  if (is.null(kind)) {
    stop(
      sprintf(
        "'%s' must be an AlignmentModel, AlignmentResult, AlignedStudy, or AlignedResampleSet",
        context
      ),
      call. = FALSE
    )
  }
  validity <- methods::validObject(x, test = TRUE)
  if (!isTRUE(validity)) {
    stop(
      sprintf("Invalid %s: %s", kind, paste(validity, collapse = "; ")),
      call. = FALSE
    )
  }
  kind
}

.alignment_object_hash <- function(x, algorithm = "sha256") {
  digest::digest(x, algo = algorithm, serialize = TRUE)
}


#' Save an Alignment Artifact
#'
#' Save an alignment artifact to a versioned RDS envelope. The artifact can be
#' loaded later with \code{\link{load_alignment}}. Reusable models can then be
#' applied to new data when their aligner capabilities permit it.
#'
#' @family serialization
#'
#' @param model An AlignmentModel, AlignmentResult, AlignedStudy, or
#'   AlignedResampleSet.
#' @param path File path for saving. If doesn't end in .rds, .rds is appended.
#' @param compress Compression type: "gzip" (default), "bzip2", "xz", or FALSE.
#' @param include_data Logical; if TRUE and model is an AlignmentResult, also
#'   save the aligned data. AlignedStudy and AlignedResampleSet objects always
#'   retain their analysis-facing data.
#'
#' @return Invisibly returns the path where the artifact was saved.
#'
#' @examples
#' set.seed(1)
#' adat <- AlignmentData(list(
#'   s1 = matrix(rnorm(30), 6, 5),
#'   s2 = matrix(rnorm(30), 6, 5)
#' ))
#' res <- fit_alignment(adat, method = "procrustes", reference = "s1", compute_quality = FALSE)
#'
#' path <- tempfile(fileext = ".rds")
#' save_alignment(res, path)
#' mdl <- load_alignment(path)
#'
#' new_res <- apply_alignment(mdl, AlignmentData(list(s3 = matrix(rnorm(30), 6, 5))))
#' unlink(path)
#'
#' @export
save_alignment <- function(model,
                           path,
                           compress = "gzip",
                           include_data = FALSE) {
  if (!is.logical(include_data) || length(include_data) != 1L || is.na(include_data)) {
    stop("'include_data' must be TRUE or FALSE", call. = FALSE)
  }

  # Preserve the historical default: saving an AlignmentResult stores its
  # reusable model unless the caller explicitly asks to retain aligned data.
  save_obj <- if (inherits(model, "AlignmentResult")) {
    if (include_data) {
      model
    } else {
      get_model(model)
    }
  } else if (inherits(model, "AlignmentModel") ||
             inherits(model, "AlignedStudy") ||
             inherits(model, "AlignedResampleSet")) {
    model
  } else {
    .validate_serializable_object(model, context = "model")
  }
  object_kind <- .validate_serializable_object(save_obj, context = "model")

  # Ensure .rds extension
  if (!grepl("\\.rds$", path, ignore.case = TRUE)) {
    path <- paste0(path, ".rds")
  }

  # Version 2 is an explicit, fail-closed envelope. The object hash is separate
  # from descriptive metadata so package/R version messages do not affect it.
  save_data <- list(
    format_id = .NEURALIGN_SERIALIZATION_FORMAT,
    format_version = .NEURALIGN_SERIALIZATION_VERSION,
    object = save_obj,
    object_kind = object_kind,
    package_version = as.character(utils::packageVersion("neuralign")),
    r_version = R.version.string,
    saved_at = Sys.time(),
    integrity = list(
      algorithm = "sha256",
      hash = .alignment_object_hash(save_obj, "sha256")
    )
  )

  saveRDS(save_data, file = path, compress = compress)

  message(sprintf("Saved %s to %s", object_kind, path))

  invisible(path)
}


#' Load an Alignment Artifact
#'
#' Load a previously saved alignment artifact from disk. Version-2 envelopes
#' are checked for object kind, supported format version, S4 validity, and
#' SHA-256 integrity. Legacy version-1 envelopes and raw RDS objects have
#' explicit compatibility paths.
#'
#' @family serialization
#'
#' @param path File path to the saved model.
#' @param verify Logical; if TRUE, verify the loaded model's integrity.
#'
#' @return The loaded AlignmentModel, AlignmentResult, AlignedStudy, or
#'   AlignedResampleSet.
#'
#' @examples
#' set.seed(1)
#' adat <- AlignmentData(list(
#'   s1 = matrix(rnorm(30), 6, 5),
#'   s2 = matrix(rnorm(30), 6, 5)
#' ))
#' res <- fit_alignment(adat, method = "procrustes", reference = "s1", compute_quality = FALSE)
#' path <- tempfile(fileext = ".rds")
#' save_alignment(res, path)
#' mdl <- load_alignment(path)
#' unlink(path)
#' inherits(mdl, "AlignmentModel") || inherits(mdl, "AlignmentResult")
#'
#' @export
load_alignment <- function(path, verify = TRUE) {
  if (!file.exists(path)) {
    stop(sprintf("File not found: %s", path), call. = FALSE)
  }

  if (!is.logical(verify) || length(verify) != 1L || is.na(verify)) {
    stop("'verify' must be TRUE or FALSE", call. = FALSE)
  }

  save_data <- readRDS(path)

  # Raw RDS objects predate the envelope. Retain this compatibility path, but
  # make the lack of integrity/version metadata explicit.
  if (!is.list(save_data) || !"object" %in% names(save_data)) {
    if (!is.null(.serialization_object_kind(save_data))) {
      warning("Loaded object was not saved with save_alignment; ",
        "cannot verify integrity",
        call. = FALSE
      )
      .validate_serializable_object(save_data, context = "loaded object")
      return(save_data)
    }
    stop("File does not contain a supported alignment object", call. = FALSE)
  }

  is_v2 <- identical(save_data$format_id %||% NULL, .NEURALIGN_SERIALIZATION_FORMAT)
  if (is_v2) {
    version <- save_data$format_version %||% NA_integer_
    if (!is.numeric(version) || length(version) != 1L || is.na(version) ||
        version != as.integer(version) || version < 1L) {
      stop("Alignment file has an invalid serialization format version", call. = FALSE)
    }
    version <- as.integer(version)
    if (version > .NEURALIGN_SERIALIZATION_VERSION) {
      stop(
        sprintf(
          "Alignment file format version %d is newer than supported version %d; update neuralign before loading it",
          version, .NEURALIGN_SERIALIZATION_VERSION
        ),
        call. = FALSE
      )
    }
    if (version < .NEURALIGN_SERIALIZATION_VERSION) {
      stop(
        sprintf("No migration is registered for alignment format version %d", version),
        call. = FALSE
      )
    }

    obj <- save_data$object
    kind <- .validate_serializable_object(obj, context = "loaded object")
    if (!identical(save_data$object_kind %||% NULL, kind)) {
      stop("Alignment file object-kind metadata does not match its object", call. = FALSE)
    }

    if (verify) {
      integrity <- save_data$integrity %||% NULL
      if (!is.list(integrity) ||
          !identical(integrity$algorithm %||% NULL, "sha256") ||
          !is.character(integrity$hash) || length(integrity$hash) != 1L) {
        stop("Alignment file is missing valid SHA-256 integrity metadata", call. = FALSE)
      }
      current_hash <- .alignment_object_hash(obj, integrity$algorithm)
      if (!identical(current_hash, integrity$hash)) {
        stop("Alignment integrity check failed; file may be corrupted", call. = FALSE)
      }
    }

    saved_version <- save_data$package_version %||% NULL
  } else {
    # Version 1 envelope written by neuralign <= 0.1.0.
    obj <- save_data$object
    .validate_serializable_object(obj, context = "loaded object")
    if (verify && !is.null(save_data$hash)) {
      current_hash <- digest::digest(obj, algo = "md5")
      if (!identical(current_hash, save_data$hash)) {
        stop("Legacy alignment integrity check failed; file may be corrupted", call. = FALSE)
      }
    }
    saved_version <- save_data$neuralign_version %||% NULL
  }

  if (!is.null(saved_version)) {
    current_version <- as.character(utils::packageVersion("neuralign"))

    if (current_version != saved_version) {
      message(sprintf(
        "Note: Object was saved with neuralign %s, current version is %s",
        saved_version, current_version
      ))
    }
  }

  obj
}


#' Export Model to Portable Format
#'
#' Export alignment transforms in a portable format (e.g., for use in
#' other software or languages).
#'
#' @family serialization
#'
#' @param model An AlignmentModel.
#' @param path Output path (without extension).
#' @param format Export format:
#'   \itemize{
#'     \item "csv" - CSV files (one per transform)
#'     \item "json" - JSON file
#'     \item "mat" - MATLAB format (requires R.matlab)
#'   }
#'
#' @return Invisibly returns list of created files.
#'
#' @examples
#' Q <- diag(3)
#' mdl <- AlignmentModel(list(s1 = Q), reference = "s1", method = "procrustes")
#' base <- tempfile("neuralign_export_")
#' export_alignment(mdl, base, format = "csv")
#' unlink(paste0(base, "_transforms"), recursive = TRUE)
#'
#' @export
export_alignment <- function(model,
                             path,
                             format = c("csv", "json", "mat")) {
  format <- match.arg(format)

  model <- .ensure_model(model, what = "model")

  files <- character(0)

  if (format == "csv") {
    # Create directory if needed
    dir <- paste0(path, "_transforms")
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
    }

    # Export each transform
    for (subj in names(model@transforms)) {
      mat <- as.matrix(model@transforms[[subj]])
      file <- file.path(dir, sprintf("%s.csv", subj))
      utils::write.csv(mat, file, row.names = FALSE)
      files <- c(files, file)
    }

    # Export metadata
    meta_file <- file.path(dir, "_metadata.csv")
    meta <- data.frame(
      subject = names(model@transforms),
      nrow = vapply(model@transforms, nrow, integer(1)),
      ncol = vapply(model@transforms, ncol, integer(1)),
      stringsAsFactors = FALSE
    )
    utils::write.csv(meta, meta_file, row.names = FALSE)
    files <- c(files, meta_file)
  }

  if (format == "json") {
    if (!requireNamespace("jsonlite", quietly = TRUE)) {
      stop("Package 'jsonlite' required for JSON export", call. = FALSE)
    }

    json_data <- list(
      method = model@method,
      subjects = names(model@transforms),
      transforms = lapply(model@transforms, function(x) {
        as.matrix(x)
      }),
      reference = if (.is_matrixish(model@reference_data)) {
        as.matrix(model@reference_data)
      } else {
        as.character(model@reference)
      },
      provenance = model@provenance
    )

    file <- paste0(path, ".json")
    json_str <- jsonlite::toJSON(json_data, auto_unbox = TRUE, pretty = TRUE)
    writeLines(json_str, file)
    files <- c(files, file)
  }

  if (format == "mat") {
    if (!requireNamespace("R.matlab", quietly = TRUE)) {
      stop("Package 'R.matlab' required for MATLAB export", call. = FALSE)
    }

    mat_data <- list(
      method = model@method,
      subjects = names(model@transforms)
    )

    # Add transforms
    for (subj in names(model@transforms)) {
      mat_data[[paste0("T_", gsub("-", "_", subj))]] <-
        as.matrix(model@transforms[[subj]])
    }

    # Add reference
    if (.is_matrixish(model@reference_data)) {
      mat_data$reference <- as.matrix(model@reference_data)
    }

    file <- paste0(path, ".mat")
    do.call(R.matlab::writeMat, c(list(con = file), mat_data))
    files <- c(files, file)
  }

  message(sprintf("Exported %d files", length(files)))
  invisible(files)
}


#' Import Transforms from External Format
#'
#' Import alignment transforms from external files.
#'
#' @family serialization
#'
#' @param path Path to import from.
#' @param format Import format (same options as export_alignment).
#' @param method Method name to assign to the imported model.
#'
#' @return An AlignmentModel with the imported transforms.
#'
#' @examples
#' Q <- diag(3)
#' mdl <- AlignmentModel(list(s1 = Q), reference = "s1", method = "procrustes")
#' base <- tempfile("neuralign_export_")
#' export_alignment(mdl, base, format = "csv")
#' mdl2 <- import_alignment(base, format = "csv", method = "imported")
#' unlink(paste0(base, "_transforms"), recursive = TRUE)
#' inherits(mdl2, "AlignmentModel")
#'
#' @export
import_alignment <- function(path,
                             format = c("csv", "json"),
                             method = "imported") {
  format <- match.arg(format)

  if (format == "csv") {
    dir <- paste0(path, "_transforms")
    if (!dir.exists(dir)) {
      stop(sprintf("Directory not found: %s", dir), call. = FALSE)
    }

    # Read metadata
    meta_file <- file.path(dir, "_metadata.csv")
    if (!file.exists(meta_file)) {
      stop("Metadata file not found", call. = FALSE)
    }
    meta <- utils::read.csv(meta_file, stringsAsFactors = FALSE)

    # Read transforms
    transforms <- list()
    for (i in seq_len(nrow(meta))) {
      subj <- meta$subject[i]
      file <- file.path(dir, sprintf("%s.csv", subj))
      if (!file.exists(file)) {
        warning(sprintf("Transform file not found for %s", subj), call. = FALSE)
        next
      }
      transforms[[subj]] <- as.matrix(utils::read.csv(file))
    }
  }

  if (format == "json") {
    if (!requireNamespace("jsonlite", quietly = TRUE)) {
      stop("Package 'jsonlite' required for JSON import", call. = FALSE)
    }

    file <- if (grepl("\\.json$", path)) path else paste0(path, ".json")
    if (!file.exists(file)) {
      stop(sprintf("File not found: %s", file), call. = FALSE)
    }

    json_data <- jsonlite::fromJSON(file)

    transforms <- lapply(json_data$transforms, as.matrix)
    names(transforms) <- json_data$subjects

    if (!is.null(json_data$method)) {
      method <- json_data$method
    }
  }

  # Validate imported transforms
  if (length(transforms) == 0) {
    stop("No transforms were imported", call. = FALSE)
  }
  for (subj in names(transforms)) {
    mat <- transforms[[subj]]
    if (!is.matrix(mat) && !inherits(mat, "Matrix")) {
      stop(sprintf("Imported transform for subject '%s' is not a matrix", subj),
        call. = FALSE)
    }
    if (!is.numeric(mat)) {
      stop(sprintf("Imported transform for subject '%s' is not numeric", subj),
        call. = FALSE)
    }
    if (any(!is.finite(mat))) {
      warning(sprintf(
        "Imported transform for subject '%s' contains non-finite values", subj
      ), call. = FALSE)
    }
  }

  AlignmentModel(
    transforms = transforms,
    reference = NULL,
    method = method,
    train_subjects = names(transforms)
  )
}
