#' Save and Load Alignment Models
#'
#' Serialize alignment models for storage and later reuse.
#'
#' @name serialize
#' @family serialization
NULL


#' Save an Alignment Model
#'
#' Save an alignment model to disk. The model can be loaded later with
#' \code{\link{load_alignment}} and applied to new data.
#'
#' @family serialization
#'
#' @param model An AlignmentModel or AlignmentResult.
#' @param path File path for saving. If doesn't end in .rds, .rds is appended.
#' @param compress Compression type: "gzip" (default), "bzip2", "xz", or FALSE.
#' @param include_data Logical; if TRUE and model is an AlignmentResult,
#'   also save the aligned data.
#'
#' @return Invisibly returns the path where the model was saved.
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
  # Extract model from result if needed
  save_obj <- if (inherits(model, "AlignmentResult")) {
    if (include_data) {
      model
    } else {
      get_model(model)
    }
  } else if (inherits(model, "AlignmentModel")) {
    model
  } else {
    stop("'model' must be an AlignmentModel or AlignmentResult", call. = FALSE)
  }

  # Ensure .rds extension
  if (!grepl("\\.rds$", path, ignore.case = TRUE)) {
    path <- paste0(path, ".rds")
  }

  # Add metadata for verification on load
  save_data <- list(
    object = save_obj,
    neuralign_version = as.character(utils::packageVersion("neuralign")),
    r_version = R.version.string,
    saved_at = Sys.time(),
    object_class = class(save_obj)[1]
  )

  # Compute hash for integrity checking
  save_data$hash <- digest::digest(save_obj, algo = "md5")

  saveRDS(save_data, file = path, compress = compress)

  message(sprintf("Saved %s to %s", class(save_obj)[1], path))

  invisible(path)
}


#' Load an Alignment Model
#'
#' Load a previously saved alignment model from disk.
#'
#' @family serialization
#'
#' @param path File path to the saved model.
#' @param verify Logical; if TRUE, verify the loaded model's integrity.
#'
#' @return The loaded AlignmentModel or AlignmentResult.
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

  save_data <- readRDS(path)

  # Check if this is a neuralign save format
  if (!is.list(save_data) || !"object" %in% names(save_data)) {
    # Maybe it's a raw model saved directly
    if (inherits(save_data, "AlignmentModel") ||
        inherits(save_data, "AlignmentResult")) {
      warning("Loaded model was not saved with save_alignment; ",
        "cannot verify integrity",
        call. = FALSE
      )
      return(save_data)
    }
    stop("File does not contain an alignment model", call. = FALSE)
  }

  obj <- save_data$object

  # Verify integrity
  if (verify && !is.null(save_data$hash)) {
    current_hash <- digest::digest(obj, algo = "md5")
    if (current_hash != save_data$hash) {
      warning("Model integrity check failed; file may be corrupted", call. = FALSE)
    }
  }

  # Version warning
  if (!is.null(save_data$neuralign_version)) {
    current_version <- as.character(utils::packageVersion("neuralign"))
    saved_version <- save_data$neuralign_version

    if (current_version != saved_version) {
      message(sprintf(
        "Note: Model was saved with neuralign %s, current version is %s",
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
