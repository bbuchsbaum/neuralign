#' Guidance Channels
#'
#' Guidance channels provide optional anatomical/geometry priors that guide
#' functional alignment. Guidance is not itself "the alignment"; it supplies
#' constraints, costs, or initialization for methods that use it.
#'
#' Guidance is stored in `AlignmentData@metadata$guidance` as a named list keyed
#' by subject. Each subject entry is a (possibly empty) named list of guidance
#' channels. A channel is a list with at least:
#' - `type`: one of `"projector"`, `"coords"`, `"geometry"` (extensible)
#' - `value`: the channel payload (often a matrix/Matrix operator)
#'
#' @name guidance
NULL

.is_nonempty_string <- function(x) {
  is.character(x) && length(x) == 1L && nzchar(x)
}

.as_guidance_channel <- function(x, name = NULL) {
  if (!is.list(x)) {
    stop("Guidance channel must be a list", call. = FALSE)
  }
  if (is.null(x$type) || !.is_nonempty_string(x$type)) {
    stop("Guidance channel must have a non-empty 'type' field", call. = FALSE)
  }
  if (is.null(x$value)) {
    stop("Guidance channel must have a 'value' field", call. = FALSE)
  }
  if (is.null(x$name) && !is.null(name) && .is_nonempty_string(name)) {
    x$name <- name
  }
  x
}

.normalize_guidance_by_subject <- function(guidance, subjects) {
  if (is.null(guidance)) {
    out <- rep(list(list()), length(subjects))
    names(out) <- subjects
    return(out)
  }

  if (!is.list(guidance)) {
    stop("'guidance' must be a list (or NULL)", call. = FALSE)
  }

  if (is.null(names(guidance))) {
    stop("'guidance' must be a named list keyed by subject", call. = FALSE)
  }

  missing <- setdiff(subjects, names(guidance))
  if (length(missing) > 0) {
    stop(
      sprintf("guidance missing subjects: %s", paste(missing, collapse = ", ")),
      call. = FALSE
    )
  }

  extra <- setdiff(names(guidance), subjects)
  if (length(extra) > 0) {
    warning(
      sprintf("Ignoring guidance for unknown subjects: %s", paste(extra, collapse = ", ")),
      call. = FALSE
    )
  }

  out <- lapply(subjects, function(subj) {
    channels <- guidance[[subj]]
    if (is.null(channels)) return(list())
    if (!is.list(channels)) {
      stop(sprintf("guidance for subject '%s' must be a list", subj), call. = FALSE)
    }
    if (is.null(names(channels))) {
      # Allow unnamed channel lists, but assign stable names.
      names(channels) <- paste0("channel-", seq_along(channels))
    }
    out_ch <- lapply(names(channels), function(nm) .as_guidance_channel(channels[[nm]], name = nm))
    names(out_ch) <- names(channels)
    out_ch
  })
  names(out) <- subjects
  out
}

.validate_guidance_dims <- function(data, guidance_by_subject) {
  data_list <- get_data_list(data)
  for (subj in names(guidance_by_subject)) {
    x <- data_list[[subj]]
    if (!.is_matrixish(x)) x <- as.matrix(x)
    n_feat <- nrow(x)

    channels <- guidance_by_subject[[subj]]
    for (ch in channels) {
      typ <- ch$type
      val <- ch$value
      if (typ %in% c("projector", "coords")) {
        if (!.is_matrixish(val)) {
          stop(
            sprintf("Guidance channel '%s' for subject '%s' must be matrix-like for type '%s'", ch$name %||% "<unnamed>", subj, typ),
            call. = FALSE
          )
        }
        if (!inherits(val, "Matrix")) {
          val <- as.matrix(val)
        }
        if (typ == "projector" && ncol(val) != n_feat) {
          stop(
            sprintf(
              "Projector guidance dimension mismatch for subject '%s': ncol(P)=%d but subject has %d features",
              subj, ncol(val), n_feat
            ),
            call. = FALSE
          )
        }
        if (typ == "coords" && nrow(val) != n_feat) {
          stop(
            sprintf(
              "Coords guidance dimension mismatch for subject '%s': nrow(Z)=%d but subject has %d features",
              subj, nrow(val), n_feat
            ),
            call. = FALSE
          )
        }
      }
    }
  }
  invisible(TRUE)
}

#' Get Guidance Channels
#'
#' @family guidance
#'
#' @param data An `AlignmentData` object.
#' @param subject Optional subject id; if supplied, returns guidance only for
#'   that subject.
#' @param type Optional type filter (e.g., `"projector"`, `"coords"`).
#'
#' @return A named list keyed by subject (or a subject's channel list).
#'
#' @examples
#' set.seed(1)
#' adat <- AlignmentData(list(
#'   s1 = matrix(rnorm(12), 3, 4),
#'   s2 = matrix(rnorm(12), 3, 4)
#' ))
#' g <- list(
#'   s1 = list(roi = guidance_channel("projector", diag(3), name = "roi")),
#'   s2 = list(roi = guidance_channel("projector", diag(3), name = "roi"))
#' )
#' adat2 <- set_guidance(adat, g)
#' get_guidance(adat2, subject = "s1")
#' @export
get_guidance <- function(data, subject = NULL, type = NULL) {
  if (!inherits(data, "AlignmentData")) {
    stop("'data' must be an AlignmentData object", call. = FALSE)
  }
  g <- data@metadata$guidance %||% NULL
  if (is.null(g)) {
    g <- rep(list(list()), length(data@subjects))
    names(g) <- data@subjects
  }

  if (!is.null(subject)) {
    if (!subject %in% data@subjects) {
      stop(sprintf("Unknown subject '%s'", subject), call. = FALSE)
    }
    g <- g[[subject]] %||% list()
    if (!is.null(type)) {
      type <- as.character(type)
      g <- Filter(function(ch) is.list(ch) && identical(ch$type, type), g)
    }
    return(g)
  }

  if (!is.null(type)) {
    type <- as.character(type)
    g <- lapply(g, function(chs) Filter(function(ch) is.list(ch) && identical(ch$type, type), chs %||% list()))
  }

  g
}

#' Set Guidance Channels
#'
#' Replace the guidance channels stored in an `AlignmentData` object.
#'
#' @family guidance
#'
#' @param data An `AlignmentData` object.
#' @param guidance A named list keyed by subject; each entry is a named list of
#'   guidance channels.
#' @param validate Logical; if TRUE, validate basic per-subject dimensional
#'   compatibility between guidance channels and the subject's data matrix.
#'
#' @return The modified `AlignmentData`.
#'
#' @examples
#' set.seed(1)
#' adat <- AlignmentData(list(
#'   s1 = matrix(rnorm(12), 3, 4),
#'   s2 = matrix(rnorm(12), 3, 4)
#' ))
#' g <- list(
#'   s1 = list(roi = guidance_channel("projector", diag(3), name = "roi")),
#'   s2 = list(roi = guidance_channel("projector", diag(3), name = "roi"))
#' )
#' adat2 <- set_guidance(adat, g)
#' names(get_guidance(adat2, subject = "s1"))
#' @export
set_guidance <- function(data, guidance, validate = TRUE) {
  if (!inherits(data, "AlignmentData")) {
    stop("'data' must be an AlignmentData object", call. = FALSE)
  }

  g_by_subj <- .normalize_guidance_by_subject(guidance, data@subjects)
  if (isTRUE(validate)) {
    .validate_guidance_dims(data, g_by_subj)
  }

  data@metadata$guidance <- g_by_subj
  if (isTRUE(validate)) {
    validObject(data)
  }
  data
}

#' Guidance Channel Constructor
#'
#' Convenience constructor for a single guidance channel.
#'
#' @family guidance
#'
#' @param type Guidance type (`"projector"`, `"coords"`, `"geometry"`, ...).
#' @param value Payload for this channel.
#' @param name Optional channel name.
#' @param ... Optional additional fields to store in the channel list.
#'
#' @return A list with elements `type`, `value`, and `name`, plus any additional
#'   fields passed via `...`.
#'
#' @examples
#' set.seed(1)
#' adat <- AlignmentData(list(
#'   s1 = matrix(rnorm(12), 3, 4),
#'   s2 = matrix(rnorm(12), 3, 4)
#' ))
#'
#' g <- list(
#'   s1 = list(roi = guidance_channel("projector", diag(3), name = "roi")),
#'   s2 = list(roi = guidance_channel("projector", diag(3), name = "roi"))
#' )
#'
#' adat2 <- set_guidance(adat, g)
#' get_guidance(adat2, subject = "s1")
#' get_guidance(adat2, type = "projector")
#' @export
guidance_channel <- function(type, value, name = NULL, ...) {
  if (!.is_nonempty_string(type)) {
    stop("'type' must be a non-empty character scalar", call. = FALSE)
  }
  if (!is.null(name) && !.is_nonempty_string(name)) {
    stop("'name' must be NULL or a non-empty character scalar", call. = FALSE)
  }
  x <- list(type = type, value = value, name = name, ...)
  .as_guidance_channel(x, name = name)
}
