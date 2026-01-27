#' multivarious interop for manifoldalign-style hyperdesign objects
#'
#' Some manifoldalign methods call `multivarious::init_transform(hyperdesign, preproc)`.
#' multivarious defines `init_transform(x, X, ...)` as an S3 generic but does not ship
#' a `hyperdesign` method. We provide one so manifoldalign can preprocess
#' lightweight hyperdesign lists constructed in neuralign.
#'
#' @keywords internal
NULL

#' @export
init_transform.hyperdesign <- function(x, X, ...) {
  # x: hyperdesign (list of domains)
  # X: pre_processor (e.g., multivarious::center())
  preproc <- X

  proc_template <- multivarious::prep(preproc)
  pdata <- x
  proclist <- list()

  for (nm in names(pdata)) {
    dom <- pdata[[nm]]
    if (!is.list(dom) || is.null(dom$x)) {
      stop("hyperdesign domains must be lists with an $x matrix", call. = FALSE)
    }
    Xmat <- dom$x
    if (!is.matrix(Xmat) && !inherits(Xmat, "Matrix")) {
      Xmat <- as.matrix(Xmat)
    }

    # IMPORTANT: use a fresh template per domain to avoid state leakage across
    # differing feature dimensions (some pre_processors store fitted params).
    proc_i <- multivarious::prep(preproc)
    transformed <- multivarious::init_transform(proc_i, Xmat)

    dom$x <- transformed
    pdata[[nm]] <- dom

    proc_attr <- attr(transformed, "preproc")
    proclist[[nm]] <- if (is.null(proc_attr)) proc_i else proc_attr
  }

  attr(pdata, "preproc") <- proclist
  pdata
}

