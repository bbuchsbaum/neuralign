#' multivarious interop for neuralign hyperdesign objects
#'
#' manifoldalign calls `multivarious::init_transform(hyperdesign, preproc)` in a
#' few methods. Rather than depending on multidesign's internal multidesign
#' constructors (which may vary across versions), neuralign constructs a
#' lightweight hyperdesign object (`neuralign_hyperdesign`) and provides an
#' S3 method for preprocessing it.
#'
#' This method is registered dynamically in `.onLoad()` when the dependency
#' stack is present.
#'
#' @name multivarious_hyperdesign
#' @keywords internal
NULL

init_transform.neuralign_hyperdesign <- function(x, X, ...) {
  # x: neuralign_hyperdesign (list of domains, each with an $x matrix)
  # X: pre_processor (e.g., multivarious::center())
  preproc <- X

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

    # Use a fresh preproc per domain to avoid leaking fitted params across blocks.
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
