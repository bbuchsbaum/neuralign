#' Package Load and Unload Hooks
#'
#' @name zzz
NULL


.onLoad <- function(libname, pkgname) {
  # Register built-in aligners
  .register_procrustes()

  # Register OT aligners if manifoldalign is available
  if (requireNamespace("manifoldalign", quietly = TRUE)) {
    tryCatch(.register_gw(), error = function(e) NULL)
    tryCatch(.register_fpgw(), error = function(e) NULL)
  }
}


.onAttach <- function(libname, pkgname) {
  n_aligners <- length(available_aligners())

  packageStartupMessage(sprintf(
    "neuralign v%s: %d alignment method(s) registered",
    utils::packageVersion(pkgname),
    n_aligners
  ))

  # Hint about available methods
  if (n_aligners > 0) {
    methods <- paste(available_aligners(), collapse = ", ")
    packageStartupMessage(sprintf("  Methods: %s", methods))
  }

  # Note about additional methods
  packageStartupMessage("  Use available_aligners(details=TRUE) for more info")
}


.onUnload <- function(libpath) {
  # Clean up registry
  .clear_registry()
}
