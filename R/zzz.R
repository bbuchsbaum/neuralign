#' Package Load and Unload Hooks
#'
#' @name zzz
NULL


.onLoad <- function(libname, pkgname) {
  # Register built-in aligners
  .register_procrustes()
  .register_procrustes_graph()
  .register_kprocrustes()

  # Register manifoldalign-backed aligners when available
  if (requireNamespace("manifoldalign", quietly = TRUE)) {
    tryCatch(.register_gw(), error = function(e) NULL)
    tryCatch(.register_fpgw(), error = function(e) NULL)
    tryCatch(.register_kema(), error = function(e) NULL)
    tryCatch(.register_coupled_diag(), error = function(e) NULL)
    tryCatch(.register_gpca(), error = function(e) NULL)
    tryCatch(.register_mma(), error = function(e) NULL)
    tryCatch(.register_grasp(), error = function(e) NULL)
    tryCatch(.register_cone(), error = function(e) NULL)
    tryCatch(.register_parrot(), error = function(e) NULL)
    tryCatch(.register_lowrank(), error = function(e) NULL)
    tryCatch(.register_lse(), error = function(e) NULL)
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
