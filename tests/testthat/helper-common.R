with_temp_registry <- function(register_fns = NULL, code) {
  registry_env <- neuralign:::.aligner_registry
  clear_registry <- neuralign:::.clear_registry

  old_registry <- as.list(registry_env)
  on.exit({
    clear_registry()
    for (nm in names(old_registry)) {
      registry_env[[nm]] <- old_registry[[nm]]
    }
  }, add = TRUE)

  clear_registry()
  if (!is.null(register_fns)) {
    if (is.function(register_fns)) {
      register_fns <- list(register_fns)
    }
    if (!is.list(register_fns)) {
      stop("'register_fns' must be a function, a list of functions, or NULL", call. = FALSE)
    }
    if (!all(vapply(register_fns, is.function, logical(1)))) {
      stop("'register_fns' must contain only functions", call. = FALSE)
    }
    for (fn in register_fns) fn()
  }

  force(code)
}

