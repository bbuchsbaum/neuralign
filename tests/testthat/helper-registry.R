ensure_test_aligner <- function(method) {
  method <- match.arg(method, c("procrustes", "procrustes_graph", "kprocrustes", "gw", "fpgw"))
  if (is_aligner_registered(method)) {
    return(invisible(TRUE))
  }

  switch(method,
    procrustes = neuralign:::.register_procrustes(),
    procrustes_graph = neuralign:::.register_procrustes_graph(),
    kprocrustes = neuralign:::.register_kprocrustes(),
    gw = neuralign:::.register_gw(),
    fpgw = neuralign:::.register_fpgw()
  )

  invisible(TRUE)
}

