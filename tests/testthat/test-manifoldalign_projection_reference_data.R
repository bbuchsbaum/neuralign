test_that("manifoldalign projection aligners store latent reference_data (k x n)", {
  skip_if_not_installed("manifoldalign")

  registry_env <- neuralign:::.aligner_registry
  clear_registry <- neuralign:::.clear_registry

  with_temp_registry <- function(register_fn, code) {
    old_registry <- as.list(registry_env)
    on.exit({
      clear_registry()
      for (nm in names(old_registry)) {
        registry_env[[nm]] <- old_registry[[nm]]
      }
    }, add = TRUE)

    clear_registry()
    register_fn()
    force(code)
  }

  set.seed(123)
  p <- 12
  n <- 20
  k <- 3

  X1 <- matrix(rnorm(p * n), p, n)
  X2 <- matrix(rnorm(p * n), p, n)
  labs <- factor(rep(c("A", "B"), length.out = n))

  data <- neuralign::AlignmentData(list(s1 = X1, s2 = X2), obs_labels = labs)

  cases <- list(
    kema = list(register = neuralign:::.register_kema, args = list(ncomp = k, knn = 5L)),
    gpca = list(register = neuralign:::.register_gpca, args = list(ncomp = k)),
    coupled_diag = list(
      register = neuralign:::.register_coupled_diag,
      args = list(ncomp = k, ncomp_per_domain = 5L, max_iter = 25L)
    ),
    lowrank = list(register = neuralign:::.register_lowrank, args = list(ncomp = k))
  )

  for (nm in names(cases)) {
    case <- cases[[nm]]
    with_temp_registry(case$register, {
      res <- do.call(
        neuralign::fit_alignment,
        c(list(data = data, method = nm, reference = "s1", compute_quality = FALSE), case$args)
      )

      model <- neuralign::get_model(res)
      aligned <- neuralign::get_aligned(res)
      ref_latent <- neuralign::get_reference(model)

      expect_equal(dim(ref_latent), dim(aligned[["s1"]]))
      expect_equal(ref_latent, aligned[["s1"]], tolerance = 1e-6)

      # reconstruction metrics should operate in the latent target space
      q <- neuralign::alignment_quality(
        res,
        metrics = c("reconstruction"),
        reference = ref_latent
      )
      expect_true("mean_reference_correlation" %in% names(q))
      expect_true(is.numeric(q$mean_reference_correlation))
    })
  }
})

