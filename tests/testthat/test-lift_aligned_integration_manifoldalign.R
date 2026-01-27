test_that("lift_aligned + gpca: lifting to reference preserves latent coordinates", {
  skip_if_not_installed("manifoldalign")
  skip_if_not_installed("multidesign")

  # Some tests clear the global aligner registry; snapshot and restore to
  # avoid ordering-dependent failures.
  registry_env <- neuralign:::.aligner_registry
  clear_registry <- neuralign:::.clear_registry
  register_gpca <- neuralign:::.register_gpca

  old_registry <- as.list(registry_env)
  on.exit({
    clear_registry()
    for (nm in names(old_registry)) {
      registry_env[[nm]] <- old_registry[[nm]]
    }
  }, add = TRUE)

  clear_registry()
  register_gpca()

  set.seed(1)
  p <- 12
  n <- 20
  k <- 3

  X1 <- matrix(rnorm(p * n), p, n)
  X2 <- matrix(rnorm(p * n), p, n)
  labs <- factor(rep(c("A", "B"), length.out = n))

  data <- neuralign::AlignmentData(
    list(s1 = X1, s2 = X2),
    obs_labels = labs
  )

  res <- neuralign::fit_alignment(
    data,
    method = "gpca",
    reference = "s1",
    ncomp = k,
    compute_quality = FALSE
  )

  model <- neuralign::get_model(res)
  aligned <- neuralign::get_aligned(res)
  A_ref <- neuralign::get_transform(model, "s1")

  # If the reference operator is rank-deficient, pinv lifting only reproduces
  # the projection onto its range; skip the exact identity check.
  sv <- svd(as.matrix(A_ref), nu = 0, nv = 0)
  r <- sum(sv$d > (max(sv$d) * 1e-8))
  if (r < nrow(A_ref)) {
    skip("reference transform is rank deficient; cannot assert exact pinv identity")
  }

  lifted_ref <- neuralign::lift_aligned(res, to = "reference", inverse = "pinv")

  for (subj in names(aligned)) {
    expect_equal(
      A_ref %*% lifted_ref[[subj]],
      aligned[[subj]],
      tolerance = 1e-6
    )
  }
})

test_that("lift_aligned + gpca: lifting into a chosen subject space works", {
  skip_if_not_installed("manifoldalign")
  skip_if_not_installed("multidesign")

  registry_env <- neuralign:::.aligner_registry
  clear_registry <- neuralign:::.clear_registry
  register_gpca <- neuralign:::.register_gpca

  old_registry <- as.list(registry_env)
  on.exit({
    clear_registry()
    for (nm in names(old_registry)) {
      registry_env[[nm]] <- old_registry[[nm]]
    }
  }, add = TRUE)

  clear_registry()
  register_gpca()

  set.seed(2)
  p <- 10
  n <- 18
  k <- 2

  X1 <- matrix(rnorm(p * n), p, n)
  X2 <- matrix(rnorm(p * n), p, n)
  labs <- factor(rep(c("A", "B", "C"), length.out = n))

  data <- neuralign::AlignmentData(
    list(s1 = X1, s2 = X2),
    obs_labels = labs
  )

  res <- neuralign::fit_alignment(
    data,
    method = "gpca",
    reference = "s1",
    ncomp = k,
    compute_quality = FALSE
  )

  model <- neuralign::get_model(res)
  aligned <- neuralign::get_aligned(res)
  A_s2 <- neuralign::get_transform(model, "s2")

  sv <- svd(as.matrix(A_s2), nu = 0, nv = 0)
  r <- sum(sv$d > (max(sv$d) * 1e-8))
  if (r < nrow(A_s2)) {
    skip("target subject transform is rank deficient; cannot assert exact pinv identity")
  }

  lifted_to_s2 <- neuralign::lift_aligned(res, to = "subject", subject = "s2", inverse = "pinv")
  for (subj in names(aligned)) {
    expect_equal(
      A_s2 %*% lifted_to_s2[[subj]],
      aligned[[subj]],
      tolerance = 1e-6
    )
  }
})
