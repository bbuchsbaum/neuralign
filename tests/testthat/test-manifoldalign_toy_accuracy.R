test_that("manifoldalign projection methods recover shared latent structure on toy data", {
  skip_if_not_installed("manifoldalign")

  set.seed(1)

  make_latent_toy <- function(subject_ids = c("s1", "s2"),
                              n_features = 40L,
                              n_obs = 40L,
                              k = 5L,
                              noise_sd = 1e-2) {
    subject_ids <- as.character(subject_ids)
    n_features <- as.integer(n_features)
    n_obs <- as.integer(n_obs)
    k <- as.integer(k)
    noise_sd <- as.numeric(noise_sd)

    Z <- matrix(stats::rnorm(k * n_obs), k, n_obs)
    data_list <- lapply(subject_ids, function(id) {
      # Orthonormal columns: V is (p x k)
      V <- qr.Q(qr(matrix(stats::rnorm(n_features * k), n_features, k)))
      X <- V %*% Z + matrix(stats::rnorm(n_features * n_obs, sd = noise_sd), n_features, n_obs)
      X
    })
    names(data_list) <- subject_ids

    list(
      data = data_list,
      obs_labels = make_test_obs_labels(n_obs, prefix = "obs")
    )
  }

  toy <- make_latent_toy()
  adat <- AlignmentData(toy$data, obs_labels = toy$obs_labels)

  cases <- list(
    gpca = list(register = neuralign:::.register_gpca, args = list(ncomp = 5L), min_abs_cor = 0.95),
    lowrank = list(register = neuralign:::.register_lowrank, args = list(ncomp = 5L), min_abs_cor = 0.95),
    coupled_diagonalization = list(
      register = neuralign:::.register_coupled_diag,
      args = list(ncomp = 5L, ncomp_per_domain = 20L),
      min_abs_cor = 0.90
    )
  )

  with_temp_registry(lapply(cases, `[[`, "register"), {
    for (nm in names(cases)) {
      spec <- cases[[nm]]
      res <- do.call(
        fit_alignment,
        c(
          list(
            data = adat,
            method = nm,
            reference = "s1",
            compute_quality = FALSE
          ),
          spec$args
        )
      )
      aligned <- get_aligned(res)
      expect_true(is.list(aligned) && length(aligned) == 2L, info = nm)
      expect_equal(dim(aligned$s1), dim(aligned$s2), info = nm)

      abs_cor <- abs(stats::cor(as.vector(aligned$s1), as.vector(aligned$s2)))
      expect_true(
        abs_cor >= spec$min_abs_cor,
        info = sprintf("%s: abs cor %.3f < %.3f", nm, abs_cor, spec$min_abs_cor)
      )
    }
  })
})


test_that("KEMA training embeddings recover labelled correspondences above chance", {
  skip_if_not_installed("manifoldalign")

  set.seed(1)
  p <- 40L
  n <- 40L
  k <- 5L
  Z <- matrix(stats::rnorm(k * n), k, n)
  xs <- lapply(seq_len(2L), function(i) {
    V <- qr.Q(qr(matrix(stats::rnorm(p * k), p, k)))
    V %*% Z + matrix(stats::rnorm(p * n, sd = 1e-2), p, n)
  })
  names(xs) <- c("s1", "s2")
  data <- AlignmentData(xs, obs_labels = sprintf("obs-%02d", seq_len(n)))

  with_temp_registry(neuralign:::.register_kema, {
    res <- fit_alignment(
      data,
      method = "kema",
      reference = "s1",
      ncomp = k,
      compute_quality = FALSE
    )
    aligned <- get_aligned(res)
    a <- t(aligned$s1)
    b <- t(aligned$s2)
    distances <- as.matrix(stats::dist(rbind(a, b)))[
      seq_len(n), n + seq_len(n), drop = FALSE
    ]
    predicted <- max.col(-distances, ties.method = "first")
    accuracy <- mean(predicted == seq_len(n))

    expect_true(
      accuracy >= 0.20,
      info = sprintf("KEMA labelled-correspondence accuracy %.3f < 0.200", accuracy)
    )
  })
})


test_that("GRASP graph accuracy meets the mandatory permutation contract", {
  skip_if_not_installed("manifoldalign")

  set.seed(1)

  # Make distinctive per-node signatures (rows = nodes/features, cols = observations)
  # so that a permutation should be recoverable.
  p <- 60L
  n_obs <- 40L
  X1 <- matrix(0, p, n_obs)
  for (i in seq_len(p)) {
    X1[i, ] <- sin(seq(0, 2 * pi, length.out = n_obs) * (i / 7)) + (i / p)
  }
  perm <- sample.int(p)
  X2 <- X1[perm, , drop = FALSE]

  adat <- AlignmentData(list(s1 = X1, s2 = X2))

  # Helper: compute recovered mapping accuracy vs true permutation.
  operator_to_perm <- function(A) {
    # A maps source->target rows; take argmax per row.
    apply(as.matrix(A), 1L, which.max)
  }
  perm_inv <- match(seq_len(p), perm)

  with_temp_registry(neuralign:::.register_grasp, {
    res <- fit_alignment(adat, method = "grasp", reference = "s1", compute_quality = FALSE)
    A <- get_transform(get_model(res), "s2")
    pred <- operator_to_perm(A)
    acc <- mean(pred == perm_inv)
    expect_true(acc >= 0.5, info = sprintf("grasp: perm accuracy %.3f", acc))
  })
})


test_that("CONE remains fail-closed after failing the graph oracle", {
  skip_if_not_installed("manifoldalign")

  with_temp_registry(NULL, {
    expect_false(is_aligner_registered("cone"))
    expect_false(is_aligner_registered("cone_align"))
    expect_error(neuralign:::.register_cone(), "accuracy contract")
  })
})
