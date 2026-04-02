test_that("manifoldalign latent aligners can return reference-space low-rank transforms", {
  skip_if_not_installed("manifoldalign")

  with_temp_registry(neuralign:::.register_gpca, {
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
      target_space = "reference",
      compute_quality = FALSE
    )

    model <- neuralign::get_model(res)
    t_s1 <- neuralign::get_transform(model, "s1")
    t_s2 <- neuralign::get_transform(model, "s2")

    expect_true(inherits(t_s1, "neuralign_low_rank_transform"))
    expect_true(inherits(t_s2, "neuralign_low_rank_transform"))
    expect_equal(nrow(t_s1$U), p)
    expect_equal(nrow(t_s1$V), p)

    aligned <- neuralign::get_aligned(res)
    expect_equal(dim(aligned$s1), c(p, n))
    expect_equal(dim(aligned$s2), c(p, n))

    # New-subject path uses apply_fn; ensure it returns a low-rank operator.
    X3 <- matrix(rnorm(p * n), p, n)
    new_data <- neuralign::AlignmentData(list(s3 = X3), obs_labels = labs)
    res2 <- neuralign::apply_alignment(res, new_data, warn_leakage = FALSE)
    aligned2 <- neuralign::get_aligned(res2)
    expect_equal(dim(aligned2$s3), c(p, n))
    t_s3 <- neuralign::get_transform(neuralign::get_model(res2), "s3")
    expect_true(inherits(t_s3, "neuralign_low_rank_transform"))
  })
})

