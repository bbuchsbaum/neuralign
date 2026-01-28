test_that("projection manifoldalign adapters fit and apply to new subjects", {
  skip_if_not_installed("manifoldalign")

  set.seed(1)
  p <- 12
  n <- 16
  k <- 3

  X1 <- matrix(rnorm(p * n), p, n)
  X2 <- matrix(rnorm(p * n), p, n)
  X3 <- matrix(rnorm(p * n), p, n)
  labs <- factor(rep(c("A", "B"), length.out = n))

  data_train <- neuralign::AlignmentData(list(s1 = X1, s2 = X2), obs_labels = labs)
  data_new <- neuralign::AlignmentData(list(s3 = X3), obs_labels = labs)

  cases <- list(
    kema = list(register = neuralign:::.register_kema, args = list(ncomp = k, knn = 5L)),
    coupled_diag = list(
      register = neuralign:::.register_coupled_diag,
      args = list(ncomp = k, ncomp_per_domain = 6L, max_iter = 25L)
    ),
    lowrank = list(register = neuralign:::.register_lowrank, args = list(ncomp = k))
  )

  for (method in names(cases)) {
    case <- cases[[method]]

    with_temp_registry(list(case$register), {
      res <- do.call(
        neuralign::fit_alignment,
        c(list(data = data_train, method = method, reference = "s1", compute_quality = FALSE), case$args)
      )

      applied <- neuralign::apply_alignment(res, data_new)
      model <- neuralign::get_model(res)
      applied_model <- neuralign::get_model(applied)

      A_new <- neuralign::get_transform(applied_model, "s3")
      expect_equal(dim(A_new), c(k, p))
      expect_equal(dim(neuralign::get_aligned(applied)[["s3"]]), c(k, n))

      # apply_alignment should use the aligner's apply_fn when present; compare
      # against a direct call with the same fit_result payload.
      aligner <- neuralign::get_aligner(method)
      expect_true(is.function(aligner$apply_fn))

      fit_payload <- list(
        transforms = neuralign::get_transforms(model),
        reference_data = neuralign::get_reference(model),
        space_from = model@space_from,
        space_to = model@space_to,
        method_state = model@method_state
      )
      direct <- aligner$apply_fn(fit_result = fit_payload, new_data = data_new)
      A_direct <- direct$transforms[[1L]]

      expect_equal(A_new, A_direct, tolerance = 1e-10)
    })
  }
})

test_that("graph manifoldalign adapters return sparse assignment operators and can apply to new subjects", {
  skip_if_not_installed("manifoldalign")

  set.seed(2)
  p <- 14
  n <- 6

  X1 <- matrix(rnorm(p * n), p, n)
  X2 <- matrix(rnorm(p * n), p, n)
  X3 <- matrix(rnorm(p * n), p, n)

  data_train <- neuralign::AlignmentData(list(s1 = X1, s2 = X2))
  data_new <- neuralign::AlignmentData(list(s3 = X3))

  cases <- list(
    grasp = neuralign:::.register_grasp,
    cone = neuralign:::.register_cone
  )

  for (method in names(cases)) {
    register_fn <- cases[[method]]

    with_temp_registry(list(register_fn), {
      res <- neuralign::fit_alignment(
        data_train,
        method = method,
        reference = "s1",
        compute_quality = FALSE
      )

      model <- neuralign::get_model(res)
      A_ref <- neuralign::get_transform(model, "s1")
      A_s2 <- neuralign::get_transform(model, "s2")

      expect_true(inherits(A_ref, "Matrix"))
      expect_equal(dim(A_ref), c(p, p))

      expect_true(inherits(A_s2, "Matrix"))
      expect_true(methods::is(A_s2, "sparseMatrix"))
      expect_equal(dim(A_s2), c(p, p))
      expect_true(all(Matrix::rowSums(A_s2) <= 1 + 1e-12))

      aligned <- neuralign::get_aligned(res)
      expect_equal(dim(aligned[["s1"]]), c(p, n))
      expect_equal(dim(aligned[["s2"]]), c(p, n))

      applied <- neuralign::apply_alignment(res, data_new)
      A_new <- neuralign::get_transform(neuralign::get_model(applied), "s3")
      expect_true(inherits(A_new, "Matrix"))
      expect_true(methods::is(A_new, "sparseMatrix"))
      expect_equal(dim(A_new), c(p, p))
      expect_true(all(Matrix::rowSums(A_new) <= 1 + 1e-12))

      aligned_new <- neuralign::get_aligned(applied)[["s3"]]
      expect_equal(dim(aligned_new), c(p, n))
    })
  }
})
