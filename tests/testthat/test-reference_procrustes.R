test_that("select_reference(distance='procrustes') returns valid subject", {
  set.seed(10)
  d <- 6
  n <- 20
  base <- matrix(rnorm(d * n), d, n)
  data_list <- list(
    s1 = base,
    s2 = base + matrix(1e-3, d, n),
    s3 = base + matrix(5, d, n)
  )
  adat <- AlignmentData(data_list)

  ref <- select_reference(adat, method = "medoid", distance = "procrustes")
  expect_true(ref %in% names(data_list))
  expect_equal(ref, "s1")
})

test_that("select_reference uses obs_labels overlap for pairwise distances", {
  set.seed(11)
  d <- 5
  labs1 <- paste0("img-", 1:10)
  labs2 <- paste0("img-", 6:15)

  X1 <- matrix(rnorm(d * length(labs1)), d, length(labs1))
  X2 <- X1[, 6:10, drop = FALSE]
  X2 <- cbind(X2, matrix(rnorm(d * 5), d, 5))

  adat <- AlignmentData(
    data = list(s1 = X1, s2 = X2),
    obs_labels = list(s1 = labs1, s2 = labs2)
  )

  ref <- select_reference(adat, method = "medoid", distance = "procrustes")
  expect_true(ref %in% c("s1", "s2"))
})

test_that("centroid reference selection errors when obs_labels differ", {
  set.seed(12)
  d <- 4
  X1 <- matrix(rnorm(d * 10), d, 10)
  X2 <- matrix(rnorm(d * 8), d, 8)
  adat <- AlignmentData(
    data = list(a = X1, b = X2),
    obs_labels = list(a = paste0("t", 1:10), b = paste0("t", 1:8))
  )

  expect_error(select_reference(adat, method = "centroid", distance = "procrustes"))
})

test_that("Procrustes reference distance honors the fitted operator class", {
  set.seed(2)
  d <- 3L
  n <- 10L
  A <- matrix(rnorm(d * n), d, n)
  H <- diag(c(-1, 1, 1))
  B <- H %*% A + matrix(rnorm(d * n, sd = 0.05), d, n)
  C <- matrix(rnorm(d * n), d, n)
  adat <- AlignmentData(list(A = A, B = B, C = C))

  reflected <- procrustes_distance(A, H %*% A, reflection = TRUE)
  proper <- procrustes_distance(A, H %*% A, reflection = FALSE)
  expect_lt(reflected, 1e-10)
  expect_gt(proper, 1)

  ref_o <- select_reference(
    adat,
    method = "medoid",
    distance = "procrustes",
    distance_args = list(reflection = TRUE)
  )
  ref_so <- select_reference(
    adat,
    method = "medoid",
    distance = "procrustes",
    distance_args = list(reflection = FALSE)
  )
  expect_identical(ref_o, "A")
  expect_identical(ref_so, "B")

  fit_o <- fit_alignment(
    adat,
    method = "procrustes",
    reference = "medoid",
    reflection = TRUE,
    compute_quality = FALSE
  )
  fit_so <- fit_alignment(
    adat,
    method = "procrustes",
    reference = "medoid",
    reflection = FALSE,
    compute_quality = FALSE
  )
  expect_identical(get_reference_spec(get_model(fit_o)), ref_o)
  expect_identical(get_reference_spec(get_model(fit_so)), ref_so)
})

test_that("built-in GPA escapes a cancelling arithmetic centroid", {
  set.seed(20260822)
  X <- matrix(rnorm(24), nrow = 3L)

  fit <- neuralign:::.gpa_builtin(
    list(a = X, b = -X),
    scale = FALSE,
    reflection = TRUE,
    rank = NULL,
    tol = 1e-10,
    max_iter = 100L
  )

  aligned <- Map(
    function(Q, Z) Q %*% Z,
    fit$transforms,
    list(a = X, b = -X)
  )
  expect_true(all(is.finite(fit$reference_data)))
  expect_identical(fit$diagnostics$initialization, "procrustes_medoid")
  expect_identical(fit$diagnostics$numerical_status, "converged")
  expect_lt(fit$diagnostics$objective, 1e-10)
  expect_equal(aligned$a, aligned$b, tolerance = 1e-10)
})
