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

