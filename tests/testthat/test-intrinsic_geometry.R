.path_adjacency <- function(n) {
  A <- matrix(0, n, n)
  for (i in seq_len(n - 1L)) {
    A[i, i + 1L] <- 1
    A[i + 1L, i] <- 1
  }
  A
}

test_that("laplacian_eigenmodes returns expected dimensions and stable signs", {
  A <- .path_adjacency(6)
  Z1 <- laplacian_eigenmodes(A, k = 2, normalized = TRUE, method = "eigen")
  Z2 <- laplacian_eigenmodes(A, k = 2, normalized = TRUE, method = "eigen")

  expect_equal(dim(Z1), c(6, 2))
  expect_false(anyNA(Z1))
  expect_equal(Z1, Z2)
})

test_that("laplacian_eigenmodes errors on non-symmetric adjacency", {
  A <- .path_adjacency(5)
  A[1, 2] <- 0
  expect_error(laplacian_eigenmodes(A, k = 2), "symmetric")
})

test_that("set_intrinsic_geometry_guidance attaches coords channels", {
  set.seed(1)
  data <- AlignmentData(list(
    s1 = matrix(rnorm(20), 5, 4),
    s2 = matrix(rnorm(24), 6, 4)
  ))

  adj <- list(
    s1 = .path_adjacency(5),
    s2 = .path_adjacency(6)
  )

  data2 <- set_intrinsic_geometry_guidance(
    data,
    adjacency_by_subject = adj,
    channel_name = "geom",
    k = 2,
    method = "eigen",
    validate = TRUE
  )

  g <- get_guidance(data2, type = "coords")
  expect_true("geom" %in% names(g$s1))
  expect_true("geom" %in% names(g$s2))
  expect_equal(dim(g$s1$geom$value), c(5, 2))
  expect_equal(dim(g$s2$geom$value), c(6, 2))
})

