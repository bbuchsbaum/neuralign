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


# ---------- Additional intrinsic geometry tests ----------

test_that("laplacian_eigenmodes unnormalized Laplacian works", {
  A <- .path_adjacency(8)
  Z <- laplacian_eigenmodes(A, k = 3, normalized = FALSE, method = "eigen")
  expect_equal(dim(Z), c(8, 3))
  expect_false(anyNA(Z))
})

test_that("laplacian_eigenmodes with return_values gives eigenvalues", {
  A <- .path_adjacency(6)
  res <- laplacian_eigenmodes(A, k = 2, method = "eigen", return_values = TRUE)
  expect_true(is.list(res))
  expect_equal(dim(res$vectors), c(6, 2))
  expect_equal(length(res$values), 2)
  # Eigenvalues should be positive (non-trivial modes)
  expect_true(all(res$values > 0))
})

test_that("laplacian_eigenmodes errors on invalid k", {
  A <- .path_adjacency(5)
  expect_error(laplacian_eigenmodes(A, k = 0), "positive integer")
  expect_error(laplacian_eigenmodes(A, k = -1), "positive integer")
})

test_that("laplacian_eigenmodes errors on negative tol", {
  A <- .path_adjacency(5)
  expect_error(laplacian_eigenmodes(A, k = 2, tol = -0.1), "non-negative")
})

test_that("laplacian_eigenmodes errors on non-square adjacency", {
  expect_error(laplacian_eigenmodes(matrix(1, 3, 4), k = 2), "square")
})

test_that("laplacian_eigenmodes with drop_zero=FALSE includes trivial modes", {
  A <- .path_adjacency(6)
  Z_no_drop <- laplacian_eigenmodes(A, k = 2, drop_zero = FALSE, method = "eigen")
  expect_equal(dim(Z_no_drop), c(6, 2))
})

test_that("laplacian_eigenmodes errors when k exceeds available modes", {
  A <- .path_adjacency(4)
  # Path of 4 nodes has 3 non-trivial modes only
  expect_error(laplacian_eigenmodes(A, k = 4, method = "eigen"), "only.*available")
})

test_that("laplacian_eigenmodes errors on negative degrees", {
  A <- diag(3)
  A[1, 2] <- -5
  A[2, 1] <- -5
  expect_error(laplacian_eigenmodes(A, k = 1), "negative degrees")
})

test_that("set_intrinsic_geometry_guidance errors on missing subjects", {
  data <- AlignmentData(list(
    s1 = matrix(rnorm(20), 5, 4),
    s2 = matrix(rnorm(24), 6, 4)
  ))

  adj <- list(s1 = .path_adjacency(5))  # missing s2
  expect_error(
    set_intrinsic_geometry_guidance(data, adjacency_by_subject = adj, k = 2, method = "eigen"),
    "missing subjects"
  )
})

test_that("set_intrinsic_geometry_guidance errors on bad channel_name", {
  data <- AlignmentData(list(s1 = matrix(1, 3, 2)))
  adj <- list(s1 = .path_adjacency(3))
  expect_error(
    set_intrinsic_geometry_guidance(data, adjacency_by_subject = adj, channel_name = ""),
    "non-empty"
  )
})

test_that("set_intrinsic_geometry_guidance errors on non-AlignmentData", {
  expect_error(
    set_intrinsic_geometry_guidance("not_adat", list()),
    "AlignmentData"
  )
})

test_that("set_intrinsic_geometry_guidance errors on unnamed adjacency list", {
  data <- AlignmentData(list(s1 = matrix(1, 3, 2)))
  expect_error(
    set_intrinsic_geometry_guidance(data, adjacency_by_subject = list(.path_adjacency(3))),
    "named list"
  )
})

