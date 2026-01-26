test_that("feature_block_diagnostics reports rank/effective rank (right convention)", {
  set.seed(1)
  P <- 10
  K <- 5
  r <- 3
  F <- matrix(rnorm(P * r), P, r) %*% matrix(rnorm(r * K), r, K)

  b <- alignment_feature_block(
    x = F,
    name = "task_beta",
    weight = 1,
    feature_names = paste0("p", seq_len(P))
  )

  diag <- feature_block_diagnostics(
    blocks = list(task_beta = b),
    convention = "right",
    tol = 1e-12
  )

  expect_s3_class(diag, "feature_block_diagnostics")
  expect_equal(diag$stacked$transform_dim, K)
  expect_equal(diag$stacked$numeric_rank, r)
  expect_true(diag$stacked$effective_rank <= r + 1e-8)
  expect_equal(diag$stacked$fraction_identified_rank, r / K)
})

test_that("feature_block_diagnostics is invariant to weight scaling (rank/effective rank)", {
  set.seed(2)
  P <- 12
  K <- 4
  r <- 2
  F <- matrix(rnorm(P * r), P, r) %*% matrix(rnorm(r * K), r, K)

  b1 <- alignment_feature_block(F, name = "blk", weight = 1)
  b2 <- alignment_feature_block(F, name = "blk", weight = 9)

  d1 <- feature_block_diagnostics(list(blk = b1), convention = "right", tol = 1e-12)
  d2 <- feature_block_diagnostics(list(blk = b2), convention = "right", tol = 1e-12)

  expect_equal(d1$stacked$numeric_rank, d2$stacked$numeric_rank)
  expect_equal(d1$stacked$effective_rank, d2$stacked$effective_rank, tolerance = 1e-10)
})

test_that("feature_block_diagnostics uses convention to define transform_dim", {
  set.seed(3)
  X <- matrix(rnorm(6 * 4), 6, 4)
  b <- alignment_feature_block(X, name = "blk", weight = 1)

  d_left <- feature_block_diagnostics(list(blk = b), convention = "left", tol = 1e-12)
  expect_equal(d_left$stacked$transform_dim, nrow(X))

  d_right <- feature_block_diagnostics(list(blk = b), convention = "right", tol = 1e-12)
  expect_equal(d_right$stacked$transform_dim, ncol(X))
})

test_that("feature_block_diagnostics supports per-subject input", {
  set.seed(4)
  b1 <- alignment_feature_block(matrix(rnorm(12), 3, 4), name = "blk")
  b2 <- alignment_feature_block(matrix(rnorm(12), 3, 4), name = "blk")

  out <- feature_block_diagnostics(
    blocks = list(s1 = list(blk = b1), s2 = list(blk = b2)),
    convention = "right"
  )

  expect_s3_class(out, "feature_block_diagnostics_by_subject")
  expect_true(all(c("s1", "s2") %in% names(out)))
  expect_s3_class(out$s1, "feature_block_diagnostics")
})

