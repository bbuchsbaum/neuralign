test_that("suggest_block_weights equalize_rows scales inversely with n_rows", {
  x_big <- matrix(1, nrow = 100, ncol = 5)
  x_small <- matrix(1, nrow = 10, ncol = 5)

  blocks_by_subject <- list(
    "sub-01" = list(
      big = alignment_feature_block(x_big, name = "big", feature_names = paste0("f", 1:100)),
      small = alignment_feature_block(x_small, name = "small", feature_names = paste0("g", 1:10))
    ),
    "sub-02" = list(
      big = alignment_feature_block(x_big, name = "big", feature_names = paste0("f", 1:100)),
      small = alignment_feature_block(x_small, name = "small", feature_names = paste0("g", 1:10))
    )
  )

  w <- suggest_block_weights(
    blocks_by_subject,
    method = "equalize_rows",
    per = "global",
    normalize = "none"
  )

  expect_true(is.numeric(w))
  expect_equal(unname(w[["small"]] / w[["big"]]), 10, tolerance = 1e-12)
})

test_that("suggest_block_weights equalize_rms downweights larger-scale blocks", {
  x1 <- matrix(1, nrow = 10, ncol = 5)
  x2 <- matrix(2, nrow = 10, ncol = 5)

  blocks_by_subject <- list(
    "sub-01" = list(
      a = alignment_feature_block(x2, name = "a", feature_names = paste0("f", 1:10)),
      b = alignment_feature_block(x1, name = "b", feature_names = paste0("g", 1:10))
    ),
    "sub-02" = list(
      a = alignment_feature_block(x2, name = "a", feature_names = paste0("f", 1:10)),
      b = alignment_feature_block(x1, name = "b", feature_names = paste0("g", 1:10))
    )
  )

  w <- suggest_block_weights(
    blocks_by_subject,
    method = "equalize_rms",
    per = "global",
    normalize = "none"
  )

  # mean(x^2) is 4 vs 1, so weight ratio is 1/4.
  expect_equal(unname(w[["a"]] / w[["b"]]), 0.25, tolerance = 1e-12)
})

test_that("suggest_block_weights per_subject returns subject-specific vectors", {
  x1 <- matrix(1, nrow = 10, ncol = 5)
  x2 <- matrix(2, nrow = 10, ncol = 5)

  blocks_by_subject <- list(
    "sub-01" = list(
      a = alignment_feature_block(x2, name = "a", feature_names = paste0("f", 1:10))
    ),
    "sub-02" = list(
      a = alignment_feature_block(x1, name = "a", feature_names = paste0("f", 1:10))
    )
  )

  w <- suggest_block_weights(
    blocks_by_subject,
    method = "equalize_rms",
    per = "per_subject",
    normalize = "none"
  )

  expect_true(is.list(w))
  expect_true(all(c("sub-01", "sub-02") %in% names(w)))
  expect_true(all(vapply(w, is.numeric, logical(1))))
  expect_equal(unname(w[["sub-01"]][["a"]] / w[["sub-02"]][["a"]]), 0.25, tolerance = 1e-12)
})

test_that("build_alignment_features accepts per-subject block_weights and suggest_weights", {
  x_a <- matrix(1, nrow = 2, ncol = 3)
  x_b <- matrix(1, nrow = 2, ncol = 3)

  blocks_by_subject <- list(
    "sub-01" = list(
      a = alignment_feature_block(x_a, name = "a", feature_names = c("f1", "f2")),
      b = alignment_feature_block(x_b, name = "b", feature_names = c("g1", "g2"))
    ),
    "sub-02" = list(
      a = alignment_feature_block(x_a, name = "a", feature_names = c("f1", "f2")),
      b = alignment_feature_block(x_b, name = "b", feature_names = c("g1", "g2"))
    )
  )

  weights_by_subject <- list(
    "sub-01" = c(a = 4),
    "sub-02" = c(a = 1)
  )

  built <- build_alignment_features(
    blocks_by_subject,
    block_weights = weights_by_subject,
    check_independence = FALSE,
    check_identifiability = FALSE
  )
  m1 <- built$matrices[["sub-01"]]
  m2 <- built$matrices[["sub-02"]]

  # Block a is first (2 rows) and gets sqrt(4)=2 scaling for sub-01 only.
  expect_equal(unname(m1[1, 1]), 2, tolerance = 1e-12)
  expect_equal(unname(m2[1, 1]), 1, tolerance = 1e-12)
  expect_equal(unname(m1[3, 1]), 1, tolerance = 1e-12) # block b unaffected

  auto <- build_alignment_features(
    blocks_by_subject,
    suggest_weights = "equalize_rows",
    check_independence = FALSE,
    check_identifiability = FALSE
  )
  manual_weights <- suggest_block_weights(blocks_by_subject, method = "equalize_rows")
  manual <- build_alignment_features(
    blocks_by_subject,
    block_weights = manual_weights,
    check_independence = FALSE,
    check_identifiability = FALSE
  )

  expect_equal(auto$matrices, manual$matrices, tolerance = 1e-12)
})

