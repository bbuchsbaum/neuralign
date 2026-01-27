test_that("alignment_feature_block validates inputs", {
  x <- matrix(1, 3, 4)
  b <- alignment_feature_block(
    x,
    name = "blk",
    weight = 2,
    feature_names = c("a", "b", "c"),
    meta = list(source_type = "supervised")
  )
  expect_s3_class(b, "alignment_feature_block")
  expect_equal(b$name, "blk")
  expect_equal(b$weight, 2)
  expect_equal(b$feature_names, c("a", "b", "c"))
  expect_equal(b$meta$source_type, "supervised")
  expect_error(alignment_feature_block(x, name = "", weight = 1))
  expect_error(alignment_feature_block(x, name = "blk", weight = -1))
  expect_error(alignment_feature_block(x, name = "blk", feature_names = c("a", "b")))
  expect_error(alignment_feature_block(x, name = "blk", meta = 1))
})

test_that("stack_feature_blocks applies sqrt(weight) scaling and rbinds", {
  x1 <- matrix(1, 2, 3)
  x2 <- matrix(2, 1, 3)
  b1 <- alignment_feature_block(x1, name = "b1", weight = 4, feature_names = c("f1", "f2"))
  b2 <- alignment_feature_block(x2, name = "b2", weight = 1, feature_names = c("g1"))

  stacked <- stack_feature_blocks(list(b1, b2))
  expect_equal(dim(stacked), c(3, 3))
  expect_equal(unname(stacked[1:2, , drop = FALSE]), 2 * x1)
  expect_equal(unname(stacked[3, , drop = FALSE]), x2)
})

test_that("stack_feature_blocks supports block_weights overrides", {
  x1 <- matrix(1, 2, 2)
  b1 <- alignment_feature_block(x1, name = "b1", weight = 1, feature_names = c("a", "b"))
  stacked <- stack_feature_blocks(list(b1), block_weights = c(b1 = 9))
  expect_equal(unname(stacked), 3 * x1)
})

test_that("harmonize_feature_blocks intersects feature names across subjects", {
  s1 <- list(
    a = alignment_feature_block(matrix(1, 3, 2), name = "a", feature_names = c("f1", "f2", "f3")),
    b = alignment_feature_block(matrix(1, 2, 2), name = "b", feature_names = c("g1", "g2"))
  )
  s2 <- list(
    a = alignment_feature_block(matrix(2, 3, 2), name = "a", feature_names = c("f2", "f3", "f4")),
    b = alignment_feature_block(matrix(2, 2, 2), name = "b", feature_names = c("g1", "g2"))
  )

  out <- harmonize_feature_blocks(list(sub1 = s1, sub2 = s2), min_features = 2)
  expect_equal(out$sub1$a$feature_names, c("f2", "f3"))
  expect_equal(dim(out$sub1$a$x), c(2, 2))
  expect_equal(out$sub2$a$feature_names, c("f2", "f3"))
})

test_that("harmonize_feature_blocks preserves feature block meta", {
  s1 <- list(
    a = alignment_feature_block(
      matrix(1, 3, 2),
      name = "a",
      feature_names = c("f1", "f2", "f3"),
      meta = list(source_type = "supervised", requires_independence = TRUE)
    )
  )
  s2 <- list(
    a = alignment_feature_block(
      matrix(2, 3, 2),
      name = "a",
      feature_names = c("f2", "f3", "f4"),
      meta = list(source_type = "supervised", requires_independence = TRUE)
    )
  )

  out <- harmonize_feature_blocks(list(sub1 = s1, sub2 = s2), min_features = 2)
  expect_equal(out$sub1$a$meta$source_type, "supervised")
  expect_true(isTRUE(out$sub1$a$meta$requires_independence))
})

test_that("harmonize_feature_blocks drops missing/low-overlap blocks", {
  s1 <- list(
    a = alignment_feature_block(matrix(1, 2, 2), name = "a", feature_names = c("f1", "f2")),
    c = alignment_feature_block(matrix(1, 2, 2), name = "c", feature_names = c("h1", "h2"))
  )
  s2 <- list(
    a = alignment_feature_block(matrix(2, 2, 2), name = "a", feature_names = c("f2", "f3"))
  )

  out <- NULL
  expect_warning(
    expect_warning(
      out <- harmonize_feature_blocks(list(sub1 = s1, sub2 = s2), min_features = 2),
      "Dropping blocks not present for all subjects"
    ),
    "Dropping block 'a'"
  )
  expect_true(is.null(out$sub1$c))
  expect_true(is.null(out$sub2$c))
  # Block a only has 1 common feature (f2) -> dropped
  expect_true(is.null(out$sub1$a))
  expect_true(is.null(out$sub2$a))
})
