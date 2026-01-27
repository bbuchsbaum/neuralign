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


# ---------- Additional feature block tests ----------

test_that("alignment_feature_block errors on non-matrix input", {
  expect_error(alignment_feature_block("text", name = "blk"), "matrix-like")
})

test_that("alignment_feature_block with NULL meta is ok", {
  b <- alignment_feature_block(matrix(1, 2, 3), name = "blk", meta = NULL)
  expect_null(b$meta)
})

test_that("stack_feature_blocks errors on empty list", {
  expect_error(stack_feature_blocks(list()), "non-empty")
})

test_that("stack_feature_blocks errors on non-block items", {
  expect_error(stack_feature_blocks(list(list(x = 1))), "alignment_feature_block")
})

test_that("stack_feature_blocks errors on unnamed block_weights", {
  b <- alignment_feature_block(matrix(1, 2, 3), name = "b")
  expect_error(stack_feature_blocks(list(b), block_weights = c(1)), "named numeric")
})

test_that("stack_feature_blocks rownames include block prefix", {
  b1 <- alignment_feature_block(matrix(1, 2, 3), name = "alpha", feature_names = c("f1", "f2"))
  b2 <- alignment_feature_block(matrix(2, 1, 3), name = "beta", feature_names = c("g1"))
  stacked <- stack_feature_blocks(list(b1, b2))
  expect_equal(rownames(stacked), c("alpha:f1", "alpha:f2", "beta:g1"))
})

test_that("stack_feature_blocks uses existing rownames when feature_names is NULL", {
  x <- matrix(1, 2, 3)
  rownames(x) <- c("r1", "r2")
  b <- alignment_feature_block(x, name = "blk")
  stacked <- stack_feature_blocks(list(b))
  expect_equal(rownames(stacked), c("blk:r1", "blk:r2"))
})

test_that("stack_feature_blocks with weight=0 zeroes the block", {
  b <- alignment_feature_block(matrix(5, 2, 3), name = "b", weight = 0)
  stacked <- stack_feature_blocks(list(b))
  expect_equal(max(abs(stacked)), 0)
})

test_that("harmonize_feature_blocks errors on unnamed subject list", {
  expect_error(harmonize_feature_blocks(list()), "non-empty")
  b <- alignment_feature_block(matrix(1, 2, 3), name = "a", feature_names = c("f1", "f2"))
  expect_error(harmonize_feature_blocks(list(list(a = b))), "named list")
})

test_that("harmonize_feature_blocks errors on invalid min_features", {
  b <- alignment_feature_block(matrix(1, 2, 3), name = "a", feature_names = c("f1", "f2"))
  expect_error(harmonize_feature_blocks(list(sub1 = list(a = b)), min_features = 0), "positive integer")
})

test_that("harmonize_feature_blocks errors when feature names are missing", {
  b1 <- alignment_feature_block(matrix(1, 2, 3), name = "a")  # no feature_names, no rownames
  b2 <- alignment_feature_block(matrix(1, 2, 3), name = "a")
  expect_error(
    harmonize_feature_blocks(list(sub1 = list(a = b1), sub2 = list(a = b2))),
    "feature names are missing"
  )
})

test_that("feature_block_diagnostics returns per_block and stacked summaries", {
  b1 <- alignment_feature_block(matrix(rnorm(12), 3, 4), name = "func", weight = 1,
                                  feature_names = c("f1", "f2", "f3"))
  b2 <- alignment_feature_block(matrix(rnorm(8), 2, 4), name = "anat", weight = 2,
                                  feature_names = c("g1", "g2"))

  diag <- feature_block_diagnostics(list(b1, b2), convention = "left")
  expect_s3_class(diag, "feature_block_diagnostics")
  expect_equal(nrow(diag$per_block), 2)
  expect_equal(diag$per_block$block, c("func", "anat"))
  expect_true(!is.null(diag$stacked$numeric_rank))
  expect_true(!is.null(diag$stacked$effective_rank))
})

test_that("feature_block_diagnostics with right convention", {
  b1 <- alignment_feature_block(matrix(rnorm(12), 3, 4), name = "func",
                                  feature_names = c("f1", "f2", "f3"))

  diag_l <- feature_block_diagnostics(list(b1), convention = "left")
  diag_r <- feature_block_diagnostics(list(b1), convention = "right")

  # Left convention: transform_dim = nrow = 3
  expect_equal(diag_l$per_block$transform_dim, 3L)
  # Right convention: transform_dim = ncol = 4
  expect_equal(diag_r$per_block$transform_dim, 4L)
})

test_that("feature_block_diagnostics include_singular_values", {
  b <- alignment_feature_block(matrix(rnorm(12), 3, 4), name = "blk",
                                feature_names = c("f1", "f2", "f3"))
  diag_no <- feature_block_diagnostics(list(b), include_singular_values = FALSE)
  expect_null(diag_no$stacked$singular_values)

  diag_yes <- feature_block_diagnostics(list(b), include_singular_values = TRUE)
  expect_true(!is.null(diag_yes$stacked$singular_values))
  expect_equal(length(diag_yes$stacked$singular_values), min(3, 4))
})

test_that("feature_block_diagnostics with per-subject input", {
  b <- alignment_feature_block(matrix(rnorm(12), 3, 4), name = "blk",
                                feature_names = c("f1", "f2", "f3"))
  per_subj <- list(
    s1 = list(b),
    s2 = list(b)
  )
  diag <- feature_block_diagnostics(per_subj, convention = "left")
  expect_s3_class(diag, "feature_block_diagnostics_by_subject")
  expect_equal(names(diag), c("s1", "s2"))
  expect_s3_class(diag$s1, "feature_block_diagnostics")
})

test_that("feature_block_diagnostics errors on empty blocks", {
  expect_error(feature_block_diagnostics(list()), "non-empty")
})

test_that("feature_block_diagnostics errors on bad tol", {
  b <- alignment_feature_block(matrix(1, 2, 3), name = "a", feature_names = c("f1", "f2"))
  expect_error(feature_block_diagnostics(list(b), tol = 0), "\\(0, 1\\)")
  expect_error(feature_block_diagnostics(list(b), tol = 1), "\\(0, 1\\)")
})

test_that("build_alignment_features harmonize=intersection stacks harmonized blocks", {
  s1 <- list(
    a = alignment_feature_block(
      matrix(1:6, 3, 2),
      name = "a",
      feature_names = c("f1", "f2", "f3")
    ),
    b = alignment_feature_block(
      matrix(c(10, 20, 30, 40), 2, 2),
      name = "b",
      feature_names = c("g1", "g2")
    )
  )
  s2 <- list(
    a = alignment_feature_block(
      matrix(101:106, 3, 2),
      name = "a",
      feature_names = c("f2", "f3", "f4")
    ),
    b = alignment_feature_block(
      matrix(c(50, 60, 70, 80), 2, 2),
      name = "b",
      feature_names = c("g1", "g2")
    )
  )

  res <- build_alignment_features(list(sub1 = s1, sub2 = s2), harmonize = "intersection", min_features = 2)
  expect_true(all(c("matrices", "blocks", "per_block", "dropped_blocks", "warnings") %in% names(res)))
  expect_equal(names(res$matrices), c("sub1", "sub2"))
  expect_equal(rownames(res$matrices$sub1), c("a:f2", "a:f3", "b:g1", "b:g2"))
  expect_equal(unname(res$matrices$sub1[1:2, , drop = FALSE]), s1$a$x[2:3, , drop = FALSE])
  expect_equal(unname(res$matrices$sub2[1:2, , drop = FALSE]), s2$a$x[1:2, , drop = FALSE])
})

test_that("build_alignment_features harmonize=union_fill uses union vocabulary with fill=0", {
  s1 <- list(
    a = alignment_feature_block(
      matrix(1:6, 3, 2),
      name = "a",
      feature_names = c("f1", "f2", "f3")
    ),
    b = alignment_feature_block(
      matrix(c(10, 20, 30, 40), 2, 2),
      name = "b",
      feature_names = c("g1", "g2")
    )
  )
  s2 <- list(
    a = alignment_feature_block(
      matrix(101:106, 3, 2),
      name = "a",
      feature_names = c("f2", "f3", "f4")
    ),
    b = alignment_feature_block(
      matrix(c(50, 60, 70, 80), 2, 2),
      name = "b",
      feature_names = c("g1", "g2")
    )
  )

  res <- build_alignment_features(
    list(sub1 = s1, sub2 = s2),
    harmonize = "union_fill",
    fill = 0,
    min_features = 2,
    union_order = "first_seen"
  )
  expect_equal(rownames(res$matrices$sub1), c("a:f1", "a:f2", "a:f3", "a:f4", "b:g1", "b:g2"))
  expect_equal(unname(res$matrices$sub2["a:f1", , drop = FALSE]), matrix(0, 1, 2))
  expect_equal(unname(res$matrices$sub1["a:f4", , drop = FALSE]), matrix(0, 1, 2))
})

test_that("build_alignment_features errors when blocks have incompatible column counts", {
  s1 <- list(
    a = alignment_feature_block(matrix(1, 2, 3), name = "a", feature_names = c("f1", "f2")),
    b = alignment_feature_block(matrix(1, 1, 2), name = "b", feature_names = "g1")
  )
  s2 <- list(
    a = alignment_feature_block(matrix(1, 2, 3), name = "a", feature_names = c("f1", "f2")),
    b = alignment_feature_block(matrix(1, 1, 2), name = "b", feature_names = "g1")
  )

  expect_error(
    build_alignment_features(list(sub1 = s1, sub2 = s2), harmonize = "intersection", min_features = 1),
    "differing column counts"
  )
})

test_that("build_alignment_features warns for requires_independence blocks unless obs_crossfit", {
  s1 <- list(
    a = alignment_feature_block(
      matrix(1:4, 2, 2),
      name = "a",
      feature_names = c("f1", "f2"),
      meta = list(requires_independence = TRUE)
    )
  )
  s2 <- list(
    a = alignment_feature_block(
      matrix(5:8, 2, 2),
      name = "a",
      feature_names = c("f1", "f2"),
      meta = list(requires_independence = TRUE)
    )
  )

  expect_equal(feature_block_requires_independence(list(sub1 = s1, sub2 = s2)), c(a = TRUE))

  expect_warning(
    res <- build_alignment_features(list(sub1 = s1, sub2 = s2), harmonize = "intersection", min_features = 2),
    "meta\\$requires_independence=TRUE: a"
  )
  expect_true(isTRUE(res$per_block$requires_independence[[1L]]))

  expect_warning(
    build_alignment_features(list(sub1 = s1, sub2 = s2),
      harmonize = "intersection",
      min_features = 2,
      obs_crossfit = TRUE
    ),
    NA
  )

  expect_warning(
    build_alignment_features(list(sub1 = s1, sub2 = s2),
      harmonize = "intersection",
      min_features = 2,
      check_independence = FALSE
    ),
    NA
  )
})
