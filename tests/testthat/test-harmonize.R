test_that("harmonize_union_fill harmonizes rows with partial overlap", {
  x1 <- matrix(1:6, 3, 2, dimnames = list(c("a", "b", "c"), c("v1", "v2")))
  x2 <- matrix(1:4, 2, 2, dimnames = list(c("b", "d"), c("v1", "v2")))

  res <- harmonize_union_fill(list(s1 = x1, s2 = x2), axis = "rows", union_order = "sorted", fill = 0)

  expect_equal(res$ids, c("a", "b", "c", "d"))
  expect_equal(rownames(res$mats$s1), res$ids)
  expect_equal(rownames(res$mats$s2), res$ids)
  expect_equal(unname(res$mats$s2["a", ]), c(0, 0))
  expect_equal(unname(res$mats$s1["d", ]), c(0, 0))

  expect_true(all(res$obs_mask$s1[c("a", "b", "c")]))
  expect_false(res$obs_mask$s1[["d"]])
  expect_true(all(res$obs_mask$s2[c("b", "d")]))
  expect_false(res$obs_mask$s2[["a"]])
})

test_that("harmonize_union_fill harmonizes columns with partial overlap", {
  x1 <- matrix(1:6, 2, 3, dimnames = list(c("r1", "r2"), c("a", "b", "c")))
  x2 <- matrix(1:4, 2, 2, dimnames = list(c("r1", "r2"), c("b", "d")))

  res <- harmonize_union_fill(list(s1 = x1, s2 = x2), axis = "cols", union_order = "sorted", fill = 0)

  expect_equal(res$ids, c("a", "b", "c", "d"))
  expect_equal(colnames(res$mats$s1), res$ids)
  expect_equal(colnames(res$mats$s2), res$ids)
  expect_equal(unname(res$mats$s2[, "a"]), c(0, 0))
  expect_equal(unname(res$mats$s1[, "d"]), c(0, 0))
})

test_that("harmonize_union_fill respects union_order and union_ids ordering", {
  x1 <- matrix(1:4, 2, 2, dimnames = list(c("b", "a"), c("v1", "v2")))
  x2 <- matrix(1:4, 2, 2, dimnames = list(c("c", "b"), c("v1", "v2")))

  res_first <- harmonize_union_fill(list(s1 = x1, s2 = x2), axis = "rows", union_order = "first_seen", fill = 0)
  expect_equal(res_first$ids, c("b", "a", "c"))

  res_sorted <- harmonize_union_fill(list(s1 = x1, s2 = x2), axis = "rows", union_order = "sorted", fill = 0)
  expect_equal(res_sorted$ids, c("a", "b", "c"))

  res_union <- harmonize_union_fill(
    list(s1 = x1, s2 = x2),
    axis = "rows",
    union_ids = c("c", "b", "a"),
    union_order = "sorted",
    fill = 0
  )
  expect_equal(res_union$ids, c("c", "b", "a"))
})

test_that("harmonize_union_fill computes pair_counts correctly", {
  x1 <- matrix(1:4, 2, 2, dimnames = list(c("a", "b"), c("v1", "v2")))
  x2 <- matrix(1:4, 2, 2, dimnames = list(c("b", "c"), c("v1", "v2")))

  res <- harmonize_union_fill(list(s1 = x1, s2 = x2), axis = "rows", union_order = "sorted", fill = 0)
  expect_equal(res$ids, c("a", "b", "c"))

  pc <- res$pair_counts
  expect_equal(unname(diag(pc)), c(1L, 2L, 1L))
  expect_equal(pc["a", "b"], 1L)
  expect_equal(pc["b", "c"], 1L)
  expect_equal(pc["a", "c"], 0L)
})

test_that("harmonize_union_fill filters by min_coverage and warns on sparse subjects", {
  x1 <- matrix(1:4, 2, 2, dimnames = list(c("a", "b"), c("v1", "v2")))
  x2 <- matrix(1:2, 1, 2, dimnames = list(c("b"), c("v1", "v2")))

  expect_warning(
    res_warn <- harmonize_union_fill(
      list(s1 = x1, s2 = x2),
      axis = "rows",
      union_order = "sorted",
      fill = 0,
      min_coverage = 1L,
      warn_sparse_below = 0.9
    ),
    "sparse coverage"
  )

  expect_equal(names(res_warn$mats), c("s1", "s2"))
  expect_equal(res_warn$dropped, character(0))

  res_drop <- harmonize_union_fill(
    list(s1 = x1, s2 = x2),
    axis = "rows",
    union_order = "sorted",
    fill = 0,
    min_coverage = 2L,
    warn_sparse_below = 0.9
  )
  expect_equal(names(res_drop$mats), "s1")
  expect_equal(res_drop$dropped, "s2")
})

test_that("harmonize_union_fill is idempotent on already-harmonized data", {
  x1 <- matrix(1:4, 2, 2, dimnames = list(c("a", "b"), c("v1", "v2")))
  x2 <- matrix(1:4, 2, 2, dimnames = list(c("b", "c"), c("v1", "v2")))

  res1 <- harmonize_union_fill(list(s1 = x1, s2 = x2), axis = "rows", union_order = "sorted", fill = 0)
  res2 <- harmonize_union_fill(res1$mats, axis = "rows", union_ids = res1$ids, fill = 0)

  expect_equal(res2$ids, res1$ids)
  expect_equal(res2$mats, res1$mats)

  # Provenance is recomputed from dimnames and thus assumes full coverage.
  expect_true(all(vapply(res2$obs_mask, all, logical(1))))
  n_subjects <- length(res2$obs_mask)
  pc_expected <- matrix(
    as.integer(n_subjects),
    nrow = length(res2$ids),
    ncol = length(res2$ids),
    dimnames = list(res2$ids, res2$ids)
  )
  expect_equal(res2$pair_counts, pc_expected)
})
