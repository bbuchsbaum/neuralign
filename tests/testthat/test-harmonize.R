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


# ---------- Edge cases ----------

test_that("harmonize_union_fill works with identical IDs (no-op)", {
  x1 <- matrix(1:6, 3, 2, dimnames = list(c("a", "b", "c"), c("v1", "v2")))
  x2 <- matrix(7:12, 3, 2, dimnames = list(c("a", "b", "c"), c("v1", "v2")))

  res <- harmonize_union_fill(list(s1 = x1, s2 = x2), axis = "rows", fill = 0)

  expect_equal(res$ids, c("a", "b", "c"))
  expect_equal(unname(res$mats$s1), unname(x1))
  expect_equal(unname(res$mats$s2), unname(x2))
  # All observed
  expect_true(all(res$obs_mask$s1))
  expect_true(all(res$obs_mask$s2))
  expect_equal(length(res$dropped), 0)
})

test_that("harmonize_union_fill with completely disjoint IDs fills everything", {
  x1 <- matrix(1:4, 2, 2, dimnames = list(c("a", "b"), c("v1", "v2")))
  x2 <- matrix(5:8, 2, 2, dimnames = list(c("c", "d"), c("v1", "v2")))

  res <- harmonize_union_fill(list(s1 = x1, s2 = x2), axis = "rows", union_order = "sorted", fill = -1)

  expect_equal(res$ids, c("a", "b", "c", "d"))
  # s1 has "c","d" filled with -1
  expect_equal(unname(res$mats$s1["c", ]), c(-1, -1))
  expect_equal(unname(res$mats$s1["d", ]), c(-1, -1))
  # s2 has "a","b" filled with -1
  expect_equal(unname(res$mats$s2["a", ]), c(-1, -1))
  expect_equal(unname(res$mats$s2["b", ]), c(-1, -1))

  # Pair counts: "a"-"b" co-occur in s1 only, "c"-"d" in s2 only
  expect_equal(res$pair_counts["a", "c"], 0L)
  expect_equal(res$pair_counts["a", "b"], 1L)
  expect_equal(res$pair_counts["c", "d"], 1L)
})

test_that("harmonize_union_fill with non-zero fill value", {
  x1 <- matrix(1:4, 2, 2, dimnames = list(c("a", "b"), c("v1", "v2")))
  x2 <- matrix(5:6, 1, 2, dimnames = list("c", c("v1", "v2")))

  res <- harmonize_union_fill(
    list(s1 = x1, s2 = x2),
    axis = "rows",
    fill = 99,
    warn_sparse_below = 0
  )

  expect_equal(unname(res$mats$s2["a", ]), c(99, 99))
  expect_equal(unname(res$mats$s2["b", ]), c(99, 99))
  expect_equal(unname(res$mats$s1["c", ]), c(99, 99))
})

test_that("harmonize_union_fill with sparse Matrix inputs preserves sparsity", {
  skip_if_not_installed("Matrix")

  x1 <- Matrix::sparseMatrix(
    i = c(1, 2), j = c(1, 2), x = c(10, 20),
    dims = c(3, 2),
    dimnames = list(c("a", "b", "c"), c("v1", "v2"))
  )
  x2 <- Matrix::sparseMatrix(
    i = c(1), j = c(1), x = 30,
    dims = c(2, 2),
    dimnames = list(c("b", "d"), c("v1", "v2"))
  )

  res <- harmonize_union_fill(list(s1 = x1, s2 = x2), axis = "rows", fill = 0)

  # With fill=0 and sparse inputs, output should be sparse
  expect_true(inherits(res$mats$s1, "Matrix"))
  expect_true(inherits(res$mats$s2, "Matrix"))
  expect_equal(nrow(res$mats$s1), 4)  # union of a,b,c,d
  expect_equal(as.matrix(res$mats$s2)[1, ], c(v1 = 0, v2 = 0))  # "a" not in s2
})

test_that("harmonize_union_fill errors on unnamed list", {
  x1 <- matrix(1:4, 2, 2, dimnames = list(c("a", "b"), c("v1", "v2")))
  expect_error(
    harmonize_union_fill(list(x1), axis = "rows"),
    "non-empty named list"
  )
})

test_that("harmonize_union_fill errors on duplicate IDs within a subject", {
  x1 <- matrix(1:4, 2, 2, dimnames = list(c("a", "a"), c("v1", "v2")))
  expect_error(
    harmonize_union_fill(list(s1 = x1), axis = "rows"),
    "Duplicate ids"
  )
})

test_that("harmonize_union_fill errors on invalid min_coverage", {
  x1 <- matrix(1:4, 2, 2, dimnames = list(c("a", "b"), c("v1", "v2")))
  expect_error(
    harmonize_union_fill(list(s1 = x1), axis = "rows", min_coverage = -1),
    "non-negative"
  )
})

test_that("harmonize_union_fill errors on invalid warn_sparse_below", {
  x1 <- matrix(1:4, 2, 2, dimnames = list(c("a", "b"), c("v1", "v2")))
  expect_error(
    harmonize_union_fill(list(s1 = x1), axis = "rows", warn_sparse_below = 2),
    "\\[0, 1\\]"
  )
})

test_that("harmonize_union_fill drops all subjects errors", {
  x1 <- matrix(1:2, 1, 2, dimnames = list("a", c("v1", "v2")))
  x2 <- matrix(3:4, 1, 2, dimnames = list("b", c("v1", "v2")))

  # Each subject has 1 ID. min_coverage = 2 drops both.
  expect_error(
    harmonize_union_fill(list(s1 = x1, s2 = x2), axis = "rows", min_coverage = 2),
    "All subjects were dropped"
  )
})

test_that("harmonize_union_fill coverage data frame is correct", {
  x1 <- matrix(1:6, 3, 2, dimnames = list(c("a", "b", "c"), c("v1", "v2")))
  x2 <- matrix(1:4, 2, 2, dimnames = list(c("b", "d"), c("v1", "v2")))
  x3 <- matrix(1:4, 2, 2, dimnames = list(c("a", "b"), c("v1", "v2")))

  res <- harmonize_union_fill(list(s1 = x1, s2 = x2, s3 = x3), axis = "rows", union_order = "sorted")

  cov <- res$coverage
  expect_equal(cov$id, c("a", "b", "c", "d"))
  # "a" in s1, s3 = 2; "b" in all = 3; "c" in s1 = 1; "d" in s2 = 1
  expect_equal(cov$n_subjects, c(2L, 3L, 1L, 1L))
})

test_that("harmonize_union_fill errors on union_ids with duplicates", {
  x1 <- matrix(1:4, 2, 2, dimnames = list(c("a", "b"), c("v1", "v2")))
  expect_error(
    harmonize_union_fill(list(s1 = x1), axis = "rows", union_ids = c("a", "a", "b")),
    "duplicates"
  )
})

test_that("harmonize_union_fill errors when union_ids misses observed IDs", {
  x1 <- matrix(1:4, 2, 2, dimnames = list(c("a", "b"), c("v1", "v2")))
  expect_error(
    harmonize_union_fill(list(s1 = x1), axis = "rows", union_ids = c("a")),
    "not present in .union_ids"
  )
})


# ---------- More harmonize coverage tests ----------

test_that("harmonize_union_fill errors on ID length mismatch", {
  x1 <- matrix(1:6, 3, 2, dimnames = list(c("a", "b", "c"), c("v1", "v2")))

  # Provide ids that don't match the matrix row count
  expect_error(
    harmonize_union_fill(
      list(s1 = x1),
      ids = list(s1 = c("a", "b")),  # 2 ids but 3 rows
      axis = "rows"
    ),
    "ID length mismatch"
  )
})

test_that("harmonize_union_fill errors on ids missing subjects", {
  x1 <- matrix(1:4, 2, 2, dimnames = list(c("a", "b"), c("v1", "v2")))
  x2 <- matrix(5:8, 2, 2, dimnames = list(c("c", "d"), c("v1", "v2")))

  # Named ids list missing s2
  expect_error(
    harmonize_union_fill(
      list(s1 = x1, s2 = x2),
      ids = list(s1 = c("a", "b")),
      axis = "rows"
    ),
    "missing subjects"
  )
})

test_that("harmonize_union_fill with named ids list reorders to match mats", {
  x1 <- matrix(1:4, 2, 2, dimnames = list(c("a", "b"), c("v1", "v2")))
  x2 <- matrix(5:8, 2, 2, dimnames = list(c("c", "d"), c("v1", "v2")))

  res <- harmonize_union_fill(
    list(s1 = x1, s2 = x2),
    ids = list(s2 = c("c", "d"), s1 = c("a", "b")),  # reversed order
    axis = "rows",
    union_order = "sorted"
  )

  expect_equal(res$ids, c("a", "b", "c", "d"))
  expect_equal(rownames(res$mats$s1), res$ids)
})

test_that("harmonize_union_fill unnamed ids list gets names from mats", {
  x1 <- matrix(1:4, 2, 2, dimnames = list(c("a", "b"), c("v1", "v2")))
  x2 <- matrix(5:8, 2, 2, dimnames = list(c("b", "c"), c("v1", "v2")))

  res <- harmonize_union_fill(
    list(s1 = x1, s2 = x2),
    ids = list(c("a", "b"), c("b", "c")),  # unnamed
    axis = "rows",
    union_order = "sorted"
  )

  expect_equal(res$ids, c("a", "b", "c"))
})

test_that("harmonize_union_fill cols axis with sparse Matrix preserves sparsity", {
  skip_if_not_installed("Matrix")

  x1 <- Matrix::sparseMatrix(
    i = c(1, 2), j = c(1, 3), x = c(10, 20),
    dims = c(2, 3),
    dimnames = list(c("r1", "r2"), c("a", "b", "c"))
  )
  x2 <- Matrix::sparseMatrix(
    i = c(1), j = c(2), x = 30,
    dims = c(2, 2),
    dimnames = list(c("r1", "r2"), c("b", "d"))
  )

  res <- harmonize_union_fill(list(s1 = x1, s2 = x2), axis = "cols", fill = 0)

  expect_true(inherits(res$mats$s1, "Matrix"))
  expect_true(inherits(res$mats$s2, "Matrix"))
  expect_equal(ncol(res$mats$s1), 4)  # union of a,b,c,d
  expect_equal(ncol(res$mats$s2), 4)
})

test_that("harmonize_union_fill cols axis with dense matrix works", {
  x1 <- matrix(1:6, 2, 3, dimnames = list(c("r1", "r2"), c("a", "b", "c")))
  x2 <- matrix(7:10, 2, 2, dimnames = list(c("r1", "r2"), c("b", "d")))

  res <- harmonize_union_fill(list(s1 = x1, s2 = x2), axis = "cols", fill = -1)

  expect_equal(colnames(res$mats$s2), colnames(res$mats$s1))
  expect_equal(unname(res$mats$s2[, "a"]), c(-1, -1))
  expect_equal(unname(res$mats$s2[, "c"]), c(-1, -1))
})

test_that("harmonize_union_fill errors on non-matrix entries", {
  expect_error(
    harmonize_union_fill(list(s1 = "not_a_matrix"), axis = "rows"),
    "must be matrices"
  )
})

test_that("harmonize_union_fill errors on empty list", {
  expect_error(
    harmonize_union_fill(list(), axis = "rows"),
    "non-empty named list"
  )
})

test_that("harmonize_union_fill with NULL ids infers from dimnames", {
  x1 <- matrix(1:4, 2, 2, dimnames = list(c("a", "b"), NULL))
  res <- harmonize_union_fill(list(s1 = x1), axis = "rows")
  expect_equal(res$ids, c("a", "b"))
})

test_that("harmonize_union_fill with warn_sparse_below=0 suppresses warnings", {
  x1 <- matrix(1:4, 2, 2, dimnames = list(c("a", "b"), c("v1", "v2")))
  x2 <- matrix(1:2, 1, 2, dimnames = list("c", c("v1", "v2")))

  # warn_sparse_below=0 should not warn
  expect_no_warning(
    harmonize_union_fill(
      list(s1 = x1, s2 = x2),
      axis = "rows",
      warn_sparse_below = 0
    )
  )
})
