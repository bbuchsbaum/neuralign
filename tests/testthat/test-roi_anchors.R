test_that("roi_anchor_projector builds normalized projector from label vector", {
  roi <- c("V1", "V1", "MT", NA, "MT", "V1")
  P <- roi_anchor_projector(roi, anchors = c("V1", "MT"), normalize = TRUE, sparse = FALSE)

  expect_equal(dim(P), c(2, 6))
  expect_equal(rownames(P), c("V1", "MT"))
  expect_equal(as.numeric(P["V1", c(1, 2, 6)]), rep(1 / 3, 3))
  expect_equal(as.numeric(P["MT", c(3, 5)]), rep(1 / 2, 2))
  expect_equal(sum(P["V1", ]), 1)
  expect_equal(sum(P["MT", ]), 1)
})

test_that("roi_anchor_projector builds projector from list of indices", {
  roi <- list(V1 = c(1, 2), MT = c(3, 5))
  P <- roi_anchor_projector(roi, n_features = 6, normalize = TRUE, sparse = FALSE)

  expect_equal(dim(P), c(2, 6))
  expect_equal(rownames(P), c("V1", "MT"))
  expect_equal(as.numeric(P["V1", c(1, 2)]), rep(1 / 2, 2))
  expect_equal(as.numeric(P["MT", c(3, 5)]), rep(1 / 2, 2))
})

test_that("roi_anchor_projectors harmonizes anchors by union", {
  roi_by_subject <- list(
    s1 = c("V1", "MT", "MT"),
    s2 = c("V1", "V1", "V1")
  )
  proj <- roi_anchor_projectors(roi_by_subject, anchor_policy = "union", sparse = FALSE)

  expect_equal(names(proj), c("s1", "s2"))
  expect_equal(rownames(proj$s1), c("MT", "V1"))
  expect_equal(rownames(proj$s2), c("MT", "V1"))

  # s2 has no MT voxels, so MT row should be all zeros
  expect_equal(sum(proj$s2["MT", ]), 0)
  expect_equal(sum(proj$s2["V1", ]), 1)
})

test_that("set_roi_guidance attaches projector channels to AlignmentData", {
  set.seed(1)
  data <- AlignmentData(list(
    s1 = matrix(rnorm(12), 3, 4),
    s2 = matrix(rnorm(15), 3, 5)
  ))

  roi <- list(
    s1 = c("V1", "V1", "MT"),
    s2 = c("V1", "V1", "V1")
  )

  data2 <- set_roi_guidance(
    data,
    roi_by_subject = roi,
    channel_name = "roi",
    anchor_policy = "union",
    sparse = TRUE,
    validate = TRUE
  )

  g <- get_guidance(data2, type = "projector")
  expect_true("roi" %in% names(g$s1))
  expect_true("roi" %in% names(g$s2))
  expect_true(inherits(g$s1$roi$value, "Matrix"))
  expect_equal(ncol(g$s1$roi$value), 3)
  expect_equal(ncol(g$s2$roi$value), 3)
})


# ---------- Additional ROI anchor tests ----------

test_that("roi_anchor_projector from matrix input with anchors subset", {
  P_full <- matrix(c(1, 0, 0, 0, 1, 0), 2, 3, dimnames = list(c("V1", "MT"), NULL))
  P_sub <- roi_anchor_projector(P_full, anchors = "V1", normalize = TRUE, sparse = FALSE)
  expect_equal(dim(P_sub), c(1, 3))
  expect_equal(rownames(P_sub), "V1")
})

test_that("roi_anchor_projector from matrix errors on missing anchors", {
  P <- matrix(1, 2, 3, dimnames = list(c("V1", "MT"), NULL))
  expect_error(roi_anchor_projector(P, anchors = "FFA"), "missing anchors")
})

test_that("roi_anchor_projector from matrix without rownames gets default names", {
  P <- matrix(c(1, 0, 0, 1), 2, 2)
  res <- roi_anchor_projector(P, normalize = FALSE, sparse = FALSE)
  expect_equal(rownames(res), c("anchor-1", "anchor-2"))
})

test_that("roi_anchor_projector from list with unnormalized weights", {
  roi <- list(V1 = c(1, 2), MT = c(3))
  P <- roi_anchor_projector(roi, n_features = 4, normalize = FALSE, sparse = FALSE)
  expect_equal(dim(P), c(2, 4))
  expect_equal(sum(P["V1", ]), 2)  # Two 1s, not normalized
  expect_equal(sum(P["MT", ]), 1)
})

test_that("roi_anchor_projector from list infers n_features from max index", {
  roi <- list(a = c(1, 5), b = c(3))
  P <- roi_anchor_projector(roi, normalize = TRUE, sparse = FALSE)
  expect_equal(ncol(P), 5)
})

test_that("roi_anchor_projector from list errors on non-positive n_features", {
  roi <- list(V1 = c(1, 2))
  expect_error(
    roi_anchor_projector(roi, n_features = 0, normalize = TRUE, sparse = FALSE),
    "positive integer"
  )
  expect_error(
    roi_anchor_projector(roi, n_features = -1, normalize = TRUE, sparse = FALSE),
    "positive integer"
  )
})

test_that("roi_anchor_projector from list errors on empty with no n_features", {
  roi <- list(a = integer(0))
  expect_error(roi_anchor_projector(roi), "infer n_features")
})

test_that("roi_anchor_projector from list with sparse=TRUE returns sparse Matrix", {
  skip_if_not_installed("Matrix")
  roi <- list(V1 = c(1, 2), MT = c(3, 4))
  P <- roi_anchor_projector(roi, n_features = 5, normalize = TRUE, sparse = TRUE)
  expect_true(inherits(P, "Matrix"))
  expect_equal(nrow(P), 2)
  expect_equal(ncol(P), 5)
})

test_that("roi_anchor_projector from vector with NAs and empty strings", {
  roi <- c("V1", NA, "", "MT", "V1")
  P <- roi_anchor_projector(roi, normalize = TRUE, sparse = FALSE)
  expect_equal(sort(rownames(P)), c("MT", "V1"))
  # V1 at indices 1 and 5
  expect_equal(sum(P["V1", ] > 0), 2)
  # MT at index 4
  expect_equal(sum(P["MT", ] > 0), 1)
})

test_that("roi_anchor_projector errors on invalid input", {
  expect_error(roi_anchor_projector(NULL), "vector, list, or matrix")
})

test_that("roi_anchor_projectors harmonizes by intersection", {
  roi_by_subject <- list(
    s1 = c("V1", "MT", "MT"),
    s2 = c("V1", "V1", "V1")
  )
  proj <- roi_anchor_projectors(roi_by_subject, anchor_policy = "intersection", sparse = FALSE)
  # Intersection is just "V1" (MT not in s2)
  expect_equal(rownames(proj$s1), "V1")
  expect_equal(rownames(proj$s2), "V1")
})

test_that("roi_anchor_projectors errors on unnamed list", {
  expect_error(
    roi_anchor_projectors(list(c("V1", "MT"), c("V1"))),
    "named list"
  )
})

test_that("set_roi_guidance errors on non-AlignmentData", {
  expect_error(set_roi_guidance("not_adat", list()), "AlignmentData")
})

test_that("set_roi_guidance errors on bad channel_name", {
  data <- AlignmentData(list(s1 = matrix(1, 3, 2)))
  expect_error(
    set_roi_guidance(data, roi_by_subject = list(s1 = c("a", "b", "c")), channel_name = ""),
    "non-empty"
  )
})

test_that("set_roi_guidance errors on missing subjects in roi_by_subject", {
  data <- AlignmentData(list(
    s1 = matrix(1, 3, 2),
    s2 = matrix(1, 4, 2)
  ))
  roi <- list(s1 = c("a", "b", "c"))  # missing s2
  expect_error(set_roi_guidance(data, roi_by_subject = roi), "missing subjects")
})

test_that("roi_anchor_projectors intersection errors when no anchors overlap", {
  roi_by_subject <- list(
    s1 = c("A", "A"),
    s2 = c("B", "B")
  )
  expect_error(
    roi_anchor_projectors(roi_by_subject, anchor_policy = "intersection", sparse = FALSE),
    "No anchors remain after harmonization"
  )
})

test_that("roi_anchor_projector errors on empty roi vector when anchors are not supplied", {
  roi <- c(NA, "", NA)
  expect_error(
    roi_anchor_projector(roi, sparse = FALSE),
    "No anchors available"
  )
})

test_that("roi_anchor_projector errors on non-positive n_features", {
  roi <- list(V1 = c(1, 2))
  expect_error(
    roi_anchor_projector(roi, n_features = 0, sparse = FALSE),
    "positive integer"
  )
})

test_that("set_roi_guidance validates projector dimension mismatch", {
  data <- AlignmentData(list(s1 = matrix(1, 3, 2)))
  roi <- list(s1 = c("a", "b", "c", "d")) # length 4, but subject has 3 features
  expect_error(
    set_roi_guidance(data, roi_by_subject = roi, sparse = FALSE, validate = TRUE),
    "dimension mismatch"
  )
})
