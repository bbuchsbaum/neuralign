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

