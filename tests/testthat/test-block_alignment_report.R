test_that("block_alignment_report combines feature building, diagnostics, and quality", {
  blocks_by_subject <- list(
    sub1 = list(
      a = alignment_feature_block(matrix(1:6, 3, 2), name = "a", feature_names = c("f1", "f2", "f3")),
      b = alignment_feature_block(matrix(c(10, 20, 30, 40), 2, 2), name = "b", feature_names = c("g1", "g2"))
    ),
    sub2 = list(
      a = alignment_feature_block(matrix(101:106, 3, 2), name = "a", feature_names = c("f1", "f2", "f3")),
      b = alignment_feature_block(matrix(c(50, 60, 70, 80), 2, 2), name = "b", feature_names = c("g1", "g2"))
    )
  )

  features <- build_alignment_features(blocks_by_subject, harmonize = "intersection", min_features = 2)
  aligned <- features$matrices
  reference <- aligned$sub1

  report <- block_alignment_report(
    blocks_by_subject,
    aligned = aligned,
    reference = reference,
    quality_metrics = c("correlation", "reconstruction"),
    per_block_quality = TRUE,
    harmonize = "intersection",
    min_features = 2
  )

  expect_s3_class(report, "block_alignment_report")
  expect_true(is.list(report$features))
  expect_true(is.data.frame(report$features$coverage))
  expect_true(is.data.frame(report$features$block_row_ranges))
  expect_true(is.list(report$diagnostics))
  expect_true(is.list(report$quality))
  expect_true(!is.null(report$quality$reconstruction_errors))
  expect_true(is.list(report$per_block_quality))
  expect_true(is.data.frame(report$per_block_quality$by_subject))
  expect_true(is.data.frame(report$per_block_quality$summary))

  # Row ranges reflect stacking order a (3 rows) then b (2 rows).
  expect_equal(report$features$block_row_ranges$row_start, c(1L, 4L))
  expect_equal(report$features$block_row_ranges$row_end, c(3L, 5L))

  # For sub2 vs sub1, block-wise differences are constant (100 for a, 40 for b).
  by_subj <- report$per_block_quality$by_subject
  sub2_a <- by_subj[by_subj$subject == "sub2" & by_subj$block == "a", ]
  sub2_b <- by_subj[by_subj$subject == "sub2" & by_subj$block == "b", ]
  expect_equal(unname(sub2_a$rmse), 100, tolerance = 1e-12)
  expect_equal(unname(sub2_a$frobenius), sqrt(60000), tolerance = 1e-12)
  expect_equal(unname(sub2_b$rmse), 40, tolerance = 1e-12)
  expect_equal(unname(sub2_b$frobenius), 80, tolerance = 1e-12)

  expect_output(print(report), "block_alignment_report")
})
