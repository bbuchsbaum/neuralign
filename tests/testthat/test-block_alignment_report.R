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
    harmonize = "intersection",
    min_features = 2
  )

  expect_s3_class(report, "block_alignment_report")
  expect_true(is.list(report$features))
  expect_true(is.data.frame(report$features$coverage))
  expect_true(is.list(report$diagnostics))
  expect_true(is.list(report$quality))
  expect_true(!is.null(report$quality$reconstruction_errors))

  expect_output(print(report), "block_alignment_report")
})

