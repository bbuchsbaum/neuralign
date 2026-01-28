test_that(".ma_get_single_subject_obs_labels handles atomic and per-subject lists", {
  X <- matrix(1, 3, 2)

  adat_atomic <- AlignmentData(list(s1 = X), obs_labels = c("a", "b"))
  expect_equal(neuralign:::.ma_get_single_subject_obs_labels(adat_atomic), c("a", "b"))

  adat_named <- AlignmentData(list(s1 = X), obs_labels = list(s1 = c("a", "b")))
  expect_equal(neuralign:::.ma_get_single_subject_obs_labels(adat_named), c("a", "b"))

  adat_bad <- AlignmentData(list(s1 = X), obs_labels = list(s1 = c("a", "b")))
  adat_bad@obs_labels <- list(s2 = c("a", "b"), s3 = c("a", "b"))
  expect_error(neuralign:::.ma_get_single_subject_obs_labels(adat_bad), "does not contain subject")
})
