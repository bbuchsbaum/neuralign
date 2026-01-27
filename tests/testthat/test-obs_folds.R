test_that("create_obs_folds(method='run') creates leave-one-run-out folds", {
  runs <- c("r1", "r1", "r1", "r2", "r2", "r3", "r3")
  spec <- create_obs_folds(runs, method = "run")
  expect_equal(spec$axis, "observation")
  expect_equal(spec$method, "run")
  expect_equal(spec$n_folds, 3)
  expect_setequal(spec$fold_ids, c("r1", "r2", "r3"))

  f1 <- spec$folds[["r1"]]
  expect_equal(f1$test_idx, 1:3)
  expect_equal(sort(c(f1$train_idx, f1$test_idx)), seq_along(runs))
  expect_equal(intersect(f1$train_idx, f1$test_idx), integer(0))
})

test_that("create_obs_folds supports per-subject run ids with union policy", {
  obs_ids <- list(
    s1 = c("A", "A", "B", "B"),
    s2 = c("C", "D", "D", "D")
  )

  spec <- NULL
  expect_warning(spec <- create_obs_folds(obs_ids, method = "run", id_policy = "union"))
  expect_equal(spec$axis, "observation")
  expect_equal(spec$method, "run")
  expect_setequal(spec$fold_ids, c("A", "B", "C", "D"))

  fold_A <- spec$folds[["A"]]
  expect_equal(fold_A$s1$test_idx, c(1, 2))
  expect_equal(fold_A$s2$test_idx, integer(0))
})

test_that("create_obs_folds(method='blocked_time') creates contiguous blocks with guard", {
  obs <- seq_len(10)
  spec <- create_obs_folds(obs, method = "blocked_time", k = 2, guard_tr = 1)
  expect_equal(spec$n_folds, 2)
  b1 <- spec$folds[["block-1"]]
  b2 <- spec$folds[["block-2"]]
  expect_equal(b1$test_idx, 1:5)
  expect_equal(b2$test_idx, 6:10)
  # Guard around block-1 removes indices 1:6 from training
  expect_equal(b1$train_idx, 7:10)
})

test_that("create_obs_folds errors when guard_tr makes training empty", {
  expect_error(create_obs_folds(seq_len(6), method = "blocked_time", k = 2, guard_tr = 3))
})


# ---------- Additional obs fold tests ----------

test_that("create_obs_folds(method='run') with guard_tr removes nearby timepoints", {
  runs <- c("r1", "r1", "r1", "r2", "r2", "r2")
  spec <- create_obs_folds(runs, method = "run", guard_tr = 1)

  # For fold "r1" (test = 1:3), guard should remove index 4 from train
  f1 <- spec$folds[["r1"]]
  expect_false(4 %in% f1$train_idx)
  expect_equal(sort(f1$test_idx), 1:3)

  # For fold "r2" (test = 4:6), guard should remove index 3 from train
  f2 <- spec$folds[["r2"]]
  expect_false(3 %in% f2$train_idx)
  expect_equal(sort(f2$test_idx), 4:6)
})

test_that("create_obs_folds(method='run') errors with < 2 unique runs", {
  expect_error(
    create_obs_folds(c("r1", "r1", "r1"), method = "run"),
    ">= 2 unique run ids"
  )
})

test_that("create_obs_folds(method='run') errors on NA ids", {
  expect_error(
    create_obs_folds(c("r1", NA, "r2"), method = "run"),
    "NA"
  )
})

test_that("create_obs_folds errors on negative guard_tr", {
  expect_error(
    create_obs_folds(c("r1", "r2"), method = "run", guard_tr = -1),
    "non-negative"
  )
})

test_that("create_obs_folds(method='blocked_time') with k=3 creates 3 blocks", {
  obs <- seq_len(12)
  spec <- create_obs_folds(obs, method = "blocked_time", k = 3, guard_tr = 0)
  expect_equal(spec$n_folds, 3)
  expect_equal(spec$axis, "observation")
  expect_equal(spec$method, "blocked_time")

  # Each block should cover 4 observations
  expect_equal(length(spec$folds[["block-1"]]$test_idx), 4)
  expect_equal(length(spec$folds[["block-2"]]$test_idx), 4)
  expect_equal(length(spec$folds[["block-3"]]$test_idx), 4)
})

test_that("create_obs_folds per-subject with intersection policy", {
  obs_ids <- list(
    s1 = c("A", "A", "B", "B"),
    s2 = c("A", "B", "B", "B")
  )

  spec <- create_obs_folds(obs_ids, method = "run", id_policy = "intersection")
  expect_equal(spec$axis, "observation")
  expect_setequal(spec$fold_ids, c("A", "B"))
  expect_equal(spec$id_policy, "intersection")
})

test_that("create_obs_folds per-subject errors on NA in per-subject ids", {
  obs_ids <- list(
    s1 = c("A", NA, "B"),
    s2 = c("A", "B", "B")
  )
  expect_error(
    create_obs_folds(obs_ids, method = "run"),
    "NA"
  )
})


# ---------- More obs fold tests for coverage ----------

test_that("create_obs_folds blocked_time per-subject creates nested folds", {
  obs_ids <- list(
    s1 = paste0("t", 1:8),
    s2 = paste0("t", 1:8)
  )
  spec <- create_obs_folds(obs_ids, method = "blocked_time", k = 2, guard_tr = 0)

  expect_equal(spec$method, "blocked_time")
  expect_equal(spec$n_folds, 2)
  # Each fold should contain per-subject entries
  for (fid in spec$fold_ids) {
    fold <- spec$folds[[fid]]
    expect_true("s1" %in% names(fold))
    expect_true("s2" %in% names(fold))
    expect_true(length(fold$s1$test_idx) > 0)
    expect_true(length(fold$s1$train_idx) > 0)
  }
})

test_that("create_obs_folds per-subject run with guard_tr removes neighbors", {
  obs_ids <- list(
    s1 = c("A", "A", "B", "B", "C", "C"),
    s2 = c("A", "A", "B", "B", "C", "C")
  )

  spec <- create_obs_folds(obs_ids, method = "run", id_policy = "intersection", guard_tr = 1)

  fold_B <- spec$folds[["B"]]
  # s1 test for B = indices 3,4; guard should remove indices 2 and 5 from train
  expect_false(2 %in% fold_B$s1$train_idx)
  expect_false(5 %in% fold_B$s1$train_idx)
})

test_that("create_obs_folds per-subject run errors with < 2 fold ids", {
  obs_ids <- list(
    s1 = c("A", "A", "A"),
    s2 = c("A", "A", "A")
  )

  expect_error(
    create_obs_folds(obs_ids, method = "run"),
    ">= 2 fold ids"
  )
})

test_that(".as_obs_ids_list returns NULL for atomic input", {
  result <- neuralign:::.as_obs_ids_list(c("a", "b", "c"))
  expect_null(result)
})

test_that(".as_obs_ids_list errors on non-list non-atomic input", {
  expect_error(
    neuralign:::.as_obs_ids_list(environment()),
    "atomic vector or a list"
  )
})

test_that(".as_obs_ids_list errors on unnamed list", {
  expect_error(
    neuralign:::.as_obs_ids_list(list(c("a", "b"))),
    "named list"
  )
})

test_that(".blocked_time_folds errors on n_obs < 2", {
  expect_error(
    neuralign:::.blocked_time_folds(1, k = 2, guard_tr = 0),
    "n_obs >= 2"
  )
})

test_that(".blocked_time_folds errors on k < 2", {
  expect_error(
    neuralign:::.blocked_time_folds(10, k = 1, guard_tr = 0),
    "k.*>= 2"
  )
})

test_that(".blocked_time_folds errors on negative guard_tr", {
  expect_error(
    neuralign:::.blocked_time_folds(10, k = 2, guard_tr = -1),
    "guard_tr.*>= 0"
  )
})

test_that(".blocked_time_folds errors on k > n_obs", {
  expect_error(
    neuralign:::.blocked_time_folds(3, k = 5, guard_tr = 0),
    "cannot exceed"
  )
})

test_that("create_obs_folds with seed parameter works", {
  runs <- c("r1", "r1", "r2", "r2")
  spec <- create_obs_folds(runs, method = "run", seed = 42)
  expect_equal(spec$n_folds, 2)
})

test_that("create_obs_folds per-subject guard_tr too large errors", {
  obs_ids <- list(
    s1 = c("A", "B"),
    s2 = c("A", "B")
  )

  expect_error(
    create_obs_folds(obs_ids, method = "run", guard_tr = 5),
    "guard_tr too large"
  )
})
