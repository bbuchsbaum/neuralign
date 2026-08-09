#!/usr/bin/env Rscript

# Reproduce neuralign-tzc against the deterministic graph oracle.
suppressPackageStartupMessages(devtools::load_all(".", quiet = TRUE))

if (!requireNamespace("manifoldalign", quietly = TRUE)) {
  stop("Install manifoldalign to run this reproduction.", call. = FALSE)
}

set.seed(1)
p <- 60L
n_obs <- 40L
target <- matrix(0, p, n_obs)
for (i in seq_len(p)) {
  target[i, ] <- sin(seq(0, 2 * pi, length.out = n_obs) * (i / 7)) + i / p
}
permutation <- sample.int(p)
source <- target[permutation, , drop = FALSE]
truth <- match(seq_len(p), permutation)

hd <- neuralign:::.ma_build_pair_hyperdesign_features(
  target,
  source,
  target_name = "target",
  source_name = "source"
)
fit <- manifoldalign::cone_align(
  hd,
  ncomp = 10L,
  sigma = 0.73,
  lambda = 0.1,
  max_iter = 30L,
  tol = 0.01
)
assignment <- as.integer(fit$assignment)

cat("target-to-source assignment accuracy:", mean(assignment == truth), "\n")
cat(
  "assignment reconstruction RMSE:",
  sqrt(mean((target - source[assignment, , drop = FALSE])^2)),
  "\n"
)
cat(
  "oracle reconstruction RMSE:",
  sqrt(mean((target - source[truth, , drop = FALSE])^2)),
  "\n"
)
cat("adapter orientation:", "target rows -> source rows", "\n")
cat("neuralign registration:")
tryCatch(
  neuralign:::.register_cone(),
  error = function(e) cat(" disabled -", conditionMessage(e), "\n")
)

sessionInfo()
