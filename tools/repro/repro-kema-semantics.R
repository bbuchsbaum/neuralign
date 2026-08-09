#!/usr/bin/env Rscript

# Reproduce neuralign-7nu.4.5 from a repository checkout.
suppressPackageStartupMessages(devtools::load_all(".", quiet = TRUE))

if (!requireNamespace("manifoldalign", quietly = TRUE)) {
  stop("Install manifoldalign to run this reproduction.", call. = FALSE)
}

set.seed(1)
p <- 40L
n <- 40L
k <- 5L
latent <- matrix(stats::rnorm(k * n), k, n)
subjects <- lapply(seq_len(2L), function(i) {
  basis <- qr.Q(qr(matrix(stats::rnorm(p * k), p, k)))
  basis %*% latent + matrix(stats::rnorm(p * n, sd = 1e-2), p, n)
})
names(subjects) <- c("s1", "s2")
labels <- sprintf("obs-%02d", seq_len(n))
data <- AlignmentData(subjects, obs_labels = labels)

hd <- neuralign:::.ma_build_hyperdesign_obs(
  data,
  labels = labels,
  label_name = "label"
)
upstream <- manifoldalign::kema(
  hd,
  y = label,
  preproc = multivarious::pass(),
  ncomp = k,
  solver = "exact",
  lambda = 1e-2
)

neuralign:::.clear_registry()
neuralign:::.register_kema()
result <- fit_alignment(
  data,
  method = "kema",
  reference = "s1",
  ncomp = k,
  compute_quality = FALSE
)
aligned <- get_aligned(result)

cat("upstream scores:", paste(dim(upstream$s), collapse = " x "), "\n")
cat("upstream primal vectors:", paste(dim(upstream$v), collapse = " x "), "\n")
cat("neuralign s1 embedding:", paste(dim(aligned$s1), collapse = " x "), "\n")
cat(
  "exact score-block match:",
  isTRUE(all.equal(t(aligned$s1), as.matrix(upstream$s)[seq_len(n), , drop = FALSE])),
  "\n"
)

sessionInfo()
