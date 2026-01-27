test_that("fit_alignment works with kprocrustes end-to-end", {
  neuralign:::.register_kprocrustes()

  set.seed(1)
  q <- 6
  r <- 3
  effects <- paste0("e", seq_len(q))
  K <- diag(q)
  dimnames(K) <- list(effects, effects)

  Uref <- k_orthonormalize(matrix(rnorm(q * r), q, r), K)
  Xref <- t(Uref) # r x q
  colnames(Xref) <- effects

  R2 <- qr.Q(qr(matrix(rnorm(r * r), r, r)))
  if (det(R2) < 0) R2[, r] <- -R2[, r]
  X2 <- t(R2) %*% Xref
  colnames(X2) <- effects

  adat <- AlignmentData(
    list(s1 = Xref, s2 = X2),
    design = list(K = K, effects = effects)
  )

  res <- fit_alignment(adat, method = "kprocrustes", reference = "s1")
  expect_s4_class(res, "AlignmentResult")

  model <- get_model(res)
  expect_equal(model@method, "kprocrustes")
  expect_equal(dim(model@transforms$s1), c(r, r))
  expect_equal(dim(model@transforms$s2), c(r, r))
  expect_equal(model@transforms$s1, diag(r), tolerance = 1e-10)
  expect_equal(model@transforms$s2 %*% X2, model@reference_data, tolerance = 1e-6)
})

test_that("kprocrustes enforces needs_design", {
  neuralign:::.register_kprocrustes()

  data_list <- list(
    s1 = matrix(rnorm(12), 3, 4),
    s2 = matrix(rnorm(12), 3, 4)
  )
  adat <- AlignmentData(data_list)

  expect_error(
    fit_alignment(adat, method = "kprocrustes", reference = "s1"),
    "requires design"
  )
})

test_that("kprocrustes validates design K dimensions", {
  neuralign:::.register_kprocrustes()

  set.seed(2)
  q <- 5
  r <- 3
  effects <- paste0("e", seq_len(q + 1L))
  K <- diag(q) # wrong q relative to effects
  dimnames(K) <- list(paste0("k", seq_len(q)), paste0("k", seq_len(q)))

  X <- matrix(rnorm(r * (q + 1L)), r, q + 1L)
  colnames(X) <- effects
  adat <- AlignmentData(list(s1 = X, s2 = X), design = list(K = K, effects = effects))

  expect_error(
    fit_alignment(adat, method = "kprocrustes", reference = "s1"),
    "design\\$effects must have length q"
  )
})

test_that("kprocrustes supports subject-axis CV", {
  neuralign:::.register_kprocrustes()

  set.seed(3)
  q <- 6
  r <- 3
  effects <- paste0("e", seq_len(q))
  K <- diag(q)
  dimnames(K) <- list(effects, effects)

  Uref <- k_orthonormalize(matrix(rnorm(q * r), q, r), K)
  Xref <- t(Uref)
  colnames(Xref) <- effects

  rot <- function() {
    R <- qr.Q(qr(matrix(rnorm(r * r), r, r)))
    if (det(R) < 0) R[, r] <- -R[, r]
    R
  }
  X2 <- t(rot()) %*% Xref
  X3 <- t(rot()) %*% Xref
  colnames(X2) <- effects
  colnames(X3) <- effects

  adat <- AlignmentData(
    list(s1 = Xref, s2 = X2, s3 = X3),
    design = list(K = K, effects = effects)
  )

  res <- fit_alignment(adat, method = "kprocrustes", reference = "consensus", cv = "loso")
  expect_s4_class(res, "AlignmentResult")
})

test_that("kprocrustes apply_alignment fits a new subject", {
  neuralign:::.register_kprocrustes()

  set.seed(4)
  q <- 6
  r <- 3
  effects <- paste0("e", seq_len(q))
  K <- diag(q)
  dimnames(K) <- list(effects, effects)

  Uref <- k_orthonormalize(matrix(rnorm(q * r), q, r), K)
  Xref <- t(Uref)
  colnames(Xref) <- effects

  R2 <- qr.Q(qr(matrix(rnorm(r * r), r, r)))
  if (det(R2) < 0) R2[, r] <- -R2[, r]
  X2 <- t(R2) %*% Xref
  colnames(X2) <- effects

  model_res <- fit_alignment(
    AlignmentData(list(s1 = Xref, s2 = X2), design = list(K = K, effects = effects)),
    method = "kprocrustes",
    reference = "s1"
  )

  R3 <- qr.Q(qr(matrix(rnorm(r * r), r, r)))
  if (det(R3) < 0) R3[, r] <- -R3[, r]
  X3 <- t(R3) %*% Xref
  colnames(X3) <- effects

  new_res <- apply_alignment(model_res, AlignmentData(list(s3 = X3), design = list(K = K, effects = effects)))
  new_model <- get_model(new_res)
  expect_true("s3" %in% names(new_model@transforms))
  expect_equal(new_model@transforms$s3 %*% X3, new_model@reference_data, tolerance = 1e-6)
})

test_that("kprocrustes supports partial effect overlap via zero-fill harmonization", {
  neuralign:::.register_kprocrustes()

  set.seed(5)
  q <- 6
  r <- 3
  effects <- paste0("e", seq_len(q))
  K <- diag(q)
  dimnames(K) <- list(effects, effects)

  X1 <- matrix(rnorm(r * q), r, q)
  colnames(X1) <- effects

  X2 <- matrix(rnorm(r * (q - 1L)), r, q - 1L)
  colnames(X2) <- effects[-1L]

  X3 <- matrix(rnorm(r * (q - 1L)), r, q - 1L)
  colnames(X3) <- effects[-q]

  adat <- AlignmentData(
    list(s1 = X1, s2 = X2, s3 = X3),
    design = list(K = K, effects = effects)
  )

  res <- fit_alignment(adat, method = "kprocrustes", reference = "consensus")
  model <- get_model(res)

  expect_equal(dim(model@reference_data), c(r, q))
  for (subj in names(model@transforms)) {
    expect_equal(dim(model@transforms[[subj]]), c(r, r))
  }
})

test_that("kprocrustes reflection control forces det(transform) >= 0", {
  neuralign:::.register_kprocrustes()

  set.seed(6)
  q <- 6
  r <- 3
  effects <- paste0("e", seq_len(q))
  K <- diag(q)
  dimnames(K) <- list(effects, effects)

  Uref <- k_orthonormalize(matrix(rnorm(q * r), q, r), K)
  Xref <- t(Uref)
  colnames(Xref) <- effects

  Rneg <- qr.Q(qr(matrix(rnorm(r * r), r, r)))
  if (det(Rneg) > 0) Rneg[, 1] <- -Rneg[, 1]
  expect_lt(det(Rneg), 0)

  X2 <- t(Rneg) %*% Xref
  colnames(X2) <- effects

  res <- fit_alignment(
    AlignmentData(list(s1 = Xref, s2 = X2), design = list(K = K, effects = effects)),
    method = "kprocrustes",
    reference = "s1",
    allow_reflection = FALSE
  )
  model <- get_model(res)
  expect_gte(det(model@transforms$s2), 0)
})
