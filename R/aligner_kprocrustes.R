#' K-metric Procrustes Alignment
#'
#' Alignment in effect space using a kernel metric `K`, producing orthogonal
#' left-multiply operators (r x r) that align each subject's rows (components)
#' to a shared reference.
#'
#' @name aligner_kprocrustes
NULL

.kprocrustes_default_effects <- function(q) {
  paste0("effect", seq_len(q))
}

.kprocrustes_extract_design <- function(design) {
  if (is.null(design) || length(design) == 0L) {
    stop("kprocrustes: AlignmentData@design must provide at least K and effects", call. = FALSE)
  }

  if (.is_matrixish(design)) {
    K <- as.matrix(design)
    if (nrow(K) != ncol(K)) stop("kprocrustes: design K must be square", call. = FALSE)
    effects <- colnames(K) %||% rownames(K) %||% .kprocrustes_default_effects(nrow(K))
    effects <- as.character(effects)
    if (is.null(rownames(K)) || is.null(colnames(K))) {
      dimnames(K) <- list(effects, effects)
    }
    return(list(K = K, effects = effects))
  }

  if (!is.list(design)) {
    stop("kprocrustes: AlignmentData@design must be a list (or a K matrix)", call. = FALSE)
  }

  K <- design$K %||% design$design_kernel %||% design$kernel
  if (is.list(K) && !is.null(K$K)) {
    K <- K$K
  }
  if (!.is_matrixish(K)) {
    stop("kprocrustes: design$K must be a matrix (or list with $K)", call. = FALSE)
  }
  K <- as.matrix(K)
  if (nrow(K) != ncol(K)) {
    stop("kprocrustes: design$K must be square", call. = FALSE)
  }

  q <- nrow(K)
  effects <- design$effects %||% colnames(K) %||% rownames(K)
  if (is.null(effects)) {
    effects <- .kprocrustes_default_effects(q)
  }
  effects <- as.character(effects)
  if (length(effects) != q) {
    stop("kprocrustes: design$effects must have length q (number of effects)", call. = FALSE)
  }

  # Ensure K has dimnames consistent with effects.
  if (is.null(rownames(K)) || is.null(colnames(K))) {
    dimnames(K) <- list(effects, effects)
  } else if (!identical(rownames(K), effects) || !identical(colnames(K), effects)) {
    if (!setequal(rownames(K), effects) || !setequal(colnames(K), effects)) {
      stop("kprocrustes: K dimnames must match design$effects", call. = FALSE)
    }
    K <- K[effects, effects, drop = FALSE]
  }

  list(K = K, effects = effects)
}

.kprocrustes_validate_data_list <- function(data_list, effects, context) {
  if (!is.list(data_list) || !length(data_list) || is.null(names(data_list))) {
    stop(context, ": expected a named list of matrices", call. = FALSE)
  }

  bad <- names(data_list)[!vapply(data_list, .is_matrixish, logical(1))]
  if (length(bad)) {
    stop(context, ": all subjects must be matrices; invalid: ", paste(bad, collapse = ", "), call. = FALSE)
  }

  # Enforce shared feature dimension (rows).
  dims_feat <- vapply(data_list, function(m) nrow(as.matrix(m)), integer(1))
  if (length(unique(dims_feat)) != 1L) {
    stop(context, ": all subjects must have the same number of rows (components). Got: ",
      paste(sprintf("%s=%d", names(dims_feat), dims_feat), collapse = ", "),
      call. = FALSE
    )
  }

  q <- length(effects)
  for (subj in names(data_list)) {
    m <- as.matrix(data_list[[subj]])
    cn <- colnames(m)
    if (!is.null(cn)) {
      cn <- as.character(cn)
      if (anyDuplicated(cn)) {
        stop(context, ": duplicated effect names for '", subj, "'", call. = FALSE)
      }
      if (!all(cn %in% effects)) {
        bad_eff <- setdiff(cn, effects)
        stop(context, ": unknown effect names for '", subj, "': ", paste(bad_eff, collapse = ", "), call. = FALSE)
      }
      # If full coverage with different ordering, reorder to effects.
      if (length(cn) == q && setequal(cn, effects) && !identical(cn, effects)) {
        m <- m[, effects, drop = FALSE]
      }
    } else {
      if (ncol(m) != q) {
        stop(
          context, ": missing colnames for '", subj, "' and ncol != q; cannot map effects",
          call. = FALSE
        )
      }
      colnames(m) <- effects
    }
    data_list[[subj]] <- m
  }

  data_list
}

.kprocrustes_reference_to_U <- function(reference,
                                        U_list,
                                        train_subjects,
                                        K,
                                        allow_reflection) {
  if (is.character(reference) && length(reference) == 1L) {
    if (reference == "consensus") {
      U_train <- U_list[train_subjects]
      cons <- k_consensus_basis(U_train, K, allow_reflection = allow_reflection)
      return(cons$U)
    }
    if (reference %in% names(U_list)) {
      if (!reference %in% train_subjects) {
        stop(
          "kprocrustes: reference subject must be in training set for fit (got: ",
          reference, ")",
          call. = FALSE
        )
      }
      return(U_list[[reference]])
    }
    stop("kprocrustes: unknown reference subject id: ", reference, call. = FALSE)
  }

  if (.is_matrixish(reference)) {
    ref <- as.matrix(reference)
    q <- nrow(K)
    if (nrow(ref) == q) {
      # q x r basis
      return(ref)
    }
    if (ncol(ref) == q) {
      # r x q "X" representation
      return(t(ref))
    }
    stop("kprocrustes: template reference has incompatible dims; expected qxr (U) or rxq (t(U))", call. = FALSE)
  }

  stop("kprocrustes: unsupported reference; expected 'consensus', a subject id, or a template matrix", call. = FALSE)
}

.kprocrustes_fit <- function(data,
                             reference = "consensus",
                             train_idx = NULL,
                             allow_reflection = TRUE,
                             min_overlap = 1L,
                             warn_sparse_below = 0.5,
                             ...) {
  if (is.null(train_idx)) train_idx <- seq_along(data@subjects)
  subjects <- data@subjects
  train_subjects <- subjects[train_idx]

  design <- .kprocrustes_extract_design(data@design)
  data_list <- .kprocrustes_validate_data_list(get_data_list(data), design$effects, "kprocrustes fit")

  min_overlap <- as.integer(min_overlap)
  if (!is.finite(min_overlap) || min_overlap < 1L) {
    stop("kprocrustes: 'min_overlap' must be a positive integer", call. = FALSE)
  }
  eff_counts <- vapply(data_list, function(m) ncol(as.matrix(m)), integer(1))
  too_small <- names(eff_counts)[eff_counts < min_overlap]
  if (length(too_small) > 0) {
    stop(
      sprintf("kprocrustes: subject(s) have < min_overlap=%d effects: %s", min_overlap, paste(too_small, collapse = ", ")),
      call. = FALSE
    )
  }

  harm <- harmonize_union_fill(
    mats = data_list,
    axis = "cols",
    union_ids = design$effects,
    fill = 0,
    min_coverage = 0L,
    warn_sparse_below = warn_sparse_below
  )

  U_list <- lapply(harm$mats, function(X) k_orthonormalize(t(as.matrix(X)), design$K))

  U_ref <- .kprocrustes_reference_to_U(
    reference = reference,
    U_list = U_list,
    train_subjects = train_subjects,
    K = design$K,
    allow_reflection = allow_reflection
  )
  U_ref <- k_orthonormalize(U_ref, design$K)

  r <- ncol(U_ref)
  transforms <- lapply(names(U_list), function(subj) {
    if (is.character(reference) && length(reference) == 1L && subj == reference) {
      diag(r)
    } else {
      pr <- k_procrustes(U_ref, U_list[[subj]], design$K, allow_reflection = allow_reflection)
      t(pr$R) # left-multiply operator on X = t(U)
    }
  })
  names(transforms) <- names(U_list)

  space_spec <- data@space %||% "kprocrustes_effect_basis"

  reference_data <- t(U_ref) # r x q
  colnames(reference_data) <- design$effects

  list(
    transforms = transforms,
    reference_data = reference_data,
    space_from = space_spec,
    space_to = space_spec,
    method_state = list(
      U_ref = U_ref,
      K = design$K,
      effects = design$effects,
      allow_reflection = isTRUE(allow_reflection),
      min_overlap = min_overlap,
      warn_sparse_below = warn_sparse_below,
      obs_mask = harm$obs_mask
    )
  )
}

.kprocrustes_apply <- function(fit_result, new_data, ...) {
  if (!inherits(new_data, "AlignmentData")) {
    stop("'new_data' must be an AlignmentData object", call. = FALSE)
  }
  if (length(new_data@subjects) != 1L) {
    stop("kprocrustes apply supports a single new subject at a time", call. = FALSE)
  }

  st <- fit_result$method_state %||% list()
  K <- st$K
  effects <- st$effects
  U_ref <- st$U_ref
  allow_reflection <- isTRUE(st$allow_reflection %||% TRUE)
  min_overlap <- as.integer(st$min_overlap %||% 1L)
  warn_sparse_below <- st$warn_sparse_below %||% 0.5

  if (is.null(K) || is.null(effects) || is.null(U_ref)) {
    stop("kprocrustes apply: missing method_state (expected K/effects/U_ref)", call. = FALSE)
  }

  data_list <- .kprocrustes_validate_data_list(get_data_list(new_data), effects, "kprocrustes apply")
  eff_counts <- vapply(data_list, function(m) ncol(as.matrix(m)), integer(1))
  if (eff_counts[[1L]] < min_overlap) {
    stop(
      sprintf("kprocrustes apply: subject has < min_overlap=%d effects", min_overlap),
      call. = FALSE
    )
  }

  harm <- harmonize_union_fill(
    mats = data_list,
    axis = "cols",
    union_ids = effects,
    fill = 0,
    min_coverage = 0L,
    warn_sparse_below = warn_sparse_below
  )

  subj_id <- new_data@subjects[[1L]]
  X <- as.matrix(harm$mats[[1L]])
  U <- k_orthonormalize(t(X), K)
  pr <- k_procrustes(U_ref, U, K, allow_reflection = allow_reflection)
  op <- t(pr$R)

  list(
    transforms = stats::setNames(list(op), subj_id),
    reference_data = fit_result$reference_data,
    space_from = fit_result$space_from,
    space_to = fit_result$space_to,
    method_state = fit_result$method_state %||% list()
  )
}

.kprocrustes_capabilities <- list(
  supports_cv = TRUE,
  cv_axes = c("subject"),
  needs_geometry = FALSE,
  needs_design = TRUE,
  requires_shared_features = TRUE,
  requires_shared_observations = FALSE,
  returns_invertible = TRUE,
  transform_type = "orthogonal",
  mass_preserving = FALSE,
  returns = "operator",
  supports_new_subject = TRUE,
  supports_new_data = TRUE,
  reference_types = c("subject", "consensus", "template")
)

.register_kprocrustes <- function() {
  register_aligner(
    name = "kprocrustes",
    fit_fn = .kprocrustes_fit,
    apply_fn = .kprocrustes_apply,
    capabilities = .kprocrustes_capabilities,
    package = "neuralign",
    description = "K-metric Procrustes alignment (design-kernel-weighted)",
    version = as.character(utils::packageVersion("neuralign"))
  )
}
