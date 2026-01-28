#' Bridge Utilities for manifoldalign Integration
#'
#' Internal helpers for converting neuralign data to manifoldalign inputs and
#' extracting operator transforms in a way that is (a) mathematically coherent
#' and (b) computationally tractable.
#'
#' Key principle:
#' - For manifoldalign's multiblock methods (KEMA/GPCA/coupled_diag/lowrank),
#'   we treat *observations* as samples: `x = t(X)` where neuralign stores
#'   `X` as `(features x observations)`.
#'   These methods naturally produce per-domain loadings `v_i (features x k)`.
#'   neuralign operators become `A_i = t(v_i)` (shape `k x features`), mapping
#'   into a shared latent space with `A_i %*% X_i -> (k x observations)`.
#' - For graph correspondence methods (CONE/GRASP), we treat *features* as nodes
#'   and use the returned assignment to build sparse permutation-style operators
#'   into a fixed reference feature space.
#'
#' @name manifoldalign_bridge
#' @keywords internal
NULL


.ma_require_manifoldalign <- function(method) {
  if (!requireNamespace("manifoldalign", quietly = TRUE)) {
    stop(
      sprintf("Package 'manifoldalign' is required for %s alignment.", method),
      call. = FALSE
    )
  }
}


.ma_require_shared_obs_labels <- function(data, method) {
  labs <- data@obs_labels
  if (is.null(labs)) {
    stop(sprintf("Method '%s' requires AlignmentData@obs_labels", method), call. = FALSE)
  }
  if (!(is.atomic(labs) || is.factor(labs))) {
    stop(
      sprintf(
        "Method '%s' currently requires shared (atomic/factor) obs_labels; per-subject obs_labels lists are not supported yet",
        method
      ),
      call. = FALSE
    )
  }
  if (any(is.na(labs))) {
    stop(sprintf("Method '%s' requires obs_labels without NA", method), call. = FALSE)
  }
  labs
}

.ma_hyperdesign_meta <- function(domains) {
  domain_names <- names(domains)
  nr <- vapply(domains, function(d) nrow(d$x), integer(1))
  nc <- vapply(domains, function(d) ncol(d$x), integer(1))
  ny <- vapply(domains, function(d) {
    if (is.null(d$design)) 0L else ncol(d$design)
  }, integer(1))

  row_end <- cumsum(nr)
  row_start <- c(1L, utils::head(row_end, -1L) + 1L)
  col_end <- cumsum(nc)
  col_start <- c(1L, utils::head(col_end, -1L) + 1L)

  data.frame(
    block = seq_along(domains),
    block_name = domain_names,
    nr = as.integer(nr),
    nxvar = as.integer(nc),
    nyvar = as.integer(ny),
    row_start = as.integer(row_start),
    row_end = as.integer(row_end),
    col_start = as.integer(col_start),
    col_end = as.integer(col_end),
    stringsAsFactors = FALSE
  )
}


.ma_build_hyperdesign_obs <- function(data, train_idx = NULL, labels, label_name = "label") {
  if (is.null(train_idx)) train_idx <- seq_along(data@subjects)
  train_data <- data[train_idx]
  data_list <- get_data_list(train_data)

  domains <- lapply(names(data_list), function(nm) {
    X <- as.matrix(data_list[[nm]])
    Xo <- t(X) # observations x features
    if (length(labels) != nrow(Xo)) {
      stop(
        sprintf(
          "obs_labels length mismatch: expected %d (n_obs), got %d",
          nrow(Xo), length(labels)
        ),
        call. = FALSE
      )
    }
    list(
      x = Xo,
      design = data.frame(stats::setNames(list(labels), label_name))
    )
  })
  names(domains) <- names(data_list)
  hdes <- .ma_hyperdesign_meta(domains)

  structure(
    domains,
    hdes = hdes,
    common_vars = label_name,
    class = c("neuralign_hyperdesign", "hyperdesign", "list")
  )
}


.ma_build_hyperdesign_features <- function(data, train_idx = NULL) {
  if (is.null(train_idx)) train_idx <- seq_along(data@subjects)
  train_data <- data[train_idx]
  data_list <- get_data_list(train_data)

  domains <- lapply(names(data_list), function(nm) {
    X <- as.matrix(data_list[[nm]]) # features x observations
    list(x = X)
  })
  names(domains) <- names(data_list)
  hdes <- .ma_hyperdesign_meta(domains)

  structure(
    domains,
    hdes = hdes,
    common_vars = character(0),
    class = c("neuralign_hyperdesign", "hyperdesign", "list")
  )
}

.ma_build_pair_hyperdesign_features <- function(target, source,
                                                target_name = "target",
                                                source_name = "source") {
  domains <- list()
  domains[[target_name]] <- list(x = as.matrix(target))
  domains[[source_name]] <- list(x = as.matrix(source))

  hdes <- .ma_hyperdesign_meta(domains)

  structure(
    domains,
    hdes = hdes,
    common_vars = character(0),
    class = c("neuralign_hyperdesign", "hyperdesign", "list")
  )
}


.ma_resolve_reference_spec <- function(data, reference, method, allow_template = FALSE) {
  if (is.character(reference) && length(reference) == 1L && reference %in% c("medoid", "centroid")) {
    ref_name <- select_reference(data, method = reference)
    ref_data <- get_subject_data(data, ref_name)
    return(list(name = ref_name, data = ref_data))
  }

  if (.is_matrixish(reference)) {
    if (!allow_template) {
      stop(sprintf("Method '%s' does not support template matrix references", method), call. = FALSE)
    }
    return(list(name = "template", data = as.matrix(reference)))
  }

  if (is.character(reference) && length(reference) == 1L) {
    if (!reference %in% data@subjects) {
      stop(sprintf("Unknown reference subject '%s'", reference), call. = FALSE)
    }
    return(list(name = reference, data = get_subject_data(data, reference)))
  }

  stop("Unsupported reference type", call. = FALSE)
}


.ma_mbp_split_loadings <- function(mbp, domain_names, feature_counts) {
  if (is.null(mbp$v)) stop("manifoldalign result missing $v loadings", call. = FALSE)
  v <- as.matrix(mbp$v)

  idx <- mbp$block_indices
  if (is.null(idx)) {
    # Fallback: contiguous blocks by feature counts
    cum <- cumsum(c(0L, feature_counts))
    idx <- lapply(seq_along(feature_counts), function(i) seq.int(cum[i] + 1L, cum[i + 1L]))
  } else if (is.matrix(idx) && ncol(idx) == 2L) {
    idx <- lapply(seq_len(nrow(idx)), function(i) seq.int(idx[i, 1], idx[i, 2]))
  } else if (!is.list(idx)) {
    stop("Unexpected block_indices format in manifoldalign result", call. = FALSE)
  }

  if (length(idx) != length(domain_names)) {
    stop("block_indices length mismatch with number of domains", call. = FALSE)
  }

  v_blocks <- lapply(seq_along(domain_names), function(i) {
    vi <- v[idx[[i]], , drop = FALSE]
    if (nrow(vi) != feature_counts[[i]]) {
      stop(
        sprintf(
          "Loadings block row mismatch for '%s': expected %d, got %d",
          domain_names[[i]], feature_counts[[i]], nrow(vi)
        ),
        call. = FALSE
      )
    }
    vi
  })
  names(v_blocks) <- domain_names
  v_blocks
}


.ma_projection_transforms_from_loadings <- function(v_blocks) {
  setNames(lapply(names(v_blocks), function(nm) {
    t(v_blocks[[nm]])
  }), names(v_blocks))
}


.ma_reference_scores <- function(A_ref, X_ref) {
  A_ref %*% X_ref
}


.ma_match_obs_indices <- function(obs_labels_ref, obs_labels_new, min_overlap = 2L) {
  obs_labels_ref <- as.character(obs_labels_ref)
  obs_labels_new <- as.character(obs_labels_new)

  common <- intersect(obs_labels_ref, obs_labels_new)
  if (length(common) < min_overlap) {
    stop(
      sprintf("Insufficient obs_labels overlap (need >= %d, got %d)", min_overlap, length(common)),
      call. = FALSE
    )
  }
  list(
    ref = match(common, obs_labels_ref),
    new = match(common, obs_labels_new)
  )
}


.ma_ridge_map_to_reference_scores <- function(X, Z_ref, lambda = 1e-2) {
  X <- as.matrix(X)
  Z_ref <- as.matrix(Z_ref)
  if (ncol(X) != ncol(Z_ref)) {
    stop("X and Z_ref must have the same number of observations (columns)", call. = FALSE)
  }
  n <- ncol(X)
  XtX <- crossprod(X) + lambda * diag(n)
  W <- solve(XtX, t(X)) # n x p
  Z_ref %*% W          # k x p
}


.ma_get_single_subject_obs_labels <- function(data) {
  labs <- data@obs_labels
  if (is.null(labs)) return(NULL)
  if (is.atomic(labs) || is.factor(labs)) return(labs)
  if (!is.list(labs)) return(NULL)

  if (length(data@subjects) != 1L) {
    stop(
      "Expected AlignmentData with exactly one subject when resolving per-subject obs_labels",
      call. = FALSE
    )
  }
  subj <- data@subjects[[1L]]

  if (!is.null(names(labs)) && subj %in% names(labs)) {
    return(labs[[subj]])
  }
  if (length(labs) == 1L) {
    return(labs[[1L]])
  }

  stop(
    sprintf(
      "AlignmentData@obs_labels is a list but does not contain subject '%s' and has length %d",
      subj, length(labs)
    ),
    call. = FALSE
  )
}


.ma_assignment_to_operator <- function(assignment, n_target, n_source) {
  if (is.null(assignment)) stop("Missing assignment for graph alignment", call. = FALSE)
  a <- as.integer(assignment)
  if (length(a) != n_target) {
    stop(
      sprintf("Assignment length mismatch: expected %d, got %d", n_target, length(a)),
      call. = FALSE
    )
  }
  ok <- !is.na(a) & a >= 1L & a <= n_source
  Matrix::sparseMatrix(
    i = which(ok),
    j = a[ok],
    x = rep(1, sum(ok)),
    dims = c(n_target, n_source)
  )
}
