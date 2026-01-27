#' Procrustes Graph Synchronization
#'
#' Align subjects when observation labels only partially overlap. Computes
#' pairwise Procrustes rotations on overlaps and synchronizes them into a
#' globally-consistent set of subject-to-reference transforms.
#'
#' @name aligner_procrustes_graph
NULL

.connected_components <- function(adj) {
  n <- nrow(adj)
  seen <- rep(FALSE, n)
  comps <- vector("list", 0)

  for (start in seq_len(n)) {
    if (seen[[start]]) next
    stack <- start
    comp <- integer(0)
    seen[[start]] <- TRUE
    while (length(stack) > 0) {
      v <- stack[[length(stack)]]
      stack <- stack[-length(stack)]
      comp <- c(comp, v)
      neigh <- which(adj[v, ])
      new <- neigh[!seen[neigh]]
      if (length(new) > 0) {
        seen[new] <- TRUE
        stack <- c(stack, new)
      }
    }
    comps[[length(comps) + 1L]] <- comp
  }

  comps
}

.project_to_orthogonal <- function(M, reflection = FALSE) {
  sv <- svd(M)
  Q <- sv$u %*% t(sv$v)
  if (!isTRUE(reflection) && det(Q) < 0) {
    sv$u[, ncol(sv$u)] <- -sv$u[, ncol(sv$u)]
    Q <- sv$u %*% t(sv$v)
  }
  Q
}

.as_obs_labels_list <- function(data) {
  labs <- data@obs_labels
  subjects <- data@subjects
  data_list <- get_data_list(data)

  if (is.null(labs)) {
    stop(
      "procrustes_graph requires obs_labels (atomic shared labels or a per-subject named list)",
      call. = FALSE
    )
  }

  if (is.atomic(labs) || is.factor(labs)) {
    labs <- as.character(labs)
    out <- rep(list(labs), length(subjects))
    names(out) <- subjects
    return(out)
  }

  if (is.list(labs)) {
    if (is.null(names(labs))) {
      stop("obs_labels list must be named by subject", call. = FALSE)
    }
    missing <- setdiff(subjects, names(labs))
    if (length(missing) > 0) {
      stop(
        sprintf("obs_labels list is missing subjects: %s", paste(missing, collapse = ", ")),
        call. = FALSE
      )
    }
    out <- lapply(subjects, function(s) as.character(labs[[s]]))
    names(out) <- subjects

    # Basic sanity: length must match ncol for each subject
    for (s in subjects) {
      x <- data_list[[s]]
      if (!.is_matrixish(x)) x <- as.matrix(x)
      if (length(out[[s]]) != ncol(x)) {
        stop(
          sprintf("obs_labels length mismatch for subject '%s'", s),
          call. = FALSE
        )
      }
    }
    return(out)
  }

  stop("obs_labels must be NULL, an atomic vector, or a named list", call. = FALSE)
}

.build_overlap_adjacency <- function(obs_labels_by_subject, min_overlap) {
  subjects <- names(obs_labels_by_subject)
  n <- length(subjects)
  adj <- matrix(FALSE, n, n, dimnames = list(subjects, subjects))
  overlap <- matrix(0L, n, n, dimnames = list(subjects, subjects))

  for (i in seq_len(n - 1L)) {
    li <- obs_labels_by_subject[[i]]
    for (j in (i + 1L):n) {
      lj <- obs_labels_by_subject[[j]]
      o <- length(intersect(li, lj))
      overlap[i, j] <- overlap[j, i] <- o
      if (o >= min_overlap) {
        adj[i, j] <- adj[j, i] <- TRUE
      }
    }
  }

  diag(adj) <- FALSE
  list(adj = adj, overlap = overlap)
}

.synchronize_rotations_spectral <- function(rotations, weights, n_nodes, d, reflection) {
  # rotations: list keyed by "i-j" with matrices mapping j -> i (d x d)
  # weights: numeric vector keyed by same "i-j" names
  S <- matrix(0, n_nodes * d, n_nodes * d)

  deg <- numeric(n_nodes)
  for (nm in names(weights)) {
    idx <- strsplit(nm, "-", fixed = TRUE)[[1]]
    i <- as.integer(idx[[1]])
    j <- as.integer(idx[[2]])
    w <- weights[[nm]]
    deg[[i]] <- deg[[i]] + w
    deg[[j]] <- deg[[j]] + w
  }
  deg[deg == 0] <- 1

  for (nm in names(rotations)) {
    idx <- strsplit(nm, "-", fixed = TRUE)[[1]]
    i <- as.integer(idx[[1]])
    j <- as.integer(idx[[2]])
    Rji <- rotations[[nm]] # maps node j -> node i
    w <- weights[[nm]] / sqrt(deg[[i]] * deg[[j]])

    ii <- ((i - 1L) * d + 1L):(i * d)
    jj <- ((j - 1L) * d + 1L):(j * d)

    S[ii, jj] <- w * Rji
    S[jj, ii] <- w * t(Rji)
  }

  eig <- eigen(S, symmetric = TRUE)
  V <- eig$vectors[, seq_len(d), drop = FALSE]

  G <- vector("list", n_nodes)
  for (i in seq_len(n_nodes)) {
    ii <- ((i - 1L) * d + 1L):(i * d)
    block <- V[ii, , drop = FALSE]
    # Important: do NOT enforce det(Q) >= 0 per-subject. Relative rotations
    # are invariant to a shared global reflection, but per-node enforcement
    # can break consistency when the leading eigenspace is sign-indeterminate.
    G[[i]] <- .project_to_orthogonal(block, reflection = TRUE)
  }
  G
}

#' Procrustes Graph Fit Function
#' @keywords internal
.procrustes_graph_fit <- function(data,
                                  reference = NULL,
                                  train_idx = NULL,
                                  min_overlap = 2L,
                                  weight = c("overlap", "uniform"),
                                  reflection = FALSE,
                                  ...) {
  weight <- match.arg(weight)
  min_overlap <- as.integer(min_overlap)
  if (!is.finite(min_overlap) || min_overlap < 1L) {
    stop("'min_overlap' must be a positive integer", call. = FALSE)
  }

  if (is.null(train_idx)) train_idx <- seq_along(data@subjects)

  # Fit on training subjects only to support subject-axis CV without leakage.
  train_data <- data[train_idx]
  subjects <- train_data@subjects
  data_list <- get_data_list(train_data)
  d <- nrow(as.matrix(data_list[[1L]]))

  obs_labels_by_subject <- .as_obs_labels_list(train_data)
  graph <- .build_overlap_adjacency(obs_labels_by_subject, min_overlap = min_overlap)
  adj <- graph$adj
  overlap <- graph$overlap

  if (!any(adj)) {
    stop(
      sprintf("No edges in overlap graph with min_overlap=%d", min_overlap),
      call. = FALSE
    )
  }

  comps <- .connected_components(adj)
  if (length(comps) > 1L) {
    comps_named <- lapply(comps, function(idx) subjects[idx])
    msg <- paste0(
      "Overlap graph is disconnected (", length(comps), " components). ",
      "Increase stimulus overlap, lower min_overlap, or provide additional anchors. Components: ",
      paste(vapply(comps_named, function(x) paste0("{", paste(x, collapse = ","), "}"), character(1)), collapse = " ")
    )
    stop(msg, call. = FALSE)
  }

  # Resolve reference subject
  if (is.null(reference)) {
    reference <- subjects[[1L]]
  }
  if (!is.character(reference) || length(reference) != 1L || !reference %in% subjects) {
    stop("procrustes_graph requires 'reference' to be a subject id", call. = FALSE)
  }
  ref_idx <- match(reference, subjects)

  rotations <- list()
  weights <- numeric(0)

  for (i in seq_len(length(subjects) - 1L)) {
    Xi <- as.matrix(data_list[[i]])
    for (j in (i + 1L):length(subjects)) {
      if (!adj[i, j]) next
      Xj <- as.matrix(data_list[[j]])

      res <- procrustes_rotation(
        source = Xj,
        target = Xi,
        convention = "left",
        scale = FALSE,
        reflection = reflection,
        obs_labels_source = obs_labels_by_subject[[j]],
        obs_labels_target = obs_labels_by_subject[[i]],
        min_overlap = min_overlap
      )
      key <- paste0(i, "-", j)
      rotations[[key]] <- res$Q
      weights[[key]] <- if (weight == "overlap") overlap[i, j] else 1
    }
  }

  if (length(rotations) < 1L) {
    stop(
      sprintf("No edges in overlap graph with min_overlap=%d", min_overlap),
      call. = FALSE
    )
  }

  G <- .synchronize_rotations_spectral(
    rotations = rotations,
    weights = weights,
    n_nodes = length(subjects),
    d = d,
    reflection = reflection
  )

  # Gauge-fix so that reference subject has identity transform.
  G_ref <- G[[ref_idx]]
  transforms_train <- lapply(G, function(Gi) G_ref %*% t(Gi))
  names(transforms_train) <- subjects

  # Fit transforms for any held-out subjects by aligning directly to the
  # reference subject using label overlap.
  transforms <- transforms_train
  heldout_subjects <- setdiff(data@subjects, subjects)
  if (length(heldout_subjects) > 0) {
    all_labels <- .as_obs_labels_list(data)
    ref_labels <- obs_labels_by_subject[[reference]]
    ref_data <- get_subject_data(train_data, reference)
    for (subj in heldout_subjects) {
      X_new <- as.matrix(get_subject_data(data, subj))
      labs_new <- all_labels[[subj]]
      res <- procrustes_rotation(
        source = X_new,
        target = as.matrix(ref_data),
        convention = "left",
        scale = FALSE,
        reflection = reflection,
        obs_labels_source = labs_new,
        obs_labels_target = ref_labels,
        min_overlap = min_overlap
      )
      transforms[[subj]] <- res$Q
    }
  }

  list(
    transforms = transforms,
    reference_data = get_subject_data(data, reference),
    space_from = data@space,
    space_to = data@space,
    method_state = list(
      reference = reference,
      reference_obs_labels = obs_labels_by_subject[[reference]],
      min_overlap = min_overlap,
      weight = weight,
      reflection = reflection
    )
  )
}

#' Procrustes Graph Apply Function (New Subjects)
#' @keywords internal
.procrustes_graph_apply <- function(fit_result, new_data, ...) {
  if (!inherits(new_data, "AlignmentData")) {
    stop("'new_data' must be an AlignmentData object", call. = FALSE)
  }
  if (length(new_data@subjects) != 1L) {
    stop("procrustes_graph apply supports a single new subject at a time", call. = FALSE)
  }

  ms <- fit_result$method_state %||% list()
  ref_subj <- ms$reference %||% NULL
  if (is.null(ref_subj) || !is.character(ref_subj) || length(ref_subj) != 1L) {
    stop("procrustes_graph apply requires method_state$reference", call. = FALSE)
  }
  min_overlap <- as.integer(ms$min_overlap %||% 2L)
  reflection <- isTRUE(ms$reflection %||% FALSE)
  ref_labels <- ms$reference_obs_labels %||% NULL

  ref_data <- fit_result$reference_data
  if (is.null(ref_data) || !.is_matrixish(ref_data)) {
    stop("procrustes_graph apply requires reference_data in fit_result", call. = FALSE)
  }

  subj <- new_data@subjects[[1L]]
  X_new <- as.matrix(get_subject_data(new_data, subj))
  labs_new <- .as_obs_labels_list(new_data)[[subj]]

  res <- procrustes_rotation(
    source = X_new,
    target = as.matrix(ref_data),
    convention = "left",
    scale = FALSE,
    reflection = reflection,
    obs_labels_source = labs_new,
    obs_labels_target = ref_labels,
    min_overlap = min_overlap
  )

  list(
    transforms = setNames(list(res$Q), subj),
    reference_data = ref_data,
    space_from = fit_result$space_from,
    space_to = fit_result$space_to,
    method_state = fit_result$method_state %||% list()
  )
}

#' Procrustes Graph Capabilities
#' @keywords internal
.procrustes_graph_capabilities <- list(
  supports_cv = TRUE,
  cv_axes = c("subject", "observation"),
  needs_geometry = FALSE,
  needs_design = FALSE,
  requires_shared_features = TRUE,
  requires_shared_observations = FALSE,
  returns_invertible = TRUE,
  transform_type = "orthogonal",
  mass_preserving = TRUE,
  returns = "operator",
  supports_new_subject = TRUE,
  supports_new_data = TRUE,
  reference_types = c("subject")
)

#' Register Procrustes Graph Aligner
#' @keywords internal
.register_procrustes_graph <- function() {
  register_aligner(
    name = "procrustes_graph",
    fit_fn = .procrustes_graph_fit,
    apply_fn = .procrustes_graph_apply,
    capabilities = .procrustes_graph_capabilities,
    package = "neuralign",
    description = "Procrustes synchronization on an overlap graph",
    version = as.character(utils::packageVersion("neuralign"))
  )
}
