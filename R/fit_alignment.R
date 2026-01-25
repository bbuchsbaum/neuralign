#' Fit an Alignment Model
#'
#' Fit an alignment model to multi-subject data. This is the main entry point
#' for alignment in neuralign. The fitted model can be saved and applied to
#' new data using \code{\link{apply_alignment}}.
#'
#' @param data AlignmentData object, or a named list of matrices that will
#'   be coerced to AlignmentData.
#' @param method Character string specifying the alignment method. Use
#'   \code{\link{available_aligners}} to see registered methods.
#' @param reference How to select the reference for alignment:
#'   \itemize{
#'     \item "medoid" - Select the medoid subject (default)
#'     \item "consensus" - Compute a consensus/mean reference
#'     \item Subject ID - Use a specific subject as reference
#'     \item Matrix - Use an external template matrix
#'   }
#' @param cv Cross-validation strategy:
#'   \itemize{
#'     \item "none" - No cross-validation (default)
#'     \item "loso" - Leave-one-subject-out
#'     \item "kfold" - K-fold cross-validation (specify k via cv_folds)
#'   }
#' @param cv_folds Number of folds for k-fold CV. Ignored if cv != "kfold".
#' @param train_idx Optional integer indices specifying which subjects to use
#'   for fitting. If NULL, all subjects are used. This enables manual CV schemes.
#' @param compute_quality Logical; if TRUE, compute quality metrics.
#' @param ... Additional arguments passed to the method's fit function.
#'
#' @return An AlignmentResult object containing the fitted model and
#'   (optionally) aligned data and quality metrics.
#'
#' @examples
#' \dontrun{
#' # Create example data
#' data_list <- list(
#'   "sub-01" = matrix(rnorm(100*50), 100, 50),
#'   "sub-02" = matrix(rnorm(100*50), 100, 50),
#'   "sub-03" = matrix(rnorm(100*50), 100, 50)
#' )
#' adat <- AlignmentData(data_list)
#'
#' # Fit alignment with default settings
#' result <- fit_alignment(adat, method = "procrustes")
#'
#' # Fit with specific reference subject
#' result <- fit_alignment(adat, method = "procrustes", reference = "sub-01")
#'
#' # Fit with leave-one-subject-out CV
#' result <- fit_alignment(adat, method = "procrustes", cv = "loso")
#' }
#'
#' @seealso \code{\link{apply_alignment}}, \code{\link{AlignmentData}},
#'   \code{\link{available_aligners}}
#'
#' @export
fit_alignment <- function(data,
                          method = "procrustes",
                          reference = "medoid",
                          cv = c("none", "loso", "kfold"),
                          cv_folds = 5,
                          train_idx = NULL,
                          compute_quality = TRUE,
                          ...) {
  cv <- match.arg(cv)

  # Coerce data to AlignmentData if needed
  if (!inherits(data, "AlignmentData")) {
    data <- as_alignment_data(data)
  }

  # Validate data
  validate_alignment_data(data, check_features = TRUE)

  # Try to load aligner if not registered
  if (!is_aligner_registered(method)) {
    if (!.try_autoload_aligner(method)) {
      stop(sprintf(
        "Unknown method '%s'. Available methods: %s",
        method,
        paste(available_aligners(), collapse = ", ")
      ))
    }
  }

  # Validate method requirements
  .validate_aligner_requirements(method, data)

  # Get aligner
  aligner <- get_aligner(method)

  # Route based on CV strategy
  if (cv == "none") {
    result <- .fit_single(data, aligner, reference, train_idx,
      compute_quality = compute_quality, ...
    )
  } else if (cv == "loso") {
    result <- .fit_cv_loso(data, aligner, reference,
      compute_quality = compute_quality, ...
    )
  } else if (cv == "kfold") {
    result <- .fit_cv_kfold(data, aligner, reference, cv_folds,
      compute_quality = compute_quality, ...
    )
  }

  result
}


#' Internal: Single Fit (No CV)
#' @keywords internal
.fit_single <- function(data, aligner, reference, train_idx,
                        compute_quality, ...) {
  # Determine training subjects
  if (is.null(train_idx)) {
    train_idx <- seq_along(data@subjects)
  }
  train_subjects <- data@subjects[train_idx]

  # Resolve reference
  ref_resolved <- .resolve_reference(data, reference, train_idx)

  # Call the fit function
  fit_result <- aligner$fit_fn(
    data = data,
    reference = ref_resolved$reference,
    train_idx = train_idx,
    ...
  )

  # Build AlignmentModel
  model <- AlignmentModel(
    transforms = fit_result$transforms,
    reference = ref_resolved$reference_spec,
    reference_data = fit_result$reference_data,
    method = aligner$name,
    space_from = fit_result$space_from,
    space_to = fit_result$space_to,
    params = list(...),
    method_state = fit_result$method_state %||% list(),
    train_subjects = train_subjects
  )

  # Apply alignment to get aligned data
  aligned <- .apply_transforms(model, data)

  # Compute quality metrics if requested
  quality <- list()
  if (compute_quality) {
    quality <- .compute_basic_quality(data, aligned, model)
  }

  AlignmentResult(
    model = model,
    aligned = aligned,
    quality = quality,
    cv_info = list(method = "none")
  )
}


#' Internal: Leave-One-Subject-Out CV
#' @keywords internal
.fit_cv_loso <- function(data, aligner, reference, compute_quality, ...) {
  n_subjects <- length(data@subjects)
  subjects <- data@subjects

  # Check CV support
  if (!isTRUE(aligner$capabilities$supports_cv)) {
    warning(sprintf(
      "Method '%s' may not fully support CV; results may have leakage",
      aligner$name
    ))
  }

  # Storage for per-fold results
  all_transforms <- list()
  all_aligned <- list()
  fold_quality <- vector("list", n_subjects)

  # Iterate over folds
  for (i in seq_len(n_subjects)) {
    test_idx <- i
    train_idx <- setdiff(seq_len(n_subjects), i)
    test_subj <- subjects[test_idx]

    # Resolve reference using only training subjects
    ref_resolved <- .resolve_reference(data, reference, train_idx)

    # Fit on training subjects
    fit_result <- aligner$fit_fn(
      data = data,
      reference = ref_resolved$reference,
      train_idx = train_idx,
      ...
    )

    # Apply to held-out subject
    test_transform <- .fit_new_subject(
      aligner, fit_result, data, test_idx, ref_resolved$reference
    )

    all_transforms[[test_subj]] <- test_transform

    # Apply transform
    test_data <- get_subject_data(data, test_subj)
    all_aligned[[test_subj]] <- test_transform %*% test_data
  }

  # Use final fold's fit for model (or could use all-subject fit)
  # For CV, we typically want the all-data fit for the model
  ref_resolved <- .resolve_reference(data, reference, seq_len(n_subjects))
  fit_result <- aligner$fit_fn(
    data = data,
    reference = ref_resolved$reference,
    train_idx = seq_len(n_subjects),
    ...
  )

  model <- AlignmentModel(
    transforms = all_transforms,
    reference = ref_resolved$reference_spec,
    reference_data = fit_result$reference_data,
    method = aligner$name,
    space_from = fit_result$space_from,
    space_to = fit_result$space_to,
    params = list(...),
    method_state = fit_result$method_state %||% list(),
    train_subjects = subjects
  )

  # Quality metrics
  quality <- list()
  if (compute_quality) {
    quality <- .compute_basic_quality(data, all_aligned, model)
  }

  AlignmentResult(
    model = model,
    aligned = all_aligned,
    quality = quality,
    cv_info = list(
      method = "loso",
      n_folds = n_subjects,
      fold_assignments = setNames(seq_len(n_subjects), subjects)
    )
  )
}


#' Internal: K-Fold CV
#' @keywords internal
.fit_cv_kfold <- function(data, aligner, reference, k, compute_quality, ...) {
  n_subjects <- length(data@subjects)
  subjects <- data@subjects

  if (k > n_subjects) {
    warning(sprintf("k=%d > n_subjects=%d; using LOSO instead", k, n_subjects))
    return(.fit_cv_loso(data, aligner, reference, compute_quality, ...))
  }

  # Create fold assignments
  fold_assignments <- sample(rep(seq_len(k), length.out = n_subjects))
  names(fold_assignments) <- subjects

  # Check CV support
  if (!isTRUE(aligner$capabilities$supports_cv)) {
    warning(sprintf(
      "Method '%s' may not fully support CV; results may have leakage",
      aligner$name
    ))
  }

  all_transforms <- list()
  all_aligned <- list()

  # Iterate over folds
  for (fold in seq_len(k)) {
    test_idx <- which(fold_assignments == fold)
    train_idx <- which(fold_assignments != fold)
    test_subjects <- subjects[test_idx]

    # Resolve reference using only training subjects
    ref_resolved <- .resolve_reference(data, reference, train_idx)

    # Fit on training subjects
    fit_result <- aligner$fit_fn(
      data = data,
      reference = ref_resolved$reference,
      train_idx = train_idx,
      ...
    )

    # Apply to held-out subjects
    for (test_subj in test_subjects) {
      test_i <- match(test_subj, subjects)
      test_transform <- .fit_new_subject(
        aligner, fit_result, data, test_i, ref_resolved$reference
      )
      all_transforms[[test_subj]] <- test_transform

      test_data <- get_subject_data(data, test_subj)
      all_aligned[[test_subj]] <- test_transform %*% test_data
    }
  }

  # All-data fit for model
  ref_resolved <- .resolve_reference(data, reference, seq_len(n_subjects))
  fit_result <- aligner$fit_fn(
    data = data,
    reference = ref_resolved$reference,
    train_idx = seq_len(n_subjects),
    ...
  )

  model <- AlignmentModel(
    transforms = all_transforms,
    reference = ref_resolved$reference_spec,
    reference_data = fit_result$reference_data,
    method = aligner$name,
    space_from = fit_result$space_from,
    space_to = fit_result$space_to,
    params = list(...),
    method_state = fit_result$method_state %||% list(),
    train_subjects = subjects
  )

  quality <- list()
  if (compute_quality) {
    quality <- .compute_basic_quality(data, all_aligned, model)
  }

  AlignmentResult(
    model = model,
    aligned = all_aligned,
    quality = quality,
    cv_info = list(
      method = "kfold",
      n_folds = k,
      fold_assignments = fold_assignments
    )
  )
}


#' Internal: Resolve Reference Specification
#' @keywords internal
.resolve_reference <- function(data, reference, train_idx) {
  subjects <- data@subjects
  train_subjects <- subjects[train_idx]

  if (is.character(reference) && length(reference) == 1) {
    if (reference == "medoid") {
      # Select medoid from training subjects only
      ref_subj <- select_reference(data[train_idx], method = "medoid")
      return(list(
        reference = ref_subj,
        reference_spec = ref_subj
      ))
    } else if (reference == "consensus") {
      return(list(
        reference = "consensus",
        reference_spec = "consensus"
      ))
    } else if (reference == "centroid") {
      ref_subj <- select_reference(data[train_idx], method = "centroid")
      return(list(
        reference = ref_subj,
        reference_spec = ref_subj
      ))
    } else {
      # Assume it's a subject ID
      if (!reference %in% train_subjects) {
        warning(sprintf(
          "Reference subject '%s' not in training set; potential leakage",
          reference
        ))
      }
      return(list(
        reference = reference,
        reference_spec = reference
      ))
    }
  } else if (is.matrix(reference)) {
    # External template
    return(list(
      reference = reference,
      reference_spec = "template"
    ))
  } else {
    stop("Invalid reference specification")
  }
}


#' Internal: Fit Transform for New Subject
#' @keywords internal
.fit_new_subject <- function(aligner, fit_result, data, subject_idx,
                             reference) {
  subject <- data@subjects[subject_idx]

  # If aligner has an apply_fn, use it
  if (!is.null(aligner$apply_fn)) {
    apply_result <- aligner$apply_fn(
      fit_result = fit_result,
      new_data = data[subject_idx]
    )
    return(apply_result$transforms[[1]])
  }

  # Otherwise, use the method's fit_fn with single subject
  # Use the reference_data (actual matrix) from fit_result, not the reference spec
  # This avoids issues when reference is a subject ID not in the single-subject data
  ref_to_use <- fit_result$reference_data
  if (is.null(ref_to_use)) {
    ref_to_use <- reference
  }

  single_fit <- aligner$fit_fn(
    data = data,
    reference = ref_to_use,
    train_idx = subject_idx
  )

  single_fit$transforms[[subject]]
}


#' Internal: Apply Transforms to Data
#' @keywords internal
.apply_transforms <- function(model, data) {
  aligned <- list()

  for (subj in names(model@transforms)) {
    if (subj %in% data@subjects) {
      transform <- model@transforms[[subj]]
      subj_data <- get_subject_data(data, subj)

      # Left-multiply: Y = A %*% X
      aligned[[subj]] <- transform %*% subj_data
    }
  }

  aligned
}


#' Internal: Compute Basic Quality Metrics
#' @keywords internal
.compute_basic_quality <- function(data, aligned, model) {
  metrics <- list()

  # Mean pairwise correlation of aligned data
  if (length(aligned) >= 2) {
    pairs <- combn(names(aligned), 2)
    cors <- apply(pairs, 2, function(p) {
      x <- aligned[[p[1]]]
      y <- aligned[[p[2]]]
      mean(diag(cor(t(x), t(y))), na.rm = TRUE)
    })
    metrics$mean_pairwise_correlation <- mean(cors, na.rm = TRUE)
    metrics$pairwise_correlations <- setNames(
      cors,
      apply(pairs, 2, paste, collapse = "-")
    )
  }

  metrics
}
