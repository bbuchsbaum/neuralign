#' Fit an Alignment Model
#'
#' Fit an alignment model to multi-subject data. This is the main entry point
#' for alignment in neuralign. The fitted model can be saved and applied to
#' new data using \code{\link{apply_alignment}}.
#'
#' @family alignment_workflow
#'
#' @details
#' \code{fit_alignment()} does not perform any centering or scaling of inputs by
#' default. If you want preprocessing, apply it explicitly to your matrices or
#' feature blocks (e.g., \code{\link{preprocess_alignment_data}} or
#' \code{\link{preprocess_matrix}}) before fitting.
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
#' @param cv_folds Either an integer number of folds for k-fold CV (used when
#'   \code{cv = "kfold"}), or a fold specification list from
#'   \code{\link{create_cv_folds}} (or compatible) to control the exact train/test
#'   splits.
#' @param train_idx Optional integer indices specifying which subjects to use
#'   for fitting. If NULL, all subjects are used. This enables manual CV schemes.
#' @param obs_labels Optional shared observation labels. If \code{data} is a
#'   list, this is passed through to \code{\link{AlignmentData}} via
#'   \code{\link{as_alignment_data}}. If \code{data} is already an
#'   \code{AlignmentData} and has no observation labels, they will be set.
#' @param compute_quality Logical; if TRUE, compute quality metrics.
#' @param return_aligned Logical; if TRUE (default), apply transforms and store
#'   aligned data in the result. Set to FALSE to return only the model (saves
#'   memory in pipeline workflows where aligned data is not needed immediately).
#' @param return_fold_transforms Logical; if TRUE and using observation-axis
#'   fold specs (see \code{\link{create_obs_folds}}), retain per-fold transforms
#'   in \code{result@cv_info$transforms_by_fold}. Default FALSE.
#' @param return_resample_artifacts Logical; if TRUE for operator-returning
#'   subject CV, retain each fold's fitted model plus aligned analysis and
#'   assessment matrices. These artifacts are required to construct an
#'   [AlignedResampleSet]. Default FALSE.
#' @param restrict_to_identified Logical; if TRUE, canonicalize orthogonal
#'   operators so their action on the orthogonal complement of the identified
#'   subspace is deterministic. This is useful when correspondence signals are
#'   rank-deficient (e.g., transform_dim > numeric_rank). Defaults to FALSE.
#'   Only supported for methods that return orthogonal operators.
#' @param restrict_tol Relative tolerance for numeric rank determination used by
#'   \code{\link{canonicalize_orthogonal_operator}}.
#' @param ... Additional arguments passed to the method's fit function.
#'
#' @return An AlignmentResult object containing the fitted model and
#'   (optionally) aligned data and quality metrics.
#'
#' @examples
#' # Create example data
#' set.seed(1)
#' data_list <- list(
#'   "sub-01" = matrix(rnorm(30), 6, 5),
#'   "sub-02" = matrix(rnorm(30), 6, 5),
#'   "sub-03" = matrix(rnorm(30), 6, 5)
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
                          obs_labels = NULL,
                          compute_quality = TRUE,
                          return_aligned = TRUE,
                          return_fold_transforms = FALSE,
                          return_resample_artifacts = FALSE,
                          restrict_to_identified = FALSE,
                          restrict_tol = sqrt(.Machine$double.eps),
                          ...) {
  cv <- match.arg(cv)
  if (!is.logical(return_resample_artifacts) ||
      length(return_resample_artifacts) != 1L ||
      is.na(return_resample_artifacts)) {
    stop("'return_resample_artifacts' must be TRUE or FALSE", call. = FALSE)
  }

  # Coerce data to AlignmentData if needed
  if (inherits(data, "AlignmentData")) {
    data <- .set_obs_labels_if_missing(data, obs_labels)
  } else {
    data <- as_alignment_data(data, obs_labels = obs_labels)
  }

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
  caps <- aligner$capabilities %||% list()
  returns <- caps$returns %||% "operator"

  restrict_to_identified <- isTRUE(restrict_to_identified)
  if (restrict_to_identified) {
    if (!identical(returns, "operator")) {
      stop(
        sprintf(
          "restrict_to_identified=TRUE requires operator-returning methods; method '%s' returns '%s'",
          method, returns
        ),
        call. = FALSE
      )
    }
    transform_type <- caps$transform_type %||% NA_character_
    if (!is.character(transform_type) || length(transform_type) != 1L ||
        !transform_type %in% c("orthogonal")) {
      stop(
        sprintf(
          "restrict_to_identified=TRUE is only supported for orthogonal operators; method '%s' has transform_type='%s'",
          method, as.character(transform_type)
        ),
        call. = FALSE
      )
    }
    if (!is.numeric(restrict_tol) || length(restrict_tol) != 1L || !is.finite(restrict_tol) ||
        restrict_tol <= 0 || restrict_tol >= 1) {
      stop("'restrict_tol' must be a single number in (0, 1)", call. = FALSE)
    }
  }

  .validate_reference_for_aligner(reference, data, method, caps)

  # Validate data dimensions based on method capabilities
  validate_alignment_data(
    data,
    check_features = isTRUE(caps$requires_shared_features %||% TRUE),
    check_observations = isTRUE(caps$requires_shared_observations %||% FALSE)
  )

  # Route based on CV strategy
  # Check for explicit fold spec first (overrides cv="none" default)
  if (.is_cv_folds_spec(cv_folds)) {
    if (identical(cv_folds$axis, "observation")) {
      if (identical(returns, "embedding")) {
        if (isTRUE(return_resample_artifacts)) {
          stop(
            "return_resample_artifacts is not yet supported for embedding-returning methods",
            call. = FALSE
          )
        }
        if (isTRUE(return_fold_transforms)) {
          stop(
            "return_fold_transforms is not supported for embedding-returning methods",
            call. = FALSE
          )
        }
        result <- .fit_cv_obs_folds_embedding(
          data = data,
          aligner = aligner,
          reference = reference,
          cv_folds = cv_folds,
          compute_quality = compute_quality,
          return_aligned = return_aligned,
          ...
        )
      } else {
        if (isTRUE(return_resample_artifacts)) {
          stop(
            "Observation-fold artifacts are retained by ",
            "run_obs_crossfit_from_data(); fit_alignment() currently retains ",
            "resample artifacts only for subject-axis CV",
            call. = FALSE
          )
        }
        result <- .fit_cv_obs_folds(
          data = data,
          aligner = aligner,
          reference = reference,
          cv_folds = cv_folds,
          compute_quality = compute_quality,
          return_aligned = return_aligned,
          return_fold_transforms = return_fold_transforms,
          restrict_to_identified = restrict_to_identified,
          restrict_tol = restrict_tol,
          ...
        )
      }
    } else {
      validate_cv_setup(cv_folds, reference = reference)
      if (identical(returns, "embedding")) {
        if (isTRUE(return_resample_artifacts)) {
          stop(
            "return_resample_artifacts is not yet supported for embedding-returning methods",
            call. = FALSE
          )
        }
        result <- .fit_cv_folds_embedding(
          data = data,
          aligner = aligner,
          reference = reference,
          cv_folds = cv_folds,
          compute_quality = compute_quality,
          return_aligned = return_aligned,
          ...
        )
      } else {
        result <- .fit_cv_folds(
          data = data,
          aligner = aligner,
          reference = reference,
          cv_folds = cv_folds,
          compute_quality = compute_quality,
          return_aligned = return_aligned,
          return_resample_artifacts = return_resample_artifacts,
          restrict_to_identified = restrict_to_identified,
          restrict_tol = restrict_tol,
          ...
        )
      }
    }
  } else if (cv == "none") {
    if (isTRUE(return_resample_artifacts)) {
      stop(
        "return_resample_artifacts=TRUE requires subject-axis cross-validation",
        call. = FALSE
      )
    }
    result <- .fit_single(data, aligner, reference, train_idx,
      compute_quality = compute_quality,
      return_aligned = return_aligned,
      restrict_to_identified = restrict_to_identified,
      restrict_tol = restrict_tol,
      ...
    )
  } else if (cv == "loso") {
    result <- .fit_cv_loso(
      data = data,
      aligner = aligner,
      reference = reference,
      compute_quality = compute_quality,
      return_aligned = return_aligned,
      return_resample_artifacts = return_resample_artifacts,
      restrict_to_identified = restrict_to_identified,
      restrict_tol = restrict_tol,
      ...
    )
  } else if (cv == "kfold") {
    result <- .fit_cv_kfold(
      data = data,
      aligner = aligner,
      reference = reference,
      k = cv_folds,
      compute_quality = compute_quality,
      return_aligned = return_aligned,
      return_resample_artifacts = return_resample_artifacts,
      restrict_to_identified = restrict_to_identified,
      restrict_tol = restrict_tol,
      ...
    )
  }

  result
}
