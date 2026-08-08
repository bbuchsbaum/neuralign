#' neuralign: General-Purpose Alignment Infrastructure
#'
#' `neuralign` provides a domain-agnostic framework for fitting and applying
#' alignment transforms across a set of subjects (or other grouped datasets).
#' It is designed to centralize common alignment operations so that downstream
#' packages can reuse consistent implementations of alignment, cross-validation,
#' feature harmonization, and provenance tracking.
#'
#' ## Core conventions
#'
#' **Data matrices.** By convention, subject data matrices are oriented as
#' **features × observations**. If your data are in the opposite orientation,
#' transpose them before creating an [AlignmentData] object.
#'
#' **Transforms.** Most aligners return linear operators that are applied by
#' **left-multiplication** (e.g., `Q %*% X`). When a right-multiplication
#' convention is needed (e.g., `X %*% Q`), use a convention-aware helper such as
#' [procrustes_rotation()] or transpose the problem explicitly.
#'
#' ## Key components
#'
#' - [AlignmentData]: container for per-subject matrices (and optional observation
#'   labels, guidance, and other metadata).
#' - [fit_alignment()] / [apply_alignment()]: fit an alignment model and apply it
#'   to new subjects or held-out data.
#' - [AlignmentModel] / [AlignmentResult]: standardized objects for storing
#'   fitted transforms, aligned outputs, and quality/CV metadata.
#' - [AlignedStudy] / [SharedFeatureSpace]: analysis-facing aligned data with
#'   observation metadata, shared-space identity, lineage, and safety records
#'   (see [align_study()], [as_aligned_study()]).
#' - Feature blocks: represent heterogeneous alignment signals and harmonize them
#'   across subjects via shared feature names (see [alignment_feature_block()],
#'   [harmonize_feature_blocks()], [stack_feature_blocks()]).
#'
#' ## Cross-validation
#'
#' `neuralign` supports both subject-axis and observation-axis cross-validation.
#' Use [create_cv_folds()] for subject CV and [create_obs_folds()] for
#' observation-axis folds (e.g., run-based or blocked-time splits).
#'
#' ## Preprocessing
#'
#' `neuralign` does not perform centering, scaling, or other preprocessing by
#' default. Apply any desired preprocessing to your matrices before fitting.
#'
#' @seealso [AlignmentData], [fit_alignment()], [apply_alignment()],
#'   [alignment_feature_block()], [available_aligners()]
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
