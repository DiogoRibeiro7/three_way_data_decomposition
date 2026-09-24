#' Construct a Centered Centroid Basis for Tucker3 SCR
#'
#' Project a candidate group-mode basis into the contrast space orthogonal to
#' the mixture-probability vector and return an orthonormal basis for that
#' subspace.
#'
#' For mixture probabilities `p` and centroid basis `A`, the constraint
#' `t(p) %*% A = 0` separates the grand mean from the centroid deviations.
#'
#' @param raw_basis Numeric G-by-P candidate basis.
#' @param probabilities Strictly positive mixture probabilities of length G.
#' @param tolerance Numerical rank tolerance.
#' @return A G-by-P matrix with orthonormal columns and
#'   `crossprod(probabilities, A) = 0` up to numerical tolerance.
#' @export
scr_centroid_basis <- function(
  raw_basis,
  probabilities,
  tolerance = 1e-10
) {
  validate_centroid_basis_inputs(
    raw_basis = raw_basis,
    probabilities = probabilities,
    tolerance = tolerance
  )

  centered <- sweep(
    raw_basis,
    MARGIN = 2L,
    STATS = as.numeric(crossprod(probabilities, raw_basis)),
    FUN = "-"
  )

  decomposition <- qr(centered, tol = tolerance)
  rank <- ncol(raw_basis)

  if (decomposition$rank < rank) {
    stop(
      "raw_basis is rank deficient after probability-weighted centering."
    )
  }

  basis <- qr.Q(decomposition, complete = FALSE)[, seq_len(rank), drop = FALSE]

  if (max(abs(crossprod(probabilities, basis))) > sqrt(tolerance)) {
    stop("failed to construct a centered centroid basis.")
  }

  basis
}


#' Reconstruct Tucker3 SCR Component Means
#'
#' Construct the G component means from the Tucker3 extension described by
#' Rocci, Vichi, and Ranalli:
#'
#' deqn{
#' \mu_{gjk} = \mu_{jk} +
#' \sum_p \sum_q \sum_r
#' a_{gp} b_{jq} c_{kr} \eta_{pqr}.
#' }
#'
#' The group-mode basis must satisfy `t(probabilities) %*% centroid_basis = 0`.
#' Consequently, the probability-weighted mean of the returned component means
#' equals `grand_mean`.
#'
#' @param grand_mean Numeric J-by-K matrix or length J*K vector.
#' @param centroid_basis Numeric G-by-P centroid loading matrix.
#' @param variable_basis Numeric J-by-Q variable loading matrix.
#' @param occasion_basis Numeric K-by-R occasion loading matrix.
#' @param core Numeric array with dimensions P-by-Q-by-R.
#' @param probabilities Strictly positive mixture probabilities of length G.
#' @param tolerance Numerical tolerance for the centering constraint.
#' @return A G-by-(J*K) matrix of vectorised component means, using the same
#'   column-major variable-within-occasion ordering as the S3 implementation.
#' @export
scr_tucker3_means <- function(
  grand_mean,
  centroid_basis,
  variable_basis,
  occasion_basis,
  core,
  probabilities,
  tolerance = 1e-8
) {
  dimensions <- validate_tucker3_mean_structure(
    grand_mean = grand_mean,
    centroid_basis = centroid_basis,
    variable_basis = variable_basis,
    occasion_basis = occasion_basis,
    core = core,
    probabilities = probabilities,
    tolerance = tolerance
  )

  grand_vector <- if (is.matrix(grand_mean)) {
    as.vector(grand_mean)
  } else {
    as.numeric(grand_mean)
  }

  core_mode1 <- matrix(
    core,
    nrow = dimensions$P,
    ncol = dimensions$Q * dimensions$R
  )

  reduced_scores <- centroid_basis %*% core_mode1
  loading_basis <- kronecker(occasion_basis, variable_basis)
  deviations <- reduced_scores %*% t(loading_basis)

  sweep(
    deviations,
    MARGIN = 2L,
    STATS = grand_vector,
    FUN = "+"
  )
}


#' Parameter Count for the Tucker3 SCR Extension
#'
#' Compute the effective number of free parameters for the Gaussian SCR model
#' with Tucker3 reduction of the centroid, variable, and occasion modes and a
#' separable covariance.
#'
#' The centroid loading subspace lives in the (G - 1)-dimensional contrast
#' space induced by the grand-mean constraint. Its Grassmann dimension is
#' `P * (G - 1 - P)`.
#'
#' When `P = G - 1`, the count reduces exactly to the published S3/Tucker2
#' parameter count.
#'
#' @param groups Number of mixture components G.
#' @param variables Number of observed variables J.
#' @param occasions Number of occasions K.
#' @param centroid_rank Centroid-mode rank P.
#' @param variable_rank Variable-mode rank Q.
#' @param occasion_rank Occasion-mode rank R.
#' @param details Logical; when `TRUE`, return the parameter-count breakdown.
#' @return An integer parameter count, or a named integer vector when
#'   `details = TRUE`.
#' @export
scr_tucker3_parameter_count <- function(
  groups,
  variables,
  occasions,
  centroid_rank,
  variable_rank,
  occasion_rank,
  details = FALSE
) {
  validate_tucker3_dimensions(
    groups = groups,
    variables = variables,
    occasions = occasions,
    centroid_rank = centroid_rank,
    variable_rank = variable_rank,
    occasion_rank = occasion_rank
  )

  if (!is.logical(details) || length(details) != 1L || is.na(details)) {
    stop("details must be TRUE or FALSE.")
  }

  components <- c(
    mixing_proportions = groups - 1L,
    grand_mean = variables * occasions,
    tucker_core = centroid_rank * variable_rank * occasion_rank,
    centroid_subspace = centroid_rank * (groups - 1L - centroid_rank),
    variable_subspace = variable_rank * (variables - variable_rank),
    occasion_subspace = occasion_rank * (occasions - occasion_rank),
    variable_covariance = variables * (variables + 1L) / 2L,
    occasion_covariance = occasions * (occasions + 1L) / 2L,
    kronecker_scale_constraint = -1L
  )

  storage.mode(components) <- "integer"

  if (details) {
    return(c(components, total = sum(components)))
  }

  as.integer(sum(components))
}


validate_tucker3_mean_structure <- function(
  grand_mean,
  centroid_basis,
  variable_basis,
  occasion_basis,
  core,
  probabilities,
  tolerance
) {
  if (
    !is.matrix(centroid_basis) ||
      !is.numeric(centroid_basis) ||
      anyNA(centroid_basis) ||
      any(!is.finite(centroid_basis))
  ) {
    stop("centroid_basis must be a finite numeric matrix.")
  }

  if (
    !is.matrix(variable_basis) ||
      !is.numeric(variable_basis) ||
      anyNA(variable_basis) ||
      any(!is.finite(variable_basis))
  ) {
    stop("variable_basis must be a finite numeric matrix.")
  }

  if (
    !is.matrix(occasion_basis) ||
      !is.numeric(occasion_basis) ||
      anyNA(occasion_basis) ||
      any(!is.finite(occasion_basis))
  ) {
    stop("occasion_basis must be a finite numeric matrix.")
  }

  G <- nrow(centroid_basis)
  P <- ncol(centroid_basis)
  J <- nrow(variable_basis)
  Q <- ncol(variable_basis)
  K <- nrow(occasion_basis)
  R <- ncol(occasion_basis)

  validate_tucker3_dimensions(
    groups = G,
    variables = J,
    occasions = K,
    centroid_rank = P,
    variable_rank = Q,
    occasion_rank = R
  )

  if (
    !is.numeric(probabilities) ||
      length(probabilities) != G ||
      anyNA(probabilities) ||
      any(!is.finite(probabilities)) ||
      any(probabilities <= 0)
  ) {
    stop(
      "probabilities must contain one strictly positive finite value per component."
    )
  }

  if (abs(sum(probabilities) - 1) > tolerance) {
    stop("probabilities must sum to one.")
  }

  if (
    !is.numeric(tolerance) ||
      length(tolerance) != 1L ||
      is.na(tolerance) ||
      !is.finite(tolerance) ||
      tolerance <= 0
  ) {
    stop("tolerance must be a positive finite number.")
  }

  if (max(abs(crossprod(probabilities, centroid_basis))) > tolerance) {
    stop(
      "centroid_basis must satisfy the probability-weighted centering constraint."
    )
  }

  if (qr(centroid_basis)$rank < P) {
    stop("centroid_basis must have full column rank.")
  }
  if (qr(variable_basis)$rank < Q) {
    stop("variable_basis must have full column rank.")
  }
  if (qr(occasion_basis)$rank < R) {
    stop("occasion_basis must have full column rank.")
  }

  if (
    !is.array(core) ||
      !is.numeric(core) ||
      !identical(dim(core), c(P, Q, R)) ||
      anyNA(core) ||
      any(!is.finite(core))
  ) {
    stop("core must be a finite numeric P-by-Q-by-R array.")
  }

  if (is.matrix(grand_mean)) {
    if (!identical(dim(grand_mean), c(J, K))) {
      stop("grand_mean matrix must have dimensions J by K.")
    }
  } else if (
    !is.numeric(grand_mean) ||
      length(grand_mean) != J * K
  ) {
    stop("grand_mean must be a J-by-K matrix or a length J*K vector.")
  }

  if (anyNA(grand_mean) || any(!is.finite(grand_mean))) {
    stop("grand_mean must contain only finite values.")
  }

  list(G = G, P = P, J = J, Q = Q, K = K, R = R)
}


validate_tucker3_dimensions <- function(
  groups,
  variables,
  occasions,
  centroid_rank,
  variable_rank,
  occasion_rank
) {
  values <- c(
    groups = groups,
    variables = variables,
    occasions = occasions,
    centroid_rank = centroid_rank,
    variable_rank = variable_rank,
    occasion_rank = occasion_rank
  )

  if (
    anyNA(values) ||
      any(!is.finite(values)) ||
      any(values <= 0) ||
      any(values %% 1 != 0)
  ) {
    stop("all Tucker3 dimensions and ranks must be positive integers.")
  }

  if (groups < 2L) {
    stop("groups must be at least two.")
  }
  if (centroid_rank > groups - 1L) {
    stop("centroid_rank cannot exceed groups - 1.")
  }
  if (variable_rank > variables) {
    stop("variable_rank cannot exceed variables.")
  }
  if (occasion_rank > occasions) {
    stop("occasion_rank cannot exceed occasions.")
  }

  invisible(TRUE)
}


validate_centroid_basis_inputs <- function(
  raw_basis,
  probabilities,
  tolerance
) {
  if (
    !is.matrix(raw_basis) ||
      !is.numeric(raw_basis) ||
      nrow(raw_basis) < 2L ||
      ncol(raw_basis) < 1L ||
      ncol(raw_basis) > nrow(raw_basis) - 1L ||
      anyNA(raw_basis) ||
      any(!is.finite(raw_basis))
  ) {
    stop(
      "raw_basis must be a finite G-by-P numeric matrix with P <= G - 1."
    )
  }

  if (
    !is.numeric(probabilities) ||
      length(probabilities) != nrow(raw_basis) ||
      anyNA(probabilities) ||
      any(!is.finite(probabilities)) ||
      any(probabilities <= 0)
  ) {
    stop(
      "probabilities must contain one strictly positive finite value per component."
    )
  }

  if (abs(sum(probabilities) - 1) > tolerance) {
    stop("probabilities must sum to one.")
  }

  if (
    !is.numeric(tolerance) ||
      length(tolerance) != 1L ||
      is.na(tolerance) ||
      !is.finite(tolerance) ||
      tolerance <= 0
  ) {
    stop("tolerance must be a positive finite number.")
  }

  invisible(TRUE)
}
