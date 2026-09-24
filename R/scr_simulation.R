#' Generate Data from a Heteroscedastic Gaussian Mixture
#'
#' Generate observations from the Gaussian mixture used in the SCR simulation
#' study of Rocci, Vichi, and Ranalli. The function also returns the posterior
#' component-membership probabilities under the generating model.
#'
#' @param n Positive integer number of observations.
#' @param probabilities Strictly positive mixing probabilities summing to one.
#' @param means Numeric matrix with one component mean per row.
#' @param covariances Either a list of covariance matrices, one per component,
#'   or the MATLAB reference representation obtained by vertically stacking the
#'   covariance matrices.
#' @return A list with:
#'   - `X`: generated observation matrix;
#'   - `U`: posterior component-membership probabilities;
#'   - `z`: latent component labels used to generate the observations.
#' @export
generate_gaussian_mixture <- function(
  n,
  probabilities,
  means,
  covariances
) {
  validate_sample_size(n)
  validate_mixture_means(means)

  groups <- nrow(means)
  dimension <- ncol(means)

  validate_mixture_probabilities(probabilities, groups)
  covariance_list <- normalize_covariances(
    covariances = covariances,
    groups = groups,
    dimension = dimension
  )

  cumulative_probability <- cumsum(probabilities)
  uniforms <- stats::runif(n)
  z <- findInterval(
    uniforms,
    cumulative_probability,
    left.open = TRUE
  ) + 1L
  z <- as.integer(z)

  X <- matrix(0, nrow = n, ncol = dimension)

  for (g in seq_len(groups)) {
    index <- which(z == g)
    group_size <- length(index)

    if (group_size == 0L) {
      next
    }

    standard_normal <- matrix(
      stats::rnorm(group_size * dimension),
      nrow = group_size,
      ncol = dimension
    )

    X[index, ] <- sweep(
      standard_normal %*% chol(covariance_list[[g]]),
      MARGIN = 2L,
      STATS = means[g, ],
      FUN = "+"
    )
  }

  log_kernel <- matrix(0, nrow = n, ncol = groups)

  for (g in seq_len(groups)) {
    decomposition <- svd(covariance_list[[g]])
    singular_values <- decomposition$d

    centered <- sweep(
      X,
      MARGIN = 2L,
      STATS = means[g, ],
      FUN = "-"
    )
    whitened <- centered %*%
      decomposition$v %*%
      diag(1 / sqrt(singular_values), nrow = dimension)

    log_kernel[, g] <- -0.5 * sum(log(singular_values)) -
      0.5 * rowSums(whitened^2)
  }

  # Preserve the underflow guard used by the MATLAB reference implementation.
  kernel <- exp(pmax(log_kernel, -700))
  weighted_kernel <- sweep(
    kernel,
    MARGIN = 2L,
    STATS = probabilities,
    FUN = "*"
  )
  U <- weighted_kernel / rowSums(weighted_kernel)

  list(X = X, U = U, z = z)
}


#' Legacy SCR Gaussian-Mixture Generator
#'
#' Compatibility wrapper for the original MATLAB/R function name.
#'
#' @param n See `generate_gaussian_mixture()`.
#' @param p Mixing probabilities.
#' @param M Component-mean matrix.
#' @param Sig Component covariance matrices, either as an R list or as the
#'   vertically stacked MATLAB representation.
#' @return See `generate_gaussian_mixture()`.
#' @export
genmixhet <- function(n, p, M, Sig) {
  generate_gaussian_mixture(
    n = n,
    probabilities = p,
    means = M,
    covariances = Sig
  )
}


validate_sample_size <- function(n) {
  if (
    !is.numeric(n) ||
      length(n) != 1L ||
      is.na(n) ||
      !is.finite(n) ||
      n <= 0 ||
      n %% 1 != 0
  ) {
    stop("n must be a positive integer.")
  }

  invisible(TRUE)
}


validate_mixture_means <- function(means) {
  if (
    !is.matrix(means) ||
      !is.numeric(means) ||
      nrow(means) == 0L ||
      ncol(means) == 0L
  ) {
    stop("means must be a non-empty numeric matrix.")
  }

  if (anyNA(means) || any(!is.finite(means))) {
    stop("means must contain only finite values.")
  }

  invisible(TRUE)
}


validate_mixture_probabilities <- function(probabilities, groups) {
  if (
    !is.numeric(probabilities) ||
      length(probabilities) != groups ||
      anyNA(probabilities) ||
      any(!is.finite(probabilities)) ||
      any(probabilities <= 0)
  ) {
    stop(
      "probabilities must contain one strictly positive finite value per component."
    )
  }

  if (abs(sum(probabilities) - 1) > 1e-10) {
    stop("probabilities must sum to one.")
  }

  invisible(TRUE)
}


normalize_covariances <- function(covariances, groups, dimension) {
  if (is.list(covariances)) {
    if (length(covariances) != groups) {
      stop("covariances must contain one matrix per component.")
    }

    covariance_list <- covariances
  } else if (is.matrix(covariances) && is.numeric(covariances)) {
    expected_dimension <- c(groups * dimension, dimension)

    if (!identical(dim(covariances), expected_dimension)) {
      stop(
        "stacked covariances must have dimensions G * J by J."
      )
    }

    covariance_list <- lapply(
      seq_len(groups),
      function(g) {
        rows <- seq.int(
          from = (g - 1L) * dimension + 1L,
          length.out = dimension
        )
        covariances[rows, , drop = FALSE]
      }
    )
  } else {
    stop(
      "covariances must be a list of matrices or a stacked numeric matrix."
    )
  }

  for (g in seq_len(groups)) {
    covariance <- covariance_list[[g]]

    if (
      !is.matrix(covariance) ||
        !is.numeric(covariance) ||
        !identical(dim(covariance), c(dimension, dimension))
    ) {
      stop("each covariance matrix must be numeric and J by J.")
    }

    if (anyNA(covariance) || any(!is.finite(covariance))) {
      stop("covariance matrices must contain only finite values.")
    }

    if (
      max(abs(covariance - t(covariance))) >
        sqrt(.Machine$double.eps)
    ) {
      stop("covariance matrices must be symmetric.")
    }

    tryCatch(
      chol(covariance),
      error = function(e) {
        stop(
          "covariance matrices must be positive definite.",
          call. = FALSE
        )
      }
    )
  }

  covariance_list
}
