#' Nearest Kronecker Covariance Approximation
#'
#' Approximate a common covariance matrix by a separable covariance
#' `Sigma_O %x% Sigma_V` using the Pitsianis--Van Loan rearrangement and a
#' leading singular-vector approximation.
#'
#' The returned covariance factors are symmetrized, projected to positive
#' definiteness when necessary, and scale-normalized so the variable covariance
#' has mean diagonal one. The occasion factor absorbs the reciprocal scale.
#'
#' @param covariance Positive-definite covariance matrix of dimension
#'   `(J*K) x (J*K)`.
#' @param variables Number of variables J.
#' @param occasions Number of occasions K.
#' @param eigen_floor Relative eigenvalue floor used only if an approximate
#'   factor requires positive-definite projection.
#' @return A list containing normalized covariance factors, the reconstructed
#'   Kronecker covariance, relative approximation errors, rank-one explained
#'   fraction, rearrangement singular values, and projection diagnostics.
#' @examples
#' variable_covariance <- matrix(c(2, 0.3, 0.3, 1), 2, 2)
#' occasion_covariance <- matrix(c(1.5, 0.2, 0.2, 0.8), 2, 2)
#' covariance <- kronecker(occasion_covariance, variable_covariance)
#' approximation <- scr_nearest_kronecker_covariance(
#'   covariance,
#'   variables = 2,
#'   occasions = 2
#' )
#' approximation$relative_error
#' @export
scr_nearest_kronecker_covariance <- function(
  covariance,
  variables,
  occasions,
  eigen_floor = 1e-10
) {
  validate_covariance_dimensions(
    covariance = covariance,
    variables = variables,
    occasions = occasions
  )
  validate_eigen_floor(eigen_floor)

  rearranged <- rearrange_kronecker_covariance(
    covariance = covariance,
    variables = variables,
    occasions = occasions
  )
  decomposition <- svd(rearranged)

  leading <- decomposition$d[1L]
  variable_vector <- sqrt(leading) * decomposition$u[, 1L]
  occasion_vector <- sqrt(leading) * decomposition$v[, 1L]

  variable_raw <- matrix(
    variable_vector,
    nrow = variables,
    ncol = variables
  )
  occasion_raw <- matrix(
    occasion_vector,
    nrow = occasions,
    ncol = occasions
  )

  variable_raw <- (variable_raw + t(variable_raw)) / 2
  occasion_raw <- (occasion_raw + t(occasion_raw)) / 2

  if (mean(diag(variable_raw)) < 0) {
    variable_raw <- -variable_raw
    occasion_raw <- -occasion_raw
  }

  raw_kronecker <- kronecker(occasion_raw, variable_raw)
  raw_error <- relative_frobenius_error(
    covariance,
    raw_kronecker
  )

  variable_projection <- project_positive_definite(
    variable_raw,
    eigen_floor = eigen_floor
  )
  occasion_projection <- project_positive_definite(
    occasion_raw,
    eigen_floor = eigen_floor
  )

  variable_covariance <- variable_projection$matrix
  occasion_covariance <- occasion_projection$matrix

  scale <- mean(diag(variable_covariance))
  if (!is.finite(scale) || scale <= 0) {
    stop("failed to normalize the variable covariance scale.")
  }

  variable_covariance <- variable_covariance / scale
  occasion_covariance <- occasion_covariance * scale

  kronecker_covariance <- kronecker(
    occasion_covariance,
    variable_covariance
  )
  projected_error <- relative_frobenius_error(
    covariance,
    kronecker_covariance
  )

  singular_energy <- sum(decomposition$d^2)
  rank_one_fraction <- if (singular_energy > 0) {
    leading^2 / singular_energy
  } else {
    1
  }

  list(
    variable_covariance = variable_covariance,
    occasion_covariance = occasion_covariance,
    covariance = kronecker_covariance,
    relative_error = projected_error,
    raw_relative_error = raw_error,
    rank_one_fraction = as.numeric(rank_one_fraction),
    singular_values = decomposition$d,
    projected = c(
      variable = variable_projection$projected,
      occasion = occasion_projection$projected
    ),
    projection_shift = c(
      variable = variable_projection$shift,
      occasion = occasion_projection$shift
    ),
    dimensions = c(
      variables = as.integer(variables),
      occasions = as.integer(occasions)
    )
  )
}


#' Construct a Kronecker-plus-Nugget Covariance
#'
#' Form the controlled covariance-departure family
#'
#' `Sigma(tau) = Sigma_O %x% Sigma_V + tau I`
#'
#' with `tau >= 0`. Since the Kronecker factor eigenvalues are products of
#' factor eigenvalues, the nugget model admits direct spectral diagnostics.
#'
#' @param variable_covariance Positive-definite J-by-J covariance.
#' @param occasion_covariance Positive-definite K-by-K covariance.
#' @param nugget Non-negative scalar tau.
#' @param details Logical; when TRUE return covariance and spectral details.
#' @return The covariance matrix, or a list of covariance and diagnostics when
#'   `details = TRUE`.
#' @examples
#' covariance <- scr_kronecker_nugget_covariance(
#'   variable_covariance = diag(2),
#'   occasion_covariance = diag(2),
#'   nugget = 0.25,
#'   details = TRUE
#' )
#' covariance$condition_number
#' @export
scr_kronecker_nugget_covariance <- function(
  variable_covariance,
  occasion_covariance,
  nugget = 0,
  details = FALSE
) {
  validate_positive_definite_matrix(
    variable_covariance,
    "variable covariance"
  )
  validate_positive_definite_matrix(
    occasion_covariance,
    "occasion covariance"
  )

  if (
    !is.numeric(nugget) ||
      length(nugget) != 1L ||
      is.na(nugget) ||
      !is.finite(nugget) ||
      nugget < 0
  ) {
    stop("nugget must be a non-negative finite number.")
  }

  if (!is.logical(details) || length(details) != 1L || is.na(details)) {
    stop("details must be TRUE or FALSE.")
  }

  base <- kronecker(
    occasion_covariance,
    variable_covariance
  )
  dimension <- nrow(base)
  covariance <- base + nugget * diag(dimension)

  if (!details) {
    return(covariance)
  }

  variable_eigenvalues <- eigen(
    variable_covariance,
    symmetric = TRUE,
    only.values = TRUE
  )$values
  occasion_eigenvalues <- eigen(
    occasion_covariance,
    symmetric = TRUE,
    only.values = TRUE
  )$values

  eigenvalues <- as.vector(
    outer(
      variable_eigenvalues,
      occasion_eigenvalues,
      FUN = "*"
    )
  ) + nugget

  list(
    covariance = covariance,
    separable_covariance = base,
    nugget = nugget,
    eigenvalues = eigenvalues,
    log_determinant = sum(log(eigenvalues)),
    condition_number = max(eigenvalues) / min(eigenvalues)
  )
}


#' Covariance Parameter Count for SCR Extensions
#'
#' Count free parameters in alternative common-covariance structures.
#'
#' @param variables Number of variables J.
#' @param occasions Number of occasions K.
#' @param model One of `"separable"`, `"nugget"`, or `"unrestricted"`.
#' @return Integer number of free covariance parameters.
#' @export
scr_covariance_parameter_count <- function(
  variables,
  occasions,
  model = c("separable", "nugget", "unrestricted")
) {
  validate_model_selection_dimension(variables, "variables")
  validate_model_selection_dimension(occasions, "occasions")
  model <- match.arg(model)

  separable <- variables * (variables + 1L) / 2L +
    occasions * (occasions + 1L) / 2L -
    1L

  count <- switch(
    model,
    separable = separable,
    nugget = separable + 1L,
    unrestricted = {
      dimension <- variables * occasions
      dimension * (dimension + 1L) / 2L
    }
  )

  as.integer(count)
}


rearrange_kronecker_covariance <- function(
  covariance,
  variables,
  occasions
) {
  rearranged <- matrix(
    0,
    nrow = variables * variables,
    ncol = occasions * occasions
  )

  column <- 1L
  for (occasion_col in seq_len(occasions)) {
    columns <- seq.int(
      from = (occasion_col - 1L) * variables + 1L,
      length.out = variables
    )

    for (occasion_row in seq_len(occasions)) {
      rows <- seq.int(
        from = (occasion_row - 1L) * variables + 1L,
        length.out = variables
      )
      block <- covariance[rows, columns, drop = FALSE]
      rearranged[, column] <- as.vector(block)
      column <- column + 1L
    }
  }

  rearranged
}


relative_frobenius_error <- function(target, approximation) {
  denominator <- sqrt(sum(target^2))
  if (!is.finite(denominator) || denominator <= 0) {
    stop("target matrix must have positive Frobenius norm.")
  }

  sqrt(sum((target - approximation)^2)) / denominator
}


project_positive_definite <- function(matrix, eigen_floor) {
  symmetric <- (matrix + t(matrix)) / 2
  decomposition <- eigen(symmetric, symmetric = TRUE)
  values <- decomposition$values

  scale <- max(abs(values), 1)
  floor_value <- eigen_floor * scale
  adjusted <- pmax(values, floor_value)
  projected <- any(values < floor_value)

  result <- decomposition$vectors %*%
    diag(adjusted, nrow = length(adjusted)) %*%
    t(decomposition$vectors)
  result <- (result + t(result)) / 2

  list(
    matrix = result,
    projected = projected,
    shift = max(adjusted - values)
  )
}


validate_covariance_dimensions <- function(
  covariance,
  variables,
  occasions
) {
  validate_model_selection_dimension(variables, "variables")
  validate_model_selection_dimension(occasions, "occasions")
  validate_positive_definite_matrix(covariance, "covariance")

  expected <- variables * occasions
  if (
    length(dim(covariance)) != 2L ||
      any(dim(covariance) != c(expected, expected))
  ) {
    stop(
      "covariance dimensions must equal (variables * occasions) squared."
    )
  }

  invisible(TRUE)
}


validate_eigen_floor <- function(eigen_floor) {
  if (
    !is.numeric(eigen_floor) ||
      length(eigen_floor) != 1L ||
      is.na(eigen_floor) ||
      !is.finite(eigen_floor) ||
      eigen_floor <= 0
  ) {
    stop("eigen_floor must be a positive finite number.")
  }

  invisible(TRUE)
}
