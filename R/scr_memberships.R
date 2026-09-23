#' Convert Fuzzy Memberships to a Hard Partition
#'
#' Convert an observation-by-group membership matrix to a one-hot partition by
#' assigning each observation to the group with the largest membership. Ties
#' are resolved in favour of the first group, matching the MATLAB reference
#' implementation used by Rocci, Vichi, and Ranalli.
#'
#' @param membership Numeric matrix with observations in rows and groups in
#'   columns.
#' @return An integer one-hot matrix with the same dimensions as `membership`.
#' @export
hard_partition <- function(membership) {
  if (!is.matrix(membership) || !is.numeric(membership)) {
    stop("membership must be a numeric matrix.")
  }
  if (nrow(membership) == 0L || ncol(membership) == 0L) {
    stop("membership must have at least one row and one column.")
  }
  if (anyNA(membership) || any(!is.finite(membership))) {
    stop("membership must contain only finite values.")
  }

  group <- max.col(membership, ties.method = "first")
  hard <- matrix(0L, nrow = nrow(membership), ncol = ncol(membership))
  hard[cbind(seq_len(nrow(membership)), group)] <- 1L
  hard
}


#' Hubert-Arabie Adjusted Rand Index from a Contingency Table
#'
#' Compute the adjusted Rand index using the contingency-table formula employed
#' by the SCR MATLAB reference code.
#'
#' @param contingency Non-negative matrix whose rows and columns represent two
#'   hard partitions.
#' @return A scalar adjusted Rand index.
#' @export
adjusted_rand_index <- function(contingency) {
  if (!is.matrix(contingency) || !is.numeric(contingency)) {
    stop("contingency must be a numeric matrix.")
  }
  if (nrow(contingency) == 0L || ncol(contingency) == 0L) {
    stop("contingency must have at least one row and one column.")
  }
  if (
    anyNA(contingency) ||
      any(!is.finite(contingency)) ||
      any(contingency < 0)
  ) {
    stop("contingency must contain finite, non-negative values.")
  }

  n <- sum(contingency)
  if (n < 2) {
    stop("contingency must describe at least two observations.")
  }

  row_pairs <- 0.5 * (sum(rowSums(contingency)^2) - n)
  col_pairs <- 0.5 * (sum(colSums(contingency)^2) - n)
  expected_pairs <- 2 * row_pairs * col_pairs / (n * (n - 1))
  observed_pairs <- 0.5 * (sum(contingency^2) - n)
  denominator <- (row_pairs + col_pairs) / 2 - expected_pairs

  if (abs(denominator) <= .Machine$double.eps) {
    if (identical(as.numeric(contingency), c(n))) {
      return(1)
    }
    stop("adjusted Rand index is undefined for this contingency table.")
  }

  (observed_pairs - expected_pairs) / denominator
}


#' Legacy SCR Name for Hard Partition Conversion
#'
#' Compatibility wrapper for the original MATLAB/R function name.
#'
#' @inheritParams hard_partition
#' @return See `hard_partition()`.
#' @export
ftoh <- function(membership) {
  hard_partition(membership)
}


#' Legacy SCR Name for the Adjusted Rand Index
#'
#' Compatibility wrapper for the original MATLAB/R function name.
#'
#' @param N See `adjusted_rand_index()`.
#' @return See `adjusted_rand_index()`.
#' @export
mrand <- function(N) {
  adjusted_rand_index(N)
}
