#' Calculate the Modified Rand Index
#'
#' Computes the modified Rand index (Hubert & Arabie, 1985), a measure of similarity
#' between two data clusterings, based on a contingency table (cross-tabulation) of the
#' cluster assignments.
#'
#' @param N A matrix representing the contingency table where the element at the i-th row
#'   and j-th column represents the number of objects that are in both cluster i of the
#'   first clustering and cluster j of the second clustering.
#'
#' @return A numeric value representing the modified Rand index, which measures the
#'   agreement of the two clusterings adjusted for chance. It ranges from -1 (no agreement)
#'   to 1 (perfect agreement).
#'
#' @examples
#' # Assuming two clustering results are given as factors
#' cluster1 <- factor(c(1, 1, 2, 2, 3, 3))
#' cluster2 <- factor(c(1, 2, 2, 1, 3, 3))
#' # Create contingency table
#' N <- table(cluster1, cluster2)
#' # Calculate modified Rand index
#' mri_value <- mrand(N)
#' print(mri_value)
#'
#' @export
mrand <- function(N) {
  n <- sum(N)
  sumi <- 0.5 * (sum(rowSums(N)^2) - n)
  sumj <- 0.5 * (sum(colSums(N)^2) - n)
  pb <- 2 * sumi * sumj / (n * (n - 1))
  mri <- (0.5 * (sum(N^2) - n) - pb) / ((sumi + sumj) / 2 - pb)
  
  return(mri)
}
