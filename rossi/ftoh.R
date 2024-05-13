#' Convert a Fuzzy Partition to a Hard Partition
#'
#' This function takes a fuzzy partition matrix where each element represents
#' the degree of belonging of a data point to a cluster and converts it into
#' a hard partition matrix where each data point is assigned to exactly one
#' cluster with a membership value of 1 (the cluster with the highest membership
#' value in the fuzzy partition), and 0 for all other clusters.
#'
#' @param Uf A numeric matrix where rows represent data points and columns represent
#' clusters. Each entry in the matrix is a fuzzy membership value between 0 and 1.
#'
#' @return A binary matrix of the same dimensions as Uf where each row has exactly one
#' '1' indicating the cluster to which the data point most strongly belongs, and '0's
#' elsewhere.
#'
#' @examples
#' Uf <- matrix(c(0.1, 0.6, 0.3, 0.7, 0.2, 0.1), nrow = 2, byrow = TRUE)
#' ftoh(Uf)
#' # Returns:
#' #      [,1] [,2] [,3]
#' # [1,]    0    1    0
#' # [2,]    1    0    0
#'
#' @export
ftoh <- function(Uf) {
  # Get the number of rows (data points) and columns (clusters)
  n <- nrow(Uf)
  nk <- ncol(Uf)
  
  # Initialize the hard partition matrix with zeros
  Uh <- matrix(0, n, nk)
  
  # Find the indices of the maximum values in each row
  ind <- max.col(Uf, ties.method = "first")
  
  # Set the corresponding entries to 1
  Uh[cbind(1:n, ind)] <- 1
  
  # Return the hard partition matrix
  return(Uh)
}
