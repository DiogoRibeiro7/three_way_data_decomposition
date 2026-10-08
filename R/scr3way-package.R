#' scr3way: Simultaneous Clustering and Reduction for Three-Way Data
#'
#' Implements simultaneous clustering and dimensionality reduction methods for
#' three-way data, including the Gaussian SCR baseline and research extensions
#' for Tucker3 centroid reduction, structural model selection, stability,
#' covariance departures, sparse discriminating subspaces, and robust
#' Student-t likelihoods.
#'
#' The SCR baseline follows Rocci, Vichi, and Ranalli (2025).
#'
#' @references
#' Rocci, R., Vichi, M., & Ranalli, M. (2025).
#' Mixture models for simultaneous classification and reduction of three-way
#' data. *Computational Statistics*, 40, 469--507.
#' \doi{10.1007/s00180-024-01478-1}
#'
#' @import rTensor
#' @importFrom methods is
#' @keywords internal
"_PACKAGE"
