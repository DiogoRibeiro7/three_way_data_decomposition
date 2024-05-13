# Define the Higher-order SVD (HOSVD) function for a tensor
# This function decomposes a tensor into a core tensor and a set of orthogonal matrices for each mode.
# It allows for dimensionality reduction via truncation, specified by the 'ranks' parameter.
#
# @param tnsr A Tensor object with K modes.
# @param ranks Optional vector specifying the dimensions of the core tensor for each mode.
#              If NULL, the original dimensions of 'tnsr' are used (no truncation).
# @return A list containing:
#   - \code{Z}: Core tensor with dimensions specified by 'ranks'.
#   - \code{U}: List of orthogonal matrices, one for each mode.
#   - \code{est}: Estimate of 'tnsr' after compression.
#   - \code{fnorm_resid}: Frobenius norm of the residual, measuring the compression error.
# @seealso \code{\link{tucker}}
# @references L. De Lathauwer, B. Moor, J. Vandewalle "A multilinear singular value decomposition".
#             SIAM Journal on Matrix Analysis and Applications, 2000.
# @note The length of 'ranks' must match the number of modes in 'tnsr' if provided.
# @examples
# tnsr <- rand_tensor(c(6, 7, 8))
# hosvd_result <- hosvd(tnsr)
# print(hosvd_result$fnorm_resid)
# hosvd_result_truncated <- hosvd(tnsr, ranks = c(3, 3, 4))
# print(hosvd_result_truncated$fnorm_resid)
hosvd <- function(tnsr, ranks = NULL) {
    requireNamespace("methods", quietly = TRUE)
    stopifnot(is(tnsr, "Tensor"))  # Ensuring input is a Tensor object
    
    num_modes <- tnsr@num_modes
    
    if (is.null(ranks)) {
        ranks <- tnsr@modes
    } else {
        stopifnot(length(ranks) == num_modes)  # Ensure ranks match the number of modes
    }
    
    # Initialize a progress bar
    pb <- txtProgressBar(min = 0, max = num_modes, style = 3)
    U_list <- vector("list", num_modes)
    
    # Loop through each mode and perform SVD on the mode-m matricization of 'tnsr'
    for (m in 1:num_modes) {
        temp_mat <- rs_unfold(tnsr, m = m)@data
        U_list[[m]] <- svd(temp_mat, nu = ranks[m])$u
        setTxtProgressBar(pb, m)
    }
    close(pb)
    
    # Compute the core tensor by applying the transpose of each U matrix to 'tnsr'
    Z <- ttl(tnsr, lapply(U_list, t), ms = 1:num_modes)
    est <- ttl(Z, U_list, ms = 1:num_modes)
    resid <- fnorm(est - tnsr)
    
    # Return a list with results
    list(Z = Z, U = U_list, est = est, fnorm_resid = resid)
}

#' Canonical Polyadic (CP) Decomposition
#'
#' Performs a CP decomposition on a given tensor using an iterative approach
#' to approximate the tensor as a sum of component rank-one tensors.
#'
#' @param tnsr A `Tensor` object on which CP decomposition is performed.
#' @param num_components Integer, the number of components for each mode in the decomposition.
#' @param max_iter Maximum number of iterations to perform (default is 25).
#' @param tol Tolerance for convergence, the algorithm stops if the change in fit 
#'        is less than this value (default is 1e-5).
#'
#' @return A list containing:
#' \itemize{
#'   \item{lambdas}{Vector of lambda values which are the weights of the rank-one tensors.}
#'   \item{U}{List of matrices, one per mode, containing the factors of the decomposition.}
#'   \item{converged}{Logical, indicating whether the algorithm converged within the given iterations.}
#'   \item{estimated}{The estimated tensor reconstructed from the CP decomposition.}
#'   \item{norm_percent}{Percentage of the norm of the original tensor explained by the model.}
#'   \item{fnorm_resid_last}{Final residual Frobenius norm.}
#'   \item{all_resids}{Vector of Frobenius norms of residuals for each iteration.}
#' }
#' @details The CP decomposition is performed using an alternating least squares approach
#' where each step attempts to minimize the Frobenius norm of the difference between the 
#' tensor and its current approximation. The factor matrices are updated iteratively.
#' A progress bar is shown in the console during the computation to indicate progress.
#'
#' @examples
#' matrixA <- matrix(rnorm(20), nrow = 5, ncol = 4)
#' matrixB <- matrix(rnorm(20), nrow = 5, ncol = 4)
#' matrixX1 <- matrix(rnorm(75), nrow = 5, ncol = 15)
#' matrixX2 <- matrix(rnorm(75), nrow = 5, ncol = 15)
#' nComp3 <- 3
#' results <- cp(tnsr = as.tensor(matrixX1), num_components = 2, max_iter = 100, tol = 1e-4)
#' print(results)
#'
#' @export
cp <- function(tnsr, num_components = NULL, max_iter = 25, tol = 1e-5) {
    if (is.null(num_components)) stop("num_components must be specified")
    stopifnot(is(tnsr, "Tensor"))

    num_modes <- tnsr@num_modes
    modes <- tnsr@modes
    U_list <- vector("list", num_modes)
    unfolded_mat <- vector("list", num_modes)
    tnsr_norm <- fnorm(tnsr)

    # Initialize factor matrices using random values
    for (m in 1:num_modes) {
        unfolded_mat[[m]] <- rs_unfold(tnsr, m = m)@data
        U_list[[m]] <- matrix(rnorm(modes[m] * num_components), nrow = modes[m], ncol = num_components)
    }

    est <- tnsr
    curr_iter <- 1
    converged <- FALSE
    fnorm_resid <- numeric(max_iter)

    # Progress bar setup
    pb <- txtProgressBar(min = 0, max = max_iter, style = 3)

    while (curr_iter <= max_iter && !converged) {
        setTxtProgressBar(pb, curr_iter)

        for (m in 1:num_modes) {
            # Updating factor matrices
            V <- Reduce(function(x, y) t(x) %*% x, U_list[-m], init = diag(num_components))
            V_inv <- solve(V)
            tmp <- unfolded_mat[[m]] %*% khatri_rao_list(U_list[-m], reverse = TRUE) %*% V_inv
            lambdas <- apply(tmp, 2, function(x) sqrt(sum(x^2)))
            U_list[[m]] <- sweep(tmp, 2, lambdas, "/")
        }

        # Construct estimated tensor and check for convergence
        est <- ttl(tensor(lambdas, c(num_components, 1, 1)), U_list, ms = 1:num_modes)
        curr_resid <- fnorm(est - tnsr)
        fnorm_resid[curr_iter] <- curr_resid

        if (curr_iter > 1 && abs(curr_resid - fnorm_resid[curr_iter - 1]) / tnsr_norm < tol) {
            converged <- TRUE
        }

        curr_iter <- curr_iter + 1
    }

    close(pb)

    # Results and diagnostics
    fnorm_resid <- fnorm_resid[fnorm_resid != 0]
    norm_percent <- (1 - (tail(fnorm_resid, 1) / tnsr_norm)) * 100

    invisible(list(
        lambdas = lambdas,
        U = U_list,
        converged = converged,
        estimated = est,
        norm_percent = norm_percent,
        fnorm_resid_last = tail(fnorm_resid, 1),
        all_resids = fnorm_resid
    ))
}

# Helper function to compute the Khatri-Rao product of a list of matrices
khatri_rao_list <- function(mat_list, reverse = FALSE) {
    if (reverse) mat_list <- rev(mat_list)
    Reduce(khatri_rao, mat_list)
}

# Helper function to compute the Khatri-Rao product of two matrices
khatri_rao <- function(mat1, mat2) {
    m <- nrow(mat1)
    n <- nrow(mat2)
    kr <- matrix(nrow = m * n, ncol = ncol(mat1))
    for (j in seq_len(ncol(mat1))) {
        kr[, j] <- rep(mat1[, j], each = n) * rep(mat2[, j], times = m)
    }
    kr
}


# Define the Tucker decomposition function for tensors
# This function performs Tucker decomposition, approximating a tensor using a smaller core tensor
# and orthogonal factor matrices for each mode. This method is similar to HOSVD but includes
# the possibility of iterative improvement via ALS.

# @param tnsr A Tensor object with K modes.
# @param ranks A vector specifying the desired modes of the output core Tensor.
# @param max_iter Maximum number of iterations for the ALS algorithm.
# @param tol Tolerance for the relative error in Frobenius norm for convergence.
# @return A list containing:
#   - \code{Z}: Core tensor with dimensions specified by 'ranks'.
#   - \code{U}: List of orthogonal factor matrices, one for each mode.
#   - \code{conv}: Boolean indicating if the algorithm converged based on the tolerance.
#   - \code{est}: Estimate of 'tnsr' after compression.
#   - \code{norm_percent}: Percentage of the Frobenius norm of the original tensor captured by the approximation.
#   - \code{fnorm_resid}: Frobenius norm of the residual error.
#   - \code{all_resids}: Vector containing the Frobenius norm of error for all iterations.
# @seealso \code{\link{hosvd}}, \code{\link{mpca}}
# @references T. Kolda, B. Bader, "Tensor Decompositions and Applications".
#             SIAM Applied Mathematics and Applications, 2009.
# @examples
# tnsr <- rand_tensor(c(6,7,8))
# tucker_result <- tucker(tnsr, ranks = c(3,3,4))
# print(tucker_result$conv)
# print(tucker_result$norm_percent)
# plot(tucker_result$all_resids)
tucker <- function(tnsr, ranks = NULL, max_iter = 25, tol = 1e-5) {
    requireNamespace("methods", quietly = TRUE)
    stopifnot(is(tnsr, "Tensor"))
    if (is.null(ranks)) {
        stop("ranks must be specified")
    }

    # Initialize using a truncated HOSVD
    num_modes <- tnsr@num_modes
    U_list <- vector("list", num_modes)
    for (m in 1:num_modes) {
        temp_mat <- rs_unfold(tnsr, m = m)@data
        U_list[[m]] <- svd(temp_mat, nu = ranks[m])$u
    }

    tnsr_norm <- fnorm(tnsr)
    curr_iter <- 1
    converged <- FALSE
    fnorm_resid <- numeric(max_iter)

    # Set up progress bar
    pb <- txtProgressBar(min = 0, max = max_iter, style = 3)

    while (curr_iter <= max_iter && !converged) {
        setTxtProgressBar(pb, curr_iter)

        for (m in 1:num_modes) {
            # Core tensor Z minus mode m
            X <- ttl(tnsr, lapply(U_list[-m], t), ms = setdiff(1:num_modes, m))

            # Truncated SVD of X
            U_list[[m]] <- svd(rs_unfold(X, m = m)@data, nu = ranks[m])$u
        }

        # Recompute core tensor Z
        Z <- ttm(X, mat = t(U_list[[num_modes]]), m = num_modes)

        # Check convergence
        est <- ttl(Z, U_list, ms = 1:num_modes)
        curr_resid <- fnorm(tnsr - est)
        fnorm_resid[curr_iter] <- curr_resid
        if (curr_iter > 1 && abs(curr_resid - fnorm_resid[curr_iter - 1]) / tnsr_norm < tol) {
            converged <- TRUE
        }

        curr_iter <- curr_iter + 1
    }
    
    close(pb)
    fnorm_resid <- fnorm_resid[fnorm_resid != 0]
    norm_percent <- (1 - (tail(fnorm_resid, 1) / tnsr_norm)) * 100

    # Return results
    list(
        Z = Z,
        U = U_list,
        conv = converged,
        est = est,
        norm_percent = norm_percent,
        fnorm_resid = tail(fnorm_resid, 1),
        all_resids = fnorm_resid
    )
}
