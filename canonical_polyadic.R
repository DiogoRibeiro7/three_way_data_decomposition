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
