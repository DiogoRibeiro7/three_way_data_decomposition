#' Multilinear Principal Component Analysis (MPCA)
#'
#' Performs MPCA on a given tensor by iteratively refining the factor matrices
#' to minimize the reconstruction error, under a rank constraint for each mode.
#'
#' @param tnsr A tensor object, usually an S4 object of class `Tensor` that the MPCA
#'        will be performed on.
#' @param ranks A numeric vector indicating the number of principal components
#'        for each mode of the tensor (excluding the last mode, which is set to 1).
#'        The length of this vector should be one less than the number of modes in the tensor.
#' @param max_iter Maximum number of iterations for the MPCA algorithm (default: 25).
#' @param tol Convergence tolerance for the algorithm, based on the relative change
#'        in Frobenius norm of the tensor approximation between iterations (default: 1e-5).
#'
#' @return A list containing the following components:
#' \itemize{
#'   \item{Z_ext}{The extended core tensor obtained after the final iteration.}
#'   \item{U}{A list of factor matrices, one for each mode, including the last mode.}
#'   \item{converged}{Logical indicating whether the algorithm converged (TRUE) or reached the maximum number of iterations without convergence (FALSE).}
#'   \item{estimated}{The estimated tensor reconstructed from the final factor matrices and core tensor.}
#'   \item{norm_percent}{The percentage of the original tensor's norm that is explained by the estimated tensor.}
#'   \item{fnorm_resid_last}{The last recorded Frobenius norm of the residual tensor (original minus estimated).}
#'   \item{all_resids}{A numeric vector containing all non-zero Frobenius norms of the residuals recorded at each iteration.}
#' }
#'
#' @details The function initializes the factor matrices using a higher-order singular value decomposition
#' (HOSVD) of the unfolded modes of the tensor. It then iteratively updates each factor matrix
#' using a truncated singular value decomposition (SVD) of the mode-wise unfolded tensor,
#' where all but the current mode's factor matrix are held constant. The process continues
#' until the change in the reconstruction error falls below a specified tolerance or the maximum number of iterations is reached.
#'
#' @examples
#' \dontrun{
#' tnsr <- as.tensor(array(rnorm(100), dim = c(5, 5, 4)))
#' results <- mpca(tnsr, ranks = c(2, 2, 1), max_iter = 100, tol = 1e-4)
#' print(results$estimated)  # Display the estimated tensor
#' }
#'
#' @export
mpca <- function(tnsr, ranks = NULL, max_iter = 25, tol = 1e-5) {
    # Validate input parameters
    if (is.null(ranks)) stop("ranks must be specified")
    stopifnot(is(tnsr, "Tensor"), length(ranks) == (tnsr@num_modes - 1))

    # Initialization of variables
    num_modes <- tnsr@num_modes
    modes <- tnsr@modes
    U_list <- vector("list", num_modes)
    unfolded_mat <- vector("list", num_modes)

    # Initialize factor matrices using HOSVD for the first M-1 modes
    for (m in 1:(num_modes - 1)) {
        unfolded_mat[[m]] <- rs_unfold(tnsr, m = m)@data
        mode_m_cov <- unfolded_mat[[m]] %*% t(unfolded_mat[[m]])
        U_list[[m]] <- svd(mode_m_cov, nu = ranks[m])$u
    }

    # Initialize extended core tensor using transposes of the factor matrices
    Z_ext <- ttl(tnsr, lapply(U_list[-num_modes], t), ms = 1:(num_modes - 1))
    tnsr_norm <- fnorm(tnsr)
    curr_iter <- 1
    converged <- FALSE
    fnorm_resid <- numeric(max_iter)

    # Setup the progress bar for visual feedback
    pb <- txtProgressBar(min = 0, max = max_iter, style = 3)

    # Main iterative process
    while (curr_iter <= max_iter && !converged) {
        setTxtProgressBar(pb, curr_iter)

        for (m in 1:(num_modes - 1)) {
            # Update factor matrices by excluding the current mode
            X <- ttl(tnsr, lapply(U_list[-c(m, num_modes)], t), ms = setdiff(1:(num_modes - 1), m))
            U_list[[m]] <- svd(rs_unfold(X, m = m)@data, nu = ranks[m])$u
        }

        # Recompute the core tensor
        Z_ext <- ttm(X, mat = t(U_list[[num_modes - 1]]), m = num_modes - 1)

        # Check for convergence
        converged <- check_convergence(tnsr, Z_ext, U_list, curr_iter, fnorm_resid, tnsr_norm, tol)
        curr_iter <- curr_iter + 1
    }

    # Clean up the progress bar
    close(pb)

    # Prepare and return the results
    est <- ttl(Z_ext, U_list[-num_modes], ms = 1:(num_modes - 1))
    fnorm_resid <- fnorm_resid[fnorm_resid != 0]
    norm_percent <- (1 - (tail(fnorm_resid, 1) / tnsr_norm)) * 100

    invisible(list(
        Z_ext = Z_ext, 
        U = U_list, 
        converged = converged, 
        estimated = est, 
        norm_percent = norm_percent, 
        fnorm_resid_last = tail(fnorm_resid, 1), 
        all_resids = fnorm_resid
    ))
}

# Function to check for convergence in MPCA
check_convergence <- function(tnsr, Z_ext, U_list, curr_iter, fnorm_resid, tnsr_norm, tol) {
    est <- ttl(Z_ext, U_list[-tnsr@num_modes], ms = 1:(tnsr@num_modes - 1))
    curr_resid <- fnorm(tnsr - est)
    fnorm_resid[curr_iter] <<- curr_resid
    if (curr_iter == 1) return(FALSE)
    return(abs(curr_resid - fnorm_resid[curr_iter - 1]) / tnsr_norm < tol)
}
