#' Perform T-SVD (Tensor Singular Value Decomposition) on a 3D tensor
#'
#' This function performs T-SVD on a 3D tensor. T-SVD is a decomposition technique
#' that decomposes a tensor into three components: U, V, and S. U and V are matrices
#' that represent the spatial modes of the tensor, and S is a tensor that represents
#' the temporal mode of the tensor.
#'
#' @param tnsr The input 3D tensor
#'
#' @return A list containing the U, V, and S components of the T-SVD decomposition
#'
#' @details This function implements the T-SVD algorithm for 3D tensors. It first checks
#' if the input tensor has exactly 3 modes. If not, an error is thrown. Then, it extracts
#' the dimensions of the tensor and initializes a progress bar. Next, it defines an inverse
#' FFT function and applies FFT along the third mode of the tensor, permuting the dimensions
#' to arrange the third mode along the third dimension. It then initializes arrays to store
#' the SVD components and computes the SVD for each face of the permuted tensor. Finally,
#' it applies inverse FFT to transform the frequency-domain results back to the spatial domain.
#' The resulting U, V, and S components are returned as a list.
#'
#' @examples
#' # Create a 3D tensor
#' tnsr <- as.tensor(array(1:27, dim = c(3, 3, 3)))
#'
#' # Perform T-SVD on the tensor
#' result <- t_svd(tnsr)
#' 
#' # Access the U, V, and S components
#' U <- result$U
#' V <- result$V
#' S <- result$S
#'
#' @export
t_svd <- function(tnsr) {
    if (tnsr@num_modes != 3) stop("T-SVD only implemented for 3d so far")
    
    # Extract tensor dimensions
    n1 <- tnsr@modes[1]
    n2 <- tnsr@modes[2]
    n3 <- tnsr@modes[3]
    
    # Initialize progress bar
    pb <- txtProgressBar(min = 0, max = n3, style = 3)
    
    # Define inverse FFT
    ifft <- function(x) suppressWarnings(as.numeric(fft(x, inverse = TRUE)) / length(x))
    
    # Apply FFT along mode 3, then permute to arrange the third mode along the third dimension
    fftz <- aperm(apply(tnsr@data, MARGIN = 1:2, fft), c(2, 3, 1))
    
    # Arrays to store the SVD components
    U_arr <- array(0, dim = c(n1, n1, n3))
    V_arr <- array(0, dim = c(n2, n2, n3))
    S_arr <- array(0, dim = c(n1, n2, n3))
    
    # Compute SVD for each face of the permuted tensor
    for (j in 1:n3) {
        setTxtProgressBar(pb, j)
        decomp <- svd(fftz[,,j], nu = n1, nv = n2)
        U_arr[,,j] <- decomp$u
        V_arr[,,j] <- decomp$v
        S_arr[,,j] <- diag(decomp$d, nrow = n1, ncol = n2)
    }
    
    close(pb)
    
    # Apply inverse FFT to transform frequency-domain results back to the spatial domain
    U <- as.tensor(aperm(apply(U_arr, MARGIN = 1:2, ifft), c(2, 3, 1)))
    V <- as.tensor(aperm(apply(V_arr, MARGIN = 1:2, ifft), c(2, 3, 1)))
    S <- as.tensor(aperm(apply(S_arr, MARGIN = 1:2, ifft), c(2, 3, 1)))
    
    invisible(list(U = U, V = V, S = S))
}
