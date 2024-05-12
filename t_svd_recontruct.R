#' Reconstruct Tensor From TSVD
#'
#' Reconstructs the original 3-Tensor after it has been decomposed into U, S, and V
#' using `t_svd`. This function performs the tensor multiplication of U, S, and the
#' transpose of V to obtain the reconstructed tensor.
#'
#' @param L A list containing the components U, S, and V as returned by `t_svd`.
#' @return A 3-Tensor reconstructed from the TSVD components.
#' @seealso \link{t_svd}
#' @examples
#' tnsr <- rand_tensor(c(10, 10, 10))  # Ensure this function generates a proper 3-tensor
#' tsvdD <- t_svd(tnsr)
#' reconstructed <- t_svd_reconstruct(tsvdD)
#' # Assess reconstruction quality
#' reconstruction_error <- 1 - fnorm(reconstructed - tnsr) / fnorm(tnsr)
#' print(reconstruction_error)
#' @export
t_svd_reconstruct <- function(L) {
    if (!is.list(L) || !all(c("U", "S", "V") %in% names(L))) {
        stop("Input must be a list containing U, S, and V components from TSVD.")
    }
    
    # Perform tensor multiplication: U * S * transpose(V)
    # Assuming t_mult and t (transpose) are available and handle tensor operations
    t_mult(t_mult(L$U, L$S), t(L$V))
}
