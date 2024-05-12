#' Define an S4 class 'Tensor' to handle tensor data
#'
#' @slot num_modes Number of dimensions in the tensor
#' @slot modes Vector of dimensions sizes for each mode
#' @slot data Actual tensor data as a multi-dimensional array
#' 
#' @return An instance of the 'Tensor' class
#' 
#' @export
setClass(
    "Tensor",
    representation(
        num_modes = "integer", # Number of dimensions in the tensor
        modes = "integer", # Vector of dimensions sizes for each mode
        data = "array" # Actual tensor data as a multi-dimensional array
    ),
    validity = function(object) {
        # Initialize a vector to store error messages
        errors <- character()

        # Check if 'num_modes' is a single positive integer
        if (length(object@num_modes) != 1 || object@num_modes <= 0) {
            errors <- c(errors, "'num_modes' must be a single, positive integer.")
        }

        # Check if 'modes' length matches 'num_modes'
        if (length(object@modes) != object@num_modes) {
            errors <- c(errors, paste(
                "Length of 'modes' must be equal to 'num_modes'. Expected",
                object@num_modes, "got", length(object@modes)
            ))
        }

        # Check if all entries in 'modes' are strictly positive integers
        if (any(object@modes <= 0)) {
            errors <- c(errors, "'modes' must contain strictly positive values; if any mode is 1, consider a smaller 'num_modes'.")
        }

        # Validate the dimensionality of 'data' matches 'modes'
        if (!identical(dim(object@data), object@modes)) {
            errors <- c(errors, "Dimensions of 'data' must exactly match the 'modes'.")
        }

        # Return TRUE if there are no errors, or return the vector of error messages
        if (length(errors) == 0) TRUE else errors
    }
)


#' Unfolds a tensor along specified dimensions
#'
#' @param tnsr The tensor to be unfolded
#' @param row_idx The indices of the rows to be unfolded
#' @param col_idx The indices of the columns to be unfolded
#'
#' @return The unfolded tensor
#'
#' @export
setGeneric(
    name = "unfold",
    def = function(tnsr, row_idx, col_idx) {
        standardGeneric("unfold")
    })


#' Unfold a Tensor Along a Specified Mode
#'
#' This function unfolds a tensor along a specified mode, turning it into a matrix. The function is generic and needs to have specific methods implemented for different tensor classes.
#'
#' @param tnsr A tensor object.
#' @param m An integer specifying the mode along which the tensor should be unfolded.
#'
#' @return A matrix obtained by unfolding the tensor along the specified mode.
#'
#' @examples
#' \dontrun{
#' tnsr <- rand_tensor(c(10, 10, 10)) # Assuming rand_tensor generates a 3D tensor
#' unfolded_matrix <- k_unfold(tnsr, 1)
#' }
#' @export
setGeneric(
    name = "k_unfold",
    def = function(tnsr, m) {
        standardGeneric("k_unfold")
    })

#' Unfold a Tensor Along a Specified Mode
#'
#' This generic function unfolds a tensor along a specified mode `m`, 
#' transforming it into a matrix. The unfolding process rearranges the tensor's elements 
#' into a matrix such that the mode `m` dimensions form the rows, and the remaining dimensions 
#' form the columns of the resulting matrix.
#'
#' @param tnsr A tensor object that needs to be unfolded.
#' @param m An integer representing the mode along which the tensor will be unfolded.
#'
#' @return A matrix obtained by unfolding the tensor along the specified mode.
#'
#' @examples
#' \dontrun{
#'   tnsr <- array(1:24, dim = c(4, 3, 2))  # Create a 3D tensor
#'   unfolded_matrix <- rs_unfold(tnsr, 2)
#'   print(unfolded_matrix)
#' }
#' @export
#' @name rs_unfold
#' @aliases rs_unfold
#' @rdname rs_unfold
setGeneric("rs_unfold",
           def = function(tnsr, m) {
               standardGeneric("rs_unfold")
           })

#' Unfold a Tensor Along a Specified Mode
#'
#' This generic function unfolds a tensor along a specified mode `m`, 
#' transforming it into a matrix. The unfolding process rearranges the tensor's elements 
#' into a matrix such that the mode `m` dimensions form the rows, and the remaining dimensions 
#' form the columns of the resulting matrix.
#'
#' @param tnsr A tensor object that needs to be unfolded.
#' @param m An integer representing the mode along which the tensor will be unfolded.
#'
#' @return A matrix obtained by unfolding the tensor along the specified mode.
#'
#' @examples
#' \dontrun{
#'   tnsr <- array(1:24, dim = c(4, 3, 2))  # Create a 3D tensor
#'   unfolded_matrix <- rs_unfold(tnsr, 2)
#'   print(unfolded_matrix)
#' }
#' @export
#' @name rs_unfold
#' @aliases rs_unfold
#' @rdname rs_unfold
setGeneric("rs_unfold",
           def = function(tnsr, m) {
               standardGeneric("rs_unfold")
           })


#' Unfold a Tensor Along a Specified Mode (cs_unfold)
#'
#' This generic function unfolds a tensor along a specified mode `m`, which reorganizes
#' the elements of the tensor into a two-dimensional matrix. The unfolding is done such
#' that one dimension of the matrix corresponds to the `m`-th mode, and the other dimension
#' corresponds to all other modes concatenated.
#'
#' @param tnsr A tensor object that needs to be unfolded. The tensor object should be 
#'        compatible with the methods that will be implemented for this generic function.
#' @param m An integer representing the mode along which the tensor will be unfolded.
#'
#' @return A matrix representing the unfolded tensor along the specified mode. 
#'         The structure of the matrix depends on the specific method implementation
#'         for the class of the input tensor.
#'
#' @examples
#' \dontrun{
#'   # Assuming 'ArrayTensor' is a class of tensor and is properly defined
#'   tnsr <- ArrayTensor(data = array(1:24, dim = c(4, 3, 2)))
#'   unfolded_matrix <- cs_unfold(tnsr, 2)
#'   print(unfolded_matrix)
#' }
#'
#' @export
#' @name cs_unfold
#' @aliases cs_unfold
#' @rdname cs_unfold
setGeneric("cs_unfold",
           def = function(tnsr, m) {
               standardGeneric("cs_unfold")
           })


#' Sum Over a Specified Mode of a Tensor
#'
#' This function computes the sum of elements along a specified mode `m` of a tensor `tnsr`.
#' If `drop` is TRUE, the result will drop the specified mode, reducing the dimension of the tensor.
#'
#' @param tnsr The tensor object over which to compute the sum.
#' @param m An integer, the mode along which the sum is computed.
#' @param drop A logical value indicating whether to drop the summed mode in the result.
#'
#' @return The result of summing the tensor along the specified mode, potentially with reduced
#'         dimension if `drop` is TRUE.
#'
#' @examples
#' \dontrun{
#'   tnsr <- array(1:24, dim = c(4, 3, 2))  # Create a 3D tensor
#'   summed_tensor <- modeSum(tnsr, 1, TRUE)
#'   print(summed_tensor)  # Display the result
#' }
#'
#' @export
#' @name modeSum
#' @aliases modeSum
#' @rdname modeSum
setGeneric("modeSum",
           def = function(tnsr, m, drop) {
               standardGeneric("modeSum")
           })


# Define a generic function 'fnorm'
# This function is intended to compute the Frobenius norm of a tensor.
# @param tnsr A tensor (multi-dimensional array) whose norm is to be calculated.
# @return Returns the Frobenius norm of the tensor.
setGeneric(
  name = "fnorm",
  def = function(tnsr) {
    # Check if 'tnsr' is a matrix or an array (which includes tensors)
    if (!is.matrix(tnsr) && !is.array(tnsr)) {
      stop("Input 'tnsr' must be a matrix or an array.")
    }
    standardGeneric("fnorm")
  }
)


# Define a generic function 'innerProd'
# This function computes the inner product between two tensors.
# @param tnsr1 A tensor (multi-dimensional array), the first tensor.
# @param tnsr2 A tensor (multi-dimensional array), the second tensor.
# @return Returns the inner product of the two tensors.
setGeneric(
  name = "innerProd",
  def = function(tnsr1, tnsr2) {
    # Check if both 'tnsr1' and 'tnsr2' are matrices or arrays
    if (!is.matrix(tnsr1) && !is.array(tnsr1)) {
      stop("Input 'tnsr1' must be a matrix or an array.")
    }
    if (!is.matrix(tnsr2) && !is.array(tnsr2)) {
      stop("Input 'tnsr2' must be a matrix or an array.")
    }
    standardGeneric("innerProd")
  }
)
