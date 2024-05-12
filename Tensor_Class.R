# Define an S4 class 'Tensor' to handle tensor data
setClass(
  "Tensor",
  representation(
    num_modes = "integer",  # Number of dimensions in the tensor
    modes = "integer",      # Vector of dimensions sizes for each mode
    data = "array"          # Actual tensor data as a multi-dimensional array
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
      errors <- c(errors, paste("Length of 'modes' must be equal to 'num_modes'. Expected", 
                                object@num_modes, "got", length(object@modes)))
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


# Define a generic function 'unfold'
setGeneric(
  name = "unfold",
  def = function(tnsr, row_idx, col_idx) {
    standardGeneric("unfold")
  },
  # Add documentation within the generic to guide developers on how to implement specific methods
  doc = "Generic function for unfolding a tensor into a matrix based on specified row and column indices. 
         This function needs specific methods to be defined for different tensor classes."
)


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
#'   tnsr <- rand_tensor(c(10, 10, 10))  # Assuming rand_tensor generates a 3D tensor
#'   unfolded_matrix <- k_unfold(tnsr, 1)
#' }
#' @export
setGeneric("k_unfold", function(tnsr, m) standardGeneric("k_unfold"))

setGeneric(
  name = "k_unfold",
  def = function(tnsr, m) {
    standardGeneric("k_unfold")
  },
)

