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
    }
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
#' tnsr <- rand_tensor(c(10, 10, 10)) # Assuming rand_tensor generates a 3D tensor
#' unfolded_matrix <- k_unfold(tnsr, 1)
#' }
#' @export
setGeneric(
    name = "k_unfold",
    def = function(tnsr, m) {
        standardGeneric("k_unfold")
    }
)

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
#' tnsr <- array(1:24, dim = c(4, 3, 2)) # Create a 3D tensor
#' unfolded_matrix <- rs_unfold(tnsr, 2)
#' print(unfolded_matrix)
#' }
#' @export
#' @name rs_unfold
#' @aliases rs_unfold
#' @rdname rs_unfold
setGeneric("rs_unfold",
    def = function(tnsr, m) {
        standardGeneric("rs_unfold")
    }
)

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
#' tnsr <- array(1:24, dim = c(4, 3, 2)) # Create a 3D tensor
#' unfolded_matrix <- rs_unfold(tnsr, 2)
#' print(unfolded_matrix)
#' }
#' @export
#' @name rs_unfold
#' @aliases rs_unfold
#' @rdname rs_unfold
setGeneric("rs_unfold",
    def = function(tnsr, m) {
        standardGeneric("rs_unfold")
    }
)


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
#' # Assuming 'ArrayTensor' is a class of tensor and is properly defined
#' tnsr <- ArrayTensor(data = array(1:24, dim = c(4, 3, 2)))
#' unfolded_matrix <- cs_unfold(tnsr, 2)
#' print(unfolded_matrix)
#' }
#'
#' @export
#' @name cs_unfold
#' @aliases cs_unfold
#' @rdname cs_unfold
setGeneric("cs_unfold",
    def = function(tnsr, m) {
        standardGeneric("cs_unfold")
    }
)


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
#' tnsr <- array(1:24, dim = c(4, 3, 2)) # Create a 3D tensor
#' summed_tensor <- modeSum(tnsr, 1, TRUE)
#' print(summed_tensor) # Display the result
#' }
#'
#' @export
#' @name modeSum
#' @aliases modeSum
#' @rdname modeSum
setGeneric("modeSum",
    def = function(tnsr, m, drop) {
        standardGeneric("modeSum")
    }
)


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


# Define the method 'initialize' for objects of class 'Tensor'
# This method initializes a Tensor object with the specified number of modes, mode dimensions, and data.
# @param .Object An instance of class 'Tensor' to be initialized.
# @param num_modes Integer, the number of modes (dimensions) of the tensor. If NULL, it will be inferred from 'data'.
# @param modes Integer vector, specifying the size of each mode (dimension). If NULL, it will be inferred from 'data'.
# @param data An array or vector representing the tensor's data. If data is a vector and no mode information is provided, it is treated as a single-mode tensor.
# @return An initialized 'Tensor' object.
setMethod(
    f = "initialize",
    signature = "Tensor",
    definition = function(.Object, num_modes = NULL, modes = NULL, data = NULL) {
        # Infer num_modes from data if not provided
        if (is.null(num_modes)) {
            if (is.vector(data)) {
                num_modes <- 1L
            } else {
                num_modes <- length(dim(data))
            }
        }

        # Infer modes from data if not provided
        if (is.null(modes)) {
            if (is.vector(data)) {
                modes <- length(data)
            } else {
                modes <- dim(data)
            }
        }

        # Assign the properties to the Tensor object
        .Object@num_modes <- num_modes
        .Object@modes <- modes
        .Object@data <- array(data, dim = modes)

        # Validate the Tensor object to ensure it meets all class requirements
        validObject(.Object)

        # Return the initialized Tensor object
        return(.Object)
    }
)

# Set options to suppress warnings globally in R
# This setting is applied to prevent warnings from being displayed.
# Use this setting cautiously, as it can hide important alerts about potential issues in your code.
options(warn = -1)


# Define the 'dim' method for objects of class 'Tensor'
# This method returns the dimensions of the tensor, specifically the size of each mode.
# @param x An object of class 'Tensor'.
# @return Returns the size of each mode of the tensor as a numeric vector.
setMethod(
    f = "dim",
    signature = "Tensor",
    definition = function(x) {
        # Ensure 'x' is a valid Tensor object
        validObject(x)
        # Access and return the 'modes' slot of the Tensor object
        return(x@modes)
    }
)

# Define the 'show' method for objects of class 'Tensor'
# This method provides a user-friendly display of the Tensor object's properties.
# @param object An instance of class 'Tensor'.
setMethod(
    f = "show",
    signature = "Tensor",
    definition = function(object) {
        # Ensure 'object' is a valid Tensor object
        validObject(object)

        # Display the basic information about the Tensor
        cat("Numeric Tensor of", object@num_modes, "Modes\n", sep = " ")

        # Display the dimensions (modes) of the Tensor
        cat("Modes: ", paste(object@modes, collapse = ", "), "\n", sep = " ")

        # Display a preview of the tensor's data
        cat("Data: \n")
        # Print only the first few elements of the tensor's data to avoid overwhelming the user
        print(head(object@data))
    }
)


# Define the 'print' method for objects of class 'Tensor'
# This method is intended to print the Tensor object in a user-readable format
# by delegating to the 'show' method.
# @param x An instance of class 'Tensor'.
# @param ... Additional arguments that can be handled by the 'print' method.
setMethod(
    f = "print",
    signature = "Tensor",
    definition = function(x, ...) {
        # Call the 'show' method for the 'Tensor' class
        # This displays the Tensor object with its structure and data
        show(x)
        # It is common to invisibly return the object itself after printing
        invisible(x)
    }
)


# Define the 'head' method for objects of class 'Tensor'
# This method displays the first few elements of the Tensor's data array.
# @param x An instance of class 'Tensor'.
# @param ... Additional arguments (like 'n' for the number of items) that can be passed to the standard 'head' function.
setMethod(
    f = "head",
    signature = "Tensor",
    definition = function(x, ...) {
        # Ensure 'x' is a valid Tensor object
        validObject(x)

        # Call the base 'head' function on the 'data' slot of the Tensor
        # Additional arguments are passed through using '...'
        result <- head(x@data, ...)

        # Return the result of the 'head' function
        return(result)
    }
)

# Define the 'tail' method for objects of class 'Tensor'
# This method displays the last few elements of the Tensor's data array.
# @param x An instance of class 'Tensor'.
# @param ... Additional arguments (like 'n' for the number of items) that can be passed to the standard 'tail' function.
setMethod(
    f = "tail",
    signature = "Tensor",
    definition = function(x, ...) {
        # Ensure 'x' is a valid Tensor object
        validObject(x)

        # Call the base 'tail' function on the 'data' slot of the Tensor
        # Additional arguments are passed through using '...'
        result <- tail(x@data, ...)

        # Return the result of the 'tail' function
        return(result)
    }
)


# Define the subsetting method '[' for objects of class 'Tensor'
# This method allows subsetting of Tensor objects, retaining or dropping dimensions based on the 'drop' argument.
# @param x A 'Tensor' object.
# @param i Indices or conditions for subsetting along the first dimension.
# @param j Indices or conditions for subsetting along the second dimension.
# @param ... Additional indices or conditions for subsetting further dimensions.
# @param drop Logical, determines whether dimensions that have only one level should be dropped.
# @return Returns a subsetted 'Tensor' object.
setMethod(
    f = "[",
    signature = "Tensor",
    definition = function(x, i, j, ..., drop = TRUE) {
        # Ensure 'x' is a valid Tensor object
        validObject(x)

        # Use base subsetting function but control the 'drop' behavior explicitly
        if (!drop) {
            # If 'drop' is FALSE, use 'as.tensor' to ensure the result is still a 'Tensor' object
            # and retain original dimensions without dropping
            new_data <- as.tensor(`[`(x@data, i, j, ..., drop = FALSE), drop = drop)
        } else {
            # If 'drop' is TRUE, allow the default behavior and convert the result back to a 'Tensor'
            new_data <- as.tensor(`[`(x@data, i, j, ...))
        }

        # Return the new, possibly subsetted 'Tensor' object
        return(new_data)
    }
)

# Define the replacement method '[<-' for objects of class 'Tensor'
# This method modifies elements of the Tensor object based on specified indices and assigns new values to these elements.
# @param x A 'Tensor' object to be modified.
# @param i Indices or conditions for selecting elements along the first dimension.
# @param j Indices or conditions for selecting elements along the second dimension.
# @param ... Additional indices or conditions for further dimensions.
# @param value Values to be assigned to the selected elements.
# @return Returns the modified 'Tensor' object with the new values assigned.
setMethod(
    f = "[<-",
    signature = "Tensor",
    definition = function(x, i, j, ..., value) {
        # Ensure 'x' is a valid Tensor object
        validObject(x)

        # Perform the replacement operation on the 'data' slot of the Tensor
        new_data <- `[<-`(x@data, i, j, ..., value = value)

        # Convert the modified data back to a 'Tensor' object
        # Assuming 'as.tensor' is a function designed to convert arrays to Tensor objects, handling Tensor-specific properties
        updated_tensor <- as.tensor(new_data)

        # Return the modified Tensor object
        return(updated_tensor)
    }
)

# Define the transpose method 't' for objects of class 'Tensor'
# This method performs a transpose operation specifically designed for 3D tensors.
# @param x A 'Tensor' object expected to be three-dimensional.
# @return Returns a new 'Tensor' object with the dimensions swapped.
setMethod(
  f = "t",
  signature = "Tensor",
  definition = function(x) {
    # Ensure 'x' is a valid Tensor object
    validObject(x)

    # Check if the tensor is 3-dimensional as required
    if (x@num_modes != 3) {
      stop("Tensor Transpose is currently only implemented for 3D tensors.")
    }

    # Retrieve the modes (dimensions) of the tensor
    modes <- x@modes

    # Perform the transpose operation on the tensor's data
    # The 'apply' function is used to apply the transpose (t) function across the 3rd dimension
    new_arr <- array(apply(x@data, MARGIN = c(1, 3), FUN = function(mat) t(mat)), dim = modes[c(2, 1, 3)])

    # Assuming 'as.tensor' properly initializes a new Tensor object with given data and updated dimensions
    as.tensor(new_arr)
  }
)

