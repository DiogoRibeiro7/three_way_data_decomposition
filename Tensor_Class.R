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

# Overload elementwise operators for the 'Tensor' class
# This functionality enhances the ability to perform basic arithmetic operations (such as addition, subtraction,
# multiplication, and division) between Tensors themselves and between Tensors and other data types (arrays and numerics).
# These methods ensure that the operations are compatible with the dimensions and structure of the Tensor objects,
# promoting data integrity and ease of mathematical operations in data analysis and scientific computing contexts.

# Tensor and Tensor operation
# Performs the specified arithmetic operation (like +, -, *, /) between two Tensors.
# Both tensors must have compatible dimensions (i.e., their structures must allow the elementwise operation).
# @param e1 A Tensor object.
# @param e2 Another Tensor object.
# @return Returns a Tensor object after performing the specified operation.
setMethod("Ops", signature(e1 = "Tensor", e2 = "Tensor"),
  definition = function(e1, e2) {
    # Apply the operation to the data of both tensors
    e1@data <- callGeneric(e1@data, e2@data)
    # Validate the result to ensure structural integrity
    validObject(e1)
    return(e1)
  }
)

# Tensor and Array operation
# Performs the specified arithmetic operation between a Tensor and an array.
# The array must be conformable to the Tensor's dimension.
# @param e1 A Tensor object.
# @param e2 An array that can be broadcast to match the Tensor's dimensions.
# @return Returns the modified Tensor object.
setMethod("Ops", signature(e1 = "Tensor", e2 = "array"),
  definition = function(e1, e2) {
    # Apply the operation between tensor data and array
    e1@data <- callGeneric(e1@data, e2)
    # Validate to ensure no corruption of the tensor's structure
    validObject(e1)
    return(e1)
  }
)

# Array and Tensor operation
# Performs the specified operation between an array and a Tensor, with the array as the left-hand operand.
# Array dimensions must be conformable to those of the Tensor.
# @param e1 An array.
# @param e2 A Tensor object.
# @return Returns the modified Tensor object.
setMethod("Ops", signature(e1 = "array", e2 = "Tensor"),
  definition = function(e1, e2) {
    # Apply the operation with the array as the left operand
    e2@data <- callGeneric(e1, e2@data)
    # Validate the Tensor after operation
    validObject(e2)
    return(e2)
  }
)

# Tensor and Numeric operation
# Performs the specified arithmetic operation between a Tensor and a numeric value.
# The operation is broadcast across all elements of the Tensor.
# @param e1 A Tensor object.
# @param e2 A numeric value.
# @return Returns the modified Tensor object.
setMethod("Ops", signature(e1 = "Tensor", e2 = "numeric"),
  definition = function(e1, e2) {
    # Apply the operation to each element of the tensor data
    e1@data <- callGeneric(e1@data, e2)
    # Validate the tensor to ensure data consistency
    validObject(e1)
    return(e1)
  }
)

# Numeric and Tensor operation
# Performs the specified arithmetic operation between a numeric value and a Tensor,
# with the numeric value as the left-hand operand.
# This operation broadcasts the numeric value across all elements of the Tensor.
# @param e1 A numeric value.
# @param e2 A Tensor object.
# @return Returns the modified Tensor object.
setMethod("Ops", signature(e1 = "numeric", e2 = "Tensor"),
  definition = function(e1, e2) {
    # Numeric value affects each element of the tensor data
    e2@data <- callGeneric(e1, e2@data)
    # Validate to check the integrity of the tensor post-modification
    validObject(e2)
    return(e2)
  }
)


# Define the 'modeSum' method for the 'Tensor' class
# This method computes the sum of elements along a specified mode (dimension) of a tensor.
# The resulting tensor has the dimension of the specified mode reduced to 1, essentially collapsing that dimension.
#
# @param tnsr A 'Tensor' object on which the sum is computed.
# @param m The mode (dimension) over which to compute the sum.
#   Must be a valid mode within the tensor's dimensions.
# @param drop Logical, indicating whether to drop singleton dimensions after the operation.
# @return Returns a new 'Tensor' object with the specified mode summed and potentially dimensions dropped.
# @examples
# tnsr <- rand_tensor(c(3, 4, 5))  # Random tensor initialization for example purpose
# mode_summed_tensor <- modeSum(tnsr, m=2)
# @rdname modeSum-methods
# @aliases modeSum,Tensor-method
setMethod("modeSum", signature = "Tensor",
  definition = function(tnsr, m = NULL, drop = FALSE) {
    # Ensure the mode 'm' is specified
    if (is.null(m)) {
      stop("must specify mode 'm'")
    }

    # Retrieve the number of modes (dimensions) of the tensor
    num_modes <- tnsr@num_modes

    # Check if 'm' is within the bounds of the tensor's dimensions
    if (m < 1 || m > num_modes) {
      stop("mode 'm' out of bounds")
    }

    # Permute dimensions so the specified mode 'm' is the first dimension for easier summation
    perm <- c(m, (1L:num_modes)[-m])
    modes <- tnsr@modes

    # Calculate the new dimensions of the tensor after summing over mode 'm'
    newmodes <- modes
    newmodes[m] <- 1

    # Perform the sum over the specified mode
    arr <- array(colSums(aperm(tnsr@data, perm), dims = 1L), dim = newmodes)

    # Convert the resulting array back to a tensor, optionally dropping singleton dimensions
    as.tensor(arr, drop = drop)
  }
)


# Define the 'modeMean' method for the 'Tensor' class
# This method calculates the mean of elements along a specified mode (dimension) of a tensor.
# The resulting tensor has the dimension of the specified mode reduced to 1,
# collapsing that dimension by averaging instead of summing.
#
# @param tnsr A 'Tensor' object on which the mean is computed.
# @param m The mode (dimension) over which to compute the mean.
#   Must be a valid mode within the tensor's dimensions.
# @param drop Logical, indicating whether to drop singleton dimensions after the operation.
# @return Returns a new 'Tensor' object with the specified mode averaged and potentially dimensions dropped.
# @examples
# tnsr <- rand_tensor(c(3, 4, 5))  # Random tensor initialization for example purpose
# mode_mean_tensor <- modeMean(tnsr, m=2)
# @rdname modeMean-methods
# @aliases modeMean,Tensor-method
setMethod("modeMean", signature = "Tensor",
  definition = function(tnsr, m = NULL, drop = FALSE) {
    # Ensure the mode 'm' is specified
    if (is.null(m)) {
      stop("must specify mode 'm'")
    }

    # Retrieve the number of modes (dimensions) of the tensor
    num_modes <- tnsr@num_modes

    # Check if 'm' is within the bounds of the tensor's dimensions
    if (m < 1 || m > num_modes) {
      stop("mode 'm' out of bounds")
    }

    # Permute dimensions so the specified mode 'm' is the first dimension for easier averaging
    perm <- c(m, (1L:num_modes)[-m])
    modes <- tnsr@modes

    # Calculate the new dimensions of the tensor after averaging over mode 'm'
    newmodes <- modes
    newmodes[m] <- 1

    # Perform the averaging over the specified mode
    sum_arr <- array(colSums(aperm(tnsr@data, perm), dims = 1L), dim = newmodes)
    mean_arr <- sum_arr / modes[m]  # Divide by the size of the dimension to get the mean

    # Convert the resulting array back to a tensor, optionally dropping singleton dimensions
    as.tensor(mean_arr, drop = drop)
  }
)


# Define the 'fnorm' method for the 'Tensor' class
# This method calculates the Frobenius norm of a tensor, which is equivalent to the Euclidean norm
# of a vector but generalized for matrices and higher-dimensional arrays.
# The Frobenius norm is computed as the square root of the sum of the absolute squares of its elements.
#
# @param tnsr A 'Tensor' object whose Frobenius norm is to be calculated.
# @return Returns the Frobenius norm as a numeric value.
# @examples
# tnsr <- rand_tensor(c(3, 4, 5))  # Assume a function to generate a random tensor
# norm_value <- fnorm(tnsr)
# @rdname fnorm-methods
# @aliases fnorm,Tensor-method
setMethod("fnorm", signature = "Tensor",
  definition = function(tnsr) {
    # Retrieve the data from the Tensor object
    arr <- tnsr@data

    # Compute the Frobenius norm: sqrt of the sum of squares of all elements
    frob_norm <- sqrt(sum(arr * arr))

    # Return the computed Frobenius norm
    return(frob_norm)
  }
)


# Define the 'innerProd' method for the 'Tensor' class
# This method calculates the inner product of two tensors. For the operation to be valid,
# both tensors must have the same dimensions. The inner product is the sum of the products
# of their corresponding elements.
#
# @param tnsr1 The first 'Tensor' object.
# @param tnsr2 The second 'Tensor' object, which must have the same dimensions as tnsr1.
# @return Returns the inner product as a numeric value.
# @examples
# tnsr1 <- rand_tensor(c(3, 4, 5))  # Assume a function to generate a random tensor
# tnsr2 <- rand_tensor(c(3, 4, 5))
# product_value <- innerProd(tnsr1, tnsr2)
# @rdname innerProd-methods
# @aliases innerProd,Tensor,Tensor-method
setMethod("innerProd", signature = c(tnsr1 = "Tensor", tnsr2 = "Tensor"),
  definition = function(tnsr1, tnsr2) {
    # Ensure that both tensors have the same dimensions
    if (!all(tnsr1@modes == tnsr2@modes)) {
      stop("The dimensions of the two tensors must match for the inner product.")
    }

    # Retrieve the data from both tensors
    arr1 <- tnsr1@data
    arr2 <- tnsr2@data

    # Calculate the inner product: sum of element-wise products
    inner_product <- sum(arr1 * arr2)

    # Return the computed inner product
    return(inner_product)
  }
)


# Define the 'unfold' method for the 'Tensor' class
# This method reorganizes a tensor into a matrix based on specified indices for the rows and columns.
# It requires the specification of indices for the rows and columns that together must cover all modes
# of the tensor without overlap and repetition.
#
# @param tnsr A 'Tensor' object to be unfolded.
# @param row_idx A vector of integers indicating which modes of the tensor should form the rows of the matrix.
# @param col_idx A vector of integers indicating which modes of the tensor should form the columns of the matrix.
# @return Returns a new 'Tensor' object that is a two-dimensional matrix representation of the original tensor.
# @examples
# tnsr <- rand_tensor(c(3, 4, 5))  # Assume a function to generate a random tensor
# unfolded_tensor <- unfold(tnsr, row_idx=c(1,2), col_idx=3)
# @rdname unfold-methods
# @aliases unfold,Tensor-method
setMethod("unfold", signature = "Tensor",
  definition = function(tnsr, row_idx = NULL, col_idx = NULL) {
    # Validation of indices
    if (is.null(row_idx) || is.null(col_idx)) {
      stop("Both row and column indices must be specified.")
    }
    num_modes <- tnsr@num_modes
    if (length(row_idx) + length(col_idx) != num_modes) {
      stop("The total number of row and column indices must match the number of modes in the tensor.")
    }
    if (any(row_idx < 1) || any(row_idx > num_modes) || any(col_idx < 1) || any(col_idx > num_modes)) {
      stop("Specified indices are out of bounds.")
    }
    perm <- c(row_idx, col_idx)
    if (any(sort(perm, decreasing = TRUE) != num_modes:1)) {
      stop("Indices are missing and/or repeated. Each index should be used exactly once.")
    }

    # Retrieve modes and reorganize tensor data
    modes <- tnsr@modes
    mat <- tnsr@data

    # Calculate new dimensions for the matrix
    new_modes <- c(prod(modes[row_idx]), prod(modes[col_idx]))

    # Rearrange the tensor into a matrix according to the specified permutation of modes
    mat <- aperm(mat, perm)
    dim(mat) <- new_modes

    # Return the newly formed matrix as a tensor
    as.tensor(mat)
  }
)


# Define the 'k_unfold' method for the 'Tensor' class
# This method unfolds the tensor along a specified mode 'm'. In mode-k unfolding,
# the mode 'm' forms the rows of the resulting matrix, and the other modes are combined to form the columns.
# This operation is commonly used in tensor decomposition techniques and multilinear algebra.
#
# @param tnsr A 'Tensor' object to be unfolded.
# @param m An integer specifying the mode that should form the rows of the unfolded matrix.
# @return Returns a new 'Tensor' object that is a two-dimensional matrix, where the specified mode is unfolded
#         to the rows and all other modes to the columns.
# @examples
# tnsr <- rand_tensor(c(3, 4, 5))  # Assume a function to generate a random tensor
# unfolded_tensor <- k_unfold(tnsr, m=2)
# @rdname k_unfold-methods
# @aliases k_unfold,Tensor-method
setMethod("k_unfold", signature = "Tensor",
  definition = function(tnsr, m = NULL) {
    # Ensure the mode 'm' is specified
    if (is.null(m)) {
      stop("Mode 'm' must be specified for k_unfolding.")
    }

    # Check if 'm' is within the bounds of tensor's dimensions
    num_modes <- tnsr@num_modes
    if (m < 1 || m > num_modes) {
      stop("Specified mode 'm' is out of bounds.")
    }

    # Define the row indices as 'm' and the column indices as all other modes
    rs <- m
    cs <- (1:num_modes)[-m]

    # Call the 'unfold' method to perform the unfolding
    unfolded_tensor <- unfold(tnsr, row_idx = rs, col_idx = cs)
    return(unfolded_tensor)
  }
)


# Define the 'matvec' method for the 'Tensor' class
# This method is specifically designed for 3D tensors where it rearranges the tensor into a matrix form.
# In this matricization, the tensor is unfolded such that the first and third modes form the rows,
# and the second mode forms the columns of the resulting matrix. This kind of unfolding is often
# used in operations that treat the tensor as a collection of vectors (mode-2 vectors) for efficient
# computation of matrix-vector products.
#
# @param tnsr A 'Tensor' object to be transformed.
# @return Returns a new 'Tensor' object that is a two-dimensional matrix, unfolded such that
#         modes 1 and 3 form the rows, and mode 2 forms the columns.
# @examples
# tnsr <- rand_tensor(c(3, 4, 5))  # Assume a function to generate a random tensor
# matrix_representation <- matvec(tnsr)
# @rdname matvec-methods
# @aliases matvec,Tensor-method matvec,Tensor,Tensor-method
setMethod('matvec', signature = "Tensor",
          definition = function(tnsr) {
            # Ensure the tensor is exactly 3-dimensional
            if (tnsr@num_modes != 3) {
              stop("Matvec currently only implemented for 3D Tensors")
            }

            # Perform the unfolding using predefined dimensions
            # This unfolds the tensor such that the rows are combinations of modes 1 and 3,
            # and the columns are formed by mode 2.
            result <- unfold(tnsr, row_idx = c(1, 3), col_idx = 2)

            # Return the resulting unfolded tensor
            return(result)
          }
)


# Define the 'rs_unfold' method for the 'Tensor' class
# This method unfolds the tensor along a specified mode 'm', arranging the tensor such that the specified
# mode forms the rows of the resulting matrix, and all other modes are combined to form the columns.
# This type of unfolding is useful in various multilinear operations and tensor decompositions where
# specific modes are treated distinctly from others.
#
# @param tnsr A 'Tensor' object to be unfolded.
# @param m An integer specifying the mode that should form the rows of the unfolded matrix.
# @return Returns a new 'Tensor' object that is a two-dimensional matrix, where the specified mode is unfolded
#         to the rows and all other modes to the columns.
# @examples
# tnsr <- rand_tensor(c(3, 4, 5))  # Assume a function to generate a random tensor
# unfolded_tensor <- rs_unfold(tnsr, m=1)
# @rdname rs_unfold-methods
# @aliases rs_unfold,Tensor-method
setMethod("rs_unfold", signature = "Tensor",
          definition = function(tnsr, m = NULL) {
              # Ensure the mode 'm' is specified
              if (is.null(m)) {
                  stop("Mode 'm' must be specified for rs_unfolding.")
              }

              # Validate mode 'm' within the bounds of tensor's dimensions
              num_modes <- tnsr@num_modes
              if (m < 1 || m > num_modes) {
                  stop("Specified mode 'm' is out of bounds.")
              }

              # Define the row indices as 'm' and the column indices as all other modes
              rs <- m
              cs <- (1:num_modes)[-m]

              # Call the 'unfold' method to perform the unfolding
              unfolded_tensor <- unfold(tnsr, row_idx = rs, col_idx = cs)
              return(unfolded_tensor)
          }
)


# Define the 'cs_unfold' method for the 'Tensor' class
# This method unfolds the tensor along a specified mode 'm', arranging the tensor such that the specified
# mode forms the columns of the resulting matrix, and all other modes are combined to form the rows.
# This type of unfolding is particularly useful for analyses that require treating one mode as the feature dimension
# in statistical models or machine learning algorithms.
#
# @param tnsr A 'Tensor' object to be unfolded.
# @param m An integer specifying the mode that should form the columns of the unfolded matrix.
# @return Returns a new 'Tensor' object that is a two-dimensional matrix, where the specified mode is unfolded
#         to the columns and all other modes to the rows.
# @examples
# tnsr <- rand_tensor(c(3, 4, 5))  # Assume a function to generate a random tensor
# unfolded_tensor <- cs_unfold(tnsr, m=3)
# @rdname cs_unfold-methods
# @aliases cs_unfold,Tensor-method
setMethod("cs_unfold", signature = "Tensor",
          definition = function(tnsr, m = NULL) {
              # Ensure the mode 'm' is specified
              if (is.null(m)) {
                  stop("Mode 'm' must be specified for cs_unfolding.")
              }

              # Validate mode 'm' within the bounds of tensor's dimensions
              num_modes <- tnsr@num_modes
              if (m < 1 || m > num_modes) {
                  stop("Specified mode 'm' is out of bounds.")
              }

              # Define the column indices as 'm' and the row indices as all other modes
              rs <- (1:num_modes)[-m]
              cs <- m

              # Call the 'unfold' method to perform the unfolding
              unfolded_tensor <- unfold(tnsr, row_idx = rs, col_idx = cs)
              return(unfolded_tensor)
          }
)


# Creation of Tensor from array, matrix, or vector

#' Tensor Conversion
#'
#' Create a \code{\link{Tensor-class}} object from an \code{array}, \code{matrix}, or \code{vector}.
#' This function allows for the seamless transformation of base R data structures into a Tensor-class object,
#' facilitating operations that are specific to tensors.
#' 
#' @param x An instance of \code{array}, \code{matrix}, or \code{vector} to be converted into a Tensor.
#' @param drop Logical; if TRUE, any singleton dimensions (dimensions of length 1) will be dropped in the resulting tensor.
#' @return A \code{\link{Tensor-class}} object representing the input data structured as a tensor.
#' @export
#' @name as.tensor
#' @rdname as.tensor
#' @aliases as.tensor
#' @examples
#' # From vector:
#' vec <- runif(100)
#' vecT <- as.tensor(vec)
#' vecT
#'
#' # From matrix:
#' mat <- matrix(runif(1000), nrow = 100, ncol = 10)
#' matT <- as.tensor(mat)
#' matT
#'
#' # From array:
#' indices <- c(10, 20, 30, 40)
#' arr <- array(runif(prod(indices)), dim = indices)
#' arrT <- as.tensor(arr)
#' arrT
as.tensor <- function(x, drop = FALSE) {
    # Validate that input x is an array, matrix, or vector
    stopifnot(is.array(x) || is.vector(x))

    # Define the modes and number of modes based on the structure of x
    if (is.vector(x)) {
        modes <- c(length(x))
        num_modes <- 1L
    } else {
        modes <- dim(x)
        num_modes <- length(modes)
        
        # Handle dropping of singleton dimensions
        if (drop) {
            dim1s <- which(modes == 1)
            if (length(dim1s) > 0) {
                modes <- modes[-dim1s]
                num_modes <- num_modes - length(dim1s)
            }
        }
    }

    # Create a new Tensor-class object with the specified modes and data
    new("Tensor", num_modes = num_modes, modes = modes, data = array(x, dim = modes))
}


# Set up a generic function for 'tperm', the tensor permutation method
setGeneric("tperm",
           def = function(tnsr, perm, ...) {
               standardGeneric("tperm")
           })

# Define the 'tperm' method for the 'Tensor' class
# This method simplifies mode permutation of a tensor by overloading the 'aperm' function,
# making it specific and convenient for tensor manipulation.
#
# @param tnsr A 'Tensor' object whose modes are to be permuted.
# @param perm A vector specifying the new ordering of the modes.
# @param ... Additional arguments to be passed to 'aperm', such as 'useNames'.
# @return Returns a new 'Tensor' object with its modes permuted according to 'perm'.
# @examples
# tnsr <- rand_tensor(c(3, 4, 5))  # Assume a function to generate a random tensor
# permuted_tensor <- tperm(tnsr, perm = c(2, 1, 3))
# permuted_tensor <- tperm(tnsr, perm = c(1, 3, 2))
# @rdname tperm-methods
# @aliases tperm tperm-methods tperm,Tensor-method
# @seealso \code{\link{aperm}}
setMethod("tperm", signature = "Tensor",
          definition = function(tnsr, perm, ...) {
              # Check if 'perm' is provided and valid
              if (is.null(perm)) {
                  stop("'perm' must be specified and cannot be NULL.")
              }
              num_modes <- length(tnsr@modes)
              if (length(perm) != num_modes || any(perm < 1) || any(perm > num_modes)) {
                  stop("Invalid permutation vector. Ensure it correctly specifies mode indices.")
              }

              # Apply permutation to the tensor's data
              permuted_data <- aperm(tnsr@data, perm, ...)

              # Return a new tensor with the permuted data
              return(as.tensor(permuted_data, modes = tnsr@modes[perm]))
          }
)


# Set up a generic function for 'vec', which flattens a tensor into a vector
setGeneric("vec", def = function(tnsr) {
    standardGeneric("vec")
})

# Define the 'vec' method for the 'Tensor' class
# This method converts a tensor into a vector, following the convention that
# earlier indices in the tensor vary slower than later indices. This is commonly used
# in operations where a tensor needs to be treated as a single sequence of elements,
# such as serialization or certain types of mathematical analysis.
#
# @param tnsr A 'Tensor' object to be flattened.
# @return Returns a vector containing all elements of the tensor.
# @examples
# tnsr <- rand_tensor(c(4, 5, 6, 7))  # Assume a function to generate a random tensor
# flattened_vector <- vec(tnsr)
# @rdname vec-methods
# @aliases vec vec,Tensor-method
# @references T. Kolda, B. Bader, "Tensor decomposition and applications".
#            SIAM Applied Mathematics and Applications, 2009.
setMethod("vec", signature = "Tensor",
          definition = function(tnsr) {
              # Flatten the tensor data into a vector
              as.vector(tnsr@data)
          }
)

