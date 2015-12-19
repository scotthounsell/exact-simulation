# Defines several coefficents for the drift and diffusion.
# Also defines functions for conversions between processes related through a lamperti transform.

#################################
## General Coefficients
#################################

make_constant_coefficient <- function(constant_value)
    # Creates a function that returns a constant.
    #
    # Args:
    #     constant_value (numeric): The value that the produced function should return
    # Returns:
    #     (function(anything, anything)): A function that will always return constant_value
    function(t,x) constant_value

make_affine_coefficient <- function(a, b)
    # Creates a affine function a*x+b.
    #
    # Args:
    #     a (numeric): The constant term of the affine function
    #     b (numeric): The first order term of the affine function
    # Returns:
    #     (function(anything, numeric)): A function that returns a*x+b
    function(t,x) a*x + b

#################################
## Drift Coefficients
#################################

make_mean_reverting_coefficient <- function(mean_reverting_rate, long_term_mean)
    # Creates a mean reverting function of the form -mean_reverting_rate(x-long_term_mean)
    #
    # Args:
    #     mean_reverting_rate (numeric): The mean reverting rate
    #     long_term_mean (numeric): The value that the process reverts to.
    # Returns:
    #     (function(anything, numeric)): A function that returns -mean_reverting_rate(x-long_term_mean)
    function(t,x) -mean_reverting_rate * (x - long_term_mean)

mu_from_paper <- function(t, y){
    # The drift coefficient. Implementing equation 5.2 from the paper.
    #
    # Args:
    #     t (numeric): The current time
    #     y (TODO): The previous value of the process for each path
    sigma <- 0.4
    x <- convert_y_to_x_1d(y)
    return(2*sigma*x/(1+x^2)^2)
}

make_mu_section_52_euler <- function(d){
    # Returns the drift coefficient function described in section 5.2 from the paper.
    #
    # Args:
    #     d (int): The current dimension
    function(t, x){
        # Args:
        #     t (numeric): The current time
        #     x (matrix): The ending value of the process for each dimension and path
        x <- matrix(x)
        x[d,]*0.1*(sqrt(x[d,]) - 1)
    }
}

convert_y_to_x_1d <- function(y){
    # Converts Y_t to X_t in equation 5.2
    #
    # Args:
    #     y (TODO): TODO
    X0 <- 1
    sigma <- 0.4
    temp1 <- 27*X0^3 + 81*X0 + 162*sigma*y
    temp2 <- ((temp1^2 + 2916)^0.5 + temp1)^(1.0/3) / (3 * 2^(1.0/3))
    return(temp2 - 1/temp2)
}

convert_x_to_y_1d <- function(x){
    # Converts X_t to Y_t in equation 5.2
    #
    # Args:
    #     x (TODO): TODO
    X0 <- 1
    sigma <- .4
    return((x-X0 + (x^3 - X0^3)/3) / (2*sigma))
}

make_mu_from_paper_multi_dimension_lamperti <- function(d){
    # Returns the drift coefficient function of the lamperti transform 
    # of the process described in section 5.2 from the paper.
    #
    # Args:
    #     d (int): The current dimension
    function(t, y){
        # Args:
        #     t (numeric): The current time
        #     y (matrix): The ending value of the process for each dimension and path
        X0 <- 1
        x <- convert_y_to_x_multi_d(y)
        return(0.2*X0*(sqrt(x[d])-1) - 0.25*X0)
    }
}

convert_y_to_x_multi_d <- function(y){
    # Converts Y_t (lamperti) to X_t from section 5.2
    #
    # Args:
    #     y (matrix or vector or numeric): Values of Y_t potentially across dimensions and timesteps
    # Returns:
    #     (matrix or vector or numeric): Values of X_t corresponding to each Y_t
    X0 <- 1.0
    return(X0*exp(0.5*y))
}

convert_x_to_y_multi_d <- function(x){
    # Converts X_t to Y_t (lamperti) from section 5.2
    #
    # Args:
    #     x (matrix or vector or numeric): Values of X_t potentially across dimensions and timesteps
    # Returns:
    #     (matrix or vector or numeric): Values of Y_t corresponding to each X_t
    X0 <- 1.0
    return(2*log(x/X0))
}

#################################
## Diffusion Coefficients
#################################

sigma_from_paper <- function(t, x){
    # The diffusion coefficient. Implementing equation 5.1 from the paper.
    #
    # Args:
    #     t (numeric): The current time
    #     x (array_like): The previous value of the process for each path
    sigma <- 0.4
    return(2*sigma/(1+x^2))
}

sigma_deriv_from_paper <- function(t, x){
    # The derivative of the diffusion coefficient. Implementing equation 5.1 from the paper.
    #
    # Args:
    #     t (numeric): The current time
    #     x (array_like): The previous value of the process for each path
    sigma <- 0.4
    return(-4*x*sigma / (1+x^2)^2)
}

make_sigma_section_52_euler <- function(d){
    function(t, x){
        x <- matrix(x)
        0.5*x[d,]
    }
}
