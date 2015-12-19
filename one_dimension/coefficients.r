# Defines several coefficents for the drift and diffusion.

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

make_mean_reverting_coefficient <- function(mean_reversion_rate, long_term_mean)
    # Creates a mean reverting function of the form -mean_reversion_rate(x-long_term_mean)
    #
    # Args:
    #     mean_reversion_rate (numeric): The mean reverting speed
    #     long_term_mean (numeric): The value that the process reverts to.
    # Returns:
    #     (function(anything, numeric)): A function that returns -mean_reversion_rate(x-long_term_mean)
    function(t,x) -mean_reversion_rate * (x - long_term_mean)

mu_sin <- function(t, x)
    # The drift coefficient. Sin function.
    # Args:
    #     t (numeric): The current time
    #     x (array-like): The previous value of the process for each path
    sin(x)

mu_from_paper <- function(t, y){
    # The drift coefficient. Implementing equation 5.2 from the paper.
    #
    # Args:
    #     t (float): The current time
    #     y (TODO): The previous value of the process for each path
    sigma <- 0.4
    x <- convert_y_to_x_1d(y)
    return(2*sigma*x/(1+x^2)^2)
}
    

convert_y_to_x_1d <- function(y){
    # This formula converts Y_t from the Lamperti transformation back to X_t
    X0 <- 1
    sigma <- .4
    temp1 <- 27*X0^3 + 81*X0 + 162*sigma*y
    temp2 <- ((temp1^2 + 2916)^0.5 + temp1)^(1.0/3) / (3 * 2^(1.0/3))
    return(temp2 - 1/temp2)
}

convert_x_to_y_1d <- function(x){
    # This formula converts X_t to Y_t via the Lamperti transformation
    X0 <- 1
    sigma <- .4
    return((x-X0 + (x^3 - X0^3)/3) / (2*sigma))
}

#################################
## Diffusion Coefficients
#################################

sigma_from_paper <- function(t, x){
    # The diffusion coefficient. Implementing equation 5.1 from the paper.
    #
    # Args:
    #     t (numeric): The current time
    #     x (array-like): The previous value of the process for each path
    sigma <- 0.4
    return(2*sigma/(1+x^2))
}

sigma_deriv_from_paper <- function(t, x){
    # The derivative of the diffusion coefficient. Implementing equation 5.1 from the paper.
    #
    # Args:
    #     t (numeric): The current time
    #     x (array-like): The previous value of the process for each path
    sigma <- 0.4
    return(-4*x*sigma / (1+x^2)^2)
}
