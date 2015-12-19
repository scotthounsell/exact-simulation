# Defines several coefficents for the drift and diffusion.
# Also defines functions for conversions between processes related through a lamperti transform.

#################################
## General Coefficients
#################################

make_constant_coefficient <- function(const_val)
    # Creates a function that returns a constant.
    #
    # Args:
    #     const_val (numeric): The value that the produced function should return
    # Returns:
    #     (function(anything, anything)): A function that will always return const_val
    function(t,x) const_val

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

make_mean_reverting_coefficient <- function(reversion_rate, long_term_mean)
    # Creates a mean reverting function of the form -reversion_rate(x-long_term_mean)
    #
    # Args:
    #     reversion_rate (numeric): The mean reverting rate
    #     long_term_mean (numeric): The value that the process reverts to.
    # Returns:
    #     (function(anything, numeric)): A function that returns -reversion_rate(x-long_term_mean)
    function(t,x) -reversion_rate*(x - long_term_mean)

mu_1d_from_paper <- function(t, y){
    # The drift coefficient. Implementing equation 5.2 from the paper.
    #
    # Args:
    #     t (numeric): The current time
    #     y (numeric): The Lamperti-transformed process
    sigma <- 0.4
    x <- convert_y_to_x_1d(y)
    denom <- 1 + x*x
    return(2*sigma*x/(denom*denom))
}

convert_y_to_x_1d <- function(y){
    # Converts Y_t to X_t in equation 5.2
    #
    # Args:
    #     y (numeric): The Lamperti-transformed process
    X0 <- 1
    sigma <- 0.4
    temp <- (27*X0*X0 + 81)*X0 + 162*sigma*y
    temp <- (.5*(sqrt(temp*temp + 2916) + temp))^(1.0/3) / 3
    return(temp - 1/temp)
}

convert_x_to_y_1d <- function(x){
    # Converts X_t to Y_t in equation 5.2
    #
    # Args:
    #     x (numeric): The original process
    X0 <- 1
    sigma <- .4
    return((x-X0 + (x*x*x - X0*X0*X0)/3) / (2*sigma))
}

make_mu_section_52_euler <- function(d)
    # Returns the drift coefficient function described in section 5.2 from the paper.
    #
    # Args:
    #     d (int): The current dimension
    function(t, x){
        # Args:
        #     t (numeric): The current time
        #     x (matrix): The current value of the process across all paths for this dimension
        x <- matrix(x)
        x[d,]*0.1*(sqrt(x[d,]) - 1)
    }

make_mu_from_paper_multi_dimension_lamperti <- function(d)
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
        return(0.2*(sqrt(x[d])-1) - 0.25)
    }

convert_y_to_x_multi_d <- function(y){
    # Converts Y_t (lamperti) to X_t from section 5.2
    #
    # Args:
    #     y (matrix or vector or numeric): Values of Y_t potentially across dimensions and timesteps
    # Returns:
    #     (matrix or vector or numeric): Values of X_t corresponding to each Y_t
    X0 <- 1
    return(X0*exp(0.5*y))
}

convert_x_to_y_multi_d <- function(x){
    # Converts X_t to Y_t (lamperti) from section 5.2
    #
    # Args:
    #     x (matrix or vector or numeric): Values of X_t potentially across dimensions and timesteps
    # Returns:
    #     (matrix or vector or numeric): Values of Y_t corresponding to each X_t
    X0 <- 1
    return(2*log(x/X0))
}

make_mu_heston_variance <- function(lambda, vbar)
    # Returns the drift coefficient function of the lamperti transform
    # of the Heston variance process
    #
    # Args:
    #     lambda (numeric): The mean reversion rate
    #     vbar (numeric): The long term mean
    function(t, v){
        # Args:
        #     t (numeric): The current time
        #     v (vector or numeric): The current variance
        if (any(v<=0)) print('(mu) v went <=0 !!!') # FIXME
        return(-lambda*((v>0)*v - vbar)) # full truncation
    }

make_mu_heston_variance_lamperti <- function(lambda, vbar, eta, v0){
    # Returns the drift coefficient function of the lamperti transform
    # of the Heston variance process
    #
    # Args:
    #     lambda (numeric): The mean reversion rate
    #     vbar (numeric): The long term mean
    #     eta (numeric): The volatility of volatility
    #     v0 (numeric): The initial variance
    convert_y_to_x_heston_variance <- make_convert_y_to_x_heston_variance(eta, v0)
    function(t, y){
        # Args:
        #     t (numeric): The current time
        #     v (vector or numeric): The current variance
        v <- convert_y_to_x_heston_variance(y)
        if (any(v<=0)) print('(mu L) v went <=0 !!!') # FIXME
        v <- (v>0)*v # full truncation
        return((-lambda*(v-vbar)/eta - 0.25*eta) / sqrt(v))
    }
}

make_convert_y_to_x_heston_variance <- function(eta, v0)
    function(y){
        # Converts Y_t (lamperti) to v_t from the Heston variance process
        #
        # Args:
        #     y (matrix or vector or numeric): Values of Y_t potentially across dimensions and timesteps
        # Returns:
        #     (matrix or vector or numeric): Values of X_t corresponding to each Y_t
        rootv <- 0.5*eta*y + sqrt(v0)
        return(rootv*rootv)
    }

make_convert_x_to_y_heston_variance <- function(eta, v0)
    function(v)
        # Converts X_t to Y_t (lamperti) from the Heston variance process
        #
        # Args:
        #     x (matrix or vector or numeric): Values of X_t potentially across dimensions and timesteps
        # Returns:
        #     (matrix or vector or numeric): Values of Y_t corresponding to each X_t
        return(2*(sqrt(v)-sqrt(v0))/eta)

make_mu_black_scholes_lamperti <- function(mu, sigma)
    # Returns the drift coefficient function of the lamperti transform
    # of the Black Scholes process
    #
    # Args:
    #     TODO
    function(t, y)
        # Args:
        #     t (numeric): The current time
        #     y (matrix): The location of the Lamperti transformed process
        return(mu/sigma - 0.5*sigma)

make_convert_y_to_x_black_scholes <- function(sigma, X0)
    function(y)
        # Converts Y_t (lamperti) to v_t from the Black Scholes process
        #
        # Args:
        #     y (matrix or vector or numeric): Values of Y_t potentially across dimensions and timesteps
        # Returns:
        #     (matrix or vector or numeric): Values of X_t corresponding to each Y_t
        return(X0*exp(sigma*y))

make_convert_x_to_y_black_scholes <- function(sigma, X0)
    function(x)
        # Converts X_t to Y_t (lamperti) from the Black Scholes process
        #
        # Args:
        #     x (matrix or vector or numeric): Values of X_t potentially across dimensions and timesteps
        # Returns:
        #     (matrix or vector or numeric): Values of Y_t corresponding to each X_t
        return(log(x/X0)/sigma)
        
make_mu_sin <- function(a, k, phi, b)
    # Returns the drift coefficient function of the sin process
    # dX_t = (b + a*sin(n*pi*x))*dt + sigma*dW_t
    #
    # Args:
    #     a (numeric): The amplitude of the oscillatory process
    #     n (numeric): The frequency of the oscillatory process
    #     b (numeric): The drift of the sine process
    function(t, x)
        # Args:
        #     t (numeric): The current time
        #     x (matrix): The location of the process
        return(a*(sin(k*t + phi) + b))


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
    return(2*sigma/(1+x*x))
}

sigma_deriv_from_paper <- function(t, x){
    # The derivative of the diffusion coefficient. Implementing equation 5.1 from the paper.
    #
    # Args:
    #     t (numeric): The current time
    #     x (array_like): The previous value of the process for each path
    sigma <- 0.4
    denom <- 1 + x*x
    return(-4*x*sigma / (denom*denom))
}

make_sigma_section_52_euler <- function(d)
    # Returns the diffusion coefficient function described in section 5.2 from the paper.
    #
    # Args:
    #     d (int): The current dimension
    function(t, x){
        # Args:
        #     t (numeric): The current time
        #     x (matrix): The current value of the process across all paths for this dimension
        x <- matrix(x)
        return(0.5*x[d,])
    }

make_sigma_cir <- function(eta)
    # Returns the CIR diffusion coefficient function eta*sqrt(X_t)
    #
    # Args:
    #     eta (numeric): The diffusion scaling parameter.
    function(t, x){
        # Args:
        #     t (numeric): The current time
        #     x (matrix): The ending value of the process for each dimension and path
        if (any(x<=0)) print('(CIR) x went <=0 !!!') # FIXME
        return(eta*sqrt((x>0)*x)) # full truncation
    }
