beta_calculator <- function(L, T, mu0, sigma0){
    # Calculates the sub-optimal beta for the poisson process.
    #
    # Args:
    #    L (float): The Lipschitz/Holder coefficient decribed in Assumption 3.1 of the reference paper
    #    T (float): The maturity of the asset
    #    mu0 (np.array): The starting value of mu for each dimension (not sure about this)
    #    sigma0 (np.array): The starting value of sigma for each dimension (not sure about this)
    # Returns:
    #     (float): The sub-optimal beta
    
    # Get the dimension of the process
    d <- length(mu0)
    
    mu_inf <- sqrt(sum(mu0^2))
    
    trace <- sum(diag(solve(sigma0 %*% t(sigma0))))
        
    L_prime <- L^2 * ( 2*(1+mu_inf*sqrt(T))^2 * trace + 2*(3*d+d*(d-1)) )
    
    return(sqrt(L_prime + 0.25*T^2) + 0.5*T)
}

test <- function(L){
    beta_calculator(L, T, mu0, sigma0)
}
