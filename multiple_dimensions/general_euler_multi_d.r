library(MASS)

General_Euler_Multi_D <- function(num_of_paths, num_of_steps, mu, sigma, cov_matrix, T=1, X0=0, convert_y_to_x_func=function(y) y){
    # Simulate a 1-dimensional process using Monte Carlo with the standard Euler scheme.
    # This can be used to simulate the payoff of an option and calculate the implied volatility smile.
    #
    # Args:
    #     num_of_paths (int): The number of Monte Carlo paths
    #     num_of_steps (int): The number of steps in each Monte Carlo path
    #     mu (vector of functions): A vector of the drift coefficient functions, a function of t (float) and X (array-like) returning a scalar
    #     sigma (vector of functions): A vector of the diffusion coefficient, a function of t (float) and X (array-like) returning a scalar
    #     cov_matrix (matrix of floats): A covariance matrix of the diffusion
    #     T (Optional[float]): Time, in years, to simulate process. Defaults to 1.
    #     X0 (vector of floats): Initial values for each of the process, Xhat. Defaults to 0.
    #     convert_y_to_x (function): A function to convert Y_t to X_t. Defaults to X_t <- Y_t
    # Returns:
    #     (list of function): A list with a single entry - the run_monte_carlo function stored in the key run_monte_carlo
    
    # Check that the inputs are correct dimensions
    if (!is.vector(mu)) stop('mu must be a vector')
    if (!is.vector(sigma)) stop('sigma must be a vector')
    num_dim <- length(mu)
    if (num_dim != length(sigma)) stop('mu and sigma must be the same length')
    if (nrow(cov_matrix) != num_dim) stop('cov_matrix must be the same size as mu and sigma')
    
    dt <- T / num_of_steps

    run_monte_carlo <- function(g_maker, strikes){
        # Runs the Monte Carlo simulation to calculate E[g(X_T)] for different strikes.
        #
        # Args:
        #     g_maker (function): A function that returns a function f(x) <- g(x,K)
        #     strikes (array-like): An array of strikes to pass in to g_maker
        # Returns:
        #     (array-like): A array of E[g(X_T)] for the range of strikes

        # Simulates the process X across time.
        # X is a 3-dimensional array-like:
        #     -The first dimension represents the dimensions of the process
        #     -The second dimension represents the timesteps.
        #     -The third dimension represents the paths. Each "sheet" is a path.
        
        # Generate dW and put it into the correct 3-dimensional form described above
        dW <- mvrnorm(n=num_of_paths*num_of_steps, mu=rep(0, num_dim), Sigma=dt*cov_matrix)
        dW <- t(dW) # Transposed so that time goes across cols instead of rows
        dim(dW) <- c(num_dim, num_of_steps, num_of_paths)
        
        # Generate X
        X <- array(rep(NA, num_dim*(num_of_steps+1)*num_of_paths), c(num_dim, num_of_steps+1, num_of_paths))
        X[,1,] <- X0
        for (t in 1:num_of_steps)
            # for (d in 1:num_dim)
            #     X[d,t+1,] <- X[d,t,] + mu[[d]]((t-1)*dt, X[,t,]) * dt + sigma[[d]]((t-1)*dt, X[,t,]) * dW[d,t,]
            # Hard code for Section 5.2
            X[,t+1,] <- X[,t,]*exp(0.5*dW[,t,] + (0.1*(sqrt(X[,t,])-1) - 1/8)*dt)
        
        #X_Ts <- sapply(X[,ncol(X),], convert_y_to_x_func)
        X_Ts <- X[,ncol(X),]

        return(sapply(strikes, function(K) mean(g_maker(K)(X_Ts))))
    }
    return(list(run_monte_carlo=run_monte_carlo)) # Return a list holding the function so we can call it oop style
}
