General_Euler <- function(num_of_paths, num_of_steps, mu, sigma, T=1, X0=0, convert_y_to_x=function(y) y, convert_x_to_y=function(x) x){
    # Simulate a 1-dimensional process using Monte Carlo with the standard Euler scheme.
    # This can be used to simulate the payoff of an option and calculate the implied volatility smile.
    #
    # Args:
    #     num_of_paths (integer): The number of Monte Carlo paths
    #     num_of_steps (integer): The number of steps in each Monte Carlo path
    #     mu (function): The drift coefficient, a function of t (numeric) and X (array-like) returning a scalar
    #     sigma (function): The diffusion coefficient, a function of t (numeric) and X (array-like) returning a scalar
    #     T (Optional[numeric]): Time, in years, to simulate process.  Defaults to 1.
    #     X0 (Optional[numeric]): Initial value for the estimator process, X. Defaults to 0.
    #     convert_y_to_x (Optional[function]): A function to convert Y_t to X_t. Defaults to X_t <- Y_t
    #     convert_x_to_y (Optional[function]): A function to convert X_t to Y_t. Defaults to Y_t <- X_t
    # Raises:
    #     Error: If mu is not a function.
    #     Error: If sigma is not a function.
    # Returns:
    #     (list of function(s)):  Currently returns a list holding a function which runs a simulation
    if (!is.function(mu)) stop('mu must be a function')
    if (!is.function(sigma)) stop('sigma must be a function')

    dt <- T / num_of_steps

    run_monte_carlo <- function(g_maker, strikes, num_paths_to_plot=FALSE){
        # Runs the Monte Carlo simulation to calculate E[g(X_T)] for different strikes.
        #
        # Args:
        #     g_maker (function): A function that returns a function f(x) <- g(x,K)
        #     strikes (array-like): An array of strikes to pass in to g_maker
        #     num_paths_to_plot (Optional[integer or logical]): How many paths should be plotted. Defaults to FALSE.  If TRUE, defaults to 100.
        # Returns:
        #     (array-like): A array of E[g(X_T)] for the range of strikes

        # Simulates the process X_t across time.
        # X is a 2-dimensional matrix with each row representing a path and each entry within it a single timestep.

        dW <- matrix(rnorm(num_of_paths*num_of_steps,sd=sqrt(dt)), nrow=num_of_paths, ncol=num_of_steps)

        # Generate the process Y using the Lamperti transform version of the process X
        Y0 <- convert_x_to_y(X0)
        Y <- cbind(Y0, matrix(nrow=num_of_paths, ncol=num_of_steps))
        for (t in 1:num_of_steps)
            Y[,t+1] <- Y[,t] + mu((t-1)*dt, Y[,t]) * dt + sigma((t-1)*dt, Y[,t]) * dW[,t]

        X <- convert_y_to_x(Y) # Convert the process Y back into the process X

        if (num_paths_to_plot!=FALSE){ # Plots the first num_paths_to_plot paths.
            if (num_paths_to_plot==TRUE) num_paths_to_plot <- 100 # Default value
            num_paths_to_plot <- min(num_paths_to_plot, num_of_paths)
            #matplot plots columns so we transpose
            matplot(t(X[1:num_paths_to_plot,]), main=sprintf('Generated paths (first %i paths) - Euler scheme',num_paths_to_plot), xlab='t', ylab='X_t', type='l')
        }

        X_Ts <- X[,ncol(X)] # Final prices for each path

        # return np.array([np.mean(g_maker(K)(X_Ts)) for K in strikes]) # For reference
        return(sapply(strikes, function(K) mean(g_maker(K)(X_Ts))))
    }
    return(list(run_monte_carlo=run_monte_carlo)) # Return a list holding the function so we can call it OO style
}