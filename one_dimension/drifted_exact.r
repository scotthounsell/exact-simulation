source("poisson_process.r")

Drifted_Exact <- function(num_of_paths, beta, mu, sigma0=1, T=1, X0=0, convert_y_to_x=function(y) y, convert_x_to_y=function(x) x){
    # Simulate a 1-dimensional process using Monte Carlo with the exact method described in section 3 of the reference paper.
    # This process has a general drift process and a constant diffusion coefficient.
    # This can be used to simulate the payoff of an option and calculate the implied volatility smile.
    #
    # Args:
    #     num_of_paths (int): The number of Monte Carlo paths
    #     beta (numeric): Arrival rate parameter passed to Poisson process
    #     mu (function): The drift coefficient, a function of t (numeric) and X (array-like) returning a scalar
    #     sigma0 (Optional[numeric]): The constant diffusion coefficient. Defaults to 1.
    #     T (Optional[numeric]): Time, in years, to simulate process.  Defaults to 1.
    #     X0 (Optional[numeric]): Initial value for the estimator process, Xhat. Defaults to 0.
    #     convert_y_to_x (Optional[function]): A function to convert Y_t to X_t. Defaults to X_t <- Y_t
    # Raises:
    #     Error: If mu is not a function.
    # Returns:
    #     list of function(s):  Currently returns a list holding a function which runs a simulation
    if (!is.function(mu)) stop('mu must be a function')

    run_monte_carlo <- function(g_maker, strikes, num_paths_to_plot=FALSE){
        # Runs the Monte Carlo simulation to calculate E[g(X_T)] for different strikes.
        #
        # Args:
        #     g_maker (function): A function that returns a function f(x) <- g(x,K)
        #     strikes (array-like): An array of strikes to pass in to g_maker
        #     num_paths_to_plot (Optional[integer or logical]): How many paths should be plotted. Defaults to FALSE.  If TRUE, defaults to 100.
        # Returns:
        #     array-like: A array of E[g(X_T)] for the range of strikes

        malliavin_weight <- function(k)
            # Generate the Malliavin weights. Needed to calculate psi.
            #
            # Args:
            #     k (int): The subscript index to identify which weight to calculate.
            return((mu(poisson$t[k], Yhat[k]) - mu(poisson$t[k-1], Yhat[k-1])) * dW[k+1] / (sigma0*poisson$dt[k+1]))

        psi <- function(g)
            # Calculate psi from the reference paper.
            # E[psi] <- E[g(X_T)] where X_T is the last step of the standard Euler scheme process.
            #
            # Args:
            #     g (function): A function taking X_T and returning a scalar
            # Returns:
            #     numeric: An estimator in expectation for g(X_T)
            return((g(Xhat[length(Xhat)]) - ifelse(N_T > 0, g(Xhat[length(Xhat)-1]), 0)) * g_multiple)

        # The simulation starts here
        totals <- rep(0,length(strikes)) # A vector of prices for each strike, summed over paths

        # For plotting
        plotted_paths <- 0
        if (num_paths_to_plot==TRUE) num_paths_to_plot <- 100 # Default value
        num_paths_to_plot <- min(num_paths_to_plot, num_of_paths)
        line_colors <- rainbow(num_paths_to_plot)

        for (path in 1:num_of_paths){
            # Generate timesteps for the process Yhat determined by a Poisson process
            poisson <- Poisson_Process(beta, T)
            N_T <- length(poisson$t) - 2 # N_T is the # of arrivals before T (so last arrival and starting time don't count towards N_T)

            # Generate dW; Including 1st element of dt (zero) ensures a zero in front to match reference paper indices
            dW <- rnorm(length(poisson$dt), sd=sqrt(poisson$dt))

            # Simulate the process Yhat using the exact method on the Lamperti transform version of the process X
            Y0 <- convert_x_to_y(X0)
            Yhat <- c(Y0, rep(NA, length(poisson$t)-1))
            for (k in 1:(length(poisson$t)-1))
                Yhat[k+1] <- Yhat[k] + mu(poisson$t[k], Yhat[k])*poisson$dt[k+1] + sigma0*dW[k+1]

            # Convert the process Yhat back into the process Xhat
            Xhat <- sapply(Yhat, convert_y_to_x)

            # Calculates the product from k=1..N_T of malliavin_weight(k)
            # Using a default of 1 to protect against the edge case of no arrivals before T
            product <- ifelse(N_T>0, prod(sapply(2:(1+N_T), malliavin_weight)), 1)

            # g_multiple represents the portion of psi that is common across g functions of different strikes
            g_multiple <- exp(beta*T) * beta^-N_T * product

            totals <- totals + sapply(strikes, function(K) psi(g_maker(K)))

            if (plotted_paths < num_paths_to_plot){ # Plots the first num_paths_to_plot paths.
                if (plotted_paths==0){ # First plot does labeling
                    plot(poisson$t, Xhat, main='Generated Paths', xlab='t', ylab=expression(hat('X')[t]), ylim=convert_y_to_x(X0)+sqrt(T)*c(-1,1), type='l', col=line_colors[plotted_paths+1])
                } else {
                    lines(poisson$t, Xhat, col=line_colors[plotted_paths+1])
                }
                plotted_paths <- plotted_paths + 1
            }
        }
        return(totals/num_of_paths)
        
    }
    return(list(run_monte_carlo=run_monte_carlo)) # Return a list holding the function so we can call it OO style
}
