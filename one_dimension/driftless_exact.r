source("poisson_process.r")

Driftless_Exact <- function(num_of_paths, beta, sigma, sigma_deriv, T=1, X0=0){
    # Simulate a 1-dimensional process using Monte Carlo with the exact method described in section 4 of the reference paper.
    # This process has zero drift and a general diffusion coefficient.
    # This can be used to simulate the payoff of an option and calculate the i}mplied volatility smile.
    #
    # Args:
    #     num_of_paths (int): The number of Monte Carlo paths
    #     beta (numeric): Arrival rate parameter passed to Poisson process
    #     sigma (function): The diffusion coefficient, a function of t (numeric) and X (np.array) returning a scalar
    #     sigma_deriv (function): The derivative of the diffusion coefficient with resect to X, a function of t (numeric) and X (np.array) returning a scalar
    #     T (Optional[numeric]): Time, in years, to simulate process.  Defaults to 1.
    #     X0 (Optional[numeric]): Initial value for the estimator process, Xhat. Defaults to 0.
    # Raises:
    #     Error: If sigma is not a function.
    #     Error: If sigma_deriv is not a function.
    if (!is.function(sigma)) stop('sigma must be a function')
    if (!is.function(sigma_deriv)) stop('sigma_deriv must be a function')

    run_monte_carlo <- function(g_maker, strikes, num_paths_to_plot=FALSE){
        # Runs the Monte Carlo simulation to calculate E[g(X_T)] for different strikes.
        #
        # Args:
        #     g_maker (function): A function that returns a function f(x) <- g(x,K)
        #     strikes (array-like): An array of strikes to pass in to g_maker
        #     num_paths_to_plot (Optional[integer or logical]): How many paths should be plotted. Defaults to FALSE.  If TRUE, defaults to 100.
        # Returns:
        #     (array-like): An array of E[g(X_T)] for the range of strikes

        malliavin_weight <- function(k){
            # Generate the Malliavin weights. Needed to calculate psi.
            #
            # Args:
            #     k (int): The subscript index to identify which weight to calculate.
            a <- .5*sigma(poisson$t[k], Xhat[k])^2
            sigma_tilde <- sigma(poisson$t[k-1], Xhat[k-1]) + sigma_deriv(poisson$t[k-1], Xhat[k-1]) * (Xhat[k] - Xhat[k-1])
            a_tilde <- .5*sigma_tilde^2
            dW <- dW[k+1]
            dt <- poisson$dt[k+1]
            (a-a_tilde) / (2*a) * (-sigma_deriv(poisson$t[k], Xhat[k]) * dW/dt + (dW^2 - dt)/dt^2)
        }

        psi <- function(g)
            # Calculate psi from the reference paper.
            # E[psi] <- E[g(X_T)] where X_T is the last step of the standard Euler scheme process.
            #
            # Args:
            #     g (function): A function taking X_T and returning a scalar
            # Returns:
            #     numeric: An estimator in expectation for g(X_T)
            (g(Xhat[length(Xhat)]) - ifelse(N_T>0, g(Xhat[length(Xhat)-1]), 0)) * g_multiple

        get_antithetic_malliavin <- function(){
            # Calculates the last malliavin weight used in the antithetic version of psi.
            T_NT <- poisson$t[length(poisson$t)-1]
            X_TNT <- Xhat[length(Xhat)-1]

            a <- .5*(sigma(T_NT, X_TNT))^2
            sigma_tilde <- sigma(poisson$t[length(poisson$t)-2], Xhat[length(Xhat)-2]) + sigma_deriv(poisson$t[length(poisson$t)-2], Xhat[length(Xhat)-2]) * (X_TNT - Xhat[length(Xhat)-2])
            a_tilde <- .5*sigma_tilde^2
            dW <- dW[length(dW)]
            dt <- poisson$dt[length(poisson$dt)]
            (a-a_tilde) / (2*a) * (sigma_deriv(T_NT, X_TNT) * dW/dt + (dW^2 - dt)/dt^2) # Sign change for sigma_deriv
        }

        psi_antithetic <- function(g){
            # Calculate psi^- from the reference paper.
            #
            # Args:
            #     g (function): A function taking X_t and returning a scalar
            # Returns:
            #     numeric: An antithetic variable for psi
            T_NT <- poisson$t[length(poisson$t)-1]
            X_TNT <- Xhat[length(Xhat)-1]
            c2 <- sigma_deriv(T_NT, X_TNT)
            if (c2 == 0){
                X_bar <- X_TNT - sigma(T_NT, X_TNT) * dW[length(dW)-1]
            } else {
                c1 <- sigma(T_NT, X_TNT) - c2 * X_TNT
                X_bar <- -c1/c2 + (c1/c2 + X_TNT) * exp(-.5*(c2^2)*poisson$dt[length(poisson$dt)] - c2*dW[length(dW)])
            }
            (g(X_bar) - ifelse(N_T>0, g(X_TNT), 0)) * g_multiple * get_antithetic_malliavin() / malliavin_weight(N_T+1) # N_T + 1 b/c R starts from 1
        }

        psi_average <- function(g)
            # Calculate average of psi and the antithetic psi.
            # Used to turn psi into finite variance.
            #
            # Args:
            #     g (function): A function taking X_T and returning a scalar
            # Returns:
            #     numeric: An estimator in expectation for g(X_T)
            ifelse(N_T>0, .5 * (psi(g) + psi_antithetic(g)), psi(g))

        # The simulation starts here
        totals <- rep(0,length(strikes)) # A vector of prices for each strike, summed over paths

        # For plotting
        plotted_paths <- 0
        if (num_paths_to_plot==TRUE) num_paths_to_plot <- 100 # Default value
        num_paths_to_plot <- min(num_paths_to_plot, num_of_paths)
        line_colors <- rainbow(num_paths_to_plot)

        for (path in 1:num_of_paths){
            # Generate timesteps for the process Xhat determined by a Poisson process
            poisson <- Poisson_Process(beta, T)
            N_T <- length(poisson$t) - 2 # N_T is the # of arrivals before T (so last arrival and starting time don't count towards N_T)

            # Generate dW; Including 1st element of dt (zero) ensures a zero in front to match reference paper indices
            dW <- rnorm(length(poisson$dt), sd=sqrt(poisson$dt))

            # Simulate the process Xhat across time
            Xhat <- c(X0, rep(NA, length(poisson$t)-1))
            for (k in 1:(length(poisson$t)-1)){
                c2 <- sigma_deriv(poisson$t[k], Xhat[k])
                if (c2 == 0){
                    Xhat[k+1] <- Xhat[k] + sigma(poisson$t[k], Xhat[k]) * dW[k+1]
                } else {
                    c1 <- sigma(poisson$t[k], Xhat[k]) - c2 * Xhat[k]
                    Xhat[k+1] <- -c1/c2 + (c1/c2 + Xhat[k]) * exp(-.5*(c2^2)*poisson$dt[k+1] + c2*dW[k+1])
                }
            }

            # Calculates the product from k=1..N_T of malliavin_weight(k)
            # Using a default of 1 to protect against the edge case of no arrivals before T
            product <- ifelse(N_T>0, prod(sapply(2:(1+N_T), malliavin_weight)), 1)

            # g_multiple represents the portion of psi that is common across g functions of different strikes
            g_multiple <- exp(beta*T) * beta^-N_T * product

            totals <- totals + sapply(strikes, function(K) psi_average(g_maker(K)))

            if (plotted_paths < num_paths_to_plot){ # Plots the first num_paths_to_plot paths.
                if (plotted_paths==0){ # First plot does labeling
                    plot(poisson$t, Xhat, main='Generated Paths', xlab='t', ylab=expression(hat('X')[t]), ylim=X0+sqrt(T)*c(-1,1), type='l', col=line_colors[plotted_paths+1])
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