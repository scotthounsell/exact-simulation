source("poisson_process.r")

Driftless_Exact <- function(num_of_paths, beta, sigma, sigma_deriv, Texp=1, X0=0){
    # Simulate a 1-dimensional process using Monte Carlo with the exact method described in section 4 of the reference paper.
    # This process has zero drift and a general diffusion coefficient.
    # This can be used to simulate the payoff of an option and calculate the i}mplied volatility smile.
    #
    # Args:
    #     num_of_paths (int): The number of Monte Carlo paths
    #     beta (numeric): Arrival rate parameter passed to Poisson process
    #     sigma (function): The diffusion coefficient, a function of t (numeric) and X (np.array) returning a scalar
    #     sigma_deriv (function): The derivative of the diffusion coefficient with resect to X, a function of t (numeric) and X (np.array) returning a scalar
    #     Texp (Optional[numeric]): Time, in years, to simulate process.  Defaults to 1.
    #     X0 (Optional[numeric]): Initial value for the estimator process, Xhat. Defaults to 0.
    # Raises:
    #     Error: If sigma is not a function.
    #     Error: If sigma_deriv is not a function.
    if (!is.function(sigma)) stop('sigma must be a function')
    if (!is.function(sigma_deriv)) stop('sigma_deriv must be a function')
    ebT <- exp(beta*Texp)

    run_monte_carlo <- function(g_maker, strikes, num_paths_to_plot=FALSE){
        # Runs the Monte Carlo simulation to calculate E[g(X_T)] for different strikes.
        #
        # Args:
        #     g_maker (function): A function that returns a function f(x) <- g(x,K)
        #     strikes (array-like): An array of strikes to pass in to g_maker
        #     num_paths_to_plot (Optional[integer or logical]): How many paths should be plotted. Defaults to FALSE.  If TRUE, defaults to 100.
        # Returns:
        #     (array-like): An array of E[g(X_T)] for the range of strikes

        malliavin_weight <- function(k,antithetic=FALSE){
            # Generate the Malliavin weights. Needed to calculate psi and psi_antithetic.
            #
            # Args:
            #     k (int): The subscript index to identify which weight to calculate.
            a <- .5*sigma(t[k], Xhat[k])^2
            sigma_tilde <- sigma(t[k-1], Xhat[k-1]) + sigma_deriv(t[k-1], Xhat[k-1]) * (Xhat[k] - Xhat[k-1])
            a_tilde <- .5*sigma_tilde*sigma_tilde
            dt <- dt[k+1]
            dW_dt <- dW[k+1]/dt
            sgn <- if(antithetic && k==N_T+1) 1 else -1 # antithetic only differs for last weight at T_NT
            (a-a_tilde) / (2*a) * (sgn*sigma_deriv(t[k], Xhat[k]) * dW_dt + (dW_dt*dW_dt - 1/dt))
        }

        psi <- function(g,antithetic=FALSE)
            # Calculate psi from the reference paper.
            # E[psi] <- E[g(X_T)] where X_T is the last step of the standard Euler scheme process.
            #
            # Args:
            #     g (function): A function taking X_T and returning a scalar
            # Returns:
            #     numeric: An estimator in expectation for g(X_T)
            (g(X_T) - if(N_T>0) g(X_TNT) else 0) * g_multiple

        psi_antithetic <- function(g)
            # Calculate psi^- from the reference paper.
            #
            # Args:
            #     g (function): A function taking X_t and returning a scalar
            # Returns:
            #     numeric: An antithetic variable for psi
            (g(X_bar) - if(N_T>0) g(X_TNT) else 0) * g_multiple * antithetic_adjustment

        psi_average <- function(g)
            # Calculate average of psi and the antithetic psi.
            # Used to turn psi into finite variance.
            #
            # Args:
            #     g (function): A function taking X_T and returning a scalar
            # Returns:
            #     numeric: An estimator in expectation for g(X_T)
            (if(N_T>0) .5 * (psi(g) + psi_antithetic(g)) else psi(g))

        # The simulation starts here
        totals <- rep(0,length(strikes)) # A vector of prices for each strike, summed over paths

        # For plotting
        plotted_paths <- 0
        if (num_paths_to_plot==TRUE) num_paths_to_plot <- 100 # Default value
        num_paths_to_plot <- min(num_paths_to_plot, num_of_paths)
        line_colors <- rainbow(num_paths_to_plot)

        for (path in 1:num_of_paths){
            # Generate timesteps for the process Xhat determined by a Poisson process
            dt <- Poisson_Process(beta, Texp)
            t <- cumsum(dt)
            N_T <- length(t) - 2 # N_T is the # of arrivals before T (so last arrival and starting time don't count towards N_T)

            # Generate dW; Including 1st element of dt (zero) ensures a zero in front to match reference paper indices
            dW <- rnorm(length(dt), sd=sqrt(dt))

            # Simulate the process Xhat across time
            Xhat <- rep(X0, length(t)) # More efficient to just initialize all values as Y0
            for (k in 1:(length(t)-1)){
                c2 <- sigma_deriv(t[k], Xhat[k])
                if (c2 == 0){
                    Xhat[k+1] <- Xhat[k] + sigma(t[k], Xhat[k]) * dW[k+1]
                } else {
                    c1 <- sigma(t[k], Xhat[k]) - c2 * Xhat[k]
                    Xhat[k+1] <- -c1/c2 + (c1/c2 + Xhat[k]) * exp(-.5*(c2^2)*dt[k+1] + c2*dW[k+1])
                }
            }

            X_T <- Xhat[N_T+2] # Save the final value for easy access

            # Calculations for psi_antithetic (independent of g)
            if(N_T>0){
                # For brevity/quick access
                T_NT <- t[N_T+1];    X_TNT <- Xhat[N_T+1]

                # Calculate Xbar for psi_antithetic
                c2 <- sigma_deriv(T_NT, X_TNT)
                if (c2 == 0){
                    X_bar <- X_TNT - sigma(T_NT, X_TNT) * dW[N_T+1]
                } else {
                    c1 <- sigma(T_NT, X_TNT) - c2 * X_TNT
                    temp <- c1/c2
                    X_bar <- -temp + (temp + X_TNT) * exp(-c2*(.5*c2*dt[N_T+2] + dW[N_T+2]))
                }

                # Calculate the adjustment ratio for g_multiple in psi_antithetic (last_antithetic_weight / last_weight)
                antithetic_adjustment <- malliavin_weight(N_T+1,antithetic=TRUE) / malliavin_weight(N_T+1) # N_T + 1 b/c R starts from 1
            }

            # Calculates the product from k=1..N_T of malliavin_weight(k)
            # Using a default of 1 to protect against the edge case of no arrivals before T
            product <- if(N_T>0) prod(malliavin_weight(2:(1+N_T))) else 1

            # g_multiple represents the portion of psi that is common across g functions of different strikes
            g_multiple <- ebT * beta^-N_T * product

            # totals <- totals + vapply(strikes, function(K) psi_average(g_maker(K)), FUN.VALUE=1)
            totals <- totals + psi_average(g_maker(strikes))

            if (plotted_paths < num_paths_to_plot){ # Plots the first num_paths_to_plot paths.
                if (plotted_paths==0){ # First plot does labeling
                    plot(t, Xhat, main='Generated Paths', xlab='t', ylab=expression(hat('X')[t]), ylim=X0+sqrt(Texp)*c(-1,1), type='l', col=line_colors[plotted_paths+1])
                } else {
                    lines(t, Xhat, col=line_colors[plotted_paths+1])
                }
                plotted_paths <- plotted_paths + 1
            }
        }
        return(totals/num_of_paths)
    }
    return(list(run_monte_carlo=run_monte_carlo)) # Return a list holding the function so we can call it OO style
}
