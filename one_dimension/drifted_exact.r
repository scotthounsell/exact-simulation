library(parallel)
source('poisson_process.r')
source('utility_functions.r')

Drifted_Exact <- function(num_of_paths, beta, mu, sigma0=1, Texp=1, X0=0, convert_y_to_x=function(y) y, convert_x_to_y=function(x) x){
    # Simulate a 1-dimensional process using Monte Carlo with the exact method described in section 3 of the reference paper.
    # This process has a general drift process and a constant diffusion coefficient.
    # This can be used to simulate the payoff of an option and calculate the implied volatility smile.
    #
    # Args:
    #     num_of_paths (int): The number of Monte Carlo paths
    #     beta (numeric): Arrival rate parameter passed to Poisson process
    #     mu (function): The drift coefficient, a function of t (numeric) and X (array-like) returning a scalar
    #     sigma0 (Optional[numeric]): The constant diffusion coefficient. Defaults to 1.
    #     Texp (Optional[numeric]): Time, in years, to simulate process.  Defaults to 1.
    #     X0 (Optional[numeric]): Initial value for the estimator process, Xhat. Defaults to 0.
    #     convert_y_to_x (Optional[function]): A function to convert Y_t to X_t. Defaults to X_t <- Y_t
    #     convert_x_to_y (Optional[function]): A function to convert X_t to Y_t. Defaults to Y_t <- X_t
    # Raises:
    #     Error: If mu is not a function.
    # Returns:
    #     list of function(s): Returns a list holding a function which runs a simulation
    if (!is.function(mu)) stop('mu must be a function')
    Y0 <- convert_x_to_y(X0)
    ebT <- exp(beta*Texp)

    run_monte_carlo <- function(g_maker, strikes, num_paths_to_plot=FALSE){
        # Runs the Monte Carlo simulation to calculate E[g(X_T)] for different strikes.
        #
        # Args:
        #     g_maker (function): A function that returns a function f(x) <- g(x,K)
        #     strikes (array-like): An array of strikes to pass in to g_maker
        #     num_paths_to_plot (Optional[integer or logical]): How many paths should be plotted. Defaults to FALSE (no plot). If TRUE, defaults to 100.
        # Returns:
        #     array-like: An array of E[g(X_T)] for the range of strikes

        # Generate the timesteps for each path of Yhat. Timesteps are determined by a Poisson process.
        all_dt <- mclapply(1:num_of_paths, function(i) Poisson_Process(beta, Texp))
        all_t <- lapply(all_dt, cumsum)

        # N_T is the number of default events (so last arrival and starting time do not count towards N_T)
        all_N_T <- as.numeric(lapply(all_dt, length)) - 2

        # For plotting
        plot_num <- 1
        if (num_paths_to_plot==TRUE) num_paths_to_plot <- 100 # Default value
        num_paths_to_plot <- min(num_paths_to_plot, num_of_paths)
        if (num_paths_to_plot) {
            plot_amount <- round(table(all_N_T)/length(all_N_T) * num_paths_to_plot) # Get the number of paths to plot for each unique N_T
            line_colors <- rainbow(num_paths_to_plot)
        }

        sum_prices <- function(N_T){
            # Returns the sum of prices for each strike based on the paths with N_T defaults
            #
            # Args:
            #    N_T (int): The number of defaults to filter the paths on
            # Returns:
            #    (vector): The sum of the prices over all paths of length N_T + 2 for each strike. Same length as "strikes".

            malliavin_weight <- function(k)
                # Generate the Malliavin weights. Needed to calculate psi.
                #
                # Args:
                #     k (int): The subscript index to identify which weight to calculate.
                (mu(t[,k], Yhat[,k]) - mu(t[,k-1], Yhat[,k-1])) * dW[,k+1] / (sigma0*dt[,k+1])

            psi <- function(g)
                # Calculate psi from the reference paper.
                # E[psi] <- E[g(X_T)] where X_T is the last step of the standard Euler scheme process.
                #
                # Args:
                #     g (function): A function taking X_T and returning a scalar
                # Returns:
                #     numeric: An estimator in expectation for g(X_T)
                (g(XhatT) - if(N_T > 0) g(XhatTNT) else 0) * g_multiple

            # Hold all dt and t paths of length N_T+2 in a matrix.
            # Each row of the matrix is a path. Time spans columns.
            dt <- matrix(unlist(all_dt[all_N_T==N_T]), nc=N_T+2, byrow=TRUE)
            t <- matrix(unlist(all_t[all_N_T==N_T]), nc=N_T+2, byrow=TRUE)

            # Generate dW; Including 1st element of dt (zero) ensures a zero in front to match reference paper indices
            dW <- matrix(rnorm(length(dt), sd=sqrt(dt)), ncol=N_T+2)

            # Simulate the process Yhat using the exact method on the Lamperti transformation of the process X

            # Initialize Y0. More efficient to just initialize all values as Y0
            Yhat <- matrix(Y0, nr=nrow(dt), nc=ncol(dt))

            # Step through time
            for (k in 1:(N_T+1))
                Yhat[,k+1] <- Yhat[,k] + mu(t[,k], Yhat[,k])*dt[,k+1] + sigma0*dW[,k+1]

            # Convert the process Yhat back into the process Xhat (we only need the last two elements at most)
            XhatT <- convert_y_to_x(Yhat[,ncol(Yhat)])
            if (N_T > 0) XhatTNT <- convert_y_to_x(Yhat[,ncol(Yhat)-1])

            # Calculates the product from k=1..N_T of malliavin_weight(k)
            # Defaults to 1 if no arrivals before T
            malliavin_product <- if (N_T>0) rowProds(malliavin_weight(2:(N_T+1)), N_T) else 1

            # g_multiple represents the portion of psi that is common across g functions of different strikes
            g_multiple <- ebT * beta^-N_T * malliavin_product

            # Plot paths
            if (num_paths_to_plot){
                plots <- plot_amount[plot_num]
                if (plot_num == 1){
                    matplot(t(t[1:plots,]), t(convert_y_to_x(Yhat[1:plots,])), main=sprintf('Generated paths (first %i paths) - Exact scheme', num_paths_to_plot), xlab='t', ylab='X_t', type='l', lty=1)
                } else {
                    matplot(t(t[1:plots,]), t(convert_y_to_x(Yhat[1:plots,])), type='l', lty=1, add=TRUE)
                }
                plot_num <<- plot_num + 1
            }

            # Returns the sum of psi across all paths of length N_T + 2 for each strike
            #return(vapply(strikes, function(K) sum(psi(g_maker(K))), FUN.VALUE=1))

            # Calculates prices for each path by row, each strike by column
            psi_matrix <- vapply(strikes, function(K) psi(g_maker(K)), FUN.VALUE=numeric(length(XhatT)))

            # The first row is the sum of the prices across paths.
            # The second row is the sum of the squared prices across paths.
            # Each column represents a strike
#~             if (is.vector(psi_matrix)) return(rbind(psi_matrix, psi_matrix*psi_matrix)) #If only one path
#~             return(rbind(colSums(psi_matrix), colSums(psi_matrix*psi_matrix))) # Sum paths for each strike

            # TEST This may be easier to understand, just return the price matrix
            return(vapply(strikes, function(K) psi(g_maker(K)), FUN.VALUE=numeric(length(XhatT))))
        }

        # Returns the mean of psi across all paths for each strike
        #return(rowSums(vapply(sort(unique(all_N_T)), sum_prices, FUN.VALUE=strikes)) / num_of_paths)

        # For each number of defaults, we have a matrix with a sum row (1) and a sum squared row (2) with a column for each strike
#~         sum_and_squared <- Reduce('+', lapply(sort(unique(all_N_T)), sum_prices))
#~         # Print the variance of the prices for each strike
#~         print(data.frame(K=strikes, Variance=1/(num_of_paths-1) *(sum_and_squared[2,] - 1/num_of_paths * (sum_and_squared[1,])^2)))
#~         return(sum_and_squared[1,]/num_of_paths) # Return the mean prices for each strike

        # TEST This may be easier to understand, if we were to just return the prices (psi_matrix)
        price_matrix <- do.call(rbind, lapply(sort(unique(all_N_T)), sum_prices)) # Concatenate the price matrices for each N_T
        print(data.frame(K=strikes, Variance=diag(var(price_matrix)))) # Diagonals of cov mtx are column variances
        return(colMeans(price_matrix))
    }

    # Return a list holding the function so we can call it object oriented style
    return(list(run_monte_carlo=run_monte_carlo))
}
