library(parallel)
source("poisson_process.r")

rowProds <- function(x){
    # More efficient row product based on r.789695.n4.nabble.com/matrix-row-product-and-cumulative-product-tp841548.html
    if(is.vector(x)) return(x)
    y <- x[,1]
    for (j in 2:ncol(x)) y <- y*x[,j]
    y
}

Drifted_Exact_Vec <- function(num_of_paths, beta, mu, sigma0=1, Texp=1, X0=0, convert_y_to_x=function(y) y, convert_x_to_y=function(x) x){
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
    # Raises:
    #     Error: If mu is not a function.
    # Returns:
    #     list of function(s):  Currently returns a list holding a function which runs a simulation
    if (!is.function(mu)) stop('mu must be a function')
    Y0 <- convert_x_to_y(X0)
    ebT <- exp(beta*Texp) # No need to recalculate 100k times

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

        # The simulation starts here
        totals <- rep(0,length(strikes)) # A vector of prices for each strike, summed over paths

        # Generate timesteps for the process Yhat determined by a Poisson process
        # Using mclapply from the 'parallel' library for multithreading (overhead might only be outweighed for all_dt, unclear)
        all_dt <- mclapply(1:num_of_paths, function(i) Poisson_Process(beta, Texp))
        all_t <- lapply(all_dt,cumsum)
        all_N_T <- as.numeric(lapply(all_dt,length)) - 2 # N_T is the # of default events (so last arrival and starting time don't count towards N_T)

        for (N_T in unique(all_N_T)){
            # Hold all dt and t paths of the same length in a matrix, by row (time spans columns)
            dt <- matrix(unlist(all_dt[all_N_T==N_T]), nc=N_T+2, byrow=TRUE) # Faster than do.call(rbind, all_dt[all_N_T==N_T])
            t <- matrix(unlist(all_t[all_N_T==N_T]), nc=N_T+2, byrow=TRUE)

            # Generate dW; Including 1st element of dt (zero) ensures a zero in front to match reference paper indices
            dW <- matrix(rnorm(length(dt), sd=sqrt(dt)), ncol=N_T+2)

            # Simulate the process Yhat using the exact method on the Lamperti transformation of the process X
            Yhat <- matrix(Y0, nr=nrow(dt), nc=ncol(dt)) # More efficient to just initialize all values as Y0
            
            for (k in 1:(N_T+1))
                Yhat[,k+1] <- Yhat[,k] + mu(t[,k], Yhat[,k])*dt[,k+1] + sigma0*dW[,k+1]

            # Convert the process Yhat back into the process Xhat (we only need the last two elements at most)
            XhatT <- convert_y_to_x(Yhat[,ncol(Yhat)])
            if(N_T > 0) XhatTNT <- convert_y_to_x(Yhat[,ncol(Yhat)-1])

            # Calculates the product from k=1..N_T of malliavin_weight(k)
            # Defaults to 1 if no arrivals before T
            product <- if(N_T>0) rowProds(malliavin_weight(2:(N_T+1))) else 1

            # g_multiple represents the portion of psi that is common across g functions of different strikes
            g_multiple <- ebT * beta^-N_T * product

            totals <- totals + vapply(strikes, function(K) sum(psi(g_maker(K))), FUN.VALUE=1)
        }
        return(totals/num_of_paths)
    }
    return(list(run_monte_carlo=run_monte_carlo)) # Return a list holding the function so we can call it OO style
}
