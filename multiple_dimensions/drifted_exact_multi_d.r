library(MASS)

source("poisson_process.r")

Drifted_Exact_Multi_D <- function(num_of_paths, beta, mu, sigma0=NA, cov_matrix=NA, T=1, X0=0, convert_y_to_x=function(y) y, convert_x_to_y=function(x) x){
    # Simulate a d-dimensional process using Monte Carlo with the exact method described in section 3 of the reference paper.
    # This process has a general drift process and a constant diffusion coefficient.
    # This can be used to simulate the payoff of an option and calculate the implied volatility smile.
    #
    # Args:
    #     num_of_paths (int): The number of Monte Carlo paths
    #     beta (float): Arrival rate parameter passed to Poisson process
    #     mu (vector of functions): A vector of the drift coefficient functions, a function of t (float) and X (array-like) returning a scalar
    #     sigma0 (vector of floats): A vector of the diffusion coefficient scalars
    #     cov_matrix (matrix of floats): A covariance matrix of the diffusion
    #     T (Optional[float]): Time, in years, to simulate process. Defaults to 1.
    #     X0 (vector of floats): Initial values for each of the process, Xhat. Defaults to 0.
    #     convert_y_to_x (function): A function to convert Y_t to X_t. Defaults to X_t <- Y_t
    #     convert_x_to_y (function): A function to convert X_t to Y_t. Defaults to Y_t <- X_t
    # Returns:
    #     (list of function): A list with a single entry - the run_monte_carlo function stored in the key run_monte_carlo

    # Check that the inputs are correct dimensions
    if (!is.vector(mu)) stop('mu must be a vector')
    num_dim <- length(mu)
    if (!is.matrix(sigma0) && is.na(sigma0)) sigma0 <- diag(num_dim)
    if (num_dim != nrow(sigma0)) stop('mu and sigma0 must be the same length')
    if (!is.matrix(cov_matrix) && is.na(cov_matrix)) cov_matrix <- diag(num_dim)
    if (nrow(cov_matrix) != num_dim) stop('cov_matrix must be the same size as mu and sigma')

    run_monte_carlo <- function(g_maker, strikes){
        # Runs the Monte Carlo simulation to calculate E[g(X_T)] for different strikes.
        #
        # Args:
        #     g_maker (function): A function g(K) that returns a function f(x,K)
        #     strikes (array-like): An array of strikes to pass in to g_maker
        # Returns:
        #     (array-like): A array of E[g(X_T)] for the range of strikes

        malliavin_weight <- function(k){
            # Generate the Malliavin weights. Needed to calculate psi.
            #
            # Args:
            #     k (int): The subscript index to identify which weight to calculate.

            # Calculate the mu difference at between step k-1 and k for each dimension d
            mu_diffs <- sapply(1:num_dim, function(d) mu[[d]](poisson$t[k], Yhat[,k]) - mu[[d]](poisson$t[k-1], Yhat[,k-1]))

            return(mu_diffs %*% solve(sigma0) %*% dW[,k+1] / poisson$dt[k+1])
        }

        psi <- function(g){
            # Calculate psi from the reference paper.
            # E[psi] <- E[g(X_T)] where X_T is the last step of the standard Euler scheme process.
            #
            # Args:
            #     g (function): A function taking X_T and returning a scalar
            # Returns:
            #     float: An estimator in expectation for g(X_T)
            return((g(Xhat[,ncol(Xhat)]) - ifelse(N_T > 0, g(Xhat[,ncol(Xhat)-1]), 0)) * g_multiple)
        }

        # Run Monte Carlo

        # totals keeps track of the sum of the g(X_T) across all paths for each strike
        totals <- rep(0,length(strikes))

        # Calculate psi for each path
        for (path in 1:num_of_paths){
            # Generate the time steps for the Xhat process
            poisson <- Poisson_Process(beta, T)
            #poisson <- data.frame(dt=c(0, 0.7971959, 0.2028041), t=c(0, 0.7971959, 1))
            #print(poisson)
            N_T <- length(poisson$t) - 2 # N_T is the number of arrivals before T (so last arrival and starting time don't count towards N_T)
#~             cat('\nN_T',N_T)
            # Generate dW for each dimension d. dW will be a d-by-length(poisson$t) matrix
            dW <- matrix(NA, nrow=num_dim, ncol=length(poisson$dt))
            for (i in 1:length(poisson$dt)){
                dW[,i] <- mvrnorm(n=1, mu=rep(0, num_dim), Sigma=poisson$dt[i]*cov_matrix)
            }
#~             cat('\n\ndW\n')
#~             print(dW)

            # Generate the process Yhat using the exact method on the Lamperti transform version of the process X
            Y0 <- convert_x_to_y(X0)
#~             cat('\n\nY0\n')
#~             print(Y0)
            Yhat <- cbind(Y0, matrix(NA, nrow=num_dim, ncol=length(poisson$t)-1))
            for (k in 1:(length(poisson$t)-1)){
                for (d in 1:num_dim){
                    Yhat[d,k+1] <- Yhat[d,k] + mu[[d]](poisson$t[k], Yhat[,k])*poisson$dt[k+1] + sigma0[d,d]*dW[d,k+1]
                }
            }
#~             cat('\n\nYhat\n')
#~             print(Yhat)

            # Convert the process Yhat back into the process Xhat
            Xhat <- convert_y_to_x(Yhat)
            
#~             cat('\n\nXhat\n')
#~             print(Xhat)

            # Calculate the product from k=1..N_T of malliavin_weight(k)
            # If there are no arrivals before T, set the product to 1
            product <- ifelse(N_T>0, prod(sapply(2:(1+N_T), malliavin_weight)), 1)
          
#~             cat('\n\nproduct\n')
#~             print(product)
            
            g_multiple <- exp(beta*T) * beta^-N_T * product
#~             cat('\n\ng_multiple\n')
#~             print(g_multiple)

            # For each strike, get the value of psi and add that to totals
            totals <- totals + sapply(strikes, function(K) psi(g_maker(K)))
#~             cat('\n\ntotals\n')
#~             print(totals)
        }
        print(totals/num_of_paths)
        return(totals/num_of_paths)
    }
    return(list(run_monte_carlo=run_monte_carlo))# Return a list holding the function so we can call it oop style
}
