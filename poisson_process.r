Poisson_Process <- function(beta, Texp){
    # Generates a Poisson process with a truncated last step.
    # Constructor. Sets dt:
    # dt (array-like): Waiting times between arrivals (exponential random variables).
    #     The last waiting time is shrunk so that its arrival coincides with Texp.
    #     A 0 is inserted in the beginning so that the indices match the reference paper.
    #
    # Args:
    #     beta (numeric): Parameter for exponential random variable.
    #     Texp (numeric): Expiry. The Poisson process is generated up to this point.

    # Generate exponential random variables until their sum adds to more than Texp
    total <- dt <- 0
    while (total <= Texp){
        random_exponential <- rexp(1,beta)
        total <- total + random_exponential
        dt <- c(dt,random_exponential)
    }
    # Change the last exponential so that the sum of all exponentials equals Texp
    dt[length(dt)] <- dt[length(dt)] - (total - Texp)

    return(dt)
}

# An example to show how the class works. Keep this commented out unless testing the function.
# print(Poisson_Process(1,3))
