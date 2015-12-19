Poisson_Process <- function(beta, T){
    # Generates a Poisson process with a truncated last step.
    # Constructor. Sets two arrays:
    # dt (array-like): Waiting times between arrivals (exponential random variables).
    #     The last waiting time is shrunk so that its arrival coincides with T.
    #     A 0 is inserted in the beginning so that the indices match the reference paper.
    # t (array-like): Arrival times.
    #     Note that a 0 is inserted to the front of the array.
    #     Also, the last arrival time will always be T.
    #
    # Args:
    #     beta (numeric): Parameter for exponential random variable.
    #     T (numeric): Expiry. The Poisson process is generated up to this point.

    # Generate exponential random variables until their sum adds to more than T
    total <- dt <- 0
    while (total <= T){
        random_exponential <- rexp(1,beta)
        total <- total + random_exponential
        dt <- c(dt,random_exponential)
    }
    # Change the last exponential so that the sum of all exponentials equals T
    dt[length(dt)] <- dt[length(dt)] - (total - T)

    # t is the cumulative sum of the exponentials
    return(data.frame(dt=dt, t=cumsum(dt)))
}

# An example to show how the class works. Keep this commented out unless testing the class.
# print(Poisson_Process(1,3))