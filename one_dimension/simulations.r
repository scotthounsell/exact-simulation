source("black_scholes.r")
source("payoffs.r")

# allImpVols is a data frame which should be filled with a strike (K) column and subsequent implied vol columns

run_euler <- function(process, strikes, T, X0, payoff_maker, num_paths_to_plot=TRUE, saveColumn=NA){
    # Runs a Monte Carlo simulation using Euler method to calculate the implied volatility smile.
    #
    # Args:
    #     process (object): An object of one of the exact method classes.
    #     strikes (array-like): Array of strikes used to price the call.
    #     T (numeric): Time, in years, to simulate process.
    #     X0 (numeric): Initial value for the estimator process, X.
    #     saveColumn (character): Column name to save implied vols in global data.frame

    # Simulate the process using Euler scheme and calculate the implied volatility for a range of strikes
    prices <- process$run_monte_carlo(function(K) payoff_maker(K), strikes, num_paths_to_plot)
    implied_vols <- calculate_smile(X0, T, 0, 0, strikes, prices, plot_smile=TRUE, title="Euler Scheme Implied Vol")
    if (!is.na(saveColumn)) allImpVols[saveColumn] <<- implied_vols

    print(data.frame(K=strikes, price=prices, implied.vol=implied_vols)) # Print a table of the results
}

run_exact <- function(process, strikes, T, X0, payoff_maker, num_paths_to_plot=TRUE, saveColumn=NA){
    # Runs a Monte Carlo simulation using Exact method to calculate the implied volatility smile.
    #
    # Args:
    #     process (object): An object of one of the exact method classes.
    #     strikes (array-like): Array of strikes used to price the call.
    #     T (numeric): Time, in years, to simulate process.
    #     X0 (numeric): Initial value for the estimator process, X.
    #     saveColumn (character): Column name to save implied vols in global data.frame

    # Simulate the process using Exact method and calculate the implied volatility for a range of strikes
    prices <- process$run_monte_carlo(function(K) payoff_maker(K), strikes, num_paths_to_plot)
    implied_vols <- calculate_smile(X0, T, 0, 0, strikes, prices, plot_smile=TRUE, title="Exact Scheme Implied Vol")
    if (!is.na(saveColumn)) allImpVols[,saveColumn] <<- implied_vols

    print(data.frame(K=strikes, price=prices, implied.vol=implied_vols)) # Print a table of the results
}