source("black_scholes.r")
source("payoffs.r")

run_simulation <- function(process, strikes, Texp, S0, payoff_maker, plot_title=NA, num_paths_to_plot=FALSE){
    # Runs a Monte Carlo simulation using Euler method to calculate the implied volatility smile.
    #
    # Args:
    #     process (function instance): An instance of one of the simulation functions.
    #     strikes (array-like): Array of strikes used to price the call.
    #     Texp (numeric): Time, in years, to simulate process.
    #     S0 (numeric): Initial value for the estimator process, X. Used for the implied vol calculation.
    #     payoff_maker (function): A function returning the 'g' functions to apply to to the final process values.
    #     plot_title (Optional[string]): Title for implied volatility plot. Defaults to NA (no plot).
    #     num_paths_to_plot (Optional[integer or logical]): How many paths should be plotted. Defaults to FALSE.  If TRUE, defaults to 100.

    # Time the function
    start_time <- Sys.time()

    # Simulate the process using Euler scheme and calculate the implied volatility for a range of strikes
    prices <- process$run_monte_carlo(function(K) payoff_maker(K), strikes, num_paths_to_plot)
    cat('Prices\n', prices)
    implied_vols <- calculate_smile(S0, Texp, 0, 0, strikes, prices, plot_smile=!is.na(plot_title), title=plot_title)

    # Print the elapsed time
    cat('\n\nTime Elapsed:', Sys.time() - start_time, '\n\n')

    # Print a table of the results
    print(data.frame(K=strikes, Price=prices, Implied.Vol=implied_vols))

    return(implied_vols)
}
