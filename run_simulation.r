source("black_scholes.r")

run_simulation <- function(process, strikes, Texp, S0, payoff_maker, plot_title=NA, num_paths_to_plot=FALSE, solver='newton', return_ivols=TRUE){
    # Runs a Monte Carlo simulation using Euler method to calculate the implied volatility smile.
    #
    # Args:
    #     process (function instance): An instance of one of the simulation functions.
    #     strikes (array-like): Array of strikes used to price the call.
    #     Texp (numeric): Time, in years, to simulate process.
    #     S0 (numeric): Initial value for the estimator process, X. Used for the implied vol calculation.
    #     payoff_maker (function): A function returning the 'g' functions to apply to to the final process values.
    #     plot_title (Optional[string]): Title for implied volatility plot. Defaults to NA (no plot).
    #     num_paths_to_plot (Optional[integer or logical]): How many paths should be plotted. Defaults to FALSE (no plot). If TRUE, defaults to 100.
    #     solver (Optional[string]): What root finding method to use for calculating the implied volatility. Can be newton, bisection,
    #         or uniroot. Defaults to newton.
    #     return_ivols (Optional[boolean]): If TRUE, calculate and return the implied volatilities.
    #         Otherwise, don't calculate the implied volatility and just return the price. Defaults to TRUE.
    # Returns:
    #     (vector): A vector of the prices or implied volatilities for a range of strikes.

    # Time the function
    start_time <- Sys.time()

    # Simulate the process using Euler scheme and calculate the implied volatility for a range of strikes. Save the results to a dataframe.
    prices <- process$run_monte_carlo(function(K) payoff_maker(K), strikes, num_paths_to_plot)

    # Print the elapsed time
    cat('\n\nTime Elapsed:', Sys.time()-start_time, '\n\n')

    df <- data.frame(K=strikes, Price=prices)

    # If we don't want to calculate the implied volatility, just print and return the prices.
    if (!return_ivols){
        print(df)
        if (!is.na(plot_title)) plot(strikes, prices, main=plot_title, xlab='K', ylab='Price', type='l', col='blue', las=1)
        return(prices)
    }

    # Otherwise, calculate the implied volatility, add it to the dataframe, and print the prices and implied volatilities
    implied_vols <- tryCatch(calculate_smile(S0, Texp, 0, 0, strikes, prices, plot_smile=!is.na(plot_title), plot_title=plot_title, solver=solver), error=function(e) {print(e); NA})
    df['Implied.Vol'] <- implied_vols
    print(df)

    return(implied_vols)
}
