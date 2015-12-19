newton <- function(f,fprime,x0,tol_approx=1e-8,tol_consec=1e-8){
    # Solves a given equation for its roots using Newton's Method given its derivative
    # and an initial guess.  Also takes in function value tolerance and consecutive value tolerance.
    #
    # Args:
    #     f (function): equation to be solved
    #     fprime (function): derivative of f
    #     x0 (float): initial guess
    #     tol_approx (float): tolerance from 0
    #     tol_consec (float): tolerance from previous value
    # Returns:
    #     (float): root of f
    x1 <- x0
    x0 <- x0 - 1 # To enter loop
    # tryCatch({
    # print(data.frame(x0=x0, x1=x1, f1=f(x1), first=abs(f(x1))>tol_approx, second=abs(x1-x0)>tol_consec))
    while (abs(f(x1))>tol_approx || abs(x1-x0)>tol_consec){
        x0 <- x1
        x1 <- x1 - f(x0)/fprime(x0)
    }
    # },error=function(e){
    #     print('Hit error')
    #     print(data.frame(x0=x0,x1=x1,f0=f(x0),f1=f(x1),fp=fprime(x0)))
    # })
    return(x1)
}

# bisection <- function(f,a,b,tol_approx=1e-8,tol_consec=1e-8){
#     # Solves a given equation for its roots using Bisection Method given an initial
#     # interval.  Also takes in function value tolerance and consecutive value tolerance.
#     xL <- a;    xR <- b;    xM <- NA
#     cat('\nxl ',xL,' fxl ',f(xL),'\n')
#     cat('\nxr ',xR,' fxr ',f(xR),'\n')
#     if (f(xL)*f(xM) > 0) stop('Initial bounds must evaluate to opposite sign.')
#     while (abs(xR-xL) > tol_consec || max(abs(c(f(xL),f(xR)))) > tol_approx){
#         xM <- .5 * (xL + xR)
#         ifelse(f(xL)*f(xM) < 0, xR <- xM, xL <- xM) # If sign changes in [xL,xM], make xM the new right bound, otherwise left
#         #print(max( abs(f(xL)) , abs(f(xR)) ),'\t',xR-xL)
#     }
#     return(xM)
# }

call_black_scholes <- function(S, K, T, sigma, r, q){
    # Returns the Black-Scholes value of a vanilla call option.
    d1 <- (log(S/K) + (r - q + 0.5*sigma^2)*T) / (sigma*sqrt(T))
    return(S*exp(-q*T)*pnorm(d1) - K*exp(-r*T)*pnorm(d1 - sigma*sqrt(T)))
}

calculate_implied_vol <- function(S, K, T, r, q, price){
    # Returns the Black-Scholes implied volatility.
    vega <- function(sigma){
        # Vega of a vanilla option in the Black-Scholes framework.
        # This is plugged into Newton's method as the first derivative.
        d1 <- (log(S/K) + (r - q + 0.5*sigma^2)*T) / (sigma*sqrt(T))
        return(S*exp(-q*T)*sqrt(T)*dnorm(d1))
    }
    # uniroot(function(sigma) call_black_scholes(S, K, T, sigma, r, q) - price, interval=c(1e-6,6))
    newton(function(sigma) call_black_scholes(S, K, T, sigma, r, q) - price, vega, x0=.25)
    # bisection(function(sigma) call_black_scholes(S, K, T, sigma, r, q) - price, a=1e-4, b=1)
}

calculate_smile <- function(S, T, r, q, strikes, prices, plot_smile=FALSE, title=''){
    # Calculates the implied volatility for a range of strikes.
    # Can plot the implied volatility smile.

    # Args:
    #     S (float): Spot price
    #     T (float): Maturity
    #     r (float): Risk-neutral return rate
    #     q (float): Continuous dividend rate
    #     strikes (array-like): A range of strike prices
    #     prices (array-like): Prices for range of strikes
    #     plot_smile (boolean): Should the function plot the implied volatility smile
    # Returns:
    #     (list): The list of implied volatilities for the range of strike prices
    implied_vols <- sapply(1:length(strikes), function(i) calculate_implied_vol(S, strikes[i], T, r, q, prices[i]))
    if (plot_smile) plot(strikes, implied_vols, main=title, xlab='K', ylab='Implied Volatility', type='l', col='blue')
    return(implied_vols)
}
