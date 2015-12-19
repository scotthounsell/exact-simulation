from numpy import exp, sqrt, log as ln
from scipy.stats import norm
from scipy.optimize import newton
import matplotlib.pyplot as plt


def call_black_scholes(S, K, T, sigma, r, q):
    """
    Returns the Black-Scholes value of a vanilla call option.
    """
    d1 = (ln(S/K) + (r - q + 0.5*sigma**2)*T) / (sigma*sqrt(T))
    return S*exp(-q*T)*norm.cdf(d1) - K*exp(-r*T)*norm.cdf(d1 - sigma*sqrt(T))

def calculate_implied_vol(S, K, T, r, q, price):
    """
    Returns the Black-Scholes implied volatility.
    """
    def vega(S, K, T, sigma, r, q):
        """
        Vega of a vanilla option in the Black-Scholes framework. 
        This is plugged into Newton's method as the first derivative.
        """
        d1 = (ln(S/K) + (r - q + 0.5*sigma**2)*T) / (sigma*sqrt(T))
        return S * exp(-q*T)*sqrt(T)*norm.pdf(d1)
    
    initial_guess = 0.25
    return newton(lambda sigma: call_black_scholes(S, K, T, sigma, r, q) - price, x0=initial_guess, fprime=lambda sigma: vega(S, K, T, sigma, r, q))\

def calculate_smile(S, T, r, q, strikes, prices, plot_smile=False, label=''):
    """
    Calculates the implied volatility a range of strikes.
    Can plot the implied volatility smile.
    
    Args:
        S (float): Spot price
        T (float): Maturity
        r (float): Risk-neutral return rate
        q (float): Continuous dividend rate
        strikes (np.array): A range of strike prices
        prices (list): Prices for range of strikes
        plot_smile (boolean): Should the function plot the implied volatility smile
    Returns:
        (list): The list of implied volatilities for the range of strike prices
    """
    implied_vols = [calculate_implied_vol(S, K, T, r, q, price) for K, price in zip(strikes, prices)]
      
    if plot_smile:
        plt.plot(strikes, implied_vols, label=label)
        plt.xlabel('K')
        plt.ylabel('Implied Volatility')
        #plt.show()
    
    return implied_vols
