import numpy as np

from general_euler import General_Euler
from driftless_exact import Driftless_Exact
from drifted_exact import Drifted_Exact

import coefficients
import payoffs
from black_scholes import calculate_smile

def run_euler(num_of_paths, num_of_steps, T, X0, strikes):
    """
    Runs a Monte Carlo simulation using Euler method to calculate the implied volatility smile.
    
    Args:
        num_of_paths (int): The number of Monte Carlo paths
        num_of_steps (int): The number of steps in each Monte Carlo path
        T (float): Time, in years, to simulate process.  Defaults to 1.
        X0 (float): Initial value for the estimator process, X. Defaults to 0.
        strikes (np.array): Array of strikes used to price the call.
    """
    # Implement equation 5.1 from the paper
    mu = coefficients.make_constant_coefficient(0)
    sigma = coefficients.sigma_from_paper
    
    # Some other mu and sigma functions to try out
    #sigma = coefficients.make_affine_coefficient(1, 0) # Black-Scholes with sigma = 1
    #mu = coefficients.make_mean_reverting_coefficient(0.8, 1.0)
    #sigma = coefficients.make_constant_coefficient(1)
    
    # Simulate the process using Euler scheme and calculate the implied volatility for a range of strikes
    process = General_Euler(num_of_paths, num_of_steps, mu, sigma, T, X0)
    prices = process.run_monte_carlo(lambda K: payoffs.make_call_payoff(K), strikes, plot_paths=True)
    implied_vols = calculate_smile(X0, T, 0, 0, strikes, prices, plot_smile=True, title="Euler Scheme Implied Vol - equation 5.1")
    
    # Print a table of the results
    print "K, price, implied vol"
    for K, price, implied_vol in zip(strikes, prices, implied_vols):
        print K, price, implied_vol
        
def run_exact(process, strikes):
    """
    Runs a Monte Carlo simulation using Exact method to calculate the implied volatility smile.
    
    Args:
        process (object): An object of one of the exact method classes.
        strikes (np.array): Array of strikes used to price the call.
    """
    # Simulate the process using Exact method and calculate the implied volatility for a range of strikes
    prices = process.run_monte_carlo(lambda K: payoffs.make_call_payoff(K), strikes, plot_paths=True)
    implied_vols = [None]*len(strikes) # Placeholder. Implied vol won't converge to a value unless the prices are reasonable
    
    # Print a table of the results
    print "K, price, implied vol"
    for K, price, implied_vol in zip(strikes, prices, implied_vols):
        print K, price, implied_vol

if __name__ == '__main__':
    """
    This runs a couple of simulations to show where we are up to. 
    
    First, it runs a simulation using the Euler scheme. The process it simulates is equation 5.1 from the reference paper.
    This method works correctly. The code outputs (some of) the paths generated. It also outputs the volatility smile.
    Note that the volatility smile matches the smile in Figure 1 on page 18 of the reference paper.
    
    Next, it runs a simulation of the same process (equation 5.1), this time using the Exact method. 
    Something is not working properly with this. The prices that it outputs are ridiculous.
    It does not converge, and sometimes gives really large numbers or negative numbers which is obviously wrong.
    Furthermore, if the poisson process has less than TODO arrivals, then the code runs into an index error which we are still working on.
    This is why we use so few paths - to avoid this issue. To temporarily resolve this, beta can be made to be small (beta=0.01)
    which will make the probability of less than TODO arrivals tiny. With beta=0.01, the number of paths can be increased.
    However, this still doesn't solve the problem. The Exact method gives wrong values and doesn't converge.
    
    Lastly, the program runs a simulation with general drift and a constant diffusion coefficient of 1.
    This is done in the spirit of equation 5.2 (and section 3). We believe that 5.2 depends on 5.1 so we will not attempt
    the process described in 5.2 until 5.1 works properly. For this example, we just used a mean reverting drift 
    coefficient of the form -b(X_t - x_bar). Once again, this does not work properly. It does not converge
    and sometimes gives very large or negative values.
    
    The output for each method is the E[max(X_T-K, 0)] for a range of K, along with the implied volatility in the Euler scheme.
    Each method also generates a plot of (some of) the simulated paths.
    """
    
    #np.random.seed(0) # For testing purposes
    
    # Parameters
    num_of_paths_euler = 1000
    num_of_steps_euler = 100
    T = 1.0
    X0 = 1.0
    strikes = np.arange(0.7, 1.6, 0.1)
    
    # Run the driftless process from equation 5.1 in the reference paper using Euler scheme
    print "\nEuler Scheme - equation 5.1"
    run_euler(num_of_paths_euler, num_of_steps_euler, T, X0, strikes)
    
    
    # Parameters
    num_of_paths_exact = 10
    beta = 0.2
    
    # Run the driftless process from equation 5.1 in the reference paper using the Exact method
    print "\nExact Method - equation 5.1"
    sigma = coefficients.sigma_from_paper
    sigma_deriv = coefficients.sigma_deriv_from_paper
    driftless_process = Driftless_Exact(num_of_paths_exact, beta, sigma, sigma_deriv, T, X0)
    run_exact(driftless_process, strikes)
    
    
    # Example of a drifted process using the Exact method
    print "\nExact Method - mean reverting process"
    mu = coefficients.make_mean_reverting_coefficient(0.8, 1.0)
    sigma0 = 1
    drifted_process = Drifted_Exact(num_of_paths_exact, beta, mu, sigma0, T, X0)
    run_exact(drifted_process, strikes)



