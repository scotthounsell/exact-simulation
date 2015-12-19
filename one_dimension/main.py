import numpy as np
import matplotlib.pyplot as plt

from general_euler import General_Euler
from driftless_exact import Driftless_Exact
from drifted_exact import Drifted_Exact

import coefficients
import payoffs
from black_scholes import calculate_smile

def run_euler(process, strikes, T, X0, convert_y_to_x_func=lambda y: y, label=''):
    """
    Runs a Monte Carlo simulation using Euler method to calculate the implied volatility smile.

    Args:
        process (object): An object of one of the exact method classes.
        strikes (np.array): Array of strikes used to price the call.
        T (float): Time, in years, to simulate process.
        X0 (float): Initial value for the estimator process, X.
        convert_y_to_x_func (Optional[function]): A function to convert Y_T to X_T. Defaults to X_T = Y_T
    """

    # Simulate the process using Euler scheme and calculate the implied volatility for a range of strikes
    prices = process.run_monte_carlo(lambda K: payoffs.make_call_payoff(K), strikes, plot_paths=False, convert_y_to_x_func=convert_y_to_x_func)
    implied_vols = calculate_smile(X0, T, 0, 0, strikes, prices, plot_smile=True, label=label)

    # Print a table of the results
    print "K, price, implied vol"
    for K, price, implied_vol in zip(strikes, prices, implied_vols):
        print K, price, implied_vol

def run_exact(process, strikes, T, X0, label=''):
    """
    Runs a Monte Carlo simulation using Exact method to calculate the implied volatility smile.

    Args:
        process (object): An object of one of the exact method classes.
        strikes (np.array): Array of strikes used to price the call.
    """
    # Simulate the process using Exact method and calculate the implied volatility for a range of strikes
    prices = process.run_monte_carlo(lambda K: payoffs.make_call_payoff(K), strikes, plot_paths=False)
    implied_vols = calculate_smile(X0, T, 0, 0, strikes, prices, plot_smile=True, label=label)

    # Print a table of the results
    print "K, price, implied vol"
    for K, price, implied_vol in zip(strikes, prices, implied_vols):
        print K, price, implied_vol

if __name__ == '__main__':
    """
    This program runs 4 different simulations successfully
    
    1) Equation 5.1 from the reference paper using standard Euler scheme.
    2) Equation 5.2 from the reference paper using standard Euler scheme.
    3) Equation 5.1 from the reference paper using standard Exact scheme.
    4) Equation 5.2 from the reference paper using standard Exact scheme.9

    The output for each method is the E[max(X_T-K, 0)] for a range of K, along with the implied volatility in the Euler scheme.
    The program also generates a plot of the implied volatility for each method.
    Note that the volatility smile matches the smile in Figure 1 on page 18 of the reference paper.
    """

    #np.random.seed(0) # For testing purposes

    # Parameters
    num_of_paths_euler = 10000
    num_of_steps_euler = 1000
    T = 1.0
    X0 = 1.0
    strikes = np.arange(0.6, 1.6, 0.1)

    # Run the driftless process from equation 5.1 in the reference paper using Euler scheme
    print "\nEuler Scheme - equation 5.1"
    mu = coefficients.make_constant_coefficient(0)
    sigma = coefficients.sigma_from_paper

    euler51_process = General_Euler(num_of_paths_euler, num_of_steps_euler, mu, sigma, T, X0)
    run_euler(euler51_process, strikes, T, X0, label="Euler - equation 5.1")


    # Run the drifted process from equation 5.2 in the reference paper using Euler scheme
    print "\nEuler Scheme - equation 5.2"
    mu = coefficients.mu_from_paper
    sigma = coefficients.make_constant_coefficient(1)
    Y0 = 0.0
    euler52_process = General_Euler(num_of_paths_euler, num_of_steps_euler, mu, sigma, T, Y0)
    run_euler(euler52_process, strikes, T, X0, coefficients.convert_y_to_x, "Euler - equation 5.2")


    # Parameters
    num_of_paths_exact = 10000
    beta = 0.2

    # Run the driftless process from equation 5.1 in the reference paper using the Exact method
    print "\nExact Method - equation 5.1"
    sigma = coefficients.sigma_from_paper
    sigma_deriv = coefficients.sigma_deriv_from_paper
    exact51_process = Driftless_Exact(num_of_paths_exact, beta, sigma, sigma_deriv, T, X0)
    run_exact(exact51_process, strikes, T, X0, "Exact - equation 5.1")

    # Run the drifted process from equation 5.2 in the reference paper using the Exact method
    print "\nExact Method - equation 5.2"
    mu = coefficients.mu_from_paper
    sigma0 = 1.0
    Y0 = 0.0
    exact52_process = Drifted_Exact(num_of_paths_exact, beta, mu, sigma0, T, Y0, coefficients.convert_y_to_x)
    run_exact(exact52_process, strikes, T, X0, "Exact - equation 5.2")

    plt.legend(loc='upper right')
    plt.savefig('results.png')
    plt.show()



