import numpy as np
from numpy import exp
import matplotlib.pyplot as plt
from operator import mul
from poisson_process import Poisson_Process

class Drifted_Exact:
    """
    Simulate a 1-dimensional process using Monte Carlo with the exact method described in section 3 of the reference paper.
    This process has a general drift process and a constant diffusion coefficient.
    This can be used to simulate the payoff of an option and calculate the implied volatility smile.
    """
    def __init__(self, num_of_paths, beta, mu, sigma0=1, T=1, X0=0):
        """
        Args:
            num_of_paths (int): The number of Monte Carlo paths
            beta (float): Arrival rate parameter passed to Poisson process
            mu (function): The drift coefficient, a function of t (float) and X (np.array) returning a scalar
            sigma0 (Optional[float]): The constant diffusion coefficient. Defaults to 1
            T (Optional[float]): Time, in years, to simulate process.  Defaults to 1.
            X0 (Optional[float]): Initial value for the estimator process, Xhat. Defaults to 0.
        Raises:
            TypeError: If mu is not a function.
        """
        if not hasattr(mu, '__call__'): raise TypeError('mu must be a function')
        
        self.num_of_paths = num_of_paths
        self.beta = beta
        self.mu = mu
        self.sigma0 = float(sigma0)
        self.T = float(T)
        self.X0 = float(X0)
        
        self.plotted_paths = 0
        
    def simulate_process(self):
        """
        This function simulates the process Xhat across time at steps determined by a Poisson process.
        """
        self.poisson = Poisson_Process(self.beta, self.T)
        self.dW = np.insert(np.random.normal(0, np.sqrt(self.poisson.dt[1:])), 0, 0) # Zero inserted in front to match reference paper indices
        self.Xhat = np.empty_like(self.poisson.t)
        self.Xhat[0] = self.X0
        
        for k in xrange(len(self.poisson.t)-1):
            self.Xhat[k+1] = self.Xhat[k] + self.mu(self.poisson.t[k], self.Xhat[k])*self.poisson.dt[k+1] + self.sigma0*self.dW[k+1]
            
    def malliavin_weight(self, k):
        """
        Generate the Malliavin weights. Needed to calculate psi.
        
        Args:
            k (int): The subscript index to identify which weight to calculate.
        """
        return self.mu(self.poisson.t[k], self.Xhat[k]) - self.mu(self.poisson.t[k-1], self.Xhat[k-1]) * self.dW[k+1] / (self.sigma0*self.poisson.dt[k+1])
        
    def generate_psi_g_multiple(self):
        """
        Calculate the and store the part of psi that does not need to be recalculated for different g functions.
        Saves computation time by doing this.
        """
        # N_T is the number of arrivals before T (so the last arrival and starting time does not count towards N_T)
        N_T = len(self.poisson.t) - 2

        # product from k=1..N_T of self.malliavin_weight(k)
        # Using initializer=1 to protect against the edge case of no arrivals before self.T
        product = reduce(mul, [self.malliavin_weight(k) for k in xrange(1, N_T+1)], 1)
        self.g_multiple = exp(self.beta*self.T) * self.beta**(-N_T) * product
        
        
    def psi(self, g):
        """
        Calculate psi from the reference paper.
        E[psi] = E[g(X_T)] where X_T is the last step of the standard Euler scheme process.
        
        Args:
            g (function): A function taking X_T and returning a scalar
        Returns:
            float: An estimator in expectation for g(X_T)
        """

        # N_T is the number of arrivals before T (so the last arrival and starting time does not count towards N_T)
        N_T = len(self.poisson.t) - 2

        return (g(self.Xhat[-1]) - (g(self.Xhat[-2]) if N_T > 0 else 0)) * self.g_multiple
        
    def plot_paths(self, paths_to_plot=100):
        """
        Plots the first paths_to_plot paths.
        """
        if self.plotted_paths < paths_to_plot:
            self.plotted_paths += 1
            plt.plot(self.poisson.t, self.Xhat)
        
    def run_monte_carlo(self, g_maker, strikes, plot_paths=False):
        """
        Runs the Monte Carlo simulation to calculate E[g(X_T)] for different strikes.
        
        Args:
            g_maker (funcion): A function that returns a function f(x) = g(x,K)
            strikes (np.array): An array of strikes to pass in to g_maker
        Returns:
            np.array: A array of E[g(X_T)] for the range of strikes
        """
        totals = np.zeros(len(strikes))
        for path in xrange(self.num_of_paths):
            self.simulate_process()
            self.generate_psi_g_multiple()
            totals += [self.psi(g_maker(K)) for K in strikes]
            if plot_paths: self.plot_paths()
            
        if plot_paths:
            plt.title('Generated paths')
            plt.xlabel('t')
            plt.ylabel('Xhat_t')
            plt.show()
        return totals/self.num_of_paths
