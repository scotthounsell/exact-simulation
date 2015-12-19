import numpy as np
from numpy import exp
import matplotlib.pyplot as plt
from operator import mul
from poisson_process import Poisson_Process

class Driftless_Exact:
    """
    Simulate a 1-dimensional process using Monte Carlo with the exact method described in section 4 of the reference paper.
    This process has zero drift and a general diffusion coefficient.
    This can be used to simulate the payoff of an option and calculate the implied volatility smile.
    """
    def __init__(self, num_of_paths, beta, sigma, sigma_deriv, T=1, X0=0):
        """
        Args:
            num_of_paths (int): The number of Monte Carlo paths
            beta (float): Arrival rate parameter passed to Poisson process
            sigma (function): The diffusion coefficient, a function of t (float) and X (np.array) returning a scalar
            sigma_deriv (function): The derivative of the diffusion coefficient with resect to X, a function of t (float) and X (np.array) returning a scalar
            T (Optional[float]): Time, in years, to simulate process.  Defaults to 1.
            X0 (Optional[float]): Initial value for the estimator process, Xhat. Defaults to 0.
        Raises:
            TypeError: If sigma is not a function.
            TypeError: If sigma_deriv is not a function.
        """
        if not hasattr(sigma, '__call__'): raise TypeError('sigma must be a function')
        if not hasattr(sigma_deriv, '__call__'): raise TypeError('sigma_deriv must be a function')
        
        self.num_of_paths = num_of_paths
        self.beta = beta
        self.sigma = sigma
        self.sigma_deriv = sigma_deriv
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
            c2 = self.sigma_deriv(self.poisson.t[k], self.Xhat[k])
            if c2 == 0:
                self.Xhat[k+1] = self.Xhat[k] + self.sigma(self.poisson.t[k], self.Xhat[k]) * self.dW[k+1]
            else:
                c1 = self.sigma(self.poisson.t[k], self.Xhat[k]) - c2 * self.Xhat[k]
                self.Xhat[k+1] = -c1/c2 + (c1/c2 + self.Xhat[k]) * exp(-0.5*(c2**2)*self.poisson.dt[k+1] + c2*self.dW[k+1])

    def malliavin_weight(self, k):
        """
        Generate the Malliavin weights. Needed to calculate psi.
        
        Args:
            k (int): The subscript index to identify which weight to calculate.
        """
        a = 0.5*(self.sigma(self.poisson.t[k], self.Xhat[k]))**2
        sigma_tilde = self.sigma(self.poisson.t[k-1], self.Xhat[k-1]) + self.sigma_deriv(self.poisson.t[k-1], self.Xhat[k-1])*(self.Xhat[k]-self.Xhat[k-1])
        a_tilde = 0.5*sigma_tilde**2
        dW = self.dW[k+1]
        dt = self.poisson.dt[k+1]
        return (a-a_tilde) / (2*a) * (-self.sigma_deriv(self.poisson.t[k], self.Xhat[k])*dW/dt + (dW**2 - dt)/dt**2)

        
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
    
    def get_antithetic_malliavin(self):
        """
        Calculates the last malliavin weight used in the antithetic version of psi.
        """
        c2 = self.sigma_deriv(self.poisson.t[-2], self.Xhat[-2])
        if c2 == 0:
            self.X_bar = self.Xhat[-2] - self.sigma(self.poisson.t[-2], self.Xhat[-2]) * self.dW[-2]
        else:
            c1 = self.sigma(self.poisson.t[-2], self.Xhat[-2]) - c2 * self.Xhat[-2]
            self.X_bar = -c1/c2 + (c1/c2 + self.Xhat[-2]) * exp(-0.5*(c2**2)*self.poisson.dt[-1] - c2*self.dW[-1])
            
        a = 0.5*(self.sigma(self.poisson.t[-2], self.Xhat[-2]))**2
        sigma_tilde = self.sigma(self.poisson.t[-3], self.Xhat[-3]) + self.sigma_deriv(self.poisson.t[-3], self.Xhat[-3])*(self.Xhat[-2]-self.Xhat[-3])
        a_tilde = 0.5*sigma_tilde**2
        dW = self.dW[-1]
        dt = self.poisson.dt[-1]
        return (a-a_tilde) / (2*a) * (self.sigma_deriv(self.poisson.t[-2], self.Xhat[-2])*dW/dt + (dW**2 - dt)/dt**2)
    
    def psi_average(self, g):
        """
        Calculate average of psi and the antithetic psi.
        Used to turn psi into finite variance.
        
        Args:
            g (function): A function taking X_T and returning a scalar
        Returns:
            float: An estimator in expectation for g(X_T)
        """

        # N_T is the number of arrivals before T (so the last arrival does not count towards N_T)
        N_T = len(self.poisson.dt) - 2
        return 0.5 * self.psi(g) * (1 + self.get_antithetic_malliavin() / self.malliavin_weight(N_T))
        
        
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
            totals += [self.psi_average(g_maker(K)) for K in strikes]
            if plot_paths: self.plot_paths()
            
        if plot_paths:
            plt.title('Generated paths')
            plt.xlabel('t')
            plt.ylabel('Xhat_t')
            plt.show()
        return totals/self.num_of_paths

