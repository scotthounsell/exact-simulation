import numpy as np
import matplotlib.pyplot as plt

class General_Euler:
    """
    Simulate a 1-dimensional process using Monte Carlo with the standard Euler scheme.
    This can be used to simulate the payoff of an option and calculate the implied volatility smile.
    """
    def __init__(self, num_of_paths, num_of_steps, mu, sigma, T=1, X0=0):
        """
        Args:
            num_of_paths (int): The number of Monte Carlo paths
            num_of_steps (int): The number of steps in each Monte Carlo path
            mu (function): The drift coefficient, a function of t (float) and X (np.array) returning a scalar
            sigma (function): The diffusion coefficient, a function of t (float) and X (np.array) returning a scalar
            T (Optional[float]): Time, in years, to simulate process.  Defaults to 1.
            X0 (Optional[float]): Initial value for the estimator process, X. Defaults to 0.
        Raises:
            TypeError: If mu is not a function.
            TypeError: If sigma is not a function.
        """
        if not hasattr(mu, '__call__'): raise TypeError('mu must be a function')
        if not hasattr(sigma, '__call__'): raise TypeError('sigma must be a function')
        
        self.num_of_paths = num_of_paths
        self.num_of_steps = num_of_steps
        self.mu = mu
        self.sigma = sigma
        self.T = float(T)
        self.X0 = float(X0)

    def simulate_process(self):
        """
        Simulates the process self.X across time.
        self.X is a 2-dimensional np.array with each row representing a single path.
        Each entry in the row represents a single timestep.
        """
        dt = self.T / self.num_of_steps
        self.dW = np.sqrt(dt)*np.random.normal(size=(self.num_of_paths, self.num_of_steps))
        self.X = np.empty((self.num_of_paths, self.num_of_steps+1))
        self.X[:,0] = self.X0
        for t in xrange(self.num_of_steps):
            self.X[:,t+1] = self.X[:,t] + self.mu(t*dt, self.X[:,t]) * dt + self.sigma(t*dt, self.X[:,t]) * self.dW[:,t]
    
    def run_monte_carlo(self, g_maker, strikes, plot_paths=False, convert_y_to_x_func=lambda y: y):
        """
        Runs the Monte Carlo simulation to calculate E[g(X_T)] for different strikes.
        
        Args:
            g_maker (funcion): A function that returns a function f(x) = g(x,K)
            strikes (np.array): An array of strikes to pass in to g_maker
            plot_paths (Optional[boolean]): Should the paths be plotted. Defaults to False
            convert_y_to_x_func (Optional[function]): A function to convert Y_T to X_T. Defaults to X_T = Y_T
        Returns:
            (np.array): A array of E[g(X_T)] for the range of strikes
        """
        self.simulate_process()
        if plot_paths: self.plot_paths()
        
        X_Ts = [convert_y_to_x_func(X_T) for X_T in self.X[:,-1]]
        
        return np.array([np.mean(g_maker(K)(X_Ts)) for K in strikes])

    def plot_paths(self, paths_to_plot=100):
        """
        Plots the first paths_to_plot paths.
        """
        paths_to_plot = min(paths_to_plot, len(self.X))
        for path in self.X[:paths_to_plot]:
            plt.plot(path)
        plt.title('Generated paths (first {0} paths) - Euler scheme'.format(paths_to_plot))
        plt.xlabel('t')
        plt.ylabel('X_t')
        plt.show()
