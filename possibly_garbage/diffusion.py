import numpy as np
import matplotlib.pyplot as plt

class Diffusion(object):
    """
    Simulates Standard Brownian Motion with the specified number of steps, paths, and starting location.
    Parent class for all other diffusion processes
    """
    def __init__(self, num_of_paths, num_of_steps, T=1, X0=0):
        """
        Constructor.
        """
        self.num_of_paths = num_of_paths
        self.num_of_steps = num_of_steps
        self.T = float(T) # cast to a float to avoid integer division when calculating dt
        self.X0 = X0

    def simulate_W(self):
        """
        Generate the Standard Brownian Motion from time 0 to T.
        Start the process at X0.
        Each path has self.num_of_steps steps.
        Generate self.num_of_paths paths.
        Each row of the matrix self.W is one such path.
        """
        dt = self.T / self.num_of_steps
        dW = np.sqrt(dt)*np.random.normal(size=(self.num_of_paths, self.num_of_steps))
        self.W = np.cumsum(dW, axis=1)
        self.W = np.insert(self.W, 0, self.X0, axis=1)
        #self.dW = dW

    def get_V0(self, func, X):
        """
        For each path, calculate func(X_T). 
        Return the mean and std of this value across all paths.
        """
        value = func(X[:,-1]) # evaluate the function at the last timestep of each path (X_T)
        V0 = np.mean(value)
        V0_std = np.std(value)
        return V0, V0_std
        
    def plot_paths(self, X):
        """
        Plots all the paths.
        """
        for path in X:
            plt.plot(path)
        plt.title('Generated paths')
        plt.xlabel('t')
        plt.ylabel('X_t')
        plt.show()
        
        
        

