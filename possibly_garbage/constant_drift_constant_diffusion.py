import numpy as np
from diffusion import Diffusion

class Constant_Drift_Constant_Diffusion(Diffusion):
    """
    Simulates the diffusion process dX_t = mu*dt + sigma * dW_t
    where mu and sigma are constant.
    """
    def __init__(self, num_of_paths, num_of_steps, T=1, X0=0, mu=0, sigma=1):
        """
        Constructor. Calls the parent constructor and sets mu and sigma.
        """
        super(Constant_Drift_Constant_Diffusion, self).__init__(num_of_paths, num_of_steps, T, X0)
        self.mu = mu
        self.sigma = sigma

    def drift(self):
        """
        Returns a vector of the drift. 
        The same drift vector can be used for each simulation since mu_t does not 
        depend on X_t (since mu is constant)
        """
        return self.mu * np.linspace(0, self.T, self.num_of_steps+1)

    def simulate_process(self):
        """
        Simulates all num_of_paths paths. 
        """
        
        # Simulate the diffusion part
        self.simulate_W()
        
        # Add the drift vector to each of the (scaled) paths 
        self.X = self.drift() + self.sigma * self.W

