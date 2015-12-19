import numpy as np
from numpy import exp, sqrt
import matplotlib.pyplot as plt
from poisson_process import Poisson_Process

class Driftless_Exact:
    """
    Simulates a driftless process with using the exact method.
    The goal is to be able to calculate E[psi] = E[g(X_T)] where X_T is the last step of the standard Euler scheme process.
    This avoids discretization errors normally associated with Monte Carlo simulations using the standard Euler scheme.
    """
    def __init__(self, num_of_paths, sigma_func, sigma_func_deriv, beta, T=1.0, X0=0.0):
        """
        Constructor.
        """
        #self.num_of_paths = num_of_paths
        self.sigma_func = sigma_func
        self.sigma_func_deriv = sigma_func_deriv
        self.beta = beta
        self.T = T
        self.X0 = X0

    def simulate_process(self):
        """
        Simulates the process X (called X hat in the paper).
        """

        # Simulate a Poisson Process
        self.poisson = Poisson_Process(self.beta, self.T)
        num_of_steps = len(self.poisson.t)

        # Simulate Brownian Motion with step size equal to the waiting time between arrivals of the Poisson Process
        self.dW = np.random.normal(0, sqrt(self.poisson.dt))
        if type(self.dW) == float: self.dW = [self.dW] # If dW is only a single value, convert it to a list to ensure iterability
        print "dw", self.dW

        # Generate self.X
        self.X = np.empty(num_of_steps)
        self.X[0] = self.X0
        for k in xrange(num_of_steps-1):
            c2 = self.sigma_func_deriv(self.X[k])
            if c2 == 0:
                self.X[k+1] = self.X[k] + self.sigma_func(self.X[k]) * self.dW[k]
            else:
                c1 = self.sigma_func(self.X[k]) - c2 * self.X[k]
                self.X[k+1] = -c1/c2 + (c1/c2 + self.X[k]) * exp(-0.5*(c2**2)*self.poisson.dt[k] + c2*self.dW[k])

        #print self.X

    def malliavin_weight(self, k):
        """
        Generate the Malliavin weights. Needed to calculate psi.
        """
        a = 0.5*(self.sigma_func(self.X[k]))**2
        sigma_tilde = self.sigma_func(self.X[k-1]) + self.sigma_func_deriv(self.X[k-1])*(self.X[k]-self.X[k-1])
        a_tilde = 0.5*sigma_tilde**2
        dW = self.dW[k] # k or k+1
        dt = self.poisson.dt[k] # k or k+1
        return (a-a_tilde) / (2*a) * (-self.sigma_func_deriv(self.X[k])*dW/dt + (dW**2 - dt)/dt**2)

    def psi(self, func):
        """
        Calculate psi.
        E[psi] = E[g(X_T)] where X_T is the last step of the standard Euler scheme process.
        """

        # N_T is the number of arrivals before T (so the last arrival does not count towards N_T)
        N_T = len(self.poisson.dt) - 1

        product = 1
        for k in xrange(1, N_T): # N_T+1 ?
            print k, self.malliavin_weight(k)
            product *= self.malliavin_weight(k)

        return exp(self.beta*self.T) * (func(self.X[-1]) - func(self.X[-2])*(N_T > 0)) * self.beta**(-N_T) * product

    def plot_paths(self):
        """
        Plots all the paths.
        """
        plt.plot(self.poisson.t, self.X, 'o-')
        plt.title('Generated paths')
        plt.xlabel('t')
        plt.ylabel('X_t')
        plt.show()






###########################
# TEST
###########################
#np.random.seed(0)

def sigma_func(x):
    return 2*0.4 / (1+x**2)

def sigma_func_deriv(x):
    return -4*0.4*x / (1+x**2)**2

K=0
#K = np.arange(0.5, 1.5, 0.1)
def g(x):
    return np.maximum(x-K, 0)

beta = 0.2
T = 1
X0 = 1

N = 1
results = np.empty(N)
for i in xrange(N):
    test = Driftless_Exact(1, sigma_func, sigma_func_deriv, beta, T, X0)
    test.simulate_process()
    #test.plot_paths()
    results[i] = test.psi(g)

print "results", results
