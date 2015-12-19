import numpy as np, matplotlib.pyplot as plt
from poisson_process import Poisson_Process
from operator import mul
from black_scholes import calculate_implied_vol


class Drifted_Exact_1D(object):
    '''Simulates a one-dimensional constant diffusion process without discretization error...hopefully.'''
    def __init__(self, mu, sigma0, beta, T=1., X0=.0):
        '''
        Args:
            mu (function): The drift coefficient, a scalar function of t_i,X_i
            sigma0 (float): The diffusion coefficient, a non-zero scalar
            beta (float): Exponential parameter passed to Poisson process
            T (Optional[float]): Time, in years, to simulate process
            X0 (Optional[float]): Initial value for the estimator process, Xhat
        Raises:
            TypeError: If mu is not a function or sigma0 is not (or is not coercible to) the proper dimension.
            ZeroDivisionError: If sigma0 is zero (non-invertible)
        '''
        if not hasattr(mu,'__call__'): raise TypeError('mu must be a function')
        if sigma0==0: raise ZeroDivisionError('sigma0 must be non-zero')
        self.mu     = mu
        self.sigma0 = float(sigma0)
        self.beta   = float(beta)
        self.T      = float(T)
        self.X0     = float(X0)

        self.poisson = self.dW = self.X = None # Placeholders (X is actually Xhat, an estimator for X)

    def simulate(self):
        '''This function evolves the one-dimensional process, X, across time at steps determined by a Poisson process.'''
        self.poisson = Poisson_Process(self.beta, self.T)
        numSteps     = len(self.poisson.dt)
        self.dW      = np.random.normal(0, np.sqrt(self.poisson.dt)) # For multidim, size=(self.d,numSteps)
        self.X       = np.empty_like(self.poisson.t) # Extra step is the initial value
        self.X[0]    = self.X0

        if numSteps==1: self.dW = [self.dW] # If dW is only a single value, convert it to a list to ensure iterability
        for k,(dt,dW) in enumerate(zip(self.poisson.dt, self.dW)):
            self.X[k+1] = self.X[k] + self.mu(self.poisson.t[k], self.X[k])*dt + self.sigma0*dW

    def Wbar(self,k):
        '''
        Args:
            k (int): The subscript index to identify which weight to calculate.
        Returns:
            float: Malliavin weight
        '''
        return ((self.mu(self.poisson.t[k], self.X[k]) - self.mu(self.poisson.t[k-1], self.X[k-1])) * self.dW[k]/self.sigma0) / self.poisson.dt[k]

    def psi(self,g):
        '''
        Args:
            g (function): A function taking X_i and returning a scalar
        Returns:
            float: An estimator in expectation for g(X_T)
        '''
        N_T = len(self.poisson.dt) - 1
        prod = reduce(mul, map(self.Wbar, range(1,N_T+1)), 1) # product from k=1..N_T of Wbar[k] using initial=1 to protect against empty range
        return np.exp(self.beta*self.T) * (g(self.X[-1]) - g(self.X[-2])*(N_T>0)) * self.beta**-N_T * prod

    def plot(self):
        plt.plot(self.poisson.t, self.X)
        plt.xlabel('t'); plt.ylabel('X_t')


class Drifted_Exact(object):
    '''Simulates a multi-dimensional constant diffusion process without discretization error...hopefully.'''
    def __init__(self, d, mu, sigma0, beta, T=1., X0=None):
        '''
        Args:
            d (int): The dimensionality of the process X.
            mu (function): The drift coefficient, a function of t_i,X_i returning a vector
            sigma0 (np.array or np.matrix): The diffusion coefficient, an invertible dxd matrix
            beta (float): Exponential parameter passed to Poisson process
            T (Optional[float]): Time, in years, to simulate process.  Defaults to 1.
            X0 (Optional[int or float or np.array or np.matrix]): dx1 matrix of initial value(s) for the estimator process, Xhat.
                A scalar will be converted to a vector of that constant.  Defaults to 0 vector.
        Raises:
            TypeError: If mu is not a function or sigma0 is not (or is not coercible to) the proper dimension.
        '''
        if not hasattr(mu,'__call__'): raise TypeError('mu must be a function')
        self.d         = d
        self.mu        = mu
        self.sigma0    = np.asmatrix(float(sigma0) if np.isscalar(sigma0) else sigma0)
        if self.sigma0.shape!=(d,d): raise TypeError('sigma0 must be {0}x{0}'.format(d))
        self.sigma0inv = np.linalg.inv(self.sigma0)
        self.beta      = float(beta)
        self.T         = float(T)
        if X0 is not None:
            if np.isscalar(X0): # If given a single initial value...
                self.X0 = np.asmatrix(np.tile(float(X0), (d,1))) # ...create a column vector of that constant
            else: # Otherwise must be given a numpy iterable of length d
                if X0.size!=d: raise TypeError('If X0 is non-scalar, it must be {}x1'.format(d))
                self.X0 = np.asmatrix(np.reshape(X0, (d,1))) # Ensure X0 is a column vector
        else:
            self.X0 = np.asmatrix(np.zeros((d,1)))

        self.poisson = self.dW = self.X = None # Placeholders (X is actually Xhat, an estimator for X)

    def simulate(self):
        '''This function evolves the d-dimensional process, X, across time at steps determined by a Poisson process.'''
        self.poisson = Poisson_Process(self.beta, self.T)
        numSteps     = len(self.poisson.dt)
        self.dW      = np.asmatrix(np.random.normal(0, np.sqrt(self.poisson.dt), size=(self.d,numSteps)))
        self.X       = np.asmatrix(np.empty((self.d,numSteps+1))) # Extra step is the initial value
        self.X[:,0]  = self.X0

        for k,(dt,dW) in enumerate(zip(self.poisson.dt, self.dW.T)): # Transpose to get columns of dW over time
            self.X[:,k+1] = self.X[:,k] + self.mu(self.poisson.t[k], self.X[:,k])*dt + self.sigma0*dW.T

    def Wbar(self,k):
        '''
        Args:
            k (int): The subscript index to identify which weight to calculate.
        Returns:
            float: Malliavin weight
        '''
        d_mu = self.mu(self.poisson.t[k], self.X[:,k]) - self.mu(self.poisson.t[k-1], self.X[:,k-1]) # Column vector
        return float((d_mu.T * self.sigma0inv*self.dW[:,k]) / self.poisson.dt[k])

    def psi(self,g):
        '''
        Args:
            g (function): A function taking X_i and returning a scalar
        Returns:
            float: An estimator in expectation for g(X_T)
        '''
        N_T = len(self.poisson.dt) - 1
        prod = reduce(mul, map(self.Wbar, range(1,N_T+1)), 1) # product from k=1..N_T of Wbar[k] using initial=1 to protect against empty range
        return np.exp(self.beta*self.T) * (g(self.X[:,-1]) - g(self.X[:,-2])*(N_T>0)) * self.beta**-N_T * prod

    def plot(self):
        '''Plots the rows of X_t'''
        plt.plot(self.poisson.t, self.X.T) # Transpose to plot rows
        plt.xlabel('t'); plt.ylabel('X_t')


if __name__ == '__main__':
    # NOTE: I think this implementation is missing something as we are given the dynamics of Y in terms of X.  Not sure how to simulate.

    K=1.

    #Testing 1-D

    def mu_scalar(t,x): return 2*.4*x / (1. + x*x)**2 # sigma=.4
    

    for K in np.arange(.5,1.5,.1):
        def call_payoff_scalar(x): return float(np.maximum(x-K, 0)) # Takes and returns a scalar

        psiVals = np.empty(200)
        for i in xrange(len(psiVals)):
            process = Drifted_Exact(d=1,mu=mu_scalar,sigma0=1.,beta=.2)
            process.simulate()
            # process.plot()
            psiVals[i] = process.psi(call_payoff_scalar)

        price = np.mean(psiVals)
        print 'K:',K,'Price:',price#,'Imp Vol:',imp_vol(1,K,1,0,0,price) # E[psi]
        # plt.show()
