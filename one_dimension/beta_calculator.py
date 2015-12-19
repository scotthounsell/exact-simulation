import numpy as np
from numpy import sqrt

def beta_calculator(L, T, mu0, sigma0):
    """
    Calculates the sub-optimal beta for the poisson process.
    
    Args:
        L (float): The Lipschitz/Holder coefficient decribed in Assumption 3.1 of the reference paper
        T (float): The maturity of the asset
        mu0 (np.array): The starting value of mu for each dimension (not sure about this)
        sigma0 (np.array): The starting value of sigma for each dimension (not sure about this)
    Returns:
        (float): The sub-optimal beta
    """
    
    # Get the dimension of the process
    d = sigma0.shape[0]
    
    mu_inf = sqrt(sum(mu**2 for mu in mu0))
    
    if d == 1:
        trace = 1.0 / sigma0[0]**2
    else:
        trace = np.trace(np.linalg.inv(np.dot(sigma0, np.transpose(sigma0))))
        
    L_prime = L**2 * ( 2*(1+mu_inf*sqrt(T))**2 * trace  + 2*(3*d+d*(d-1)) )
    
    return sqrt(L_prime + 0.25*T**2) + 0.5*T


# TEST
L = 0.0
T = 1.0
mu0 = np.array([0.0])
sigma0 = np.array([0.4])

print beta_calculator(L, T, mu0, sigma0)
