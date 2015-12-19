import numpy as np

"""
Defines several coefficents for the drift and diffusion.
"""

#################################
## General Coefficients
#################################

def make_constant_coefficient(constant_value):
    """
    Creates a function that returns a constant.
    
    Args:
        constant_value (float): The value that the produced function should return
    Returns:
        (function(anything, anything)): A function that will always return constant_value
    """ 
    return lambda t, x: constant_value
    
def make_affine_coefficient(a, b):
    """
    Creates a affine function a*x+b.
    
    Args:
        a (float): The constant term of the affine function
        b (float): The first order term of the affine function
    Returns:
        (function(anything, float)): A function that returns a*x+b
    """ 
    return lambda t, x: a*x + b

#################################
## Drift Coefficients
#################################

def make_mean_reverting_coefficient(mean_reverting_speed, long_term_mean):
    """
    Creates a mean reverting function of the form -mean_reverting_speed(x-long_term_mean)
    
    Args:
        mean_reverting_speed (float): The mean reverting speed
        long_term_mean (float): The value that the process reverts to.
    Returns:
        (function(anything, float)): A function that returns -mean_reverting_speed(x-long_term_mean)
    """ 
    return lambda t, x: -mean_reverting_speed * (x - long_term_mean)

def mu_sin(t, x):
    """
    The drift coefficient. Sin function.
    Args:
        t (float): The current time
        x (np.array): The previous value of the process for each path
    """
    return np.sin(x)
    

#################################
## Diffusion Coefficients
#################################
    
def sigma_from_paper(t, x):
    """
    The diffusion coefficient. Implementing equation 5.1 from the paper.
    Args:
        t (float): The current time
        x (np.array): The previous value of the process for each path
    """
    sigma = 0.4
    return 2*sigma/(1+x**2)
    
def sigma_deriv_from_paper(t, x):
    """
    The derivitive of the diffusion coefficient. Implementing equation 5.1 from the paper.
    Args:
        t (float): The current time
        x (np.array): The previous value of the process for each path
    """
    sigma = 0.4
    return -4*x*sigma/(1+x**2)**2
    
