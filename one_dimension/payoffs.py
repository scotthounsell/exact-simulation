import numpy as np

"""
Defines several payoff functions.
"""

def make_call_payoff(K):
    """
    Returns a call payoff function with strike K.
    The call payoff works for a scalar value as wel as an np.array of values.
    """
    return lambda x: np.maximum(np.array(x)-K, 0)

def identity_payoff(X_T):
    """
    Identity payoff function. This can be used to reduce the complexity of psi
    in the Exact method for testing purposes.
    
    Args:
        X_T (int or float or list or np.array): The ending value of the process for each path
    Returns:
        (int or float or list or np.array): X_T
    """
    return X_T
