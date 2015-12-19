# Defines several payoff functions.

make_call_payoff <- function(K)
    # Returns a call payoff function with strike K.
    # The call payoff works for a scalar value as well as an array of values.
    function(x) pmax(x-K, 0)

identity_payoff <- function(X_T)
    # Identity payoff function. This can be used to reduce the complexity of psi
    # in the Exact method for testing purposes.
    #
    # Args:
    #     X_T (numeric or array-like): The ending value of the process for each path
    # Returns:
    #     (int or float or array-like): X_T
    return(X_T)