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
    #     X_T (int or float or array-like): The ending value of the process for each path
    # Returns:
    #     (int or float or array-like): X_T
    return(X_T)
    
make_basket_call_payoff <- function(K)
    # Returns a basket call payoff function with strike K.
    # The returned function can be applied to many paths at once
    function(x) {
        # If there is only one asset in the basket (d=1), change x from a vector to a matrix.
        # This will ensure that colMeans works appropriately.
        if (is.vector(x)) x <- matrix(x, nrow=1)
        
        # Calculate the basket call payoff for each path
        return(pmax(colMeans(x)-K, 0))
    }
    
make_basket_call_payoff_exact <- function(K)
    # Returns a basket call payoff function with strike K.
    # The returned function only be applied to a single path.
    function(x) {
        return(max(mean(x)-K, 0))
    }
