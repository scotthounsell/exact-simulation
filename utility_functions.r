# This file contains some utility functions.

plot_values <- function(df, title=NA, ylab='Implied Volatility'){
    # Draw a plot that compares the implied volatility
    #
    # Args:
    #     df (data.frame): A dataframe containing the strikes (in column 'K') and either prices or implied volatilities.
    #     title (Optional[string]): The title for the plot and file name. Defaults to NA (no file saved in this case).
    #     ylab (Optional[string]): The label for the y-axis. Defaults to 'Implied Volatility'.

    # Draw the plot
    matplot(df$K, df[names(df)!='K'], main=title, xlab='K', ylab=ylab, type='l', lty=1)
    legend('topright', inset=.05, legend=names(df)[names(df)!='K'], lty=1, col=1:(ncol(df)-1))

    # Save plot to a pdf file in the "results" folder. Replace spaces in file name with underscores
    if (!is.na(title)) dev.copy2pdf(file=paste('results/', gsub(" ", "_", title), '.pdf', sep=''))
}

rowProds <- function(x, N_T){
    # More efficient row product based on r.789695.n4.nabble.com/matrix-row-product-and-cumulative-product-tp841548.html
    #
    # Args:
    #     x (array-like): A matrix or vector or scalar to do a row product over.
    #     N_T (int): The number of time steps not including 0 and T.
    # Returns:
    #     (scalar or vector): The product of each row.

    # If N_T==1, then there is only one malliavin weight per row.  To avoid mistakenly doing a product over the column, return x as is.
    if (N_T == 1) return(x)

    # If there is only one one row, return the product of that row.
    if(is.vector(x)) return(prod(x))

    # Otherwise, return the product of each row in the matrix.
    y <- x[,1]
    for (j in 2:ncol(x)) y <- y*x[,j]
    return(y)
}

feller <- function(params)
    # Args:
    #     params (list): A list containing the parameters to plug into the Heston model
    # Returns:
    #     (boolean): Do the parameters meet the Feller Condition?
    2*params$lambda*params$vbar > params$eta^2