#' Computes grid of predictions for ridge regression for a specific value of lambda
#' 
#' @name order_ridge_lambda
#' 
#' @description Computes prediction errors for a grid of ridge regression models for a specified value of lambda
#' 
#' @param x design matrix for training
#' @param z design matrix for testing
#' @param y response for training
#' @param yz response for testing
#' @param var_order For user-specified ordering of variables. Indices start at 0, start with least important variable and end with most. By default order will be induced from scaling of columns in design matrix
#' @param lambda Value of lambda to be used for the calculation
#' @param x1x1inv If this is already computed then it can be passed here. Only for use within grid_ridge
#' @param grid.size Number of subsets of variables for which a solution path will be computed for
#' @param errors_mean Controls whether MSPE or SPE (without dividing by sample size) is returned. Used only for within cv_grid_ridge
#' 
#' @return A vector of errors of length grid.size

order_ridge_lambda = function( x, z, y, yz, var_order = NULL, lambda, x1x1inv = NULL, grid.size = p, errors_mean = TRUE ) {
  # This is for a specific value of lambda, not a path of solutions. However this does appear to be relatively fast
  p = dim(x)[ 2 ]
  n = length(y)
  mu = mean(y)
  y = y - mean(y)
  y = as.vector(y)
  null_dev = mean(y^2)
  
  if (is.null(var_order)) {
    if (is.null(x)) stop("x must be specified")
    var_order = order(apply(x, 2, var)) 
  }
  
  
  # before we scale the design matrix we should store the column means and variances
  x = scale_mle(x)
  col.means = attr(x, "scaled:center")
  col.vars = attr(x, "scaled:scale") # bad naming, this is the square root MLE variance estimate
  
  grid_seq = ( exp(log(p) / ( grid.size - 1 )) )^(seq.int(from = 0, to = grid.size - 1L, by = 1L))
  grid_seq = pmax(floor(grid_seq + 1e-8), rank(grid_seq, ties.method = "first")) # added rounding factor to ensure final term 
  # includes all variables
  # grid_seq = p + 1 - grid_seq
  grid_seq = rev(p - grid_seq)
  # print(grid_seq)
  # print(var_order)
  
  z = scale(z, scale = F)
  z = t(t(z) / col.vars)
  nz = dim(z)[ 1 ]
  
  z1x1 = z %*% t(x)
  if ( is.null(x1x1inv) )  x1x1inv = solve(x %*% t(x) + lambda * diag(n))
  
  predictions = matrix(0, nz, grid.size)
  predictions[ , 1 ] = z1x1 %*% x1x1inv %*% y
  
  if ( grid.size > 1 ) {
    for ( l in 2:grid.size ) {
      # print(var_order[ (grid_seq[ l - 1 ] + 1):grid_seq[ l ] ])
      x2 = x[ , var_order[ (grid_seq[ l - 1 ] + 1):grid_seq[ l ] ], drop = F ]
      z2 = z[ , var_order[ (grid_seq[ l - 1 ] + 1):grid_seq[ l ] ], drop = F ]
      updated_predictions = update_predictions(z1x1, z2, x2, x1x1inv, y)
      predictions[ , l ] = updated_predictions[[ 1 ]]
      z1x1 = updated_predictions[[ 2 ]]
      x1x1inv = updated_predictions[[ 3 ]]
    }
  }
  
  yz = yz - mu
  errors = yz - predictions
  
  if ( errors_mean ) {
    return(apply(errors^2, 2, mean)) # returns the mean squared error for the models
  } else {
    return(apply(errors^2, 2, sum))
  }
  
}
