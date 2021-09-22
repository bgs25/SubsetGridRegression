#' Computes grid of predictions for ridge regression
#' 
#' @name order_ridge
#' 
#' @description Computes prediction errors for a grid (over lambda and variable subsets) of ridge regression models
#' 
#' @param x design matrix for training
#' @param z design matrix for testing
#' @param y response for training
#' @param yz response for testing
#' @param var_order For user-specified ordering of variables. Indices start at 0, start with least important variable and end with most. By default order will be induced from scaling of columns in design matrix
#' @param lambda For user-specified sequence of tuning parameter lambda
#' @param nlambda Length of automatically generated sequence of tuning parameters lambda
#' @param lambda.min.ratio Ratio of max/min lambda for automatically generated sequence of tuning parameters lambda
#' @param grid.size Number of subsets of variables for which a solution path will be computed for
#' @param lambda.mult Scales the sequence of lambda by a constant
#' @param errors_mean Controls whether MSPE or SPE (without dividing by sample size) is returned. Used only for within cv_grid_ridge
#' 
#' @return A list of objects:
#' \itemize{
#' \item errors -- a matrix of errors for the models (grid.size times nlambda)
#' \item lambda -- sequence of lambda values used
#' }
#' 
#' @examples 
#' set.seed(1)
#' X = matrix(0, 50, 500)
#' Z = matrix(0, 10, 500)
#' betavec = c(rep(1,5),rep(0,495))
#' X[ , 1:5 ] = matrix(rnorm(250), 50, 5)
#' Z[ , 1:5 ] = matrix(rnorm(50), 10, 5)
#' YZ = as.vector(Z %*% betavec)
#' Y = as.vector(X %*% betavec)
#' Y = Y + rnorm(50)
#' X = X + matrix(rnorm(50*500), 50, 500)
#' mod1 = order_ridge(X, Z, Y, YZ, grid.size = 50)
#' 
#' @export



order_ridge = function( x, z, y, yz, var_order = NULL, lambda = NULL, nlambda = 100, lambda.min.ratio = 1e-5, grid.size = p, lambda.mult = 1e5, errors_mean = TRUE ) {

  p = dim(x)[ 2 ]
  if (is.null(lambda)) {
    lambda_max = max(abs(colMeans(scale_mle(x) * (y - mean(y)))))
    lambda_min = lambda_max * lambda.min.ratio
    lambda = exp(seq(from=log(lambda_max), to=log(lambda_min), length.out=nlambda))
  } else {
    nlambda = length(lambda)
  }
  lambda = lambda * lambda.mult
  xs = scale(x, scale = FALSE)
  # xs = scale_mle(x) # do I need to be careful with the scaling of this?
  svd1 = svd(xs %*% t(xs)) # doing this and just alterning the diagonal entries will allow fast setting of initial matrices
  
  
  
  
  errors = sapply(lambda, function(ll) grid_ridge_lambda(x, z, y, yz, lambda = ll, x1x1inv = svd1$u %*% diag(1/(svd1$d + ll)) %*% t(svd1$v), 
                                                         grid.size = grid.size, var_order = var_order, errors_mean = errors_mean))
  sols = list()
  sols$errors = errors
  sols$lambda = lambda
  attr(sols,"class")<-"order_ridge" 
  return(sols)
}




update_predictions = function( z1x1, z2, x2, x1x1inv, y ) {
  # This is using the woodbury idenity
  x2 = as.matrix(x2)
  z2 = as.matrix(z2)
  delta_p = dim(x2)[ 2 ]
  mat1 = z1x1 - z2 %*% t(x2)
  if ( delta_p == 1 ) {
    mat2 = x1x1inv + ( x1x1inv %*% x2 %*% t(x2) %*% x1x1inv ) / as.double( 1 - t(x2) %*% x1x1inv %*% x2 )
  } else {
    mat2 = x1x1inv + ( x1x1inv %*% x2 %*% solve(diag(delta_p) - t(x2) %*% x1x1inv %*% x2 ) %*% t(x2) %*% x1x1inv )
  }
  predictions = mat1 %*% mat2 %*% y
  return(list(predictions, mat1, mat2))
} 
