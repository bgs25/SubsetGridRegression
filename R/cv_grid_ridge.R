#' Computes cross-validated solution for ridge regression
#' 
#' @name cv_grid_ridge
#' 
#' @description Cross-validation wrapper for grid_ridge that computes solutions, selects and fits the optimal model.
#' 
#' @param x Design matrix, n x p. 
#' @param y Vector of responses, length n.
#' @param K Number of folds for cross-validation. Must be at least 2
#' @param var_order For user-specified ordering of variables. Indices start at 0, start with least important variable and end with most. By default order will be induced from scaling of columns in design matrix
#' @param lambda For user-specified sequence of tuning parameter lambda
#' @param nlambda Length of automatically generated sequence of tuning parameters lambda
#' @param lambda.min.ratio Ratio of max/min lambda for automatically generated sequence of tuning parameters lambda
#' @param grid.size Number of subsets of variables for which a solution path will be computed for
#' @param lambda.mult Scales the sequence of lambda by a constant
#' @param fold_assign For user-specified vector of assignment of folds for cross-validation. Must be of the form of integer vector with entries in 1 , ... , K.
#' 
#' @return A list of objects:
#' \itemize{
#' \item mu -- estimated intercept
#' \item beta -- coefficient estimate
#' \item cv -- a matrix of errors for the models (grid.size times nlambda)
#' \item lambda -- sequence of lambda values used
#' }
#' 
#' @examples 
#' set.seed(1)
#' X = matrix(0, 50, 100)
#' betavec = c(rep(1,5),rep(0,95))
#' X[ , 1:5 ] = matrix(rnorm(250), 50, 5)
#' Y = as.vector(X %*% betavec)
#' Y = Y + rnorm(50)
#' X = X + matrix(rnorm(50*100), 50, 100)
#' mod1 = cv_grid_ridge(X, Y, grid.size = 50)
#' 
#' @export

cv_grid_ridge = function( x, y, K = 5, var_order = NULL, lambda = NULL, nlambda = 100L, lambda.min.ratio = ifelse(n<p, 0.01, 0.0001),
                          grid.size = p, lambda.mult = 1e5, fold_assign = NULL ) {
  # First setup the path of lambdas and the fold assignments, as with grid lasso method. Then compute the predictions
  # and fit the final model on all the data for the best model found
  n = length(y)
  cv_err = list()
  p = dim(x)[ 2 ]
  if ( is.null(var_order) ) {
    var_order = order(apply(x, 2, var))
  }
  if ( is.null(fold_assign) ) fold_assign = ceiling(K * (sample(1:n) / n))
  
  for ( k in 1:K ) {
    print(paste0("Computing solutions for fold ", k))
    cv_err[[ k ]] = grid_ridge(x[ fold_assign != k, ], x[ fold_assign == k, ], y[ fold_assign != k ], 
                                              y[ fold_assign == k ], lambda = lambda, nlambda = nlambda, lambda.min.ratio = lambda.min.ratio,
                                              grid.size = grid.size, var_order = var_order, lambda.mult = lambda.mult, errors_mean = FALSE )
    cv_err[[ k ]] = cv_err[[ k ]]$errors
  }


  cv_err = Reduce('+', cv_err) / n

  best_mod = which(cv_err == min(cv_err), arr.ind = T)[ 1, ]
  
  if (is.null(lambda)) {
    lambda_max = max(abs(colMeans(scale_mle(x) * (y - mean(y)))))
    lambda_min = lambda_max * lambda.min.ratio
    lambda = exp(seq(from=log(lambda_max), to=log(lambda_min), length.out=nlambda))
  }
  
  
  grid_seq = ( exp(log(p) / ( grid.size - 1 )) )^(seq.int(from = 0, to = grid.size - 1L, by = 1L))
  grid_seq = pmax(floor(grid_seq + 1e-8), rank(grid_seq, ties.method = "first")) # added rounding factor to ensure final term 
  grid_seq = rev(p - grid_seq)
  
  x = x[ , var_order[ (grid_seq[ best_mod[ 1 ] ] + 1):p ] ]
  lambda_path = lambda
  lambda = lambda[ best_mod[ 2 ] ]

  obj1 = svd(x)
  mu = mean(y)
  y = y - mu
  fit = list()
  fit$mu = mu
  fit$beta = rep(0, p)
  fit$beta[ var_order[ (grid_seq[ best_mod[ 1 ] ] + 1):p ] ] = obj1$v %*% diag(1/(obj1$d^2 + lambda)) %*% diag(obj1$d) %*% t(obj1$u) %*% y
  fit$cv = cv_err
  fit$lambda = lambda_path
  attr(fit,"class")<-"cv_grid_ridge" 
  return(fit)
}